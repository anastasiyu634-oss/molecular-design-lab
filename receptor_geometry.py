"""Conservative receptor coordinate handling for exploratory docking.

Never delete failed residues, infer corrected chemistry from an error alone, or
silently report modelled atoms as experimental. This module provides a coherent
altloc selection and a restrained minimization for PDBFixer-built atoms.
"""
from __future__ import annotations

from collections import defaultdict
from io import StringIO
from math import dist, sqrt


def choose_altlocs(pdb_text: str) -> tuple[str, dict]:
    """Select ONE occupancy-supported alternative per residue, preserving coordinates.

    Blank altloc atoms are shared between configurations. In a tie choose 'A'
    before other labels; PDB remains unchanged if there are no alternatives.
    """
    scores: dict[tuple[str, str, str], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for line in pdb_text.splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')) or len(line) < 60:
            continue
        alt = line[16:17].strip()
        if not alt:
            continue
        key = (line[21:22], line[22:26], line[26:27])
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = 0.0
        scores[key][alt] += max(0., occupancy)
    chosen = {key: sorted(opts, key=lambda label: (-opts[label], label != 'A', label))[0]
              for key, opts in scores.items()}
    removed = 0
    output = []
    for line in pdb_text.splitlines(keepends=True):
        if line.startswith(('ATOM  ', 'HETATM')) and len(line) >= 27:
            alt = line[16:17].strip()
            if alt:
                key = (line[21:22], line[22:26], line[26:27])
                if alt != chosen.get(key):
                    removed += 1
                    continue
                line = line[:16] + ' ' + line[17:]
        output.append(line)
    report = {
        'altloc_choices': [f'{key[0].strip() or "?"}:{key[1].strip()}{key[2].strip()}={alt}'
                           for key, alt in sorted(chosen.items())],
        'excluded_altloc_atoms': removed,
        'coordinate_changes': False,
    }
    return ''.join(output), report


def _atom_key(atom):
    res = atom.residue
    # Identity must include the ATOM NAME, not just the residue. Without this,
    # all atoms within one residue were incorrectly restrained to the same
    # coordinate (the last atom inserted into the observation dictionary).
    return (res.chain.id, str(res.id), getattr(res, 'insertionCode', '') or '',
            res.name, atom.name.strip())


def _heavy_atom_coordinates(pdb_text: str) -> dict:
    """Read observed heavy atoms by stable PDB residue IDs."""
    from openmm.app import PDBFile
    from openmm import unit
    source = PDBFile(StringIO(pdb_text))
    positions = source.positions.value_in_unit(unit.nanometer)
    seen = {}
    for atom in source.topology.atoms():
        if not atom.element or atom.element.symbol == 'H':
            continue
        key = _atom_key(atom)
        if key in seen:
            raise ValueError(f'В исходном PDB повторяется идентификатор атома: {key}. '
                             'Проверь альтернативные положения и нумерацию.')
        seen[key] = tuple(float(x) for x in positions[atom.index])
    return seen


def _anchor_plan(atoms, observed: dict) -> tuple[dict[int, tuple[float, float, float]], list[int]]:
    """Determine *observed atom* indices that must remain fixed during repair.

    Unlike per-residue matching, this uses chain, residue ID, insertion code,
    residue name and atom name.  Missing observed atoms are a hard failure.
    """
    anchored = {}
    matched = set()
    new_heavy = []
    for atom in atoms:
        if atom.element is None or atom.element.symbol == 'H':
            continue
        key = _atom_key(atom)
        if key in observed:
            if key in matched:
                raise ValueError(f'При восстановлении один атом встречается дважды: {key}')
            matched.add(key)
            anchored[atom.index] = observed[key]
        else:
            new_heavy.append(atom.index)
    lost = set(observed).difference(matched)
    if lost:
        raise ValueError('При восстановлении потеряны исходные тяжёлые атомы; подготовка остановлена.')
    return anchored, new_heavy


def minimize_repaired_receptor(repaired_pdb: str, observed_pdb: str, *,
                               max_iterations: int = 250) -> tuple[str, dict]:
    """Minimize *only reconstructed atoms and hydrogens*, preserving observed heavy atoms.

    Existing experimental heavy atoms are truly fixed (zero-mass particles in
    OpenMM), not harmonically restrained.  Corrects the previous issue where
    soft positional restraints let existing coordinates drift >0.5 A.

    This is still a provisional repaired structure, not a validated receptor.
    Failure to parameterize, optimize, or pass clash checks is fatal.
    """
    from openmm import LangevinMiddleIntegrator, Vec3
    from openmm.app import ForceField, PDBFile, Simulation, NoCutoff
    from openmm import unit

    observed = _heavy_atom_coordinates(observed_pdb)
    if not observed:
        raise ValueError('Нет наблюдаемых тяжёлых атомов для контроля геометрии.')
    pdb = PDBFile(StringIO(repaired_pdb))
    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    # Do not create bond-length constraints here.  Original heavy atoms are
    # frozen below by giving them zero mass; an OpenMM constraint involving a
    # zero-mass particle causes Context initialization to fail.  For this
    # minimization we need unconstrained bonds, including H bonds and water.
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff,
                                     constraints=None, rigidWater=False)
    if system.getNumConstraints() != 0:
        raise ValueError('В модели остались ограничения связей OpenMM, несовместимые '
                         'с фиксацией экспериментальных атомов. Оптимизация отменена.')
    atoms = list(pdb.topology.atoms())
    anchors, reconstructed = _anchor_plan(atoms, observed)
    initial = list(pdb.positions.value_in_unit(unit.nanometer))
    for idx, (x, y, z) in anchors.items():
        # OpenMM explicitly treats zero-mass particles as fixed during minimization.
        system.setParticleMass(idx, 0 * unit.dalton)
        # PDBFixer may have slightly moved an observed atom even before OpenMM.
        # Restore its exact observed coordinate before freezing it.
        initial[idx] = Vec3(float(x), float(y), float(z))

    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond,
                                         0.002*unit.picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(unit.Quantity(initial, unit.nanometer))
    try:
        simulation.minimizeEnergy(maxIterations=max_iterations)
        final_positions = simulation.context.getState(getPositions=True).getPositions()
        final = final_positions.value_in_unit(unit.nanometer)
    finally:
        del simulation, integrator

    displacements = [dist(tuple(float(v) for v in final[idx]), coords)
                     for idx, coords in anchors.items()]
    max_d = max(displacements, default=0.) * 10.0
    rms_d = sqrt(sum(x*x for x in displacements) / max(len(displacements), 1)) * 10.0
    # This small tolerance verifies actual immobility, not a relaxed "good enough" threshold.
    if max_d > 0.01 or rms_d > 0.005:
        raise ValueError(f'Зафиксированные исходные атомы неожиданно сместились: '
                         f'максимум {max_d:.3f} Å, RMSD {rms_d:.3f} Å.')

    # New heavy atoms must not overlap observed atoms from other residues.
    heavy_original = [a for a in atoms if a.index in anchors]
    n_clashes = 0
    for idx in reconstructed:
        atom = atoms[idx]
        p = tuple(float(x) for x in final[idx])
        for other in heavy_original:
            if atom.residue == other.residue:
                continue
            if (atom.residue.chain.id == other.residue.chain.id and
                    abs(int(atom.residue.id) - int(other.residue.id)) <= 1):
                continue
            if dist(p, tuple(float(x) for x in final[other.index])) < 0.125:
                n_clashes += 1
                break
    if n_clashes:
        raise ValueError(f'После оптимизации осталось {n_clashes} столкновений '
                         'новых тяжёлых атомов с исходными (менее 1,25 Å).')
    out = StringIO()
    PDBFile.writeFile(pdb.topology, final_positions, out, keepIds=True)
    return out.getvalue(), {
        'restrained_minimization': False,
        'fixed_experimental_heavy_atoms': True,
        'fixed_original_heavy_atoms': len(anchors),
        'new_heavy_atoms': len(reconstructed),
        'original_heavy_rmsd_a': round(rms_d, 3),
        'original_heavy_max_shift_a': round(max_d, 3),
        'remaining_short_contacts': n_clashes,
        'validated_experimentally': False,
    }
