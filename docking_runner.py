"""Optional AutoDock Vina workflow for Molecular Design Lab.

The routine is exploratory, not a validated docking protocol.  It never
labels predicted affinities as experimental binding energies.
"""
from __future__ import annotations

import math
import os
import re
from pathlib import Path
from io import StringIO
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from typing import Any


def binding_site_box(pdb_id: str, cif_text: str, ligand: dict,
                     padding: float = 8.0, min_edge: float = 18.0,
                     max_edge: float = 35.0) -> dict[str, list[float]]:
    """Define a search box from the observed ligand in the selected complex.

    This is suitable as a starting point for *focused* docking only.  Users can
    edit center and dimensions before performing the calculation.
    """
    from analysis import load_model

    model = load_model(pdb_id, cif_text)
    coords: list[tuple[float, float, float]] = []
    for chain in model:
        if chain.id != ligand.get('Цепь'):
            continue
        for residue in chain:
            if (residue.resname != ligand.get('Код компонента') or
                    residue.id[1] != int(ligand.get('Номер')) or
                    str(residue.id[2]).strip() != ligand.get('Вставка', '')):
                continue
            if not str(residue.id[0]).startswith('H_'):
                continue
            for atom in residue:
                if atom.element.upper() not in ('H', 'D'):
                    coords.append(tuple(float(v) for v in atom.coord))
    if not coords:
        raise ValueError('Нет координат выбранного лиганда для автоматического определения области связывания.')
    if not 3 <= padding <= 15:
        raise ValueError('Отступ области поиска должен быть от 3 до 15 Å.')
    lows = [min(row[i] for row in coords) for i in range(3)]
    highs = [max(row[i] for row in coords) for i in range(3)]
    center = [round((lo+hi)/2, 2) for lo, hi in zip(lows, highs)]
    size = [round(max(min_edge, hi-lo+2*padding),2) for lo,hi in zip(lows,highs)]
    if any(edge > max_edge for edge in size):
        raise ValueError('Исходный лиганд слишком велик для стандартной области поиска. Нужна ручная настройка или другой комплекс.')
    return {'center':center, 'size':size}


def validate_box(center: list[float], size: list[float]):
    if len(center) != 3 or len(size) != 3:
        raise ValueError('Область поиска задаётся тремя координатами центра и тремя размерами.')
    if not all(math.isfinite(float(x)) for x in center+size):
        raise ValueError('Все значения координат и размеров должны быть конечными числами.')
    if not all(10 <= float(x) <= 40 for x in size):
        raise ValueError('Для первого запуска размеры каждой стороны области должны быть от 10 до 40 Å.')


def ligand_pdbqt(sdf_text: str) -> str:
    """Parameterize RDKit 3D, explicit-H ligand via Meeko."""
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromMolBlock(sdf_text.split('$$$$',1)[0], removeHs=False)
    if mol is None or mol.GetNumConformers() < 1:
        raise ValueError('Сначала построй 3D-конформер для корректного SMILES.')
    if not any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        raise ValueError('В SDF отсутствуют явные атомы водорода. Перестрой конформер.')
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if len(setups) != 1:
        raise RuntimeError('Meeko не подготовил единственный вариант лиганда.')
    output = PDBQTWriterLegacy.write_string(setups[0])
    if isinstance(output, tuple):
        pdbqt, ok, err = output[0], output[1], output[2] if len(output)>2 else ''
        if not ok:
            raise ValueError(f'Meeko не смог подготовить лиганд: {err}')
    else:
        pdbqt = output
    if not isinstance(pdbqt,str) or 'ATOM' not in pdbqt:
        raise ValueError('Meeko не вернул атомы лиганда в PDBQT.')
    return pdbqt


def _receptor_program() -> list[str]:
    """Resolve Meeko CLI from the exact Python environment running Streamlit."""
    scripts = Path(sysconfig.get_path('scripts'))
    for filename in ('mk_prepare_receptor.exe', 'mk_prepare_receptor.py', 'mk_prepare_receptor'):
        path = scripts / filename
        if path.exists():
            return ([sys.executable, str(path)] if path.suffix == '.py' else [str(path)])
    for name in ('mk_prepare_receptor.py','mk_prepare_receptor'):
        tool=shutil.which(name)
        if tool:
            return [sys.executable,tool] if tool.endswith('.py') else [tool]
    raise RuntimeError('Не найдена команда Meeko mk_prepare_receptor. Установи requirements-docking.txt в том же окружении, где запущен сайт.')



class ReceptorPreparationError(RuntimeError):
    """Краткая научно понятная ошибка + отдельные подробности Meeko."""
    def __init__(self, message: str, details: str = '', residues: list[str] | None = None):
        super().__init__(message)
        self.details = details
        self.residues = residues or []


def diagnose_meeko_failure(output: str) -> dict[str, object]:
    """Extract grounded issues from a real Meeko log; never guess structural fixes."""
    output = (output or '').strip()
    positions = re.findall(r"(?:residue[_ ]key\s*[=:]?\s*['\"]?|ключу\s*)([A-Za-z0-9]:-?\d+[A-Za-z]?)", output, re.I)
    positions += re.findall(r'\b([A-Za-z0-9]:-?\d+[A-Za-z]?)\b', output)
    positions = list(dict.fromkeys(positions))[:16]
    valence = 'AtomValenceException' in output or 'valence' in output.lower() or 'валентност' in output.lower()
    if valence:
        issue = ('Meeko/RDKit обнаружил недопустимую валентность: программа считает, '
                 'что у одного из атомов слишком много химических связей.')
    elif 'template' in output.lower() or 'шаблон' in output.lower():
        if 'heavy_miss=' in output:
            issue = 'У некоторых аминокислот отсутствуют тяжёлые атомы, необходимые для химического шаблона.'
        else:
            issue = 'Meeko не смог сопоставить один или несколько остатков с химическими шаблонами.'
        if 'alternate location' in output.lower() or 'altloc' in output.lower():
            issue += ' Также обнаружены альтернативные положения атомов.'
    else:
        issue = 'Meeko не смог подготовить химически корректную модель рецептора.'
    if positions:
        issue += ' Проверьте остатки: ' + ', '.join(positions) + '.'
    return {'summary': issue, 'residues': positions, 'valence': valence}


def inspect_receptor_pdb(pdb_text: str) -> dict[str, object]:
    """Lightweight, conservative check of backbone chain continuity; no repair."""
    from collections import OrderedDict
    chains = OrderedDict()
    for line in (pdb_text or '').splitlines():
        if not line.startswith('ATOM  ') or len(line) < 54:
            continue
        try:
            chain = line[21:22].strip() or '?'
            ident = (int(line[22:26]), line[26].strip())
            atom = line[12:16].strip()
            xyz = tuple(float(line[i:i+8]) for i in (30,38,46))
            name = line[17:20].strip()
        except (ValueError,IndexError):
            continue
        residue = chains.setdefault(chain,OrderedDict()).setdefault(ident,{'name':name,'atoms':{}})
        # Keep first alt-loc; ambiguous alternatives need separate validation.
        residue['atoms'].setdefault(atom, xyz)
    breaks=[]; missing=[]; gaps=[]
    for chain,residues in chains.items():
        data=list(residues.items())
        for (num,ins), residue in data:
            absent=[x for x in ('N','CA','C','O') if x not in residue['atoms']]
            if absent:
                missing.append(f'{chain}:{num}{ins} ({residue["name"]}; отсутствуют {",".join(absent)})')
        for (prev_id,prev),(cur_id,cur) in zip(data,data[1:]):
            if cur_id[0]-prev_id[0] > 1:
                gaps.append(f'{chain}:{prev_id[0]}–{cur_id[0]}')
            ca, na = prev['atoms'].get('C'), cur['atoms'].get('N')
            if ca and na and cur_id[0]-prev_id[0] in (0,1):
                dist=math.dist(ca,na)
                if dist > 2.1 or dist < 1.0:
                    breaks.append(f'{chain}:{prev_id[0]}→{cur_id[0]} (C–N {dist:.2f} Å)')
    return {'chains':len(chains), 'residues':sum(len(x) for x in chains.values()),
            'missing_backbone': missing[:20], 'sequence_gaps':gaps[:20],
            'suspect_peptide_geometry':breaks[:20],
            'terminal_residues':[
                f'{ch}:{next(reversed(rows))[0]} ({next(reversed(rows.values()))["name"]})'
                for ch,rows in chains.items() if rows]}


def prepare_receptor_pdbqt(protein_pdb: str, directory: Path, parser: str = 'rdkit') -> str:
    """Try Meeko; alternative ProDy parser only on explicit user selection."""
    if parser not in {'rdkit','prody'}:
        raise ValueError('Неизвестный способ чтения рецептора.')
    if parser == 'prody':
        import importlib.util
        if importlib.util.find_spec('prody') is None:
            raise RuntimeError('Для альтернативного чтения установите ProDy в том же Python: '
                               r'.\.venv\Scripts\python.exe -m pip install prody')
    directory.mkdir(parents=True, exist_ok=True)
    source = directory / 'protein.pdb'
    target = directory / 'receptor.pdbqt'
    source.write_text(protein_pdb,encoding='utf-8')
    option = '--read_with_prody' if parser == 'prody' else '--read_pdb'
    command = _receptor_program() + [option, str(source), '--write_pdbqt', str(target)]
    try:
        completed = subprocess.run(command,capture_output=True,text=True,timeout=120,
                                   cwd=directory,encoding='utf-8',errors='replace')
    except subprocess.TimeoutExpired as exc:
        raise ReceptorPreparationError('Meeko превысил ограничение времени подготовки (120 секунд).') from exc
    if completed.returncode != 0:
        details=((completed.stderr or '') + '\n' + (completed.stdout or '')).strip()
        issue=diagnose_meeko_failure(details)
        raise ReceptorPreparationError(str(issue['summary']),details,issue['residues'])
    candidates=[target, directory/'receptor_rigid.pdbqt']
    candidates += list(directory.glob('*rigid*.pdbqt'))
    for path in candidates:
        if path.exists() and path.stat().st_size > 80:
            data=path.read_text(encoding='utf-8',errors='replace')
            if any(x.startswith(('ATOM ','HETATM')) for x in data.splitlines()):
                return data
    raise ReceptorPreparationError('Meeko завершил работу, но не создал читаемый рецептор PDBQT.')


def repair_receptor_pdb(protein_pdb: str, directory: Path, ph: float = 7.0) -> tuple[str, dict[str, int]]:
    """Explicit, opt-in preliminary repair using PDBFixer.

    Only complete *existing* residues.  Do not fabricate absent loop residues.
    Return a report so restored coordinates are never mistaken for experimental data.
    The result is NOT automatically a validated receptor for production docking.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise RuntimeError(
            'Для восстановления рецептора установи PDBFixer и OpenMM в окружение сайта: '
            r'.\.venv\Scripts\python.exe -m pip install "pdbfixer>=1.11,<2"'
        ) from exc
    directory.mkdir(parents=True, exist_ok=True)
    inp = directory / 'unrepaired_receptor.pdb'
    out = directory / 'repaired_receptor.pdb'
    inp.write_text(protein_pdb, encoding='utf-8')
    try:
        fixer = PDBFixer(filename=str(inp))
        fixer.findMissingResidues()
        # Only fix atoms in existing residues. Filling missing internal loops or
        # missing termini requires separate structural validation.
        missing_residues = sum(len(x) for x in fixer.missingResidues.values())
        fixer.missingResidues = {}
        fixer.findNonstandardResidues()
        nonstandard = list(fixer.nonstandardResidues)
        if nonstandard:
            names = ', '.join(f'{res.chain.id}:{res.id}({res.name})' for res, _ in nonstandard[:5])
            raise RuntimeError(
                f'Нестандартные остатки ({names}) требуют ручной проверки; '
                'автоматическая замена аминокислот отключена.')
        fixer.findMissingAtoms()
        missing_atoms = sum(len(value) for value in fixer.missingAtoms.values())
        missing_term_atoms = sum(len(value) for value in fixer.missingTerminals.values())
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=float(ph))
        # PDBFixer may put rebuilt side chains on top of existing atoms.
        # Meeko then mistakenly infers an extra covalent bond and RDKit reports
        # "carbon valence 5". Relax ONLY rebuilt atoms/hydrogens freely; restrain
        # all observed heavy atoms and check their displacement afterwards.
        raw = StringIO()
        PDBFile.writeFile(fixer.topology, fixer.positions, raw, keepIds=True)
        from receptor_geometry import minimize_repaired_receptor
        optimized, geometry_report = minimize_repaired_receptor(raw.getvalue(), protein_pdb)
        out.write_text(optimized, encoding='utf-8')
    except RuntimeError:
        raise
    except Exception as exc:
        raise RuntimeError(f'Не удалось достроить и геометрически проверить рецептор: '
                           f'{type(exc).__name__}: {exc}') from exc
    result = out.read_text(encoding='utf-8')
    if not any(line.startswith('ATOM  ') for line in result.splitlines()):
        raise RuntimeError('PDBFixer не создал координаты атомов белка.')
    return result, {
        'added_heavy_atoms': missing_atoms,
        'missing_terminal_atoms': missing_term_atoms,
        'unrestored_missing_residues': missing_residues,
        **geometry_report,
    }

def read_pose_first(pdbqt: str) -> str:
    """Extract top pose PDBQT for 3Dmol; leave full multiple-pose output intact."""
    if 'MODEL' in pdbqt:
        lines = pdbqt.splitlines(keepends=True)
        started = False
        result=[]
        for line in lines:
            if line.startswith('MODEL'):
                if started:break
                started=True
            if started:
                result.append(line)
                if line.startswith('ENDMDL'):break
        return ''.join(result)
    return pdbqt


def find_vina_executable() -> str | None:
    """Prefer an official local Vina CLI binary; do not accept an uploaded executable."""
    root = Path(__file__).resolve().parent
    scripts = Path(sysconfig.get_path('scripts'))
    candidates = [root/'bin'/'vina.exe', root/'bin'/'vina',
                  scripts/'vina.exe', scripts/'vina']
    for item in candidates:
        if item.exists() and item.is_file():
            return str(item)
    return shutil.which('vina')


def _affinities_from_pdbqt(pdbqt: str) -> list[float]:
    values=[]
    for line in pdbqt.splitlines():
        match=re.match(r'^REMARK\s+VINA RESULT:\s*(-?\d+(?:\.\d+)?)', line)
        if match:
            values.append(float(match.group(1)))
    return values


def _run_vina_cli(executable: str,lig_pdbqt: str,receptor_pdbqt: str,
                  center: list[float],box_size: list[float],
                  exhaustiveness: int,n_poses: int) -> dict[str,Any]:
    with tempfile.TemporaryDirectory(prefix='mdl_vina_cli_') as work:
        wd=Path(work)
        rec=wd/'receptor.pdbqt'; lig=wd/'ligand.pdbqt'; out=wd/'docked_poses.pdbqt'
        rec.write_text(receptor_pdbqt,encoding='utf-8')
        lig.write_text(lig_pdbqt,encoding='utf-8')
        command=[executable,'--receptor',str(rec),'--ligand',str(lig),
                 '--center_x',str(center[0]),'--center_y',str(center[1]),'--center_z',str(center[2]),
                 '--size_x',str(box_size[0]),'--size_y',str(box_size[1]),'--size_z',str(box_size[2]),
                 '--out',str(out),'--cpu','2','--seed','2026',
                 '--exhaustiveness',str(exhaustiveness),'--num_modes',str(n_poses)]
        try:
            process=subprocess.run(command,cwd=wd,timeout=300,stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,text=True,encoding='utf-8',errors='replace')
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError('Vina не завершил расчёт за 5 минут. Уменьши область поиска или сложность расчёта.') from exc
        if process.returncode != 0:
            detail=(process.stderr or process.stdout or 'Нет сообщения')[-2000:]
            raise RuntimeError('AutoDock Vina завершил работу с ошибкой: '+detail)
        if not out.exists():
            raise RuntimeError('Vina не создал файл рассчитанных поз. '+(process.stdout or '')[-1500:])
        poses=out.read_text(encoding='utf-8',errors='replace')
    affinities=_affinities_from_pdbqt(poses)
    if not affinities:
        # Some formats print the score only in standard output. No invented score.
        raise RuntimeError('В выводе Vina отсутствуют строки с энергией «REMARK VINA RESULT». Расчёт нельзя представить как проверенный.')
    return {'affinities':[round(x,3) for x in affinities],
            'poses_pdbqt':poses,'first_pose_pdbqt':read_pose_first(poses),
            'center':center,'box_size':box_size,'exhaustiveness':exhaustiveness,
            'n_poses':len(affinities),'engine':'Vina CLI'}


def run_focused_docking(sdf_text: str, receptor_pdbqt: str,
                        center: list[float], box_size: list[float],
                        exhaustiveness: int = 8, n_poses: int = 5) -> dict[str,Any]:
    """Perform rigid Vina docking via local Windows CLI or optional Python API."""
    validate_box(center,box_size)
    if not 1 <= exhaustiveness <= 32 or not 1 <= n_poses <= 10:
        raise ValueError('Неправильное число проходов поиска или поз.')
    if not any(line.startswith(('ATOM ','HETATM')) for line in receptor_pdbqt.splitlines()):
        raise ValueError('Файл рецептора PDBQT не содержит координат атомов.')
    lig_pdbqt=ligand_pdbqt(sdf_text)
    executable=find_vina_executable()
    if executable:
        return _run_vina_cli(executable,lig_pdbqt,receptor_pdbqt,center,box_size,exhaustiveness,n_poses)
    # Non-Windows users may install Vina Python bindings directly.
    try:
        from vina import Vina
    except ImportError as exc:
        raise RuntimeError('Не найден AutoDock Vina. Для Windows скачай официальный vina.exe и положи в папку bin проекта.') from exc
    with tempfile.TemporaryDirectory(prefix='mdl_vina_') as work:
        workdir=Path(work)
        receptor=workdir/'receptor.pdbqt'
        poses_path=workdir/'docked_poses.pdbqt'
        receptor.write_text(receptor_pdbqt,encoding='utf-8')
        vina=Vina(sf_name='vina',cpu=2,seed=2026,verbosity=0)
        vina.set_receptor(str(receptor))
        vina.set_ligand_from_string(lig_pdbqt)
        vina.compute_vina_maps(center=list(map(float,center)),box_size=list(map(float,box_size)))
        vina.dock(exhaustiveness=exhaustiveness,n_poses=n_poses)
        scores=vina.energies(n_poses=n_poses)
        vina.write_poses(str(poses_path),n_poses=n_poses,overwrite=True)
        poses = poses_path.read_text(encoding='utf-8',errors='replace')
    affinities=[round(float(row[0]),3) for row in scores]
    if not affinities or not all(math.isfinite(x) for x in affinities):
        raise RuntimeError('Vina завершил расчёт, но не вернул пригодные оценки.')
    if not any(line.startswith(('ATOM ','HETATM')) for line in poses.splitlines()):
        raise RuntimeError('Vina не сохранил координаты рассчитанных положений.')
    return {'affinities':affinities,'poses_pdbqt':poses,
            'first_pose_pdbqt':read_pose_first(poses),
            'center':center,'box_size':box_size,
            'exhaustiveness':exhaustiveness,'n_poses':len(affinities),'engine':'Vina Python API'}
