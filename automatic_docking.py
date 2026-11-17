"""One-click exploratory docking orchestration.

The code attempts non-destructive preparation first. Only as a last resort,
PDBFixer may construct missing atoms; any such coordinates are explicitly
flagged as *modelled*, never experimental. No residue is silently removed.
"""
from __future__ import annotations
from pathlib import Path
import importlib.util
import tempfile

from docking_runner import (ReceptorPreparationError, repair_receptor_pdb,
                            prepare_receptor_pdbqt, run_focused_docking)
from receptor_workflow import attempt_original_receptor, AutoPreparationError
from receptor_geometry import choose_altlocs


class AutomaticDockingError(RuntimeError):
    def __init__(self, summary: str, attempts: list[dict]):
        super().__init__(summary)
        self.attempts = attempts


def automatic_receptor(protein_pdb: str, temp: Path, *, original_prep=None,
                       repair=None, next_prep=None, prody_available=None,
                       fixer_available=None) -> tuple[str, dict]:
    """Prepare rigid receptor with provenance; never claim validation."""
    original_prep = original_prep or attempt_original_receptor
    repair = repair or repair_receptor_pdb
    next_prep = next_prep or prepare_receptor_pdbqt
    if prody_available is None: prody_available = importlib.util.find_spec('prody') is not None
    if fixer_available is None: fixer_available = importlib.util.find_spec('pdbfixer') is not None
    # Resolve crystallographic alternatives once, consistently per residue,
    # using the reported occupancy; keep the original input unchanged outside
    # this in-memory preparation attempt.
    selected_pdb, altloc_report = choose_altlocs(protein_pdb)
    failures: list[dict] = []
    try:
        output, attempts = original_prep(selected_pdb, temp/'original', prody_available=prody_available)
        return output, {'route':'original','attempts':attempts,'repaired_pdb':None,
                        'repair_report':None,'reconstruction':False,
                        'altloc_report':altloc_report}
    except AutoPreparationError as exc:
        failures.extend(exc.attempts)
    except Exception as exc:
        failures.append({'method':'Оригинальная PDB','status':'Ошибка','message':str(exc),
                         'details':str(getattr(exc,'details','')),'residues':getattr(exc,'residues',[])})
    if not fixer_available:
        raise AutomaticDockingError('Не удалось подготовить исходный белок автоматически. Для возможного восстановления недостающих атомов нужен PDBFixer; другой вариант — проверенный PDBQT.', failures)
    try:
        repaired, report = repair(selected_pdb, temp/'repair')
        if not repaired or not report:
            raise RuntimeError('Восстановление не вернуло готовую модель.')
        # Retry standard parser, then alternate only if necessary. Output atom
        # changes are retained and prominently disclosed to the user.
        try:
            pdbqt = next_prep(repaired,temp/'repaired_rdkit',parser='rdkit')
            parser='Meeko'
        except Exception as exc:
            failures.append({'method':'После PDBFixer → Meeko','status':'Ошибка',
                             'message':str(exc),'details':str(getattr(exc,'details','')),
                             'residues':getattr(exc,'residues',[])})
            if not prody_available: raise
            pdbqt=next_prep(repaired,temp/'repaired_prody',parser='prody')
            parser='ProDy → Meeko'
        if not any(s.startswith(('ATOM  ','HETATM')) for s in pdbqt.splitlines()):
            raise RuntimeError('Подготовленный PDBQT не содержит атомов.')
        return pdbqt, {'route':'PDBFixer → '+parser,'attempts':failures,
                       'repaired_pdb':repaired,'repair_report':report,
                       'reconstruction': bool(report.get('added_heavy_atoms') or report.get('missing_terminal_atoms')),
                       'altloc_report':altloc_report}
    except Exception as exc:
        failures.append({'method':'Восстановление PDBFixer','status':'Ошибка','message':str(exc),
                         'details':str(getattr(exc,'details','')),
                         'residues':getattr(exc,'residues',[])})
        raise AutomaticDockingError('Эту структуру белка не удалось безопасно подготовить автоматически. Попробуй другую экспериментальную структуру той же мишени или проверенный PDBQT.',failures) from exc


def run_automatic_docking(smiles: str, protein_pdb: str,
                          center: list[float], size: list[float],
                          exhaustiveness: int = 8, poses: int = 5,
                          receptor_pdbqt: str | None = None,
                          sdf_text: str | None = None) -> dict:
    """Full preparation and Vina run upon one user click, no prior SDF action."""
    from chem_design import conformer_sdf
    from docking_runner import validate_box
    validate_box(center,size)
    with tempfile.TemporaryDirectory(prefix='mdl_oneclick_') as directory:
        if receptor_pdbqt:
            rec = receptor_pdbqt
            report={'route':'загруженный PDBQT','attempts':[], 'reconstruction':False,
                    'repaired_pdb':None,'repair_report':None}
        else:
            rec,report=automatic_receptor(protein_pdb,Path(directory))
        sdf = sdf_text or conformer_sdf(smiles)
        result=run_focused_docking(sdf,rec,center,size,exhaustiveness,poses)
        result.update(receptor_pdbqt=rec,sdf_text=sdf,receptor_route=report)
        return result
