"""Conservative receptor preparation workflow; no deletion or silent rebuilding.

Tries ordinary Meeko and then ProDy reading of the IDENTICAL original PDB. If both
fail, caller may explicitly choose PDBFixer (unvalidated) or upload reviewed PDBQT.
"""
from __future__ import annotations
from pathlib import Path
from docking_runner import prepare_receptor_pdbqt, ReceptorPreparationError

class AutoPreparationError(ReceptorPreparationError):
    def __init__(self, summary, attempts):
        self.attempts=attempts
        details='\n\n'.join(f'{a["method"]}: {a["message"]}\n{a["details"]}' for a in attempts)
        residues=list(dict.fromkeys(r for a in attempts for r in a['residues']))
        super().__init__(summary,details,residues)


def attempt_original_receptor(protein_pdb:str, root:Path, prepare=None, prody_available=True):
    """Return (pdbqt, attempts), or raise with exact diagnostics.

    Does not call any atom repair. Does not perform docking.
    """
    prepare=prepare or prepare_receptor_pdbqt
    trials=[('Meeko (обычное чтение)','rdkit')]
    if prody_available:
        trials.append(('Meeko (чтение ProDy)','prody'))
    failures=[]
    for label,parser in trials:
        try:
            result=prepare(protein_pdb, root/parser, parser=parser)
            if not result or not any(line.startswith(('ATOM  ','HETATM')) for line in result.splitlines()):
                raise ReceptorPreparationError('Выходной PDBQT не содержит координат атомов.')
            return result,failures+[{'method':label,'status':'OK','message':'Подготовлено из исходного PDB', 'details':'', 'residues':[]}]
        except Exception as exc:
            failures.append({'method':label,'status':'Ошибка','message':str(exc),
                             'details':str(getattr(exc,'details','')),
                             'residues':list(getattr(exc,'residues',[]))})
    raise AutoPreparationError('Оба неразрушающих способа подготовки не прошли. Для продолжения проверьте структуру или загрузите проверенный PDBQT.',failures)
