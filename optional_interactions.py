"""Экспериментальная интеграция ProLIF; не заменяет контроль препаративных данных.

Требует отдельной установки `pip install prolif` и химического шаблона
лиганда из RCSB CCD. Если преобразование невозможно, возвращает причину,
а не выдумывает наличие водородных или других связей.
"""
from __future__ import annotations


def run_prolif(pdb_text:str, ligand:dict, smiles:str)->list[dict]:
    if not smiles:
        raise RuntimeError('Для ProLIF нужен SMILES лиганда из CCD; он не получен.')
    try:
        from rdkit import Chem
        import prolif as plf
        from prolif.io import MoleculeStandardizer
    except ImportError as exc:
        raise RuntimeError('ProLIF не установлен. Отдельная команда: py -m pip install prolif') from exc
    template=Chem.MolFromSmiles(smiles)
    if template is None:raise RuntimeError('SMILES лиганда не удалось прочитать RDKit.')
    mol=Chem.MolFromPDBBlock(pdb_text,removeHs=False,sanitize=False)
    if mol is None:raise RuntimeError('Невозможно преобразовать комплекс из PDB в RDKit.')
    try:
        standardizer=MoleculeStandardizer(templates=[(ligand['Код компонента'],template)])
        system=standardizer(mol)
        def matches(r):
            return (r.name==ligand['Код компонента'] and r.number==ligand['Номер']
                    and r.chain==ligand['Цепь'])
        ligand_mol, protein_mol=plf.split_molecule(system,matches)
        if not ligand_mol.n_residues or not protein_mol.n_residues:
            raise RuntimeError('Не удалось однозначно разделить лиганд и белок.')
        fp=plf.Fingerprint(implicit_hydrogens=True)
        fp.run_from_iterable([ligand_mol],protein_mol, n_jobs=1)
    except Exception as exc:
        raise RuntimeError(f'ProLIF не смог подготовить этот комплекс: {exc}') from exc
    names={
        'Hydrophobic':'Гидрофобный контакт','HBDonor':'Водородная связь (донор лиганда)',
        'HBAcceptor':'Водородная связь (акцептор лиганда)',
        'ImplicitHBDonor':'Водородная связь (донор лиганда)',
        'ImplicitHBAcceptor':'Водородная связь (акцептор лиганда)',
        'PiStacking':'Ароматическое π-взаимодействие',
        'Anionic':'Ионное (анион)','Cationic':'Ионное (катион)',
        'CationPi':'Катион–π', 'PiCation':'π–катион', 'VdWContact':'Ван-дер-Ваальсов контакт',
    }
    out=[]
    for (_l,protein_res), kinds in fp.ifp.get(0,{}).items():
        for kind, occurrences in kinds.items():
            dist = None
            if occurrences and isinstance(occurrences[0],dict):
                dist=occurrences[0].get('distance')
            out.append({'Остаток':str(protein_res),'Взаимодействие':names.get(kind,kind),
                        'Расстояние (Å)':round(float(dist),2) if dist is not None else None,
                        'Метод':'ProLIF (расчётное)'} )
    return out
