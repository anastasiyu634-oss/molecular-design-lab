"""Convert Vina/Meeko PDBQT poses to bond-complete SDF conformers.

PDBQT contains a torsion tree, not a complete bond graph. Meeko writes a
SMILES and an atom-index mapping in REMARK lines. Use those to preserve the
intended chemistry instead of guessing bonds from distances.
"""
from __future__ import annotations

import io
import math
import re
from rdkit import Chem
from rdkit.Geometry import Point3D


class PoseFormatError(ValueError):
    pass


def _blocks(text: str) -> list[str]:
    clean=text.replace('\x00','').replace('\r','')
    model_blocks=re.findall(r'(?ms)^MODEL\s+\d+[^\n]*\n(.*?)^ENDMDL[^\n]*$',clean)
    if model_blocks:
        return model_blocks
    return [clean] if any(line.startswith(('ATOM','HETATM')) for line in clean.splitlines()) else []


def _pose_to_mol(block: str, index: int):
    smi_line=next((s for s in block.splitlines() if s.startswith('REMARK SMILES ')
                   and not s.startswith('REMARK SMILES IDX ')),None)
    if not smi_line:
        raise PoseFormatError('В PDBQT нет REMARK SMILES: нельзя восстановить правильные связи автоматически.')
    smiles=smi_line[len('REMARK SMILES '):].strip()
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:
        raise PoseFormatError('Невозможно разобрать SMILES рассчитанной позы.')
    mapping={}
    for line in block.splitlines():
        if line.startswith('REMARK SMILES IDX '):
            try:
                numbers=[int(x) for x in line[len('REMARK SMILES IDX '):].split()]
            except ValueError as exc:
                raise PoseFormatError('Неверная карта атомов SMILES IDX.') from exc
            if len(numbers)%2:
                raise PoseFormatError('Неполная карта атомов SMILES IDX.')
            mapping.update(zip(numbers[::2],numbers[1::2]))
    if set(mapping)!=set(range(1,mol.GetNumAtoms()+1)) or len(set(mapping.values()))!=len(mapping):
        raise PoseFormatError('Карта атомов неполна или содержит повторы.')
    coord_map={}
    elements={}
    for line in block.splitlines():
        if not line.startswith(('ATOM  ','HETATM')):
            continue
        try:
            pdbqt_atom_index=int(line[6:11])
            x,y,z=(float(line[start:start+8]) for start in (30,38,46))
        except (ValueError,IndexError) as exc:
            raise PoseFormatError('Невозможно прочитать координаты атомов PDBQT.') from exc
        if not all(map(math.isfinite,(x,y,z))):
            raise PoseFormatError('Найдены нечисловые координаты PDBQT.')
        coord_map[pdbqt_atom_index]=(x,y,z)
        elements[pdbqt_atom_index]=line[12:16].strip()[0:1].upper()
    if not set(mapping.values()).issubset(coord_map):
        raise PoseFormatError('Координаты некоторых тяжёлых атомов отсутствуют.')
    conf=Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for smiles_atom_index,pdbqt_atom_index in mapping.items():
        atom=mol.GetAtomWithIdx(smiles_atom_index-1)
        if atom.GetSymbol().upper()!=elements[pdbqt_atom_index]:
            raise PoseFormatError('SMILES и PDBQT не совпадают по химическим элементам.')
        conf.SetAtomPosition(smiles_atom_index-1,Point3D(*coord_map[pdbqt_atom_index]))
    mol.RemoveAllConformers()
    mol.AddConformer(conf,assignId=True)
    # Protect against a corrupt atom map even if all indices were present.
    for bond in mol.GetBonds():
        a=conf.GetAtomPosition(bond.GetBeginAtomIdx())
        b=conf.GetAtomPosition(bond.GetEndAtomIdx())
        d=(a-b).Length()
        if not 0.85 <= d <= 2.30:
            raise PoseFormatError(f'Обнаружена некорректная длина связи {d:.2f} Å в позе {index}.')
    score_match=re.search(r'(?m)^REMARK VINA RESULT:\s*(-?\d+(?:\.\d+)?)',block)
    score=float(score_match.group(1)) if score_match else None
    mol.SetProp('_Name',f'Docking pose {index}')
    mol.SetProp('Pose',str(index))
    mol.SetProp('Source_SMILES',smiles)
    if score is not None:
        mol.SetProp('Vina_score_kcal_mol',str(score))
    return mol, score


def pdbqt_poses_to_mols(text: str):
    blocks=_blocks(text)
    if not blocks:
        raise PoseFormatError('В файле нет рассчитанных поз.')
    return [_pose_to_mol(block,i+1) for i,block in enumerate(blocks)]


def pdbqt_poses_to_sdf(text: str) -> list[dict]:
    """One full SDF molblock per pose, including bonds and score metadata."""
    results=[]
    for i,(mol,score) in enumerate(pdbqt_poses_to_mols(text),1):
        out=io.StringIO()
        w=Chem.SDWriter(out)
        w.write(mol)
        w.flush()
        w.close()
        results.append({'index':i,'score':score,'sdf':out.getvalue(),'heavy_atoms':mol.GetNumHeavyAtoms()})
    return results
