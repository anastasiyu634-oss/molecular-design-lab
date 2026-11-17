"""Reference-ligand redocking and *in-place* heavy-atom RMSD validation.

Do not superimpose predicted and crystallographic ligands: they already share
protein coordinates. Symmetry-related atom assignments may be exchanged.
"""
from __future__ import annotations

import math


class RedockingValidationError(ValueError):
    pass


def experimental_ligand(complex_pdb: str, ligand: dict, smiles: str):
    """Use the same CCD-backed chemical reconstruction as working ProLIF."""
    from rdkit import Chem
    from pharmacophore import _extract_ligand
    mol = Chem.RemoveHs(Chem.Mol(_extract_ligand(complex_pdb, ligand, smiles)))
    if not mol.GetNumConformers() or mol.GetNumHeavyAtoms() < 3:
        raise RedockingValidationError('В экспериментальном комплексе нет полноценной 3D-структуры лиганда.')
    return mol


def crystallographic_sdf(reference_mol) -> str:
    """Add hydrogens without changing observed heavy-atom coordinates."""
    from rdkit import Chem
    mol = Chem.AddHs(Chem.Mol(reference_mol), addCoords=True)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            raise RedockingValidationError('Неизвестный химический элемент экспериментального лиганда.')
    mol.SetProp('_Name', 'Experimental ligand for redocking')
    return Chem.MolToMolBlock(mol) + '\n$$$$\n'


def rmsd_against_reference(reference_mol, predicted_mols: list) -> list[float]:
    """Symmetry-aware, NO molecular alignment, heavy atom RMSD in receptor frame."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    reference=Chem.RemoveHs(Chem.Mol(reference_mol))
    if reference.GetNumConformers()!=1:
        raise RedockingValidationError('Экспериментальный лиганд должен иметь один 3D-конформер.')
    values=[]
    for i,pose in enumerate(predicted_mols,1):
        probe=Chem.RemoveHs(Chem.Mol(pose))
        if probe.GetNumAtoms()!=reference.GetNumAtoms():
            raise RedockingValidationError(f'Поза {i}: различается число тяжёлых атомов с экспериментальным лигандом.')
        if not probe.HasSubstructMatch(reference, useChirality=False) or not reference.HasSubstructMatch(probe, useChirality=False):
            raise RedockingValidationError(f'Поза {i}: не совпадает граф химических связей с эталоном.')
        try:
            value=rdMolAlign.CalcRMS(probe, reference, maxMatches=10000, symmetrizeConjugatedTerminalGroups=True)
        except Exception as exc:
            raise RedockingValidationError(f'Поза {i}: невозможно рассчитать RMSD без дополнительного выравнивания: {exc}') from exc
        if not math.isfinite(value):
            raise RedockingValidationError(f'Поза {i}: RMSD не является конечным числом.')
        values.append(round(float(value),3))
    if not values:
        raise RedockingValidationError('В результатах нет поз для сравнения.')
    return values


def run_control_redocking(reference_mol, protein_pdb, center, size,
                          *, exhaustiveness=8, poses=5, receptor_pdbqt=None):
    """Run Vina on the *co-crystal ligand*, not the user's test molecule."""
    from automatic_docking import run_automatic_docking
    from docking_pose_format import pdbqt_poses_to_mols
    from rdkit import Chem
    source_smiles=Chem.MolToSmiles(Chem.RemoveHs(reference_mol), isomericSmiles=True)
    output=run_automatic_docking(source_smiles, protein_pdb, center, size,
                                 exhaustiveness=exhaustiveness,poses=poses,
                                 receptor_pdbqt=receptor_pdbqt,
                                 sdf_text=crystallographic_sdf(reference_mol))
    predicted=[mol for mol,_ in pdbqt_poses_to_mols(output['poses_pdbqt'])]
    rmsds=rmsd_against_reference(reference_mol,predicted)
    scores=output.get('affinities',[])
    if len(scores)!=len(rmsds):
        raise RedockingValidationError('Количество RMSD и оценок Vina не совпадает.')
    output['reference_rmsds']=rmsds
    output['best_rmsd']=min(rmsds)
    output['first_pose_rmsd']=rmsds[0]
    return output
