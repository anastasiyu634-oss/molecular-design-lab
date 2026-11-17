"""Optional ProLIF fingerprints for actual predicted Vina poses.

Interprets chemical contacts in a *prepared/modelled* structure. Does not
validate docking nor transform a Vina score into experimental affinity.
"""
from __future__ import annotations


def pose_contact_rows(pdbqt_poses: str, receptor_pdb: str, pose_index: int = 1):
    from rdkit import Chem
    from docking_pose_format import pdbqt_poses_to_mols
    from optional_interactions import run_prolif

    poses=pdbqt_poses_to_mols(pdbqt_poses)
    if not 1<=pose_index<=len(poses):
        raise ValueError('Выбранной позы нет в рассчитанных данных.')
    ligand=Chem.Mol(poses[pose_index-1][0])
    smi=ligand.GetProp('Source_SMILES')
    # ProLIF's working pipeline accepts a combined, residue-tagged PDB. Keep
    # the prepared protein coordinates and add ONLY the selected docked ligand.
    protein_lines=[line for line in receptor_pdb.splitlines()
                   if line.startswith(('ATOM  ','HETATM')) and line[17:20].strip() not in ('HOH','WAT')]
    if not protein_lines:
        raise ValueError('У модели рецептора нет PDB-координат для ProLIF.')
    # Ligand must have an explicit, unique residue label for splitting.
    ligand_block=Chem.MolToPDBBlock(ligand, flavor=0)
    ligand_lines=[]
    serial=80000
    for row in ligand_block.splitlines():
        if not row.startswith(('ATOM  ','HETATM')):continue
        serial+=1
        if serial>99999:
            raise ValueError('Слишком много атомов для PDB-формата.')
        row='HETATM'+f'{serial:5d}'+row[11:17]+'LIG'+' Z'+f'{901:4d}'+row[26:]
        ligand_lines.append(row)
    if len(ligand_lines)!=ligand.GetNumAtoms():
        raise ValueError('Не удалось перенести все атомы позы в PDB.')
    joined='\n'.join(protein_lines+['TER']+ligand_lines+['END'])+'\n'
    return run_prolif(joined,{'Код компонента':'LIG','Цепь':'Z','Номер':901},smi)
