"""Exploratory structure-based ligand pharmacophore features in the protein reference frame.

Source: RDKit BaseFeatures.fdef on experimentally positioned, chemically repaired ligand
via ProLIF. Features and their spatial coincidences are hypotheses, not verified activity.
"""
from __future__ import annotations

from collections import defaultdict
import os

from Bio.PDB import Superimposer
from analysis import load_model, sequence_alignment

FEATURE_NAMES = {
    'Donor': 'Донор H-связи',
    'Acceptor': 'Акцептор H-связи',
    'Aromatic': 'Ароматический центр',
    'Hydrophobe': 'Гидрофобный центр',
    'LumpedHydrophobe': 'Гидрофобная область',
    'PosIonizable': 'Положительно ионизируемый',
    'NegIonizable': 'Отрицательно ионизируемый',
}
FEATURE_COLORS = {
    'Donor':'#6A9EC9','Acceptor':'#CF7F83','Aromatic':'#B59B5F',
    'Hydrophobe':'#7C9E80','LumpedHydrophobe':'#7C9E80',
    'PosIonizable':'#967AB4','NegIonizable':'#D28C5E',
}


def _extract_ligand(pdb_text, ligand, smiles):
    """Identical molecular reconstruction route to the working ProLIF module."""
    from rdkit import Chem
    import prolif as plf
    from prolif.io import MoleculeStandardizer
    if not smiles:
        raise ValueError('Не найден химический шаблон (SMILES) в базе PDB.')
    template = Chem.MolFromSmiles(smiles)
    mol = Chem.MolFromPDBBlock(pdb_text, removeHs=False, sanitize=False)
    if template is None or mol is None:
        raise ValueError('Не удалось прочитать молекулу/комплекс в RDKit.')
    system = MoleculeStandardizer(templates=[(ligand['Код компонента'], template)])(mol)
    def matches(residue):
        return (residue.name == ligand['Код компонента'] and
                residue.number == ligand['Номер'] and residue.chain == ligand['Цепь'])
    ligand_mol, _ = plf.split_molecule(system, matches)
    if not ligand_mol.n_residues or not ligand_mol.GetNumConformers():
        raise ValueError('Химически подготовленный лиганд не имеет 3D-координат.')
    # ProLIF can successfully use a split ligand without initializing RDKit's
    # RingInfo cache. Its chemical-feature factory cannot: sanitize a COPY
    # after splitting, preserving the crystallographic conformer coordinates.
    # Do not silently invent ring/aromatic properties if sanitization fails.
    return prepare_feature_molecule(ligand_mol)


def prepare_feature_molecule(mol):
    """Return a sanitized RDKit copy with 3D geometry retained for featurization.

    Submolecule extraction may invalidate/omit the ring cache even when the
    parent complex had been standardized using a CCD chemical template.
    """
    from rdkit import Chem
    repaired = Chem.Mol(mol)
    if repaired.GetNumConformers() == 0:
        raise ValueError('У лиганда нет 3D-координат, необходимые признаки не определить.')
    try:
        Chem.SanitizeMol(repaired)
        # Explicitly require initialized ring data before SMARTS-based features.
        repaired.GetRingInfo().NumRings()
    except Exception as exc:
        raise ValueError(
            'Не удалось подготовить кольца, ароматичность и валентности '
            'лиганда по химической модели PDB; проверьте CCD-шаблон и связи: '
            f'{exc}'
        ) from exc
    return repaired


def ligand_feature_points(ligand_mol, label):
    """RDKit molecular-feature positions in the ligand's crystallographic frame."""
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    factory = ChemicalFeatures.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
    # Also protect direct callers; extraction is not the only possible source.
    clean_mol = prepare_feature_molecule(ligand_mol)
    points=[]
    for feat in factory.GetFeaturesForMol(clean_mol):
        family=feat.GetFamily()
        if family not in FEATURE_NAMES:
            continue
        xyz=feat.GetPos()
        points.append({'family':family,'name':FEATURE_NAMES[family],
                       'color':FEATURE_COLORS[family], 'xyz':[float(xyz.x),float(xyz.y),float(xyz.z)],
                       'ligand':label})
    return points


def fit_to_reference(reference, other, cifs, chains):
    """Other -> reference 3D coordinate transform; reject unrelated protein targets."""
    import numpy as np
    if reference['_uid'] == other['_uid']:
        return np.eye(3),np.zeros(3)
    ref_chain=reference.get('_primary_chain','')
    other_chain=other.get('_primary_chain','')
    ref_seq=chains.get(reference['PDB'],{}).get(ref_chain,[])
    other_seq=chains.get(other['PDB'],{}).get(other_chain,[])
    align=sequence_alignment(ref_seq,other_seq)
    if not align['ok'] or align['identity']<0.9 or align['coverage']<0.8:
        raise ValueError('Белковые цепи не относятся к достаточно похожей мишени.')
    fixed_model=load_model(reference['PDB'],cifs[reference['PDB']])
    moving_model=load_model(other['PDB'],cifs[other['PDB']])
    fixed=[];moving=[]
    for i,j in align['pairs']:
        r1=ref_seq[i];r2=other_seq[j]
        rid=(' ',r1['number'],r1['insertion'] or ' ')
        oid=(' ',r2['number'],r2['insertion'] or ' ')
        if rid in fixed_model[ref_chain] and oid in moving_model[other_chain]:
            a=fixed_model[ref_chain][rid];b=moving_model[other_chain][oid]
            if 'CA' in a and 'CA' in b:
                fixed.append(a['CA']);moving.append(b['CA'])
    if len(fixed)<3:
        raise ValueError('Недостаточно совпадающих атомов Cα для совмещения.')
    s=Superimposer();s.set_atoms(fixed,moving)
    return s.rotran


def cluster_features(points, tolerance=2.0, min_ligands=2):
    """Greedy same-feature centroid clusters with distinct source ligands.

    Deterministic exploratory clustering; avoids counting multiple ligand atoms from
    a single ligand as replication across compounds.
    """
    import numpy as np
    groups=defaultdict(list)
    for pt in points:groups[pt['family']].append(pt)
    clusters=[]
    for family in sorted(groups):
        members=sorted(groups[family],key=lambda v:(v['ligand'], tuple(round(x,3) for x in v['xyz'])))
        grouped=[]
        for point in members:
            xyz=np.asarray(point['xyz'],dtype=float)
            options=[(np.linalg.norm(np.mean([np.array(x['xyz']) for x in cluster],axis=0)-xyz),i)
                     for i,cluster in enumerate(grouped) if point['ligand'] not in {x['ligand'] for x in cluster}]
            if options and min(options)[0]<=tolerance:
                grouped[min(options)[1]].append(point)
            else:
                grouped.append([point])
        for members in grouped:
            ligands=sorted({x['ligand'] for x in members})
            if len(ligands)<min_ligands:
                continue
            center=np.mean([x['xyz'] for x in members],axis=0)
            clusters.append({'Признак':FEATURE_NAMES[family], 'Семейство':family,
               'Лигандов':len(ligands),'Комплексы':', '.join(ligands),
               'X (Å)':round(float(center[0]),2),'Y (Å)':round(float(center[1]),2),
               'Z (Å)':round(float(center[2]),2),'_xyz':[float(v) for v in center],
               '_color':FEATURE_COLORS[family]})
    return sorted(clusters,key=lambda c:(-c['Лигандов'],c['Признак'],c['X (Å)']))


def compute_pharmacophore(selected, cifs, chains, complex_func, info_func, tolerance=2.0):
    """Returns per-ligand feature points and recurring consensuses, with errors."""
    if not selected:return [],[],[]
    ref=selected[0]
    points=[];errors=[];valid=[]
    for lig in selected:
        label=f"{lig['PDB']} · {lig['Код компонента']}"
        try:
            rot,trans=fit_to_reference(ref,lig,cifs,chains)
            mol=_extract_ligand(complex_func(lig),lig,info_func(lig).get('smiles',''))
            feature_points=ligand_feature_points(mol,label)
            for p in feature_points:
                import numpy as np
                p['xyz']=(np.array(p['xyz'])@rot+trans).tolist()
            points.extend(feature_points)
            valid.append(label)
        except Exception as ex:
            errors.append(f'{label}: {ex}')
    # When the user requests consensus across multiple complexes, do not report
    # features from a single successfully processed ligand as common features.
    min_num=2 if len(selected)>=2 else 1
    clusters=cluster_features(points,tolerance=tolerance,min_ligands=min_num)
    return points,clusters,errors
