"""Exploratory BRICS recombination and candidate structure preparation.

This is fragment recombination, not validated de novo optimization or pharmacophore matching.
"""
from __future__ import annotations

from itertools import islice


def properties(smiles):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Неверный SMILES: RDKit не смог прочитать структуру.')
    return {'SMILES':Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True),
            'Масса (Да)':round(Descriptors.MolWt(mol),1),
            'cLogP':round(Crippen.MolLogP(mol),2),
            'HBA':rdMolDescriptors.CalcNumHBA(mol),
            'HBD':rdMolDescriptors.CalcNumHBD(mol),
            'Тяжёлых атомов':mol.GetNumHeavyAtoms()}


def generate_candidates(smiles_list, limit=8):
    from rdkit import Chem
    from rdkit.Chem import BRICS
    parents=[]
    for sm in smiles_list:
        mol=Chem.MolFromSmiles(sm or '')
        if mol is not None:
            parents.append(mol)
    if not parents:
        raise ValueError('Для генерации нужны хотя бы одна корректная исходная молекула.')
    originals={Chem.MolToSmiles(m) for m in parents}
    frags=set()
    for parent in parents:
        frags.update(BRICS.BRICSDecompose(parent))
    pool=[Chem.MolFromSmiles(s) for s in sorted(frags)[:60]]
    pool=[x for x in pool if x is not None]
    if len(pool)<2:
        return [],len(pool)
    found=[];seen=set()
    products=BRICS.BRICSBuild(pool, scrambleReagents=False, maxDepth=2)
    # Bound evaluation to avoid spending unbounded CPU time on large fragment pools.
    for mol in islice(products,250):
        try:
            Chem.SanitizeMol(mol)
            sm=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
            if sm in seen or sm in originals:continue
            data=properties(sm)
            if not (10<=data['Тяжёлых атомов']<=45 and data['Масса (Да)']<=550 and
                    data['cLogP']<=6 and data['HBA']<=10 and data['HBD']<=5):
                continue
            found.append(data);seen.add(sm)
            if len(found)>=limit:break
        except (ValueError,RuntimeError):
            continue
    return found,len(pool)


def conformer_sdf(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('Не удалось прочитать SMILES.')
    mol=Chem.AddHs(mol)
    status=AllChem.EmbedMolecule(mol,randomSeed=42,maxAttempts=80)
    if status!=0:raise ValueError('RDKit не смог построить 3D-конформер для этой молекулы.')
    try:
        props=AllChem.MMFFGetMoleculeProperties(mol)
        if props is not None:AllChem.MMFFOptimizeMolecule(mol,maxIters=250)
        else:AllChem.UFFOptimizeMolecule(mol,maxIters=250)
    except (ValueError,RuntimeError):
        pass
    mol.SetProp('_Name','Molecular Design Lab — candidate (not docked)')
    return Chem.MolToMolBlock(mol)+'\n$$$$\n'


def protein_only_pdb(pdb_id,cif_text,chain_id):
    """Extract protein coordinates, *not* docking-ready protonated PDBQT."""
    from io import StringIO
    from Bio.PDB import PDBIO,Select,is_aa
    from analysis import load_model
    class ProteinChain(Select):
        def accept_chain(self,chain):
            return int(chain.id==chain_id)
        def accept_residue(self,res):
            return int(res.id[0]==' ' and is_aa(res,standard=False))
    buf=StringIO();io=PDBIO();io.set_structure(load_model(pdb_id,cif_text))
    io.save(buf,select=ProteinChain())
    return buf.getvalue()


def generate_analogs(smiles_list, limit=8):
    """Small, single-site structural analogues with original scaffold retained.

    Experimental chemical enumeration, *not* activity prediction, not a
    validated pharmacophore filter and not a synthesis recommendation.
    """
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions, rdFMCS
    rules=(
        ('Ароматическое фторирование', '[cH:1]>>[c:1]F'),
        ('Ароматическое метилирование', '[cH:1]>>[c:1]C'),
        ('Ароматическое хлорирование', '[cH:1]>>[c:1]Cl'),
        ('Метилирование O–H', '[O;H1:1]>>[O:1]C'),
        ('Замена F на Cl', '[F:1]>>[Cl:1]'),
        ('Замена Cl на F', '[Cl:1]>>[F:1]'),
    )
    parents=[]
    for smi in smiles_list:
        m=Chem.MolFromSmiles(smi or '')
        if m is not None:parents.append(m)
    if not parents:raise ValueError('Нет корректных исходных SMILES.')
    originals={Chem.MolToSmiles(m,isomericSmiles=True) for m in parents}
    out=[];seen=set();generated=0
    for parent in parents:
        for label,smarts in rules:
            reaction=rdChemReactions.ReactionFromSmarts(smarts)
            for product_tuple in reaction.RunReactants((parent,))[:100]:
                generated+=1
                if generated>300:return out,generated-1
                product=product_tuple[0]
                try:
                    Chem.SanitizeMol(product)
                    sm=Chem.MolToSmiles(product,canonical=True,isomericSmiles=True)
                    if sm in seen or sm in originals:continue
                    data=properties(sm)
                    if not (8<=data['Тяжёлых атомов']<=45 and data['Масса (Да)']<=550 and
                            -2<=data['cLogP']<=6 and data['HBA']<=10 and data['HBD']<=5):continue
                    # A single local substitution must preserve almost all
                    # heavy atoms of the original chemical graph.
                    mcs=rdFMCS.FindMCS([parent,product],timeout=1,ringMatchesRingOnly=True,
                                      completeRingsOnly=True)
                    if mcs.canceled or mcs.numAtoms/max(parent.GetNumHeavyAtoms(),1)<0.8:continue
                    data.update({'Способ':label,'Сохранение каркаса (%)':round(100*mcs.numAtoms/parent.GetNumHeavyAtoms()),
                                 'Фармакофор':'Не проверен','Докинг':'Не проведён'})
                    out.append(data);seen.add(sm)
                    if len(out)>=limit:return out,generated
                except (ValueError,RuntimeError):continue
    return out,generated
