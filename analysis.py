"""Molecular Design Lab: mmCIF, выбор небелковых молекул и геометрия.

Контакты по расстоянию НЕ являются доказательством H-связей, ингибирования
или фармакофорной модели. Дополнительная ProLIF-проверка находится в
optional_interactions.py и требует корректных химических шаблонов.
"""
from __future__ import annotations

from collections import Counter
from io import StringIO
import re

import requests
from Bio.PDB import MMCIFParser, NeighborSearch, Superimposer, is_aa
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner

SERVICE_CODES = {'HOH','DOD','SO4','PO4','GOL','EDO','DMS','ACT','PEG','FMT','CL','NA','K','MG','CA','ZN','IOD'}
KNOWN_ROLES = {'HEM': 'Кофактор (гем)', 'PO4': 'Фосфат, не ингибитор', 'SB2':'Ингибитор SB203580 (для 1A9U)'}


def parse_pdb_ids(raw: str, max_count: int = 10) -> list[str]:
    ids = list(dict.fromkeys(p.strip().upper() for p in re.split(r'[,;\s]+', raw) if p.strip()))
    if not ids:
        raise ValueError('Введите хотя бы один код PDB.')
    bad = [p for p in ids if not re.fullmatch(r'[A-Z0-9]{4}', p)]
    if bad:
        raise ValueError('Некорректные PDB ID: ' + ', '.join(bad))
    if len(ids) > max_count:
        raise ValueError(f'Можно загрузить не более {max_count} структур.')
    return ids


def download_structure(pdb_id: str) -> str:
    r = requests.get(f'https://files.rcsb.org/download/{pdb_id}.cif', timeout=40)
    r.raise_for_status()
    return r.text


def download_chemcomp(comp_id: str) -> dict:
    """Данные CCD: название, формула и SMILES, если они опубликованы."""
    comp_id = comp_id.upper()
    if not re.fullmatch(r'[A-Z0-9]{1,5}', comp_id):
        return {}
    try:
        resp = requests.get(f'https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}', timeout=12)
        resp.raise_for_status()
        data = resp.json()
        chem = data.get('chem_comp') or {}
        desc = data.get('rcsb_chem_comp_descriptor') or {}
        smiles = (desc.get('SMILES_CANONICAL') or desc.get('SMILES') or desc.get('smiles')
                  or desc.get('smiles_stereo') or '')
        if not smiles:
            for d in data.get('pdbx_chem_comp_descriptor', []):
                if 'SMILES' in d.get('type', '') and d.get('descriptor'):
                    smiles = d['descriptor']
                    break
        return {'name': chem.get('name') or comp_id, 'formula': chem.get('formula') or '', 'smiles': smiles}
    except (requests.RequestException, ValueError, KeyError):
        return {}


def load_model(pdb_id: str, cif_text: str):
    return next(MMCIFParser(QUIET=True).get_structure(pdb_id, StringIO(cif_text)).get_models())


def protein_chains(model) -> dict[str, list[dict]]:
    result = {}
    for ch in model:
        records = []
        for r in ch:
            if r.id[0] != ' ' or not is_aa(r, standard=False):
                continue
            if 'CA' not in r:
                continue
            ins = str(r.id[2]).strip()
            records.append({'code': seq1(r.resname, custom_map={'MSE':'M'}),
                            'resname':r.resname, 'number':r.id[1],
                            'insertion':ins, 'key':f'{r.resname}:{r.id[1]}{ins}',
                            'chain':ch.id})
        if records:
            result[ch.id] = records
    return result


def classify_atom_pair(ligand_element: str, protein_element: str, distance: float) -> str:
    l, p = ligand_element.upper(), protein_element.upper()
    if l in {'N','O'} and p in {'N','O'} and 2.4 <= distance <= 3.5:
        return 'Полярное сближение N/O'
    if l == 'C' and p == 'C' and 3.3 <= distance <= 4.0:
        return 'Неполярное сближение C/C'
    return 'Прочая близость'


def get_ligand_contacts(ligand_atoms: list, neighbor_search: NeighborSearch | None,
                        radius: float = 4.0) -> list[dict]:
    if neighbor_search is None:
        return []
    hits = {}
    for a in ligand_atoms:
        for pa in neighbor_search.search(a.coord, radius, level='A'):
            r = pa.get_parent(); ch = r.get_parent()
            ins = str(r.id[2]).strip(); resid = f'{r.id[1]}{ins}'
            ref = f'{r.resname} {ch.id}:{resid}'
            d = float(a-pa)
            typ = classify_atom_pair(a.element, pa.element, d)
            row = hits.setdefault(ref, {
                'Остаток':ref, 'Цепь белка':ch.id,
                'Минимальная дистанция (Å)':float('inf'),
                'Ближайший атом лиганда':'', 'Ближайший атом белка':'',
                'Полярные N/O':False, 'Неполярные C/C':False,
                'Пар атомов ≤ 4 Å':0, '_residue_key':f'{r.resname}:{resid}',
            })
            if d < row['Минимальная дистанция (Å)']:
                row['Минимальная дистанция (Å)'] = d
                row['Ближайший атом лиганда'] = a.name
                row['Ближайший атом белка'] = pa.name
            row['Пар атомов ≤ 4 Å'] += 1
            if typ == 'Полярное сближение N/O': row['Полярные N/O'] = True
            if typ == 'Неполярное сближение C/C': row['Неполярные C/C'] = True
    for row in hits.values():
        row['Минимальная дистанция (Å)'] = round(row['Минимальная дистанция (Å)'], 2)
        row['Тип предварительного контакта'] = contact_label(row)
    return sorted(hits.values(), key=lambda x:(x['Минимальная дистанция (Å)'],x['Остаток']))


def contact_label(row:dict) -> str:
    pol, non = row['Полярные N/O'], row['Неполярные C/C']
    if pol and non: return 'Полярное N/O + неполярное C/C'
    if pol: return 'Полярное N/O'
    if non: return 'Неполярное C/C'
    return 'Геометрическая близость (тип не установлен)'


def analyze_structure(pdb_id: str, cif_text: str, radius: float=4.0) -> tuple[list[dict], dict]:
    model = load_model(pdb_id, cif_text)
    protein_res = [r for ch in model for r in ch if r.id[0]==' ' and is_aa(r, standard=False)]
    protein_atoms = [a for r in protein_res for a in r if a.element not in {'H','D'}]
    near = NeighborSearch(protein_atoms) if protein_atoms else None
    seqs = protein_chains(model)
    ligs = []
    for ch in model:
        for r in ch:
            if not r.id[0].startswith('H_') or r.resname in {'HOH','DOD'}:
                continue
            atoms = [a for a in r if a.element not in {'H','D'}]
            if not atoms:
                continue
            contacts = get_ligand_contacts(atoms,near,radius)
            counts = Counter(x['Цепь белка'] for x in contacts)
            primary = counts.most_common(1)[0][0] if counts else ''
            code=r.resname
            large = len(atoms)>=6 and len(contacts)>0 and code not in SERVICE_CODES
            ligs.append({
                'PDB':pdb_id,'Код компонента':code,'Цепь':ch.id,
                'Номер':r.id[1], 'Вставка':str(r.id[2]).strip(),
                'Тяжёлых атомов':len(atoms),
                'Категория':'Лиганд-кандидат' if large else 'Другой небелковый компонент',
                'Близких остатков (≤ 4 Å)':len(contacts),
                'Полярных N/O остатков':sum(x['Полярные N/O'] for x in contacts),
                'Неполярных C/C остатков':sum(x['Неполярные C/C'] for x in contacts),
                'Ближайшие остатки':', '.join(x['Остаток'] for x in contacts[:8]),
                '_residue_keys':[x['_residue_key'] for x in contacts],
                '_contact_rows': contacts, '_primary_chain':primary,
                '_uid':f'{pdb_id}:{code}:{ch.id}:{r.id[1]}:{str(r.id[2]).strip()}',
            })
    ligs.sort(key=lambda r:(r['Категория']=='Лиганд-кандидат',r['Близких остатков (≤ 4 Å)'],r['Тяжёлых атомов']),reverse=True)
    return ligs, {'Белковых остатков':len(protein_res),
                  'Белковых цепей':len({r.get_parent().id for r in protein_res}),
                  '_chains':seqs}


def sequence_alignment(ref:list[dict], other:list[dict]):
    """Соотнести номера остатков по последовательностям; не по совпадению чисел."""
    if not ref or not other:
        return {'ok':False,'reason':'Нет аминокислотной последовательности цепи'}
    refseq = ''.join(x['code'] for x in ref)
    othseq = ''.join(x['code'] for x in other)
    aligner=PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                            open_gap_score=-6, extend_gap_score=-0.5)
    align = aligner.align(refseq, othseq)[0]
    pairs=[(int(i),int(j)) for i,j in zip(*align.indices) if i>=0 and j>=0]
    matches=sum(refseq[i] == othseq[j] for i,j in pairs)
    identity=matches/max(len(pairs),1)
    cov=min(len(pairs)/len(ref),len(pairs)/len(other))
    ok = identity>=0.55 and cov>=0.65
    mapping={other[j]['key']:ref[i]['key'] for i,j in pairs} if ok else {}
    return {'ok':ok,'identity':identity,'coverage':cov,'mapping':mapping,
            'pairs':pairs,'reason':'' if ok else 'Цепи недостаточно похожи для корректного сравнения'}


def compare_selected_ligands(selected:list[dict], chains:dict, cif_texts:dict|None=None) -> tuple[list[dict],list[dict],list[str]]:
    """Матрица контактов после выравнивания цепей и RMSD Cα по возможности."""
    if len(selected)<2:
        return [],[],['Для сравнения нужны как минимум две выбранные молекулы.']
    ref=selected[0]
    refchain=ref['_primary_chain']
    ref_records=chains.get(ref['PDB'],{}).get(refchain,[])
    if not ref_records:
        return [],[],['Не удалось определить цепь белка для первого лиганда.']
    refs=set(x['_residue_key'] for x in ref['_contact_rows'] if x['Цепь белка']==refchain)
    mapping_by_uid={ref['_uid']:{r['key']:r['key'] for r in ref_records}}
    report=[{'Структура':ref['PDB'], 'Цепь':refchain, 'Идентичность последовательности':'100%',
            'Покрытие':'100%', 'RMSD Cα (Å)':0.0, 'Статус':'Опорная структура'}]
    problems=[]
    for l in selected[1:]:
        ochain=l['_primary_chain']
        other=chains.get(l['PDB'],{}).get(ochain,[])
        align=sequence_alignment(ref_records,other)
        if not align['ok']:
            problems.append(f"{l['PDB']}: {align['reason']}. Матрица для этой структуры отключена.")
            continue
        mapping_by_uid[l['_uid']]=align['mapping']
        score=None
        if cif_texts and ref['PDB'] in cif_texts and l['PDB'] in cif_texts:
            try:
                rm=load_model(ref['PDB'],cif_texts[ref['PDB']]); om=load_model(l['PDB'],cif_texts[l['PDB']])
                ra=list(r for r in rm[refchain] if r.id[0]==' ' and is_aa(r,standard=False) and 'CA' in r)
                oa=list(r for r in om[ochain] if r.id[0]==' ' and is_aa(r,standard=False) and 'CA' in r)
                a,b=[],[]
                for i,j in align['pairs']:
                    if i<len(ra) and j<len(oa): a.append(ra[i]['CA']);b.append(oa[j]['CA'])
                if len(a)>=3:
                    sup=Superimposer();sup.set_atoms(a,b);score=round(sup.rms,2)
            except (KeyError,ValueError,IndexError):
                pass
        report.append({'Структура':l['PDB'],'Цепь':ochain,
                       'Идентичность последовательности':f"{align['identity']:.0%}",
                       'Покрытие':f"{align['coverage']:.0%}",
                       'RMSD Cα (Å)':score if score is not None else 'н/д',
                       'Статус':'Сопоставлена с опорной структурой'})
    if len(mapping_by_uid)<2:
        return [],report,problems
    for l in selected[1:]:
        if l['_uid'] not in mapping_by_uid:continue
        mapping=mapping_by_uid[l['_uid']]
        refs.update(mapping[x['_residue_key']] for x in l['_contact_rows']
                    if x['Цепь белка']==l['_primary_chain'] and x['_residue_key'] in mapping)
    result=[]
    for k in sorted(refs, key=lambda x: (int(re.search(r'(-?\d+)',x).group(1)) if re.search(r'(-?\d+)',x) else 0,x)):
        row={'Остаток опорной цепи':k}
        for l in selected:
            if l['_uid'] not in mapping_by_uid:continue
            m=mapping_by_uid[l['_uid']]
            hits=[x for x in l['_contact_rows'] if x['Цепь белка']==l['_primary_chain'] and m.get(x['_residue_key'])==k]
            title=f"{l['PDB']} · {l['Код компонента']} ({l['Цепь']}:{l['Номер']})"
            row[title]= contact_label(hits[0]) if hits else '—'
        result.append(row)
    return result,report,problems


# Старые функции сохранены для совместимости с v3/тестами

def common_residues(ligands:list[dict])->list[str]:
    if len(ligands)<2:return []
    return sorted(set.intersection(*(set(x.get('_residue_keys',[])) for x in ligands)))


def make_contact_matrix(ligands:list[dict])->list[dict]:
    keys=sorted({r['_residue_key'] for l in ligands for r in l.get('_contact_rows',[])})
    result=[]
    for key in keys:
        row={'Остаток':key}
        for l in ligands:
            hits=[r for r in l['_contact_rows'] if r['_residue_key']==key]
            if not hits: cell='—'
            else:
                tags=[]
                if any(x['Полярные N/O'] for x in hits):tags.append('N/O')
                if any(x['Неполярные C/C'] for x in hits):tags.append('C/C')
                cell='+'.join(tags) if tags else 'рядом'
            row[f"{l['PDB']} · {l['Код компонента']}"]=cell
        result.append(row)
    return result
