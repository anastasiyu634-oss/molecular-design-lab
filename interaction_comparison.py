"""Сравнение *рассчитанных* взаимодействий ProLIF после выравнивания цепей.

Никаких заранее вписанных PDB, типов связей или результатов.
Не рассчитанные лиганды не принимаются за отсутствие контактов.
"""
from __future__ import annotations

from collections import defaultdict
import re

from analysis import sequence_alignment

# Более короткие обозначения для таблицы, не меняющие определения ProLIF.
TYPES = (
    "H-связь", "Гидрофобный", "π-взаимодействие", "Ионное",
    "Катион–π", "Ван-дер-Ваальса",
)


def classify_interaction(text: str) -> str:
    """Унифицировать подпись русской/английской категории ProLIF."""
    v = text.lower()
    if 'водородн' in v or 'hbacceptor' in v or 'hbdonor' in v:
        return 'H-связь'
    if 'гидрофобн' in v or 'hydrophobic' in v:
        return 'Гидрофобный'
    if ('ароматич' in v or 'pistacking' in v or 'edgefacetoface' in v or 'facetoface' in v):
        return 'π-взаимодействие'
    if 'катион–π' in v or 'π–катион' in v or 'cationpi' in v or 'pication' in v:
        return 'Катион–π'
    if 'ионн' in v or v in {'cationic', 'anionic'}:
        return 'Ионное'
    if 'ван-дер' in v or 'vdwcontact' in v:
        return 'Ван-дер-Ваальса'
    return text


def parse_residue_id(value: str) -> tuple[str, str] | None:
    """`TYR35.A` → (`TYR:35`, `A`). Вставочные коды допустимы: TYR35B.A."""
    match = re.fullmatch(r'([A-Za-z]{3})(-?\d+)([A-Za-z]?)\.([A-Za-z0-9]+)', value.strip())
    if not match:
        return None
    name, number, ins, chain = match.groups()
    return f'{name.upper()}:{number}{ins.upper()}', chain


def ligand_summary(records: list[dict]) -> dict[str, int]:
    """Число уникальных остатков по каждому типу, а не число атомных пар."""
    groups = defaultdict(set)
    for record in records:
        kind = classify_interaction(str(record.get('Взаимодействие', '')))
        groups[kind].add(str(record.get('Остаток', '')))
    return {k:len(v) for k,v in groups.items()}


def readable_residue(ref_key: str) -> str:
    """`MET:109` -> `MET109`."""
    return ref_key.replace(':', '')


def compare_prolif(selected: list[dict], records_by_uid: dict, chain_data: dict):
    """Вернуть матрицу и общие типы *только* по успешно посчитанным комплексам.

    records_by_uid: {ligand['_uid']: list[ProLIF record]}.
    Возвращает (matrix, shared, issues, used_ligands), где:
      matrix = [{"Остаток":..., "<PDB> · <код>": "H-связь, ..."}, ...]
      shared = [{"Остаток": ..., "Общие взаимодействия": ...}, ...].
    Несопоставимые цепи не включаются.
    """
    calculated = [x for x in selected if x['_uid'] in records_by_uid]
    if not calculated:
        return [], [], [], []
    reference = calculated[0]
    ref_chain = reference.get('_primary_chain', '')
    ref_seq = chain_data.get(reference['PDB'], {}).get(ref_chain, [])
    if not ref_seq:
        return [], [], ['Не удалось прочитать последовательность опорной цепи.'], []

    refkeys = {r['key'] for r in ref_seq}
    maps = {}
    included = []
    errors = []
    for lig in calculated:
        chain = lig.get('_primary_chain', '')
        if lig['_uid'] == reference['_uid']:
            mapping = {r['key']:r['key'] for r in ref_seq}
        else:
            other_seq = chain_data.get(lig['PDB'],{}).get(chain,[])
            aln = sequence_alignment(ref_seq, other_seq)
            if not aln['ok'] or aln['identity']<0.9 or aln['coverage']<0.8:
                errors.append(f"{lig['PDB']}: белок недостаточно похож на опорную мишень для строгого сравнения (идентичность {aln['identity']:.0%}, покрытие {aln['coverage']:.0%}).")
                continue
            mapping = aln['mapping']
        maps[lig['_uid']] = mapping
        included.append(lig)

    # ProLIF возвращает конкретные взаимодействия в формате RESNUM.CHAIN.
    mapped = {}
    for lig in included:
        contacts = defaultdict(set)
        unparsed = 0
        off_chain = 0
        for rec in records_by_uid[lig['_uid']]:
            parsed = parse_residue_id(str(rec.get('Остаток','')))
            if parsed is None:
                unparsed += 1
                continue
            key, chain = parsed
            if chain != lig.get('_primary_chain', ''):
                off_chain += 1
                continue
            ref_key = maps[lig['_uid']].get(key)
            if ref_key in refkeys:
                contacts[ref_key].add(classify_interaction(str(rec.get('Взаимодействие',''))))
        if unparsed:
            errors.append(f"{lig['PDB']}: {unparsed} записей ProLIF с неизвестным обозначением остатка пропущены.")
        if off_chain:
            errors.append(f"{lig['PDB']}: {off_chain} записей другой белковой цепи не включены в сравнение.")
        mapped[lig['_uid']] = contacts

    all_keys = set().union(*(m.keys() for m in mapped.values())) if mapped else set()
    ref_idx = {r['key']:i for i,r in enumerate(ref_seq)}
    labels = {lig['_uid']:f"{lig['PDB']} · {lig['Код компонента']}" for lig in included}
    table, shared = [], []
    for key in sorted(all_keys,key=lambda k:(ref_idx.get(k,100000),k)):
        row = {'Остаток': readable_residue(key)}
        interactions = []
        for lig in included:
            found = mapped[lig['_uid']].get(key,set())
            interactions.append(found)
            row[labels[lig['_uid']]] = ', '.join(sorted(found,key=lambda k:(TYPES.index(k) if k in TYPES else 999,k))) if found else '—'
        table.append(row)
        if len(included) >= 2 and all(interactions):
            common = set.intersection(*interactions)
            if common:
                shared.append({'Остаток':readable_residue(key),
                               'Общие взаимодействия':', '.join(sorted(common,key=lambda k:(TYPES.index(k) if k in TYPES else 999,k)))})
    return table,shared,errors,included
