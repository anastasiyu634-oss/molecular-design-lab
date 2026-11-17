"""Select a *binding region*, not merely a protein chain.

Sites are anchored either to an experimentally observed hetero molecule
from this exact coordinate model, or to a user-specified residue set.
No site discovery is claimed when the PDB lacks a reference ligand.
"""
from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any


def observed_sites(pdb_id: str, ligands: list[dict], primary_chain: str = '') -> list[dict]:
    """Stable, transparent choices backed by actual PDB ligand coordinates."""
    rows = []
    for item in ligands:
        if item.get('PDB') != pdb_id or item.get('Категория') != 'Лиганд-кандидат':
            continue
        if not item.get('_primary_chain'):
            continue
        if primary_chain and str(item['_primary_chain']) != str(primary_chain):
            continue
        # Multiple sites in the same PDB can target different chains.
        number = int(item['Номер'])
        rows.append({
            'key': str(item.get('_uid') or f"{pdb_id}:{item['Код компонента']}:{item['Цепь']}:{number}"),
            'ligand': item,
            'chain': str(item['_primary_chain']),
            'label': (f"Область около {item['Код компонента']} · "
                      f"цепь белка {item['_primary_chain']} · молекула №{number}"),
            'type': 'observed',
        })
    # Keep all real sites; no arbitrary selection based solely on PDB chain.
    return sorted(rows, key=lambda x: (x['chain'] != primary_chain, x['label']))


def _parse_site_residues(query: str, default_chain: str) -> list[tuple[str, int, str]]:
    """Accept 'GLU353, HIS524', 'A:353, A:524', or bare numbered residues."""
    tokens = [v.strip() for v in re.split(r'[,;\n]+', query.strip()) if v.strip()]
    if not tokens or len(tokens) > 20:
        raise ValueError('Укажи от 1 до 20 остатков, например GLU353, HIS524.')
    result = []
    for item in tokens:
        m = re.fullmatch(r'(?:(?P<chain>[A-Za-z0-9]):)?(?:[A-Za-z]{3}\s*)?(?P<num>-?\d+)(?P<ins>[A-Za-z]?)', item)
        if not m:
            raise ValueError(f'Не удалось прочитать остаток «{item}». Пример: A:353, HIS524.')
        result.append((m.group('chain') or default_chain, int(m.group('num')), (m.group('ins') or '').upper()))
    return list(dict.fromkeys(result))


def box_from_residues(protein_pdb: str, residues: str, chain: str,
                      size: float = 22.0) -> dict[str, Any]:
    """A conservative local box centred on existing atoms of named residues.

    These are NOT automatically discovered catalytic residues. The user
    explicitly selects a region; absent residues are a hard error.
    """
    ids = _parse_site_residues(residues, chain)
    atoms: dict[tuple[str, int, str], list[tuple[float, float, float]]] = defaultdict(list)
    for line in protein_pdb.splitlines():
        if not line.startswith('ATOM  ') or len(line) < 54:
            continue
        try:
            key = (line[21].strip() or '?', int(line[22:26]), line[26].strip().upper())
            element = line[76:78].strip().upper() if len(line) >= 78 else ''
            if element in {'H', 'D'}:
                continue
            atoms[key].append(tuple(float(line[p:p+8]) for p in (30, 38, 46)))
        except ValueError:
            continue
    absent = [f'{ch}:{num}{ins}' for ch,num,ins in ids if not atoms.get((ch,num,ins))]
    if absent:
        raise ValueError('Эти остатки отсутствуют в выбранной модели белка: ' + ', '.join(absent))
    points = [v for key in ids for v in atoms[key]]
    center = [round(sum(p[i] for p in points)/len(points), 2) for i in range(3)]
    edges = [max(p[i] for p in points)-min(p[i] for p in points)+12.0 for i in range(3)]
    if any(e > 40.0 for e in edges):
        raise ValueError('Выбранные остатки расположены слишком далеко друг от друга для одной области докинга. Укажи более близкие остатки.')
    side = [round(max(18.0, float(size), e), 2) for e in edges]
    return {'center':center,'size':side,'source':'residues','residues':ids}
