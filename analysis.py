"""Первичный анализ PDB/mmCIF: загрузка структур и близкие остатки белка.

Не является фармакофорным анализом, докингом или подтверждением активности.
"""
from __future__ import annotations

import re
from io import StringIO

import requests
from Bio.PDB import MMCIFParser, NeighborSearch, is_aa


def parse_pdb_ids(raw: str, max_count: int = 10) -> list[str]:
    """Проверить пользовательский список PDB ID; порядок сохраняется."""
    ids = list(dict.fromkeys(p.strip().upper() for p in re.split(r"[,;\s]+", raw) if p.strip()))
    if not ids:
        raise ValueError("Введите хотя бы один идентификатор PDB.")
    invalid = [p for p in ids if not re.fullmatch(r"[A-Z0-9]{4}", p)]
    if invalid:
        raise ValueError("Проверьте идентификаторы PDB: " + ", ".join(invalid) + ". Нужны четырёхсимвольные коды.")
    if len(ids) > max_count:
        raise ValueError(f"Пока можно загрузить не более {max_count} структур за раз.")
    return ids


def download_structure(pdb_id: str) -> str:
    """Скачать структуру mmCIF из RCSB PDB."""
    resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=35)
    resp.raise_for_status()
    return resp.text


def analyze_structure(pdb_id: str, cif_text: str, radius: float = 4.0) -> tuple[list[dict], dict]:
    """Крупные неводные HET-компоненты и близкие аминокислотные остатки.

    Расстояние отсчитывается между тяжёлыми атомами; один остаток учитывается
    один раз, даже если найдено много близких атомов. При нескольких моделях
    анализируется только первая. Данные не классифицируют химические связи.
    """
    structure = MMCIFParser(QUIET=True).get_structure(pdb_id, StringIO(cif_text))
    try:
        model = next(structure.get_models())
    except StopIteration as exc:
        raise ValueError("В структуре нет координатных моделей.") from exc

    protein_residues = [
        residue for chain in model for residue in chain
        if residue.id[0] == " " and is_aa(residue, standard=False)
    ]
    protein_atoms = [
        atom for residue in protein_residues for atom in residue
        if atom.element not in {"H", "D"}
    ]
    neighbor_search = NeighborSearch(protein_atoms) if protein_atoms else None

    ligands = []
    for chain in model:
        for residue in chain:
            if not residue.id[0].startswith("H_"):
                continue
            atoms = [a for a in residue if a.element not in {"H", "D"}]
            if len(atoms) < 6:
                continue  # Эвристика: отсеиваем воду, ионы, небольшие добавки.

            neighbors: dict[str, str] = {}
            if neighbor_search is not None:
                for atom in atoms:
                    for nearby in neighbor_search.search(atom.coord, radius, level="R"):
                        chain_id = nearby.get_parent().id
                        seq_id = nearby.id[1]
                        ins = nearby.id[2].strip()
                        key = f"{nearby.resname}:{seq_id}{ins}"
                        label = f"{nearby.resname} {chain_id}:{seq_id}{ins}"
                        neighbors[label] = key

            ligands.append({
                "PDB": pdb_id,
                "Код компонента": residue.resname,
                "Цепь": chain.id,
                "Номер": residue.id[1],
                "Тяжёлых атомов": len(atoms),
                "Близких остатков (≤ 4 Å)": len(neighbors),
                "Ближайшие остатки": ", ".join(sorted(neighbors)),
                "_residue_keys": sorted(set(neighbors.values())),
            })

    ligands.sort(key=lambda r: (r["Тяжёлых атомов"], r["Близких остатков (≤ 4 Å)"]), reverse=True)
    return ligands, {
        "Белковых остатков": len(protein_residues),
        "Белковых цепей": len({r.get_parent().id for r in protein_residues}),
    }


def common_residues(ligands: list[dict]) -> list[str]:
    """Грубое совпадение по имени/номеру остатка БЕЗ выравнивания белков."""
    if len(ligands) < 2:
        return []
    sets = [set(ligand.get("_residue_keys", [])) for ligand in ligands]
    return sorted(set.intersection(*sets))
