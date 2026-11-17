"""Combine geometric proximity and confirmed-by-fingerprint ProLIF predictions.

These are distinct computed measures, not proof of experimental chemical bonding.
ProLIF not attempted is *never* interpreted as zero interactions.
"""
from __future__ import annotations
from collections import defaultdict
import re
from interaction_comparison import parse_residue_id


def proximity_key(row):
    chain=str(row.get('Цепь белка','')).strip()
    name=str(row.get('Остаток','')).strip()
    # GLU A:353 or GLU A:353B
    m=re.fullmatch(r'([A-Za-z]{3})\s+([A-Za-z0-9]+):(-?\d+[A-Za-z]?)',name)
    if m:
        return f'{m[1].upper()}:{m[3].upper()}',m[2]
    # internal canonical key from analysis.py: GLU:353
    k=str(row.get('_residue_key','')).upper()
    return (k,chain) if k and chain else None


def unified_contact_rows(proximity:list[dict], records:list[dict]|None, chain:str|None=None):
    """Union all residue positions, including ProLIF-only >4A contacts.

    records = None => not calculated. [] => calculated, none found.
    """
    contacts={}
    for g in proximity:
        k=proximity_key(g)
        if not k or (chain and k[1]!=chain):
            continue
        contacts.setdefault(k,{'name':k[0].replace(':',''),'chain':k[1], 'geo':None,'prolif':[]})['geo']=g
    if records is not None:
        for rec in records:
            parsed=parse_residue_id(str(rec.get('Остаток','')))
            if not parsed:
                continue
            key,ch=parsed
            if chain and chain!=ch:
                continue
            d=contacts.setdefault((key,ch),{'name':key.replace(':',''),'chain':ch,'geo':None,'prolif':[]})
            d['prolif'].append(rec)
    out=[]
    for _k,r in contacts.items():
        geo=r['geo']; found=r['prolif']; types=list(dict.fromkeys(str(x['Взаимодействие']) for x in found))
        dist=[x.get('Расстояние (Å)') for x in found if isinstance(x.get('Расстояние (Å)'),(int,float))]
        out.append({
            'Остаток':r['name']+' · цепь '+r['chain'],
            'Ближайшие атомы, Å':geo.get('Минимальная дистанция (Å)') if geo else None,
            'Контакт по ProLIF':('Не рассчитано' if records is None else ', '.join(types) if types else 'Не обнаружено'),
            'Расстояние ProLIF, Å':min(dist) if dist else None,
            'Атомы (лиганд / белок)':(
                str(geo.get('Ближайший атом лиганда',''))+' / '+str(geo.get('Ближайший атом белка',''))
                if geo else '—')
        })
    def sorter(v):
        m=re.search(r'(-?\d+)',v['Остаток'])
        return (int(m.group()) if m else 999999,v['Остаток'])
    return sorted(out,key=sorter)
