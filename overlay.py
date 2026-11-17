"""Visual superposition of selected ligands after fitting protein C-alpha atoms."""
from __future__ import annotations
from html import escape
import json
from Bio.PDB import Superimposer
from analysis import load_model, sequence_alignment
from viewer import complex_pdb

COLORS=['#3D7395','#BC8B55','#9577B2','#71998C','#AE7780']


def _hetero(model, ligand):
    for chain in model:
        if chain.id != ligand['Цепь']:
            continue
        for res in chain:
            if res.id[0].startswith('H_') and res.resname == ligand['Код компонента'] and res.id[1] == ligand['Номер'] and str(res.id[2]).strip() == ligand.get('Вставка',''):
                return res
    raise ValueError('лиганд отсутствует в данной модели')


def _transformed_pdb(model, ligand, transform):
    rot,trans=transform
    res=_hetero(model,ligand)
    rows=[]
    for idx,atom in enumerate(res,1):
        if atom.element in ('H','D'):continue
        x,y,z=atom.coord @ rot + trans
        rows.append(f'HETATM{idx:5d} {atom.name[:4]:4s} {res.resname:>3s} {ligand["Цепь"][:1]:1s}'
                    f'{res.id[1]:4d}{str(res.id[2]):1s}   {x:8.3f}{y:8.3f}{z:8.3f}'
                    f'{1.00:6.2f}{0.00:6.2f}          {atom.element:>2s}')
    return '\n'.join(rows)+'\nEND\n'


def prepare_overlay(selected, cifs, chain_data):
    """Returns actual ligand coordinates transformed to the first protein frame."""
    ref=selected[0]
    rm=load_model(ref['PDB'],cifs[ref['PDB']])
    rc=ref['_primary_chain']
    ref_seq=chain_data[ref['PDB']][rc]
    ref_pdb=complex_pdb(ref['PDB'],cifs[ref['PDB']],ref)
    layers=[dict(label=ref['PDB']+' · '+ref['Код компонента'],color=COLORS[0],
                 reference=True,code=ref['Код компонента'],chain=ref['Цепь'],number=ref['Номер'],pdb='')]
    errors=[]
    for lig in selected[1:]:
        try:
            chain=lig['_primary_chain']
            other_seq=chain_data[lig['PDB']][chain]
            aln=sequence_alignment(ref_seq,other_seq)
            if not aln['ok']:raise ValueError('цепи недостаточно похожи')
            other=load_model(lig['PDB'],cifs[lig['PDB']])
            fixed,moving=[],[]
            for i,j in aln['pairs']:
                rr=ref_seq[i]; orr=other_seq[j]
                rid=(' ',rr['number'],rr['insertion'] or ' ')
                oid=(' ',orr['number'],orr['insertion'] or ' ')
                if rid in rm[rc] and oid in other[chain]:
                    if 'CA' in rm[rc][rid] and 'CA' in other[chain][oid]:
                        fixed.append(rm[rc][rid]['CA']);moving.append(other[chain][oid]['CA'])
            if len(fixed)<3:raise ValueError('мало соответствующих атомов Cα')
            sup=Superimposer();sup.set_atoms(fixed,moving)
            layers.append(dict(label=lig['PDB']+' · '+lig['Код компонента'],
                color=COLORS[(len(layers))%len(COLORS)],reference=False,
                pdb=_transformed_pdb(other,lig,sup.rotran)))
        except (KeyError,ValueError,IndexError) as exc:
            errors.append(f"{lig['PDB']}: не удалось наложить лиганд: {exc}")
    return ref_pdb,layers,errors


def overlay_html(reference_pdb, layers, common_residues, height=580):
    """HTML 3D view: one reference protein, aligned inhibitor molecules, shared contact positions."""
    pdb=json.dumps(reference_pdb,ensure_ascii=True).replace('<','\\u003c')
    layers_j=json.dumps(layers,ensure_ascii=True).replace('<','\\u003c')
    common_j=json.dumps(common_residues,ensure_ascii=True)
    view_h=height-83
    return f'''<!DOCTYPE html><html lang="ru"><head><meta charset="utf-8"/>
      <meta name="viewport" content="width=device-width, initial-scale=1"/>
      <style>
       html,body{{margin:0;padding:0;background:#FBFAF7;color:#334C3B;font:12px system-ui,sans-serif}}
       #legend{{min-height:38px;display:flex;gap:12px;flex-wrap:wrap;align-items:center;padding:5px}}
       #legend label{{display:flex;gap:5px;align-items:center;cursor:pointer}}
       #mol{{height:{view_h}px;width:100%;border:1px solid #DDE6DA;border-radius:11px;box-sizing:border-box;position:relative}}
       #status{{position:absolute;top:10px;left:10px;background:#F4E9D8;padding:8px;z-index:1;border-radius:7px}}
       #hint{{font-size:12px;color:#687D6E;padding:8px 4px}}
      </style></head><body><div id="legend"></div>
      <div id="mol"><div id="status">Загружаем сравнение…</div></div>
      <div id="hint">Левая кнопка — вращать · колесо — масштаб · правая кнопка — перемещать.</div>
      <script>
       const pdb={pdb}, layers={layers_j}, common={common_j};
       function launch(){{
        if(!window.$3Dmol){{document.getElementById('status').textContent='3D-библиотека не загрузилась';return;}}
        try{{
         const view=$3Dmol.createViewer('mol',{{backgroundColor:'#FBFAF7'}});
         const base=view.addModel(pdb,'pdb');
         base.setStyle({{hetflag:false}},{{}});
         common.forEach(r=>view.addStyle({{hetflag:false,chain:r.chain,resi:r.resi}},
                          {{stick:{{radius:0.15,color:'#52795B'}}}}));
         const proteinBox=document.createElement('label');
         const proteinCb=document.createElement('input');proteinCb.type='checkbox';proteinCb.checked=false;
         proteinCb.onchange=()=>{{base.setStyle({{hetflag:false}},proteinCb.checked?{{cartoon:{{color:'#BFD1C2',opacity:0.25}}}}:{{}});
            common.forEach(r=>view.addStyle({{hetflag:false,chain:r.chain,resi:r.resi}},{{stick:{{radius:0.15,color:'#52795B'}}}}));view.render();}};
         proteinBox.append(proteinCb,document.createTextNode('Показать белок'));document.getElementById('legend').append(proteinBox);
         const refsel={{hetflag:true,resn:layers[0].code,chain:layers[0].chain,resi:layers[0].number}};
         view.addStyle(refsel,{{stick:{{radius:0.23,color:layers[0].color}}}});
         const models=[base];
         layers.slice(1).forEach(layer=>{{const model=view.addModel(layer.pdb,'pdb');
                 model.setStyle({{}},{{stick:{{radius:0.23,color:layer.color}}}});models.push(model);}});
         layers.forEach((layer,i)=>{{
          const label=document.createElement('label');const cb=document.createElement('input');
          cb.type='checkbox';cb.checked=true;
          cb.onchange=()=>{{if(i===0)base.setStyle(refsel,cb.checked?{{stick:{{radius:0.23,color:layer.color}}}}:{{}});
                         else models[i].setStyle({{}},cb.checked?{{stick:{{radius:0.23,color:layer.color}}}}:{{}});
                         view.render();}};
          const dot=document.createElement('span');dot.style.cssText='width:10px;height:10px;border-radius:50%;background:'+layer.color;
          label.append(cb,dot,document.createTextNode(layer.label));document.getElementById('legend').append(label);
         }});
         view.zoomTo(refsel);view.zoom(0.58);view.render();
         document.getElementById('status').style.display='none';
        }}catch(e){{document.getElementById('status').textContent='Не удалось построить наложение: '+e.message;}}
       }}
       function load(src,backup){{const s=document.createElement('script');s.src=src;s.onload=launch;
        s.onerror=()=>backup?load(backup,null):document.getElementById('status').textContent='Нет доступа к 3Dmol.js';
        document.head.appendChild(s);}}
       load('https://3Dmol.org/build/3Dmol-min.js','https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js');
      </script></body></html>'''
