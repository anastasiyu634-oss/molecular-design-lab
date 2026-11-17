"""3Dmol viewer for the selected real PDB complex; 2D structures from RDKit."""
from __future__ import annotations

from html import escape
from io import StringIO
import json
import re

from Bio.PDB import PDBIO, Select, is_aa
from analysis import load_model


class SelectedComplex(Select):
    def __init__(self, ligand):
        self.ligand = ligand

    def accept_residue(self, residue):
        if residue.id[0] == ' ' and is_aa(residue, standard=False):
            return True
        ligand = self.ligand
        return int(residue.resname == ligand['Код компонента'] and
                   residue.get_parent().id == ligand['Цепь'] and
                   residue.id[1] == ligand['Номер'] and
                   str(residue.id[2]).strip() == ligand.get('Вставка', ''))


def complex_pdb(pdb_id: str, cif_text: str, ligand: dict) -> str:
    model = load_model(pdb_id, cif_text)
    io = PDBIO()
    io.set_structure(model)
    buf = StringIO()
    io.save(buf, select=SelectedComplex(ligand))
    return buf.getvalue()


def _pdb_selection(residue_row):
    key = residue_row['_residue_key']
    match = re.fullmatch(r'([A-Z]{3}):(-?\d+)([A-Z]?)', key)
    if match is None:
        return None
    return {'chain': residue_row['Цепь белка'], 'resi': int(match.group(2)),
            'shared_key': key}


def viewer_html(pdb_text: str, ligand: dict, height: int = 420,
                mode: str = 'full', highlight_keys: set[str] | None = None) -> str:
    """Interactive 3D viewer with a ligand-focused pocket mode and highlighted residues."""
    safe_pdb = json.dumps(pdb_text, ensure_ascii=True).replace('<', '\\u003c').replace('>', '\\u003e')
    code = json.dumps(ligand['Код компонента'])
    chain = json.dumps(ligand['Цепь'])
    num = json.dumps(int(ligand['Номер']))
    highlights = set(highlight_keys or [])
    nearby = []
    for row in ligand.get('_contact_rows', [])[:100]:
        obj = _pdb_selection(row)
        if obj:
            obj['shared'] = row['_residue_key'] in highlights
            nearby.append(obj)
    nearby_json = json.dumps(nearby, ensure_ascii=True)
    start_pocket = json.dumps(mode == 'pocket')
    view_height = max(360, height - 130)
    return f'''<!doctype html><html lang="ru"><head><meta charset="utf-8" />
    <style>
       html,body{{margin:0;padding:0;background:#F8F7F3;color:#344A3A;font:12px system-ui,sans-serif;}}
       .tools{{min-height:58px;display:flex;flex-wrap:wrap;align-items:center;justify-content:space-between;
           padding:5px 6px;gap:9px;box-sizing:border-box;}}
       .tools button{{font:inherit;cursor:pointer;border:1px solid #D8E1D7;border-radius:7px;
           background:#fff;color:#3E5546;padding:6px 10px;}}
       .gesture{{font-size:12px;color:#596D60;line-height:1.5;padding:4px 7px;}}
       #view{{height:{view_height}px;width:100%;border-radius:10px;
           border:1px solid #DFE7DE;box-sizing:border-box;position:relative;overflow:hidden;background:#F8F7F3;}}
       #label{{color:#3F5947;height:46px;padding:8px 8px 3px;line-height:1.35;box-sizing:border-box;}}
       #status{{position:absolute;left:10px;top:8px;z-index:3;background:#fffdd9;padding:5px;border-radius:6px;}}
    </style></head><body>
      <div class="tools"><div><button id="pocket">Область связывания</button> <button id="full">Белок целиком</button></div>
        <span class="gesture">Левая кнопка + движение — вращать · колесо — масштаб · правая кнопка + движение — сдвигать</span>
      </div>
      <div id="view"><span id="status">Загружаем 3D…</span></div>
      <div id="label">Нажмите на атом белка, чтобы увидеть аминокислоту и расстояние до ближайшего атома лиганда.</div>
      <script>
      const modelText = {safe_pdb};
      const near = {nearby_json};
      const ligandSel = {{hetflag:true,resn:{code},chain:{chain},resi:{num}}};
      let viewer, model;
      function launch() {{
        if(!window.$3Dmol) {{
          document.getElementById('status').textContent = '3Dmol.js не загрузился. Проверь доступ к интернету.';
          return;
        }}
        try {{
          viewer=$3Dmol.createViewer('view',{{backgroundColor:'#F8F7F3'}});
          model=viewer.addModel(modelText,'pdb');
          document.getElementById('status').style.display='none';
          const atoms=model.selectedAtoms(ligandSel);
          if(!atoms.length) {{
            document.getElementById('label').textContent = 'Не удалось выделить лиганд. Проверь выбранную молекулу.';
          }}
          function apply(pocket) {{
            viewer.removeAllSurfaces();
            viewer.setStyle({{}},{{}});
            viewer.setStyle({{hetflag:false}},{{cartoon:{{color:'#9CB4A2',opacity:pocket?0.22:0.9}}}});
            near.forEach(r => viewer.addStyle(
              {{hetflag:false,chain:r.chain,resi:r.resi}},
              {{stick:{{radius:pocket?0.16:0.13,color:r.shared?'#6B9C76':'#C5A16D'}}}}
            ));
            viewer.addStyle(ligandSel,{{stick:{{radius:0.22,colorscheme:'Jmol'}}}});
            viewer.setClickable({{hetflag:false}},true,function(atom){{
              if(!atom)return;
              let dist=Infinity,closest=null;
              atoms.forEach(a => {{const d=Math.hypot(atom.x-a.x,atom.y-a.y,atom.z-a.z);
                if(d<dist){{dist=d;closest=a;}}}});
              document.getElementById('label').textContent =
                atom.resn+' '+atom.resi+' · цепь '+atom.chain+' · атом '+atom.atom+
                (closest?' · ближайший атом лиганда: '+closest.atom+' · '+dist.toFixed(2)+' Å':'');
              viewer.removeAllLabels();
              viewer.addLabel(atom.resn+' '+atom.resi,{{position:atom,backgroundColor:'#3E5546',fontColor:'white',fontSize:11}});
              viewer.render();
            }});
            if(pocket) {{viewer.zoomTo(ligandSel);viewer.zoom(0.60);}}
            else {{viewer.zoomTo({{hetflag:false}});viewer.zoom(0.83);}}
            viewer.render();
          }}
          document.getElementById('pocket').addEventListener('click',()=>apply(true));
          document.getElementById('full').addEventListener('click',()=>apply(false));
          apply({start_pocket});
        }} catch(err) {{document.getElementById('status').textContent='Ошибка 3D: '+err.message;}}
      }}
      function loadLibrary(url, fallback) {{
         const s=document.createElement('script');s.src=url;s.onload=launch;
         s.onerror=()=>fallback?loadLibrary(fallback,null):
           document.getElementById('status').textContent='Библиотека 3D недоступна. Проверь интернет.';
         document.head.appendChild(s);
      }}
      loadLibrary('https://3Dmol.org/build/3Dmol-min.js',
                  'https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js');
      </script></body></html>'''


def ligand_svg(smiles: str) -> str | None:
    if not smiles:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(640, 360)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


def responsive_ligand_html(svg: str, expanded: bool = False) -> str:
    if not svg or '<svg' not in svg:
        return ''
    svg = svg[svg.find('<svg'):]
    svg = re.sub(r"width=['\"]\d+px['\"]", "width='100%'", svg, count=1)
    svg = re.sub(r"height=['\"]\d+px['\"]", "height='100%'", svg, count=1)
    return ('<!doctype html><html><head><meta charset="utf-8" />'
            '<style>html,body{height:100%;margin:0;background:#fff;overflow:hidden}'
            '.chemical{height:100vh;width:100%;display:flex;align-items:center;justify-content:center}'
            'svg{display:block;width:100% !important;height:95% !important;max-height:100%}'
            '</style></head><body><div class="chemical">'+svg+'</div></body></html>')


def sdf_viewer_html(sdf_text: str, height: int = 450) -> str:
    """Display a generated SDF conformer interactively; this is not a docked pose."""
    safe_sdf = (json.dumps(sdf_text, ensure_ascii=True)
                .replace('<', r'\u003c').replace('>', r'\u003e').replace('&', r'\u0026'))
    view_h = max(300, int(height) - 60)
    return f'''<!doctype html><html lang="ru"><head><meta charset="utf-8" />
    <style>
      html,body{{margin:0;padding:0;background:#F8F7F3;font:13px system-ui,sans-serif;color:#3E5546}}
      #controls{{display:flex;gap:8px;align-items:center;flex-wrap:wrap;padding:7px 4px}}
      button{{font:inherit;border:1px solid #D9E3D9;background:#fff;color:#385242;padding:7px 11px;border-radius:7px;cursor:pointer}}
      #view{{height:{view_h}px;width:100%;background:#F8F7F3;border:1px solid #DFE7DE;border-radius:10px;position:relative;box-sizing:border-box;overflow:hidden}}
      #status{{position:absolute;top:10px;left:10px;background:#F9F0DD;padding:7px;border-radius:6px;z-index:2}}
    </style></head><body>
    <div id="controls">
      <button type="button" id="sticks">Палочки</button>
      <button type="button" id="balls">Шары и палочки</button>
      <button type="button" id="fit">Показать целиком</button>
      <span style="color:#607566">Мышь — вращать · колесо — масштаб · правая кнопка — двигать</span>
    </div>
    <div id="view"><span id="status">Загружаем 3D-молекулу…</span></div>
    <script>
      const sdf={safe_sdf};
      function showError(msg){{document.getElementById('status').textContent=msg;}}
      function launch(){{
        if(!window.$3Dmol){{showError('Не удалось загрузить 3Dmol.js. Проверь интернет.');return;}}
        try{{
          const view=$3Dmol.createViewer('view',{{backgroundColor:'#F8F7F3'}});
          const model=view.addModel(sdf,'sdf');
          if(!model || !model.selectedAtoms({{}}).length){{showError('Не удалось прочитать 3D-координаты SDF.');return;}}
          function style(balls){{
            view.setStyle({{}}, balls
              ? {{stick:{{radius:0.13}}, sphere:{{scale:0.28}}}}
              : {{stick:{{radius:0.2}}}});
            view.render();
          }}
          style(true);
          view.zoomTo();view.zoom(0.86);view.render();
          document.getElementById('status').style.display='none';
          document.getElementById('sticks').onclick=()=>style(false);
          document.getElementById('balls').onclick=()=>style(true);
          document.getElementById('fit').onclick=()=>{{view.zoomTo();view.zoom(0.86);view.render();}};
        }}catch(e){{showError('Ошибка отображения 3D: '+e.message);}}
      }}
      function load(url,fallback){{
        const s=document.createElement('script');s.src=url;s.onload=launch;
        s.onerror=()=>fallback?load(fallback,null):showError('Библиотека 3Dmol.js не загрузилась. Проверь интернет.');
        document.head.appendChild(s);
      }}
      load('https://3Dmol.org/build/3Dmol-min.js','https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js');
    </script></body></html>'''


def docked_pose_viewer_html(protein_pdb: str, poses_pdbqt: str) -> str:
    """Render full docked ligands (SDF bond graph) and let users select poses."""
    from docking_pose_format import pdbqt_poses_to_sdf, PoseFormatError
    try:
        payload=pdbqt_poses_to_sdf(poses_pdbqt)
        error=''
    except PoseFormatError as exc:
        payload=[]
        error=str(exc)
    safe_protein=json.dumps(protein_pdb,ensure_ascii=True).replace('<',r'\u003c').replace('>',r'\u003e').replace('&',r'\u0026')
    safe_payload=json.dumps(payload,ensure_ascii=True).replace('<',r'\u003c').replace('>',r'\u003e').replace('&',r'\u0026')
    safe_error=json.dumps(error,ensure_ascii=True).replace('<',r'\u003c').replace('>',r'\u003e')
    return r'''<!doctype html><html lang="ru"><head><meta charset="utf-8"><style>
    html,body{margin:0;font:13px system-ui,sans-serif;color:#3E5546;background:#F8F7F3}
    .tools{min-height:43px;display:flex;gap:8px;align-items:center;flex-wrap:wrap;padding:5px 3px}
    select,button{font:inherit;padding:6px 10px;border:1px solid #CCD8CC;border-radius:7px;background:white;color:#385242;cursor:pointer}
    #score{font-size:12px;color:#566A5D}
    #view{height:430px;width:100%;border:1px solid #DFE7DE;border-radius:10px;position:relative;overflow:hidden;box-sizing:border-box}
    #status{position:absolute;top:10px;left:10px;background:#FFF4E2;padding:7px;border-radius:6px;z-index:4;max-width:95%}
    .help{padding:6px;color:#637365;font-size:12px}
    </style></head><body>
    <div class="tools"><label for="poses">Положение молекулы:</label>
      <select id="poses"></select>
      <button id="fit" type="button">Показать молекулу целиком</button>
      <button id="protein" type="button">Скрыть белок</button>
      <span id="score"></span>
    </div>
    <div id="view"><div id="status">Загружаем рассчитанную позу…</div></div>
    <div class="help">Левая кнопка — вращать · колёсико — приближать · правая кнопка — сдвигать. Это расчётные позы, не экспериментальные.</div>
    <script>
    // The input is a PDBQT file, but SDF is used for bond-complete rendering.
    const rec=__REC__,poses=__POSES__,initialError=__ERROR__;
    const status=document.getElementById('status');
    function fail(msg){status.textContent=msg;status.style.display='block';}
    function start(){
      if(!poses.length){fail('Не удалось восстановить 3D-лиганд: '+initialError);return;}
      try{
        const v=$3Dmol.createViewer('view',{backgroundColor:'#F8F7F3'});
        const receptor=v.addModel(rec,'pdb');
        receptor.setStyle({},{cartoon:{color:'#B0C0B2',opacity:.5}});
        const selector=document.getElementById('poses');
        poses.forEach((pose,i)=>{
          const opt=document.createElement('option');opt.value=String(i);
          opt.textContent='Поза '+pose.index;selector.appendChild(opt);
        });
        let ligand=null,proteinVisible=true;
        function paint(i,zoom=true){
          if(ligand)v.removeModel(ligand);
          const pose=poses[i];
          ligand=v.addModel(pose.sdf,'sdf');
          const atoms=ligand.selectedAtoms({});
          if(atoms.length!==pose.heavy_atoms)throw new Error('Просмотрщик отобразил не все атомы лиганда.');
          ligand.setStyle({},{stick:{radius:.22,colorscheme:'Jmol'}});
          if(zoom){v.zoomTo({model:ligand.getID()});v.zoom(.75);}
          v.render();
          document.getElementById('score').textContent=pose.score===null?'':'Оценка Vina: '+pose.score.toFixed(3)+' ккал/моль';
          status.style.display='none';
        }
        selector.onchange=()=>{try{paint(Number(selector.value));}catch(e){fail('Ошибка отображения: '+e.message);}};
        document.getElementById('fit').onclick=()=>{v.zoomTo({model:ligand.getID()});v.zoom(.75);v.render();};
        document.getElementById('protein').onclick=()=>{
          proteinVisible=!proteinVisible;
          receptor.setStyle({},proteinVisible?{cartoon:{color:'#B0C0B2',opacity:.5}}:{});
          document.getElementById('protein').textContent=proteinVisible?'Скрыть белок':'Показать белок';
          v.render();
        };
        paint(0);
      }catch(e){fail('Ошибка 3D: '+e.message);}
    }
    function load(url,fallback){const script=document.createElement('script');script.src=url;script.onload=start;
      script.onerror=()=>fallback?load(fallback,null):fail('3Dmol.js не загрузился; проверь интернет.');
      document.head.appendChild(script);}
    load('https://3Dmol.org/build/3Dmol-min.js','https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js');
    </script></body></html>'''.replace('__REC__',safe_protein).replace('__POSES__',safe_payload).replace('__ERROR__',safe_error)
