"""Stage 03: exploratory pharmacophore and fragment-design, stage 04: export for docking."""
from __future__ import annotations
import json
from html import escape
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from chem_design import conformer_sdf, generate_candidates, generate_analogs, properties, protein_only_pdb
from pharmacophore import compute_pharmacophore, FEATURE_COLORS, FEATURE_NAMES
from viewer import sdf_viewer_html


def _pharmacophore_html(complex_pdb, features, height=510):
    """Selected ligand + pharmacophore hypothesis in experimentally aligned frame."""
    pdb=json.dumps(complex_pdb).replace('<','\\u003c')
    fs=json.dumps([{'name':v['Признак'],'family':v['Семейство'],
                    'xyz':v['_xyz'],'color':v['_color']} for v in features],ensure_ascii=False).replace('<','\\u003c')
    html = r'''<!doctype html><html lang="ru"><head><meta charset="utf-8"/><style>
    body{margin:0;background:#F8F7F3;font:13px system-ui;color:#3E5546}
    #m{height:440px;width:100%;position:relative;border:1px solid #DFE7DE;border-radius:10px;overflow:hidden}
    #status{position:absolute;z-index:2;top:10px;left:10px;background:#F4E9D8;padding:7px;border-radius:6px}
    #legend{display:flex;gap:10px;flex-wrap:wrap;padding:8px 0}
    .item{display:flex;align-items:center;gap:5px;font-size:12px}
    .dot{height:10px;width:10px;border-radius:50%;display:inline-block}
    .help{padding:4px 2px;font-size:12px;color:#617367}
    </style></head><body><div id="m"><span id="status">Загрузка 3D-модели…</span></div>
    <div id="legend"></div><div class="help">Мышь: левая кнопка — вращение; колесо — масштаб; правая кнопка — перемещение. Сферы — расчетные центры признаков.</div>
    <script>const pdb=__PDB__,features=__FEATURES__;
    function launch(){try{const v=$3Dmol.createViewer('m',{backgroundColor:'#F8F7F3'});
    v.addModel(pdb,'pdb');v.setStyle({hetflag:false},{cartoon:{color:'#B5C7B7',opacity:0.27}});
    v.addStyle({hetflag:true},{stick:{radius:.19,colorscheme:'Jmol'}});
    const legend=document.getElementById('legend');let unique={};
    for(const f of features){const c=f.xyz;
      v.addSphere({center:{x:c[0],y:c[1],z:c[2]},radius:.78,color:f.color,alpha:.69});
      unique[f.family]={label:f.name,color:f.color};
    }
    for(const v of Object.values(unique)){
      const e=document.createElement('span');e.className='item';
      const d=document.createElement('span');d.className='dot';d.style.background=v.color;
      const t=document.createElement('span');t.textContent=v.label;e.append(d,t);legend.appendChild(e);
    }
    v.zoomTo({hetflag:true});v.zoom(.82);v.render();
    document.getElementById('status').style.display='none';
    }catch(e){document.getElementById('status').textContent='Ошибка 3D: '+e.message;}}
    function load(url,backup){const s=document.createElement('script');s.src=url;s.onload=launch;
    s.onerror=()=>backup?load(backup,null):document.getElementById('status').textContent='Библиотека 3Dmol недоступна';document.head.appendChild(s);}
    load('https://3Dmol.org/build/3Dmol-min.js','https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js');
    </script></body></html>'''
    return html.replace('__PDB__',pdb).replace('__FEATURES__',fs)


def render_design(selected, complex_fn, chem_info_fn, session_cifs, chain_data):
    st.markdown('<div class="lab-kicker">Этап 03 / Фармакофор и проектирование</div>', unsafe_allow_html=True)
    st.subheader('Поиск общих химических признаков')
    st.write('Исследуем, какие типы химических групп занимают сходные места **относительно одной белковой мишени**. Это предварительная гипотеза для проектирования, а не подтверждение активности.')
    if not selected:
        st.info('Сначала загрузите структуру белка с лигандом на этапе 01.')
        return
    if len(selected) == 1:
        st.info('Загружен один комплекс: покажем его химические признаки, но искать повторения пока не с чем.')
    else:
        st.caption(f'В анализе участвуют выбранные комплексы: {", ".join(x["PDB"] for x in selected)}.')
    dist = st.slider('Допуск для совпадения однотипных признаков (Å)',
                     min_value=1.0,max_value=3.0,value=2.0,step=0.25,
                     help='Максимальное расстояние между центрами признаков после выравнивания белковых структур. Это настроечный параметр, не доказанный порог фармакофора.')
    if st.button('Найти пространственные признаки', type='primary'):
        st.session_state['feature_attempted'] = True
        st.session_state['feature_points'] = []
        st.session_state['pharma_clusters'] = []
        st.session_state['pharma_errors'] = []
        try:
            with st.spinner('Сопоставляем структуры и ищем химические признаки…'):
                points, groups, errors = compute_pharmacophore(
                    selected,session_cifs,chain_data,complex_fn,chem_info_fn,dist)
            st.session_state['feature_points'] = points
            st.session_state['pharma_clusters'] = groups
            st.session_state['pharma_errors'] = errors
            st.session_state['pharma_tolerance'] = dist
        except Exception as exc:
            st.session_state['pharma_errors'] = [f'Ошибка расчёта: {type(exc).__name__}: {exc}']

    if st.session_state.get('feature_attempted',False):
        errors=st.session_state.get('pharma_errors',[])
        for error in errors:
            st.error(error)
        points=st.session_state.get('feature_points',[])
        groups=st.session_state.get('pharma_clusters',[])
        if not points:
            st.warning('Не удалось получить пространственные признаки. Это не означает, что у молекулы их нет. Проверьте ошибки выше и доступность химических данных PDB.')
            st.caption('Если хочешь, пришли мне текст сообщения об ошибке — по нему можно исправить конкретный этап подготовки.')
        else:
            a,b=st.columns(2)
            a.metric('Химических признаков в обработанных лигандах',len(points))
            b.metric('Пространственных совпадений',len(groups))
            if groups:
                st.markdown('**Фармакофорные точки в 3D**')
                st.caption('Сферы показывают пространственные центры химических признаков. Это предварительная гипотеза, а не проверенный фармакофор.')
                families=sorted({v['Семейство'] for v in groups})
                labels={x:FEATURE_NAMES.get(x,x) for x in families}
                shown=st.multiselect('Типы признаков на модели',families,
                                     default=families,format_func=lambda x:labels[x],
                                     key='pharma_visible_families')
                visible=[v for v in groups if v['Семейство'] in shown]
                if visible:
                    try:
                        components.html(_pharmacophore_html(complex_fn(selected[0]),visible),height=522,scrolling=False)
                    except Exception as exc:
                        st.warning(f'Не удалось показать 3D-признаки: {exc}')
                else:
                    st.info('Выбери хотя бы один тип признаков для просмотра.')
                if len(selected)>1:
                    st.caption('Повторяемость считается по успешно обработанным лигандам; одиночные признаки не выдаются за общие.')
                else:
                    st.caption('Загружен один комплекс: это признаки одной молекулы, а не общий фармакофор.')
                with st.expander('Подробные данные · координаты и экспорт'):
                    df=pd.DataFrame([{k:v for k,v in item.items() if not k.startswith('_')} for item in groups])
                    st.dataframe(df,hide_index=True,use_container_width=True)
                    st.download_button('Скачать признаки · CSV',df.to_csv(index=False).encode('utf-8-sig'),file_name='pharmacophore_features.csv',mime='text/csv')
                    st.download_button('Скачать координаты · JSON',json.dumps(
                        [{k:v for k,v in r.items() if k != '_color'} for r in groups],ensure_ascii=False,indent=2).encode('utf-8'),
                        file_name='pharmacophore_features.json',mime='application/json')
            else:
                st.info('Повторяющихся признаков между разными лигандами пока не найдено. Ниже можно просмотреть отдельные признаки, но это НЕ общий фармакофор.')
                items=[{'Признак':v['name'],'Семейство':v['family'],
                        '_xyz':v['xyz'],'_color':v['color']} for v in points]
                families=sorted({v['Семейство'] for v in items})
                included=st.multiselect('Какие отдельные признаки показать',families,
                    default=families,format_func=lambda x:FEATURE_NAMES.get(x,x),key='individual_visible_families')
                individual=[v for v in items if v['Семейство'] in included]
                if individual:
                    components.html(_pharmacophore_html(complex_fn(selected[0]),individual),
                                    height=522,scrolling=False)
            with st.expander('Технические данные · все исходные координаты'):
                details=pd.DataFrame([{
                    k:v for k,v in item.items() if k != 'xyz'
                }|{'X':round(item['xyz'][0],2),'Y':round(item['xyz'][1],2),'Z':round(item['xyz'][2],2)} for item in points])
                st.dataframe(details,hide_index=True,use_container_width=True)
            st.caption('Это предварительное сопоставление 3D-признаков. Для фармакофорной модели дополнительно проверяются взаимодействия, геометрия, активность и отрицательные примеры.')
    else:
        st.caption('Нажмите «Найти пространственные признаки», чтобы начать расчёт.')

    st.divider()
    st.subheader('Дополнительно · создать молекулярные варианты')
    st.write('**Что такое варианты?** Можно комбинировать фрагменты (BRICS) или создавать ограниченные структурные аналоги исходных молекул. Это только гипотезы: фармакофор, активность, синтезируемость и докинг ещё не проверены.')
    with st.expander('Экспериментальная генерация структур',expanded=False):
        n=st.select_slider('Сколько вариантов показать',options=[5,10,15,20],value=5)
        if st.button('Собрать варианты из фрагментов'):
            st.session_state['candidates_attempted']=True
            st.session_state['generated_candidates']=[]
            st.session_state['candidates_error']=None
            inputs=[]
            for ligand in selected:
                try:
                    value=chem_info_fn(ligand).get('smiles','')
                    if value:inputs.append(value)
                except Exception as exc:
                    st.warning(f"{ligand['PDB']}: не удалось получить SMILES ({exc})")
            if not inputs:
                st.session_state['candidates_error']='Для генерации не удалось получить ни одного SMILES исходных лигандов из PDB.'
            else:
                try:
                    with st.spinner('Собираем химически допустимые комбинации фрагментов…'):
                        rows,num_frags=generate_candidates(inputs,n)
                    st.session_state['generated_candidates']=rows
                    st.session_state['fragments_count']=num_frags
                    st.session_state['candidate_generation_method']='BRICS'
                except Exception as exc:
                    st.session_state['candidates_error']=f'{type(exc).__name__}: {exc}'
        if st.button('Создать аналоги исходного лиганда',key='run_analogs'):
            st.session_state['candidates_attempted']=True
            st.session_state['generated_candidates']=[]
            st.session_state['candidates_error']=None
            st.session_state['candidate_generation_method']='аналоги'
            inputs=[]
            for ligand in selected:
                try:
                    sm=chem_info_fn(ligand).get('smiles','')
                    if sm:inputs.append(sm)
                except Exception as exc:
                    st.warning(f"{ligand['PDB']}: SMILES недоступен ({exc})")
            try:
                rows,n_checked=generate_analogs(inputs,n)
                st.session_state['generated_candidates']=rows
                st.session_state['fragments_count']=n_checked
            except Exception as exc:
                st.session_state['candidates_error']=f'{type(exc).__name__}: {exc}'
        if st.session_state.get('candidates_error'):
            st.error(st.session_state['candidates_error'])
        candidates=st.session_state.get('generated_candidates',[])
        if candidates:
            st.success((f'Получено {len(candidates)} аналогов исходных молекул.' if st.session_state.get('candidate_generation_method')=='аналоги' else f'Получено {len(candidates)} вариантов из {st.session_state.get("fragments_count","—")} фрагментов.'))
            try:
                from rdkit import Chem
                from rdkit.Chem import Draw
                for start in range(0,len(candidates),3):
                    columns=st.columns(3)
                    for index,c in enumerate(candidates[start:start+3],start=start+1):
                        with columns[index-start-1]:
                            st.markdown(f'**Вариант {index}**')
                            mol=Chem.MolFromSmiles(c['SMILES'])
                            if mol is not None:st.image(Draw.MolToImage(mol,size=(350,220)),use_container_width=True)
                            st.caption(f"Масса: {c['Масса (Да)']} Да · cLogP: {c['cLogP']}")
                            st.caption(f"Метод: {c.get('Способ','BRICS')} · фармакофор: не проверен · докинг: не проведён")
                            if st.button('Выбрать для докинга →',key=f'candidate_to_dock_{index}',use_container_width=True):
                                st.session_state['selected_candidate_smiles']=c['SMILES']
                                st.session_state['candidate_sdf']=None
                                st.session_state['candidate_sdf_smiles']=None
                                st.session_state['stage']='docking'
                                st.rerun()
                            with st.popover('SMILES и характеристики'):
                                st.code(c['SMILES'],language=None)
                                st.json({k:v for k,v in c.items() if k!='SMILES'})
            except Exception as exc:
                st.warning(f'Не удалось нарисовать 2D-формулы: {exc}. Химические данные доступны ниже.')
            with st.expander('Все характеристики вариантов · таблица'):
                st.dataframe(pd.DataFrame(candidates),hide_index=True,use_container_width=True)
            st.download_button('Скачать варианты · CSV',pd.DataFrame(candidates).to_csv(index=False).encode('utf-8-sig'),file_name='molecular_variants.csv',mime='text/csv')
            st.caption('«Выбрать для докинга» перенесёт SMILES в следующий раздел. Эти методы пока не используют найденную фармакофорную модель как фильтр.')
        elif st.session_state.get('candidates_attempted',False) and not st.session_state.get('candidates_error'):
            st.info('Варианты не получены этим методом. Можно попробовать вторую кнопку — «Создать аналоги исходного лиганда». Если и она не даст результатов, структура не подходит для наших ограниченных правил перечисления.')
        else:
            st.caption('Генерация ещё не запускалась. Этот этап необязателен и не влияет на поиск фармакофора выше.')

def render_docking_prep(selected,session_cifs,all_ligands=None,*,complex_fn=None,chem_info_fn=None):
    """Simple user path: choose molecule, target structure, binding region, Dock."""
    st.markdown('<div class="lab-kicker">Этап 03 / Докинг</div>',unsafe_allow_html=True)
    st.subheader('Проверить молекулу в выбранной области белка')
    st.caption('Выбери молекулу и белковую структуру. Область связывания выберем по экспериментальному лиганду, если он есть. Подготовку файлов сайт выполнит самостоятельно.')
    if not selected:
        st.info('Сначала загрузи хотя бы один комплекс белок–лиганд в разделе «Исследование комплексов».')
        return

    candidates=st.session_state.get('generated_candidates',[])
    previously=st.session_state.get('selected_candidate_smiles','')
    choices={'Ввести свой SMILES':None}
    for idx,c in enumerate(candidates):
        choices[f"Созданная молекула {idx+1}"]=c['SMILES']
    if previously and previously not in choices.values():
        choices={'Перенесённый кандидат':previously,**choices}
    selected_label=st.selectbox('Молекула для докинга',list(choices),index=0,key='dock_candidate_choice')
    selected_sm=choices[selected_label]
    if selected_sm:
        raw=selected_sm
        st.code(raw,language=None)
    else:
        raw=st.text_input('SMILES молекулы',value=previously,placeholder='Введи SMILES',key='dock_manual_smiles')
    smiles_value=''
    sdf_value=None
    if raw.strip():
        try:
            data=properties(raw.strip())
            smiles_value=data['SMILES']
            with st.expander('Посмотреть свойства и 3D-конформер молекулы',expanded=False):
                st.dataframe(pd.DataFrame([data]),hide_index=True,use_container_width=True)
                if st.button('Построить отдельную 3D-модель',key='dock_preview'):
                    st.session_state['candidate_sdf']=conformer_sdf(smiles_value)
                    st.session_state['candidate_sdf_smiles']=smiles_value
                if st.session_state.get('candidate_sdf_smiles')==smiles_value and st.session_state.get('candidate_sdf'):
                    sdf_value=st.session_state['candidate_sdf']
                    components.html(sdf_viewer_html(st.session_state['candidate_sdf']),height=430,scrolling=False)
                    st.download_button('Скачать SDF',sdf_value,file_name='candidate_3D.sdf',mime='chemical/x-mdl-sdfile')
        except ValueError as exc:
            st.error('Некорректная структура: '+str(exc))

    targets={f"{x['PDB']} · белковая цепь {x.get('_primary_chain','?')}":x for x in selected}
    chosen=st.selectbox('Структура белка для докинга',list(targets),index=0,key='dock_receptor_choice')
    ligand=targets[chosen]
    try:
        protein=protein_only_pdb(ligand['PDB'],session_cifs[ligand['PDB']],ligand['_primary_chain'])
    except Exception as exc:
        st.error('Не удалось прочитать белок из PDB: '+str(exc))
        return
    from docking_ui import render_vina_panel
    def reference_provider(anchor):
        if not complex_fn or not chem_info_fn:
            raise ValueError('Не подключено чтение экспериментального лиганда из PDB.')
        return complex_fn(anchor),chem_info_fn(anchor).get('smiles','')
    render_vina_panel(ligand,session_cifs[ligand['PDB']],protein,
                      sdf_value,smiles_value,available_sites=all_ligands or selected,
                      reference_provider=reference_provider)
    with st.expander('Скачать исходный белок для проверки',expanded=False):
        st.download_button('Белок · PDB',protein,file_name=f"{ligand['PDB']}_protein.pdb",mime='chemical/x-pdb')
