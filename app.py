"""Molecular Design Lab v4.7 — исследование комплексов, ProLIF и предварительный фармакофор.

На страницах отображаются только данные из загруженных пользователем PDB.
Таблица химических взаимодействий строится из успешно проведённых расчётов.
"""
from __future__ import annotations

from collections import defaultdict
from html import escape
import importlib.util
import re

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components

from analysis import (analyze_structure, compare_selected_ligands,
                      download_chemcomp, download_structure, parse_pdb_ids)
from interaction_comparison import compare_prolif, ligand_summary, readable_residue
from optional_interactions import run_prolif
from unified_contacts import unified_contact_rows
from pharma_ui import render_design, render_docking_prep
from overlay import prepare_overlay, overlay_html
from ui_styles import hero, inject_styles, metric_card, render_header
from translation_support import install_browser_translation_guard
from viewer import complex_pdb, ligand_svg, responsive_ligand_html, viewer_html

st.set_page_config(page_title='Molecular Design Lab', page_icon='🧬', layout='wide',
                   initial_sidebar_state='collapsed')
install_browser_translation_guard()
inject_styles()

STAGES = [('research','01 · Исследование комплексов'),
          ('design','02 · Фармакофор и дизайн'),('docking','03 · Докинг')]
DEFAULTS = {'stage':'research','summary':[],'ligands':[],'errors':[],
            'cifs':{},'chain_data':{},'loaded_ids':[],'prolif_results':{},'feature_points':[],'pharma_clusters':[],
            'pharma_errors':[],'generated_candidates':[], 'feature_attempted':False,'candidates_attempted':False, 'candidates_error':None}
for key,value in DEFAULTS.items():
    if key not in st.session_state:
        st.session_state[key]=value


@st.cache_data(ttl=86400, show_spinner=False)
def cached_download(pdb_id):
    return download_structure(pdb_id)


@st.cache_data(ttl=86400, show_spinner=False)
def cached_chemcomp(code):
    return download_chemcomp(code)


@st.cache_data(ttl=86400, show_spinner=False)
def cached_complex(pdb_id, cif_text, uid, code, chain, number, ins):
    return complex_pdb(pdb_id, cif_text, {'Код компонента':code,'Цепь':chain,
                                         'Номер':number,'Вставка':ins})


def pdb_text(lig):
    return cached_complex(lig['PDB'],st.session_state.cifs[lig['PDB']],lig['_uid'],
                          lig['Код компонента'],lig['Цепь'],lig['Номер'],lig.get('Вставка',''))


def chemical_info(lig):
    return cached_chemcomp(lig['Код компонента'])


def chemical_name(lig):
    return chemical_info(lig).get('name') or lig['Код компонента']


def short_name(lig):
    name=chemical_name(lig)
    return name if len(name)<=68 else name[:65].rstrip()+'…'


def public_rows(rows):
    return pd.DataFrame([{k:v for k,v in r.items() if not k.startswith('_')} for r in rows])


def select_ligands():
    """Выбор появляется только если в структуре больше одного кандидата."""
    grouped=defaultdict(list)
    for l in st.session_state.ligands:
        grouped[l['PDB']].append(l)
    picked=[]
    for pdb in st.session_state.loaded_ids:
        available=grouped[pdb]
        preferred=[x for x in available if x['Категория']=='Лиганд-кандидат']
        additional=[x for x in available if x['Категория']!='Лиганд-кандидат']
        if not available:
            st.info(f'{pdb}: небелковые молекулы не найдены.')
            continue
        if preferred and additional:
            with st.expander(f'{pdb} · другие компоненты ({len(additional)})'):
                st.checkbox('Разрешить выбирать кофакторы, ионы и малые компоненты',
                            key=f'extra_{pdb}',value=False)
        choices=(available if st.session_state.get(f'extra_{pdb}',False) or not preferred
                 else preferred)
        if len(choices)==1:
            chosen=choices[0]
        else:
            keys=[l['_uid'] for l in choices]
            uid=st.selectbox(f'{pdb} · молекула для исследования',keys,
                key=f'choose_{pdb}_{bool(st.session_state.get(f"extra_{pdb}",False))}',
                format_func=lambda u:next(
                    f"{short_name(l)} · {l['Код компонента']}, цепь {l['Цепь']}, №{l['Номер']}"
                    for l in choices if l['_uid']==u))
            chosen=next(l for l in choices if l['_uid']==uid)
        picked.append(chosen)
    return picked


def render_complex(lig):
    pdb=lig['PDB']; code=lig['Код компонента']
    fullname=chemical_name(lig)
    displayed=short_name(lig)
    st.markdown(
        f'<div class="lab-complex-kicker">КОМПЛЕКС {escape(pdb)} · {escape(code)}</div>'
        f'<div class="lab-compound-name">{escape(displayed)}</div>'
        f'<div class="lab-intro">Лиганд: {escape(code)} · цепь {escape(lig["Цепь"])} · '
        f'номер в структуре {escape(str(lig["Номер"]))}</div>',
        unsafe_allow_html=True)
    if fullname != displayed:
        with st.expander('Полное название вещества'):
            st.write(fullname)
    visual,formula=st.tabs(['3D · белок и лиганд','2D · структурная формула'])
    with visual:
        try:
            components.html(viewer_html(pdb_text(lig),lig,height=585,mode='full'),
                            height=585,scrolling=False)
        except Exception as exc:
            st.warning(f'Не удалось открыть 3D: {exc}')
    with formula:
        info=chemical_info(lig)
        svg=ligand_svg(info.get('smiles',''))
        if svg:
            components.html(responsive_ligand_html(svg),height=380,scrolling=False)
        else:
            st.info('В каталоге химических компонентов нет доступной 2D-формулы.')
        if info.get('formula'):
            st.caption('Молекулярная формула: '+info['formula'])
    c1,c2=st.columns(2)
    with c1:
        st.link_button('Открыть комплекс в RCSB PDB',
                       f'https://www.rcsb.org/structure/{pdb}',use_container_width=True)
    with c2:
        st.link_button('Открыть сведения о молекуле',
                       f'https://www.rcsb.org/ligand/{code}',use_container_width=True)
    # Geometrical contacts are now displayed WITH ProLIF outcomes in a single
    # understandable residue table, after the interaction calculation widget.


def render_research():
    """Единый исследовательский этап: один комплекс или сравнение одной мишени."""
    st.markdown('<div class="lab-kicker">Этап 01 / Исследование комплексов</div>', unsafe_allow_html=True)
    st.subheader('Белок и взаимодействующие с ним молекулы')
    st.caption('Загрузите один комплекс, чтобы изучить лиганд. Несколько комплексов той же мишени можно сравнить между собой.')
    with st.form('pdb_entry'):
        raw = st.text_input('Идентификаторы PDB', placeholder='Один или несколько PDB ID через запятую')
        go = st.form_submit_button('Загрузить структуры', type='primary')
    if go:
        try:
            ids = parse_pdb_ids(raw)
        except ValueError as exc:
            st.error(str(exc)); ids = []
        if ids:
            summaries=[]; ligands=[]; errors=[]; cifs={}; chains={}
            for code in ids:
                try:
                    with st.spinner(f'Загрузка {code}…'):
                        cif=cached_download(code)
                        found,meta=analyze_structure(code,cif)
                    cifs[code]=cif
                    chains[code]=meta.pop('_chains')
                    summaries.append({'PDB':code,**meta,'Небелковых молекул':len(found)})
                    ligands.extend(found)
                except Exception as exc:
                    errors.append(f'{code}: {exc}')
            st.session_state.update(
                summary=summaries,ligands=ligands,errors=errors,cifs=cifs,
                chain_data=chains,loaded_ids=list(cifs),prolif_results={},
                feature_points=[],pharma_clusters=[],pharma_errors=[],
                generated_candidates=[],feature_attempted=False,
                candidates_attempted=False,candidates_error=None,
                candidate_sdf=None,candidate_sdf_smiles=None,vina_result=None,
                vina_signature=None,receptor_repair_report=None,
                repaired_receptor_pdb=None)
    for error in st.session_state.errors:
        st.warning(error)
    if not st.session_state.summary:
        st.caption('Сайт не содержит заранее выбранных белков: результаты появятся после загрузки ваших PDB.')
        return
    cols=st.columns(3)
    with cols[0]:metric_card('Комплексов',len(st.session_state.summary))
    with cols[1]:metric_card('Небелковых молекул',len(st.session_state.ligands))
    with cols[2]:metric_card('Цепей белка',sum(m['Белковых цепей'] for m in st.session_state.summary))
    with st.expander('Состав загруженных структур · технические данные'):
        st.dataframe(pd.DataFrame(st.session_state.summary),hide_index=True,use_container_width=True)
    selected=select_ligands()
    if not selected:
        st.info('Для анализа нужно выбрать хотя бы одну небелковую молекулу.')
        return
    if len(selected)>1:
        mode=st.radio('Режим исследования', ['Изучить один комплекс','Сравнить комплексы'],
                      horizontal=True, key='research_mode')
    else:
        mode='Изучить один комплекс'
    if mode=='Изучить один комплекс':
        if len(selected)>1:
            options={lig['_uid']:lig for lig in selected}
            choice=st.selectbox('Комплекс для просмотра',list(options),
                format_func=lambda key:f"{options[key]['PDB']} · {short_name(options[key])} ({options[key]['Код компонента']})",
                key='research_active_complex')
            active=options[choice]
        else:
            active=selected[0]
        render_complex(active)
        st.markdown('<div class="lab-section">Химические взаимодействия</div>',unsafe_allow_html=True)
        st.caption('Для одного лиганда ProLIF покажет, с какими остатками он взаимодействует и каким образом.')
        has_prolif=importlib.util.find_spec('prolif') is not None
        if st.button('Рассчитать взаимодействия этого лиганда',key='one_prolif',
                     disabled=not has_prolif,type='primary'):
            compute_all([active])
        if not has_prolif:
            st.info('Для химического анализа установите ProLIF в окружение сайта.')
        result=get_result(active['_uid'])
        if result is not None:
            if result[1]:st.error(f'ProLIF: {result[1]}')
            else:
                st.success(f'ProLIF обнаружил {len(result[0])} записей взаимодействий.')
                show_one_records(active)
        show_unified_contacts(active)
        if len(selected)>1:
            st.caption('Чтобы увидеть общие взаимодействия нескольких лигандов, выберите режим «Сравнить комплексы» выше.')
    else:
        render_analysis(selected)
    with st.expander('Все небелковые компоненты · техническая таблица'):
        st.dataframe(public_rows(st.session_state.ligands),hide_index=True,use_container_width=True)
    st.divider()
    if st.button('Перейти к фармакофору и проектированию →',type='primary',key='research_next'):
        st.session_state.stage='design'
        st.rerun()


def get_result(uid):
    entry=st.session_state.prolif_results.get(uid)
    if entry is None:return None
    return entry if isinstance(entry,tuple) else (entry.get('records',[]),entry.get('error'))


def compute_all(selected):
    if importlib.util.find_spec('prolif') is None:
        st.error('ProLIF отсутствует в том Python, которым запущен сайт.')
        st.code('.\\.venv\\Scripts\\python.exe -m pip install "prolif>=2.2.2,<3"',language='powershell')
        return
    progress=st.progress(0,text='Подготовка расчётов')
    for i,lig in enumerate(selected):
        uid=lig['_uid']
        try:
            with st.spinner(f"ProLIF: {lig['PDB']} · {lig['Код компонента']}…"):
                records=run_prolif(pdb_text(lig),lig,chemical_info(lig).get('smiles',''))
            st.session_state.prolif_results[uid]=(records,None)
        except Exception as exc:
            st.session_state.prolif_results[uid]=([],str(exc))
        progress.progress((i+1)/len(selected),text=f'Рассчитано {i+1} из {len(selected)}')
    progress.empty()


def show_unified_contacts(lig):
    """One table merges close atom distances and predicted interaction kinds."""
    result=get_result(lig['_uid'])
    records=None if result is None or result[1] else result[0]
    rows=unified_contact_rows(lig['_contact_rows'],records,chain=lig.get('_primary_chain'))
    st.markdown('**Контакты лиганда с белком**')
    st.caption('Расстояние до 4 Å показывает геометрическую близость. ProLIF отдельно оценивает тип взаимодействия и может найти контакты за пределами 4 Å.')
    if result and result[1]:
        for row in rows:
            row['Контакт по ProLIF']='Ошибка расчёта'
        st.warning('ProLIF завершился с ошибкой. В таблице показаны только геометрические данные.')
    if rows:
        df=pd.DataFrame(rows)
        for col in ('Ближайшие атомы, Å','Расстояние ProLIF, Å'):
            df[col]=df[col].map(lambda v: f'{v:.2f}' if pd.notna(v) else '—')
        st.dataframe(df,hide_index=True,use_container_width=True)
        st.download_button('Скачать общую таблицу · CSV',
            df.to_csv(index=False).encode('utf-8-sig'),
            file_name=f"{lig['PDB']}_{lig['Код компонента']}_unified_contacts.csv",
            mime='text/csv',key=f"unified_contacts_{lig['_uid']}")
    else:
        st.info('Близких остатков нет; если ProLIF ещё не запускался, рассчитай взаимодействия выше.')
    with st.expander('Технические подробности · атомные пары и исходные записи ProLIF'):
        st.dataframe(public_rows(lig['_contact_rows']),hide_index=True,use_container_width=True)
        if records is not None:
            st.markdown('**Исходные записи ProLIF**')
            st.dataframe(pd.DataFrame(records),hide_index=True,use_container_width=True)


def show_one_records(lig):
    result=get_result(lig['_uid'])
    if result is None:return
    records,error=result
    if error:
        st.error(f"{lig['PDB']} · {lig['Код компонента']}: {error}")
        return
    with st.expander(f"{lig['PDB']} · {lig['Код компонента']} — все взаимодействия ({len(records)})"):
        if records:
            df=pd.DataFrame(records)
            st.dataframe(df,hide_index=True,use_container_width=True)
            st.download_button('Скачать CSV',df.to_csv(index=False).encode('utf-8-sig'),
                file_name=f"prolif_{lig['PDB']}_{lig['Код компонента']}.csv",
                mime='text/csv',key='prolif_csv_'+lig['_uid'])
        else:
            st.info('В выбранной подготовке взаимодействия не распознаны.')


def summary_table(selected):
    rows=[]
    for lig in selected:
        res=get_result(lig['_uid'])
        item={'Комплекс':f"{lig['PDB']} · {lig['Код компонента']}"}
        if res is None:
            item.update({'Статус':'Не рассчитан','H-связи':None,'Гидрофобные':None,
                         'π-взаимодействия':None,'Ионные':None,'Ван-дер-Ваальса':None})
        elif res[1]:
            item.update({'Статус':'Ошибка расчёта','H-связи':None,'Гидрофобные':None,
                         'π-взаимодействия':None,'Ионные':None,'Ван-дер-Ваальса':None})
        else:
            counts=ligand_summary(res[0])
            item.update({'Статус':'Рассчитан','H-связи':counts.get('H-связь',0),
                         'Гидрофобные':counts.get('Гидрофобный',0),
                         'π-взаимодействия':counts.get('π-взаимодействие',0),
                         'Ионные':counts.get('Ионное',0),
                         'Ван-дер-Ваальса':counts.get('Ван-дер-Ваальса',0)})
        rows.append(item)
    return pd.DataFrame(rows)


def show_alignment_and_geometry(selected):
    with st.expander('Технические данные · выравнивание и расстояния'):
        if len(selected)<2:
            st.caption('Для выравнивания нужны две структуры.')
            return
        matrix,report,issues=compare_selected_ligands(
            selected,st.session_state.chain_data,st.session_state.cifs)
        for issue in issues:st.warning(issue)
        st.caption('RMSD Cα описывает различие координат сопоставленных белков после наложения, не силу ингибирования.')
        if report:st.dataframe(pd.DataFrame(report),hide_index=True,use_container_width=True)
        if matrix:
            df=pd.DataFrame(matrix)
            st.dataframe(df,hide_index=True,use_container_width=True)
            st.download_button('Скачать геометрическую матрицу',
                df.to_csv(index=False).encode('utf-8-sig'),
                file_name='geometry_contacts.csv',mime='text/csv')


def overlay_ref_highlights(included,shared,selected_label=None):
    if not included:return []
    reference=included[0]
    names={r['Остаток'] for r in shared}
    if selected_label:names.add(selected_label)
    seq=st.session_state.chain_data.get(reference['PDB'],{}).get(reference['_primary_chain'],[])
    return [{'chain':reference['_primary_chain'],'resi':r['number']}
            for r in seq if readable_residue(r['key']) in names]


def render_analysis(selected):
    st.subheader('Что общего у взаимодействий разных лигандов?')
    if not selected:return
    st.caption('Типы химических взаимодействий определяются ProLIF. Структуры должны относиться к одной или сопоставимой белковой мишени.')
    has_prolif=importlib.util.find_spec('prolif') is not None
    if not has_prolif:
        st.warning('ProLIF не найден в текущем окружении. Установите библиотеку и перезапустите Streamlit.')
        with st.expander('Как установить ProLIF в локальное окружение'):
            st.code('.\\.venv\\Scripts\\python.exe -m pip install "prolif>=2.2.2,<3"',language='powershell')
            st.code('.\\.venv\\Scripts\\python.exe -m streamlit run app.py',language='powershell')
    if st.button('Рассчитать взаимодействия для всех выбранных молекул',type='primary',
                 disabled=not has_prolif,use_container_width=True):
        compute_all(selected)
        st.rerun()

    if any(get_result(l['_uid']) is not None for l in selected):
        st.markdown('<div class="lab-section">Типы взаимодействий</div>',unsafe_allow_html=True)
        st.dataframe(summary_table(selected),hide_index=True,use_container_width=True)
        st.caption('Количество — это число различных аминокислотных остатков, для которых ProLIF обнаружил данный тип взаимодействия, а не число связей или оценка силы ингибирования.')
        for lig in selected:
            res=get_result(lig['_uid'])
            if res and res[1]:
                st.error(f"{lig['PDB']} · {lig['Код компонента']}: {res[1]}")
        st.caption('Подробные записи ProLIF доступны отдельно для каждого комплекса:')
        for lig in selected:
            show_one_records(lig)
    else:
        st.info('Нажмите «Рассчитать взаимодействия» — после расчёта здесь появятся данные, полученные из загруженных структур.')

    successful={lig['_uid']:get_result(lig['_uid'])[0]
                for lig in selected if get_result(lig['_uid']) is not None
                and get_result(lig['_uid'])[1] is None}
    matrix,shared,issues,included=compare_prolif(selected,successful,st.session_state.chain_data)
    for issue in issues:
        st.warning(issue)
    if len(included)>=2:
        st.markdown('<div class="lab-section">Сравнение по аминокислотам</div>',unsafe_allow_html=True)
        st.caption('В строке — соответствующая позиция белка. В столбцах — взаимодействия с каждым рассчитанным лигандом. Прочерк означает, что ProLIF не обнаружил взаимодействие этого типа.')
        interaction_df=pd.DataFrame(matrix)
        choice=st.dataframe(interaction_df,hide_index=True,use_container_width=True,
                            height=min(500,90+35*len(matrix)),
                            on_select='rerun',selection_mode='single-row',key='matrix_selection')
        chosen_residue=None
        try:
            row_ids=choice.selection.rows
            if row_ids:
                chosen_residue=interaction_df.iloc[row_ids[0]]['Остаток']
        except (AttributeError,IndexError,KeyError,TypeError):
            pass
        if chosen_residue:
            st.caption(f'Выбрана позиция {chosen_residue}. Она будет дополнительно выделена в 3D-наложении ниже.')
        st.download_button('Скачать сравнительную таблицу CSV',
            interaction_df.to_csv(index=False).encode('utf-8-sig'),
            file_name='prolif_comparison.csv',mime='text/csv')
        st.markdown('<div class="lab-section">Что действительно совпадает</div>',unsafe_allow_html=True)
        if shared:
            st.dataframe(pd.DataFrame(shared),hide_index=True,use_container_width=True)
            st.caption('Здесь указаны только те типы взаимодействий, которые ProLIF обнаружил с соответствующим остатком у каждого включённого в сравнение лиганда.')
        else:
            st.info('Одинаковых типов взаимодействий с одними и теми же позициями пока не найдено.')
        if len(included)!=len(selected):
            st.warning('Сравнение неполное: часть выбранных комплексов ещё не рассчитана или не сопоставляется с опорным белком.')
    elif len(selected)>1:
        st.caption('Сравнительная таблица появится, когда хотя бы две молекулы будут успешно рассчитаны и белковые цепи сопоставлены.')
    elif included:
        st.caption('Вы загрузили один комплекс. Его взаимодействия доступны выше; для поиска общих признаков понадобятся дополнительные комплексы той же мишени.')

    if len(included)>=2:
        with st.expander('3D · сравнить положения лигандов'):
            st.caption('Здесь показаны несколько разных лигандов после совмещения соответствующих белковых цепей. Каждый лиганд имеет свой цвет.')
            try:
                base,layers,problems=prepare_overlay(included,st.session_state.cifs,st.session_state.chain_data)
                for problem in problems:st.warning(problem)
                if len(layers)>=2:
                    components.html(overlay_html(base,layers,overlay_ref_highlights(included,shared,chosen_residue)),
                                    height=590,scrolling=False)
                else:
                    st.info('Для наложения не удалось подготовить две совместимые структуры.')
            except Exception as exc:
                st.warning(f'Не удалось построить совмещение: {exc}')
    show_alignment_and_geometry(selected)
    st.caption('Расчёты ProLIF — компьютерная интерпретация структуры, а не экспериментальное доказательство активности.')


def render_future(stage):
    selected=select_ligands() if st.session_state.summary else []
    chosen_key=tuple(l['_uid'] for l in selected)
    if st.session_state.get('_design_source')!=chosen_key:
        st.session_state.update(_design_source=chosen_key,feature_points=[],pharma_clusters=[],
                                pharma_errors=[],generated_candidates=[],feature_attempted=False,candidates_attempted=False,candidates_error=None,candidate_sdf=None,
                                candidate_sdf_smiles=None)
    if stage=='design':
        render_design(selected,pdb_text,chemical_info,st.session_state.cifs,st.session_state.chain_data)
    else:
        render_docking_prep(selected,st.session_state.cifs,st.session_state.ligands,
                            complex_fn=pdb_text,chem_info_fn=chemical_info)


def navigation():
    st.markdown('<div class="lab-section">Ход исследования</div>',unsafe_allow_html=True)
    cols=st.columns(3,gap='small')
    for col,(stage,label) in zip(cols,STAGES):
        with col:
            st.button(label,key='nav_'+stage,
                      type='primary' if st.session_state.stage==stage else 'secondary',
                      on_click=lambda p=stage:st.session_state.update(stage=p),
                      use_container_width=True)


render_header()
hero('От белковой структуры к новой молекуле','Molecular Design Lab',
     'Исследование комплексов → фармакофор и дизайн → молекулярный докинг.')
navigation()
st.divider()
if st.session_state.stage=='research':
    render_research()
else:
    render_future(st.session_state.stage)
