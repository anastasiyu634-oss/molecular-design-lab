"""Minimal one-click exploratory docking screen.

The user chooses a scientifically grounded binding region, not just a chain.
SMILES -> SDF, receptor preparation, Meeko, and Vina are internal steps.
Failures never generate invented affinity scores.
"""
from __future__ import annotations

from html import escape
import importlib.util
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from binding_sites import observed_sites, box_from_residues
from docking_runner import binding_site_box, find_vina_executable
from automatic_docking import run_automatic_docking, AutomaticDockingError
from viewer import docked_pose_viewer_html, sdf_viewer_html


def render_vina_panel(ligand: dict, cif_text: str, protein_pdb: str,
                      sdf_text: str | None, candidate_smiles: str,
                      available_sites: list[dict] | None = None,
                      reference_provider=None):
    """Render one-button flow; optional expert controls cannot obstruct default."""
    st.subheader('Докинг')
    if not candidate_smiles:
        st.info('Для проверки собственной молекулы введи SMILES. Контрольный редокинг известного лиганда доступен ниже и без собственного SMILES.')
    has_meeko=importlib.util.find_spec('meeko') is not None
    has_vina=(find_vina_executable() is not None or importlib.util.find_spec('vina') is not None)
    if not (has_meeko and has_vina):
        missing=[]
        if not has_meeko: missing.append('Meeko')
        if not has_vina: missing.append('AutoDock Vina')
        st.error('Для докинга пока недоступно: '+', '.join(missing)+'.')
        with st.expander('Как подключить недостающие компоненты'):
            if not has_meeko:
                st.code(r'.\.venv\Scripts\python.exe -m pip install -r requirements-docking.txt',language='powershell')
            if not has_vina:
                st.write('Положи официальный Windows-файл AutoDock Vina под именем `bin/vina.exe`.')
                st.link_button('Официальные выпуски Vina','https://github.com/forlilab/AutoDock-Vina/releases')
        return

    st.caption('Сайт самостоятельно построит 3D-лиганд, подготовит белок и выполнит расчёт Vina. Это компьютерная оценка, не экспериментальное подтверждение активности.')
    with st.expander('Посмотреть ранее сохранённые позы · PDBQT',expanded=False):
        saved=st.file_uploader('Загрузить vina_poses.pdbqt',type=['pdbqt'],key='saved_vina_poses_495')
        if saved is not None:
            raw=saved.getvalue()
            if len(raw)>5_000_000:
                st.error('Файл слишком большой для интерактивного просмотра (максимум 5 МБ).')
            else:
                from docking_pose_format import pdbqt_poses_to_sdf, PoseFormatError
                poses_text=raw.decode('utf-8',errors='replace')
                try:
                    preview=pdbqt_poses_to_sdf(poses_text)
                    st.caption(f"Найдено {len(preview)} поз. Визуализация не запускает новый докинг. Используйте исходную структуру белка из того же проекта.")
                    components.html(docked_pose_viewer_html(protein_pdb,poses_text),height=535,scrolling=False)
                    st.download_button('Скачать все позы · SDF',
                                       ''.join(pose['sdf'] for pose in preview),
                                       file_name='vina_poses_complete.sdf',key='saved_vina_sdf_495')
                except PoseFormatError as exc:
                    st.error('Не удалось восстановить полную структуру лиганда: '+str(exc))

    sites=observed_sites(ligand['PDB'],available_sites or [ligand],ligand.get('_primary_chain',''))
    selected_site=None
    site_key=''
    site_option=None
    if sites:
        # The default can be inferred ONLY if a reference ligand is present.
        values=[x['key'] for x in sites]+['__residues__']
        label_lookup={x['key']:x['label'] for x in sites}
        label_lookup['__residues__']='Указать другую область по аминокислотам'
        preferred=next((x['key'] for x in sites if x['ligand'].get('_uid')==ligand.get('_uid')),sites[0]['key'])
        idx=values.index(preferred)
        site_option=st.selectbox('Где искать связывание?',values,index=idx,
                                 format_func=lambda v:label_lookup[v],
                                 key=f"binding_site_{ligand['PDB']}_{ligand.get('_uid','')}")
        selected_site=next((s for s in sites if s['key']==site_option),None)
    else:
        st.info('В этой структуре нет подходящего экспериментального лиганда. Для докинга укажи область белка по аминокислотам.')
        site_option='__residues__'
    residues=''
    if site_option=='__residues__':
        st.caption('Выбери участок, например `GLU353, HIS524`. Это не поиск кармана по всей поверхности, а центрирование области около указанных остатков.')
        residues=st.text_input('Аминокислоты в нужной области',value='',placeholder='GLU353, HIS524',
                               key=f"dock_residues_{ligand['PDB']}")
    try:
        if selected_site:
            anchor=selected_site['ligand']
            # the ligand coordinates must come from the SAME PDB model
            proposed=binding_site_box(anchor['PDB'],cif_text,anchor)
            site_key=selected_site['key']
            chain=selected_site['chain']
        elif residues.strip():
            chain=ligand.get('_primary_chain') or ligand.get('Цепь','A')
            proposed=box_from_residues(protein_pdb,residues,chain)
            site_key='residues:'+residues.strip()+':'+chain
        else:
            proposed=None
    except Exception as exc:
        proposed=None
        st.error('Не удалось определить выбранную область: '+str(exc))

    # Specialist settings are truly optional, with sensible defaults.
    center=(proposed or {}).get('center',[0.0,0.0,0.0]); size=(proposed or {}).get('size',[22.0]*3)
    exhaust=8; num_poses=5; receptor_override=None
    with st.expander('Дополнительно · настройки для исследователей',expanded=False):
        st.caption('Не обязательно менять для первого расчёта. Координаты относятся к модели выбранного белка.')
        if proposed:
            cols=st.columns(3)
            center=[cols[i].number_input(f'Центр {axis}, Å',value=float(center[i]),format='%.2f',step=0.5,
                        key=f"oneclick_center_{ligand['PDB']}_{site_key}_{i}") for i,axis in enumerate('XYZ')]
            cols=st.columns(3)
            size=[cols[i].number_input(f'Размер {axis}, Å',value=float(size[i]),min_value=10.0,max_value=40.0,
                        format='%.2f',step=1.0,key=f"oneclick_size_{ligand['PDB']}_{site_key}_{i}") for i,axis in enumerate('XYZ')]
        exhaust=st.slider('Тщательность поиска',min_value=4,max_value=24,value=8,step=4,key='oneclick_exhaust')
        num_poses=st.slider('Число сохраняемых поз',min_value=1,max_value=10,value=5,key='oneclick_poses')
        upload=st.file_uploader('Если есть проверенный рецептор PDBQT, можно загрузить его',type=['pdbqt'],
                                key=f"oneclick_verified_{ligand['PDB']}")
        if upload is not None:
            receptor_override=upload.getvalue().decode('utf-8',errors='replace')
            if not any(x.startswith(('ATOM ','HETATM')) for x in receptor_override.splitlines()):
                st.error('Загруженный PDBQT не содержит строк с атомами.')
                receptor_override=None
        st.caption('Автоматическая подготовка пробует Meeko и ProDy, при необходимости PDBFixer. Неизвестные или химически некорректные остатки не удаляются молча.')

    if st.button('Запустить докинг',key='oneclick_run',type='primary',disabled=proposed is None or not bool(candidate_smiles),use_container_width=False):
        st.session_state['oneclick_result']=None
        st.session_state['oneclick_error']=None
        st.session_state['oneclick_signature']=None
        st.session_state['pose_prolif_result']=None
        st.session_state['pose_prolif_signature']=None
        st.session_state['all_pose_prolif_signature']=None
        try:
            with st.spinner('Готовим молекулу и рецептор, затем выполняем докинг…'):
                result=run_automatic_docking(candidate_smiles,protein_pdb,center,size,
                                             exhaustiveness=exhaust,poses=num_poses,
                                             receptor_pdbqt=receptor_override,sdf_text=sdf_text)
            st.session_state['oneclick_result']=result
            st.session_state['oneclick_signature']=(ligand.get('_uid'),candidate_smiles,site_key)
        except AutomaticDockingError as exc:
            st.session_state['oneclick_error']={'summary':str(exc),'attempts':exc.attempts}
        except Exception as exc:
            st.session_state['oneclick_error']={'summary':f'{type(exc).__name__}: {exc}','attempts':[]}

    error=st.session_state.get('oneclick_error')
    if error:
        st.error('Докинг не запущен. '+error['summary'])
        st.caption('Исходная модель белка сохранена. Результаты Vina не выдумываются; можно попробовать другую PDB-структуру той же мишени.')
        if error['attempts']:
            with st.expander('Технические подробности подготовки белка'):
                for attempt in error['attempts']:
                    st.write(f"**{attempt.get('method','Подготовка')}** — {attempt.get('message','ошибка')}")
                    if attempt.get('residues'):
                        st.caption('Проблемные остатки: '+', '.join(attempt['residues']))
                    if attempt.get('details'):
                        st.code(str(attempt['details'])[-10000:],language='text')

    # Control redocking is a separate validation calculation using the
    # crystallographically observed ligand and its 3D coordinates. It is not
    # a claim that a user's arbitrary candidate is the reference ligand.
    with st.expander('Проверка метода · контрольный редокинг',expanded=False):
        st.write('Программа извлекает экспериментальный лиганд выбранной области, докирует ЕГО обратно и сравнивает положения атомов без дополнительного совмещения. Это проверка протокола, а не доказательство активности.')
        if not selected_site or not reference_provider:
            st.info('Для редокинга выбери область вокруг экспериментального лиганда из той же PDB-структуры.')
        else:
            anchor=selected_site['ligand']
            control_sig=(ligand['PDB'],anchor.get('_uid'),site_key,tuple(center),tuple(size))
            if st.button('Запустить контрольный редокинг',key='start_redocking',disabled=proposed is None):
                st.session_state['redock_result']=None
                st.session_state['redock_error']=None
                st.session_state['redock_signature']=control_sig
                try:
                    from redocking import experimental_ligand, run_control_redocking
                    ref_pdb,ref_smi=reference_provider(anchor)
                    with st.spinner('Извлекаем кристаллический лиганд, выполняем редокинг и считаем RMSD…'):
                        experimental=experimental_ligand(ref_pdb,anchor,ref_smi)
                        validation=run_control_redocking(experimental,protein,center,size,
                                  exhaustiveness=exhaust,poses=num_poses,receptor_pdbqt=receptor_override)
                    st.session_state['redock_result']=validation
                except Exception as exc:
                    st.session_state['redock_error']=f'{type(exc).__name__}: {exc}'
            if st.session_state.get('redock_signature')==control_sig:
                if st.session_state.get('redock_error'):
                    st.error('Редокинг не выполнен: '+st.session_state['redock_error'])
                validation=st.session_state.get('redock_result')
                if validation:
                    st.success('Контрольный расчёт выполнен; RMSD рассчитан без совмещения лигандов.')
                    a,b=st.columns(2)
                    a.metric('RMSD позы с лучшей оценкой Vina',f"{validation['first_pose_rmsd']:.2f} Å")
                    b.metric('Минимальный RMSD среди поз',f"{validation['best_rmsd']:.2f} Å")
                    st.caption('Это РАЗНЫЕ показатели: поза с лучшей оценкой энергии не обязана быть ближе всех к экспериментальной. Значение ~2 Å часто используют как ориентир, но оно не заменяет научную валидацию.')
                    quality=pd.DataFrame([{'Поза':i+1,'Vina (ккал/моль)':score,'RMSD (Å)':rmsd} for i,(score,rmsd) in enumerate(zip(validation['affinities'],validation['reference_rmsds']))])
                    st.dataframe(quality,hide_index=True,use_container_width=True)
                    st.download_button('Скачать таблицу редокинга · CSV',quality.to_csv(index=False).encode('utf-8-sig'),file_name='redocking_validation.csv')
                    st.caption('Если атомы рецептора были достроены, RMSD относится к моделированному рецептору. Для строгой проверки требуется оценить качество подготовки и сравнить другие настройки докинга.')

    result=st.session_state.get('oneclick_result')
    if not result or st.session_state.get('oneclick_signature')!=(ligand.get('_uid'),candidate_smiles,site_key):
        return
    st.success('Расчёт завершён — доступны рассчитанные положения молекулы в белке.')
    route=result.get('receptor_route') or {}
    if route.get('reconstruction'):
        st.warning('Часть атомов белка достроена PDBFixer. Исходные тяжёлые атомы были зафиксированы во время оптимизации; новые атомы остаются модельными. Необходимы проверка геометрии и контрольный редокинг.')
    if route.get('altloc_report', {}).get('altloc_choices'):
        st.caption('Для альтернативных атомных положений выбраны согласованные варианты по заселённости (occupancy). Подробности доступны в технических данных.')
    st.metric('Лучшая оценка Vina',f"{result['affinities'][0]:.2f} ккал/моль")
    st.caption('Число относится к скоринговой функции Vina, а не измеренной энергии связывания. Оно не доказывает ингибирующую активность.')
    frame=pd.DataFrame({'Поза':list(range(1,len(result['affinities'])+1)),
                        'Оценка Vina (ккал/моль)':result['affinities']})
    st.dataframe(frame,hide_index=True,use_container_width=True)
    with st.expander('Контакты рассчитанных поз · ProLIF',expanded=False):
        st.caption('ProLIF анализирует выбранную рассчитанную позу в подготовленной модели белка. Это не измеренная активность.')
        pose_id=st.selectbox('Выбери позу',range(1,len(result['affinities'])+1),key='prolif_pose_id')
        pose_sig=(st.session_state.get('oneclick_signature'),pose_id)
        if st.button('Рассчитать взаимодействия этой позы',key='pose_prolif_calc'):
            st.session_state['pose_prolif_signature']=pose_sig
            st.session_state['pose_prolif_result']=None
            st.session_state['pose_prolif_error']=None
            try:
                from docked_interactions import pose_contact_rows
                prepared=(result.get('receptor_route') or {}).get('repaired_pdb') or protein_pdb
                with st.spinner('Определяем типы взаимодействий по координатам рассчитанной позы…'):
                    rows=pose_contact_rows(result['poses_pdbqt'],prepared,pose_id)
                st.session_state['pose_prolif_result']=rows
            except Exception as exc:
                st.session_state['pose_prolif_error']=f'{type(exc).__name__}: {exc}'
        if st.button('Рассчитать ProLIF для всех рассчитанных поз',key='all_prolif_calc'):
            st.session_state['all_pose_prolif_signature']=st.session_state.get('oneclick_signature')
            st.session_state['all_pose_prolif_rows']=[]
            st.session_state['all_pose_prolif_errors']=[]
            try:
                from docked_interactions import pose_contact_rows
                prepared=(result.get('receptor_route') or {}).get('repaired_pdb') or protein_pdb
                with st.spinner('Определяем взаимодействия для каждой позы…'):
                    for i in range(1,len(result['affinities'])+1):
                        try:
                            for r in pose_contact_rows(result['poses_pdbqt'],prepared,i):
                                st.session_state['all_pose_prolif_rows'].append({'Поза':i,**r})
                        except Exception as exc:
                            st.session_state['all_pose_prolif_errors'].append(f'Поза {i}: {type(exc).__name__}: {exc}')
            except Exception as exc:
                st.session_state['all_pose_prolif_errors'].append(f'{type(exc).__name__}: {exc}')
        if st.session_state.get('all_pose_prolif_signature')==st.session_state.get('oneclick_signature'):
            all_rows=st.session_state.get('all_pose_prolif_rows') or []
            all_errors=st.session_state.get('all_pose_prolif_errors') or []
            if all_rows:
                df=pd.DataFrame(all_rows)
                st.dataframe(df,hide_index=True,use_container_width=True)
                st.download_button('Скачать все взаимодействия · CSV',df.to_csv(index=False).encode('utf-8-sig'),file_name='all_docked_prolif.csv')
            for message in all_errors:
                st.warning('Некоторые позы не удалось проанализировать: '+message)
            if not all_rows and not all_errors:st.info('В рассчитанных позах не обнаружены взаимодействия выбранных типов.')
        if st.session_state.get('pose_prolif_signature')==pose_sig:
            if st.session_state.get('pose_prolif_error'):
                st.warning('Этот анализ пока недоступен: '+st.session_state['pose_prolif_error'])
            else:
                rows=st.session_state.get('pose_prolif_result')
                if rows is not None:
                    st.write(f'Получено записей взаимодействий: {len(rows)}')
                    if rows:
                        contacts=pd.DataFrame(rows)
                        st.dataframe(contacts,hide_index=True,use_container_width=True)
                        st.download_button('Скачать контакты · CSV',contacts.to_csv(index=False).encode('utf-8-sig'),file_name=f'prolif_docked_pose_{pose_id}.csv')
                    else:st.info('ProLIF не обнаружил контактов среди настроенных типов для этой позы.')

    st.markdown('**Предсказанное положение молекулы в белке**')
    display_pdb=route.get('repaired_pdb') or protein_pdb
    components.html(docked_pose_viewer_html(display_pdb,result['poses_pdbqt']),height=500,scrolling=False)
    with st.expander('Скачать результаты и технические данные'):
        st.write('Рецептор подготовлен способом: '+str(route.get('route','неизвестно')))
        if route.get('altloc_report', {}).get('altloc_choices'):
            st.write('Выбор альтернативных положений: '+', '.join(route['altloc_report']['altloc_choices']))
        report = route.get('repair_report') or {}
        if report.get('restrained_minimization') or report.get('fixed_experimental_heavy_atoms'):
            st.write(f"Достроено тяжёлых атомов: {report.get('added_heavy_atoms', 0)}; "
                     f"RMSD исходных атомов после оптимизации: {report.get('original_heavy_rmsd_a', 0):.3f} Å; "
                     f"максимальное смещение: {report.get('original_heavy_max_shift_a', 0):.3f} Å.")
        st.download_button('Скачать рассчитанные позы · PDBQT',result['poses_pdbqt'],file_name='vina_poses.pdbqt')
        try:
            from docking_pose_format import pdbqt_poses_to_sdf
            sdf_poses=pdbqt_poses_to_sdf(result['poses_pdbqt'])
            st.download_button('Скачать все позы с химическими связями · SDF',
                               ''.join(p['sdf'] for p in sdf_poses),
                               file_name='vina_poses_complete.sdf')
        except ValueError as exc:
            st.caption('SDF с позами пока недоступен: '+str(exc))
        st.download_button('Скачать рецептор · PDBQT',result['receptor_pdbqt'],file_name='receptor.pdbqt')
        st.download_button('Скачать таблицу оценок · CSV',frame.to_csv(index=False).encode('utf-8-sig'),file_name='vina_scores.csv')
        if route.get('repaired_pdb'):
            st.download_button('Скачать восстановленный белок · PDB',route['repaired_pdb'],file_name='receptor_modelled.pdb')
        st.caption('Для проверки качества модели необходим контрольный редокинг с известным лигандом и анализ RMSD.')
