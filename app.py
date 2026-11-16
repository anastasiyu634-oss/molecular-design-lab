"""Molecular Design Lab — спокойный интерфейс с пошаговым исследованием.

Реализованы: шаг 1 (PDB + геометрические контакты), шаг 2 (первичное
сравнение по номерам остатков). Шаги 3 и 4 — честные экраны-планы.
"""
from __future__ import annotations

from collections import defaultdict
from html import escape

import pandas as pd
import requests
import streamlit as st

from analysis import analyze_structure, common_residues, download_structure, parse_pdb_ids
from ui_styles import hero, hint, inject_styles, metric_card, render_header


st.set_page_config(
    page_title="Molecular Design Lab",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)
inject_styles()

STAGES = [
    ("structures", "01 · Структуры"),
    ("compare", "02 · Анализ"),
    ("design", "03 · Дизайн"),
    ("docking", "04 · Докинг"),
]

if "stage" not in st.session_state:
    st.session_state.stage = "structures"
if "summary" not in st.session_state:
    st.session_state.summary = []
if "ligands" not in st.session_state:
    st.session_state.ligands = []
if "load_errors" not in st.session_state:
    st.session_state.load_errors = []
if "loaded_ids" not in st.session_state:
    st.session_state.loaded_ids = []


@st.cache_data(ttl=24 * 60 * 60, show_spinner=False)
def cached_download(pdb_id: str) -> str:
    return download_structure(pdb_id)


def choose_stage(stage: str) -> None:
    st.session_state.stage = stage


def public_table(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame([{k: v for k, v in row.items() if not k.startswith("_")} for row in rows])


def render_navigation() -> None:
    st.markdown("<div class='lab-section'>Ход исследования</div>", unsafe_allow_html=True)
    current = st.session_state.stage
    cols = st.columns(4, gap="small")
    for col, (key, label) in zip(cols, STAGES):
        with col:
            st.button(
                label,
                key=f"nav_{key}",
                type="primary" if current == key else "secondary",
                use_container_width=True,
                on_click=choose_stage,
                args=(key,),
            )
    st.caption("Выберите этап. Сейчас работают загрузка структур и первичное сравнение контактов; генерация и докинг пока в разработке.")


def render_summary() -> None:
    summary = st.session_state.summary
    ligands = st.session_state.ligands
    if not summary:
        return
    st.markdown("<div class='lab-section'>Результат загрузки</div>", unsafe_allow_html=True)
    a, b, c = st.columns(3, gap="small")
    with a:
        metric_card("Структур загружено", len(summary), "Из базы RCSB PDB")
    with b:
        metric_card("Кандидатов в лиганды", len(ligands), "Крупные неводные компоненты")
    with c:
        metric_card("Белковых цепей", sum(s["Белковых цепей"] for s in summary), "Суммарно во всех структурах")

    if summary:
        st.markdown("**Загруженные структуры**")
        st.dataframe(pd.DataFrame(summary), hide_index=True, use_container_width=True)
    if ligands:
        st.markdown("**Найденные компоненты рядом с белком**")
        for i, ligand in enumerate(ligands):
            st.markdown(
                f"<span class='lab-chip'>{escape(ligand['PDB'])} · {escape(ligand['Код компонента'])}</span>",
                unsafe_allow_html=True,
            )
            if i >= 6:
                break
        with st.expander("Показать подробную таблицу кандидатов и контактов"):
            df = public_table(ligands)
            st.dataframe(df, hide_index=True, use_container_width=True)
            st.download_button(
                "Скачать CSV",
                data=df.to_csv(index=False).encode("utf-8-sig"),
                file_name="pdb_nearby_residues.csv",
                mime="text/csv",
            )
        hint(
            "Мы нашли атомно-близкие аминокислотные остатки (до 4 Å), но это ещё не доказанные химические связи. "
            "Обнаруженные компоненты — кандидаты в лиганды; их биологическую роль нужно проверять отдельно."
        )
        if st.button("Перейти к сравнению →", type="primary", key="go_compare"):
            st.session_state.stage = "compare"
            st.rerun()
    else:
        st.info("Крупных неводных компонентов в загруженных моделях не найдено. Можно попробовать другие PDB ID.")


def render_structures() -> None:
    st.markdown("<div class='lab-kicker'>Этап 01 / Подготовка</div>", unsafe_allow_html=True)
    st.subheader("Выбор белковых структур")
    st.write("Введите один или несколько PDB ID. Сайт загрузит структуры и найдёт кандидатов в лиганды и ближайшие аминокислотные остатки белка.")

    with st.form("pdb_form"):
        raw = st.text_input(
            "Идентификаторы PDB",
            value="1A9U, 1OUY, 1W7H",
            help="Разделяйте коды запятой, пробелом или точкой с запятой. Максимум 10 структур.",
        )
        submitted = st.form_submit_button("Загрузить структуры", type="primary", use_container_width=False)

    if submitted:
        try:
            pdb_ids = parse_pdb_ids(raw)
        except ValueError as exc:
            st.error(str(exc))
        else:
            summary, ligands, errors = [], [], []
            for pdb_id in pdb_ids:
                try:
                    with st.spinner(f"Загрузка и анализ {pdb_id}…"):
                        cif = cached_download(pdb_id)
                        found, metadata = analyze_structure(pdb_id, cif)
                    summary.append({"PDB": pdb_id, **metadata, "Кандидатов в лиганды": len(found)})
                    ligands.extend(found)
                except requests.HTTPError as exc:
                    status = exc.response.status_code if exc.response is not None else "?"
                    errors.append(f"{pdb_id}: сервер PDB вернул HTTP {status}. Проверьте идентификатор.")
                except requests.RequestException:
                    errors.append(f"{pdb_id}: не удалось загрузить структуру. Проверьте соединение с интернетом.")
                except Exception as exc:
                    errors.append(f"{pdb_id}: ошибка разбора структуры — {exc}")

            st.session_state.summary = summary
            st.session_state.ligands = ligands
            st.session_state.loaded_ids = [s["PDB"] for s in summary]
            st.session_state.load_errors = errors

    for err in st.session_state.load_errors:
        st.warning(err)

    if st.session_state.summary:
        render_summary()
    else:
        st.markdown(
            "<div class='lab-placeholder'><b>Пока нет загруженных структур.</b> "
            "Введи PDB ID и нажми «Загрузить структуры». Здесь появятся понятные результаты и ссылка на подробные данные.</div>",
            unsafe_allow_html=True,
        )
        hint("Для первого теста используем структуры 1A9U, 1OUY и 1W7H из нашей учебной работы с p38α MAPK.")


def render_compare() -> None:
    st.markdown("<div class='lab-kicker'>Этап 02 / Первичный анализ</div>", unsafe_allow_html=True)
    st.subheader("Сравнение лигандов")
    st.write("Выберем по одному компоненту из каждой структуры и посмотрим, какие аминокислотные остатки находятся поблизости.")
    if not st.session_state.summary:
        st.info("Сначала загрузите структуры на этапе 01.")
        if st.button("Открыть загрузку", type="primary"):
            st.session_state.stage = "structures"
            st.rerun()
        return

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in st.session_state.ligands:
        grouped[row["PDB"]].append(row)

    selected = []
    for pdb_id in st.session_state.loaded_ids:
        variants = grouped[pdb_id]
        if not variants:
            st.warning(f"{pdb_id}: подходящих крупных компонентов не обнаружено.")
            continue
        index = st.selectbox(
            f"Компонент для структуры {pdb_id}",
            options=list(range(len(variants))),
            format_func=lambda i, variants=variants: (
                f"{variants[i]['Код компонента']} · цепь {variants[i]['Цепь']} · "
                f"{variants[i]['Тяжёлых атомов']} тяжёлых атомов"
            ),
            key=f"select_{pdb_id}",
        )
        selected.append(variants[index])

    if not selected:
        st.info("Нет выбранных компонентов для сравнения.")
        return

    st.markdown("**Обзор выбранных компонентов**")
    overview = [
        {
            "Структура": row["PDB"],
            "Компонент": row["Код компонента"],
            "Тяжёлых атомов": row["Тяжёлых атомов"],
            "Близких остатков (≤ 4 Å)": row["Близких остатков (≤ 4 Å)"],
        }
        for row in selected
    ]
    st.dataframe(pd.DataFrame(overview), hide_index=True, use_container_width=True)

    common = common_residues(selected)
    if len(selected) >= 2:
        st.markdown("**Предварительно совпадающие остатки**")
        if common:
            st.markdown(" ".join(f"<span class='lab-chip'>{escape(item)}</span>" for item in common), unsafe_allow_html=True)
            st.caption(f"Совпали {len(common)} обозначений по названию и номеру аминокислотного остатка.")
        else:
            st.info("Общих обозначений близких остатков пока не обнаружено.")
    else:
        st.info("Для сравнения нужно выбрать компоненты хотя бы из двух структур.")

    with st.expander("Посмотреть ближайшие остатки для каждой структуры"):
        for row in selected:
            st.markdown(f"**{row['PDB']} · {row['Код компонента']}**")
            st.write(row["Ближайшие остатки"] or "Близких остатков не найдено")
        st.download_button(
            "Скачать выбранные данные (CSV)",
            data=public_table(selected).to_csv(index=False).encode("utf-8-sig"),
            file_name="selected_pdb_contacts.csv",
            mime="text/csv",
        )

    hint(
        "Это предварительное сравнение по номерам аминокислот, без пространственного выравнивания белков. "
        "Одинаковые номера не гарантируют эквивалентное положение в разных структурах. "
        "Следующая научная доработка — выравнивание белков и распознавание типов взаимодействий "
        "(водородных, гидрофобных и других); полноценного фармакофора здесь пока нет."
    )


def render_design() -> None:
    st.markdown("<div class='lab-kicker'>Этап 03 / В разработке</div>", unsafe_allow_html=True)
    st.subheader("Проектирование новых молекул")
    st.write("Здесь будет модуль de novo design. Сначала нужно корректно построить фармакофор и проверить его на известных лигандах.")
    st.markdown(
        "<div class='lab-placeholder'><b>План модуля:</b><br>"
        "1. Выбор подтверждённых фармакофорных признаков.<br>"
        "2. Генерация нескольких химически допустимых структур.<br>"
        "3. Проверка валентности, основных свойств и построение 3D-конформеров.<br>"
        "4. Сохранение выбранных кандидатов для последующей проверки.</div>",
        unsafe_allow_html=True,
    )
    hint("Пока генерация молекул НЕ запускается. Здесь нет выдуманных результатов — этот экран показывает будущий этап исследования.")


def render_docking() -> None:
    st.markdown("<div class='lab-kicker'>Этап 04 / В разработке</div>", unsafe_allow_html=True)
    st.subheader("Молекулярный докинг")
    st.write("На этом этапе мы добавим автоматическую подготовку белка и молекулы, расчёты докинга и сравнение результатов по одному протоколу.")
    st.markdown(
        "<div class='lab-placeholder'><b>План модуля:</b><br>"
        "1. Проверка и подготовка структур белка и лиганда.<br>"
        "2. Выбор центра и размера области докинга.<br>"
        "3. Запуск AutoDock Vina на вычислительном сервере.<br>"
        "4. Просмотр предсказанных поз и энергий связывания.<br>"
        "5. Экспорт параметров, результатов и изображений.</div>",
        unsafe_allow_html=True,
    )
    hint("Докинг пока НЕ подключён. Энергии связывания не рассчитываются и не выводятся как реальные результаты.")


render_header()
hero(
    "От структуры белка к новой молекуле",
    "Исследование / p38α MAPK",
    "Спокойное рабочее пространство: загружаем структуры, изучаем взаимодействия и постепенно переходим к фармакофору, дизайну и докингу.",
)
render_navigation()
st.divider()
if st.session_state.stage == "structures":
    render_structures()
elif st.session_state.stage == "compare":
    render_compare()
elif st.session_state.stage == "design":
    render_design()
else:
    render_docking()

st.divider()
st.caption(
    "Molecular Design Lab · локальная учебно-исследовательская версия. "
    "Геометрические контакты и компьютерные прогнозы не подтверждают биологическую активность вещества."
)
