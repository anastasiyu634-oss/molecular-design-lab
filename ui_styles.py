"""Общий визуальный стиль. Его можно менять отдельно от научных расчётов."""
from __future__ import annotations

from html import escape

import streamlit as st


CSS = """
<style>
:root { --lab-olive:#3E5546; --lab-ink:#26372E; --lab-muted:#718076;
 --lab-page:#F8F7F3; --lab-line:#E1E6DF; --lab-mist:#EFF3ED; }
html, body, [class*="css"], .stApp {font-family: 'Segoe UI', Inter, system-ui, sans-serif;}
.stApp, [data-testid="stAppViewContainer"] {background:var(--lab-page); color:var(--lab-ink);}
[data-testid="stHeader"] {background:transparent;}
.block-container {max-width:1160px; padding-top:2.2rem; padding-bottom:5rem;}
h1,h2,h3 {letter-spacing:-0.025em; color:var(--lab-ink);}
h1 {font-weight:700 !important;}
[data-testid="stMarkdownContainer"] p {line-height:1.55;}
[data-testid="stMetric"] {background:#FFF;border:1px solid var(--lab-line);border-radius:14px;padding:14px;}
div.stButton > button[kind="primary"] {background:#3E5546; border:1px solid #3E5546; color:white; font-weight:600;}
div.stButton > button[kind="primary"]:hover {background:#304538; border-color:#304538; color:#fff;}
div.stButton > button[kind="secondary"] {background:white; border:1px solid var(--lab-line); color:#45594B;}
div.stButton > button[kind="secondary"]:hover {border-color:#9FB2A3;color:#304538;}
[data-testid="stExpander"], [data-testid="stDataFrame"] {border-radius:12px;}
.lab-brand {font-size:1.04rem;font-weight:750;letter-spacing:-.02em;color:#334B3C;}
.lab-kicker {font-size:.74rem;color:#647E69;font-weight:700;letter-spacing:.10em;text-transform:uppercase;}
.lab-hero {background:#fff;border:1px solid var(--lab-line);border-radius:18px;padding:24px 26px;margin-bottom:18px;}
.lab-hero h1 {font-size:2rem;line-height:1.18;margin:8px 0 8px;color:#293B30;}
.lab-hero p {font-size:.98rem;color:#65766B;max-width:760px;margin:0;}
.lab-intro {font-size:.88rem;color:#708274;margin-top:5px;}
.lab-card {border:1px solid var(--lab-line);border-radius:14px;background:#fff;padding:18px;margin:7px 0;}
.lab-card-label {font-size:.77rem;color:#78867C;margin-bottom:5px;}
.lab-card-value {font-weight:750;font-size:1.6rem;color:var(--lab-olive);line-height:1.3;}
.lab-card-foot {font-size:.78rem;color:#809083;margin-top:5px;}
.lab-section {font-weight:750;color:#2D4434;font-size:1.24rem;margin:14px 0 8px;}
.lab-note {background:#EDF3EA;color:#4D6753;border-radius:12px;padding:14px 16px;line-height:1.5;font-size:.88rem;margin-top:12px;}
.lab-placeholder {background:#fff;border:1px dashed #C3D1C5;border-radius:14px;padding:22px;color:#627668;}
.lab-chip {display:inline-block;background:#EAF0E7;padding:5px 9px;color:#45614A;border-radius:6px;font-size:.81rem;font-weight:650;margin-right:6px;}
.lab-rule {border-top:1px solid var(--lab-line);margin:14px 0;}
@media (max-width:700px){.block-container{padding-top:1.2rem}.lab-hero{padding:18px}.lab-hero h1{font-size:1.5rem}}
</style>
"""


def inject_styles() -> None:
    st.markdown(CSS, unsafe_allow_html=True)


def render_header() -> None:
    st.markdown("<div class='lab-brand'>◈ &nbsp; Molecular Design Lab</div>", unsafe_allow_html=True)
    st.caption("Рабочее пространство для исследований молекулярного связывания")
    st.markdown("<div class='lab-rule'></div>", unsafe_allow_html=True)


def hero(title: str, kicker: str, body: str) -> None:
    st.markdown(
        "<div class='lab-hero'>"
        f"<div class='lab-kicker'>{escape(kicker)}</div>"
        f"<h1>{escape(title)}</h1><p>{escape(body)}</p></div>",
        unsafe_allow_html=True,
    )


def metric_card(label: str, value: str | int, foot: str = "") -> None:
    st.markdown(
        "<div class='lab-card'>"
        f"<div class='lab-card-label'>{escape(str(label))}</div>"
        f"<div class='lab-card-value'>{escape(str(value))}</div>"
        f"<div class='lab-card-foot'>{escape(str(foot))}</div></div>",
        unsafe_allow_html=True,
    )


def hint(text: str) -> None:
    st.markdown(f"<div class='lab-note'>{escape(text)}</div>", unsafe_allow_html=True)
