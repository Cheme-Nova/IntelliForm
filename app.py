"""
IntelliForm v2.1 — ChemeNova LLC
AI-Powered Specialty Chemical Formulation Intelligence

v2.1 Changelog (reconciliation build):
  - Keeps v2.0 mobile-first UX: no sidebar, full tab navigation, responsive CSS
  - Restores ALL real modules from v1.5 (no _try_import fallback stubs):
      modules/verticals.py           → vertical selector + ingredient filtering
      modules/optimizer.py           → vertical-aware LP with max_conc
      modules/pareto_optimizer.py    → full Pareto frontier tab
      modules/bayesian_optimizer.py  → Bayesian GP optimizer
      modules/ecometrics.py          → real radar vs petrochemical baseline
      modules/qsar.py                → Mordred/GBR QSAR with active learning
      modules/regulatory.py          → REACH/EPA/COSMOS per-ingredient
      modules/vertical_regulatory.py → per-vertical regulatory reports
      modules/agents.py              → 4-agent swarm commentary
      modules/chem_utils.py          → RDKit molecular structure
      modules/certification_oracle.py→ CertificationOracle™
      modules/carbon_passport.py     → ISO 14067 Carbon Passport
      modules/carbon_credits.py      → carbon credit calculator
      modules/reformulation_intelligence.py → Reformulation Lab
      modules/memory_network.py      → Formulation Memory Network
      modules/pdf_proposal.py        → branded PDF proposals
      modules/stability.py           → stability & viscosity
      modules/notifications.py       → email (SendGrid)
      modules/persistence.py         → Supabase save/load
      modules/analytics.py           → PostHog tracking
      modules/pharma.py              → full ICH/BCS pharma deep dive
  - New v2.0 features kept: CBAM calculator, supply chain AI, Pro CTA, pilot mailto
  - Pharma Deep Dive tab only visible when Pharmaceutical vertical selected
  - Citations: ChemRxiv DOI + NJIT Showcase 2025 + UIC Indigo 2025
"""

import os
import json
import re
from datetime import datetime

import pandas as pd
import numpy as np
import streamlit as st

# ── PAGE CONFIG (must be first Streamlit call) ────────────────────────────────
st.set_page_config(
    page_title="IntelliForm | ChemeNova",
    page_icon="⚗",
    layout="wide",
    initial_sidebar_state="collapsed",
    menu_items={
        "Get Help": "mailto:shehan@chemenova.com",
        "Report a bug": "https://github.com/Cheme-Nova/IntelliForm/issues",
        "About": (
            "**IntelliForm v2.1** — AI-powered green chemistry formulation.\n"
            "© 2026 ChemeNova LLC · MIT License\n"
            "DOI: 10.26434/chemrxiv.15000857"
        ),
    },
)

# ── ENVIRONMENT ───────────────────────────────────────────────────────────────
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

GROQ_API_KEY      = os.getenv("GROQ_API_KEY", "")
ANTHROPIC_API_KEY = os.getenv("ANTHROPIC_API_KEY", "")
OPENAI_API_KEY    = os.getenv("OPENAI_API_KEY", "")
SUPABASE_URL      = os.getenv("SUPABASE_URL", "")
SUPABASE_KEY      = os.getenv("SUPABASE_ANON_KEY", "")
POSTHOG_KEY       = os.getenv("POSTHOG_API_KEY", "")
DEBUG             = os.getenv("DEBUG", "").lower() in ("1", "true", "yes")

# EU ETS reference price — update quarterly
CBAM_RATE_EUR_PER_TON = 65.0  # €/tCO₂e, March 2026

# ── MOBILE-FIRST CSS (from v2.0, kept verbatim) ───────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@300;400;500;700&family=DM+Sans:wght@300;400;500;600;700&family=Syne:wght@700;800&display=swap');

html, body, [data-testid="stAppViewContainer"] {
    background: #050E1F !important;
    color: #f1f5f9 !important;
    font-family: 'DM Sans', sans-serif !important;
}
[data-testid="stHeader"]  { background: rgba(5,14,31,0.95) !important; }
footer                    { visibility: hidden !important; }
.block-container          { padding: 1rem 1.5rem 2rem !important; max-width: 1400px !important; }

/* Hide sidebar — tabs only */
[data-testid="stSidebar"]        { display: none !important; }
[data-testid="collapsedControl"] { display: none !important; }

/* ── TABS ── */
.stTabs [data-baseweb="tab-list"] {
    gap: 4px;
    background: rgba(10,22,40,0.85);
    padding: 6px;
    border-radius: 10px;
    border: 1px solid rgba(30,58,95,0.55);
    overflow-x: auto;
    -webkit-overflow-scrolling: touch;
    flex-wrap: nowrap;
    scrollbar-width: none;
}
.stTabs [data-baseweb="tab-list"]::-webkit-scrollbar { display: none; }
.stTabs [data-baseweb="tab"] {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.72rem !important;
    font-weight: 600 !important;
    letter-spacing: 0.05em !important;
    text-transform: uppercase !important;
    color: #94a3b8 !important;
    background: transparent !important;
    border-radius: 7px !important;
    padding: 8px 16px !important;
    white-space: nowrap !important;
    border: none !important;
    transition: all 0.2s !important;
    min-width: fit-content !important;
}
.stTabs [aria-selected="true"] {
    background: linear-gradient(135deg, #0D9488, #0f766e) !important;
    color: #fff !important;
    box-shadow: 0 2px 12px rgba(13,148,136,0.4) !important;
}
[data-testid="stTabContent"] { padding-top: 24px !important; }

/* ── INPUTS ── */
.stTextInput > div > div > input,
.stTextArea > div > div > textarea,
.stNumberInput > div > div > input {
    background: rgba(10,22,40,0.9) !important;
    border: 1px solid rgba(30,58,95,0.7) !important;
    border-radius: 8px !important;
    color: #f1f5f9 !important;
    font-family: 'DM Sans', sans-serif !important;
}
.stTextInput > div > div > input:focus,
.stTextArea > div > div > textarea:focus {
    border-color: #0D9488 !important;
    box-shadow: 0 0 0 2px rgba(13,148,136,0.2) !important;
}
label { color: #94a3b8 !important; font-size: 0.85rem !important; }

/* ── SELECTBOX ── */
.stSelectbox > div > div {
    background: rgba(10,22,40,0.9) !important;
    border: 1px solid rgba(30,58,95,0.7) !important;
    border-radius: 8px !important;
    color: #f1f5f9 !important;
}

/* ── BUTTONS ── */
.stButton > button {
    background: linear-gradient(135deg, #0D9488, #0f766e) !important;
    color: #fff !important;
    border: none !important;
    border-radius: 8px !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-weight: 700 !important;
    font-size: 0.78rem !important;
    letter-spacing: 0.08em !important;
    text-transform: uppercase !important;
    padding: 10px 20px !important;
    transition: all 0.2s !important;
    width: 100% !important;
}
.stButton > button:hover {
    transform: translateY(-1px) !important;
    box-shadow: 0 8px 24px rgba(13,148,136,0.4) !important;
}

/* ── METRICS ── */
[data-testid="metric-container"] {
    background: rgba(10,22,40,0.8) !important;
    border: 1px solid rgba(30,58,95,0.55) !important;
    border-radius: 10px !important;
    padding: 16px !important;
    transition: border-color 0.2s !important;
}
[data-testid="metric-container"]:hover {
    border-color: #0D9488 !important;
}
[data-testid="stMetricLabel"] {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.65rem !important;
    font-weight: 600 !important;
    letter-spacing: 0.1em !important;
    text-transform: uppercase !important;
    color: #64748b !important;
    white-space: nowrap !important;
    overflow: hidden !important;
    text-overflow: ellipsis !important;
}
[data-testid="stMetricValue"] {
    font-family: 'Syne', sans-serif !important;
    font-size: 1.3rem !important;
    font-weight: 800 !important;
    color: #f1f5f9 !important;
    white-space: nowrap !important;
    overflow: hidden !important;
    text-overflow: ellipsis !important;
}

/* ── EXPANDERS ── */
[data-testid="stExpander"] {
    background: rgba(10,22,40,0.7) !important;
    border: 1px solid rgba(30,58,95,0.5) !important;
    border-radius: 8px !important;
}
[data-testid="stExpander"] summary {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.8rem !important;
    font-weight: 600 !important;
    color: #94a3b8 !important;
    padding: 10px 14px !important;
}

/* ── ALERTS ── */
[data-testid="stSuccess"] { background: rgba(13,148,136,0.08) !important; border-left: 3px solid #0D9488 !important; }
[data-testid="stWarning"] { background: rgba(217,119,6,0.08) !important; border-left: 3px solid #D97706 !important; }
[data-testid="stError"]   { background: rgba(239,68,68,0.08) !important; border-left: 3px solid #ef4444 !important; }
[data-testid="stInfo"]    { background: rgba(59,130,246,0.08) !important; border-left: 3px solid #3b82f6 !important; }

/* ── DATAFRAMES ── */
[data-testid="stDataFrame"] {
    border: 1px solid rgba(30,58,95,0.5) !important;
    border-radius: 8px !important;
    overflow: hidden !important;
}
[data-testid="stDataFrame"] th {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.68rem !important;
    letter-spacing: 0.06em !important;
    text-transform: uppercase !important;
    background: #0A1628 !important;
    color: #64748b !important;
}

/* ── CUSTOM CARDS ── */
.cn-card { background: rgba(10,22,40,0.8); border: 1px solid rgba(30,58,95,0.55); border-radius: 10px; padding: 20px 24px; margin-bottom: 16px; }
.cn-result-card { background: rgba(10,22,40,0.9); border: 1px solid rgba(13,148,136,0.35); border-radius: 10px; padding: 24px; margin-bottom: 16px; }
.cn-pilot-cta { background: linear-gradient(135deg, #0D9488, #0f766e); border-radius: 10px; padding: 28px 24px; text-align: center; margin: 20px 0; }
.cn-cbam-card { background: rgba(217,119,6,0.06); border: 1px solid rgba(217,119,6,0.25); border-radius: 10px; padding: 20px 24px; margin-bottom: 16px; }
.cert-pass { background: rgba(13,148,136,0.12); border: 1px solid rgba(13,148,136,0.3); border-radius: 6px; padding: 6px 12px; display: inline-block; color: #14b8a8; font-family: 'IBM Plex Mono', monospace; font-size: 0.74rem; margin: 3px; }
.cert-warn { background: rgba(245,158,11,0.1); border: 1px solid rgba(245,158,11,0.25); border-radius: 6px; padding: 6px 12px; display: inline-block; color: #f59e0b; font-family: 'IBM Plex Mono', monospace; font-size: 0.74rem; margin: 3px; }
.cert-fail { background: rgba(239,68,68,0.08); border: 1px solid rgba(239,68,68,0.22); border-radius: 6px; padding: 6px 12px; display: inline-block; color: #f87171; font-family: 'IBM Plex Mono', monospace; font-size: 0.74rem; margin: 3px; }

/* ── SCROLLBAR ── */
::-webkit-scrollbar { width: 5px; height: 5px; }
::-webkit-scrollbar-track { background: #050E1F; }
::-webkit-scrollbar-thumb { background: #1e3a5f; border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: #0D9488; }

@media (max-width: 768px) {
    .block-container { padding: 0.75rem 0.8rem 2rem !important; }
    .stTabs [data-baseweb="tab"] { font-size: 0.62rem !important; padding: 7px 10px !important; }
}
</style>
""", unsafe_allow_html=True)

# ── REAL MODULE IMPORTS ───────────────────────────────────────────────────────
# All modules are imported directly — no _try_import fallback stubs.
# If a module fails to import, a clear error is shown rather than silent fallback.

from modules.analytics        import track, identify_user, get_session_id
from modules.tiers             import detect_tier
from modules.llm_parser        import parse_request
from modules.optimizer         import run_optimization
from modules.pareto_optimizer  import run_pareto_optimization, pareto_frontier_dataframe
from modules.bayesian_optimizer import run_bayesian_optimization, BAYES_OK
from modules.agents            import run_agent_swarm
from modules.chem_utils        import draw_mol, enrich_db
from modules.ecometrics        import compute_ecometrics, ecometrics_radar_data
from modules.qsar              import initialize_models, predict_properties, submit_feedback
from modules.regulatory        import get_blend_report, regulatory_table_df
from modules.vertical_regulatory import generate_vertical_regulatory_report
from modules.verticals         import (VERTICAL_OPTIONS, get_profile,
                                        filter_db_by_vertical, get_vertical_constraints)
from modules.certification_oracle import run_certification_oracle, CERTIFICATION_PROFILES
from modules.carbon_passport   import generate_carbon_passport, passport_to_json
from modules.carbon_credits    import calculate_carbon_credits
from modules.reformulation_intelligence import (
    run_reformulation_intelligence, FAILURE_TYPES, ReformulationReport)
from modules.memory_network    import get_memory_network, FormulationMemoryNetwork
from modules.pdf_proposal      import generate_proposal_pdf
from modules.stability         import predict_stability
from modules.notifications     import send_pilot_booking_confirmation
from modules.persistence       import (save_project, load_projects, save_booking,
                                        save_feedback as db_save_feedback,
                                        is_connected, MIGRATION_SQL)
from modules.pharma            import (run_pharma_deep_dive, BCS_STRATEGIES,
                                        DOSAGE_FORMS, REGULATORY_PATHWAYS,
                                        ICH_STABILITY_ZONES)

# ── LLM AVAILABILITY ──────────────────────────────────────────────────────────
_available_llms = []
if GROQ_API_KEY:      _available_llms.append("Groq (llama-3.1-8b) — Free")
if ANTHROPIC_API_KEY: _available_llms.append("Anthropic (claude-sonnet) — Requires credits")
if OPENAI_API_KEY:    _available_llms.append("OpenAI (gpt-4o-mini) — Requires credits")
_available_llms.append("Regex fallback — No API needed")

_llm_to_provider = {
    "Groq (llama-3.1-8b) — Free": "groq",
    "Anthropic (claude-sonnet) — Requires credits": "anthropic",
    "OpenAI (gpt-4o-mini) — Requires credits": "openai",
    "Regex fallback — No API needed": "regex",
}

# ── DB LOAD ───────────────────────────────────────────────────────────────────
@st.cache_data(ttl=3600, show_spinner=False)
def load_db() -> pd.DataFrame:
    df = pd.read_csv("data/ingredients_db.csv")
    if "Vertical" not in df.columns:
        df["Vertical"] = "personal_care"
    return enrich_db(df)

@st.cache_resource
def load_models(n: int, db_hash: str):
    db = load_db()
    return initialize_models(db)

# ── SESSION STATE ─────────────────────────────────────────────────────────────
def _init_state():
    defaults = {
        "last_result":      None,
        "last_parsed":      None,
        "last_eco":         None,
        "last_pareto":      None,
        "last_reg":         None,
        "last_stability":   None,
        "last_carbon":      None,
        "last_v_reg":       None,
        "last_cert":        None,
        "last_passport":    None,
        "last_refo":        None,
        "pharma_result":    None,
        "bayes_state":      None,
        "model_card":       None,
        "memory_net":       None,
        "projects":         [],
        "blend_history":    [],
        "pilot_sent":       False,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

_init_state()

# ── BOOTSTRAP ─────────────────────────────────────────────────────────────────
ingredients_db = load_db()
_db_hash = str(len(ingredients_db)) + str(
    ingredients_db["Ingredient"].iloc[0] if len(ingredients_db) > 0 else "")

if st.session_state.model_card is None:
    st.session_state.model_card = load_models(len(ingredients_db), _db_hash)

if st.session_state.memory_net is None:
    st.session_state.memory_net = get_memory_network()

mc = st.session_state.model_card
memory = st.session_state.memory_net
qsar_ok = mc and mc.sklearn_version != "unavailable"
_mordred_on = mc and getattr(mc, "mordred_active", False)

# ── TIER DETECTION ─────────────────────────────────────────────────────────────
current_tier = detect_tier()

# ── POSTHOG SESSION START (fires once per session) ────────────────────────────
if "session_tracked" not in st.session_state:
    session_id = get_session_id()
    identify_user(session_id, {"tier": current_tier, "vertical": "unset",
                                "app_version": "2.1", "db_size": len(ingredients_db)})
    track("session_started", {"tier": current_tier, "app_version": "2.1",
                               "db_size": len(ingredients_db),
                               "mordred": bool(_mordred_on)})
    st.session_state["session_tracked"] = True

# Load persisted projects from Supabase
if not st.session_state.projects and is_connected():
    try:
        stored = load_projects(get_session_id(), limit=20)
        if stored:
            st.session_state.projects = stored
    except Exception:
        pass

# ── HEADER ────────────────────────────────────────────────────────────────────
def _render_header():
    c1, c2, c3 = st.columns([2, 5, 2])
    with c1:
        st.markdown(
            '<div style="padding:10px 0;font-family:\'IBM Plex Mono\',monospace;font-size:1rem;'
            'font-weight:700;color:#f1f5f9">Cheme<span style="color:#0D9488">Nova</span>'
            '<span style="font-size:0.6rem;color:#64748b;background:rgba(13,148,136,0.1);'
            'border:1px solid rgba(13,148,136,0.2);padding:2px 8px;border-radius:3px;'
            'margin-left:10px;vertical-align:middle">v2.1</span></div>',
            unsafe_allow_html=True)
    with c2:
        st.markdown(
            f'<div style="text-align:center;padding:10px 0;font-family:\'IBM Plex Mono\',monospace;'
            f'font-size:0.82rem;font-weight:700;color:#f1f5f9">'
            f'IntelliForm &nbsp;·&nbsp; '
            f'<span style="color:#0D9488">{len(ingredients_db):,}</span> ingredients &nbsp;·&nbsp; '
            f'<span style="color:#0D9488">{"Mordred" if _mordred_on else "GBR"}</span> QSAR &nbsp;·&nbsp; '
            f'<span style="color:{"#0D9488" if qsar_ok else "#64748b"}">'
            f'{"ML Active" if qsar_ok else "ML Offline"}</span> &nbsp;·&nbsp; '
            f'<span style="color:{"#0D9488" if is_connected() else "#64748b"}">'
            f'{"Supabase ✓" if is_connected() else "Session only"}</span></div>',
            unsafe_allow_html=True)
    with c3:
        st.markdown(
            '<div style="text-align:right;padding:10px 0">'
            '<a href="https://chemenova.com" target="_blank" style="font-family:\'IBM Plex Mono\','
            'monospace;font-size:0.65rem;color:#0D9488;text-decoration:none;border:1px solid '
            'rgba(13,148,136,0.3);padding:5px 12px;border-radius:4px">chemenova.com →</a></div>',
            unsafe_allow_html=True)
    st.divider()

_render_header()

# ── PILOT CTA HTML ────────────────────────────────────────────────────────────
def _pilot_mailto(result=None) -> str:
    if result:
        blend = result.blend if hasattr(result, "blend") else {}
        lines = "\n".join(f"  - {k}: {v}%" for k, v in blend.items())
        cost  = result.cost_per_kg if hasattr(result, "cost_per_kg") else 0
        bio   = result.bio_pct if hasattr(result, "bio_pct") else 0
    else:
        lines, cost, bio = "", 0, 0
    subj = "Pilot Batch Request — IntelliForm Formulation"
    body = (
        f"Hello ChemeNova/ChemRich team,\n\nI would like to request a pilot batch.\n\n"
        f"Cost: ${cost}/kg\nBio: {bio}%\n\nBlend:\n{lines}\n\nBest regards"
    )
    return (f"mailto:shehan@chemenova.com"
            f"?subject={subj.replace(' ','%20')}"
            f"&body={body.replace(chr(10),'%0A').replace(' ','%20')}")

def _pilot_cta_html(result=None) -> str:
    mailto = _pilot_mailto(result)
    return f"""
    <div class="cn-pilot-cta">
        <div style="color:rgba(255,255,255,0.7);font-family:'IBM Plex Mono',monospace;
             font-size:0.65rem;letter-spacing:0.18em;text-transform:uppercase;margin-bottom:6px">
            Ready to manufacture this formulation?
        </div>
        <div style="color:#fff;font-size:1.05rem;font-weight:500;margin-bottom:4px">
            200 kg pilot batch via ChemRich Global
        </div>
        <div style="color:rgba(255,255,255,0.65);font-size:0.84rem;margin-bottom:18px">
            No-cost reformulation guarantee if certification criteria aren't met.
        </div>
        <a href="{mailto}" style="background:#fff;color:#0D9488;padding:13px 32px;
           border-radius:6px;font-family:'IBM Plex Mono',monospace;font-weight:700;
           font-size:0.8rem;text-transform:uppercase;letter-spacing:0.08em;
           text-decoration:none;display:inline-block">Request Pilot Batch →</a>
        <div style="color:rgba(255,255,255,0.4);font-size:0.7rem;margin-top:14px">
            NJ / NY / PA manufacturing · 21-day lead time · $6,500–$9,500
        </div>
    </div>"""

def _cert_badge(cert: str, prob: float) -> str:
    cls  = "cert-pass" if prob >= 0.8 else "cert-warn" if prob >= 0.6 else "cert-fail"
    icon = "✓" if prob >= 0.8 else "~" if prob >= 0.6 else "✗"
    return f'<span class="{cls}">{icon} {cert} {prob*100:.0f}%</span>'

# ── TAB DEFINITIONS ───────────────────────────────────────────────────────────
# Determine vertical selection (needed for pharma tab visibility)
# We need to read vertical_sel before building tabs, so we use a top-level selectbox
# hidden inside an invisible container and then re-use it everywhere.

_VERTICAL_DISPLAY = list(VERTICAL_OPTIONS.keys())

# ── TOP-LEVEL VERTICAL + LLM CONTROLS ─────────────────────────────────────────
_ctrl_cols = st.columns([3, 2, 2])
with _ctrl_cols[0]:
    selected_vertical = st.selectbox(
        "🎯 Industry Vertical",
        options=_VERTICAL_DISPLAY,
        format_func=lambda x: VERTICAL_OPTIONS[x],
        key="global_vertical",
        help="Filters ingredients and constraints for your specific industry."
    )
with _ctrl_cols[1]:
    selected_llm_label = st.selectbox(
        "🤖 AI Model",
        options=_available_llms,
        key="global_llm",
    )
    os.environ["LLM_PROVIDER"] = _llm_to_provider.get(selected_llm_label, "groq")
with _ctrl_cols[2]:
    _opt_options = ["LP (fast)", "Pareto"]
    if BAYES_OK:
        _opt_options.append("Bayesian")
    _opt_mode = st.radio(
        "⚙️ Optimizer",
        _opt_options,
        horizontal=True,
        key="global_opt",
    )
    if not BAYES_OK:
        st.caption("⚠ Bayesian unavailable — install scikit-optimize")

vprofile = get_profile(selected_vertical)
if vprofile:
    st.caption(f"{vprofile.description} · Frameworks: {' · '.join(vprofile.regulatory_frameworks[:3])}")

st.divider()

# Build tab list — pharma deep dive only when pharma vertical selected
_base_tabs = [
    "⚗ Formulate",
    "🏆 Certify",
    "🌿 EcoMetrics",
    "📋 Regulatory",
    "📈 Pareto",
    "🔬 QSAR",
    "🧪 Stability",
    "🌍 Carbon",
    "🔁 Reformulation",
    "🧠 Memory",
    "📄 Proposal",
    "📊 History",
    "⚡ Pro",
]
_show_pharma = (selected_vertical == "pharmaceutical")
if _show_pharma:
    _all_tabs = _base_tabs + ["💊 Pharma"]
    t_form, t_cert, t_eco, t_reg, t_pareto, t_qsar, t_stab, t_carbon, t_refo, t_mem, t_prop, t_hist, t_pro, t_pharma = st.tabs(_all_tabs)
else:
    t_form, t_cert, t_eco, t_reg, t_pareto, t_qsar, t_stab, t_carbon, t_refo, t_mem, t_prop, t_hist, t_pro = st.tabs(_base_tabs)
    t_pharma = None

# ─────────────────────────────────────────────────────────────────────────────
# TAB 1 — FORMULATE
# ─────────────────────────────────────────────────────────────────────────────
with t_form:
    st.subheader("⚗ Describe Your Formulation")
    st.caption("Plain English → optimized, regulatory-cleared blend in seconds.")

    # Advanced controls
    with st.expander("⚙ Formulation controls"):
        fc1, fc2, fc3 = st.columns(3)
        with fc1:
            max_conc   = st.slider("Max single ingredient %", 30, 100, 70, 5)
            batch_size = st.selectbox("Pilot batch size (kg)", [200, 500, 1000, 2000, 5000], index=1)
        with fc2:
            customer_company = st.text_input("Company name (for PDF)", placeholder="Your Company Ltd.")
            uname = st.text_input("Your name", placeholder="Your Name")
        with fc3:
            uemail = st.text_input("Email (for notifications)", placeholder="you@company.com")
            customer_logo = st.file_uploader("White-label logo (PNG/JPG)", type=["png","jpg","jpeg"])

    vprofile_sel = get_profile(selected_vertical)
    if vprofile_sel and vprofile_sel.example_prompts:
        with st.expander("💡 Example prompts for this vertical"):
            for ep in vprofile_sel.example_prompts:
                if st.button(f"📋 {ep[:60]}…" if len(ep)>60 else f"📋 {ep}", key=f"ep_{ep[:20]}"):
                    st.session_state["prefill_prompt"] = ep

    default_prompt = st.session_state.get("prefill_prompt",
        "I need a mild foaming green surfactant for personal care, under $4/kg, "
        "at least 95% bio-based, EPA Safer Choice compatible")
    nl_input = st.text_area("What do you need?", value=default_prompt, height=120, key="nl_input")

    run_btn = st.button("⚗  GENERATE FORMULATION", type="primary", use_container_width=True)

    if run_btn and nl_input.strip():
        track("formulation_run", {"vertical": selected_vertical, "optimizer": _opt_mode})

        with st.spinner("🧠 Parsing request…"):
            parsed = parse_request(nl_input)
        st.session_state.last_parsed = parsed
        track("llm_parser_used", {"backend": parsed.parser_backend,
                                   "application": parsed.application_type,
                                   "vertical": selected_vertical})

        with st.expander("🧠 Parse Result", expanded=False):
            pc1, pc2, pc3, pc4 = st.columns(4)
            pc1.metric("Max Cost",   f"${parsed.max_cost}/kg")
            pc2.metric("Min Bio%",   f"{parsed.min_bio}%")
            pc3.metric("Min Perf",   str(parsed.min_perf))
            pc4.metric("Application",parsed.application_type.replace("_"," ").title()[:14])
            st.caption(f"**{parsed.parser_backend.upper()}**: {parsed.reasoning}")

        # Vertical-filtered DB + vertical-adjusted constraints
        v_db = filter_db_by_vertical(ingredients_db, selected_vertical)
        v_cost, v_bio, v_perf = get_vertical_constraints(
            selected_vertical, parsed.max_cost, parsed.min_bio, parsed.min_perf)

        # Run optimizer
        if _opt_mode == "Pareto":
            with st.spinner("📈 Pareto optimization…"):
                pareto = run_pareto_optimization(
                    v_db, max_cost=v_cost, min_bio=v_bio, min_perf=v_perf)
            st.session_state.last_pareto = pareto
            if not pareto.success:
                st.error(pareto.error_msg)
                st.stop()
            rec = pareto.recommended
            from modules.optimizer import OptResult
            result = OptResult(success=True, blend=rec.blend, cost_per_kg=rec.cost_per_kg,
                               bio_pct=rec.bio_pct, perf_score=rec.perf_score, status="Optimal")
            st.info(f"📈 {pareto.n_solutions} Pareto solutions · backend: {pareto.backend}")

        elif _opt_mode == "Bayesian":
            with st.spinner("🧪 Bayesian GP optimization…"):
                bayes_result, new_state = run_bayesian_optimization(
                    v_db, max_cost=v_cost, min_bio=v_bio, min_perf=v_perf,
                    state=st.session_state.bayes_state,
                    max_conc=max_conc/100, vertical=selected_vertical)
            st.session_state.bayes_state = new_state
            if not bayes_result.success:
                st.error(bayes_result.error_msg)
                st.stop()
            from modules.optimizer import OptResult
            result = OptResult(success=True, blend=bayes_result.blend,
                               cost_per_kg=bayes_result.cost_per_kg,
                               bio_pct=bayes_result.bio_pct, perf_score=bayes_result.perf_score,
                               status="Optimal")
            st.info(f"🧪 Bayesian iter {bayes_result.n_observations} · "
                    f"EI={bayes_result.expected_improvement:.4f} · σ={bayes_result.uncertainty:.4f}")
        else:
            with st.spinner("⚗ LP optimization…"):
                result = run_optimization(
                    v_db, max_cost=v_cost, min_bio=v_bio, min_perf=v_perf,
                    max_concentration=max_conc/100, vertical=selected_vertical)
            if not result.success:
                st.error(result.error_msg)
                st.stop()
            if result.relaxed:
                st.warning(f"⚠ Constraints relaxed × {result.relaxation_rounds}")
                track("constraints_relaxed", {"rounds": result.relaxation_rounds,
                                               "vertical": selected_vertical,
                                               "original_cost": v_cost,
                                               "original_bio": v_bio})

        st.session_state.last_result = result

        # Agent swarm
        with st.spinner("🤖 Agent swarm…"):
            for comment in run_agent_swarm(result, parsed):
                st.info(comment)

        # Compute downstream
        eco      = compute_ecometrics(result.blend, ingredients_db)
        reg      = get_blend_report(result.blend)
        v_reg    = generate_vertical_regulatory_report(result.blend, v_db, selected_vertical)
        stability= predict_stability(result.blend, ingredients_db)
        carbon   = calculate_carbon_credits(result.blend, ingredients_db, batch_kg=batch_size)
        cert     = run_certification_oracle(result.blend, v_db, selected_vertical, result.bio_pct)

        # Store all
        st.session_state.last_eco      = eco
        st.session_state.last_reg      = reg
        st.session_state.last_v_reg    = v_reg
        st.session_state.last_stability= stability
        st.session_state.last_carbon   = carbon
        st.session_state.last_cert     = cert

        # Memory network
        try:
            memory.record("formulation_generated", result.blend, selected_vertical,
                          metadata={"cost": result.cost_per_kg, "bio": result.bio_pct},
                          outcome=None)
        except Exception:
            pass

        # Blend history for comparison
        st.session_state.blend_history.append({
            "label": f"Run {len(st.session_state.blend_history)+1} — "
                     f"{parsed.application_type.title()} @ ${result.cost_per_kg}/kg",
            "blend": result.blend, "cost": result.cost_per_kg,
            "bio": result.bio_pct, "perf": result.perf_score, "input": nl_input,
        })
        if len(st.session_state.blend_history) > 10:
            st.session_state.blend_history.pop(0)

        # Results header
        rb1, rb2, rb3, rb4 = st.columns(4)
        rb1.metric("Cost/kg",      f"${result.cost_per_kg}")
        rb2.metric("Bio-based",    f"{result.bio_pct}%")
        rb3.metric("EcoScore",     f"{eco.eco_score:.0f}/100" if eco else "—")
        rb4.metric("Regulatory",   reg.overall_status if reg else "—")

        if stability and carbon:
            sb1, sb2, sb3, sb4 = st.columns(4)
            sb1.metric("Shelf Life",    stability.shelf_life_range)
            sb2.metric("Viscosity",     f"{stability.viscosity_cp:,.0f} cP")
            sb3.metric("CO2 Displaced", f"{carbon.co2_displaced_kg:.1f} kg/batch")
            sb4.metric("Carbon Value",  f"${carbon.credit_value_mid:.2f}/batch")

        # Certification quick badges
        if cert:
            badges = " ".join(_cert_badge(c, p.pass_probability)
                              for c, p in list(cert.predictions.items())[:6]
                              if "N/A" not in p.verdict)
            st.markdown(badges, unsafe_allow_html=True)

        # Ingredient cards
        with st.expander("📋 Blend composition + molecular structures", expanded=True):
            for ing, pct in result.blend.items():
                rows = ingredients_db[ingredients_db["Ingredient"] == ing]
                if rows.empty:
                    continue
                row = rows.iloc[0]
                ct, ci = st.columns([2, 1])
                with ct:
                    st.markdown(f"### {ing} — {pct}%")
                    st.caption(
                        f"**{row.get('Function','—')}** · ${row['Cost_USD_kg']}/kg · "
                        f"{row['Bio_based_pct']}% bio · Stock {row['Stock_kg']} kg")
                    qp = predict_properties(row["SMILES"])
                    st.caption(
                        f"🔬 QSAR: Bio {qp.biodegradability:.0f}% · "
                        f"Ecotox {qp.ecotoxicity:.1f}/10 · "
                        f"Perf {qp.performance:.0f} · "
                        f"{'ML' if qp.used_ml else 'rules'} · {qp.confidence} confidence")
                with ci:
                    img = draw_mol(row["SMILES"])
                    if img:
                        st.image(img, width=180)

        # Pilot CTA
        st.markdown(_pilot_cta_html(result), unsafe_allow_html=True)
        if st.button("📤 Book ChemRich NJ Pilot (confirm)", type="primary",
                     use_container_width=True, key="pilot_btn"):
            quote = round(result.cost_per_kg * batch_size * 1.12, 0)
            save_booking(get_session_id(), result.blend, result.cost_per_kg,
                         batch_size, quote, parsed.application_type)
            track("pilot_button_clicked",
                  {"cost_per_kg": result.cost_per_kg, "quote_usd": quote})
            st.balloons()
            st.success(f"✅ Booking submitted! Quote: **${quote:,.0f}** · 5 days · shehan@chemenova.com")
            if uemail:
                send_pilot_booking_confirmation(
                    customer_email=uemail,
                    customer_name=uname or "Valued Customer",
                    blend=result.blend,
                    cost_per_kg=result.cost_per_kg,
                    batch_kg=batch_size,
                    quote_usd=quote,
                    application=parsed.application_type,
                )

        # Save project
        project = {
            "timestamp":   datetime.now().strftime("%b %d %H:%M"),
            "input":       nl_input,
            "application": parsed.application_type,
            "blend":       result.blend,
            "cost":        result.cost_per_kg,
            "bio":         result.bio_pct,
            "perf":        result.perf_score,
            "eco_score":   eco.eco_score if eco else None,
            "eco_grade":   eco.grade if eco else None,
            "relaxed":     result.relaxed,
            "savings":     round((result.cost_per_kg * 1.28 - result.cost_per_kg) * 500, 0),
            "co2_kg":      round(500 * 0.75, 0),
            "parser":      parsed.parser_backend,
            "optimizer":   _opt_mode.lower(),
        }
        st.session_state.projects.append(project)
        save_project(project, get_session_id())

# ─────────────────────────────────────────────────────────────────────────────
# TAB 2 — CERTIFY (CertificationOracle™)
# ─────────────────────────────────────────────────────────────────────────────
with t_cert:
    st.subheader("🏆 CertificationOracle™")
    st.caption("Predicts pass probability for 8 major green chemistry standards before $2k–$80k of testing.")

    cert = st.session_state.last_cert
    if not cert:
        st.info("Run a formulation first.", icon="⚗")
    else:
        ca1, ca2, ca3, ca4 = st.columns(4)
        ca1.metric("Green Score",    f"{cert.overall_green_score:.0f}/100")
        ca2.metric("Natural Origin", f"{cert.natural_origin_index:.2f}")
        ca3.metric("Top Cert",       cert.top_certification[:20])
        ca4.metric("Quick Wins",     len(cert.quick_wins))

        if cert.quick_wins:
            st.success(f"✅ Quick wins: {', '.join(cert.quick_wins)}")
        if cert.recommended_certs:
            st.info(f"💡 Best ROI: {', '.join(cert.recommended_certs)}")

        st.divider()
        for cert_name, pred in cert.predictions.items():
            if "N/A" in pred.verdict:
                continue
            with st.expander(
                f"**{cert_name}** — {pred.verdict} ({pred.pass_probability*100:.0f}%)",
                expanded=pred.pass_probability > 0.65
            ):
                p1, p2, p3, p4 = st.columns(4)
                p1.metric("Pass Prob.",  f"{pred.pass_probability*100:.0f}%")
                p2.metric("Score",       f"{pred.score:.0f}/100")
                p3.metric("Cost",        pred.estimated_cost[:18])
                p4.metric("Timeline",    pred.estimated_timeline[:18])
                if pred.blocking_issues:
                    for b in pred.blocking_issues: st.error(b)
                if pred.gap_analysis:
                    for g in pred.gap_analysis: st.warning(g)
                if pred.strengths:
                    for s in pred.strengths: st.success(s)

        st.divider()
        cb_data = [{
            "Certification": c,
            "Pass Prob": f"{p.pass_probability*100:.0f}%",
            "Cost": p.estimated_cost,
            "Timeline": p.estimated_timeline,
            "Verdict": p.verdict,
        } for c, p in cert.predictions.items()
          if p.pass_probability > 0 and "N/A" not in p.verdict]
        if cb_data:
            st.dataframe(pd.DataFrame(cb_data), use_container_width=True, hide_index=True)

        # ── CERTIFICATION_PROFILES reference table (v1.5 feature) ──────────────
        with st.expander("📚 Certification Profiles Reference — All 8 Standards"):
            ref_rows = []
            for cert_key, profile in CERTIFICATION_PROFILES.items():
                ref_rows.append({
                    "Certification":  cert_key,
                    "Min Bio%":       profile.get("min_bio_based", "—"),
                    "Max Petrochem%": profile.get("max_petrochemical", "—"),
                    "Testing Cost":   profile.get("testing_cost", "—"),
                    "Timeline":       profile.get("timeline", "—"),
                    "Body":           profile.get("certifying_body", "—"),
                    "Scope":          str(profile.get("scope", "—"))[:50],
                })
            if ref_rows:
                st.dataframe(pd.DataFrame(ref_rows), use_container_width=True, hide_index=True)

# ─────────────────────────────────────────────────────────────────────────────
# TAB 3 — ECOMETRICS
# ─────────────────────────────────────────────────────────────────────────────
with t_eco:
    st.subheader("🌿 EcoMetrics™ Sustainability Scoring")
    eco = st.session_state.last_eco
    if not eco:
        st.info("Run a formulation first.", icon="🌿")
    else:
        import plotly.graph_objects as go
        c1, c2, c3, c4, c5, c6 = st.columns(6)
        c1.metric("EcoScore™", f"{eco.eco_score:.0f}/100")
        c2.metric("Grade",     eco.grade)
        c3.metric("Biodeg.",   f"{eco.biodegradability:.0f}")
        c4.metric("Carbon FP", f"{eco.carbon_footprint:.0f}")
        c5.metric("Ecotox",    f"{eco.ecotoxicity:.0f}")
        c6.metric("Renewable", f"{eco.renewability:.0f}")

        radar = ecometrics_radar_data(eco)
        cats  = radar["categories"]
        fig   = go.Figure()
        fig.add_trace(go.Scatterpolar(
            r=radar["intelliform"] + [radar["intelliform"][0]],
            theta=cats + [cats[0]], fill="toself", name="IntelliForm",
            line_color="#0D9488", fillcolor="rgba(13,148,136,0.25)"))
        fig.add_trace(go.Scatterpolar(
            r=radar["baseline"] + [radar["baseline"][0]],
            theta=cats + [cats[0]], fill="toself", name="Petrochemical Baseline",
            line_color="#ef4444", fillcolor="rgba(239,68,68,0.1)"))
        fig.update_layout(
            polar=dict(radialaxis=dict(range=[0,100], gridcolor="#1e3a5f"),
                       angularaxis=dict(gridcolor="#1e3a5f"), bgcolor="#0a1628"),
            paper_bgcolor="#050E1F", font_color="#f1f5f9", height=420,
            title="EcoMetrics™ vs Petrochemical Baseline",
            legend=dict(bgcolor="rgba(0,0,0,0)"))
        st.plotly_chart(fig, use_container_width=True)

        baselines = {"Biodegradability":52,"Carbon Footprint":38,"Ecotoxicity":41,
                     "Renewability":25,"Regulatory":60}
        score_map = {"Biodegradability":eco.biodegradability,"Carbon Footprint":eco.carbon_footprint,
                     "Ecotoxicity":eco.ecotoxicity,"Renewability":eco.renewability,"Regulatory":eco.regulatory}
        delta_rows = [{"Axis":k,"IntelliForm":f"{score_map[k]:.1f}","Baseline":v,
                       "Δ":f"{'▲' if (score_map[k]-v)>0 else '▼'} {abs(score_map[k]-v):.1f}"}
                      for k,v in baselines.items()]
        st.dataframe(pd.DataFrame(delta_rows), use_container_width=True, hide_index=True)
        st.caption("Makani S.S., ChemRxiv 2026 · DOI: 10.26434/chemrxiv.15000857 · NJIT Showcase 2025 · UIC Indigo 2025")

# ─────────────────────────────────────────────────────────────────────────────
# TAB 4 — REGULATORY
# ─────────────────────────────────────────────────────────────────────────────
with t_reg:
    st.subheader("📋 Regulatory Intelligence")
    v_reg = st.session_state.last_v_reg
    reg   = st.session_state.last_reg

    if v_reg and v_reg.vertical != "all":
        st.caption(f"**Framework:** {v_reg.framework}")
        _status_map = {"✅ Clear":"success","⚠️ Review Required":"warning","❌ Blocked":"error"}
        getattr(st, _status_map.get(v_reg.overall_status,"info"))(
            f"**{v_reg.overall_status}** — {v_reg.vertical.replace('_',' ').title()}")
        if v_reg.certifications:
            st.subheader("🏆 Eligible Certifications")
            for c in v_reg.certifications: st.success(c)
        if v_reg.flags:
            st.subheader("🚨 Blocking Issues")
            for f in v_reg.flags: st.error(f)
        if v_reg.warnings:
            for w in v_reg.warnings: st.warning(w)
        if v_reg.passes:
            for p in v_reg.passes[:5]: st.success(p)
        st.divider()

    if not reg:
        st.info("Run a formulation first.", icon="📋")
    else:
        _sm = {"✅ Clear":"success","⚠️ Review Required":"warning","❌ Blocked":"error"}
        getattr(st, _sm.get(reg.overall_status,"info"))(
            f"**{reg.overall_status}** · EU Ecolabel: {'✅' if reg.eu_ecolabel_eligible else '❌'} · "
            f"COSMOS: {'✅' if reg.cosmos_eligible else '❌'} · "
            f"EPA SC: {'✅' if reg.epa_safer_choice_eligible else '❌'}")
        st.subheader("Per-Ingredient Detail")
        st.dataframe(regulatory_table_df(reg.blend), use_container_width=True, hide_index=True)
        if reg.certification_pathways:
            st.subheader("🏆 Certification Pathways")
            for p in reg.certification_pathways: st.success(p)
        if reg.amber_flags:
            for f in reg.amber_flags: st.warning(f)
        if reg.red_flags:
            for f in reg.red_flags: st.error(f)

# ─────────────────────────────────────────────────────────────────────────────
# TAB 5 — PARETO FRONTIER
# ─────────────────────────────────────────────────────────────────────────────
with t_pareto:
    st.subheader("📈 Multi-Objective Pareto Frontier")
    pareto = st.session_state.last_pareto
    if not pareto:
        st.info("Select **Pareto** optimizer mode and run a formulation.", icon="📈")
    elif not pareto.success:
        st.error(pareto.error_msg)
    else:
        import plotly.express as px
        pa1, pa2, pa3, pa4 = st.columns(4)
        pa1.metric("Solutions",  pareto.n_solutions)
        pa2.metric("Backend",    pareto.backend.upper())
        pa3.metric("Best Cost",  f"${min(s.cost_per_kg for s in pareto.frontier):.2f}/kg")
        pa4.metric("Best Bio%",  f"{max(s.bio_pct for s in pareto.frontier):.1f}%")

        df_p = pareto_frontier_dataframe(pareto)
        if not df_p.empty:
            rec_id = pareto.recommended.solution_id if pareto.recommended else -1
            df_p["Type"] = df_p["ID"].apply(
                lambda x: "⭐ Recommended" if x == rec_id else "Frontier")
            import plotly.graph_objects as go
            fig3 = px.scatter_3d(df_p, x="Cost ($/kg)", y="Bio-based (%)",
                                  z="Perf Score", color="Type", height=520,
                                  color_discrete_map={"⭐ Recommended":"#D97706",
                                                      "Frontier":"#0D9488"},
                                  title="Pareto Frontier: Cost vs Bio% vs Performance")
            fig3.update_layout(paper_bgcolor="#050E1F", font_color="#f1f5f9",
                scene=dict(bgcolor="#0a1628",
                           xaxis=dict(gridcolor="#1e3a5f", color="#94a3b8"),
                           yaxis=dict(gridcolor="#1e3a5f", color="#94a3b8"),
                           zaxis=dict(gridcolor="#1e3a5f", color="#94a3b8")))
            st.plotly_chart(fig3, use_container_width=True)

        if pareto.recommended:
            rec = pareto.recommended
            st.markdown(
                f'<div class="cn-result-card"><b>⭐ TOPSIS Recommended</b> — '
                f'Cost: <b>${rec.cost_per_kg}/kg</b> · Bio: <b>{rec.bio_pct}%</b> · '
                f'Perf: <b>{rec.perf_score}/100</b><br>'
                f'{", ".join(f"<b>{k}</b> ({v}%)" for k,v in rec.blend.items())}</div>',
                unsafe_allow_html=True)

        st.dataframe(df_p.drop(columns=["Type"], errors="ignore"),
                     use_container_width=True, hide_index=True)
        st.download_button("📥 Download Frontier CSV", df_p.to_csv(index=False),
                           "IntelliForm_Pareto.csv", "text/csv")

# ─────────────────────────────────────────────────────────────────────────────
# TAB 6 — QSAR MODEL CARD
# ─────────────────────────────────────────────────────────────────────────────
with t_qsar:
    st.subheader("🔬 QSAR Engine + Model Card")
    st.caption("Makani S.S., ChemRxiv 2026 · DOI: 10.26434/chemrxiv.15000857 · NJIT Showcase 2025 · UIC Indigo 2025")

    mc2 = st.session_state.model_card
    if mc2:
        qi1, qi2, qi3, qi4 = st.columns(4)
        qi1.metric("N Ingredients", mc2.n_training)
        qi2.metric("DB Hash",       mc2.training_hash[:8] if mc2.training_hash else "—")
        qi3.metric("sklearn",       mc2.sklearn_version[:8] if mc2.sklearn_version else "—")
        qi4.metric("AL Rounds",     mc2.active_learning_rounds)
        st.divider()
        for target, bench in mc2.benchmarks.items():
            with st.expander(f"**{target}** — R²={bench['cv_r2']} · RMSE={bench['cv_rmse']}", expanded=True):
                bc1, bc2, bc3, bc4 = st.columns(4)
                bc1.metric("5-fold R²",  bench["cv_r2"])
                bc2.metric("CV RMSE",    f"{bench['cv_rmse']}")
                bc3.metric("Unit",       bench["unit"][:12] if bench.get("unit") else "—")
                bc4.metric("Algorithm",  "GBR")
                st.caption(f"**Model**: {bench['model']} · **Descriptors**: {bench['descriptor']}")

    st.divider()
    st.subheader("🔮 Live SMILES Prediction")
    test_smi = st.text_input("SMILES", value="CCCCCCCCCCCCOC1OC(CO)C(O)C(O)C1O")
    if st.button("Predict", type="primary"):
        qp = predict_properties(test_smi)
        qc1, qc2, qc3 = st.columns(3)
        qc1.metric("Biodegradability", f"{qp.biodegradability:.1f}%")
        qc2.metric("Ecotoxicity",      f"{qp.ecotoxicity:.1f}/10")
        qc3.metric("Performance",      f"{qp.performance:.1f}/100")
        st.caption(f"{'ML model' if qp.used_ml else 'Rule-based'} · Confidence: {qp.confidence}")
        for w in qp.warnings: st.warning(w)

    st.divider()
    st.subheader("📬 Active Learning — Submit Validated Data")
    al_smi = st.text_input("SMILES (validated compound)", key="al_smi")
    al_tgt = st.selectbox("Property", ["Biodegradability","Ecotoxicity","Performance"])
    al_val = st.number_input("Measured value", 0.0, 100.0, 90.0)
    if st.button("Submit Feedback"):
        if al_smi:
            qp2  = predict_properties(al_smi)
            pred = {"Biodegradability": qp2.biodegradability,
                    "Ecotoxicity": qp2.ecotoxicity,
                    "Performance": qp2.performance}[al_tgt]
            db_save_feedback(get_session_id(), al_smi, al_tgt, pred, al_val)
            st.success(submit_feedback(al_smi, al_tgt, al_val, ingredients_db))

# ─────────────────────────────────────────────────────────────────────────────
# TAB 7 — STABILITY & VISCOSITY
# ─────────────────────────────────────────────────────────────────────────────
with t_stab:
    st.subheader("🧪 Stability & Viscosity Prediction")
    stab = st.session_state.last_stability
    if not stab:
        st.info("Run a formulation first.", icon="🧪")
    else:
        s1, s2, s3, s4 = st.columns(4)
        s1.metric("Shelf Life",  stab.shelf_life_range)
        s2.metric("Viscosity",   stab.viscosity_range[:15] if stab.viscosity_range else "—")
        s3.metric("pH Range",    f"{stab.ph_min:.1f} – {stab.ph_max:.1f}")
        s4.metric("Stability",   stab.overall_rating[:15] if stab.overall_rating else "—")
        st.divider()
        sr1, sr2 = st.columns(2)
        with sr1:
            st.subheader("Stability Risks")
            for risk in stab.stability_risks: st.warning(risk)
        with sr2:
            st.subheader("Stability Boosters")
            for boost in stab.stability_boosters: st.success(boost)
        st.divider()
        st.subheader("Packaging Recommendation")
        st.info(stab.recommended_packaging)

# ─────────────────────────────────────────────────────────────────────────────
# TAB 8 — CARBON (Passport + CBAM + Credits)
# ─────────────────────────────────────────────────────────────────────────────
with t_carbon:
    st.subheader("🌍 Carbon Intelligence")
    ct1, ct2, ct3 = st.tabs(["🛂 Carbon Passport", "🇪🇺 CBAM Calculator", "💰 Carbon Credits"])

    # ── Carbon Passport ──────────────────────────────────────────────────────
    with ct1:
        st.caption("ISO 14067:2018 product carbon footprint · EU CBAM-ready JSON · blockchain audit hash.")
        result = st.session_state.last_result
        if not result:
            st.info("Run a formulation first.")
        else:
            with st.expander("⚙ Configure Passport", expanded=True):
                pp1, pp2 = st.columns(2)
                with pp1:
                    product_name = st.text_input("Product name", value="IntelliForm Formulation")
                    batch_id_in  = st.text_input("Batch ID", value=f"CR-{datetime.now().strftime('%Y%m%d')}-0001")
                with pp2:
                    process_in   = st.selectbox("Process", ["blending","emulsification","granulation","spray_drying","fermentation"])
                    grid_in      = st.selectbox("Grid region", ["US","EU","UK","China","India","Global"])
                    batch_kg_pp  = st.number_input("Batch size (kg)", value=500, min_value=1)

            if st.button("🛂 Generate Carbon Passport", type="primary", use_container_width=True):
                v_db2 = filter_db_by_vertical(ingredients_db, selected_vertical)
                passport = generate_carbon_passport(
                    blend=result.blend, db=v_db2,
                    product_name=product_name, batch_id=batch_id_in,
                    batch_kg=batch_kg_pp, manufacturing_process=process_in,
                    grid_region=grid_in)
                st.session_state.last_passport = passport

            passport = st.session_state.last_passport
            if passport:
                pg1, pg2, pg3, pg4, pg5 = st.columns(5)
                pg1.metric("Total PCF",    f"{passport.total_pcf:.3f} kgCO₂eq/kg")
                pg2.metric("Net PCF",      f"{passport.net_pcf:.3f}")
                pg3.metric("Grade",        passport.carbon_intensity_grade)
                pg4.metric("vs Industry",  f"{passport.vs_industry_average:.2f}×")
                pg5.metric("Avoided",      f"{passport.avoided_emissions:.3f}")
                st.caption(f"🔐 Audit hash: {passport.blockchain_hash[:32]}…")
                st.download_button("📥 Download CBAM-Ready JSON",
                    passport_to_json(passport),
                    f"carbon_passport_{passport.batch_id}.json", "application/json",
                    use_container_width=True)

    # ── CBAM Calculator ───────────────────────────────────────────────────────
    with ct2:
        st.markdown("#### EU CBAM Liability Calculator")
        st.caption("Regulation 2023/956 — mandatory from January 2026 for chemical importers into the EU.")
        cl1, cl2 = st.columns(2)
        with cl1:
            pcf_in  = st.number_input("Product PCF (kg CO₂e/kg)", value=2.5, min_value=0.01, step=0.01)
            vol_kg  = st.number_input("Annual import volume (kg)", value=10_000, min_value=100, step=1_000)
            ets_pr  = st.number_input("EU ETS price (€/tCO₂e)", value=CBAM_RATE_EUR_PER_TON, min_value=10.0, step=1.0)
        with cl2:
            total_co2t = pcf_in * vol_kg / 1000
            cbam_total = total_co2t * ets_pr
            cbam_per_kg= cbam_total / vol_kg
            saved_eur  = cbam_total * 0.40
            st.markdown(f"""
            <div class="cn-cbam-card">
            <div style="font-family:monospace;font-size:0.6rem;color:#D97706;letter-spacing:0.14em;text-transform:uppercase;margin-bottom:14px">CBAM Liability</div>
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:14px">
            <div><div style="color:#64748b;font-size:0.75rem">Total CO₂</div><div style="color:#f1f5f9;font-size:1.2rem;font-weight:700">{total_co2t:.2f} t CO₂e</div></div>
            <div><div style="color:#64748b;font-size:0.75rem">CBAM cost</div><div style="color:#D97706;font-size:1.2rem;font-weight:700">€{cbam_total:,.0f}</div></div>
            <div><div style="color:#64748b;font-size:0.75rem">Per kg</div><div style="color:#D97706;font-size:1.2rem;font-weight:700">€{cbam_per_kg:.4f}</div></div>
            <div><div style="color:#64748b;font-size:0.75rem">Saved vs unoptimized</div><div style="color:#14b8a8;font-size:1.2rem;font-weight:700">€{saved_eur:,.0f}</div></div>
            </div></div>""", unsafe_allow_html=True)

    # ── Carbon Credits ────────────────────────────────────────────────────────
    with ct3:
        carbon = st.session_state.last_carbon
        if not carbon:
            st.info("Run a formulation first.")
        else:
            cr1, cr2, cr3, cr4 = st.columns(4)
            cr1.metric("Green PCF",      f"{carbon.green_co2_per_kg:.2f} kgCO2/kg")
            cr2.metric("Baseline PCF",   f"{carbon.baseline_co2_per_kg:.1f}")
            cr3.metric("CO2 Displaced",  f"{carbon.co2_displaced_kg:.1f} kg/batch")
            cr4.metric("Credits/batch",  f"{carbon.credits_per_batch:.4f}")
            st.divider()
            st.markdown(f'<div class="cn-card">{carbon.summary}</div>', unsafe_allow_html=True)
            a1, a2 = st.columns(2)
            a1.metric("Floor ($15/t)",   f"${carbon.credit_value_low:.2f}")
            a2.metric("Mid ($45/t)",     f"${carbon.credit_value_mid:.2f}")
            st.metric("Annual Value (12 batches)", f"${carbon.annual_credit_value_mid:,.0f}")

# ─────────────────────────────────────────────────────────────────────────────
# TAB 9 — REFORMULATION LAB
# ─────────────────────────────────────────────────────────────────────────────
with t_refo:
    st.subheader("🔁 Reformulation Lab™")
    st.caption("Pilot batch failed? Enter test results → root cause diagnosis + ranked minimal-change fixes.")

    result_refo = st.session_state.last_result
    if not result_refo:
        st.info("Run a formulation first.", icon="🔁")
    else:
        with st.expander("📋 Enter Test Results", expanded=True):
            rf1, rf2 = st.columns(2)
            with rf1:
                failure_type_input = st.selectbox(
                    "What failed?",
                    list(FAILURE_TYPES.keys()),
                    format_func=lambda x: FAILURE_TYPES[x]["label"])
                failure_info = FAILURE_TYPES[failure_type_input]
                st.caption(f"ℹ️ {failure_info['description']}")
            with rf2:
                batch_ref = st.text_input("Batch reference",
                    value=f"CR-{datetime.now().strftime('%Y%m%d')}-0001")

            test_data = {}
            required_inputs = failure_info.get("inputs", [])
            ic = st.columns(min(len(required_inputs), 3)) if required_inputs else []
            for i, inp in enumerate(required_inputs):
                col = ic[i % len(ic)] if ic else st
                label = inp.replace("_", " ").title()
                if "ph" in inp and "target" not in inp:
                    test_data[inp] = col.number_input(label, value=7.0, min_value=0.0, max_value=14.0, step=0.1, key=f"rf_{inp}")
                elif "viscosity" in inp:
                    test_data[inp] = col.number_input(label, value=1000, min_value=0, step=100, key=f"rf_{inp}")
                elif "cost" in inp:
                    test_data[inp] = col.number_input(label, value=5.0, min_value=0.0, step=0.1, key=f"rf_{inp}")
                else:
                    test_data[inp] = col.text_input(label, key=f"rf_{inp}")

        if st.button("🔬 Diagnose & Get Fixes", type="primary", use_container_width=True):
            with st.spinner("Analysing failure…"):
                v_db_r = filter_db_by_vertical(ingredients_db, selected_vertical)
                refo   = run_reformulation_intelligence(
                    blend=result_refo.blend, db=v_db_r,
                    failure_type=failure_type_input,
                    test_data=test_data, batch_id=batch_ref)
                st.session_state.last_refo = refo

        refo = st.session_state.last_refo
        if refo:
            rca = refo.root_cause
            sev_colors = {"Minor":"#0D9488","Moderate":"#D97706","Severe":"#f97316","Critical":"#ef4444"}
            sev_color  = sev_colors.get(rca.severity, "#64748b")
            st.markdown(
                f'<div style="background:{sev_color}22;border-left:3px solid {sev_color};'
                f'border-radius:0 8px 8px 0;padding:12px 16px;margin:8px 0">'
                f'<strong style="color:{sev_color}">Root Cause ({rca.severity})</strong><br>'
                f'{rca.root_cause}</div>', unsafe_allow_html=True)
            for e in rca.evidence: st.caption(f"📌 {e}")

            sp = refo.predicted_success_probability
            sf1, sf2, sf3 = st.columns(3)
            sf1.metric("Fix Success Prob.", f"{sp*100:.0f}%")
            sf2.metric("Iterations Needed", refo.iterations_to_fix)
            sf3.metric("Confidence",        f"{rca.confidence*100:.0f}%")

            st.subheader(f"🔧 Ranked Fix Suggestions ({len(refo.suggestions)})")
            for sug in refo.suggestions:
                with st.expander(
                    f"{'⭐ ' if sug.rank==1 else ''}#{sug.rank}: {sug.action_type.upper()} "
                    f"{sug.ingredient} ({sug.current_pct:.1f}% → {sug.suggested_pct:.1f}%)",
                    expanded=sug.rank==1
                ):
                    sg1, sg2, sg3, sg4 = st.columns(4)
                    sg1.metric("Confidence", f"{sug.confidence*100:.0f}%")
                    sg2.metric("Risk",       sug.risk_level)
                    sg3.metric("Cost Δ",     f"${sug.cost_delta_per_kg:+.2f}/kg")
                    sg4.metric("Action",     sug.action_type.title())
                    st.info(sug.rationale)
                    st.caption(f"💡 {sug.implementation_notes}")

            if refo.do_not_change:
                st.success(f"✅ Keep unchanged: {', '.join(refo.do_not_change)}")
            st.info(f"📚 **Learning Note:** {refo.learning_note}")

# ─────────────────────────────────────────────────────────────────────────────
# TAB 10 — MEMORY NETWORK
# ─────────────────────────────────────────────────────────────────────────────
with t_mem:
    st.subheader("🧠 Formulation Memory Network™")
    st.caption("Living knowledge graph — learns from every run, rejection, and pilot batch outcome.")

    stats = memory.get_summary_stats()
    ms1, ms2, ms3, ms4, ms5, ms6 = st.columns(6)
    ms1.metric("Total Events",    stats["total_events"])
    ms2.metric("Formulations",    stats["total_formulations"])
    ms3.metric("Accepted",        stats["total_acceptances"])
    ms4.metric("Rejected",        stats["total_rejections"])
    ms5.metric("Acceptance Rate", f"{stats['acceptance_rate']*100:.0f}%" if stats["total_acceptances"]+stats["total_rejections"] > 0 else "—")
    ms6.metric("Negative Rules",  stats["negative_rules"])

    st.divider()
    st.subheader(f"💡 Insights — {selected_vertical.replace('_',' ').title()}")
    insights = memory.get_insights(selected_vertical)
    for ins in insights:
        icon  = {"preference":"✅","warning":"⚠️","opportunity":"💡","trend":"📈","info":"ℹ️"}.get(ins.insight_type,"•")
        color = {"preference":"success","warning":"warning","opportunity":"info","trend":"info","info":"info"}.get(ins.insight_type,"info")
        getattr(st, color)(
            f"{icon} **{ins.title}** — conf:{ins.confidence*100:.0f}% n={ins.evidence_count} — {ins.description[:80]}")

    neg_rules = memory.get_negative_rules(selected_vertical)
    if neg_rules:
        st.divider()
        st.subheader(f"🚫 Negative Knowledge DB — {len(neg_rules)} rules")
        nk_data = [{"Ingredient":r["ingredient"],"Reason":r["reason"][:60],
                    "Confidence":f"{r.get('confidence',0)*100:.0f}%","Evidence":r.get("evidence_count",1)}
                   for r in sorted(neg_rules, key=lambda x: -x.get("confidence",0))[:20]]
        st.dataframe(pd.DataFrame(nk_data), use_container_width=True, hide_index=True)

    st.divider()
    st.subheader("📥 Record Feedback")
    result_mem = st.session_state.last_result
    if result_mem:
        fb1, fb2 = st.columns(2)
        with fb1:
            if st.button("✅ Accept this blend", use_container_width=True, key="mem_accept"):
                memory.record("blend_accepted", result_mem.blend, selected_vertical, outcome="positive")
                st.success("✅ Acceptance recorded")
        with fb2:
            if st.button("❌ Reject this blend", use_container_width=True, key="mem_reject"):
                memory.record("blend_rejected", result_mem.blend, selected_vertical,
                              metadata={"reason":"user rejection"}, outcome="negative")
                st.error("❌ Rejection recorded — negative knowledge updated")
        pp1, pp2 = st.columns(2)
        with pp1:
            if st.button("✅ Pilot PASSED", use_container_width=True, key="mem_pass"):
                memory.record_pilot_outcome(result_mem.blend, selected_vertical, passed=True)
                st.success("🎉 Pilot pass recorded")
        with pp2:
            if st.button("❌ Pilot FAILED", use_container_width=True, key="mem_fail"):
                memory.record_pilot_outcome(result_mem.blend, selected_vertical, passed=False,
                                            failure_type="unspecified",
                                            failure_ingredients=list(result_mem.blend.keys())[:2])
                st.error("📋 Failure recorded — use Reformulation Lab to fix")

    kb_json = memory.export_knowledge_base()
    st.download_button("📥 Export Knowledge Base (JSON)", kb_json,
                       f"intelliform_memory_{datetime.now().strftime('%Y%m%d')}.json",
                       "application/json", use_container_width=True)

# ─────────────────────────────────────────────────────────────────────────────
# TAB 11 — PROPOSAL (PDF)
# ─────────────────────────────────────────────────────────────────────────────
with t_prop:
    st.subheader("📄 Proposal Generator")
    st.caption("Branded 4-page PDF proposal with EcoMetrics, regulatory table, pilot quote, and next steps.")

    if not st.session_state.projects:
        st.info("Run a formulation first.", icon="📄")
    else:
        latest  = st.session_state.projects[-1]
        eco_res = st.session_state.last_eco
        reg_res = st.session_state.last_reg
        fmt     = st.radio("Format", ["📄 PDF (branded)", "📝 Markdown"], horizontal=True)

        if "PDF" in fmt:
            if st.button("Generate PDF", type="primary", use_container_width=True):
                with st.spinner("Generating branded PDF…"):
                    try:
                        logo_bytes = customer_logo.read() if customer_logo else None
                        pdf_bytes = generate_proposal_pdf(
                            latest, eco_res, reg_res, ingredients_db,
                            customer_name=uname or None,
                            customer_company=customer_company or None,
                            logo_bytes=logo_bytes)
                        st.download_button(
                            "📥 Download PDF", data=pdf_bytes,
                            file_name=f"IntelliForm_Proposal_{datetime.now().strftime('%Y%m%d')}.pdf",
                            mime="application/pdf", use_container_width=True)
                        track("export_proposal", {"format": "pdf",
                                                   "white_label": logo_bytes is not None})
                        st.success("✅ PDF ready — click above to download.")
                    except Exception as e:
                        st.error(f"PDF failed: {e}")
        else:
            blend_lines = "\n".join(f"  - {k}: {v}%" for k,v in latest["blend"].items())
            eco_sec = f"\n## EcoMetrics\nEcoScore: {eco_res.eco_score:.1f}/100 · Grade: {eco_res.grade}\n" if eco_res else ""
            reg_sec = f"\n## Regulatory\n{reg_res.overall_status}\n" if reg_res else ""
            md = (f"# IntelliForm Proposal — {datetime.now().strftime('%B %d, %Y')}\n"
                  f"ChemeNova LLC x ChemRich Global · {latest['application'].replace('_',' ').title()}\n"
                  f"## Blend\n{blend_lines}\n## Metrics\n"
                  f"Cost: ${latest['cost']}/kg · Bio: {latest['bio']}% · Perf: {latest['perf']}/100\n"
                  f"{eco_sec}{reg_sec}\n"
                  f"*Makani S.S., ChemRxiv 2026 · NJIT & UIC · shehan@chemenova.com*")
            st.download_button("📥 Download Markdown", md, "IntelliForm_Proposal.md", "text/markdown",
                               use_container_width=True)
            with st.expander("Preview"): st.markdown(md)

# ─────────────────────────────────────────────────────────────────────────────
# TAB 12 — HISTORY & ROI
# ─────────────────────────────────────────────────────────────────────────────
with t_hist:
    st.subheader("📊 ROI & Formulation History")
    st.caption(f"Storage: {'Supabase (persistent)' if is_connected() else 'Session memory only'}")

    if not st.session_state.projects:
        st.info("Run a formulation first.", icon="📊")
    else:
        df = pd.DataFrame(st.session_state.projects)
        h1, h2, h3, h4 = st.columns(4)
        h1.metric("Runs",         len(df))
        h2.metric("Total Savings",f"${df['savings'].sum():,.0f}")
        h3.metric("CO2 Avoided",  f"{df['co2_kg'].sum():,.0f} kg")
        avg_eco = df["eco_score"].dropna().mean() if "eco_score" in df.columns else None
        h4.metric("Avg EcoScore", f"{avg_eco:.1f}" if avg_eco else "—")

        cols = [c for c in ["timestamp","application","cost","bio","perf","eco_score",
                            "eco_grade","optimizer","parser"] if c in df.columns]
        st.dataframe(df[cols], use_container_width=True)

        # Blend Comparison
        if len(st.session_state.blend_history) >= 2:
            st.divider()
            st.subheader("🔄 Blend Comparison")
            history  = st.session_state.blend_history
            labels   = [h["label"] for h in history]
            bc1, bc2 = st.columns(2)
            with bc1:
                sel_a = st.selectbox("Blend A", labels, index=0)
            with bc2:
                sel_b = st.selectbox("Blend B", labels, index=min(1, len(labels)-1))
            blend_a = next(h for h in history if h["label"] == sel_a)
            blend_b = next(h for h in history if h["label"] == sel_b)
            m1, m2, m3 = st.columns(3)
            m1.metric("Cost/kg",   f"${blend_b['cost']}",
                      f"{round(blend_b['cost']-blend_a['cost'],2):+.2f} vs A", delta_color="inverse")
            m2.metric("Bio-based", f"{blend_b['bio']}%",
                      f"{round(blend_b['bio']-blend_a['bio'],1):+.1f}% vs A")
            m3.metric("Perf",      f"{blend_b['perf']}/100",
                      f"{round(blend_b['perf']-blend_a['perf'],1):+.1f} vs A")
            import plotly.graph_objects as go
            all_ings = sorted(set(list(blend_a["blend"].keys()) + list(blend_b["blend"].keys())))
            vals_a   = [blend_a["blend"].get(i,0) for i in all_ings]
            vals_b   = [blend_b["blend"].get(i,0) for i in all_ings]
            fig_comp = go.Figure()
            fig_comp.add_trace(go.Bar(name="Blend A", x=all_ings, y=vals_a, marker_color="#0D9488"))
            fig_comp.add_trace(go.Bar(name="Blend B", x=all_ings, y=vals_b, marker_color="#D97706"))
            fig_comp.update_layout(barmode="group", paper_bgcolor="#050E1F",
                                   plot_bgcolor="#0a1628", font_color="#f1f5f9",
                                   title="Ingredient Comparison (%)", xaxis_tickangle=-30)
            st.plotly_chart(fig_comp, use_container_width=True)
            comp_df = pd.DataFrame({"Ingredient":all_ings,"Blend A %":vals_a,"Blend B %":vals_b})
            st.download_button("📥 Download Comparison CSV", comp_df.to_csv(index=False),
                               "IntelliForm_Comparison.csv", "text/csv")

    if not is_connected():
        st.divider()
        with st.expander("🗄 Set Up Supabase — View migration SQL"):
            st.code(MIGRATION_SQL, language="sql")

# ─────────────────────────────────────────────────────────────────────────────
# TAB 13 — PRO
# ─────────────────────────────────────────────────────────────────────────────
with t_pro:
    st.subheader("⚡ Pro Access & Pilot Manufacturing")
    pr1, pr2 = st.columns([3, 2])
    with pr1:
        st.markdown(
            '<div class="cn-result-card">'
            '<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.68rem;'
            'letter-spacing:0.15em;text-transform:uppercase;color:#0D9488;margin-bottom:10px">'
            'IntelliForm Pro — $99/month</div>'
            '<div style="color:#f1f5f9;font-size:1.2rem;font-weight:500;margin-bottom:20px">'
            'Claude 3.5 Sonnet · GPT-4o · Custom databases</div>'
            '</div>', unsafe_allow_html=True)
        for feat, desc in [
            ("Claude 3.5 Sonnet + GPT-4o", "Full AI council for complex pharma/specialty formulations"),
            ("Custom ingredient DB upload", "Upload proprietary ingredient list"),
            ("White-label PDF reports",    "Branded proposals for clients"),
            ("Priority support",           "24h response from ChemeNova team"),
        ]:
            st.markdown(f"✅ **{feat}** — {desc}")

        st.divider()
        with st.form("pro_form", clear_on_submit=True):
            st.markdown("**Request Pro access**")
            p_email   = st.text_input("Work email",  placeholder="you@company.com")
            p_company = st.text_input("Company",     placeholder="Acme Chemical Corp.")
            p_use     = st.text_area("Use case",     placeholder="e.g. COSMOS-certified personal care for EU", height=80)
            sub       = st.form_submit_button("Request Pro Access →", use_container_width=True)
            if sub and p_email:
                track("pro_inquiry", {"email": p_email, "company": p_company})
                st.success(f"✅ Request received! We'll reach out to **{p_email}** within 24h.")

    with pr2:
        st.markdown(
            '<div class="cn-card">'
            '<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.6rem;color:#0D9488;'
            'letter-spacing:0.14em;text-transform:uppercase;margin-bottom:8px">Enterprise</div>'
            '<div style="color:#f1f5f9;font-size:1.4rem;font-weight:400;margin-bottom:8px">Custom</div>'
            '<div style="color:#94a3b8;font-size:0.85rem">REST API · SSO/SAML · ERP/LIMS · '
            'on-premise · 21 CFR Part 11 audit log · unlimited seats</div>'
            '</div>', unsafe_allow_html=True)
        st.markdown(
            '<a href="mailto:shehan@chemenova.com?subject=IntelliForm%20Enterprise%20Demo" '
            'style="display:block;background:linear-gradient(135deg,#D97706,#b45309);'
            'color:#fff;padding:12px;border-radius:8px;font-family:\'IBM Plex Mono\','
            'monospace;font-size:0.76rem;font-weight:700;text-transform:uppercase;'
            'text-align:center;text-decoration:none;letter-spacing:0.08em">'
            'Contact Enterprise Sales →</a>', unsafe_allow_html=True)

    st.divider()
    st.markdown("### 🏭 Pilot Batch Manufacturing")
    st.markdown(_pilot_cta_html(st.session_state.last_result), unsafe_allow_html=True)

    st.divider()
    st.markdown("### 🔗 Supply Chain Resilience")
    st.caption("Find drop-in analogs when an ingredient is out of stock or delisted. "
               "Powered by SMILES similarity search across the full ingredient DB.")
    sc1, sc2 = st.columns([3, 1])
    with sc1:
        sc_ingredient = st.selectbox(
            "Ingredient at risk",
            options=sorted(ingredients_db["Ingredient"].tolist()) if not ingredients_db.empty else [],
            key="sc_ingredient",
        )
    with sc2:
        sc_top_n = st.number_input("Top N analogs", min_value=1, max_value=10, value=5, key="sc_top_n")
    if st.button("🔍 Find Analogs", use_container_width=True, key="sc_run"):
        sc_rows = ingredients_db[ingredients_db["Ingredient"] == sc_ingredient]
        if sc_rows.empty:
            st.warning("Ingredient not found in DB.")
        else:
            sc_row = sc_rows.iloc[0]
            sc_smiles = sc_row.get("SMILES", "")
            if not sc_smiles:
                st.warning("No SMILES available for this ingredient — cannot compute similarity.")
            else:
                from modules.chem_utils import tanimoto_similar
                analogs = tanimoto_similar(sc_smiles, ingredients_db,
                                           exclude=sc_ingredient, top_n=int(sc_top_n))
                if not analogs:
                    st.info("No close analogs found in the current DB.")
                else:
                    st.success(f"Found {len(analogs)} analogs for **{sc_ingredient}**")
                    analog_df = pd.DataFrame(analogs)
                    st.dataframe(analog_df, use_container_width=True, hide_index=True)
                    track("supply_chain_search", {"ingredient": sc_ingredient, "n_found": len(analogs)})

    st.divider()
    st.markdown("### 📄 Validation & IP")
    vc = st.columns(4)
    vc[0].markdown("**ChemRxiv Preprint**\nDOI: 10.26434/chemrxiv.15000857 · 2026")
    vc[1].markdown("**USPTO Provisional**\nTwo-stage NSGA-III surrogate pipeline filed")
    vc[2].markdown("**Academic Partners**\nNJIT Showcase 2025 · UIC Indigo 2025")
    vc[3].markdown("**MIT License**\ngithub.com/Cheme-Nova/IntelliForm")

# ─────────────────────────────────────────────────────────────────────────────
# TAB 14 — PHARMA DEEP DIVE (only when pharmaceutical vertical selected)
# ─────────────────────────────────────────────────────────────────────────────
if t_pharma:
    with t_pharma:
        st.subheader("💊 Pharmaceutical Deep Dive")
        st.caption("ICH Q8/Q9/Q10 · BCS Classification · API-excipient compatibility · "
                   "ICH Q1A stability zones · Dosage form engineering · Regulatory pathway guidance.")

        result_ph = st.session_state.last_result
        if not result_ph:
            st.info("Run a pharmaceutical formulation first using the Pharmaceutical vertical.")
        else:
            with st.expander("⚙ Configure Deep Dive", expanded=True):
                pd1, pd2, pd3 = st.columns(3)
                with pd1:
                    bcs_input = st.selectbox("BCS Classification",
                        ["I — High Sol / High Perm","II — Low Sol / High Perm",
                         "III — High Sol / Low Perm","IV — Low Sol / Low Perm"], index=0)
                    bcs_class = bcs_input.split(" ")[0]
                with pd2:
                    dosage_form_input = st.selectbox("Dosage Form",
                        list(DOSAGE_FORMS.keys()),
                        format_func=lambda x: DOSAGE_FORMS[x].form, index=0)
                with pd3:
                    is_generic    = st.checkbox("Generic (ANDA) pathway")
                    is_pediatric  = st.checkbox("Pediatric formulation")
                    target_markets= st.multiselect("Target markets",
                        ["USA","EU","India","Japan","Southeast Asia"], default=["USA","EU"])

            if st.button("🔬 Run Pharma Deep Dive", type="primary", use_container_width=True):
                with st.spinner("Running ICH analysis…"):
                    v_db_ph = filter_db_by_vertical(ingredients_db, "pharmaceutical")
                    pharma_result = run_pharma_deep_dive(
                        blend=result_ph.blend, db=v_db_ph, bcs_class=bcs_class,
                        dosage_form=dosage_form_input, target_markets=target_markets,
                        is_generic=is_generic, is_pediatric=is_pediatric)
                    st.session_state.pharma_result = pharma_result

            pharma = st.session_state.pharma_result
            if pharma:
                st.divider()
                st.subheader(f"📊 BCS Class {pharma.bcs_class} — "
                             f"{pharma.bcs_profile.solubility} Sol / {pharma.bcs_profile.permeability} Perm")
                bcs1, bcs2 = st.columns([2,1])
                with bcs1:
                    st.info(pharma.bcs_profile.formulation_strategy)
                with bcs2:
                    st.metric("Bioavailability Risk", pharma.bcs_profile.bioavailability_risk)
                    st.metric("ICM Score", f"{pharma.icm_score:.0f}/100")
                st.caption("Enabling excipients: " + " · ".join(pharma.bcs_profile.enabling_excipients[:4]))

                st.divider()
                st.subheader(f"💊 {pharma.recommended_dosage_form}")
                df_prof = pharma.dosage_form_profile
                if df_prof:
                    st.caption(df_prof.description)
                    dfd1, dfd2 = st.columns(2)
                    with dfd1:
                        st.write("**Typical Composition**")
                        for func, (lo, hi) in df_prof.typical_composition.items():
                            st.caption(f"• {func}: {lo:.0f}–{hi:.0f}%")
                    with dfd2:
                        st.write("**Critical Quality Attributes**")
                        for cqa in df_prof.critical_quality_attributes:
                            st.caption(f"• {cqa}")
                    st.info(f"**Process:** {df_prof.manufacturing_process}")

                st.divider()
                st.subheader(f"🏭 {pharma.manufacturing_route}")
                st.info(pharma.manufacturing_rationale)

                st.divider()
                st.subheader("⚗ API-Excipient Compatibility")
                st.caption(pharma.compatibility_summary)
                if pharma.compatibility_results:
                    compat_df = pd.DataFrame([{
                        "Ingredient": c.ingredient,
                        "Severity":   c.severity,
                        "Interaction":c.interaction_type,
                        "Mechanism":  c.mechanism[:70],
                        "Mitigation": c.mitigation[:70],
                    } for c in pharma.compatibility_results])
                    st.dataframe(compat_df, use_container_width=True, hide_index=True)
                    for c in pharma.compatibility_results:
                        if c.severity == "Severe":
                            st.error(f"🚨 {c.ingredient}: {c.mechanism}")
                        elif c.severity == "Moderate":
                            st.warning(f"⚠ {c.ingredient}: {c.mechanism[:100]}")
                else:
                    st.success("✅ No significant incompatibilities detected.")

                st.divider()
                st.subheader(f"🌡 ICH Q1A — Zone {pharma.stability_zone}")
                if pharma.stability_profile:
                    sp = pharma.stability_profile
                    sp1, sp2, sp3 = st.columns(3)
                    sp1.metric("Long-term",   sp.long_term.split(" for")[0])
                    sp2.metric("Accelerated", sp.accelerated.split(" for")[0])
                    sp3.metric("Regions",     sp.regions[0] if sp.regions else "—")
                    st.info(f"**Packaging:** {pharma.packaging_recommendation}")
                for concern in pharma.stability_concerns:
                    st.warning(concern)

                st.divider()
                risk_col, rec_col = st.columns(2)
                with risk_col:
                    st.subheader("⚠ Development Risks")
                    for risk in pharma.development_risks: st.error(risk)
                    if not pharma.development_risks: st.success("Low risk profile")
                with rec_col:
                    st.subheader("✅ Recommendations")
                    for rec in pharma.development_recommendations[:5]: st.success(rec)

                # ── REGULATORY_PATHWAYS reference table (v1.5 feature) ─────────
                st.divider()
                with st.expander("📋 Regulatory Pathway Guide — USA / EU / Japan"):
                    rp_rows = []
                    for pathway_key, pathway in REGULATORY_PATHWAYS.items():
                        rp_rows.append({
                            "Pathway":    pathway_key,
                            "Agency":     pathway.get("agency", "—"),
                            "Timeline":   pathway.get("timeline", "—"),
                            "Cost Est.":  pathway.get("cost_estimate", "—"),
                            "Key Stds":   ", ".join(pathway.get("key_standards", []))[:60],
                            "Notes":      str(pathway.get("notes", "—"))[:80],
                        })
                    if rp_rows:
                        st.dataframe(pd.DataFrame(rp_rows), use_container_width=True, hide_index=True)

                # ── ICH_STABILITY_ZONES reference table (v1.5 feature) ─────────
                with st.expander("🌡 ICH Q1A Stability Zones — All 5 Zones"):
                    sz_rows = []
                    for zone_key, zone in ICH_STABILITY_ZONES.items():
                        sz_rows.append({
                            "Zone":        zone_key,
                            "Long-term":   zone.get("long_term", "—"),
                            "Accelerated": zone.get("accelerated", "—"),
                            "Regions":     ", ".join(zone.get("regions", []))[:60],
                            "Notes":       str(zone.get("notes", "—"))[:80],
                        })
                    if sz_rows:
                        st.dataframe(pd.DataFrame(sz_rows), use_container_width=True, hide_index=True)

# ── FOOTER ────────────────────────────────────────────────────────────────────
st.markdown("""
<hr style="border-color:rgba(30,58,95,0.5)">
<div style="text-align:center;font-family:'IBM Plex Mono',monospace;font-size:0.62rem;
     color:#64748b;padding:16px 0">
IntelliForm v2.1 ·
<a href="https://chemenova.com" style="color:#0D9488;text-decoration:none">ChemeNova LLC</a>
 ×
<a href="https://chemrichgroup.com" style="color:#0D9488;text-decoration:none">ChemRich Global</a>
 ·
<a href="https://github.com/Cheme-Nova/IntelliForm" style="color:#0D9488;text-decoration:none">MIT License</a>
 · DOI:
<a href="https://doi.org/10.26434/chemrxiv.15000857" style="color:#0D9488;text-decoration:none">10.26434/chemrxiv.15000857</a>
 · NJIT Showcase 2025 · UIC Indigo 2025
 · <a href="mailto:shehan@chemenova.com" style="color:#0D9488;text-decoration:none">shehan@chemenova.com</a>
</div>
""", unsafe_allow_html=True)
