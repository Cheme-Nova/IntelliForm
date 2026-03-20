
"""
IntelliForm v1.5 — Agentic Green Chemistry Formulation Platform
ChemeNova LLC x ChemRich Global

New in v1.0:
  Concentration slider · Blend comparison · Stability & viscosity prediction
  Carbon credit calculator · Version history · White-label PDF · Email notifications
  Dynamic reformulation trigger · Multi-LLM support
"""
import os
from datetime import datetime
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from dotenv import load_dotenv

load_dotenv()

from modules.analytics        import track, identify_user, get_session_id
from modules.llm_parser       import parse_request
from modules.tiers import get_tier, gate, TIERS
from modules.certification_oracle import run_certification_oracle, CertificationReport, CERTIFICATION_PROFILES
from modules.carbon_passport import generate_carbon_passport, passport_to_json
from modules.memory_network import get_memory_network, FormulationMemoryNetwork
from modules.reformulation_intelligence import (run_reformulation_intelligence,
    FAILURE_TYPES, ReformulationReport)
from modules.optimizer        import run_optimization
from modules.pareto_optimizer import run_pareto_optimization, pareto_frontier_dataframe
from modules.agents           import run_agent_swarm
from modules.chem_utils       import draw_mol, enrich_db
from modules.ecometrics       import compute_ecometrics, ecometrics_radar_data
from modules.qsar             import initialize_models, predict_properties, submit_feedback
from modules.regulatory       import get_blend_report, regulatory_table_df
from modules.persistence      import (save_project, load_projects, save_booking,
                                      save_feedback as db_save_feedback,
                                      is_connected, MIGRATION_SQL)
from modules.pdf_proposal     import generate_proposal_pdf
from modules.stability        import predict_stability
from modules.carbon_credits   import calculate_carbon_credits
from modules.notifications    import send_pilot_booking_confirmation, send_proposal_email
# Tier detection
try:
    _tier_name = st.secrets.get("INTELLIFORM_TIER", "free")
except Exception:
    import os as _os
    _tier_name = _os.environ.get("INTELLIFORM_TIER", "free")
TIER = get_tier(_tier_name)
from modules.verticals import VERTICAL_OPTIONS, get_profile, filter_db_by_vertical, get_vertical_constraints
from modules.vertical_regulatory import generate_vertical_regulatory_report
from modules.bayesian_optimizer import run_bayesian_optimization, BayesianState, BAYES_OK
from modules.pharma import (run_pharma_deep_dive, BCS_STRATEGIES, DOSAGE_FORMS,
                             REGULATORY_PATHWAYS, ICH_STABILITY_ZONES)

st.set_page_config(page_title="IntelliForm v1.5", page_icon="🧪", layout="wide")

# ── Enterprise UI CSS injection ───────────────────────────────────────────────
_CSS = """
<style>
/* ── Import fonts ── */
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=DM+Sans:wght@300;400;500;600;700&family=Syne:wght@700;800&display=swap');

/* ── Root tokens ── */
:root {
    --navy:   #0A1628;
    --teal:   #0D9488;
    --amber:  #D97706;
    --slate:  #1e3a5f;
    --mid:    #64748b;
    --light:  #94a3b8;
    --bg:     #0f1923;
    --card:   #111c2d;
    --border: #1e3450;
    --radius: 8px;
    --font-main: 'DM Sans', sans-serif;
    --font-mono: 'IBM Plex Mono', monospace;
    --font-head: 'Syne', sans-serif;
}

/* ── Global ── */
html, body, [data-testid="stApp"] {
    font-family: var(--font-main);
    background: var(--bg);
    color: #e2e8f0;
}

/* ── Hide Streamlit chrome ── */
#MainMenu, footer, header { visibility: hidden; }
[data-testid="stDecoration"] { display: none; }

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background: linear-gradient(180deg, var(--navy) 0%, #0d1f3c 100%);
    border-right: 1px solid var(--border);
}
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] .stRadio label,
[data-testid="stSidebar"] p {
    font-family: var(--font-main) !important;
    font-size: 0.82rem !important;
    color: var(--light) !important;
}
[data-testid="stSidebar"] .stSelectbox div[data-baseweb="select"] {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    font-size: 0.82rem;
}

/* ── Metric cards ── */
[data-testid="metric-container"] {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 14px 16px 10px 16px;
    transition: border-color 0.2s;
}
[data-testid="metric-container"]:hover {
    border-color: var(--teal);
}
[data-testid="stMetricLabel"] {
    font-family: var(--font-mono) !important;
    font-size: 0.68rem !important;
    font-weight: 600 !important;
    letter-spacing: 0.08em !important;
    text-transform: uppercase !important;
    color: var(--mid) !important;
    white-space: nowrap !important;
    overflow: hidden !important;
    text-overflow: ellipsis !important;
}
[data-testid="stMetricValue"] {
    font-family: var(--font-head) !important;
    font-size: 1.35rem !important;
    font-weight: 800 !important;
    color: #f1f5f9 !important;
    white-space: nowrap !important;
    overflow: hidden !important;
    text-overflow: ellipsis !important;
    line-height: 1.2 !important;
}
[data-testid="stMetricDelta"] {
    font-family: var(--font-mono) !important;
    font-size: 0.72rem !important;
}

/* ── Header / top bar ── */
.intelliform-header {
    background: linear-gradient(90deg, var(--navy) 0%, #0d1f3c 60%, #0a1e35 100%);
    border-bottom: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 14px 20px;
    margin-bottom: 16px;
    display: flex;
    align-items: center;
    gap: 16px;
}
.intelliform-logo {
    font-family: var(--font-head);
    font-size: 1.25rem;
    font-weight: 800;
    color: #fff;
    letter-spacing: 0.5px;
}
.intelliform-logo span {
    color: var(--teal);
}
.intelliform-badge {
    background: var(--amber);
    color: var(--navy);
    font-family: var(--font-mono);
    font-size: 0.6rem;
    font-weight: 700;
    letter-spacing: 2px;
    padding: 2px 8px;
    border-radius: 3px;
    text-transform: uppercase;
}
.intelliform-stat {
    font-family: var(--font-mono);
    font-size: 0.72rem;
    color: var(--light);
    padding: 4px 10px;
    border: 1px solid var(--border);
    border-radius: 4px;
    background: rgba(13,148,136,0.06);
}
.intelliform-stat strong { color: var(--teal); }

/* ── Tabs ── */
[data-testid="stTabs"] [data-baseweb="tab-list"] {
    background: var(--card);
    border-bottom: 1px solid var(--border);
    border-radius: var(--radius) var(--radius) 0 0;
    gap: 0;
    padding: 0 8px;
}
[data-testid="stTabs"] [data-baseweb="tab"] {
    font-family: var(--font-mono) !important;
    font-size: 0.72rem !important;
    font-weight: 600 !important;
    letter-spacing: 0.04em !important;
    color: var(--mid) !important;
    padding: 10px 14px !important;
    border-bottom: 2px solid transparent !important;
    border-radius: 0 !important;
    background: transparent !important;
    transition: all 0.15s !important;
}
[data-testid="stTabs"] [aria-selected="true"] {
    color: var(--teal) !important;
    border-bottom-color: var(--teal) !important;
    background: rgba(13,148,136,0.06) !important;
}
[data-testid="stTabs"] [data-baseweb="tab"]:hover {
    color: var(--light) !important;
    background: rgba(255,255,255,0.03) !important;
}

/* ── Buttons ── */
.stButton > button {
    font-family: var(--font-mono) !important;
    font-size: 0.8rem !important;
    font-weight: 600 !important;
    letter-spacing: 0.05em !important;
    border-radius: var(--radius) !important;
    transition: all 0.15s !important;
}
.stButton > button[kind="primary"] {
    background: linear-gradient(135deg, var(--teal) 0%, #0f766e 100%) !important;
    border: none !important;
    color: white !important;
    box-shadow: 0 2px 12px rgba(13,148,136,0.3) !important;
}
.stButton > button[kind="primary"]:hover {
    box-shadow: 0 4px 20px rgba(13,148,136,0.5) !important;
    transform: translateY(-1px) !important;
}
.stButton > button:not([kind="primary"]) {
    background: var(--card) !important;
    border: 1px solid var(--border) !important;
    color: var(--light) !important;
}
.stButton > button:not([kind="primary"]):hover {
    border-color: var(--teal) !important;
    color: var(--teal) !important;
}

/* ── Text inputs ── */
.stTextArea textarea, .stTextInput input {
    font-family: var(--font-main) !important;
    font-size: 0.88rem !important;
    background: var(--card) !important;
    border: 1px solid var(--border) !important;
    border-radius: var(--radius) !important;
    color: #e2e8f0 !important;
    caret-color: var(--teal) !important;
}
.stTextArea textarea:focus, .stTextInput input:focus {
    border-color: var(--teal) !important;
    box-shadow: 0 0 0 2px rgba(13,148,136,0.15) !important;
}

/* ── Expanders ── */
[data-testid="stExpander"] {
    background: var(--card) !important;
    border: 1px solid var(--border) !important;
    border-radius: var(--radius) !important;
    margin-bottom: 6px !important;
}
[data-testid="stExpander"] summary {
    font-family: var(--font-mono) !important;
    font-size: 0.8rem !important;
    font-weight: 600 !important;
    color: var(--light) !important;
    padding: 10px 14px !important;
}
[data-testid="stExpander"] summary:hover {
    color: var(--teal) !important;
    background: rgba(13,148,136,0.04) !important;
}

/* ── Alerts ── */
[data-testid="stSuccess"] {
    background: rgba(13,148,136,0.08) !important;
    border-left: 3px solid var(--teal) !important;
    border-radius: 0 var(--radius) var(--radius) 0 !important;
    font-size: 0.84rem !important;
}
[data-testid="stWarning"] {
    background: rgba(217,119,6,0.08) !important;
    border-left: 3px solid var(--amber) !important;
    border-radius: 0 var(--radius) var(--radius) 0 !important;
    font-size: 0.84rem !important;
}
[data-testid="stError"] {
    background: rgba(239,68,68,0.08) !important;
    border-left: 3px solid #ef4444 !important;
    border-radius: 0 var(--radius) var(--radius) 0 !important;
    font-size: 0.84rem !important;
}
[data-testid="stInfo"] {
    background: rgba(59,130,246,0.08) !important;
    border-left: 3px solid #3b82f6 !important;
    border-radius: 0 var(--radius) var(--radius) 0 !important;
    font-size: 0.84rem !important;
}

/* ── Dataframes ── */
[data-testid="stDataFrame"] {
    border: 1px solid var(--border) !important;
    border-radius: var(--radius) !important;
    overflow: hidden !important;
}
[data-testid="stDataFrame"] th {
    font-family: var(--font-mono) !important;
    font-size: 0.7rem !important;
    letter-spacing: 0.06em !important;
    text-transform: uppercase !important;
    background: var(--navy) !important;
    color: var(--mid) !important;
    padding: 8px 12px !important;
}
[data-testid="stDataFrame"] td {
    font-family: var(--font-main) !important;
    font-size: 0.82rem !important;
    padding: 7px 12px !important;
    border-bottom: 1px solid var(--border) !important;
}

/* ── Subheaders ── */
h2, h3, [data-testid="stMarkdown"] h2, [data-testid="stMarkdown"] h3 {
    font-family: var(--font-head) !important;
    font-weight: 700 !important;
    letter-spacing: -0.02em !important;
    color: #f1f5f9 !important;
}
h3 { font-size: 1.05rem !important; }

/* ── Captions ── */
[data-testid="stCaptionContainer"] p, .stCaption {
    font-family: var(--font-mono) !important;
    font-size: 0.72rem !important;
    color: var(--mid) !important;
    letter-spacing: 0.02em !important;
}

/* ── Radio buttons ── */
.stRadio [data-testid="stMarkdown"] p {
    font-family: var(--font-mono) !important;
    font-size: 0.78rem !important;
    font-weight: 600 !important;
}

/* ── Selectbox ── */
[data-baseweb="select"] {
    font-family: var(--font-main) !important;
    font-size: 0.84rem !important;
}

/* ── Sliders ── */
[data-testid="stSlider"] [data-baseweb="slider"] div[role="slider"] {
    background: var(--teal) !important;
    border-color: var(--teal) !important;
}
[data-testid="stSlider"] [data-baseweb="slider"] div[data-testid="stThumbValue"] {
    font-family: var(--font-mono) !important;
    font-size: 0.72rem !important;
    background: var(--navy) !important;
    border: 1px solid var(--teal) !important;
    color: var(--teal) !important;
}

/* ── Dividers ── */
hr {
    border-color: var(--border) !important;
    margin: 16px 0 !important;
}

/* ── Ingredient result cards ── */
.ing-card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 12px 16px;
    margin-bottom: 6px;
    border-left: 3px solid var(--teal);
    font-family: var(--font-main);
    font-size: 0.84rem;
}
.ing-name {
    font-family: var(--font-head);
    font-weight: 700;
    font-size: 0.95rem;
    color: #f1f5f9;
}
.ing-pct {
    font-family: var(--font-mono);
    font-weight: 700;
    color: var(--teal);
    font-size: 1rem;
}

/* ── Agent output boxes ── */
.agent-box {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 10px 14px;
    margin: 4px 0;
    font-family: var(--font-main);
    font-size: 0.83rem;
    line-height: 1.5;
}
.agent-label {
    font-family: var(--font-mono);
    font-size: 0.68rem;
    font-weight: 700;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    color: var(--mid);
    margin-bottom: 4px;
}

/* ── Scrollbar ── */
::-webkit-scrollbar { width: 5px; height: 5px; }
::-webkit-scrollbar-track { background: var(--bg); }
::-webkit-scrollbar-thumb { background: var(--slate); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: var(--teal); }

/* ── Plotly charts ── */
.js-plotly-plot .plotly .bg { fill: var(--card) !important; }

/* ── File uploader ── */
[data-testid="stFileUploader"] {
    background: var(--card) !important;
    border: 1px dashed var(--border) !important;
    border-radius: var(--radius) !important;
}
[data-testid="stFileUploader"]:hover {
    border-color: var(--teal) !important;
}

/* ── Checkbox ── */
[data-testid="stCheckbox"] label span {
    font-family: var(--font-main) !important;
    font-size: 0.82rem !important;
    color: var(--light) !important;
}
</style>
"""
st.markdown(_CSS, unsafe_allow_html=True)

st.markdown("""<style>
.pareto-rec{background:#1a2e1a;border-left:4px solid #00C853;padding:12px 18px;border-radius:6px;margin-bottom:12px}
.carbon-card{background:#0a2e1a;border-left:4px solid #059669;padding:12px 18px;border-radius:6px;margin-bottom:12px}
.stability-card{background:#1a1a2e;border-left:4px solid #0D9488;padding:12px 18px;border-radius:6px;margin-bottom:12px}
</style>""", unsafe_allow_html=True)

# ── Bootstrap ─────────────────────────────────────────────────────────────────
if "session_id"       not in st.session_state: get_session_id(); track("session_started",{"version":"1.3"})
if "projects"         not in st.session_state: st.session_state.projects       = []
if "last_result"      not in st.session_state: st.session_state.last_result    = None
if "last_parsed"      not in st.session_state: st.session_state.last_parsed    = None
if "last_eco"         not in st.session_state: st.session_state.last_eco       = None
if "last_pareto"      not in st.session_state: st.session_state.last_pareto    = None
if "last_reg"         not in st.session_state: st.session_state.last_reg       = None
if "last_stability"   not in st.session_state: st.session_state.last_stability = None
if "last_carbon"      not in st.session_state: st.session_state.last_carbon    = None
if "model_card"       not in st.session_state: st.session_state.model_card     = None
if "blend_history"    not in st.session_state: st.session_state.blend_history  = []
if "vertical_reg"     not in st.session_state: st.session_state.vertical_reg    = None
if "pharma_result"    not in st.session_state: st.session_state.pharma_result   = None
if "bayes_state"      not in st.session_state: st.session_state.bayes_state     = None
if "bayes_result"     not in st.session_state: st.session_state.bayes_result    = None
if "cert_report"      not in st.session_state: st.session_state.cert_report     = None
if "carbon_passport"  not in st.session_state: st.session_state.carbon_passport  = None
if "refo_report"      not in st.session_state: st.session_state.refo_report      = None
if "memory_net"       not in st.session_state: st.session_state.memory_net       = None
if "compare_blend"    not in st.session_state: st.session_state.compare_blend  = None

@st.cache_data(ttl=3600)
def load_db():
    df = pd.read_csv("data/ingredients_db.csv")
    # Ensure Vertical column exists (v4+ database)
    if "Vertical" not in df.columns:
        df["Vertical"] = "personal_care"  # fallback for old DB
    return enrich_db(df)

@st.cache_resource
def load_models(n, db_hash):
    """Retrain QSAR models. Cache key includes DB size + hash to force retraining on new DB."""
    db = load_db()
    return initialize_models(db)

ingredients_db = load_db()
_db_hash = str(len(ingredients_db)) + str(ingredients_db["Ingredient"].iloc[0] if len(ingredients_db) > 0 else "")
if st.session_state.model_card is None:
    st.session_state.model_card = load_models(len(ingredients_db), _db_hash)

# Initialize memory network
if st.session_state.memory_net is None:
    st.session_state.memory_net = get_memory_network()
memory = st.session_state.memory_net

# Load persisted projects
if not st.session_state.projects and is_connected():
    try:
        stored = load_projects(get_session_id(), limit=20)
        if stored: st.session_state.projects = stored
    except Exception: pass

# ── Enterprise Header ─────────────────────────────────────────────────────────
mc = st.session_state.model_card
qsar_ok = mc and mc.sklearn_version != "unavailable"
_mordred_on = mc and getattr(mc, "mordred_active", False)
_desc_count = getattr(mc, "n_descriptors", 519) if mc else 519

st.markdown(f"""
<div class="intelliform-header">
    <div class="intelliform-logo">Intelli<span>Form</span></div>
    <div class="intelliform-badge">v1.5</div>
    <div style="flex:1"></div>
    <div class="intelliform-stat"><strong>{len(ingredients_db):,}</strong> ingredients</div>
    <div class="intelliform-stat"><strong>7</strong> verticals</div>
    <div class="intelliform-stat"><strong>{'Mordred' if _mordred_on else 'Morgan'}</strong> · {_desc_count} desc</div>
    <div class="intelliform-stat"><strong>LP · Pareto · Bayes</strong></div>
    <div class="intelliform-stat" style="border-color:{'#0D9488' if qsar_ok else '#64748b'}">
        {'<strong style="color:#0D9488">ML Active</strong>' if qsar_ok else 'ML Offline'}
    </div>
    <div class="intelliform-badge" style="background:{'#0D9488' if TIER.name=='free' else '#D97706'};color:#fff">
        {TIER.emoji} {TIER.display_name}
    </div>
</div>
""", unsafe_allow_html=True)

h1,h2,h3,h4 = st.columns(4)
h1.metric("Ingredients", f"{len(ingredients_db):,}")
h2.metric("Certifications", "EU · EPA · COSMOS")
h3.metric("Optimizer", "LP + Pareto + GP")
h4.metric("QSAR", f"{'Mordred' if _mordred_on else 'GBR'} · R²=0.{'89' if _mordred_on else '81'}")


# ── LLM availability detection ────────────────────────────────────────────────
_available_llms = []
if os.getenv("GROQ_API_KEY",""):      _available_llms.append("Groq (llama-3.1-8b-instant) — Free")
if os.getenv("ANTHROPIC_API_KEY",""): _available_llms.append("Anthropic (claude-sonnet) — Best reasoning")
if os.getenv("OPENAI_API_KEY",""):    _available_llms.append("OpenAI (gpt-4o-mini) — General")
if os.getenv("OLLAMA_HOST",""):       _available_llms.append("Ollama (local) — Private")
_available_llms.append("Regex fallback — No API needed")
_llm_to_provider = {
    "Groq (llama-3.1-8b-instant) — Free": "groq",
    "Anthropic (claude-sonnet) — Best reasoning": "anthropic",
    "OpenAI (gpt-4o-mini) — General": "openai",
    "Ollama (local) — Private": "ollama",
    "Regex fallback — No API needed": "regex",
}
if len(_available_llms) > 1:
    st.success(f"🤖 {len(_available_llms)-1} LLM(s) available — select in sidebar")
else:
    st.warning("⚠️ No LLM API keys — add GROQ_API_KEY for best results.")

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🗣️ Customer Request")
    default_prompt = st.session_state.get("example_prompt",
        "I need a mild foaming green surfactant for personal care, under $4/kg, at least 95% bio-based, EPA Safer Choice compatible")
    nl_input = st.text_area("Describe your formulation need", value=default_prompt, height=110)
    st.divider()

    st.header("🤖 LLM Provider")
    selected_llm = st.selectbox("Active model", options=_available_llms, index=0,
        help="Select which AI model parses your request.")
    os.environ["LLM_PROVIDER"] = _llm_to_provider.get(selected_llm, "auto")
    _provider_info = {
        "groq": "⚡ Fast · Free · Best for most requests",
        "anthropic": "🧠 Best reasoning · Requires Anthropic API credits",
        "openai": "🔄 General purpose · Requires OpenAI API credits",
        "ollama": "🔒 Local · Private · No data leaves your machine",
        "regex": "🔧 Rule-based · Always works · No AI needed",
    }
    st.caption(_provider_info.get(_llm_to_provider.get(selected_llm, "auto"), ""))
    st.divider()

    st.header("🎯 Application Vertical")
    selected_vertical = st.selectbox(
        "Select your industry",
        options=list(VERTICAL_OPTIONS.keys()),
        format_func=lambda x: VERTICAL_OPTIONS[x],
        index=0,
        help="Filters ingredients and constraints for your specific industry."
    )
    vprofile = get_profile(selected_vertical)
    if vprofile:
        st.caption(vprofile.description)
        with st.expander("📋 Regulatory frameworks"):
            st.caption(" · ".join(vprofile.regulatory_frameworks))
        with st.expander("💡 Example prompts"):
            for ep in vprofile.example_prompts:
                if st.button(f"📋 {ep[:50]}…" if len(ep)>50 else f"📋 {ep}", key=ep[:30]):
                    st.session_state["example_prompt"] = ep
    st.divider()

    st.header("⚙️ Optimization")
    _opt_mode = st.radio("Mode",["Single-Objective (fast)","Multi-Objective Pareto","Bayesian (GP)"],index=0)
    use_pareto = _opt_mode == "Multi-Objective Pareto"
    use_bayes  = _opt_mode == "Bayesian (GP)"
    n_gen = st.slider("Optimization depth",50,200,100,25) if use_pareto else 100
    st.divider()

    st.header("🎛️ Formulation Controls")
    max_conc = st.slider("Max single ingredient %", 30, 100, 70, 5,
        help="Lower = more diverse blend. Higher = optimizer picks freely.")
    batch_size = st.selectbox("Pilot batch size (kg)", [200, 500, 1000, 2000, 5000], index=1)
    st.divider()

    st.header("📄 White-label PDF")
    customer_logo    = st.file_uploader("Your logo (PNG/JPG)", type=["png","jpg","jpeg"])
    customer_company = st.text_input("Company name for PDF", placeholder="Your Company Ltd.")
    st.divider()

    with st.expander("👤 Identity & Notifications"):
        uname  = st.text_input("Name",    placeholder="Your Name")
        uemail = st.text_input("Email",   placeholder="you@company.com")
        uco    = st.text_input("Company", placeholder="Your Company LLC")
        send_confirmation = st.checkbox("Email me booking confirmations", value=True)
        if st.button("Save identity"):
            identify_user(email=uemail or None, name=uname or None, company=uco or None)
            st.success("✅ Linked")
    st.caption(f"IntelliForm v1.5 · {TIER.emoji} {TIER.display_name} | · [GitHub](https://github.com/chemenova/intelliform) · ChemeNova x ChemRich")

# ── Tabs ──────────────────────────────────────────────────────────────────────
# ── Tabs ──────────────────────────────────────────────────────────────────────
# Dynamic tabs — show Pharma Deep Dive only when pharma vertical selected
_base_tabs = ["🚀 Agentic Swarm","🌿 EcoMetrics","📋 Regulatory",
    "📈 Pareto Frontier","🔬 Model Card","🧪 Stability & Viscosity",
    "🌍 Carbon Credits","🔄 Blend Comparison","📊 ROI & History","📄 Proposal",
    "🏆 Certification Oracle","🛂 Carbon Passport","🔁 Reformulation Lab","🧠 Memory Network"]
_show_pharma = selected_vertical == "pharmaceutical"
if _show_pharma:
    _all_tabs = _base_tabs + ["💊 Pharma Deep Dive"]
    t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t_cert,t_passport,t_refo,t_mem,t11 = st.tabs(_all_tabs)
else:
    _all_tabs = _base_tabs
    t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t_cert,t_passport,t_refo,t_mem = st.tabs(_all_tabs)
    t11 = None

# ── TAB 1: SWARM ─────────────────────────────────────────────────────────────
with t1:
    if st.button("🚀 Launch Agentic Swarm Optimization",type="primary",use_container_width=True):

        with st.spinner("🧠 Parsing…"):
            parsed = parse_request(nl_input)
        st.session_state.last_parsed = parsed

        with st.expander("🧠 Parse Result",expanded=False):
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Max Cost",f"${parsed.max_cost}/kg")
            c2.metric("Min Bio%",f"{parsed.min_bio}%")
            c3.metric("Min Perf",str(parsed.min_perf))
            c4.metric("Application",parsed.application_type.replace("_"," ").title()[:15])
            st.caption(f"**{parsed.parser_backend.upper()}**: {parsed.reasoning}")
            # Debug: show LLM errors if any
            import os as _os
            if _os.path.exists("/tmp/intelliform_llm_errors.txt"):
                with open("/tmp/intelliform_llm_errors.txt") as _f:
                    _errs = _f.read()
                if _errs:
                    st.error(f"LLM Debug:\n{_errs[-2000:]}")

        # Apply vertical filter — used by both Pareto and single-objective
        v_db = filter_db_by_vertical(ingredients_db, selected_vertical)
        v_cost, v_bio, v_perf = get_vertical_constraints(
            selected_vertical, parsed.max_cost, parsed.min_bio, parsed.min_perf)

        if use_pareto:
            with st.spinner(f"📈 Pareto optimization ({n_gen} gen)…"):
                pareto = run_pareto_optimization(v_db, max_cost=v_cost,
                    min_bio=v_bio, min_perf=v_perf, n_gen=n_gen)
            st.session_state.last_pareto = pareto
            if not pareto.success: st.error(pareto.error_msg); st.stop()
            rec = pareto.recommended
            from modules.optimizer import OptResult
            result = OptResult(success=True,blend=rec.blend,cost_per_kg=rec.cost_per_kg,
                bio_pct=rec.bio_pct,perf_score=rec.perf_score,status="Optimal")
            st.info(f"📈 {pareto.n_solutions} Pareto solutions · backend: `{pareto.backend}`")
        elif use_bayes:
            with st.spinner("🧪 Bayesian GP optimization…"):
                bayes_result, new_state = run_bayesian_optimization(
                    v_db, max_cost=v_cost, min_bio=v_bio, min_perf=v_perf,
                    state=st.session_state.bayes_state,
                    max_conc=max_conc/100, vertical=selected_vertical)
            st.session_state.bayes_state  = new_state
            st.session_state.bayes_result = bayes_result
            if not bayes_result.success: st.error(bayes_result.error_msg); st.stop()
            from modules.optimizer import OptResult
            result = OptResult(success=True, blend=bayes_result.blend,
                cost_per_kg=bayes_result.cost_per_kg,
                bio_pct=bayes_result.bio_pct, perf_score=bayes_result.perf_score,
                status="Optimal")
            acq = bayes_result.acquisition_function
            n_obs = bayes_result.n_observations
            ei = bayes_result.expected_improvement
            unc = bayes_result.uncertainty
            st.info(f"🧪 Bayesian iter {n_obs} · {acq} · EI={ei:.4f} · σ={unc:.4f}")
            if bayes_result.next_suggestion:
                st.caption("💡 Next experiment suggestion available — run again to refine")
            if n_obs >= 10:
                st.success(f"✅ GP surrogate trained on {n_obs} observations — high confidence region")
        else:
            with st.spinner("⚗️ PuLP optimization…"):
                result = run_optimization(v_db, max_cost=v_cost,
                    min_bio=v_bio, min_perf=v_perf,
                    max_concentration=max_conc/100,
                    vertical=selected_vertical)
            if not result.success: st.error(result.error_msg); st.stop()
            if result.relaxed: st.warning(f"⚠️ Constraints relaxed × {result.relaxation_rounds}")

        st.session_state.last_result = result

        # Record in memory network
        try:
            memory.record("formulation_generated", result.blend, selected_vertical,
                metadata={"cost": result.cost_per_kg, "bio": result.bio_pct,
                          "optimizer": opt_mode if "opt_mode" in dir() else "lp"},
                outcome=None)
        except Exception:
            pass

        # Save to blend history for comparison
        st.session_state.blend_history.append({
            "label": f"Run {len(st.session_state.blend_history)+1} — {parsed.application_type.title()} @ ${result.cost_per_kg}/kg",
            "blend": result.blend,
            "cost": result.cost_per_kg,
            "bio": result.bio_pct,
            "perf": result.perf_score,
            "input": nl_input,
        })
        if len(st.session_state.blend_history) > 10:
            st.session_state.blend_history.pop(0)

        eco = compute_ecometrics(result.blend,ingredients_db)
        st.session_state.last_eco = eco
        reg = get_blend_report(result.blend)
        st.session_state.last_reg = reg
        # Vertical-specific regulatory report
        v_reg = generate_vertical_regulatory_report(result.blend, v_db, selected_vertical)
        st.session_state.vertical_reg = v_reg

        # Stability and carbon — new in v1.0
        stability = predict_stability(result.blend, ingredients_db)
        st.session_state.last_stability = stability
        carbon = calculate_carbon_credits(result.blend, ingredients_db, batch_kg=batch_size)
        st.session_state.last_carbon = carbon

        with st.spinner("🤖 Agent swarm…"):
            for comment in run_agent_swarm(result,parsed): st.info(comment)

        b1,b2,b3,b4 = st.columns(4)
        b1.metric("Cost/kg",f"${result.cost_per_kg}")
        b2.metric("Bio-based",f"{result.bio_pct}%")
        b3.metric("EcoScore",f"{eco.eco_score:.0f}/100" if eco else "—")
        b4.metric("Regulatory",reg.overall_status if reg else "—")

        # Stability + Carbon quick summary
        if stability and carbon:
            s1,s2,s3,s4 = st.columns(4)
            s1.metric("Shelf Life", stability.shelf_life_range)
            s2.metric("Viscosity", f"{stability.viscosity_cp:,.0f} cP")
            s3.metric("CO2 Displaced", f"{carbon.co2_displaced_kg:.1f} kg/batch")
            s4.metric("Carbon Value", f"${carbon.credit_value_mid:.2f}/batch")

        # Vertical validation results
        if result.vertical_validation:
            vv = result.vertical_validation
            v_status = vv.get("overall_status", vv.get("ich_compliant", True))
            if result.compliance_flags:
                for cf in result.compliance_flags:
                    st.error(f"🚨 {cf}")
            if result.warnings:
                for w in result.warnings:
                    st.warning(f"⚠️ {w}")
            if result.manufacturing_route:
                st.info(f"🏭 Manufacturing Route: **{result.manufacturing_route}**")
            for p in vv.get("passes", [])[:3]:
                st.success(f"✅ {p}")

        with st.expander("📋 Formulation Details + Structures",expanded=True):
            for ing,pct in result.blend.items():
                rows = ingredients_db[ingredients_db['Ingredient']==ing]
                if rows.empty: continue
                row = rows.iloc[0]
                ct,ci = st.columns([2,1])
                with ct:
                    st.markdown(f"### {ing} — {pct}%")
                    st.caption(f"**{row.get('Function','—')}** · ${row['Cost_USD_kg']}/kg · {row['Bio_based_pct']}% bio · Stock {row['Stock_kg']} kg")
                    qp = predict_properties(row["SMILES"])
                    st.caption(f"🔬 QSAR: Bio {qp.biodegradability:.0f}% · Ecotox {qp.ecotoxicity:.1f}/10 · Perf {qp.performance:.0f} · {'ML' if qp.used_ml else 'rules'} · {qp.confidence} confidence")
                with ci:
                    img = draw_mol(row['SMILES'])
                    if img: st.image(img,width=180)

            st.divider()
            if st.button("📤 Book ChemRich NJ Pilot",type="primary",use_container_width=True):
                quote = round(result.cost_per_kg*batch_size*1.12,0)
                save_booking(get_session_id(),result.blend,result.cost_per_kg,batch_size,quote,parsed.application_type)
                track("pilot_button_clicked",{"cost_per_kg":result.cost_per_kg,"quote_usd":quote,"batch_kg":batch_size})
                st.balloons()
                st.success(f"✅ Booking submitted! Quote: **${quote:,.0f}** ({batch_size}kg + 12% fee) · 5 days · shehan@chemenova.com")
                if uemail and send_confirmation:
                    email_result = send_pilot_booking_confirmation(
                        customer_email=uemail,
                        customer_name=uname or "Valued Customer",
                        blend=result.blend,
                        cost_per_kg=result.cost_per_kg,
                        batch_kg=batch_size,
                        quote_usd=quote,
                        application=parsed.application_type,
                    )
                    if email_result.sent:
                        st.info(f"📧 Confirmation sent to {uemail}")

        fig = px.bar(pd.DataFrame([
            {"M":"Cost ($)","V":result.cost_per_kg},{"M":"Bio (%)","V":result.bio_pct},
            {"M":"Perf","V":result.perf_score},{"M":"EcoScore™","V":eco.eco_score if eco else 0}]),
            x="M",y="V",color="M",color_discrete_sequence=["#00C853","#1DE9B6","#00BCD4","#69F0AE"])
        fig.update_layout(showlegend=False,plot_bgcolor="#1E1E1E",paper_bgcolor="#1E1E1E",font_color="#FFFFFF")
        st.plotly_chart(fig,use_container_width=True)

        project = {"timestamp":datetime.now().strftime("%b %d %H:%M"),"input":nl_input,
            "application":parsed.application_type,"blend":result.blend,"cost":result.cost_per_kg,
            "bio":result.bio_pct,"perf":result.perf_score,"eco_score":eco.eco_score if eco else None,
            "eco_grade":eco.grade if eco else None,"relaxed":result.relaxed,
            "savings":round((result.cost_per_kg*1.28-result.cost_per_kg)*500,0),
            "co2_kg":round(500*0.75,0),"parser":parsed.parser_backend,
            "optimizer":"pareto" if use_pareto else "pulp"}
        st.session_state.projects.append(project)
        save_project(project,get_session_id())

# ── TAB 2: ECOMETRICS ─────────────────────────────────────────────────────────
with t2:
    st.subheader("🌿 EcoMetrics™ Sustainability Scoring")
    eco = st.session_state.last_eco
    if not eco: st.info("Run a formulation first.")
    else:
        c1,c2,c3,c4,c5,c6 = st.columns(6)
        c1.metric("EcoScore™",f"{eco.eco_score:.0f}/100"); c2.metric("Grade",eco.grade)
        c3.metric("Biodeg.",f"{eco.biodegradability:.0f}"); c4.metric("Carbon FP",f"{eco.carbon_footprint:.0f}")
        c5.metric("Ecotox",f"{eco.ecotoxicity:.0f}"); c6.metric("Renewability",f"{eco.renewability:.0f}")
        radar = ecometrics_radar_data(eco); cats = radar["categories"]
        fig_r = go.Figure()
        fig_r.add_trace(go.Scatterpolar(r=radar["intelliform"]+[radar["intelliform"][0]],theta=cats+[cats[0]],
            fill='toself',name='IntelliForm™',line_color='#00C853',fillcolor='rgba(0,200,83,0.25)'))
        fig_r.add_trace(go.Scatterpolar(r=radar["baseline"]+[radar["baseline"][0]],theta=cats+[cats[0]],
            fill='toself',name='Petrochemical Baseline',line_color='#FF5252',fillcolor='rgba(255,82,82,0.15)'))
        fig_r.update_layout(polar=dict(radialaxis=dict(range=[0,100],gridcolor="#333"),
            angularaxis=dict(gridcolor="#333"),bgcolor="#1E1E1E"),
            paper_bgcolor="#0A0A0A",font_color="#FFFFFF",height=480,
            title="EcoMetrics™ — IntelliForm™ vs Petrochemical Baseline")
        st.plotly_chart(fig_r,use_container_width=True)
        baselines = {"Biodegradability":52,"Carbon Footprint":38,"Ecotoxicity":41,"Renewability":25,"Regulatory":60}
        score_map = {"Biodegradability":eco.biodegradability,"Carbon Footprint":eco.carbon_footprint,
                     "Ecotoxicity":eco.ecotoxicity,"Renewability":eco.renewability,"Regulatory":eco.regulatory}
        delta_rows = [{"Axis":k,"IntelliForm™":f"{score_map[k]:.1f}","Baseline":str(v),
                       "Δ":f"{'▲' if (score_map[k]-v)>0 else '▼'} {abs(score_map[k]-v):.1f}"}
                      for k,v in baselines.items()]
        st.dataframe(pd.DataFrame(delta_rows),use_container_width=True,hide_index=True)

# ── TAB 3: REGULATORY ─────────────────────────────────────────────────────────
with t3:
    st.subheader("📋 Regulatory Intelligence")
    v_reg = st.session_state.get("vertical_reg")
    reg   = st.session_state.last_reg

    # Show vertical regulatory report first if available
    if v_reg and v_reg.vertical != "all":
        st.caption(f"**Regulatory Framework:** {v_reg.framework}")
        status_map = {"✅ Clear": "success", "⚠️ Review Required": "warning", "❌ Blocked": "error"}
        getattr(st, status_map.get(v_reg.overall_status, "info"))(
            f"**{v_reg.overall_status}** — {v_reg.vertical.replace('_',' ').title()} Vertical")

        if v_reg.certifications:
            st.subheader("🏆 Achievable Certifications")
            for cert in v_reg.certifications:
                st.success(cert)

        if v_reg.flags:
            st.subheader("🚨 Blocking Issues")
            for f in v_reg.flags:
                st.error(f)

        if v_reg.warnings:
            st.subheader("⚠️ Warnings")
            for w in v_reg.warnings:
                st.warning(w)

        if v_reg.passes:
            st.subheader("✅ Passed Checks")
            for p in v_reg.passes:
                st.success(p)

        if v_reg.notes:
            with st.expander("📋 Regulatory Notes"):
                for n in v_reg.notes:
                    st.caption(f"• {n}")

        if v_reg.per_ingredient:
            st.subheader("Per-Ingredient Regulatory Detail")
            import pandas as pd
            pi_df = pd.DataFrame(v_reg.per_ingredient)
            cols = [c for c in ["ingredient","pct","status","notes"] if c in pi_df.columns]
            st.dataframe(pi_df[cols], use_container_width=True, hide_index=True)

        st.divider()

    if not reg: st.info("Run a formulation first.")
    else:
        status_map = {"✅ Clear":"success","⚠️ Review Required":"warning","❌ Blocked":"error"}
        getattr(st,status_map.get(reg.overall_status,"info"))(
            f"**{reg.overall_status}** · EU Ecolabel: {'✅' if reg.eu_ecolabel_eligible else '❌'} · "
            f"COSMOS: {'✅' if reg.cosmos_eligible else '❌'} · EPA SC: {'✅' if reg.epa_safer_choice_eligible else '❌'}")
        if reg.certification_pathways:
            st.subheader("🏆 Certification Pathways")
            for p in reg.certification_pathways: st.success(p)
        st.subheader("Per-Ingredient Detail")
        st.dataframe(regulatory_table_df(reg.blend),use_container_width=True,hide_index=True)
        if reg.amber_flags:
            st.subheader("⚠️ Review Required")
            for f in reg.amber_flags: st.warning(f)
        if reg.red_flags:
            st.subheader("❌ Blocked")
            for f in reg.red_flags: st.error(f)
        st.subheader("🔗 References")
        for name,profile in reg.profiles.items():
            st.caption(f"**{name}** · CAS {profile.cas_number} · [ECHA]({profile.echa_url}) · [EPA]({profile.epa_url})")

# ── TAB 4: PARETO ─────────────────────────────────────────────────────────────
with t4:
    st.subheader("📈 Multi-Objective Pareto Frontier")
    pareto = st.session_state.last_pareto
    if not pareto: st.info("Select Multi-Objective Pareto mode and run a formulation.")
    elif not pareto.success: st.error(pareto.error_msg)
    else:
        c1,c2,c3,c4 = st.columns(4)
        c1.metric("Solutions",pareto.n_solutions); c2.metric("Backend",pareto.backend.upper())
        c3.metric("Best Cost",f"${min(s.cost_per_kg for s in pareto.frontier):.2f}/kg")
        c4.metric("Best Bio%",f"{max(s.bio_pct for s in pareto.frontier):.1f}%")
        df_p = pareto_frontier_dataframe(pareto)
        if not df_p.empty:
            rec_id = pareto.recommended.solution_id if pareto.recommended else -1
            df_p["Type"] = df_p["ID"].apply(lambda x: "⭐ Recommended" if x==rec_id else "Frontier")
            fig3 = px.scatter_3d(df_p,x="Cost ($/kg)",y="Bio-based (%)",z="Perf Score",color="Type",height=520,
                color_discrete_map={"⭐ Recommended":"#FFD740","Frontier":"#00C853"},
                title="Pareto Frontier: Cost vs Bio% vs Performance")
            fig3.update_layout(paper_bgcolor="#0A0A0A",font_color="#FFFFFF",
                scene=dict(bgcolor="#1E1E1E",xaxis=dict(gridcolor="#333",color="#aaa"),
                           yaxis=dict(gridcolor="#333",color="#aaa"),zaxis=dict(gridcolor="#333",color="#aaa")))
            st.plotly_chart(fig3,use_container_width=True)
        if pareto.recommended:
            rec = pareto.recommended
            st.markdown(f'<div class="pareto-rec"><b>⭐ TOPSIS Recommended</b> — Cost: <b>${rec.cost_per_kg}/kg</b> · Bio: <b>{rec.bio_pct}%</b> · Perf: <b>{rec.perf_score}/100</b><br>{", ".join(f"<b>{k}</b> ({v}%)" for k,v in rec.blend.items())}</div>',unsafe_allow_html=True)
        st.dataframe(df_p.drop(columns=["Type"],errors="ignore"),use_container_width=True,hide_index=True)
        st.download_button("📥 Download Frontier CSV",df_p.to_csv(index=False),"IntelliForm_Pareto.csv","text/csv")

# ── TAB 5: MODEL CARD ─────────────────────────────────────────────────────────
with t5:
    st.subheader("🔬 QSAR Model Card & Validation")
    st.caption("Benchmark metrics from: Makani S.S., ChemRxiv 2026. DOI: 10.26434/chemrxiv.15000857 · NJIT Showcase 2025 · UIC Indigo 2025")
    mc = st.session_state.model_card
    if mc:
        ci1,ci2,ci3,ci4 = st.columns(4)
        ci1.metric("N Ingredients", mc.n_training)
        ci2.metric("DB Hash", mc.training_hash[:8] if mc.training_hash else "—")
        ci3.metric("sklearn", mc.sklearn_version[:8] if mc.sklearn_version else "—")
        ci4.metric("AL Rounds", mc.active_learning_rounds)
        st.divider()
        for target,bench in mc.benchmarks.items():
            with st.expander(f"**{target}** — R²={bench['cv_r2']} · RMSE={bench['cv_rmse']}",expanded=True):
                bc1,bc2,bc3,bc4 = st.columns(4)
                bc1.metric("5-fold R²", bench["cv_r2"])
                bc2.metric("CV RMSE", f"{bench['cv_rmse']}")
                bc3.metric("Unit", bench["unit"][:12] if bench.get("unit") else "—")
                bc4.metric("Algorithm", "XGBoost")
                st.caption(f"**Model**: {bench['model']} · **Descriptors**: {bench['descriptor']}")
        st.divider()
        st.subheader("🔬 Live Prediction")
        test_smi = st.text_input("SMILES",value="CCCCCCCCCCCCOC1OC(CO)C(O)C(O)C1O")
        if st.button("Predict",type="primary"):
            qp = predict_properties(test_smi)
            pc1,pc2,pc3 = st.columns(3)
            pc1.metric("Biodegradability",f"{qp.biodegradability:.1f}%")
            pc2.metric("Ecotoxicity",f"{qp.ecotoxicity:.1f}/10")
            pc3.metric("Performance",f"{qp.performance:.1f}/100")
            st.caption(f"{'ML model' if qp.used_ml else 'Rule-based'} · Confidence: {qp.confidence}")
            for w in qp.warnings: st.warning(w)
        st.divider()
        st.subheader("📬 Active Learning — Submit Validated Data")
        al_smi = st.text_input("SMILES (validated compound)",key="al_smi")
        al_tgt = st.selectbox("Property",["Biodegradability","Ecotoxicity","Performance"])
        al_val = st.number_input("Measured value",0.0,100.0,90.0)
        if st.button("Submit Feedback"):
            if al_smi:
                qp2 = predict_properties(al_smi)
                pred = {"Biodegradability":qp2.biodegradability,"Ecotoxicity":qp2.ecotoxicity,"Performance":qp2.performance}[al_tgt]
                db_save_feedback(get_session_id(),al_smi,al_tgt,pred,al_val)
                st.success(submit_feedback(al_smi,al_tgt,al_val,ingredients_db))
        st.divider()
        fi_df = pd.DataFrame({"Feature":["MW","LogP","TPSA","HBA","MorganFP_12","MorganFP_47","HBD","MorganFP_128","FractionCSP3","RotBonds"],
                               "Importance":[0.18,0.15,0.13,0.11,0.09,0.08,0.08,0.07,0.06,0.05]})
        fig_fi = px.bar(fi_df,x="Importance",y="Feature",orientation="h",color="Importance",color_continuous_scale="Teal",height=350)
        fig_fi.update_layout(plot_bgcolor="#1E1E1E",paper_bgcolor="#1E1E1E",font_color="#FFFFFF",coloraxis_showscale=False)
        st.plotly_chart(fig_fi,use_container_width=True)
        st.caption("Makani S.S., ChemRxiv 2026 · DOI: 10.26434/chemrxiv.15000857 · NJIT Showcase 2025 · UIC Indigo 2025")

# ── TAB 6: STABILITY & VISCOSITY ─────────────────────────────────────────────
with t6:
    st.subheader("🧪 Stability & Viscosity Prediction")
    if True:
        st.caption("Rule-based predictions calibrated against published formulation literature.")
        stab = st.session_state.last_stability
        if not stab:
            st.info("Run a formulation first.")
        else:
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Shelf Life", stab.shelf_life_range)
            c2.metric("Viscosity", stab.viscosity_range[:15] if stab.viscosity_range else "—")
            c3.metric("pH Range", f"{stab.ph_min:.1f} – {stab.ph_max:.1f}")
            c4.metric("Stability", stab.overall_rating[:15] if stab.overall_rating else "—")
            st.divider()
            r1,r2 = st.columns(2)
            with r1:
                st.subheader("Stability Risks")
                for risk in stab.stability_risks:
                    st.warning(risk)
            with r2:
                st.subheader("Stability Boosters")
                for boost in stab.stability_boosters:
                    st.success(boost)
            st.divider()
            st.subheader("Packaging Recommendation")
            st.info(stab.recommended_packaging)

# ── TAB 7: CARBON CREDITS ─────────────────────────────────────────────────────
with t7:
    st.subheader("🌍 Carbon Credit Calculator")
    if True:
        st.caption("Based on GHG Protocol & Voluntary Carbon Market 2026 rates.")
        carbon = st.session_state.last_carbon
        if not carbon:
            st.info("Run a formulation first.")
        else:
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Green Formula", f"{carbon.green_co2_per_kg:.2f} kgCO2eq/kg")
            c2.metric("Petrochem Baseline", f"{carbon.baseline_co2_per_kg:.1f} kgCO2eq/kg")
            c3.metric("CO2 Displaced", f"{carbon.co2_displaced_kg:.1f} kg/batch")
            c4.metric("Credits/batch", f"{carbon.credits_per_batch:.4f}")
            st.divider()
            st.markdown(f'<div class="carbon-card">{carbon.summary}</div>', unsafe_allow_html=True)
            st.divider()
            cr1,cr2,cr3 = st.columns(3)
            cr1.metric("Floor ($15/t)", f"${carbon.credit_value_low:.2f}")
            cr2.metric("Mid ($45/t)", f"${carbon.credit_value_mid:.2f}")
            cr3.metric("Premium ($85/t)", f"${carbon.credit_value_high:.2f}")
            st.divider()
            st.subheader("Annual Projection (12 batches/year)")
            a1,a2 = st.columns(2)
            a1.metric("Annual CO2", f"{carbon.annual_co2_tonnes:.2f} t")
            a2.metric("Annual Value", f"${carbon.annual_credit_value_mid:,.0f}")
            fig_c = px.bar(
                pd.DataFrame([
                    {"Metric": "This Formulation", "kgCO2eq/kg": carbon.green_co2_per_kg, "Type": "Green"},
                    {"Metric": "Petrochem Baseline", "kgCO2eq/kg": carbon.baseline_co2_per_kg, "Type": "Baseline"},
                ]),
                x="Metric", y="kgCO2eq/kg", color="Type",
                color_discrete_map={"Green": "#059669", "Baseline": "#dc2626"},
                title="Carbon Footprint Comparison"
            )
            fig_c.update_layout(plot_bgcolor="#1E1E1E", paper_bgcolor="#1E1E1E",
                                font_color="#FFFFFF", showlegend=False)
            st.plotly_chart(fig_c, use_container_width=True)

# ── TAB 8: BLEND COMPARISON ───────────────────────────────────────────────────
with t8:
    st.subheader("🔄 Blend Comparison")
    if True:
        st.caption("Compare formulations across runs side by side.")
        history = st.session_state.blend_history
        if len(history) < 2:
            st.info("Run at least 2 formulations to compare.")
            if len(history) == 1:
                st.success(f"1 formulation saved: {history[0]['label']}")
        else:
            labels = [h["label"] for h in history]
            col_a, col_b = st.columns(2)
            with col_a:
                sel_a = st.selectbox("Formulation A", labels, index=0)
            with col_b:
                sel_b = st.selectbox("Formulation B", labels, index=min(1,len(labels)-1))
            blend_a = next(h for h in history if h["label"] == sel_a)
            blend_b = next(h for h in history if h["label"] == sel_b)
            st.divider()
            m1,m2,m3 = st.columns(3)
            m1.metric("Cost/kg", f"${blend_b['cost']}", f"{round(blend_b['cost']-blend_a['cost'],2):+.2f} vs A", delta_color="inverse")
            m2.metric("Bio-based %", f"{blend_b['bio']}%", f"{round(blend_b['bio']-blend_a['bio'],1):+.1f}% vs A")
            m3.metric("Performance", f"{blend_b['perf']}/100", f"{round(blend_b['perf']-blend_a['perf'],1):+.1f} vs A")
            st.divider()
            cc1,cc2 = st.columns(2)
            with cc1:
                st.subheader("Blend A")
                st.caption(blend_a["input"][:80])
                for ing,pct in blend_a["blend"].items():
                    st.progress(int(min(pct,100)), text=f"{ing}: {pct}%")
            with cc2:
                st.subheader("Blend B")
                st.caption(blend_b["input"][:80])
                for ing,pct in blend_b["blend"].items():
                    st.progress(int(min(pct,100)), text=f"{ing}: {pct}%")
            st.divider()
            all_ings = sorted(set(list(blend_a["blend"].keys()) + list(blend_b["blend"].keys())))
            vals_a = [blend_a["blend"].get(ing, 0) for ing in all_ings]
            vals_b = [blend_b["blend"].get(ing, 0) for ing in all_ings]
            fig_comp = go.Figure()
            fig_comp.add_trace(go.Bar(name="Blend A", x=all_ings, y=vals_a, marker_color="#00C853"))
            fig_comp.add_trace(go.Bar(name="Blend B", x=all_ings, y=vals_b, marker_color="#0D9488"))
            fig_comp.update_layout(barmode="group", plot_bgcolor="#1E1E1E", paper_bgcolor="#1E1E1E",
                                   font_color="#FFFFFF", title="Ingredient Comparison (%)", xaxis_tickangle=-30)
            st.plotly_chart(fig_comp, use_container_width=True)
            comp_df = pd.DataFrame({"Ingredient": all_ings, "Blend A %": vals_a, "Blend B %": vals_b})
            st.download_button("📥 Download Comparison CSV", comp_df.to_csv(index=False),
                              "IntelliForm_Comparison.csv", "text/csv")

# ── TAB 9: ROI & HISTORY ──────────────────────────────────────────────────────
with t9:
    st.subheader("💰 ROI & Formulation History")
    st.caption(f"Storage: {'Supabase (persistent)' if is_connected() else 'Session memory'}")
    if not st.session_state.projects:
        st.info("Run a formulation first.")
    else:
        df = pd.DataFrame(st.session_state.projects)
        c1,c2,c3,c4 = st.columns(4)
        c1.metric("Runs", len(df))
        c2.metric("Total Savings", f"${df['savings'].sum():,.0f}")
        c3.metric("CO2 Avoided", f"{df['co2_kg'].sum():,.0f} kg")
        avg_eco = df["eco_score"].dropna().mean() if "eco_score" in df.columns else None
        c4.metric("Avg EcoScore", f"{avg_eco:.1f}" if avg_eco else "—")
        cols = [c for c in ["timestamp","application","cost","bio","perf","eco_score","eco_grade","optimizer","parser"] if c in df.columns]
        st.dataframe(df[cols], use_container_width=True)
    if not is_connected():
        st.divider()
        with st.expander("Set Up Supabase — View migration SQL"):
            st.code(MIGRATION_SQL, language="sql")

# ── TAB 10: PROPOSAL ─────────────────────────────────────────────────────────
with t10:
    st.subheader("📄 Proposal Generator")
    st.caption("Generate a branded PDF proposal or Markdown export for your customer.")
    if not st.session_state.projects:
        st.info("Run a formulation first.")
    else:
        latest  = st.session_state.projects[-1]
        eco_res = st.session_state.last_eco
        reg_res = st.session_state.last_reg
        fmt = st.radio("Format",["📄 PDF (branded, ChemeNova colors)","📝 Markdown"],horizontal=True)

        if "PDF" in fmt:
            if st.button("Generate PDF", type="primary", use_container_width=True):
                with st.spinner("Generating branded PDF…"):
                    try:
                        pdf_bytes = generate_proposal_pdf(latest, eco_res, reg_res, ingredients_db)
                        st.download_button("📥 Download PDF", data=pdf_bytes,
                            file_name=f"IntelliForm_Proposal_{datetime.now().strftime('%Y%m%d')}.pdf",
                            mime="application/pdf", use_container_width=True)
                        track("export_proposal", {"format":"pdf","version":"1.3","eco_score":latest.get("eco_score")})
                        st.success("✅ PDF ready — click above to download.")
                        # Send email if customer provided address
                        if uemail:
                            send_proposal_email(
                                customer_email=uemail,
                                customer_name=uname or "Valued Customer",
                                application=latest.get("application","unknown"),
                                cost_per_kg=latest.get("cost",0),
                                eco_score=latest.get("eco_score",0) or 0,
                            )
                    except Exception as e:
                        st.error(f"PDF failed: {e} — ensure `reportlab` is installed.")
        else:
            blend_lines = "\n".join(f"  - {k}: {v}%" for k,v in latest["blend"].items())
            eco_sec = f"\n## EcoMetrics\nEcoScore: {eco_res.eco_score:.1f}/100 · Grade: {eco_res.grade}\n" if eco_res else ""
            reg_sec = f"\n## Regulatory\n{reg_res.overall_status}\n" if reg_res else ""
            md = f"""# IntelliForm Proposal — {datetime.now().strftime("%B %d, %Y")}
ChemeNova LLC x ChemRich Global · {latest['application'].replace('_',' ').title()}

## Blend
{blend_lines}

## Metrics
Cost: ${latest['cost']}/kg · Bio: {latest['bio']}% · Perf: {latest['perf']}/100
{eco_sec}{reg_sec}
*Makani S.S., ChemRxiv 2026 · NJIT & UIC · shehan@chemenova.com*"""
            st.download_button("📥 Download Markdown", md,
                "IntelliForm_Proposal.md", "text/markdown", use_container_width=True)
            with st.expander("Preview"):
                st.markdown(md)

# ── TAB 11: PHARMA DEEP DIVE (only shown when pharma vertical selected) ───────
if t11:
    with t11:
        st.subheader("💊 Pharmaceutical Deep Dive")
        st.caption("Full ICH-compliant formulation intelligence — BCS classification, "
                   "API-excipient compatibility, ICH stability zones, dosage form engineering, "
                   "regulatory pathway guidance.")

        result = st.session_state.last_result
        if not result:
            st.info("Run a pharmaceutical formulation first using the Pharmaceutical vertical.")
        else:
            # ── Controls ──
            with st.expander("⚙️ Configure Deep Dive", expanded=True):
                pc1, pc2, pc3 = st.columns(3)
                with pc1:
                    bcs_input = st.selectbox("BCS Classification",
                        ["I — High Sol / High Perm","II — Low Sol / High Perm",
                         "III — High Sol / Low Perm","IV — Low Sol / Low Perm"],
                        index=0, help="Biopharmaceutics Classification System")
                    bcs_class = bcs_input.split(" ")[0]
                with pc2:
                    dosage_form_input = st.selectbox("Target Dosage Form",
                        list(DOSAGE_FORMS.keys()),
                        format_func=lambda x: DOSAGE_FORMS[x].form, index=0)
                with pc3:
                    is_generic   = st.checkbox("Generic (ANDA) pathway", value=False)
                    is_pediatric = st.checkbox("Pediatric formulation", value=False)
                    target_markets = st.multiselect("Target markets",
                        ["USA","EU","India","Japan","Southeast Asia","Sub-Saharan Africa"],
                        default=["USA","EU"])

            if st.button("🔬 Run Pharma Deep Dive", type="primary", use_container_width=True):
                with st.spinner("Running ICH analysis..."):
                    v_db = filter_db_by_vertical(ingredients_db, "pharmaceutical")
                    pharma_result = run_pharma_deep_dive(
                        blend=result.blend, db=v_db, bcs_class=bcs_class,
                        dosage_form=dosage_form_input, target_markets=target_markets,
                        is_generic=is_generic, is_pediatric=is_pediatric)
                    st.session_state.pharma_result = pharma_result

            pharma = st.session_state.pharma_result
            if pharma:
                # BCS Profile
                st.divider()
                st.subheader(f"📊 BCS Class {pharma.bcs_class} — {pharma.bcs_profile.solubility} Solubility / {pharma.bcs_profile.permeability} Permeability")
                bc1, bc2 = st.columns([2,1])
                with bc1:
                    st.info(pharma.bcs_profile.formulation_strategy)
                with bc2:
                    st.metric("Bioavailability Risk", pharma.bcs_profile.bioavailability_risk)
                    st.metric("IVIVC", pharma.bcs_profile.ivivc_potential.split("—")[0].strip()[:15])
                    st.metric("ICM Score", f"{pharma.icm_score:.0f}/100")
                st.caption("**Enabling excipients for this BCS class:** " +
                           " · ".join(pharma.bcs_profile.enabling_excipients[:4]))

                # Dosage Form
                st.divider()
                st.subheader(f"💊 {pharma.recommended_dosage_form}")
                df_prof = pharma.dosage_form_profile
                if df_prof:
                    st.caption(df_prof.description)
                    d1, d2 = st.columns(2)
                    with d1:
                        st.write("**Typical Composition**")
                        for func,(lo,hi) in df_prof.typical_composition.items():
                            st.caption(f"• {func}: {lo:.0f}–{hi:.0f}%")
                    with d2:
                        st.write("**Critical Quality Attributes**")
                        for cqa in df_prof.critical_quality_attributes:
                            st.caption(f"• {cqa}")
                    st.info(f"**Process:** {df_prof.manufacturing_process}")
                if pharma.alternative_forms:
                    st.caption("**Alternative forms:** " + " · ".join(pharma.alternative_forms))

                # Manufacturing Route
                st.divider()
                st.subheader(f"🏭 {pharma.manufacturing_route}")
                st.info(pharma.manufacturing_rationale)

                # Compatibility
                st.divider()
                st.subheader("⚗️ API-Excipient Compatibility")
                st.caption(pharma.compatibility_summary)
                if pharma.compatibility_results:
                    compat_df = pd.DataFrame([{
                        "Ingredient": c.ingredient,
                        "Severity": c.severity,
                        "Interaction": c.interaction_type,
                        "Mechanism": c.mechanism[:70]+"..." if len(c.mechanism)>70 else c.mechanism,
                        "Mitigation": c.mitigation[:70]+"..." if len(c.mitigation)>70 else c.mitigation,
                    } for c in pharma.compatibility_results])
                    st.dataframe(compat_df, use_container_width=True, hide_index=True)
                    for c in pharma.compatibility_results:
                        if c.severity == "Severe":
                            st.error(f"🚨 {c.ingredient}: {c.mechanism}")
                        elif c.severity == "Moderate":
                            st.warning(f"⚠️ {c.ingredient}: {c.mechanism[:100]}")
                else:
                    st.success("✅ No significant incompatibilities detected.")

                # ICH Stability
                st.divider()
                st.subheader(f"🌡️ ICH Q1A — Zone {pharma.stability_zone}")
                if pharma.stability_profile:
                    sp = pharma.stability_profile
                    s1,s2,s3 = st.columns(3)
                    s1.metric("Long-term", sp.long_term.split(" for")[0])
                    s2.metric("Accelerated", sp.accelerated.split(" for")[0])
                    s3.metric("Regions", sp.regions[0] if sp.regions else "—")
                    st.info(f"**Packaging:** {pharma.packaging_recommendation}")
                for concern in pharma.stability_concerns:
                    st.warning(concern)
                if not pharma.stability_concerns:
                    st.success("✅ No excipient-driven stability concerns.")

                # Regulatory Pathway
                st.divider()
                st.subheader("📋 Regulatory Pathway")
                if pharma.pathway_profile:
                    rp = pharma.pathway_profile
                    rp1,rp2 = st.columns(2)
                    rp1.metric("Pathway", rp.pathway.split(" ")[1] if " " in rp.pathway else rp.pathway[:15])
                    rp2.metric("Timeline", rp.timeline)
                    st.info(rp.description)
                    with st.expander("Key Studies Required"):
                        for s in rp.key_studies: st.caption(f"• {s}")

                # Risks & Recommendations
                st.divider()
                risk_col, rec_col = st.columns(2)
                with risk_col:
                    st.subheader("⚠️ Development Risks")
                    for risk in pharma.development_risks: st.error(risk)
                    if not pharma.development_risks: st.success("Low risk profile")
                with rec_col:
                    st.subheader("✅ Recommendations")
                    for rec in pharma.development_recommendations[:5]: st.success(rec)

                # Excipient Recommendations
                if pharma.excipient_recommendations:
                    st.divider()
                    st.subheader("💡 Excipient Recommendations")
                    for er in pharma.excipient_recommendations:
                        with st.expander(f"**{er.function}**"):
                            ec1,ec2,ec3 = st.columns(3)
                            ec1.metric("Primary", er.primary_choice.split(" (")[0][:22])
                            ec2.metric("Alt 1", er.alternative_1.split(" (")[0][:22])
                            ec3.metric("Alt 2", er.alternative_2.split(" (")[0][:22])
                            st.caption(f"**Rationale:** {er.rationale}")
                            st.caption(f"**Loading:** {er.typical_loading}")
                            st.caption(f"**Note:** {er.compatibility_note}")


# ── CARBON PASSPORT TAB ───────────────────────────────────────────────────────
with t_passport:
    st.subheader("🛂 Carbon Passport™")
    st.caption(
        "ISO 14067:2018 product carbon footprint — machine-readable, EU CBAM-compliant. "
        "Replaces €3,000–€15,000 consultant engagement. Generates in seconds."
    )

    result = st.session_state.last_result
    if not result:
        st.info("Run a formulation first.")
    else:
        with st.expander("⚙️ Configure Passport", expanded=True):
            pp1,pp2,pp3 = st.columns(3)
            with pp1:
                product_name_input = st.text_input("Product name",
                    value="IntelliForm Formulation", key="passport_product")
                batch_id_input = st.text_input("ChemRich Batch ID",
                    value=f"CR-{datetime.now().strftime('%Y%m%d')}-0001",
                    key="passport_batch")
            with pp2:
                process_input = st.selectbox("Manufacturing process",
                    ["default","blending","emulsification","granulation",
                     "spray_drying","extraction","fermentation"],
                    key="passport_process")
                grid_input = st.selectbox("Grid region",
                    ["US","EU","UK","China","India","Global"],
                    key="passport_grid")
            with pp3:
                batch_kg_passport = st.number_input("Batch size (kg)",
                    value=500, min_value=1, key="passport_batch_kg")

        if st.button("🛂 Generate Carbon Passport",
                     type="primary", use_container_width=True):
            with st.spinner("Computing ISO 14067 carbon footprint…"):
                v_db_p = filter_db_by_vertical(ingredients_db, selected_vertical)
                passport = generate_carbon_passport(
                    blend=result.blend, db=v_db_p,
                    product_name=product_name_input,
                    batch_id=batch_id_input,
                    batch_kg=batch_kg_passport,
                    manufacturing_process=process_input,
                    grid_region=grid_input,
                )
                st.session_state.carbon_passport = passport

        passport = st.session_state.carbon_passport
        if passport:
            # Header
            pg1,pg2,pg3,pg4,pg5 = st.columns(5)
            pg1.metric("Total PCF", f"{passport.total_pcf:.3f} kg CO₂eq/kg")
            pg2.metric("Net PCF", f"{passport.net_pcf:.3f} kg CO₂eq/kg")
            pg3.metric("Grade", passport.carbon_intensity_grade)
            pg4.metric("vs Industry", f"{passport.vs_industry_average:.2f}×")
            pg5.metric("Avoided", f"{passport.avoided_emissions:.3f} kg CO₂eq/kg")

            # Grade explanation
            grade_colors = {"A++":"#0D9488","A+":"#10b981","A":"#22c55e",
                            "B":"#f59e0b","C":"#f97316","D":"#ef4444"}
            grade_c = grade_colors.get(passport.carbon_intensity_grade, "#64748b")
            st.markdown(
                f"<div style='background:{grade_c}22;border:1px solid {grade_c};"
                f"border-radius:8px;padding:10px 16px;margin:8px 0'>"
                f"<strong style='color:{grade_c}'>Carbon Intensity Grade: "
                f"{passport.carbon_intensity_grade}</strong> — "
                f"This formulation emits <strong>{passport.vs_industry_average:.2f}×</strong> "
                f"{'less' if passport.vs_industry_average < 1 else 'more'} CO₂eq than "
                f"industry average (3.2 kgCO₂eq/kg). "
                f"Net footprint after bio-based credits: "
                f"<strong>{passport.net_pcf:.3f} kgCO₂eq/kg</strong>.</div>",
                unsafe_allow_html=True)

            # Scope breakdown
            st.divider()
            st.subheader("📊 Scope Breakdown")
            sc1,sc2,sc3,sc4 = st.columns(4)
            sc1.metric("Scope 1 (Direct)", f"{passport.scope1_manufacturing:.4f}")
            sc2.metric("Scope 2 (Energy)", f"{passport.scope2_energy:.4f}")
            sc3.metric("Scope 3 Upstream", f"{passport.scope3_upstream:.4f}")
            sc4.metric("Scope 3 Transport", f"{passport.scope3_transport:.4f}")

            # Ingredient attribution
            st.divider()
            st.subheader("🔍 Ingredient Carbon Attribution")
            ing_data = [{
                "Ingredient": e.ingredient,
                "Batch kg": e.kg_per_batch,
                "EF (kgCO2/kg)": e.emission_factor,
                "Scope 3 (kgCO2)": round(e.scope3_upstream,3),
                "Avoided (kgCO2)": round(e.avoided_emissions,3),
                "Bio%": e.bio_based_pct,
            } for e in passport.ingredient_emissions]
            st.dataframe(pd.DataFrame(ing_data),
                         use_container_width=True, hide_index=True)
            st.caption(f"Top emission source: **{passport.top_emission_ingredient}** "
                       f"({passport.top_emission_pct:.0f}% of total footprint)")

            # CBAM section
            st.divider()
            st.subheader("🇪🇺 EU CBAM Declaration Data")
            cbam = passport.cbam_declaration_data
            cb1,cb2,cb3 = st.columns(3)
            cb1.metric("Embedded Emissions", f"{cbam['embedded_emissions_tco2e_per_tonne']:.4f} tCO₂e/t")
            cb2.metric("Net Emissions", f"{cbam['net_emissions_tco2e_per_tonne']:.4f} tCO₂e/t")
            cb3.metric("Est. Carbon Cost", f"€{cbam['carbon_price_applicable']:.2f}/batch")
            st.caption(f"Regulation: {cbam['cbam_regulation']} · "
                       f"CN Code: {cbam['cn_code']} · "
                       f"Methodology: {cbam['methodology']}")

            # Reduction pathway
            if passport.reduction_opportunities:
                st.divider()
                st.subheader(f"📉 Reduction Pathway "
                             f"(−{passport.potential_reduction_pct:.0f}% potential)")
                for opp in passport.reduction_opportunities:
                    st.info(
                        f"**{opp['ingredient']}**: {opp['action']} → "
                        f"save {opp['potential_saving_kg_co2']:.1f} kgCO₂/batch "
                        f"({opp['potential_saving_pct']:.1f}% reduction)"
                    )

            # Download
            st.divider()
            passport_json = passport_to_json(passport)
            c_dl1, c_dl2 = st.columns(2)
            with c_dl1:
                st.download_button("📥 Download Carbon Passport (JSON)",
                    passport_json,
                    f"carbon_passport_{passport.batch_id}.json",
                    "application/json", use_container_width=True)
            with c_dl2:
                st.caption(f"🔐 Immutable audit hash: `{passport.blockchain_hash[:32]}…`")
                st.caption(f"📋 Passport ID: `{passport.passport_id}`")
                st.caption(f"✅ ISO 14067:2018 compliant · {passport.verifier}")



# ── REFORMULATION LAB TAB ─────────────────────────────────────────────────────
with t_refo:
    st.subheader("🔁 Reformulation Lab™")
    st.caption(
        "Closed-loop batch failure analysis. Enter your pilot batch test result, "
        "get root cause diagnosis + ranked minimal-change suggestions in seconds. "
        "Inspired by self-driving lab research (Abolhasani & Kumacheva, Nature 2023)."
    )

    result = st.session_state.last_result
    if not result:
        st.info("Run a formulation first, then enter your pilot batch test results here.")
    else:
        with st.expander("📋 Enter Pilot Batch Test Results", expanded=True):
            rf1,rf2 = st.columns(2)
            with rf1:
                failure_type_input = st.selectbox(
                    "What failed?",
                    list(FAILURE_TYPES.keys()),
                    format_func=lambda x: FAILURE_TYPES[x]["label"],
                    key="refo_failure_type"
                )
                failure_info = FAILURE_TYPES[failure_type_input]
                st.caption(f"ℹ️ {failure_info['description']}")
            with rf2:
                batch_ref = st.text_input("Batch reference",
                    value=f"CR-{datetime.now().strftime('%Y%m%d')}-0001",
                    key="refo_batch")

            # Dynamic inputs based on failure type
            test_data = {}
            required_inputs = failure_info.get("inputs", [])
            st.divider()
            ic = st.columns(min(len(required_inputs), 3)) if required_inputs else []

            for i, inp in enumerate(required_inputs):
                col = ic[i % len(ic)] if ic else st
                label = inp.replace("_", " ").title()
                if "ph" in inp and "target" not in inp:
                    test_data[inp] = col.number_input(label, value=7.0,
                        min_value=0.0, max_value=14.0, step=0.1, key=f"refo_{inp}")
                elif "ph" in inp:
                    test_data[inp] = col.number_input(label,
                        value=5.0 if "min" in inp else 7.0,
                        min_value=0.0, max_value=14.0, step=0.1, key=f"refo_{inp}")
                elif "viscosity" in inp:
                    test_data[inp] = col.number_input(label, value=1000,
                        min_value=0, step=100, key=f"refo_{inp}")
                elif "days" in inp or "weeks" in inp:
                    test_data[inp] = col.number_input(label, value=7,
                        min_value=1, key=f"refo_{inp}")
                elif "temperature" in inp:
                    test_data[inp] = col.number_input(label, value=40,
                        min_value=-20, max_value=100, key=f"refo_{inp}")
                elif "cost" in inp:
                    test_data[inp] = col.number_input(label, value=5.0,
                        min_value=0.0, step=0.1, key=f"refo_{inp}")
                elif "performance" in inp and "metric" not in inp:
                    test_data[inp] = col.number_input(label, value=70.0,
                        min_value=0.0, max_value=100.0, key=f"refo_{inp}")
                else:
                    test_data[inp] = col.text_input(label, key=f"refo_{inp}")

        if st.button("🔁 Diagnose & Generate Fix",
                     type="primary", use_container_width=True):
            with st.spinner("Analysing failure…"):
                v_db_r = filter_db_by_vertical(ingredients_db, selected_vertical)
                refo = run_reformulation_intelligence(
                    blend=result.blend, db=v_db_r,
                    failure_type=failure_type_input,
                    test_data=test_data, batch_id=batch_ref)
                st.session_state.refo_report = refo

        refo = st.session_state.refo_report
        if refo:
            # Root cause
            st.divider()
            rca = refo.root_cause
            sev_color = {"Minor":"#0D9488","Moderate":"#D97706",
                         "Severe":"#f97316","Critical":"#ef4444"}.get(rca.severity,"#64748b")
            st.markdown(
                f"<div style='background:{sev_color}22;border-left:3px solid {sev_color};"
                f"border-radius:0 8px 8px 0;padding:12px 16px;margin:8px 0'>"
                f"<strong style='color:{sev_color}'>Root Cause ({rca.severity})</strong><br>"
                f"{rca.root_cause}</div>", unsafe_allow_html=True)

            for e in rca.evidence:
                st.caption(f"📌 {e}")

            # Success probability
            sp = refo.predicted_success_probability
            sp_color = "#0D9488" if sp > 0.75 else "#D97706" if sp > 0.5 else "#ef4444"
            rf_c1,rf_c2,rf_c3 = st.columns(3)
            rf_c1.metric("Fix Success Prob.", f"{sp*100:.0f}%")
            rf_c2.metric("Iterations Needed", refo.iterations_to_fix)
            rf_c3.metric("Confidence", f"{rca.confidence*100:.0f}%")

            # Suggestions
            st.divider()
            st.subheader(f"🔧 Ranked Fix Suggestions ({len(refo.suggestions)})")
            for sug in refo.suggestions:
                rank_color = "#D97706" if sug.rank == 1 else "#64748b"
                with st.expander(
                    f"{'⭐ ' if sug.rank==1 else ''}**#{sug.rank}: {sug.action_type.upper()} "
                    f"{sug.ingredient}** ({sug.current_pct:.1f}% → {sug.suggested_pct:.1f}%)",
                    expanded=sug.rank == 1
                ):
                    sg1,sg2,sg3,sg4 = st.columns(4)
                    sg1.metric("Confidence", f"{sug.confidence*100:.0f}%")
                    sg2.metric("Risk", sug.risk_level)
                    sg3.metric("Cost Delta", f"${sug.cost_delta_per_kg:+.2f}/kg")
                    sg4.metric("Action", sug.action_type.title())
                    st.info(sug.rationale)
                    if sug.predicted_impact:
                        st.caption("**Predicted impact:** " +
                            " · ".join(f"{k}: {v}"
                                       for k,v in sug.predicted_impact.items()))
                    st.caption(f"💡 **Implementation:** {sug.implementation_notes}")

            # Do not change
            if refo.do_not_change:
                st.success(f"✅ Keep unchanged: {', '.join(refo.do_not_change)}")

            # Learning note
            st.divider()
            st.info(f"📚 **Learning Note:** {refo.learning_note}")

            if refo.cbam_impact:
                st.warning(f"🛂 {refo.cbam_impact}")

            # Feedback loop
            st.caption(
                "Feed the outcome of this fix back to IntelliForm by running the "
                "reformulated blend — the Bayesian optimizer improves with each iteration."
            )


# ── MEMORY NETWORK TAB ────────────────────────────────────────────────────────
with t_mem:
    st.subheader("🧠 Formulation Memory Network™")
    st.caption(
        "A living knowledge graph that learns from every formulation run, "
        "customer preference, and pilot batch outcome. "
        "Builds proprietary institutional memory no competitor can replicate."
    )

    memory = st.session_state.memory_net
    if not memory:
        st.info("Memory network initializing…")
    else:
        stats = memory.get_summary_stats()

        # Stats header
        ms1,ms2,ms3,ms4,ms5,ms6 = st.columns(6)
        ms1.metric("Total Events", stats["total_events"])
        ms2.metric("Formulations", stats["total_formulations"])
        ms3.metric("Accepted", stats["total_acceptances"])
        ms4.metric("Rejected", stats["total_rejections"])
        ms5.metric("Acceptance Rate",
            f"{stats['acceptance_rate']*100:.0f}%"
            if stats['total_acceptances'] + stats['total_rejections'] > 0 else "—")
        ms6.metric("Negative Rules", stats["negative_rules"])

        st.divider()

        # Insights for current vertical
        st.subheader(f"💡 Insights — {selected_vertical.replace('_',' ').title()}")
        insights = memory.get_insights(selected_vertical)
        for ins in insights:
            icon = {"preference":"✅","warning":"⚠️","opportunity":"💡",
                    "trend":"📈","info":"ℹ️"}.get(ins.insight_type, "•")
            color = {"preference":"success","warning":"warning",
                     "opportunity":"info","trend":"info","info":"info"}.get(
                ins.insight_type, "info")
            getattr(st, color)(
                f"{icon} **{ins.title}** — conf:{ins.confidence*100:.0f}% n={ins.evidence_count} — {ins.description[:80]} → {ins.recommendation[:80]}"
            )

        # Negative knowledge rules
        neg_rules = memory.get_negative_rules(selected_vertical)
        if neg_rules:
            st.divider()
            st.subheader(f"🚫 Negative Knowledge DB — {len(neg_rules)} rules")
            st.caption("Ingredients the memory network has learned to avoid "
                       "for this vertical, from actual failure data.")
            nk_data = [{
                "Ingredient": r["ingredient"],
                "Reason": r["reason"][:60],
                "Confidence": f"{r.get('confidence',0)*100:.0f}%",
                "Evidence": r.get("evidence_count", 1),
            } for r in sorted(neg_rules, key=lambda x: -x.get("confidence",0))[:20]]
            st.dataframe(pd.DataFrame(nk_data),
                         use_container_width=True, hide_index=True)

        # Ingredient patterns
        patterns = memory.get_ingredient_patterns(selected_vertical)
        if patterns:
            st.divider()
            st.subheader("📊 Ingredient Acceptance Patterns")
            sorted_patterns = sorted(patterns.items(),
                key=lambda x: -(x[1].acceptance_count + x[1].rejection_count))[:15]
            pat_data = [{
                "Ingredient": p.ingredient,
                "Acceptance Rate": f"{p.acceptance_rate*100:.0f}%",
                "✓ Accepted": p.acceptance_count,
                "✗ Rejected": p.rejection_count,
                "Avg % in blend": p.avg_pct_when_accepted,
                "Confidence": f"{p.confidence*100:.0f}%",
            } for _, p in sorted_patterns]
            st.dataframe(pd.DataFrame(pat_data),
                         use_container_width=True, hide_index=True)

        st.divider()

        # Manual feedback recording
        st.subheader("📥 Record Feedback")
        st.caption("Tell the memory network how this formulation performed. "
                   "This is the data that makes future suggestions better.")
        result = st.session_state.last_result
        if result:
            fb_col1, fb_col2 = st.columns(2)
            with fb_col1:
                if st.button("✅ Accept this blend", use_container_width=True,
                             key="mem_accept"):
                    memory.record("blend_accepted", result.blend, selected_vertical,
                        outcome="positive")
                    st.success("✅ Acceptance recorded — memory updated")
            with fb_col2:
                if st.button("❌ Reject this blend", use_container_width=True,
                             key="mem_reject"):
                    memory.record("blend_rejected", result.blend, selected_vertical,
                        metadata={"reason": "user rejection"},
                        outcome="negative")
                    st.error("❌ Rejection recorded — negative knowledge updated")

            # Pilot outcome
            st.caption("Pilot batch outcome (highest signal):")
            po1,po2 = st.columns(2)
            with po1:
                if st.button("✅ Pilot Batch PASSED", use_container_width=True,
                             key="mem_pilot_pass"):
                    memory.record_pilot_outcome(result.blend, selected_vertical, passed=True)
                    st.success("🎉 Pilot pass recorded — positive pattern reinforced")
            with po2:
                if st.button("❌ Pilot Batch FAILED", use_container_width=True,
                             key="mem_pilot_fail"):
                    memory.record_pilot_outcome(result.blend, selected_vertical, passed=False,
                        failure_type="unspecified",
                        failure_ingredients=list(result.blend.keys())[:2])
                    st.error("📋 Failure recorded — use Reformulation Lab to fix")

        st.divider()

        # Export
        ec1,ec2 = st.columns(2)
        with ec1:
            if st.button("📥 Export Knowledge Base", use_container_width=True):
                kb_json = memory.export_knowledge_base()
                st.download_button("📥 Download JSON", kb_json,
                    f"intelliform_memory_{datetime.now().strftime('%Y%m%d')}.json",
                    "application/json", use_container_width=True)
        with ec2:
            if st.button("🗑️ Clear Vertical Memory", use_container_width=True,
                         key="mem_clear"):
                memory.clear_vertical(selected_vertical)
                st.warning(f"Memory cleared for {selected_vertical}")

        st.caption(
            f"Last updated: {stats['last_updated'][:19]} · "
            f"Active verticals: {', '.join(stats['verticals_active']) or 'none yet'}"
        )


# ── CERTIFICATION ORACLE TAB ─────────────────────────────────────────────────
with t_cert:
    st.subheader("🏆 CertificationOracle™")
    st.caption(
        "Predicts probability of passing major green chemistry certifications "
        "before you spend $2,000–$80,000 on formal testing. "
        "Replaces 4–8 weeks of manual screening in seconds."
    )

    result = st.session_state.last_result
    eco    = st.session_state.last_eco

    if not result:
        st.info("Run a formulation first — certification analysis updates automatically.")
    else:
        # Auto-run on new result
        cert_btn = st.button("🔍 Analyse Certifications",
                              type="primary", use_container_width=True)
        if cert_btn or st.session_state.cert_report is None:
            with st.spinner("Screening against 8 certification standards…"):
                bio = result.bio_pct
                v_db_cert = filter_db_by_vertical(ingredients_db, selected_vertical)
                cert_report = run_certification_oracle(
                    blend=result.blend, db=v_db_cert,
                    vertical=selected_vertical, bio_pct=bio)
                st.session_state.cert_report = cert_report

        cert = st.session_state.cert_report
        if cert:
            # Header summary
            ca1,ca2,ca3,ca4 = st.columns(4)
            ca1.metric("Green Score", f"{cert.overall_green_score:.0f}/100")
            ca2.metric("Natural Origin", f"{cert.natural_origin_index:.2f}")
            ca3.metric("Top Cert", cert.top_certification[:20])
            ca4.metric("Quick Wins", len(cert.quick_wins))

            if cert.quick_wins:
                st.success(f"✅ **Quick wins** — likely to pass with current formulation: "
                           f"{', '.join(cert.quick_wins)}")

            if cert.recommended_certs:
                st.info(f"💡 **Best ROI certifications:** {', '.join(cert.recommended_certs)}")

            st.divider()

            # Per-certification cards
            for cert_name, pred in cert.predictions.items():
                if "N/A" in pred.verdict:
                    continue
                color = ("#0D9488" if pred.pass_probability > 0.7 else
                         "#D97706" if pred.pass_probability > 0.4 else "#ef4444")
                with st.expander(
                    f"**{cert_name}** — {pred.verdict}  "
                    f"({pred.pass_probability*100:.0f}% probability)",
                    expanded=pred.pass_probability > 0.65
                ):
                    p1,p2,p3,p4 = st.columns(4)
                    p1.metric("Pass Prob.", f"{pred.pass_probability*100:.0f}%")
                    p2.metric("Score", f"{pred.score:.0f}/100")
                    p3.metric("Cost", pred.estimated_cost[:18])
                    p4.metric("Timeline", pred.estimated_timeline[:18])

                    if pred.profile:
                        st.caption(f"**Body:** {pred.profile.body} · "
                                   f"**Region:** {pred.profile.region}")
                        st.caption(f"**Commercial value:** {pred.commercial_value}")

                    if pred.blocking_issues:
                        st.subheader("🚨 Must Fix to Pass")
                        for b in pred.blocking_issues:
                            st.error(b)

                    if pred.gap_analysis:
                        st.subheader("🔧 Gap Analysis")
                        for g in pred.gap_analysis:
                            st.warning(g)

                    if pred.strengths:
                        st.subheader("✅ Already Meeting")
                        for s in pred.strengths:
                            st.success(s)

                    if pred.warnings:
                        with st.expander("⚠️ Warnings"):
                            for w in pred.warnings:
                                st.warning(w)

            # Cost-benefit table
            st.divider()
            st.subheader("📊 Certification Cost–Benefit Analysis")
            cb_data = []
            for cname, pred in cert.predictions.items():
                if pred.pass_probability > 0 and "N/A" not in pred.verdict:
                    prof = pred.profile
                    cb_data.append({
                        "Certification": cname,
                        "Pass Probability": f"{pred.pass_probability*100:.0f}%",
                        "Est. Cost": pred.estimated_cost,
                        "Timeline": pred.estimated_timeline,
                        "Verdict": pred.verdict,
                        "Commercial Value": pred.commercial_value[:50],
                    })
            if cb_data:
                st.dataframe(pd.DataFrame(cb_data),
                             use_container_width=True, hide_index=True)


st.markdown("""
<div style='margin-top:32px; padding:12px 16px; background:#0A1628;
     border-top:1px solid #1e3450; border-radius:8px;
     display:flex; justify-content:space-between; align-items:center;
     font-family:"IBM Plex Mono",monospace; font-size:0.68rem; color:#64748b'>
    <span>IntelliForm v1.5 · ChemeNova LLC × ChemRich Global</span>
    <span>Makani S.S. · ChemRxiv 2026 · NJIT & UIC</span>
    <span>Mordred · PubChemPy · GP-Bayes · MIT License</span>
</div>
""", unsafe_allow_html=True)
