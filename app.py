"""
IntelliForm™ v0.7 — Agentic Green Chemistry Formulation Platform
ChemeNova LLC × ChemRich Global

Architecture:
  app.py → thin UI orchestrator
  modules/llm_parser.py  → Groq / Ollama / regex NL parser
  modules/optimizer.py   → PuLP LP solver with constraint relaxation
  modules/agents.py      → LLM-powered or template agent swarm
  modules/chem_utils.py  → RDKit molecule rendering (LRU-cached)
  modules/analytics.py   → PostHog event tracking (never crashes app)
"""
import os
from datetime import datetime

import pandas as pd
import plotly.express as px
import streamlit as st
from dotenv import load_dotenv

load_dotenv()  # reads .env for GROQ_API_KEY, POSTHOG_API_KEY, etc.

from modules.analytics  import track, identify_user, get_session_id
from modules.llm_parser import parse_request
from modules.optimizer  import run_optimization
from modules.agents     import run_agent_swarm
from modules.chem_utils import draw_mol, enrich_db

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(page_title="IntelliForm™ v0.7", page_icon="🧪", layout="wide")

# ── Session bootstrap ─────────────────────────────────────────────────────────
if "session_id" not in st.session_state:
    # First load — fire session_started event
    get_session_id()
    track("session_started", {"version": "0.7"})

if "projects" not in st.session_state:
    st.session_state.projects = []

# ── Data loading (cached) ─────────────────────────────────────────────────────
@st.cache_data
def load_db():
    df = pd.read_csv("data/ingredients_db.csv")
    return enrich_db(df)

ingredients_db = load_db()

# ── Header ────────────────────────────────────────────────────────────────────
st.title("🧪 IntelliForm™ v0.7")
st.subheader("ChemeNova LLC × ChemRich Global — Agentic Green Chemistry")

# Show which LLM backend is active
groq_key = os.getenv("GROQ_API_KEY", "")
ollama_host = os.getenv("OLLAMA_HOST", "")
if groq_key:
    st.success("🤖 LLM: Groq (llama-3.1-8b-instant) — real language understanding active")
elif ollama_host:
    st.info(f"🤖 LLM: Ollama @ {ollama_host} — local inference active")
else:
    st.warning("⚠️ LLM: Regex fallback mode — add GROQ_API_KEY to .env for full NL understanding")

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🗣️ Customer Request")
    nl_input = st.text_area(
        "Describe your formulation need",
        value="I need a mild foaming green surfactant for cosmetics "
              "under $4/kg, 98% bio-based, high skin compatibility",
        height=120,
        help="Plain English — describe application, budget, sustainability target, and any constraints."
    )

    st.divider()
    st.header("👤 Session Identity (optional)")
    with st.expander("Link this session to your name/email"):
        user_name  = st.text_input("Name",    placeholder="Shehan Makani")
        user_email = st.text_input("Email",   placeholder="shehan@chemenova.com")
        user_co    = st.text_input("Company", placeholder="ChemeNova LLC")
        if st.button("Save identity"):
            identify_user(email=user_email or None,
                          name=user_name   or None,
                          company=user_co  or None)
            st.success("✅ Session linked to PostHog identity")

    st.divider()
    st.caption(f"Session: `{get_session_id()[:8]}…`")
    st.caption("IntelliForm™ v0.7 • [GitHub](https://github.com/chemenova/intelliform) • ChemeNova")

# ── Main tabs ─────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4 = st.tabs([
    "🚀 Agentic Swarm",
    "📊 ROI & Impact",
    "📜 My Projects",
    "📄 Proposal"
])

# ── TAB 1: Optimization ───────────────────────────────────────────────────────
with tab1:
    if st.button("🚀 Launch Agentic Swarm Optimization", type="primary", use_container_width=True):

        # Step 1: Parse
        with st.spinner("🧠 Parsing your request…"):
            parsed = parse_request(nl_input)

        # Show LLM reasoning
        with st.expander("🧠 LLM Parsing Result", expanded=False):
            cols = st.columns(4)
            cols[0].metric("Max Cost/kg", f"${parsed.max_cost}")
            cols[1].metric("Min Bio %",   f"{parsed.min_bio}%")
            cols[2].metric("Min Perf",    f"{parsed.min_perf}")
            cols[3].metric("Application", parsed.application_type.replace("_"," ").title())
            st.caption(f"**{parsed.parser_backend.upper()} reasoning:** {parsed.reasoning}")

        # Step 2: Optimize
        with st.spinner("⚗️ Running PuLP optimization…"):
            result = run_optimization(
                ingredients_db,
                max_cost=parsed.max_cost,
                min_bio=parsed.min_bio,
                min_perf=parsed.min_perf
            )

        if not result.success:
            st.error(f"❌ {result.error_msg}")
            st.stop()

        if result.relaxed:
            st.warning(
                f"⚠️ Your original constraints were too tight. "
                f"IntelliForm auto-relaxed them over {result.relaxation_rounds} round(s) "
                f"to find a feasible green blend. Review the metrics below."
            )

        st.success("✅ Swarm complete — real ChemRich ingredients, inventory-checked")

        # Step 3: Agent commentary
        with st.spinner("🤖 Running agent swarm analysis…"):
            agent_comments = run_agent_swarm(result, parsed)

        for comment in agent_comments:
            st.info(comment)

        # Step 4: Formulation details
        with st.expander("**📋 Formulation Details + Molecular Structures**", expanded=True):
            for ing, pct in result.blend.items():
                col_text, col_img = st.columns([2, 1])
                with col_text:
                    row = ingredients_db[ingredients_db['Ingredient'] == ing].iloc[0]
                    st.markdown(f"### {ing} — {pct}%")
                    st.caption(
                        f"Cost: ${row['Cost_USD_kg']}/kg | "
                        f"Bio: {row['Bio_based_pct']}% | "
                        f"Perf: {row['Performance_Score']} | "
                        f"Stock: {row['Stock_kg']} kg | "
                        f"MW: {row.get('MW','—')} | LogP: {row.get('LogP','—')}"
                    )
                with col_img:
                    img = draw_mol(row['SMILES'])
                    if img:
                        st.image(img, width=180)

            st.divider()
            # Pilot booking
            if st.button("📤 Book ChemRich NJ Pilot (500 kg batch)", use_container_width=True, type="primary"):
                track("pilot_button_clicked", {
                    "cost_per_kg":       result.cost_per_kg,
                    "bio_based_pct":     result.bio_pct,
                    "perf_score":        result.perf_score,
                    "batch_kg":          500,
                    "quote_usd":         round(result.cost_per_kg * 500 * 1.12, 0),
                    "application_type":  parsed.application_type,
                    "relaxed":           result.relaxed
                })
                st.balloons()
                st.success(
                    f"✅ Pilot booking submitted! "
                    f"Quote: **${round(result.cost_per_kg * 500 * 1.12, 0):,.0f}** (500 kg + 12% pilot fee). "
                    f"Lead time: 5 business days. "
                    f"Contact: shehan@chemenova.com"
                )

        # Step 5: Metrics chart
        fig = px.bar(
            pd.DataFrame([
                {"Metric": "Cost/kg ($)",    "Value": result.cost_per_kg},
                {"Metric": "Bio-based (%)",  "Value": result.bio_pct},
                {"Metric": "Perf Score",     "Value": result.perf_score}
            ]),
            x="Metric", y="Value",
            color="Metric",
            color_discrete_sequence=["#00C853", "#1DE9B6", "#00BCD4"]
        )
        fig.update_layout(showlegend=False, plot_bgcolor="#1E1E1E", paper_bgcolor="#1E1E1E",
                          font_color="#FFFFFF")
        st.plotly_chart(fig, use_container_width=True)

        # Save to session
        st.session_state.projects.append({
            "timestamp":   datetime.now().strftime("%b %d %H:%M"),
            "input":       nl_input,
            "application": parsed.application_type,
            "blend":       result.blend,
            "cost":        result.cost_per_kg,
            "bio":         result.bio_pct,
            "perf":        result.perf_score,
            "relaxed":     result.relaxed,
            "savings":     round((result.cost_per_kg * 1.28 - result.cost_per_kg) * 500, 0),
            "co2_kg":      round(500 * 0.75, 0),
            "parser":      parsed.parser_backend
        })

# ── TAB 2: ROI ────────────────────────────────────────────────────────────────
with tab2:
    st.subheader("💰 ROI & Impact Dashboard")
    if not st.session_state.projects:
        st.info("Run a formulation first to see ROI metrics.")
    else:
        df = pd.DataFrame(st.session_state.projects)
        col1, col2, col3 = st.columns(3)
        col1.metric("Formulations Run",     len(df))
        col2.metric("Total Projected Savings", f"${df['savings'].sum():,.0f}")
        col3.metric("CO₂ Avoided (total)",  f"{df['co2_kg'].sum():,.0f} kg")
        st.dataframe(
            df[["timestamp", "application", "cost", "bio", "perf", "savings", "co2_kg", "parser"]],
            use_container_width=True
        )

# ── TAB 3: Projects ───────────────────────────────────────────────────────────
with tab3:
    st.subheader("📜 Saved Formulations (this session)")
    if not st.session_state.projects:
        st.info("No projects yet.")
    else:
        for p in reversed(st.session_state.projects):
            label = f"{p['timestamp']} — {p['application'].title()} — ${p['cost']}/kg"
            with st.expander(label):
                st.json(p["blend"])
                if p["relaxed"]:
                    st.warning("⚠️ Constraints were relaxed for this formulation.")

# ── TAB 4: Proposal ───────────────────────────────────────────────────────────
with tab4:
    st.subheader("📄 Ready-to-Send Proposal")
    if not st.session_state.projects:
        st.info("Run a formulation first.")
    else:
        latest = st.session_state.projects[-1]
        blend_lines = "\n".join(f"  - {k}: {v}%" for k, v in latest["blend"].items())
        md = f"""# IntelliForm™ Green Formulation Proposal
**Date:** {datetime.now().strftime("%B %d, %Y")}  
**Prepared by:** ChemeNova LLC × ChemRich Global  
**Application:** {latest['application'].replace('_',' ').title()}

## Optimized Blend
{blend_lines}

## Key Metrics
| Metric | Value |
|--------|-------|
| Cost | ${latest['cost']}/kg |
| Bio-based | {latest['bio']}% |
| Performance Score | {latest['perf']}/100 |
| Projected savings (500 kg batch) | ${latest['savings']:,.0f} |
| CO₂ avoided | {latest['co2_kg']} kg/batch |

## Notes
{"⚠️ Constraints were auto-relaxed to find this blend — please review before filing." if latest['relaxed'] else "✅ All original constraints satisfied."}

## Next Step
Book your NJ pilot line slot — 5-day turnaround.  
Contact: shehan@chemenova.com | chemrichgroup.com  

---
*Powered by IntelliForm™ v0.7 Agentic AI | Parser: {latest['parser'].upper()}*
"""
        if st.download_button(
            "📥 Download Proposal (Markdown → PDF)",
            data=md,
            file_name="IntelliForm_Proposal.md",
            mime="text/markdown",
            use_container_width=True
        ):
            track("export_proposal", {
                "cost_per_kg":    latest["cost"],
                "bio_based_pct":  latest["bio"],
                "application":    latest["application"],
                "parser_backend": latest["parser"],
                "format":         "markdown"
            })

        st.info("💡 Open in any markdown viewer → Print → Save as PDF for a professional deliverable.")
        with st.expander("Preview"):
            st.markdown(md)

st.caption("IntelliForm™ v0.7 • github.com/chemenova/intelliform • ChemeNova × ChemRich Global")
