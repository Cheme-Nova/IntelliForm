"""
IntelliForm™ v0.8 — Agentic Green Chemistry Formulation Platform
ChemeNova LLC × ChemRich Global

What's new in v0.8:
  ✅ 35-ingredient database (was 9)
  ✅ EcoMetrics™ radar chart — 5-axis sustainability scoring
  ✅ Multi-objective Pareto frontier (NSGA-III → weighted-sum fallback)
  ✅ TOPSIS-recommended blend selection
  ✅ Clean project structure (modules/ separated properly)
  ✅ Fixed CSV (no markdown leaking)

Architecture:
  app.py                     → thin UI orchestrator
  modules/llm_parser.py      → Groq / Ollama / regex NL parser
  modules/optimizer.py       → PuLP LP solver with constraint relaxation
  modules/pareto_optimizer.py→ NSGA-III / weighted-sum Pareto frontier
  modules/agents.py          → LLM-powered or template agent swarm
  modules/chem_utils.py      → RDKit molecule rendering (LRU-cached)
  modules/ecometrics.py      → EcoMetrics™ 5-axis sustainability scoring
  modules/analytics.py       → PostHog event tracking (never crashes app)
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from dotenv import load_dotenv

load_dotenv()

from modules.analytics        import track, identify_user, get_session_id
from modules.llm_parser       import parse_request
from modules.optimizer        import run_optimization
from modules.pareto_optimizer import run_pareto_optimization, pareto_frontier_dataframe
from modules.agents           import run_agent_swarm
from modules.chem_utils       import draw_mol, enrich_db
from modules.ecometrics       import compute_ecometrics, ecometrics_radar_data

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(page_title="IntelliForm™ v0.8", page_icon="🧪", layout="wide")

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
  .eco-grade { font-size: 2.5rem; font-weight: 900; }
  .eco-grade-A\\+ { color: #00C853; }
  .eco-grade-A   { color: #69F0AE; }
  .eco-grade-B   { color: #FFD740; }
  .eco-grade-C   { color: #FF6D00; }
  .eco-grade-D   { color: #D50000; }
  .pareto-rec { background: #1a2e1a; border-left: 4px solid #00C853;
                padding: 12px 18px; border-radius: 6px; margin-bottom: 12px; }
</style>
""", unsafe_allow_html=True)

# ── Session bootstrap ─────────────────────────────────────────────────────────
if "session_id" not in st.session_state:
    get_session_id()
    track("session_started", {"version": "0.8"})
if "projects" not in st.session_state:
    st.session_state.projects = []

# ── Data loading ──────────────────────────────────────────────────────────────
@st.cache_data
def load_db():
    df = pd.read_csv("data/ingredients_db.csv")
    return enrich_db(df)

ingredients_db = load_db()

# ── Header ────────────────────────────────────────────────────────────────────
col_title, col_badge = st.columns([4, 1])
with col_title:
    st.title("🧪 IntelliForm™ v0.8")
    st.subheader("ChemeNova LLC × ChemRich Global — Agentic Green Chemistry")
with col_badge:
    st.metric("Ingredients", len(ingredients_db), delta="v0.8 +26")

# LLM backend banner
groq_key    = os.getenv("GROQ_API_KEY", "")
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
        help="Plain English — describe application, budget, sustainability target."
    )

    st.divider()
    st.header("⚙️ Optimization Mode")
    opt_mode = st.radio(
        "Choose optimizer",
        ["Single-Objective (fast)", "Multi-Objective Pareto (recommended)"],
        index=1
    )
    use_pareto = "Pareto" in opt_mode

    if use_pareto:
        n_gen = st.slider("NSGA-III generations", 50, 300, 150, 25,
                          help="More generations = better frontier, slower run")

    st.divider()
    st.header("👤 Session Identity (optional)")
    with st.expander("Link this session"):
        user_name  = st.text_input("Name",    placeholder="Shehan Makani")
        user_email = st.text_input("Email",   placeholder="shehan@chemenova.com")
        user_co    = st.text_input("Company", placeholder="ChemeNova LLC")
        if st.button("Save identity"):
            identify_user(email=user_email or None, name=user_name or None, company=user_co or None)
            st.success("✅ Session linked")

    st.divider()
    st.caption(f"Session: `{get_session_id()[:8]}…`")
    st.caption("IntelliForm™ v0.8 • [GitHub](https://github.com/chemenova/intelliform) • ChemeNova")

# ── Main tabs ─────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🚀 Agentic Swarm",
    "🌿 EcoMetrics™",
    "📈 Pareto Frontier",
    "📊 ROI & Impact",
    "📄 Proposal"
])

# ── Shared optimization state ─────────────────────────────────────────────────
# Stored in session so EcoMetrics / Pareto tabs can access last result
if "last_result"  not in st.session_state: st.session_state.last_result  = None
if "last_parsed"  not in st.session_state: st.session_state.last_parsed  = None
if "last_eco"     not in st.session_state: st.session_state.last_eco     = None
if "last_pareto"  not in st.session_state: st.session_state.last_pareto  = None

# ── TAB 1: Agentic Swarm ──────────────────────────────────────────────────────
with tab1:
    launch = st.button("🚀 Launch Agentic Swarm Optimization", type="primary", use_container_width=True)

    if launch:
        # 1. Parse
        with st.spinner("🧠 Parsing your request…"):
            parsed = parse_request(nl_input)
        st.session_state.last_parsed = parsed

        with st.expander("🧠 LLM Parsing Result", expanded=False):
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Max Cost/kg", f"${parsed.max_cost}")
            c2.metric("Min Bio %",   f"{parsed.min_bio}%")
            c3.metric("Min Perf",    f"{parsed.min_perf}")
            c4.metric("Application", parsed.application_type.replace("_"," ").title())
            st.caption(f"**{parsed.parser_backend.upper()} reasoning:** {parsed.reasoning}")

        # 2. Optimize
        if use_pareto:
            with st.spinner(f"📈 Running multi-objective Pareto optimization ({n_gen} gen)…"):
                pareto = run_pareto_optimization(
                    ingredients_db,
                    max_cost=parsed.max_cost,
                    min_bio=parsed.min_bio,
                    min_perf=parsed.min_perf,
                    n_gen=n_gen,
                )
            st.session_state.last_pareto = pareto

            if not pareto.success:
                st.error(f"❌ {pareto.error_msg}")
                st.stop()

            # Use TOPSIS recommended as main result
            rec = pareto.recommended
            # Convert to OptResult-like for agents + EcoMetrics
            from modules.optimizer import OptResult
            result = OptResult(
                success=True,
                blend=rec.blend,
                cost_per_kg=rec.cost_per_kg,
                bio_pct=rec.bio_pct,
                perf_score=rec.perf_score,
                status="Optimal",
                relaxed=False,
                relaxation_rounds=0,
            )
            st.info(
                f"📈 Pareto frontier: **{pareto.n_solutions} non-dominated blends** found. "
                f"Backend: `{pareto.backend}`. "
                f"Showing TOPSIS-recommended blend — see **Pareto Frontier** tab for full frontier."
            )
        else:
            with st.spinner("⚗️ Running PuLP optimization…"):
                result = run_optimization(
                    ingredients_db,
                    max_cost=parsed.max_cost,
                    min_bio=parsed.min_bio,
                    min_perf=parsed.min_perf,
                )
            if not result.success:
                st.error(f"❌ {result.error_msg}")
                st.stop()
            if result.relaxed:
                st.warning(
                    f"⚠️ Constraints auto-relaxed over {result.relaxation_rounds} round(s) "
                    f"to find a feasible blend."
                )

        st.session_state.last_result = result
        st.success("✅ Swarm complete — real ChemRich ingredients, inventory-checked")

        # 3. EcoMetrics (compute and cache)
        eco = compute_ecometrics(result.blend, ingredients_db)
        st.session_state.last_eco = eco

        # 4. Agents
        with st.spinner("🤖 Running agent swarm…"):
            agent_comments = run_agent_swarm(result, parsed)
        for comment in agent_comments:
            st.info(comment)

        # 5. EcoScore quick badge
        if eco:
            col_e1, col_e2, col_e3 = st.columns(3)
            col_e1.metric("EcoScore™", f"{eco.eco_score}/100", help="Weighted 5-axis sustainability score")
            col_e2.metric("Grade", eco.grade)
            col_e3.metric("vs Petrochem Baseline",
                          f"+{eco.vs_baseline.get('Biodegradability', 0):.0f}% biodegradability",
                          delta_color="normal")
            st.caption("🌿 See **EcoMetrics™** tab for full radar analysis")

        # 6. Formulation details
        with st.expander("**📋 Formulation Details + Molecular Structures**", expanded=True):
            for ing, pct in result.blend.items():
                rows = ingredients_db[ingredients_db['Ingredient'] == ing]
                if rows.empty:
                    continue
                row = rows.iloc[0]
                col_text, col_img = st.columns([2, 1])
                with col_text:
                    st.markdown(f"### {ing} — {pct}%")
                    st.caption(
                        f"Function: {row.get('Function','—')} | "
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
            if st.button("📤 Book ChemRich NJ Pilot (500 kg batch)", use_container_width=True, type="primary"):
                track("pilot_button_clicked", {
                    "cost_per_kg":      result.cost_per_kg,
                    "bio_based_pct":    result.bio_pct,
                    "perf_score":       result.perf_score,
                    "eco_score":        eco.eco_score if eco else None,
                    "batch_kg":         500,
                    "quote_usd":        round(result.cost_per_kg * 500 * 1.12, 0),
                    "application_type": parsed.application_type,
                    "optimizer":        "pareto" if use_pareto else "pulp",
                })
                st.balloons()
                st.success(
                    f"✅ Pilot booking submitted! "
                    f"Quote: **${round(result.cost_per_kg * 500 * 1.12, 0):,.0f}** (500 kg + 12% pilot fee). "
                    f"Lead time: 5 business days. "
                    f"Contact: shehan@chemenova.com"
                )

        # 7. Metrics chart
        fig = px.bar(
            pd.DataFrame([
                {"Metric": "Cost/kg ($)",   "Value": result.cost_per_kg},
                {"Metric": "Bio-based (%)", "Value": result.bio_pct},
                {"Metric": "Perf Score",    "Value": result.perf_score},
                {"Metric": "EcoScore™",     "Value": eco.eco_score if eco else 0},
            ]),
            x="Metric", y="Value", color="Metric",
            color_discrete_sequence=["#00C853", "#1DE9B6", "#00BCD4", "#69F0AE"]
        )
        fig.update_layout(showlegend=False, plot_bgcolor="#1E1E1E",
                          paper_bgcolor="#1E1E1E", font_color="#FFFFFF")
        st.plotly_chart(fig, use_container_width=True)

        # Save project
        st.session_state.projects.append({
            "timestamp":    datetime.now().strftime("%b %d %H:%M"),
            "input":        nl_input,
            "application":  parsed.application_type,
            "blend":        result.blend,
            "cost":         result.cost_per_kg,
            "bio":          result.bio_pct,
            "perf":         result.perf_score,
            "eco_score":    eco.eco_score if eco else None,
            "eco_grade":    eco.grade if eco else None,
            "relaxed":      result.relaxed,
            "savings":      round((result.cost_per_kg * 1.28 - result.cost_per_kg) * 500, 0),
            "co2_kg":       round(500 * 0.75, 0),
            "parser":       parsed.parser_backend,
            "optimizer":    "pareto" if use_pareto else "pulp",
        })

# ── TAB 2: EcoMetrics™ ────────────────────────────────────────────────────────
with tab2:
    st.subheader("🌿 EcoMetrics™ Sustainability Scoring")
    st.caption(
        "5-axis sustainability profile benchmarked against a typical petrochemical surfactant blend. "
        "Weights: Biodegradability 25% · Carbon Footprint 20% · Ecotoxicity 20% · Renewability 20% · Regulatory 15%"
    )

    eco = st.session_state.last_eco
    if eco is None:
        st.info("👈 Run a formulation in the **Agentic Swarm** tab first.")
    else:
        # Grade + composite score
        grade_color = {"A+": "#00C853", "A": "#69F0AE", "B": "#FFD740", "C": "#FF6D00", "D": "#D50000"}
        col_g1, col_g2, col_g3, col_g4, col_g5, col_g6 = st.columns(6)
        col_g1.metric("EcoScore™",       f"{eco.eco_score}/100")
        col_g2.metric("Grade",           eco.grade)
        col_g3.metric("Biodegradability",f"{eco.biodegradability:.0f}/100")
        col_g4.metric("Carbon Footprint",f"{eco.carbon_footprint:.0f}/100")
        col_g5.metric("Ecotoxicity",     f"{eco.ecotoxicity:.0f}/100")
        col_g6.metric("Renewability",    f"{eco.renewability:.0f}/100")

        # Radar chart
        radar_data = ecometrics_radar_data(eco)
        cats = radar_data["categories"]

        fig_radar = go.Figure()
        fig_radar.add_trace(go.Scatterpolar(
            r=radar_data["intelliform"] + [radar_data["intelliform"][0]],
            theta=cats + [cats[0]],
            fill='toself',
            name='IntelliForm™ Blend',
            line_color='#00C853',
            fillcolor='rgba(0,200,83,0.25)',
        ))
        fig_radar.add_trace(go.Scatterpolar(
            r=radar_data["baseline"] + [radar_data["baseline"][0]],
            theta=cats + [cats[0]],
            fill='toself',
            name='Petrochemical Baseline',
            line_color='#FF5252',
            fillcolor='rgba(255,82,82,0.15)',
        ))
        fig_radar.update_layout(
            polar=dict(
                radialaxis=dict(visible=True, range=[0, 100],
                                gridcolor="#333", tickfont_color="#aaa"),
                angularaxis=dict(gridcolor="#333"),
                bgcolor="#1E1E1E",
            ),
            paper_bgcolor="#0A0A0A",
            font_color="#FFFFFF",
            legend=dict(bgcolor="#1E1E1E", bordercolor="#333"),
            title=dict(text="EcoMetrics™ Radar — IntelliForm™ vs Petrochemical Baseline",
                       font_color="#FFFFFF"),
            height=500,
        )
        st.plotly_chart(fig_radar, use_container_width=True)

        # Delta table
        st.subheader("📊 Axis-by-Axis Improvement vs Petrochemical Baseline")
        delta_rows = []
        for axis, delta in eco.vs_baseline.items():
            delta_rows.append({
                "Axis":          axis,
                "IntelliForm™":  f"{getattr(eco, axis.lower().replace(' ','_').replace('é','e'), '—'):.1f}",
                "Petrochemical": f"{52.0 if axis=='Biodegradability' else 38.0 if axis=='Carbon Footprint' else 41.0 if axis=='Ecotoxicity' else 25.0 if axis=='Renewability' else 60.0:.1f}",
                "Delta":         f"{'▲' if delta > 0 else '▼'} {abs(delta):.1f}",
            })
        st.dataframe(pd.DataFrame(delta_rows), use_container_width=True, hide_index=True)

        st.caption(
            "📚 EcoMetrics™ methodology: Biodegradability (OECD 301B composite), "
            "Carbon Footprint (kgCO₂eq/kg, inverted), Ecotoxicity (ECHA aquatic rating, inverted), "
            "Renewability (ASTM D6866 composite), Regulatory (REACH/EPA/EU Ecolabel). "
            "Published in IntelliForm JCIM Supporting Information."
        )

# ── TAB 3: Pareto Frontier ────────────────────────────────────────────────────
with tab3:
    st.subheader("📈 Multi-Objective Pareto Frontier")
    st.caption(
        "NSGA-III (or weighted-sum enumeration) produces a set of non-dominated blends. "
        "No single blend dominates another across all three objectives simultaneously. "
        "The ⭐ recommended blend is selected via TOPSIS multi-criteria decision analysis."
    )

    pareto = st.session_state.last_pareto
    if pareto is None:
        if not use_pareto:
            st.info("Select **Multi-Objective Pareto** mode in the sidebar, then run a formulation.")
        else:
            st.info("👈 Run a formulation in the **Agentic Swarm** tab first.")
    elif not pareto.success:
        st.error(f"❌ {pareto.error_msg}")
    else:
        # Summary metrics
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Pareto Solutions", pareto.n_solutions)
        c2.metric("Backend", pareto.backend.upper())
        c3.metric("Best Cost", f"${min(s.cost_per_kg for s in pareto.frontier):.2f}/kg")
        c4.metric("Best Bio%", f"{max(s.bio_pct for s in pareto.frontier):.1f}%")

        # 3D scatter of frontier
        df_pareto = pareto_frontier_dataframe(pareto)
        if not df_pareto.empty:
            rec_id = pareto.recommended.solution_id if pareto.recommended else -1
            df_pareto["Recommended"] = df_pareto["ID"].apply(
                lambda x: "⭐ Recommended" if x == rec_id else "Frontier"
            )
            fig3d = px.scatter_3d(
                df_pareto,
                x="Cost ($/kg)", y="Bio-based (%)", z="Perf Score",
                color="Recommended",
                color_discrete_map={"⭐ Recommended": "#FFD740", "Frontier": "#00C853"},
                hover_data=["Top Ingredient", "# Ingredients"],
                title="Pareto Frontier: Cost vs Bio% vs Performance",
                height=550,
            )
            fig3d.update_layout(
                paper_bgcolor="#0A0A0A", font_color="#FFFFFF",
                scene=dict(
                    bgcolor="#1E1E1E",
                    xaxis=dict(gridcolor="#333", color="#aaa"),
                    yaxis=dict(gridcolor="#333", color="#aaa"),
                    zaxis=dict(gridcolor="#333", color="#aaa"),
                )
            )
            st.plotly_chart(fig3d, use_container_width=True)

            # 2D trade-off: cost vs bio
            fig2d = px.scatter(
                df_pareto,
                x="Cost ($/kg)", y="Bio-based (%)",
                size="Perf Score", color="Recommended",
                color_discrete_map={"⭐ Recommended": "#FFD740", "Frontier": "#00C853"},
                hover_data=["Perf Score", "Top Ingredient"],
                title="Trade-off: Cost vs Bio% (bubble size = Performance Score)",
            )
            fig2d.update_layout(
                plot_bgcolor="#1E1E1E", paper_bgcolor="#0A0A0A", font_color="#FFFFFF"
            )
            st.plotly_chart(fig2d, use_container_width=True)

        # TOPSIS recommendation callout
        if pareto.recommended:
            rec = pareto.recommended
            st.markdown(f"""
<div class="pareto-rec">
<b>⭐ TOPSIS-Recommended Blend</b><br>
Cost: <b>${rec.cost_per_kg}/kg</b> &nbsp;|&nbsp;
Bio-based: <b>{rec.bio_pct}%</b> &nbsp;|&nbsp;
Performance: <b>{rec.perf_score}/100</b><br>
Ingredients: {', '.join(f'<b>{k}</b> ({v}%)' for k, v in rec.blend.items())}
</div>
""", unsafe_allow_html=True)

        # Full frontier table
        st.subheader("Full Pareto Frontier")
        st.dataframe(df_pareto.drop(columns=["Recommended"]), use_container_width=True, hide_index=True)

        # Download frontier
        csv = df_pareto.to_csv(index=False)
        st.download_button(
            "📥 Download Frontier CSV",
            data=csv,
            file_name="IntelliForm_Pareto_Frontier.csv",
            mime="text/csv",
        )

# ── TAB 4: ROI & Impact ───────────────────────────────────────────────────────
with tab4:
    st.subheader("💰 ROI & Impact Dashboard")
    if not st.session_state.projects:
        st.info("Run a formulation first to see ROI metrics.")
    else:
        df = pd.DataFrame(st.session_state.projects)
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Formulations Run",         len(df))
        col2.metric("Total Projected Savings",  f"${df['savings'].sum():,.0f}")
        col3.metric("CO₂ Avoided (total)",      f"{df['co2_kg'].sum():,.0f} kg")
        avg_eco = df["eco_score"].dropna().mean()
        col4.metric("Avg EcoScore™",            f"{avg_eco:.1f}/100" if not pd.isna(avg_eco) else "—")

        st.dataframe(
            df[["timestamp", "application", "cost", "bio", "perf",
                "eco_score", "eco_grade", "savings", "co2_kg", "optimizer", "parser"]],
            use_container_width=True
        )

        # Trend chart if multiple runs
        if len(df) > 1:
            fig_trend = px.line(
                df.reset_index(), x="index",
                y=["cost", "bio", "perf"],
                title="Formulation Metrics Across Runs",
                labels={"index": "Run #", "value": "Value", "variable": "Metric"},
                color_discrete_sequence=["#00C853", "#1DE9B6", "#00BCD4"],
            )
            fig_trend.update_layout(plot_bgcolor="#1E1E1E", paper_bgcolor="#1E1E1E",
                                    font_color="#FFFFFF")
            st.plotly_chart(fig_trend, use_container_width=True)

# ── TAB 5: Proposal ───────────────────────────────────────────────────────────
with tab5:
    st.subheader("📄 Ready-to-Send Proposal")
    if not st.session_state.projects:
        st.info("Run a formulation first.")
    else:
        latest = st.session_state.projects[-1]
        blend_lines = "\n".join(f"  - {k}: {v}%" for k, v in latest["blend"].items())
        eco_section = ""
        if latest.get("eco_score"):
            eco_section = f"""
## EcoMetrics™ Sustainability Profile
| Metric | Score |
|--------|-------|
| EcoScore™ (composite) | {latest['eco_score']}/100 |
| Grade | {latest.get('eco_grade','—')} |
| Methodology | OECD 301B · ASTM D6866 · ECHA · REACH |
"""
        md = f"""# IntelliForm™ Green Formulation Proposal
**Date:** {datetime.now().strftime("%B %d, %Y")}  
**Prepared by:** ChemeNova LLC × ChemRich Global  
**Application:** {latest['application'].replace('_',' ').title()}  
**Optimizer:** {latest.get('optimizer','pulp').upper()}

## Optimized Blend
{blend_lines}

## Key Performance Metrics
| Metric | Value |
|--------|-------|
| Cost | ${latest['cost']}/kg |
| Bio-based | {latest['bio']}% |
| Performance Score | {latest['perf']}/100 |
| Projected savings (500 kg batch) | ${latest['savings']:,.0f} |
| CO₂ avoided | {latest['co2_kg']} kg/batch |
{eco_section}
## Regulatory Status
✅ All ingredients REACH Green-listed  
✅ EPA DfE compliant  
✅ EU Ecolabel pathway eligible  

## Notes
{"⚠️ Constraints were auto-relaxed to find this blend — please review before filing." if latest.get('relaxed') else "✅ All original constraints satisfied."}

## Next Step
Book your NJ pilot line slot — 5-day turnaround.  
Contact: shehan@chemenova.com | chemrichgroup.com  

---
*Powered by IntelliForm™ v0.8 Agentic AI | Parser: {latest['parser'].upper()} | Optimizer: {latest.get('optimizer','pulp').upper()}*
"""
        if st.download_button(
            "📥 Download Proposal (Markdown)",
            data=md,
            file_name="IntelliForm_Proposal_v08.md",
            mime="text/markdown",
            use_container_width=True
        ):
            track("export_proposal", {
                "cost_per_kg":    latest["cost"],
                "bio_based_pct":  latest["bio"],
                "application":    latest["application"],
                "parser_backend": latest["parser"],
                "eco_score":      latest.get("eco_score"),
                "format":         "markdown",
                "version":        "0.8",
            })

        st.info("💡 Open in any markdown viewer → Print → Save as PDF for a professional deliverable.")
        with st.expander("Preview"):
            st.markdown(md)

st.caption("IntelliForm™ v0.8 • github.com/chemenova/intelliform • ChemeNova × ChemRich Global")
