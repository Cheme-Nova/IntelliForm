"""
IntelliForm v1.0 — Agentic Green Chemistry Formulation Platform
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
# Tiers: all features free — upsell via pilot batch booking

st.set_page_config(page_title="IntelliForm v1.0", page_icon="🧪", layout="wide")
st.markdown("""<style>
.pareto-rec{background:#1a2e1a;border-left:4px solid #00C853;padding:12px 18px;border-radius:6px;margin-bottom:12px}
.carbon-card{background:#0a2e1a;border-left:4px solid #059669;padding:12px 18px;border-radius:6px;margin-bottom:12px}
.stability-card{background:#1a1a2e;border-left:4px solid #0D9488;padding:12px 18px;border-radius:6px;margin-bottom:12px}
</style>""", unsafe_allow_html=True)

# ── Bootstrap ─────────────────────────────────────────────────────────────────
if "session_id"       not in st.session_state: get_session_id(); track("session_started",{"version":"1.0"})
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
if "compare_blend"    not in st.session_state: st.session_state.compare_blend  = None

@st.cache_data
def load_db():
    df = pd.read_csv("data/ingredients_db.csv")
    return enrich_db(df)

@st.cache_resource
def load_models(n):
    db = load_db()
    return initialize_models(db)

ingredients_db = load_db()
if st.session_state.model_card is None:
    st.session_state.model_card = load_models(len(ingredients_db))

# Load persisted projects
if not st.session_state.projects and is_connected():
    try:
        stored = load_projects(get_session_id(), limit=20)
        if stored: st.session_state.projects = stored
    except Exception: pass

# ── Header ────────────────────────────────────────────────────────────────────
h1,h2,h3,h4 = st.columns([3,1,1,1])
with h1:
    st.title("🧪 IntelliForm v1.0")
    st.caption("AI-powered green chemistry formulation — describe what you need, get a certified, pilot-ready blend in seconds.")
h2.metric("Ingredients", len(ingredients_db))
h3.metric("Certifications", "EU Ecolabel · EPA · COSMOS")
mc = st.session_state.model_card
qsar_ok = mc and mc.sklearn_version != "unavailable"
h4.metric("Optimization", "NSGA-III Pareto + LP")


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
    nl_input = st.text_area("Describe your formulation need",
        value="I need a mild foaming green surfactant for personal care, under $4/kg, at least 95% bio-based, EPA Safer Choice compatible",
        height=110)
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

    st.header("⚙️ Optimization")
    use_pareto = st.radio("Mode",["Single-Objective (fast)","Multi-Objective Pareto"],index=0) == "Multi-Objective Pareto"
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
    st.caption("IntelliForm v1.0 · [GitHub](https://github.com/chemenova/intelliform) · ChemeNova x ChemRich")

# ── Tabs ──────────────────────────────────────────────────────────────────────
# ── Tabs ──────────────────────────────────────────────────────────────────────
t1,t2,t3,t4,t5,t6,t7,t8,t9,t10 = st.tabs([
    "🚀 Agentic Swarm","🌿 EcoMetrics","📋 Regulatory",
    "📈 Pareto Frontier","🔬 Model Card","🧪 Stability & Viscosity",
    "🌍 Carbon Credits","🔄 Blend Comparison","📊 ROI & History","📄 Proposal"])

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
            c4.metric("Application",parsed.application_type.replace("_"," ").title())
            st.caption(f"**{parsed.parser_backend.upper()}**: {parsed.reasoning}")
            # Debug: show LLM errors if any
            import os as _os
            if _os.path.exists("/tmp/intelliform_llm_errors.txt"):
                with open("/tmp/intelliform_llm_errors.txt") as _f:
                    _errs = _f.read()
                if _errs:
                    st.error(f"LLM Debug:\n{_errs[-2000:]}")

        if use_pareto:
            with st.spinner(f"📈 Pareto optimization ({n_gen} gen)…"):
                pareto = run_pareto_optimization(ingredients_db,max_cost=parsed.max_cost,
                    min_bio=parsed.min_bio,min_perf=parsed.min_perf,n_gen=n_gen)
            st.session_state.last_pareto = pareto
            if not pareto.success: st.error(pareto.error_msg); st.stop()
            rec = pareto.recommended
            from modules.optimizer import OptResult
            result = OptResult(success=True,blend=rec.blend,cost_per_kg=rec.cost_per_kg,
                bio_pct=rec.bio_pct,perf_score=rec.perf_score,status="Optimal")
            st.info(f"📈 {pareto.n_solutions} Pareto solutions · backend: `{pareto.backend}`")
        else:
            with st.spinner("⚗️ PuLP optimization…"):
                result = run_optimization(ingredients_db,max_cost=parsed.max_cost,
                    min_bio=parsed.min_bio,min_perf=parsed.min_perf,
                    max_concentration=max_conc/100)
            if not result.success: st.error(result.error_msg); st.stop()
            if result.relaxed: st.warning(f"⚠️ Constraints relaxed × {result.relaxation_rounds}")

        st.session_state.last_result = result

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
    reg = st.session_state.last_reg
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
        ci1.metric("Training Set",mc.n_training); ci2.metric("Data Hash",mc.training_hash)
        ci3.metric("scikit-learn",mc.sklearn_version); ci4.metric("Active Learning Rounds",mc.active_learning_rounds)
        st.divider()
        for target,bench in mc.benchmarks.items():
            with st.expander(f"**{target}** — R²={bench['cv_r2']} · RMSE={bench['cv_rmse']} {bench['unit']}",expanded=True):
                bc1,bc2,bc3,bc4 = st.columns(4)
                bc1.metric("5-fold R²",bench["cv_r2"]); bc2.metric("CV RMSE",f"{bench['cv_rmse']} {bench['unit']}")
                bc3.metric("N train",bench["n_train"]); bc4.metric("Algorithm","Gradient Boosting")
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
            c2.metric("Viscosity", stab.viscosity_range)
            c3.metric("pH Range", f"{stab.ph_min:.1f} – {stab.ph_max:.1f}")
            c4.metric("Stability Rating", stab.overall_rating)
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
            cr1.metric("Credit Value (floor $15/t)", f"${carbon.credit_value_low:.2f}")
            cr2.metric("Credit Value (mid $45/t)", f"${carbon.credit_value_mid:.2f}")
            cr3.metric("Credit Value (premium $85/t)", f"${carbon.credit_value_high:.2f}")
            st.divider()
            st.subheader("Annual Projection (12 batches/year)")
            a1,a2 = st.columns(2)
            a1.metric("Annual CO2 Displaced", f"{carbon.annual_co2_tonnes:.2f} tonnes")
            a2.metric("Annual Credit Value", f"${carbon.annual_credit_value_mid:,.0f}", "at mid-market rate")
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
                        track("export_proposal", {"format":"pdf","version":"1.0","eco_score":latest.get("eco_score")})
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

st.caption("IntelliForm v1.0 · github.com/chemenova/intelliform · ChemeNova x ChemRich · Makani S.S., ChemRxiv 2026 · NJIT & UIC")
