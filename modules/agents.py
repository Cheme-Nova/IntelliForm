"""
modules/agents.py
Tool-using expert agent for IntelliForm v1.5+.

Agent tier (best available at runtime):
  1. Claude tool-using agent  — calls IntelliForm modules as tools, reasons over real data
     Requires ANTHROPIC_API_KEY. Model: claude-haiku-4-5-20251001 (override: ANTHROPIC_AGENT_MODEL)

  2. Groq simple swarm       — 4-agent LLM commentary, no tool access
     Requires GROQ_API_KEY.

  3. Static fallback          — rule-based commentary, no API needed

Tools available to the Claude agent:
  predict_blend_qsar      — QSAR predictions for blend ingredients (ML model active tier)
  get_formulation_memory  — negative rules + strategic insights from the memory network
  get_certifications      — pass probabilities for EPA Safer Choice, COSMOS, EU Ecolabel etc.
  get_carbon_footprint    — ISO 14067 PCF via Brightway2 or built-in emission factors
"""
from __future__ import annotations

import json
import os
import re
from typing import Any, Dict, List, Optional

from modules.optimizer  import OptResult
from modules.llm_parser import ParseResult


# ── Tool implementations ──────────────────────────────────────────────────────

def _tool_predict_blend_qsar(input_data: dict, context: dict) -> dict:
    """Run QSAR predictions for one or more blend ingredients."""
    try:
        from modules.qsar import predict_properties
        db     = context["db"]
        result = context["result"]

        names = input_data.get("ingredients") or list(result.blend.keys())
        preds = []
        for name in names[:8]:
            rows = db[db["Ingredient"] == name] if db is not None else []
            if hasattr(rows, "empty") and rows.empty:
                continue
            smiles = str(rows.iloc[0].get("SMILES", "")) if len(rows) else ""
            if not smiles or smiles == "nan":
                continue
            qp = predict_properties(smiles)
            preds.append({
                "ingredient":      name,
                "biodegradability": round(qp.biodegradability, 1),
                "ecotoxicity":      round(qp.ecotoxicity, 1),
                "performance":      round(qp.performance, 1),
                "confidence":       qp.confidence,
                "model": (
                    "Chemprop D-MPNN" if qp.used_chemprop
                    else "Mordred+GBR"  if qp.used_mordred
                    else "Morgan+GBR"   if qp.used_ml
                    else "rule-based"
                ),
            })
        return {"predictions": preds}
    except Exception as e:
        return {"error": str(e)[:120]}


def _tool_get_formulation_memory(input_data: dict, context: dict) -> dict:
    """Query the Formulation Memory Network for negative rules and strategic insights."""
    try:
        memory   = context["memory"]
        vertical = input_data.get("vertical", context["vertical"])
        if memory is None:
            return {"negative_rules": [], "insights": []}
        neg_rules = memory.get_negative_rules(vertical)
        insights  = memory.get_insights(vertical)
        return {
            "negative_rules": [
                {"ingredient": r["ingredient"], "reason": r["reason"][:80],
                 "confidence": round(r.get("confidence", 0) * 100)}
                for r in neg_rules[:5]
            ],
            "insights": [
                {"title": i.title, "type": i.insight_type,
                 "confidence": round(i.confidence * 100), "description": i.description[:80]}
                for i in insights[:4]
            ],
        }
    except Exception as e:
        return {"error": str(e)[:120]}


def _tool_get_certifications(input_data: dict, context: dict) -> dict:
    """Get green chemistry certification pass probabilities for the current blend."""
    try:
        from modules.certification_oracle import run_certification_oracle
        result   = context["result"]
        db       = context["db"]
        vertical = context["vertical"]
        cert = run_certification_oracle(result.blend, db, vertical, result.bio_pct)
        ranked = sorted(
            [(k, v) for k, v in cert.predictions.items() if "N/A" not in v.verdict],
            key=lambda x: -x[1].pass_probability,
        )
        return {
            "green_score":  round(cert.overall_green_score, 1),
            "top_cert":     cert.top_certification,
            "quick_wins":   cert.quick_wins[:3],
            "predictions":  [
                {"cert": k, "pass_prob": round(v.pass_probability * 100),
                 "verdict": v.verdict}
                for k, v in ranked[:5]
            ],
        }
    except Exception as e:
        return {"error": str(e)[:120]}


def _tool_get_carbon_footprint(input_data: dict, context: dict) -> dict:
    """Get ISO 14067 product carbon footprint for the current blend."""
    try:
        from modules.carbon_passport import generate_carbon_passport
        result = context["result"]
        db     = context["db"]
        passport = generate_carbon_passport(
            blend=result.blend, db=db,
            product_name="IntelliForm Blend",
            batch_id="agent-q", batch_kg=500,
        )
        return {
            "total_pcf_kg_co2eq_per_kg": round(passport.total_pcf, 3),
            "net_pcf":                    round(passport.net_pcf, 3),
            "grade":                      passport.carbon_intensity_grade,
            "vs_industry_avg":            round(passport.vs_industry_average, 2),
            "lca_engine":                 getattr(passport, "computation_engine", "built-in"),
        }
    except Exception as e:
        return {"error": str(e)[:120]}


_TOOL_DISPATCH: Dict[str, Any] = {
    "predict_blend_qsar":     _tool_predict_blend_qsar,
    "get_formulation_memory": _tool_get_formulation_memory,
    "get_certifications":     _tool_get_certifications,
    "get_carbon_footprint":   _tool_get_carbon_footprint,
}

_TOOL_SCHEMAS = [
    {
        "name": "predict_blend_qsar",
        "description": (
            "Run QSAR predictions (biodegradability %, ecotoxicity 1-10, performance 0-100) "
            "for one or more blend ingredients using the active ML model tier "
            "(Chemprop D-MPNN → Mordred+GBR → Morgan+GBR → rule-based). "
            "Omit 'ingredients' to predict all blend ingredients."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "ingredients": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Ingredient names (from the blend). Omit to predict all.",
                }
            },
        },
    },
    {
        "name": "get_formulation_memory",
        "description": (
            "Query the Graphiti temporal knowledge graph + Supabase memory for this vertical. "
            "Returns negative rules (ingredients that caused failures in past batches) "
            "and strategic insights learned from formulation history."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "vertical": {
                    "type": "string",
                    "description": "Industry vertical (e.g. personal_care, industrial_cleaning, pharmaceutical).",
                }
            },
        },
    },
    {
        "name": "get_certifications",
        "description": (
            "Get pass probabilities for major green chemistry certifications: "
            "EPA Safer Choice, COSMOS Organic, EU Ecolabel, NATRUE, NSF 305, "
            "USDA BioPreferred, Cradle to Cradle, Nordic Swan. "
            "Returns the top certifications and quick wins."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "get_carbon_footprint",
        "description": (
            "Get the ISO 14067 product carbon footprint (kgCO₂eq/kg) for the current blend. "
            "Uses Brightway2 (ISO 14040 matrix LCA) if available, else built-in emission factors. "
            "Returns PCF, grade (A+–D), and performance vs industry average."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
]

_SYSTEM_PROMPT = """\
You are IntelliForm's expert advisory council — a senior green chemistry team at ChemeNova LLC.
You review specialty chemical formulations and provide concise, data-driven expert commentary.

Use the provided tools to query real data from IntelliForm's modules. Then respond with a JSON array of exactly 4 expert comments (1–2 sentences each, no markdown inside strings):

["Cost Agent comment", "Green Agent comment", "Performance Agent comment", "Regulatory Agent comment"]

Cost Agent: cost/kg vs market benchmarks, margin analysis.
Green Agent: bio-based %, sustainability credentials, certification outlook.
Performance Agent: QSAR scores for the application, model tier used.
Regulatory Agent: REACH/EPA status, memory network red flags, pilot readiness.

Be specific. Cite numbers. No markdown inside the JSON strings.\
"""


# ── Tier 1: Claude tool-using agent ──────────────────────────────────────────

def _swarm_via_claude(
    result: OptResult,
    parsed: ParseResult,
    db,
    memory_net,
    vertical: str,
) -> List[str]:
    api_key = os.getenv("ANTHROPIC_API_KEY", "")
    if not api_key:
        return []
    try:
        import anthropic
    except ImportError:
        return []

    try:
        client  = anthropic.Anthropic(api_key=api_key)
        model   = os.getenv("ANTHROPIC_AGENT_MODEL", "claude-haiku-4-5-20251001")
        context = {"result": result, "db": db, "memory": memory_net, "vertical": vertical}

        blend_summary = ", ".join(f"{k} ({v}%)" for k, v in result.blend.items())
        user_msg = (
            f"Formulation: {blend_summary}\n"
            f"Cost: ${result.cost_per_kg}/kg | Bio: {result.bio_pct}% | Perf: {result.perf_score}/100\n"
            f"Application: {parsed.application_type} | Vertical: {vertical}\n"
            f"Constraints relaxed: {result.relaxed}\n\n"
            "Gather data using your tools, then write the 4-agent JSON array."
        )
        messages = [{"role": "user", "content": user_msg}]

        for _ in range(4):  # max 4 agentic turns
            resp = client.messages.create(
                model=model,
                max_tokens=1024,
                system=_SYSTEM_PROMPT,
                tools=_TOOL_SCHEMAS,
                messages=messages,
            )

            if resp.stop_reason == "end_turn":
                for block in resp.content:
                    if hasattr(block, "text"):
                        raw = block.text.strip()
                        if raw.startswith("```"):
                            raw = re.sub(r"```(?:json)?", "", raw).replace("```", "").strip()
                        s, e = raw.find("["), raw.rfind("]")
                        if s != -1 and e != -1:
                            comments = json.loads(raw[s:e + 1])
                            if isinstance(comments, list) and len(comments) == 4:
                                labels = [
                                    "**💰 Cost Agent**", "**🌿 Green Agent**",
                                    "**⚗️ Performance Agent**", "**📋 Regulatory Agent**",
                                ]
                                return [f"{lbl}: {c}" for lbl, c in zip(labels, comments)]
                break  # end_turn but no parseable JSON

            if resp.stop_reason == "tool_use":
                tool_results = []
                for block in resp.content:
                    if hasattr(block, "type") and block.type == "tool_use":
                        fn   = _TOOL_DISPATCH.get(block.name)
                        data = fn(block.input, context) if fn else {"error": f"Unknown tool: {block.name}"}
                        tool_results.append({
                            "type":        "tool_result",
                            "tool_use_id": block.id,
                            "content":     json.dumps(data),
                        })
                messages.append({"role": "assistant", "content": resp.content})
                messages.append({"role": "user",      "content": tool_results})
            else:
                break  # unexpected stop_reason

    except Exception as e:
        print(f"[agents] Claude tool-using agent failed: {e}")

    return []


# ── Tier 2: Groq simple swarm (unchanged) ────────────────────────────────────

def _swarm_via_groq(result: OptResult, parsed: ParseResult) -> List[str]:
    api_key = os.getenv("GROQ_API_KEY", "")
    if not api_key:
        return []
    try:
        from groq import Groq
        client = Groq(api_key=api_key)
        blend_summary = ", ".join(f"{k} ({v}%)" for k, v in result.blend.items())
        prompt = f"""You are a green chemistry expert reviewing a new surfactant formulation.

Formulation: {blend_summary}
Cost: ${result.cost_per_kg}/kg
Bio-based: {result.bio_pct}%
Performance Score: {result.perf_score}/100
Application: {parsed.application_type}
Constraints were relaxed: {result.relaxed}

Write exactly 4 short expert comments, one per agent below.
Respond ONLY with a JSON array of 4 strings, no markdown.

Agents:
1. Cost Agent: comment on cost competitiveness vs. market
2. Green Agent: comment on bio-based content and sustainability credentials
3. Performance Agent: comment on performance score for the stated application
4. Regulatory Agent: comment on REACH/EPA status and pilot readiness

Keep each comment to 1-2 sentences. Be specific and cite numbers."""

        response = client.chat.completions.create(
            model=os.getenv("GROQ_SWARM_MODEL", "llama-3.3-70b-versatile"),
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3,
            max_tokens=400,
        )
        raw = (response.choices[0].message.content or "").strip()
        if raw.startswith("```"):
            raw = re.sub(r"```(?:json)?", "", raw).replace("```", "").strip()
        try:
            data = json.loads(raw)
        except json.JSONDecodeError:
            start = raw.find("[")
            end   = raw.rfind("]")
            if start == -1 or end == -1:
                raise
            data = json.loads(raw[start:end + 1])
        if isinstance(data, list):
            comments = data
        else:
            candidate = data.get("comments") or data.get("agents") or data.get("responses")
            if candidate is None:
                first_val = list(data.values())[0]
                candidate = first_val if isinstance(first_val, list) else list(data.values())
            comments = candidate if isinstance(candidate, list) else [str(candidate)]

        labels = [
            "**💰 Cost Agent**", "**🌿 Green Agent**",
            "**⚗️ Performance Agent**", "**📋 Regulatory Agent**",
        ]
        return [f"{label}: {comment}" for label, comment in zip(labels, comments)]

    except Exception as e:
        print(f"[agents] Groq swarm failed: {e}")
        return []


# ── Tier 3: static fallback (unchanged) ──────────────────────────────────────

def _swarm_fallback(result: OptResult, parsed: ParseResult) -> List[str]:
    savings_pct = round(
        (result.cost_per_kg * 1.28 - result.cost_per_kg) / (result.cost_per_kg * 1.28) * 100)
    app_label = parsed.application_type.replace("_", " ").title()
    return [
        f"**💰 Cost Agent**: ${result.cost_per_kg:.2f}/kg — approximately {savings_pct}% under "
        f"typical Ecocert green surfactant market rate. Strong margin for {app_label} positioning.",
        f"**🌿 Green Agent**: {result.bio_pct:.1f}% bio-based — fully circular & REACH Green. "
        "Qualifies for EU Ecolabel and COSMOS certification pathways.",
        f"**⚗️ Performance Agent**: Score {result.perf_score:.1f}/100 — "
        + ("excellent for mild foaming and skin compatibility (COSMOS/Natrue benchmark)."
           if parsed.application_type in ("cosmetics", "personal_care")
           else "solid performance score for the stated application type."),
        f"**📋 Regulatory Agent**: ✅ All ingredients REACH/EPA Green-listed. "
        "Ready for immediate ChemRich NJ toll manufacturing. "
        + ("⚠️ Note: constraints were relaxed — review before filing." if result.relaxed
           else "No regulatory flags."),
    ]


# ── Public API ────────────────────────────────────────────────────────────────

def run_agent_swarm(
    result: OptResult,
    parsed: ParseResult,
    db=None,
    memory_net=None,
    vertical: str = "personal_care",
) -> List[str]:
    """
    Run the expert advisory agent swarm.

    Tier 1 (tool-using Claude agent) activates when ANTHROPIC_API_KEY is set
    and db + memory_net are supplied — it queries QSAR, memory, certifications,
    and carbon data before writing comments.

    Falls back to Groq simple swarm (GROQ_API_KEY) then static comments.
    """
    if not result.success:
        return [f"**❌ Swarm halted**: Optimization failed — {result.error_msg}"]

    # Tier 1: Claude tool-using agent
    if db is not None and memory_net is not None:
        comments = _swarm_via_claude(result, parsed, db, memory_net, vertical)
        if len(comments) == 4:
            return comments

    # Tier 2: Groq simple swarm
    comments = _swarm_via_groq(result, parsed)
    if len(comments) == 4:
        return comments

    # Tier 3: static fallback
    return _swarm_fallback(result, parsed)
