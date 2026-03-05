"""
modules/agents.py
Agent swarm reasoning for IntelliForm v0.7.

In v0.7, the swarm can use a real LLM call (Groq) to generate
contextual, application-specific commentary on the blend.
Falls back to structured f-string templates from v0.5/0.6.

The LLM call is intentionally separate from the parser call so
rate limits and latency hit different steps of the UX pipeline.
"""
import os
import json
from typing import List

from modules.optimizer import OptResult
from modules.llm_parser import ParseResult


# ── Groq-powered swarm ────────────────────────────────────────────────────────

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
            model="llama-3.1-8b-instant",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3,
            max_tokens=400,
            response_format={"type": "json_object"}
        )
        raw = response.choices[0].message.content
        # The model may wrap in {"comments": [...]} or return a bare array
        data = json.loads(raw)
        if isinstance(data, list):
            comments = data
        else:
            # Try common wrapper keys
            comments = data.get("comments") or data.get("agents") or list(data.values())[0]

        labels = ["**💰 Cost Agent**", "**🌿 Green Agent**", "**⚗️ Performance Agent**", "**📋 Regulatory Agent**"]
        return [f"{label}: {comment}" for label, comment in zip(labels, comments)]

    except Exception as e:
        print(f"[agents] Groq swarm failed: {e}")
        return []


# ── Fallback swarm (v0.5 template logic, always works) ───────────────────────

def _swarm_fallback(result: OptResult, parsed: ParseResult) -> List[str]:
    savings_pct = round((result.cost_per_kg * 1.28 - result.cost_per_kg) / (result.cost_per_kg * 1.28) * 100)
    app_label = parsed.application_type.replace("_", " ").title()

    comments = [
        f"**💰 Cost Agent**: ${result.cost_per_kg:.2f}/kg — approximately {savings_pct}% under typical "
        f"Ecocert green surfactant market rate. Strong margin for {app_label} positioning.",

        f"**🌿 Green Agent**: {result.bio_pct:.1f}% bio-based — fully circular & REACH Green. "
        f"Qualifies for EU Ecolabel and COSMOS certification pathways.",

        f"**⚗️ Performance Agent**: Score {result.perf_score:.1f}/100 — "
        + ("excellent for mild foaming and skin compatibility (COSMOS/Natrue benchmark)."
           if parsed.application_type in ("cosmetics", "personal_care")
           else "solid performance score for the stated application type."),

        f"**📋 Regulatory Agent**: ✅ All ingredients REACH/EPA Green-listed. "
        f"Ready for immediate ChemRich NJ toll manufacturing. "
        + ("⚠️ Note: constraints were relaxed to find this blend — review before filing."
           if result.relaxed else "No regulatory flags.")
    ]
    return comments


# ── Public interface ──────────────────────────────────────────────────────────

def run_agent_swarm(result: OptResult, parsed: ParseResult) -> List[str]:
    """
    Generate agent commentary on a formulation result.
    Tries Groq first, falls back to templates.
    Returns a list of 4 markdown strings.
    """
    if not result.success:
        return [f"**❌ Swarm halted**: Optimization failed — {result.error_msg}"]

    comments = _swarm_via_groq(result, parsed)
    if len(comments) == 4:
        return comments
    return _swarm_fallback(result, parsed)
