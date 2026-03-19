"""
modules/llm_parser.py
Natural language → formulation constraints parser for IntelliForm.

Priority chain:
  1. Groq  (llama-3.1-8b-instant)
  2. Ollama local
  3. Regex fallback (always works)
"""
import os, re, json, requests
from dataclasses import dataclass
from typing import Optional
from modules.analytics import track


@dataclass
class ParseResult:
    max_cost: float
    min_bio: float
    min_perf: float
    application_type: str
    reasoning: str
    parser_backend: str
    raw_input: str


_SYSTEM_PROMPT = """You are a green chemistry formulation assistant for IntelliForm.
Extract formulation constraints from a customer request.
Respond ONLY with valid JSON — no markdown, no preamble.

JSON schema:
{
  "max_cost":        <float, USD/kg, default 5.0>,
  "min_bio":         <float, bio-based %, default 95>,
  "min_perf":        <float, performance score 0-100, default 82>,
  "application_type":<string: cosmetics|personal_care|industrial|food_safe|textile|unknown>,
  "reasoning":       <string, ≤2 sentences>
}

Rules:
- "mild" or "skin" → min_perf ≥ 85
- "foaming" or "surfactant" → min_perf ≥ 88
- "98%","99%","100%","fully bio" → min_bio ≥ 98
- Extract cost from "under $4/kg", "max $3.50", "budget $5"
- Never set min_bio < 80 or min_perf < 65"""

_USER_PROMPT = 'Customer request: """{text}"""'


def _parse_with_groq(text: str) -> Optional[ParseResult]:
    api_key = os.getenv("GROQ_API_KEY", "")
    if not api_key:
        return None
    try:
        from groq import Groq
        client = Groq(api_key=api_key)
        response = client.chat.completions.create(
            model="llama-3.1-8b-instant",
            messages=[{"role": "system", "content": _SYSTEM_PROMPT},
                      {"role": "user", "content": _USER_PROMPT.format(text=text)}],
            temperature=0.1, max_tokens=300,
            response_format={"type": "json_object"}
        )
        data = json.loads(response.choices[0].message.content)
        return ParseResult(max_cost=float(data.get("max_cost", 5.0)),
                           min_bio=float(data.get("min_bio", 95)),
                           min_perf=float(data.get("min_perf", 82)),
                           application_type=str(data.get("application_type", "unknown")),
                           reasoning=str(data.get("reasoning", "")),
                           parser_backend="groq", raw_input=text)
    except Exception as e:
        print(f"[llm_parser] Groq failed: {e}")
        return None


def _parse_with_ollama(text: str) -> Optional[ParseResult]:
    host = os.getenv("OLLAMA_HOST", "http://localhost:11434")
    model = os.getenv("OLLAMA_MODEL", "llama3")
    try:
        resp = requests.post(f"{host}/api/generate",
                             json={"model": model, "prompt": f"{_SYSTEM_PROMPT}\n\n{_USER_PROMPT.format(text=text)}", "stream": False, "format": "json"},
                             timeout=15)
        resp.raise_for_status()
        data = json.loads(resp.json().get("response", ""))
        return ParseResult(max_cost=float(data.get("max_cost", 5.0)),
                           min_bio=float(data.get("min_bio", 95)),
                           min_perf=float(data.get("min_perf", 82)),
                           application_type=str(data.get("application_type", "unknown")),
                           reasoning=str(data.get("reasoning", "")),
                           parser_backend="ollama", raw_input=text)
    except Exception as e:
        print(f"[llm_parser] Ollama failed: {e}")
        return None


def _parse_with_regex(text: str) -> ParseResult:
    t = text.lower()
    m = re.search(r'(?:under|max|budget|below)\s*\$?\s*(\d+\.?\d*)\s*(?:/kg)?', t)
    max_cost = float(m.group(1)) if m else 5.0
    min_bio = 98.0 if any(x in t for x in ['100%','99%','98%','fully bio','fully natural']) else 95.0
    if any(x in t for x in ['foaming','high foam','surfactant']): min_perf = 88.0
    elif any(x in t for x in ['skin','cosmetics','personal care','mild']): min_perf = 85.0
    else: min_perf = 82.0
    if any(x in t for x in ['cosmetics','skin','personal care','lotion','cream']): app = "cosmetics"
    elif any(x in t for x in ['food','edible','fda']): app = "food_safe"
    elif any(x in t for x in ['textile','fabric','laundry']): app = "textile"
    elif any(x in t for x in ['industrial','cleaning','degreaser']): app = "industrial"
    else: app = "unknown"
    return ParseResult(max_cost=max_cost, min_bio=min_bio, min_perf=min_perf,
                       application_type=app, reasoning="Regex fallback (no LLM API key).",
                       parser_backend="regex", raw_input=text)


def parse_request(text: str) -> ParseResult:
    result = _parse_with_groq(text) or _parse_with_ollama(text) or _parse_with_regex(text)
    track("llm_parser_used", {"parser_backend": result.parser_backend,
                               "application_type": result.application_type,
                               "max_cost": result.max_cost, "min_bio": result.min_bio,
                               "min_perf": result.min_perf, "input_length": len(text)})
    return result
