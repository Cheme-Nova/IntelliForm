# IntelliForm Upgrade Guide: v0.7 → v0.9

This document explains every file that changed or was added between the GitHub repo (v0.7) and the current codebase (v0.9).

---

## Files to REPLACE (changed significantly)

| File | What changed |
|------|-------------|
| `app.py` | Completely rewritten — 4 tabs → 7 tabs, all v0.8/v0.9 modules wired in |
| `requirements.txt` | Added: pymoo, reportlab, supabase, xgboost, scikit-learn |
| `modules/optimizer.py` | Refactored for cleaner module imports |
| `modules/agents.py` | Minor cleanup, better error handling |
| `modules/chem_utils.py` | Added graceful fallback when RDKit not installed |
| `modules/llm_parser.py` | Minor cleanup |
| `modules/analytics.py` | Minor cleanup |

---

## Files to ADD (new — did not exist in v0.7)

### v0.8 additions
| File | Description |
|------|-------------|
| `modules/ecometrics.py` | EcoMetrics™ — 5-axis sustainability scoring, radar chart, petrochemical baseline comparison |
| `modules/pareto_optimizer.py` | NSGA-III multi-objective Pareto optimizer (pymoo) with PuLP fallback and TOPSIS selection |
| `data/ingredients_db.csv` | Expanded from 9 → 35 ingredients with Function, Biodegradability, CarbonFootprint, Ecotoxicity, Renewability columns |

### v0.9 additions
| File | Description |
|------|-------------|
| `modules/qsar.py` | QSAR/QSPR models — Morgan fingerprints + XGBoost for biodegradability, ecotoxicity, performance prediction |
| `modules/regulatory.py` | Regulatory intelligence — REACH, EPA Safer Choice, EU Ecolabel, COSMOS per ingredient |
| `modules/persistence.py` | Supabase persistence with in-memory fallback |
| `modules/pdf_proposal.py` | Branded PDF proposals (ReportLab) — ChemeNova navy/teal/amber |
| `migrations/001_create_tables.sql` | Supabase schema — run once in SQL editor |

---

## Installation changes

```bash
# New dependencies for v0.8+
pip install pymoo xgboost scikit-learn

# New dependencies for v0.9+
pip install reportlab supabase

# Full install
pip install -r requirements.txt
```

---

## Environment variables added

```bash
# v0.9 — Supabase persistence (optional)
SUPABASE_URL=https://yourproject.supabase.co
SUPABASE_ANON_KEY=your-anon-key
```

---

## Breaking changes

**None.** All new modules degrade gracefully:
- No Supabase credentials → falls back to in-memory storage
- No RDKit → molecular structures hidden, optimizer still works  
- No Groq key → regex fallback parser
- No pymoo → PuLP fallback optimizer
- No ReportLab → Markdown download instead of PDF

---

## Version history

| Version | Key additions |
|---------|--------------|
| v0.5 | Initial release — PuLP optimizer, agent swarm, basic UI |
| v0.6 | Regex NL parser, PostHog analytics |
| v0.7 | Groq LLM parser, Ollama support, modular codebase |
| v0.8 | EcoMetrics™, NSGA-III Pareto optimizer, expanded ingredient DB (35 ingredients) |
| v0.9 | QSAR models, regulatory intelligence, Supabase persistence, PDF proposals, model card tab |
