# Contributing to IntelliForm™

**IntelliForm is open-source green chemistry infrastructure.** We welcome chemists, ML engineers, sustainability researchers, and developers. Every contribution — from correcting a SMILES string to adding a new optimization objective — makes sustainable formulation more accessible.

---

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Ways to Contribute](#ways-to-contribute)
- [Development Setup](#development-setup)
- [Project Structure](#project-structure)
- [Submitting Changes](#submitting-changes)
- [Ingredient Data Guidelines](#ingredient-data-guidelines)
- [Testing](#testing)
- [Style Guide](#style-guide)

---

## Code of Conduct

Be kind. Be precise. Cite your sources (especially for chemical data). We're building tools for a sustainable future — let's work that way too.

---

## Ways to Contribute

### 🧪 Chemistry Contributions (no coding needed)
- **Add ingredients** to `data/ingredients_db.csv` — glucosides, biosurfactants, esters, humectants
- **Fix SMILES strings** — verify against PubChem or ChemSpider
- **Add performance data** — lab-tested bio%, cost ranges, application suitability
- **Flag regulatory issues** — REACH restrictions, EPA DfE status, EU Ecolabel conflicts

### 🤖 ML / AI Contributions
- Improve the LLM parser prompts (`modules/llm_parser.py`)
- Add new optimization objectives (e.g. HLB targeting, foam index)
- Add alternative LLM backends (OpenAI, Anthropic, Mistral)
- Build a retrieval-augmented ingredient lookup

### 🛠️ Engineering Contributions
- Add unit tests (`tests/`)
- Improve PuLP model with multi-objective optimization
- Add API endpoint (FastAPI wrapper around optimizer)
- Add persistent storage (SQLite or Supabase for project history)
- Docker / docker-compose setup

### 📊 Data & Analytics
- PostHog dashboard improvements
- Add funnel analysis for pilot conversion
- Add A/B test infrastructure via PostHog feature flags

---

## Development Setup

```bash
# 1. Clone
git clone https://github.com/chemenova/intelliform.git
cd intelliform

# 2. Create virtual environment
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Configure environment
cp .env.example .env
# Edit .env and add GROQ_API_KEY (free at console.groq.com)

# 5. Run
streamlit run app.py
```

The app runs without a Groq key — it falls back to regex parsing.

---

## Project Structure

```
intelliform/
├── app.py                  # UI layer — Streamlit tabs and layout only
├── modules/
│   ├── llm_parser.py       # NL request → constraints (Groq/Ollama/regex)
│   ├── optimizer.py        # PuLP LP solver + constraint relaxation
│   ├── agents.py           # Agent swarm commentary
│   ├── chem_utils.py       # RDKit helpers (mol rendering, property calc)
│   └── analytics.py        # PostHog wrapper
├── data/
│   └── ingredients_db.csv  # Ingredient database (community-maintained)
├── tests/                  # pytest tests (help us expand coverage!)
└── .github/
    └── ISSUE_TEMPLATE/     # Structured issue forms
```

**Key rule:** `app.py` should contain only UI code. Business logic goes in `modules/`.

---

## Submitting Changes

1. **Fork** the repository
2. **Create a branch**: `git checkout -b feat/add-lauryl-betaine` or `fix/smiles-decyl-glucoside`
3. **Make your changes** (see guidelines below)
4. **Test**: `pytest tests/` and run the app locally
5. **Open a Pull Request** — use the PR template

### Branch naming conventions
- `feat/` — new features
- `fix/` — bug fixes
- `data/` — ingredient data additions or corrections
- `docs/` — documentation only

---

## Ingredient Data Guidelines

When adding rows to `data/ingredients_db.csv`:

| Column | Rules |
|--------|-------|
| `Ingredient` | INCI name preferred. Include "(bio)" if bio-derived variant. |
| `SMILES` | Verify on [PubChem](https://pubchem.ncbi.nlm.nih.gov) or [ChemSpider](https://www.chemspider.com). Use canonical SMILES. |
| `Cost_USD_kg` | Use publicly available bulk market price (≥1 MT). Add source in PR description. |
| `Bio_based_pct` | ASTM D6866 or supplier SDS value. 0–100. |
| `Performance_Score` | 0–100 composite score (foaming 40% + mildness 40% + compatibility 20%). Document how you scored it. |
| `Stock_kg` | ChemRich NJ current stock or 0 if unknown. |
| `REACH_Flag` | Green / Amber / Red based on ECHA database status. |

**Do not add** petroleum-derived ingredients, ingredients on SVHC lists, or ingredients with Amber/Red REACH status without flagging them clearly.

---

## Testing

```bash
pytest tests/ -v
```

We need tests for:
- `modules/optimizer.py` — feasibility, infeasibility, relaxation logic
- `modules/llm_parser.py` — regex fallback parsing
- `modules/chem_utils.py` — SMILES validity, mol property ranges
- CSV data integrity (valid SMILES, numeric ranges)

If you add a feature, add a test.

---

## Style Guide

- **Python**: PEP8, `ruff` for linting (`ruff check .`)
- **Docstrings**: Google-style
- **Type hints**: required for all new `modules/` functions
- **No hardcoded secrets** — use `os.getenv()` and `.env.example`
- **Analytics**: all new trackable events must be added to the catalogue in `modules/analytics.py`

---

## Questions?

Open a [Discussion](https://github.com/chemenova/intelliform/discussions) or email shehan@chemenova.com.

**Thank you for making green chemistry more accessible. 🌿**
