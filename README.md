# IntelliForm™ v0.9

**Agentic AI Green Chemistry Formulation Platform**  
ChemeNova LLC × ChemRich Global

> Natural language → RDKit cheminformatics → Multi-objective Pareto optimization → EcoMetrics™ sustainability scoring → Regulatory intelligence → ChemRich NJ pilot batch

[![MIT License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue)](https://python.org)
[![Streamlit](https://img.shields.io/badge/streamlit-1.32+-red)](https://streamlit.io)

---

## What IntelliForm Does

IntelliForm is an open-source agentic AI platform for designing sustainable chemical formulations. It combines:

- **AI-powered NL parsing** — describe your formulation need in plain English
- **PuLP LP optimizer** — finds the optimal blend from 35+ bio-based ingredients
- **NSGA-III Pareto optimizer** — multi-objective optimization across cost, sustainability, and performance
- **EcoMetrics™ scoring** — 5-axis sustainability radar vs. petrochemical baseline
- **QSAR/QSPR models** — XGBoost-based property prediction from molecular fingerprints
- **Regulatory intelligence** — REACH, EPA Safer Choice, EU Ecolabel, COSMOS status per ingredient
- **Branded PDF proposals** — one-click professional deliverable
- **Supabase persistence** — project history and active learning feedback loop
- **PostHog analytics** — usage tracking and funnel analysis

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/chemenova/intelliform.git
cd intelliform

# 2. Create environment (RDKit requires conda)
conda create -n intelliform python=3.11 -y
conda activate intelliform
conda install -c conda-forge rdkit -y

# 3. Install dependencies
pip install -r requirements.txt

# 4. Configure (all optional for first run)
cp .env.example .env
# Add GROQ_API_KEY for full NL understanding (free at console.groq.com)

# 5. Run
streamlit run app.py
```

---

## Repo Structure

```
intelliform/
├── app.py                          ← Main UI (7 tabs)
├── requirements.txt
├── .env.example
├── .streamlit/
│   └── config.toml                 ← Dark theme
├── data/
│   └── ingredients_db.csv          ← 35 bio-based ingredients
├── modules/
│   ├── llm_parser.py               ← Groq / Ollama / regex NL parser
│   ├── optimizer.py                ← PuLP LP solver
│   ├── pareto_optimizer.py         ← NSGA-III multi-objective optimizer
│   ├── ecometrics.py               ← EcoMetrics™ 5-axis scoring
│   ├── qsar.py                     ← QSAR/QSPR ML models (XGBoost)
│   ├── regulatory.py               ← Regulatory intelligence engine
│   ├── persistence.py              ← Supabase / in-memory storage
│   ├── pdf_proposal.py             ← Branded PDF generation (ReportLab)
│   ├── agents.py                   ← LLM agent swarm
│   ├── chem_utils.py               ← RDKit helpers
│   └── analytics.py                ← PostHog wrapper
├── migrations/
│   └── 001_create_tables.sql       ← Run in Supabase SQL editor
└── tests/
    └── test_optimizer.py
```

---

## Tabs

| Tab | Description |
|-----|-------------|
| 🚀 Agentic Swarm | NL parse → optimize → EcoMetrics → Regulatory → QSAR |
| 🌿 EcoMetrics™ | Sustainability radar + petrochemical baseline |
| 📋 Regulatory | REACH, EPA Safer Choice, COSMOS, EU Ecolabel per ingredient |
| 📈 Pareto Frontier | 3D multi-objective scatter + TOPSIS recommendation |
| 🔬 Model Card | QSAR benchmarks, live prediction, active learning |
| 📊 ROI & History | Project history + Supabase migration SQL |
| 📄 Proposal | Branded PDF or Markdown download |

---

## Supabase Setup (optional)

1. Create a project at [supabase.com](https://supabase.com)
2. Run `migrations/001_create_tables.sql` in the SQL editor
3. Add `SUPABASE_URL` and `SUPABASE_ANON_KEY` to `.env`

Without Supabase, IntelliForm falls back to in-memory storage (session-only).

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Chemistry contributions (new ingredients, SMILES corrections) are as welcome as code contributions.

---

## Citation

If you use IntelliForm in research:

```
Makani, S. et al. "IntelliForm: An Agentic AI Platform for Green Chemistry Formulation."
ChemRxiv (2026). DOI: 10.26434/chemrxiv.15000857
```

---

## License

MIT — see [LICENSE](LICENSE)

**ChemeNova LLC × ChemRich Global | shehan@chemenova.com | Pearl River, NJ**
