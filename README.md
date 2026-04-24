# IntelliForm: AI-Powered Chemical Formulation & Process Optimization

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://intelliform.streamlit.app/)

**IntelliForm** is an open-source agentic AI framework designed to accelerate specialty chemical R&D. It translates natural language requirements into high-precision, RDKit-validated chemical formulations and manufacturing protocols.

---

## 🌐 Ecosystem & Live Links

* **[Full Platform (Streamlit Lab)](https://intelliform.streamlit.app/):** The core R&D cockpit for internal research, lab-scale prototyping, and deep optimization.
* **[IntelliForm Free (Web)](https://chemenova.com/intelliform/free):** Public production app featuring Google Sign-in and rate-limited QSAR generation.
* **[Interactive Demo](https://chemenova.com/intelliform/demo):** A guided showcase of agentic reasoning in chemical manufacturing.
* **[Official Project Site](https://chemenova.com/intelliform/):** Documentation and ecosystem overview by ChemeNova.

---

## 🚀 Key Capabilities

* **Natural Language → Chemistry:** Convert briefs into structured SMILES and formulation tables using LLM-based parsing.
* **Multi-Objective Optimization:** Balance assay purity (e.g., target 98–99% $CaCl_2·2H_2O$), cost, and sustainability metrics.
* **Regulatory Intelligence:** Integrated QSAR modeling and regulatory screening (v0.9+).
* **Agentic Reliability:** Local JSON normalization and contradiction filtering to ensure "pilot-ready" outputs.

---

## 🛠 Architecture: "Lab vs. Product"

This repository supports a dual-stack architecture to separate research from public access:

| Feature | **Streamlit Lab** (`app.py`) | **Public Web/API** (`web/` + `api/`) |
| :--- | :--- | :--- |
| **Primary Use** | Internal R&D / Heavy Computation | Public Showcase / User Onboarding |
| **Tech Stack** | Python + Streamlit + RDKit | React (Vite) + FastAPI + Supabase |
| **UX Focus** | Complexity & Scientific Tooling | Performance, Auth, & Responsiveness |
| **Auth** | Local/Session | Google Sign-In |

---

## 📖 Quick Start

### 1. Streamlit Lab (The Chemistry Workbench)
```bash
conda create -n intelliform python=3.11 -y
conda activate intelliform
conda install -c conda-forge rdkit -y
pip install -r requirements.txt
streamlit run app.py
