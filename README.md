# IntelliForm

Open-source AI formulation intelligence for specialty chemical R&D teams.

This repo now supports two modes:

- `web/` + `api/`: the public free product shape
- `app.py`: the original Streamlit lab app for demos, research, and internal workflows

## What To Deploy

If you want a public version that people can actually use, deploy:

- `web/` to Vercel
- `api/` to Render or Railway

Do not treat Streamlit as the long-term public product surface. Keep Streamlit for:

- internal chemistry workbench usage
- demos and investor screenshots
- rapid module prototyping

## Why This Split

The React + FastAPI stack is a better public free product because it gives you:

- a cleaner product UX than a monolithic Streamlit app
- easier auth and rate limiting
- better control over API costs
- a path to separate public and enterprise editions

## Public Free Tier Features

The API now includes a lightweight free-tier protection layer:

- hourly request limits for public generation endpoints
- separate QSAR request budget
- optional public API key support
- CORS configured for local dev plus deploy-time custom origins

Environment variables:

- `INTELLIFORM_FREE_TIER=1`
- `INTELLIFORM_FREE_TIER_MAX_REQUESTS_PER_HOUR=8`
- `INTELLIFORM_FREE_TIER_MAX_QSAR_PER_HOUR=20`
- `INTELLIFORM_FREE_TIER_REQUIRE_API_KEY=0`
- `INTELLIFORM_PUBLIC_API_KEY=...` (optional)
- `ALLOWED_ORIGINS=https://your-vercel-app.vercel.app,https://yourdomain.com`

Frontend environment:

- `VITE_API_URL=https://your-intelliform-api.onrender.com`
- `VITE_PUBLIC_MODE=1`
- `VITE_PUBLIC_API_KEY=...` (optional, if you turn on public API key mode)

## Quick Start

### 1. Streamlit Lab App

```bash
conda create -n intelliform python=3.11 -y
conda activate intelliform
conda install -c conda-forge rdkit -y
pip install -r requirements.txt
streamlit run app.py
```

### 2. API

```bash
pip install -r requirements.txt
uvicorn api.main:app --reload
```

API runs at `http://localhost:8000`.

### 3. Web

```bash
cd web
npm install
cp .env.example .env.local
npm run dev
```

Web runs at `http://localhost:5173`.

## Deployment

### Backend on Render

This repo includes:

- [render.yaml](/Users/makani/IntelliForm/render.yaml)
- [Procfile](/Users/makani/IntelliForm/Procfile)
- [runtime.txt](/Users/makani/IntelliForm/runtime.txt)

Recommended first deployment:

1. Create a new Render Web Service from this repo.
2. Use the existing `render.yaml`.
3. Add your env vars:
   - `GROQ_API_KEY`
   - `ALLOWED_ORIGINS`
   - optional Supabase/PostHog keys
4. Deploy and verify `/health`.

### Frontend on Vercel

This repo includes:

- [web/vercel.json](/Users/makani/IntelliForm/web/vercel.json)
- [web/.env.example](/Users/makani/IntelliForm/web/.env.example)

Recommended first deployment:

1. Import the `web/` directory as a Vercel project.
2. Set `VITE_API_URL` to your Render API URL.
3. Deploy.

## Product Recommendation

Best structure going forward:

- `IntelliForm Public`: free, rate-limited, showcase-focused, lighter feature set
- `IntelliForm Streamlit Lab`: internal R&D cockpit
- `IntelliForm Enterprise`: paid product with auth, persistence, workflows, and higher-trust outputs

That lets you grow without forcing one codebase surface to do everything.

## Reliability Upgrades

This repo already includes the public-controller hardening work:

- local JSON normalization instead of brittle provider-only enforcement
- canonical vertical mapping across parser, controller, and UI
- brief-aware contradiction filtering
- proof-style showcase examples under [examples/showcase/README.md](/Users/makani/IntelliForm/examples/showcase/README.md)

## License

MIT — see [LICENSE](/Users/makani/IntelliForm/LICENSE)
