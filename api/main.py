# redeploy: 2026-05-25 — add pharma/carbon-passport/pubchem endpoints
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import asyncio
import json
import queue
import threading
from typing import Optional, Dict, Any
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel
from api.models import (
    FormulateRequest, ParetoRequest, BayesianRequest,
    QSARRequest, ReformulateRequest, HealthResponse,
    ALFeedbackRequest, ALFeedbackBatchRequest,
)

class RefineRequest(BaseModel):
    instruction: str
    current_result: Dict[str, Any]
    vertical: str
    batch_size: float = 1000.0
    opt_mode: str = "auto"
from api.memory import memory
from api.public_access import (
    FREE_TIER_ENABLED,
    get_client_id,
    validate_public_access,
)
from api.auth import REQUIRE_SIGNIN, extract_bearer_token, is_auth_enabled, verify_supabase_user
from modules.persistence import load_projects_for_user, load_recent_usage_count, record_usage, save_project
from api.controller import _canonicalize_vertical, _merge_constraints, _serialize
from modules.verticals import filter_db_by_vertical

DB_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "ingredients_db.csv")
_EXTRA_ORIGINS = [origin.strip() for origin in os.getenv("ALLOWED_ORIGINS", "").split(",") if origin.strip()]

app = FastAPI(title="IntelliForm API", version="2.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://localhost:5174",
        "https://chemenova.com",
        "https://www.chemenova.com",
        *_EXTRA_ORIGINS,
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.middleware("http")
async def public_quota_middleware(request: Request, call_next):
    protected_routes = {
        "/api/v1/formulate": "formulate",
        "/api/v1/optimize/pareto": "optimize",
        "/api/v1/optimize/bayesian": "optimize",
        "/api/v1/predict/qsar": "qsar",
        "/api/v1/reformulate": "reformulate",
    }
    bucket = protected_routes.get(request.url.path)
    if bucket and request.method.upper() == "POST":
        allowed, retry_after, error = validate_public_access(request, bucket=bucket)
        if not allowed:
            payload = {
                "detail": error,
                "retry_after_seconds": retry_after,
                "free_tier_enabled": FREE_TIER_ENABLED,
            }
            headers = {"Retry-After": str(retry_after)} if retry_after else {}
            return JSONResponse(status_code=429, content=payload, headers=headers)
    response = await call_next(request)
    response.headers["X-IntelliForm-Free-Tier"] = "true" if FREE_TIER_ENABLED else "false"
    response.headers["X-Client-Id"] = get_client_id(request)
    return response

@app.get("/")
def root():
    return {
        "name": "IntelliForm API",
        "version": "2.1.0",
        "free_tier_enabled": FREE_TIER_ENABLED,
        "docs": "/docs",
        "health": "/health",
    }


def _require_user(request: Request):
    if not REQUIRE_SIGNIN or not is_auth_enabled():
        return None
    token = extract_bearer_token(request)
    user = verify_supabase_user(token)
    if not user:
        raise HTTPException(status_code=401, detail="Sign in required to generate formulations.")
    return user

@app.get("/health", response_model=HealthResponse)
def health():
    return {
        "status": "ok",
        "version": "2.1.0",
        "modules": [
            "optimizer", "pareto_optimizer", "bayesian_optimizer", "ecometrics",
            "qsar", "regulatory", "vertical_regulatory", "stability",
            "carbon_credits", "certification_oracle", "agents", "llm_parser"
        ]
    }

@app.get("/api/v1/verticals")
def get_verticals():
    return [
        "personal_care",
        "industrial",
        "agricultural",
        "pharmaceutical",
        "food",
        "fabric_laundry",
        "paint_coatings",
    ]

@app.get("/api/v1/failure-types")
def get_failure_types():
    return ["viscosity", "stability", "pH", "color", "odor", "certification", "eco_score"]

@app.get("/api/v1/memory")
def get_memory(n: int = 10):
    return memory.recent(n)


@app.get("/api/v1/me")
def get_me(request: Request):
    user = _require_user(request)
    if not user:
        return {"auth_enabled": is_auth_enabled(), "signed_in": False}
    usage_count = load_recent_usage_count(user.get("id", ""))
    return {
        "auth_enabled": is_auth_enabled(),
        "signed_in": True,
        "user": {
            "id": user.get("id"),
            "email": user.get("email"),
            "name": (user.get("user_metadata") or {}).get("full_name") or user.get("email"),
            "avatar_url": (user.get("user_metadata") or {}).get("avatar_url"),
        },
        "usage": {
            "last_24h_formulations": usage_count,
        },
    }


@app.get("/api/v1/projects")
def get_projects(request: Request, limit: int = 25):
    user = _require_user(request)
    if not user:
        return []
    return load_projects_for_user(user.get("email", ""), limit=limit)

@app.post("/api/v1/formulate")
async def formulate(req: FormulateRequest, request: Request):
    user = _require_user(request)
    try:
        from api.controller import controller
        result = controller.run(
            req.input_text, req.vertical,
            req.batch_size, req.opt_mode,
            req.constraints
        )
        if user and result.get("result", {}).get("success"):
            record_usage(user.get("id", ""), user.get("email", ""), "formulate")
            save_project({
                "application": req.vertical,
                "blend": result.get("result", {}).get("blend", {}),
                "cost": result.get("result", {}).get("cost_per_kg"),
                "bio": result.get("result", {}).get("bio_pct"),
                "perf": result.get("result", {}).get("perf_score"),
                "eco_score": (result.get("eco") or {}).get("eco_score") if isinstance(result.get("eco"), dict) else None,
                "eco_grade": (result.get("eco") or {}).get("grade") if isinstance(result.get("eco"), dict) else None,
                "optimizer": req.opt_mode,
                "parser": (result.get("parsed") or {}).get("parser_backend"),
                "relaxed": (result.get("result") or {}).get("relaxed", False),
                "input": req.input_text,
            }, session_id=user.get("id", ""), user_email=user.get("email", ""))
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/optimize/pareto")
async def optimize_pareto(req: ParetoRequest):
    try:
        from modules.pareto_optimizer import run_pareto_optimization
        import pandas as pd
        db = pd.read_csv(DB_PATH)
        vertical = _canonicalize_vertical(req.vertical)
        filtered_db = filter_db_by_vertical(db, vertical)
        max_cost, min_bio, min_perf = _merge_constraints(
            type("Parsed", (), {"max_cost": 999.0, "min_bio": 0.0, "min_perf": 0.0})(),
            req.constraints,
            vertical,
        )
        result = run_pareto_optimization(filtered_db, max_cost, min_bio, min_perf, req.n_solutions)
        return _serialize(result)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/optimize/bayesian")
async def optimize_bayesian(req: BayesianRequest):
    try:
        from modules.bayesian_optimizer import run_bayesian_optimization
        import pandas as pd
        db = pd.read_csv(DB_PATH)
        vertical = _canonicalize_vertical(req.vertical)
        filtered_db = filter_db_by_vertical(db, vertical)
        max_cost, min_bio, min_perf = _merge_constraints(
            type("Parsed", (), {"max_cost": 999.0, "min_bio": 0.0, "min_perf": 0.0})(),
            req.constraints,
            vertical,
        )
        result, state = run_bayesian_optimization(
            filtered_db,
            max_cost,
            min_bio,
            min_perf,
            state=req.state,
            n_random_init=max(5, req.n_iterations),
            vertical=vertical,
        )
        return {
            "result": _serialize(result),
            "state": _serialize(state),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# ── QSAR lazy-init state ─────────────────────────────────────────────────────
_qsar_initialized = False
_qsar_db = None

def _ensure_qsar_initialized():
    global _qsar_initialized, _qsar_db
    if not _qsar_initialized:
        import pandas as pd
        from modules.qsar import initialize_models
        _qsar_db = pd.read_csv(DB_PATH)
        initialize_models(_qsar_db)
        _qsar_initialized = True
    return _qsar_db


def _check_lims_key(request: Request):
    """Verify X-Api-Key header against LIMS_API_KEY env var. Open if env var unset."""
    required = os.getenv("LIMS_API_KEY", "")
    if not required:
        return
    provided = request.headers.get("X-Api-Key", "")
    if provided != required:
        raise HTTPException(status_code=401, detail="Invalid or missing X-Api-Key header.")


@app.post("/api/v1/predict/qsar")
async def predict_qsar(req: QSARRequest):
    try:
        from modules.qsar import predict_properties
        _ensure_qsar_initialized()
        predictions = [predict_properties(smiles) for smiles in req.smiles]
        return {
            "predictions": _serialize(predictions),
            "properties": req.properties,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── Active Learning / LIMS endpoints ─────────────────────────────────────────

@app.get("/api/v1/al/stats")
async def al_stats_endpoint(request: Request):
    """Return active learning model card stats (feedback count, rounds, conformal status)."""
    _check_lims_key(request)
    try:
        from modules.qsar import al_stats
        _ensure_qsar_initialized()
        return al_stats()
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/v1/al/feedback")
async def al_feedback_list(request: Request, limit: int = 100):
    """Return the most recent lab feedback records (newest first)."""
    _check_lims_key(request)
    try:
        from modules.qsar import get_feedback_records
        _ensure_qsar_initialized()
        records = get_feedback_records()
        trimmed = records[-limit:][::-1]
        return {"records": [r.__dict__ for r in trimmed], "total": len(records)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/v1/al/experiment-queue")
async def al_experiment_queue(request: Request, top_k: int = 10, vertical: str = "personal_care"):
    """Return top-k high-uncertainty ingredient candidates for the next lab round."""
    _check_lims_key(request)
    try:
        import pandas as pd
        from modules.qsar import query_uncertain_candidates
        from modules.verticals import filter_db_by_vertical
        db = pd.read_csv(DB_PATH)
        filtered = filter_db_by_vertical(db, _canonicalize_vertical(vertical))
        _ensure_qsar_initialized()
        queue_df = query_uncertain_candidates(filtered, top_k=top_k)
        if queue_df is None or queue_df.empty:
            return {"candidates": [], "note": "ensemble not trained — submit feedback first"}
        return {
            "candidates": queue_df.to_dict(orient="records"),
            "vertical": vertical,
            "top_k": top_k,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/al/feedback")
async def al_submit_feedback(req: ALFeedbackRequest, request: Request):
    """Submit a single lab result. Triggers immediate model retrain."""
    _check_lims_key(request)
    try:
        import pandas as pd
        from modules.qsar import submit_feedback
        db = _ensure_qsar_initialized()
        result = submit_feedback(
            smiles=req.smiles,
            target=req.target,
            actual_value=req.actual_value,
            db=db,
        )
        return {
            "status": "ok",
            "al_round": result.get("al_round"),
            "total_feedback": result.get("total_feedback"),
            "retrained": result.get("retrained", False),
            "source_system": req.source_system,
            "batch_id": req.batch_id,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/al/feedback/batch")
async def al_submit_feedback_batch(req: ALFeedbackBatchRequest, request: Request):
    """Submit up to 200 lab results in one call. Single retrain after all records are stored."""
    _check_lims_key(request)
    try:
        from modules.qsar import submit_feedback_batch
        db = _ensure_qsar_initialized()
        records = [r.model_dump() for r in req.records]
        result = submit_feedback_batch(records, db)
        return {
            **result,
            "source_system": req.source_system,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/stream/formulate")
async def stream_formulate(req: FormulateRequest, request: Request):
    """SSE streaming endpoint — yields progress events then the complete result."""
    user = _require_user(request)
    event_queue: queue.Queue = queue.Queue()

    def progress_cb(event):
        event_queue.put(event)

    def run_in_thread():
        try:
            from api.controller import controller
            controller.run(
                req.input_text, req.vertical,
                req.batch_size, req.opt_mode,
                req.constraints,
                progress_cb=progress_cb,
            )
        except Exception as e:
            event_queue.put({"step": "error", "message": str(e)})
        finally:
            event_queue.put(None)  # sentinel

    thread = threading.Thread(target=run_in_thread, daemon=True)
    thread.start()

    async def event_stream():
        loop = asyncio.get_event_loop()
        while True:
            event = await loop.run_in_executor(None, event_queue.get)
            if event is None:
                break
            yield f"data: {json.dumps(event)}\n\n"
            if event.get("step") in ("complete", "error"):
                break

    return StreamingResponse(
        event_stream(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
            "Connection": "keep-alive",
        },
    )


@app.post("/api/v1/refine")
async def refine_formulation(req: RefineRequest, request: Request):
    """Refine an existing formulation with a natural-language instruction."""
    _require_user(request)
    try:
        from api.controller import controller
        result = controller.refine(
            req.current_result,
            req.instruction,
            req.vertical,
            req.batch_size,
            req.opt_mode,
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


class PubChemRequest(BaseModel):
    names: list[str]

class PharmaRequest(BaseModel):
    blend: Dict[str, Any]
    bcs_class: str = "I"
    dosage_form: str = "immediate_release_tablet"
    target_markets: list[str] = ["USA", "EU"]
    is_generic: bool = False
    is_pediatric: bool = False

class CarbonPassportRequest(BaseModel):
    blend: Dict[str, Any]
    batch_kg: float = 500.0
    product_name: str = "IntelliForm Formulation"

@app.post("/api/v1/pharma/deep-dive")
async def pharma_deep_dive(req: PharmaRequest):
    """Run ICH/BCS pharmaceutical deep-dive analysis on a blend."""
    try:
        import pandas as pd
        from modules.pharma import run_pharma_deep_dive
        db = pd.read_csv(DB_PATH)
        result = run_pharma_deep_dive(
            req.blend, db,
            bcs_class=req.bcs_class,
            dosage_form=req.dosage_form,
            target_markets=req.target_markets,
            is_generic=req.is_generic,
            is_pediatric=req.is_pediatric,
        )
        return _serialize(result)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/export/carbon-passport")
async def export_carbon_passport(req: CarbonPassportRequest):
    """Generate an ISO 14067 Carbon Passport JSON for a formulation."""
    try:
        import pandas as pd
        from modules.carbon_passport import generate_carbon_passport, passport_to_json
        db = pd.read_csv(DB_PATH)
        passport = generate_carbon_passport(
            req.blend, db,
            product_name=req.product_name,
            batch_kg=req.batch_kg,
        )
        return {"json": passport_to_json(passport)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/pubchem/enrich")
async def pubchem_enrich(req: PubChemRequest):
    """Enrich a list of ingredient names with PubChem identity data (CAS, SMILES, MW)."""
    try:
        from modules.pubchem import enrich_blend
        if not req.names:
            return {}
        blend = {n: 0 for n in req.names[:30]}  # cap at 30 names
        return enrich_blend(blend)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/reformulate")
async def reformulate(req: ReformulateRequest):
    try:
        import pandas as pd
        from modules.reformulation_intelligence import run_reformulation_intelligence
        db = pd.read_csv(DB_PATH)
        failure_map = {
            "viscosity": "viscosity_too_high",
            "stability": "stability_failure",
            "pH": "ph_too_high",
            "color": "colour_complaint",
            "odor": "odour_complaint",
            "certification": "certification_rejection",
            "eco_score": "performance_shortfall",
        }
        normalized_failure = failure_map.get(req.failure_type, req.failure_type)
        default_test_data = {
            "viscosity_too_high": {"measured_viscosity_cP": 8500, "target_viscosity_max_cP": 5000},
            "stability_failure": {"stability_condition": "40C / 75% RH", "failure_timepoint_weeks": 4, "failure_observation": "phase drift"},
            "ph_too_high": {"measured_ph": 8.8, "target_ph_min": 5.0, "target_ph_max": 6.5},
            "colour_complaint": {"observed_colour": "yellow shift", "target_colour": "clear"},
            "odour_complaint": {"odour_description": "sharp solvent note"},
            "certification_rejection": {"certification_name": "Public demo certification screen", "rejected_ingredient": next(iter(req.blend.keys()), "unknown"), "rejection_reason": "Requires cleaner input profile"},
            "performance_shortfall": {"measured_performance": 62, "target_performance": 80, "performance_metric": "eco score proxy"},
        }
        result = run_reformulation_intelligence(
            req.blend,
            db,
            normalized_failure,
            default_test_data.get(normalized_failure, {}),
        )
        return _serialize(result)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
