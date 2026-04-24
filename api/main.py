import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from api.models import (
    FormulateRequest, ParetoRequest, BayesianRequest,
    QSARRequest, ReformulateRequest, HealthResponse
)
from api.memory import memory
from api.public_access import (
    FREE_TIER_ENABLED,
    get_client_id,
    validate_public_access,
)
from api.auth import REQUIRE_SIGNIN, extract_bearer_token, is_auth_enabled, verify_supabase_user
from modules.persistence import load_projects_for_user, load_recent_usage_count, record_usage, save_project

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
    try:
        user = _require_user(request)
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
        import pandas as pd
        from modules.pareto_optimizer import run_pareto_optimization
        db = pd.read_csv(DB_PATH)
        result = run_pareto_optimization(db, req.constraints, req.n_solutions)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/optimize/bayesian")
async def optimize_bayesian(req: BayesianRequest):
    try:
        import pandas as pd
        from modules.bayesian_optimizer import run_bayesian_optimization
        db = pd.read_csv(DB_PATH)
        result = run_bayesian_optimization(db, req.constraints, req.n_iterations)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/predict/qsar")
async def predict_qsar(req: QSARRequest):
    try:
        from modules.qsar import predict_properties
        result = predict_properties(req.smiles, req.properties)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/reformulate")
async def reformulate(req: ReformulateRequest):
    try:
        import pandas as pd
        from modules.reformulation_intelligence import run_reformulation
        db = pd.read_csv(DB_PATH)
        result = run_reformulation(req.blend, req.failure_type, req.vertical, db, req.constraints)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
