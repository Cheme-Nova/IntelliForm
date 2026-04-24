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

@app.post("/api/v1/formulate")
async def formulate(req: FormulateRequest):
    try:
        from api.controller import controller
        return controller.run(
            req.input_text, req.vertical,
            req.batch_size, req.opt_mode,
            req.constraints
        )
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
