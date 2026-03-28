from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from api.models import (
    FormulateRequest, ParetoRequest, BayesianRequest,
    QSARRequest, ReformulateRequest, HealthResponse
)
from api.memory import memory
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

app = FastAPI(title="IntelliForm API", version="2.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "https://chemenova.com"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/health", response_model=HealthResponse)
def health():
    return {
        "status": "ok",
        "version": "2.1.0",
        "modules": [
            "optimizer", "pareto", "bayesian", "eco",
            "qsar", "regulatory", "stability", "carbon",
            "certifications", "agents", "nlp"
        ]
    }

@app.get("/api/v1/verticals")
def get_verticals():
    return ["personal_care", "home_care", "industrial", "pharma", "food", "agriculture"]

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
        import pandas as pd
        db = pd.read_csv("data/ingredients.csv")
        return controller.run(
            req.input_text, req.vertical,
            req.batch_size, req.opt_mode,
            req.constraints, db
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/optimize/pareto")
async def optimize_pareto(req: ParetoRequest):
    try:
        from modules.pareto import run_pareto_optimization
        import pandas as pd
        db = pd.read_csv("data/ingredients.csv")
        result = run_pareto_optimization(db, req.constraints, req.n_solutions)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/optimize/bayesian")
async def optimize_bayesian(req: BayesianRequest):
    try:
        from modules.bayesian import run_bayesian_optimization
        import pandas as pd
        db = pd.read_csv("data/ingredients.csv")
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
        from modules.reformulation import run_reformulation
        import pandas as pd
        db = pd.read_csv("data/ingredients.csv")
        result = run_reformulation(req.blend, req.failure_type, req.vertical, db, req.constraints)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
