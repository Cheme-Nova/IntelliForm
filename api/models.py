from pydantic import BaseModel
from typing import Optional, List, Dict, Any

class FormulateRequest(BaseModel):
    input_text: str
    vertical: str
    batch_size: float = 1000.0
    opt_mode: str = "auto"
    constraints: Optional[Dict[str, Any]] = {}

class ParetoRequest(BaseModel):
    vertical: str
    constraints: Optional[Dict[str, Any]] = {}
    n_solutions: int = 10

class BayesianRequest(BaseModel):
    vertical: str
    constraints: Optional[Dict[str, Any]] = {}
    n_iterations: int = 5
    state: Optional[Dict[str, Any]] = None

class QSARRequest(BaseModel):
    smiles: List[str]
    properties: Optional[List[str]] = ["biodegradability", "ecotox", "performance"]

class ReformulateRequest(BaseModel):
    blend: Dict[str, float]
    failure_type: str
    vertical: str
    constraints: Optional[Dict[str, Any]] = {}

class HealthResponse(BaseModel):
    status: str
    version: str
    modules: List[str]
