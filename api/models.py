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

class DigitalTwinRequest(BaseModel):
    blend: Dict[str, float]
    batch_kg: float = 500.0
    manufacturing_process: Optional[str] = None
    manufacturing_route: Optional[str] = None
    vertical: str = "all"
    grid_region: str = "US"
    compare_scales: bool = False

class SDLProtocolRequest(BaseModel):
    blend: Dict[str, float]
    batch_scale_g: float = 10.0
    platform: str = "opentrons_ot2"
    target_specs: Optional[Dict[str, float]] = {}

class SDLIngestRequest(BaseModel):
    blend: Dict[str, float]
    measured: Dict[str, float]
    target_specs: Dict[str, float]
    vertical: str = "all"
    iteration: int = 1
    batch_id: Optional[str] = None

class HealthResponse(BaseModel):
    status: str
    version: str
    modules: List[str]
