from pydantic import BaseModel, Field
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


# ── Active Learning / LIMS endpoints ─────────────────────────────────────────

class ALFeedbackRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string of the tested compound")
    target: str = Field(..., description="'Biodegradability' | 'Ecotoxicity' | 'Performance'")
    actual_value: float = Field(..., description="Lab-measured value (Bio: 0-100%, Etox: 1-10, Perf: 0-100)")
    source_system: Optional[str] = Field(None, description="LIMS system identifier e.g. 'LabVantage', 'SAP-QM'")
    batch_id: Optional[str] = Field(None, description="Lab batch or experiment reference ID")

class ALFeedbackRecord(BaseModel):
    smiles: str
    target: str
    actual_value: float
    batch_id: Optional[str] = None

class ALFeedbackBatchRequest(BaseModel):
    records: List[ALFeedbackRecord] = Field(..., description="List of lab results to submit (max 200)")
    source_system: Optional[str] = Field(None, description="LIMS system identifier")


# ── Webhook endpoints ─────────────────────────────────────────────────────────

class WebhookRegisterRequest(BaseModel):
    url: str = Field(..., description="HTTPS endpoint that will receive POST events")
    events: Optional[List[str]] = Field(None, description="Event names to subscribe to (default: all)")
    description: Optional[str] = Field("", description="Human-readable label for this webhook")

class WebhookResponse(BaseModel):
    id: str
    url: str
    events: List[str]
    created_at: str
    active: bool
    description: str
