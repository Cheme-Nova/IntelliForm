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

class MOBORequest(BaseModel):
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


# ── Admin API key endpoints ───────────────────────────────────────────────────

# ── Supply chain / supplier portal ───────────────────────────────────────────

class SupplierRegisterRequest(BaseModel):
    name: str = Field(..., description="Company or individual supplier name")
    email: str = Field(..., description="Primary contact email")
    country: str = Field(..., description="Country of operation (ISO alpha-2 or full name)")
    certifications: Optional[List[str]] = Field([], description="e.g. ISO9001, Ecocert, REACH")
    description: str = Field("", description="Short company bio / product focus")

class SupplierResponse(BaseModel):
    id: str
    name: str
    email: str
    country: str
    certifications: List[str]
    status: str
    api_key_prefix: str
    joined_at: str
    description: str

class SupplierListingRequest(BaseModel):
    ingredient_name: str
    price_per_kg: float = Field(..., gt=0)
    currency: str = "USD"
    moq_kg: float = Field(25.0, gt=0, description="Minimum order quantity in kg")
    lead_time_days: int = Field(14, ge=1, description="Typical lead time in days")
    availability: str = Field("in_stock", description="in_stock | limited | out_of_stock")
    certifications: Optional[List[str]] = []
    valid_until: Optional[str] = Field(None, description="ISO date until which price is valid")

class SupplierListingResponse(BaseModel):
    id: str
    supplier_id: str
    supplier_name: str
    ingredient_name: str
    price_per_kg: float
    currency: str
    moq_kg: float
    lead_time_days: int
    availability: str
    certifications: List[str]
    valid_until: Optional[str]
    created_at: str
    updated_at: str

class SupplyChainBlendRequest(BaseModel):
    blend: Dict[str, Any] = Field(..., description="ingredient → weight %")


# ── Admin API key endpoints ───────────────────────────────────────────────────

class APIKeyCreateRequest(BaseModel):
    org_name: str = Field(..., description="Organisation or customer name")
    scopes: Optional[List[str]] = Field(None, description="Subset of ALL_SCOPES; defaults to all")
    rate_limit_per_hour: int = Field(200, description="Max requests per sliding 1-hour window")
    description: str = Field("", description="Human-readable label for this key")

class APIKeyResponse(BaseModel):
    id: str
    key_prefix: str
    org_name: str
    scopes: List[str]
    rate_limit_per_hour: int
    created_at: str
    last_used_at: Optional[str]
    active: bool
    description: str
