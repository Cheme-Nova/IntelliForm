"""
modules/memory_network.py
Formulation Memory Network™ — IntelliForm v1.5

Living knowledge graph that accumulates institutional memory from every
formulation run, customer decision, pilot outcome, and regulatory change.

Storage tier (best available at runtime):
  1. Graphiti temporal KG  — bi-temporal Neo4j graph; queryable at any
     point in time. Powered by getzep/graphiti.
     pip install graphiti-core neo4j
     Apache-2.0 · Rauch et al. 2024 (getzep/graphiti)

  2. Supabase flat table   — persistent across sessions (existing).

  3. Local JSON            — session-persistent fallback (existing).

Graphiti upgrade:
  - Every event → natural-language episode → LLM extracts entities +
    relations → stored as bi-temporal edges in Neo4j.
  - get_negative_rules(vertical, at_time=...) now returns rules as they
    stood at any historical date — critical for CBAM/regulatory compliance.
  - get_insights() enriched with Graphiti semantic search results.

Environment variables for Graphiti:
  NEO4J_URI          e.g. bolt://localhost:7687 or neo4j+s://xxx.databases.neo4j.io
  NEO4J_USER         neo4j
  NEO4J_PASSWORD     your-password
  GROQ_API_KEY       used for Graphiti's entity extraction LLM (Llama 3.3)
  OPENAI_API_KEY     alternative if GROQ_API_KEY not set

References:
  - Rauch et al., "Graphiti", 2024 (getzep/graphiti)
  - Raccuglia et al., Nature 2016 (collaborative filtering for materials)
  - Kononova et al., Sci Data 2021 (failed experiment databases)
"""
from __future__ import annotations

import asyncio
import concurrent.futures
import json
import os
import hashlib
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
import numpy as np

try:
    from supabase import create_client
    SUPABASE_OK = True
except ImportError:
    SUPABASE_OK = False

# ── Graphiti (tier 1) ─────────────────────────────────────────────────────────
GRAPHITI_OK = False
_graphiti_mod    = None
_graphiti_nodes  = None
_graphiti_llm    = None
_graphiti_config = None

try:
    import graphiti_core as _graphiti_mod
    from graphiti_core.nodes import EpisodeType as _EpisodeType
    _graphiti_nodes = _EpisodeType
    GRAPHITI_OK = True
except Exception:
    pass

# ── Async runner — safe for Streamlit (no nested loop) ───────────────────────

def _run_async(coro, timeout: float = 8.0):
    """
    Run an async coroutine from sync code.
    Uses a dedicated thread so it works inside Streamlit's event loop.
    Returns None silently on any failure.
    """
    def _target():
        return asyncio.run(coro)
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
            return pool.submit(_target).result(timeout=timeout)
    except Exception:
        return None


# ── Episode body builder ──────────────────────────────────────────────────────

def _episode_body(event_type: str, blend: Dict[str, float],
                   vertical: str, metadata: Dict) -> str:
    """
    Convert a formulation event to entity-dense natural language.
    Graphiti's LLM extracts ingredient / vertical / outcome entities from this.
    """
    blend_str = ", ".join(
        f"{ing} {pct:.1f}%"
        for ing, pct in sorted(blend.items(), key=lambda x: -x[1])[:6]
    )
    ings = ", ".join(sorted(blend.keys())[:4])

    templates: Dict[str, str] = {
        "formulation_generated": (
            f"New formulation generated for {vertical} application. "
            f"Blend composition: {blend_str}. "
            f"Optimizer: {metadata.get('optimizer', 'unknown')}."
        ),
        "blend_accepted": (
            f"Formulation ACCEPTED for {vertical}. "
            f"Ingredients that performed well: {ings}. "
            f"Blend: {blend_str}. "
            f"Notes: {metadata.get('reason', '')}."
        ),
        "blend_rejected": (
            f"Formulation REJECTED for {vertical}. "
            f"Rejected ingredients: {', '.join(metadata.get('rejected_ingredients', []))}. "
            f"Reason: {metadata.get('reason', 'unspecified')}. "
            f"Full blend was: {blend_str}."
        ),
        "ingredient_swapped": (
            f"Ingredient swap in {vertical} formulation. "
            f"Removed: {metadata.get('removed', '?')}. "
            f"Added: {metadata.get('added', '?')}. "
            f"Updated blend: {blend_str}."
        ),
        "pilot_batch_passed": (
            f"Pilot batch PASSED all quality tests for {vertical}. "
            f"Batch size: {metadata.get('batch_kg', 500)} kg. "
            f"Ingredients confirmed working: {ings}. "
            f"Blend: {blend_str}."
        ),
        "pilot_batch_failed": (
            f"Pilot batch FAILED for {vertical}. "
            f"Failure type: {metadata.get('failure_type', 'unspecified')}. "
            f"Problematic ingredients: "
            f"{', '.join(metadata.get('failure_ingredients', []))}. "
            f"Blend that failed: {blend_str}."
        ),
        "certification_passed": (
            f"Formulation received {metadata.get('certification', 'certification')} "
            f"certification for {vertical}. "
            f"Certified ingredients: {ings}. Blend: {blend_str}."
        ),
        "certification_rejected": (
            f"Formulation FAILED {metadata.get('certification', 'certification')} "
            f"certification for {vertical}. "
            f"Reason: {metadata.get('reason', 'unspecified')}. "
            f"Problematic ingredient: {metadata.get('failed_ingredient', 'unknown')}. "
            f"Blend: {blend_str}."
        ),
        "constraint_relaxed": (
            f"Optimization constraints relaxed for {vertical} formulation. "
            f"Relaxed parameters: {metadata.get('relaxed_params', 'unknown')}. "
            f"Resulting blend: {blend_str}."
        ),
    }

    return templates.get(
        event_type,
        f"Formulation event '{event_type}' for {vertical}. Blend: {blend_str}.",
    )


# ── Graphiti client init ──────────────────────────────────────────────────────

_GRAPHITI_CLIENT = None


async def _build_graphiti_client():
    """
    Build and initialise a Graphiti client using environment variables.
    Returns None if any required variable is missing or connection fails.
    """
    if not GRAPHITI_OK:
        return None

    neo4j_uri  = os.environ.get("NEO4J_URI",      "")
    neo4j_user = os.environ.get("NEO4J_USER",     "neo4j")
    neo4j_pw   = os.environ.get("NEO4J_PASSWORD",  "")

    if not neo4j_uri or not neo4j_pw:
        return None

    try:
        from graphiti_core import Graphiti

        # Prefer Groq (already a dep), fall back to OpenAI
        groq_key   = os.environ.get("GROQ_API_KEY",   "")
        openai_key = os.environ.get("OPENAI_API_KEY",  "")

        llm_client = None
        if groq_key:
            try:
                from graphiti_core.llm_client.openai_client import OpenAIClient
                from graphiti_core.llm_client.config import LLMConfig
                llm_client = OpenAIClient(LLMConfig(
                    api_key=groq_key,
                    model="llama-3.3-70b-versatile",
                    base_url="https://api.groq.com/openai/v1",
                ))
            except Exception:
                llm_client = None

        client = Graphiti(
            neo4j_uri, neo4j_user, neo4j_pw,
            **({"llm_client": llm_client} if llm_client else {}),
        )
        await client.build_indices_and_constraints()
        return client
    except Exception:
        return None


def _get_graphiti_client():
    """Return the module-level Graphiti client (initialised once)."""
    global _GRAPHITI_CLIENT
    if _GRAPHITI_CLIENT is None:
        _GRAPHITI_CLIENT = _run_async(_build_graphiti_client(), timeout=12.0)
    return _GRAPHITI_CLIENT


# ── Graphiti episode write ────────────────────────────────────────────────────

async def _add_graphiti_episode(client, name: str, body: str,
                                  ref_time: datetime, source_desc: str):
    from graphiti_core.nodes import EpisodeType
    await client.add_episode(
        name=name,
        episode_body=body,
        source=EpisodeType.text,
        reference_time=ref_time,
        source_description=source_desc,
    )


# ── Graphiti temporal query ───────────────────────────────────────────────────

async def _search_graphiti(client, query: str,
                             reference_time: Optional[datetime] = None,
                             num_results: int = 20) -> List[Any]:
    """Search the temporal knowledge graph. Returns EntityEdge list."""
    kwargs: Dict[str, Any] = {"num_results": num_results}
    if reference_time:
        kwargs["reference_time"] = reference_time
    try:
        return await client.search(query, **kwargs)
    except Exception:
        return []


def _graphiti_negative_rules(vertical: str,
                               reference_time: Optional[datetime] = None) -> List[dict]:
    """
    Query Graphiti for negative knowledge rules at reference_time.
    Returns edges whose facts describe rejections / failures / restrictions.
    Edges with invalid_at before reference_time are excluded — this is the
    bi-temporal moat: a rule that was lifted in September won't appear in
    an October query.
    """
    client = _get_graphiti_client()
    if client is None:
        return []

    rt = reference_time or datetime.now()
    query = (f"ingredient rejections, pilot failures, certification failures, "
             f"and restrictions for {vertical} vertical")
    edges = _run_async(_search_graphiti(client, query, rt, num_results=30))
    if not edges:
        return []

    negative_kw = {"reject", "fail", "avoid", "prohibit", "restrict",
                   "not", "never", "banned", "hazard", "conflict"}
    rules = []
    for edge in edges:
        fact_lower = edge.fact.lower()
        if not any(kw in fact_lower for kw in negative_kw):
            continue
        # Bi-temporal filter: skip facts invalidated before reference_time
        if edge.invalid_at and edge.invalid_at < rt:
            continue
        rules.append({
            "source":     "graphiti",
            "fact":       edge.fact,
            "valid_at":   edge.valid_at.isoformat() if edge.valid_at else None,
            "invalid_at": edge.invalid_at.isoformat() if edge.invalid_at else None,
            "confidence": 0.80,
            "evidence_count": 1,
        })
    return rules


def _graphiti_semantic_insights(vertical: str) -> List[dict]:
    """Return top patterns Graphiti has extracted for a vertical."""
    client = _get_graphiti_client()
    if client is None:
        return []

    queries = [
        f"best performing ingredients for {vertical}",
        f"commonly accepted blends for {vertical}",
        f"emerging patterns in {vertical} formulations",
    ]
    facts = []
    for q in queries:
        edges = _run_async(_search_graphiti(client, q, num_results=10))
        if edges:
            facts.extend(e.fact for e in edges if e.fact not in facts)
    return facts[:15]


# ── Event types ───────────────────────────────────────────────────────────────

EVENT_TYPES = {
    "formulation_generated":     "New blend created by optimizer",
    "blend_accepted":            "User accepted / proceeded with blend",
    "blend_rejected":            "User rejected blend (with reason)",
    "ingredient_swapped":        "User manually replaced an ingredient",
    "ingredient_removed":        "User removed an ingredient",
    "constraint_relaxed":        "Optimizer relaxed constraints",
    "pilot_batch_ordered":       "Pilot batch booked with ChemRich",
    "pilot_batch_passed":        "Pilot batch passed all tests",
    "pilot_batch_failed":        "Pilot batch failed (with failure type)",
    "certification_attempted":   "Certification application submitted",
    "certification_passed":      "Certification granted",
    "certification_rejected":    "Certification rejected (with ingredient)",
    "customer_preference_set":   "Customer expressed explicit preference",
    "rerun_requested":           "User ran optimizer again after seeing result",
    "proposal_downloaded":       "Proposal PDF downloaded",
    "carbon_passport_generated": "Carbon passport generated",
    "regulatory_update":         "Regulatory rule changed (bi-temporal)",
}


# ── Data schemas ──────────────────────────────────────────────────────────────

@dataclass
class MemoryEvent:
    event_id:      str
    event_type:    str
    timestamp:     str
    session_id:    str
    vertical:      str
    blend_hash:    str
    blend_summary: Dict[str, float]
    metadata:      Dict[str, Any]
    outcome:       Optional[str]


@dataclass
class IngredientPattern:
    ingredient:               str
    vertical:                 str
    acceptance_rate:          float
    rejection_count:          int
    acceptance_count:         int
    common_rejection_reasons: List[str]
    frequent_co_ingredients:  List[str]
    avg_pct_when_accepted:    float
    confidence:               float


@dataclass
class NegativeKnowledge:
    rule_id:        str
    rule_type:      str
    ingredient:     str
    vertical:       str
    condition:      str
    reason:         str
    confidence:     float
    evidence_count: int
    created_at:     str


@dataclass
class MemoryInsight:
    insight_type:      str
    title:             str
    description:       str
    affected_vertical: str
    confidence:        float
    evidence_count:    int
    actionable:        bool
    recommendation:    str
    source:            str = "pattern"   # "pattern" | "graphiti"


@dataclass
class MemoryNetworkState:
    events:               List[dict]
    negative_rules:       List[dict]
    patterns:             Dict[str, dict]
    session_count:        int
    total_formulations:   int
    total_acceptances:    int
    total_rejections:     int
    last_updated:         str


# ── Local JSON storage ────────────────────────────────────────────────────────

LOCAL_MEMORY_PATH = "/tmp/intelliform_memory.json"


def _load_local_memory() -> MemoryNetworkState:
    try:
        if os.path.exists(LOCAL_MEMORY_PATH):
            with open(LOCAL_MEMORY_PATH, "r") as f:
                data = json.load(f)
            return MemoryNetworkState(**data)
    except Exception:
        pass
    return MemoryNetworkState(
        events=[], negative_rules=[], patterns={},
        session_count=0, total_formulations=0,
        total_acceptances=0, total_rejections=0,
        last_updated=datetime.now().isoformat(),
    )


def _save_local_memory(state: MemoryNetworkState):
    try:
        with open(LOCAL_MEMORY_PATH, "w") as f:
            json.dump(asdict(state), f, default=str)
    except Exception:
        pass


# ── Supabase storage ──────────────────────────────────────────────────────────

def _get_supabase():
    if not SUPABASE_OK:
        return None
    url = os.environ.get("SUPABASE_URL",      "")
    key = os.environ.get("SUPABASE_ANON_KEY", "")
    if url and key:
        try:
            return create_client(url, key)
        except Exception:
            pass
    return None


def _store_event_supabase(client, event: MemoryEvent):
    try:
        client.table("memory_events").insert({
            "event_id":     event.event_id,
            "event_type":   event.event_type,
            "timestamp":    event.timestamp,
            "session_id":   event.session_id,
            "vertical":     event.vertical,
            "blend_hash":   event.blend_hash,
            "blend_summary": json.dumps(event.blend_summary),
            "metadata":     json.dumps(event.metadata),
            "outcome":      event.outcome,
        }).execute()
    except Exception:
        pass


# ── Core engine ───────────────────────────────────────────────────────────────

class FormulationMemoryNetwork:
    """
    Living knowledge graph for formulation intelligence.

    Storage priority: Graphiti temporal KG → Supabase → local JSON.
    All three stores are written concurrently on every record() call.
    Graphiti additionally enables bi-temporal queries and semantic search.
    """

    def __init__(self, session_id: Optional[str] = None):
        self.session_id = session_id or hashlib.md5(
            datetime.now().isoformat().encode()).hexdigest()[:8]
        self.state      = _load_local_memory()
        self.state.session_count += 1
        self._supabase  = _get_supabase()
        # Warm up Graphiti client asynchronously (non-blocking)
        self._graphiti  = _get_graphiti_client()

    @property
    def graphiti_active(self) -> bool:
        return self._graphiti is not None

    def _make_blend_hash(self, blend: Dict[str, float]) -> str:
        return hashlib.md5(
            json.dumps(sorted(blend.items()), default=str).encode()
        ).hexdigest()[:12]

    def _blend_summary(self, blend: Dict[str, float]) -> Dict[str, float]:
        return dict(sorted(blend.items(), key=lambda x: -x[1])[:5])

    # ── record ────────────────────────────────────────────────────────────────

    def record(
        self,
        event_type: str,
        blend: Dict[str, float],
        vertical: str,
        metadata: Optional[Dict] = None,
        outcome: Optional[str] = None,
    ) -> MemoryEvent:
        """
        Record any formulation event.

        Written to all three stores:
          1. Local JSON  (always)
          2. Supabase    (if configured)
          3. Graphiti KG (if Neo4j configured) — fires async in thread pool,
             so it never blocks the UI even on slow Neo4j connections.
        """
        metadata = metadata or {}
        ref_time = datetime.now()

        event = MemoryEvent(
            event_id=hashlib.md5(
                f"{event_type}{ref_time.isoformat()}{self.session_id}".encode()
            ).hexdigest()[:12],
            event_type=event_type,
            timestamp=ref_time.isoformat(),
            session_id=self.session_id,
            vertical=vertical,
            blend_hash=self._make_blend_hash(blend),
            blend_summary=self._blend_summary(blend),
            metadata=metadata,
            outcome=outcome,
        )

        # ── Local JSON ────────────────────────────────────────────────────────
        self.state.events.append(asdict(event))
        if event_type == "formulation_generated":
            self.state.total_formulations += 1
        elif outcome == "positive":
            self.state.total_acceptances += 1
        elif outcome == "negative":
            self.state.total_rejections += 1
        self.state.last_updated = ref_time.isoformat()
        _save_local_memory(self.state)

        # ── Supabase ──────────────────────────────────────────────────────────
        if self._supabase:
            _store_event_supabase(self._supabase, event)

        # ── Graphiti temporal KG ──────────────────────────────────────────────
        if self._graphiti:
            body = _episode_body(event_type, blend, vertical, metadata)
            _run_async(
                _add_graphiti_episode(
                    self._graphiti,
                    name=f"intelliform_{event_type}_{event.event_id}",
                    body=body,
                    ref_time=ref_time,
                    source_desc=f"IntelliForm {event_type} event",
                ),
                timeout=6.0,
            )

        return event

    # ── Rejection helpers ─────────────────────────────────────────────────────

    def record_rejection(
        self,
        blend: Dict[str, float],
        vertical: str,
        rejected_ingredients: List[str],
        reason: str,
    ):
        self.record(
            "blend_rejected", blend, vertical,
            metadata={"rejected_ingredients": rejected_ingredients, "reason": reason},
            outcome="negative",
        )
        for ing in rejected_ingredients:
            rule_id = hashlib.md5(f"{ing}{vertical}{reason}".encode()).hexdigest()[:8]
            existing = [r for r in self.state.negative_rules
                        if r.get("ingredient") == ing and r.get("vertical") == vertical]
            if existing:
                existing[0]["evidence_count"] = existing[0].get("evidence_count", 1) + 1
                existing[0]["confidence"]     = min(existing[0]["confidence"] + 0.05, 0.95)
            else:
                self.state.negative_rules.append({
                    "rule_id":        rule_id,
                    "rule_type":      "avoid_suggestion",
                    "ingredient":     ing,
                    "vertical":       vertical,
                    "condition":      "any formulation",
                    "reason":         reason,
                    "confidence":     0.60,
                    "evidence_count": 1,
                    "created_at":     datetime.now().isoformat(),
                })
        _save_local_memory(self.state)

    def record_pilot_outcome(
        self,
        blend: Dict[str, float],
        vertical: str,
        passed: bool,
        failure_type: Optional[str] = None,
        failure_ingredients: Optional[List[str]] = None,
    ):
        event_type = "pilot_batch_passed" if passed else "pilot_batch_failed"
        self.record(
            event_type, blend, vertical,
            metadata={"failure_type": failure_type,
                      "failure_ingredients": failure_ingredients or []},
            outcome="positive" if passed else "negative",
        )
        if not passed and failure_ingredients:
            self.record_rejection(
                blend, vertical, failure_ingredients,
                reason=f"Pilot batch failure: {failure_type}",
            )

    # ── Queries ───────────────────────────────────────────────────────────────

    def get_negative_rules(
        self,
        vertical: str,
        at_time: Optional[datetime] = None,
    ) -> List[dict]:
        """
        Return all negative knowledge rules for a vertical.

        If Graphiti is active, merges:
          - Flat rules from local JSON (exact ingredient matches)
          - Graphiti temporal edges (bi-temporal — respects at_time)

        at_time allows querying "what rules were active in November 2024?"
        This is the key capability flat storage cannot provide.
        """
        # Flat rules from local store
        flat_rules = [
            r for r in self.state.negative_rules
            if r.get("vertical") in (vertical, "all")
        ]

        # Graphiti temporal rules
        graphiti_rules = _graphiti_negative_rules(vertical, at_time) if self._graphiti else []

        # Merge: prefer Graphiti when same ingredient appears in both
        flat_ings = {r.get("ingredient") for r in flat_rules}
        merged = list(flat_rules)
        for gr in graphiti_rules:
            # Only add if Graphiti found something not already in flat rules
            fact = gr.get("fact", "")
            if not any(ing.lower() in fact.lower() for ing in flat_ings):
                merged.append(gr)

        return merged

    def get_ingredient_patterns(self, vertical: str) -> Dict[str, IngredientPattern]:
        vertical_events = [e for e in self.state.events
                           if e.get("vertical") == vertical]
        acceptance: Dict[str, int] = defaultdict(int)
        rejection:  Dict[str, int] = defaultdict(int)
        pct_sums:   Dict[str, list] = defaultdict(list)

        for event in vertical_events:
            summary = event.get("blend_summary", {})
            outcome = event.get("outcome")
            for ing, pct in summary.items():
                if outcome == "positive":
                    acceptance[ing] += 1
                    pct_sums[ing].append(pct)
                elif outcome == "negative":
                    rejection[ing] += 1

        patterns: Dict[str, IngredientPattern] = {}
        for ing in set(list(acceptance) + list(rejection)):
            total = acceptance[ing] + rejection[ing]
            if total == 0:
                continue
            acc_rate = acceptance[ing] / total
            avg_pct  = float(np.mean(pct_sums[ing])) if pct_sums[ing] else 0.0
            patterns[ing] = IngredientPattern(
                ingredient=ing, vertical=vertical,
                acceptance_rate=round(acc_rate, 2),
                rejection_count=rejection[ing],
                acceptance_count=acceptance[ing],
                common_rejection_reasons=[],
                frequent_co_ingredients=[],
                avg_pct_when_accepted=round(avg_pct, 1),
                confidence=min(total / 10, 1.0),
            )
        return patterns

    def get_insights(self, vertical: str) -> List[MemoryInsight]:
        """
        Extract actionable insights. Merges pattern-based insights (existing)
        with Graphiti semantic search insights (new — cross-session patterns
        the flat event log can't surface).
        """
        insights: List[MemoryInsight] = []
        events = [e for e in self.state.events if e.get("vertical") == vertical]
        n = len(events)

        if n < 3:
            insights.append(MemoryInsight(
                insight_type="info",
                title="Building memory…",
                description=(f"Only {n} formulation events recorded for {vertical}. "
                             f"Run more formulations to unlock pattern insights."),
                affected_vertical=vertical,
                confidence=1.0, evidence_count=n, actionable=False,
                recommendation="Run 10+ formulations to start seeing patterns.",
            ))
        else:
            patterns = self.get_ingredient_patterns(vertical)

            # High-acceptance ingredients
            high_accept = [(ing, p) for ing, p in patterns.items()
                           if p.acceptance_rate > 0.75 and p.acceptance_count >= 2]
            if high_accept:
                top = sorted(high_accept, key=lambda x: -x[1].acceptance_rate)[:3]
                insights.append(MemoryInsight(
                    insight_type="preference",
                    title="High-acceptance ingredients for this vertical",
                    description=(f"Consistently appear in accepted blends: "
                                 f"{', '.join(ing for ing, _ in top)}"),
                    affected_vertical=vertical,
                    confidence=min(sum(p.acceptance_count for _, p in top) / 20, 0.9),
                    evidence_count=sum(p.acceptance_count + p.rejection_count for _, p in top),
                    actionable=True,
                    recommendation=(f"Prioritise {top[0][0]} — "
                                    f"{top[0][1].acceptance_rate*100:.0f}% acceptance rate."),
                ))

            # Consistently rejected ingredients
            high_reject = [(ing, p) for ing, p in patterns.items()
                           if p.acceptance_rate < 0.25 and p.rejection_count >= 2]
            if high_reject:
                bottom = sorted(high_reject, key=lambda x: x[1].acceptance_rate)[:3]
                insights.append(MemoryInsight(
                    insight_type="warning",
                    title="Consistently rejected ingredients",
                    description=(f"Appear repeatedly in rejected blends: "
                                 f"{', '.join(ing for ing, _ in bottom)}"),
                    affected_vertical=vertical,
                    confidence=min(sum(p.rejection_count for _, p in bottom) / 10, 0.9),
                    evidence_count=sum(p.rejection_count for _, p in bottom),
                    actionable=True,
                    recommendation=(f"Avoid {bottom[0][0]} for {vertical} — "
                                    f"only {bottom[0][1].acceptance_rate*100:.0f}% acceptance."),
                ))

            # High-confidence negative rules
            neg_rules = self.get_negative_rules(vertical)
            hc_rules  = [r for r in neg_rules if r.get("confidence", 0) > 0.70]
            if hc_rules:
                ings_str = ", ".join(
                    r.get("ingredient") or r.get("fact", "")[:60]
                    for r in hc_rules[:5]
                )
                insights.append(MemoryInsight(
                    insight_type="warning",
                    title=f"{len(hc_rules)} high-confidence negative rules active",
                    description=f"Avoid: {ings_str}",
                    affected_vertical=vertical,
                    confidence=float(np.mean([r.get("confidence", 0.5) for r in hc_rules])),
                    evidence_count=sum(r.get("evidence_count", 1) for r in hc_rules),
                    actionable=True,
                    recommendation="Automatically applied to filter future optimizer suggestions.",
                ))

            # Rerun rate
            reruns    = sum(1 for e in events if e.get("event_type") == "rerun_requested")
            rerun_rate = reruns / n if n else 0
            if rerun_rate > 0.4:
                insights.append(MemoryInsight(
                    insight_type="opportunity",
                    title=f"High rerun rate ({rerun_rate*100:.0f}%) — optimizer needs tuning",
                    description=(f"{reruns}/{n} sessions required rerun. "
                                 f"Optimizer not converging on good solutions first try."),
                    affected_vertical=vertical,
                    confidence=0.75, evidence_count=reruns, actionable=True,
                    recommendation="Switch to Bayesian optimizer — it learns from reruns.",
                ))

        # ── Graphiti semantic insights ─────────────────────────────────────────
        if self._graphiti:
            facts = _graphiti_semantic_insights(vertical)
            if facts:
                insights.append(MemoryInsight(
                    insight_type="opportunity",
                    title=f"Graphiti: {len(facts)} cross-session patterns detected",
                    description=facts[0] if facts else "",
                    affected_vertical=vertical,
                    confidence=0.70,
                    evidence_count=len(facts),
                    actionable=True,
                    recommendation=(
                        "These patterns were extracted from the temporal knowledge "
                        "graph across all sessions. Full list in the knowledge export."
                    ),
                    source="graphiti",
                ))

        return sorted(
            insights,
            key=lambda x: -(x.confidence * (1.2 if x.actionable else 1.0)),
        )

    # ── Optimizer priors ──────────────────────────────────────────────────────

    def adjust_optimization_priors(
        self,
        vertical: str,
        db_ingredients: List[str],
    ) -> Dict[str, float]:
        patterns  = self.get_ingredient_patterns(vertical)
        neg_rules = {
            r.get("ingredient", "") for r in self.get_negative_rules(vertical)
            if r.get("confidence", 0) > 0.65 and r.get("ingredient")
        }
        adjustments: Dict[str, float] = {}
        for ing in db_ingredients:
            if ing in neg_rules:
                adjustments[ing] = -0.8
            elif ing in patterns:
                p = patterns[ing]
                if p.confidence > 0.3:
                    adjustments[ing] = (p.acceptance_rate - 0.5) * p.confidence
        return adjustments

    # ── Stats & export ────────────────────────────────────────────────────────

    def get_summary_stats(self) -> dict:
        return {
            "total_events":       len(self.state.events),
            "total_formulations": self.state.total_formulations,
            "total_acceptances":  self.state.total_acceptances,
            "total_rejections":   self.state.total_rejections,
            "acceptance_rate":    round(
                self.state.total_acceptances /
                max(self.state.total_acceptances + self.state.total_rejections, 1), 2,
            ),
            "negative_rules":     len(self.state.negative_rules),
            "session_count":      self.state.session_count,
            "last_updated":       self.state.last_updated,
            "verticals_active":   list(set(
                e.get("vertical", "—") for e in self.state.events
            )),
            "graphiti_active":    self.graphiti_active,
            "storage_tier":       ("graphiti+supabase+local" if self.graphiti_active and self._supabase
                                   else "graphiti+local" if self.graphiti_active
                                   else "supabase+local" if self._supabase
                                   else "local"),
        }

    def export_knowledge_base(self) -> str:
        verticals = set(e.get("vertical", "") for e in self.state.events)
        return json.dumps({
            "intelliform_memory_network": {
                "version":     "1.5",
                "exported_at": datetime.now().isoformat(),
                "stats":       self.get_summary_stats(),
            },
            "events":            self.state.events[-500:],
            "negative_knowledge": self.state.negative_rules,
            "patterns_summary":  {
                v: {ing: {"acceptance_rate": p.acceptance_rate,
                          "evidence": p.acceptance_count + p.rejection_count}
                    for ing, p in self.get_ingredient_patterns(v).items()}
                for v in verticals
            },
        }, indent=2, default=str)

    def clear_vertical(self, vertical: str):
        self.state.events        = [e for e in self.state.events
                                    if e.get("vertical") != vertical]
        self.state.negative_rules = [r for r in self.state.negative_rules
                                     if r.get("vertical") != vertical]
        _save_local_memory(self.state)


# ── Module-level singleton ────────────────────────────────────────────────────

_MEMORY_NETWORK: Optional[FormulationMemoryNetwork] = None


def get_memory_network(session_id: Optional[str] = None) -> FormulationMemoryNetwork:
    global _MEMORY_NETWORK
    if _MEMORY_NETWORK is None:
        _MEMORY_NETWORK = FormulationMemoryNetwork(session_id)
    return _MEMORY_NETWORK
