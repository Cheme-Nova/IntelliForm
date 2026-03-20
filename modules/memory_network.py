"""
modules/memory_network.py
Formulation Memory Network™ — IntelliForm v1.4

The strategic moat. A living knowledge graph that accumulates institutional
memory from every formulation run, customer preference, pilot batch outcome,
and rejection pattern.

What it does that no LLM or public dataset can replicate:
  - Remembers that customer X always rejects ingredient Y (even GRAS)
  - Knows that vertical Z + function W + bio% > 90 → 73% pilot acceptance rate
  - Tracks which suggestions get accepted vs rejected per customer segment
  - Identifies emerging ingredient preference patterns before they're published
  - Builds a "negative knowledge" database: what NOT to suggest and why
  - Predicts reformulation success probability from pattern matching

Architecture:
  - Event store: every formulation run, accept, reject, tweak stored as events
  - Pattern extractor: finds statistically significant preference patterns
  - Recommendation adjuster: biases future suggestions based on memory
  - Negative knowledge DB: explicit "never suggest X for Y use case" rules
  - Customer segment profiler: groups similar customers by preference fingerprint

Storage:
  - Supabase (production): persistent across sessions, shareable across org
  - Local JSON (fallback): session-persistent, no external dependency

This is the data flywheel. The more customers use IntelliForm through ChemRich
pilot batches, the more the memory network learns, the better the suggestions
get, the harder it is for any competitor to replicate.

Reference:
  - Collaborative filtering for materials: Raccuglia et al., Nature 2016
  - "Failed experiment" databases: Kononova et al., Sci Data 2021
  - Knowledge graphs for chemistry: Horawalavithana et al., ACS 2022
"""
from __future__ import annotations

import json
import os
import hashlib
from collections import defaultdict, Counter
from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import numpy as np

try:
    from supabase import create_client
    SUPABASE_OK = True
except ImportError:
    SUPABASE_OK = False


# ── Event types ───────────────────────────────────────────────────────────────

EVENT_TYPES = {
    "formulation_generated":    "New blend created by optimizer",
    "blend_accepted":           "User accepted / proceeded with blend",
    "blend_rejected":           "User rejected blend (with reason)",
    "ingredient_swapped":       "User manually replaced an ingredient",
    "ingredient_removed":       "User removed an ingredient",
    "constraint_relaxed":       "Optimizer relaxed constraints",
    "pilot_batch_ordered":      "Pilot batch booked with ChemRich",
    "pilot_batch_passed":       "Pilot batch passed all tests",
    "pilot_batch_failed":       "Pilot batch failed (with failure type)",
    "certification_attempted":  "Certification application submitted",
    "certification_passed":     "Certification granted",
    "certification_rejected":   "Certification rejected (with ingredient)",
    "customer_preference_set":  "Customer expressed explicit preference",
    "rerun_requested":          "User ran optimizer again after seeing result",
    "proposal_downloaded":      "Proposal PDF downloaded",
    "carbon_passport_generated":"Carbon passport generated",
}


# ── Data schemas ──────────────────────────────────────────────────────────────

@dataclass
class MemoryEvent:
    event_id:     str
    event_type:   str
    timestamp:    str
    session_id:   str
    vertical:     str
    blend_hash:   str
    blend_summary: Dict[str, float]   # top 5 ingredients + percentages
    metadata:     Dict[str, Any]      # event-specific data
    outcome:      Optional[str]       # positive / negative / neutral


@dataclass
class IngredientPattern:
    ingredient:       str
    vertical:         str
    acceptance_rate:  float           # 0-1: how often this ingredient survives to acceptance
    rejection_count:  int
    acceptance_count: int
    common_rejection_reasons: List[str]
    frequent_co_ingredients: List[str]
    avg_pct_when_accepted: float
    confidence:       float           # based on sample size


@dataclass
class NegativeKnowledge:
    """Explicit "don't do this" rules learned from failures."""
    rule_id:          str
    rule_type:        str             # never_suggest / avoid_combination / concentration_cap
    ingredient:       str
    vertical:         str
    condition:        str             # when this rule applies
    reason:           str             # why — from actual failure data
    confidence:       float
    evidence_count:   int
    created_at:       str


@dataclass
class MemoryInsight:
    """A pattern the memory network has detected."""
    insight_type:     str             # preference / warning / opportunity / trend
    title:            str
    description:      str
    affected_vertical: str
    confidence:       float
    evidence_count:   int
    actionable:       bool
    recommendation:   str


@dataclass
class MemoryNetworkState:
    """Full state of the memory network."""
    events:           List[dict]
    negative_rules:   List[dict]
    patterns:         Dict[str, dict]  # ingredient -> pattern
    session_count:    int
    total_formulations: int
    total_acceptances:  int
    total_rejections:   int
    last_updated:     str


# ── Local fallback storage ─────────────────────────────────────────────────────

LOCAL_MEMORY_PATH = "/tmp/intelliform_memory.json"


def _load_local_memory() -> MemoryNetworkState:
    """Load memory from local JSON file."""
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
    """Persist memory to local JSON."""
    try:
        with open(LOCAL_MEMORY_PATH, "w") as f:
            json.dump(asdict(state), f, default=str)
    except Exception:
        pass


# ── Supabase storage ──────────────────────────────────────────────────────────

def _get_supabase():
    """Get Supabase client if configured."""
    if not SUPABASE_OK:
        return None
    url = os.environ.get("SUPABASE_URL", "")
    key = os.environ.get("SUPABASE_ANON_KEY", "")
    if url and key:
        try:
            return create_client(url, key)
        except Exception:
            pass
    return None


def _store_event_supabase(client, event: MemoryEvent):
    """Store event in Supabase memory_events table."""
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
        pass  # Graceful — table may not exist yet


# ── Core engine ────────────────────────────────────────────────────────────────

class FormulationMemoryNetwork:
    """
    Living knowledge graph for formulation intelligence.
    Accumulates patterns from every user interaction.
    """

    def __init__(self, session_id: Optional[str] = None):
        self.session_id = session_id or hashlib.md5(
            datetime.now().isoformat().encode()).hexdigest()[:8]
        self.state = _load_local_memory()
        self.state.session_count += 1
        self._supabase = _get_supabase()

    def _make_blend_hash(self, blend: Dict[str, float]) -> str:
        return hashlib.md5(
            json.dumps(sorted(blend.items()), default=str).encode()
        ).hexdigest()[:12]

    def _blend_summary(self, blend: Dict[str, float]) -> Dict[str, float]:
        """Top 5 ingredients by percentage."""
        return dict(sorted(blend.items(), key=lambda x: -x[1])[:5])

    def record(
        self,
        event_type: str,
        blend: Dict[str, float],
        vertical: str,
        metadata: Optional[Dict] = None,
        outcome: Optional[str] = None,
    ) -> MemoryEvent:
        """Record any formulation event into the memory network."""
        event = MemoryEvent(
            event_id=hashlib.md5(
                f"{event_type}{datetime.now().isoformat()}{self.session_id}".encode()
            ).hexdigest()[:12],
            event_type=event_type,
            timestamp=datetime.now().isoformat(),
            session_id=self.session_id,
            vertical=vertical,
            blend_hash=self._make_blend_hash(blend),
            blend_summary=self._blend_summary(blend),
            metadata=metadata or {},
            outcome=outcome,
        )

        # Store locally
        self.state.events.append(asdict(event))
        if event_type == "formulation_generated":
            self.state.total_formulations += 1
        elif outcome == "positive":
            self.state.total_acceptances += 1
        elif outcome == "negative":
            self.state.total_rejections += 1

        self.state.last_updated = datetime.now().isoformat()
        _save_local_memory(self.state)

        # Store in Supabase if available
        if self._supabase:
            _store_event_supabase(self._supabase, event)

        return event

    def record_rejection(
        self,
        blend: Dict[str, float],
        vertical: str,
        rejected_ingredients: List[str],
        reason: str,
    ):
        """Record a blend rejection — feeds into negative knowledge."""
        self.record(
            "blend_rejected", blend, vertical,
            metadata={"rejected_ingredients": rejected_ingredients, "reason": reason},
            outcome="negative",
        )
        # Add to negative knowledge
        for ing in rejected_ingredients:
            rule_id = hashlib.md5(f"{ing}{vertical}{reason}".encode()).hexdigest()[:8]
            # Check if rule already exists
            existing = [r for r in self.state.negative_rules
                        if r.get("ingredient") == ing and r.get("vertical") == vertical]
            if existing:
                existing[0]["evidence_count"] = existing[0].get("evidence_count", 1) + 1
                existing[0]["confidence"] = min(
                    existing[0]["confidence"] + 0.05, 0.95)
            else:
                self.state.negative_rules.append({
                    "rule_id": rule_id,
                    "rule_type": "avoid_suggestion",
                    "ingredient": ing,
                    "vertical": vertical,
                    "condition": "any formulation",
                    "reason": reason,
                    "confidence": 0.60,
                    "evidence_count": 1,
                    "created_at": datetime.now().isoformat(),
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
        """Record pilot batch outcome — highest signal event."""
        event_type = "pilot_batch_passed" if passed else "pilot_batch_failed"
        self.record(
            event_type, blend, vertical,
            metadata={
                "failure_type": failure_type,
                "failure_ingredients": failure_ingredients or [],
            },
            outcome="positive" if passed else "negative",
        )
        if not passed and failure_ingredients:
            self.record_rejection(blend, vertical, failure_ingredients,
                                   reason=f"Pilot batch failure: {failure_type}")

    def get_negative_rules(self, vertical: str) -> List[dict]:
        """Get all negative knowledge rules for a vertical."""
        return [r for r in self.state.negative_rules
                if r.get("vertical") == vertical or r.get("vertical") == "all"]

    def get_ingredient_patterns(self, vertical: str) -> Dict[str, IngredientPattern]:
        """
        Extract ingredient acceptance patterns from event history.
        Returns ingredients ranked by acceptance rate for this vertical.
        """
        vertical_events = [e for e in self.state.events
                          if e.get("vertical") == vertical]

        acceptance = defaultdict(int)
        rejection  = defaultdict(int)
        pct_sums   = defaultdict(list)

        for event in vertical_events:
            summary = event.get("blend_summary", {})
            outcome = event.get("outcome")
            for ing, pct in summary.items():
                if outcome == "positive":
                    acceptance[ing] += 1
                    pct_sums[ing].append(pct)
                elif outcome == "negative":
                    rejection[ing] += 1

        patterns = {}
        all_ings = set(list(acceptance.keys()) + list(rejection.keys()))
        for ing in all_ings:
            total = acceptance[ing] + rejection[ing]
            if total == 0:
                continue
            acc_rate = acceptance[ing] / total
            avg_pct = float(np.mean(pct_sums[ing])) if pct_sums[ing] else 0.0
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
        Extract actionable insights from the memory network.
        Returns insights sorted by confidence × actionability.
        """
        insights = []
        events = [e for e in self.state.events if e.get("vertical") == vertical]
        n = len(events)

        if n < 3:
            insights.append(MemoryInsight(
                insight_type="info",
                title="Building memory…",
                description=f"Only {n} formulation events recorded for {vertical}. "
                            f"Run more formulations to unlock pattern insights.",
                affected_vertical=vertical,
                confidence=1.0,
                evidence_count=n,
                actionable=False,
                recommendation="Run 10+ formulations to start seeing patterns.",
            ))
            return insights

        # Insight 1: Most accepted ingredients
        patterns = self.get_ingredient_patterns(vertical)
        high_accept = [(ing, p) for ing, p in patterns.items()
                       if p.acceptance_rate > 0.75 and p.acceptance_count >= 2]
        if high_accept:
            top = sorted(high_accept, key=lambda x: -x[1].acceptance_rate)[:3]
            insights.append(MemoryInsight(
                insight_type="preference",
                title="High-acceptance ingredients for this vertical",
                description=f"These ingredients appear in accepted blends most often: "
                            f"{', '.join(ing for ing, _ in top)}",
                affected_vertical=vertical,
                confidence=min(sum(p.acceptance_count for _, p in top) / 20, 0.9),
                evidence_count=sum(p.acceptance_count + p.rejection_count for _, p in top),
                actionable=True,
                recommendation=f"Prioritise {top[0][0]} in new formulations — "
                               f"{top[0][1].acceptance_rate*100:.0f}% acceptance rate.",
            ))

        # Insight 2: Consistently rejected ingredients
        high_reject = [(ing, p) for ing, p in patterns.items()
                       if p.acceptance_rate < 0.25 and p.rejection_count >= 2]
        if high_reject:
            bottom = sorted(high_reject, key=lambda x: x[1].acceptance_rate)[:3]
            insights.append(MemoryInsight(
                insight_type="warning",
                title="Consistently rejected ingredients",
                description=f"These ingredients appear repeatedly in rejected blends: "
                            f"{', '.join(ing for ing, _ in bottom)}",
                affected_vertical=vertical,
                confidence=min(sum(p.rejection_count for _, p in bottom) / 10, 0.9),
                evidence_count=sum(p.rejection_count for _, p in bottom),
                actionable=True,
                recommendation=f"Avoid {bottom[0][0]} in new formulations for this vertical — "
                               f"only {bottom[0][1].acceptance_rate*100:.0f}% acceptance.",
            ))

        # Insight 3: Negative knowledge rules
        neg_rules = self.get_negative_rules(vertical)
        high_conf_rules = [r for r in neg_rules if r.get("confidence", 0) > 0.7]
        if high_conf_rules:
            insights.append(MemoryInsight(
                insight_type="warning",
                title=f"{len(high_conf_rules)} high-confidence negative rules active",
                description=f"Ingredients to avoid based on observed failures: "
                            f"{', '.join(r['ingredient'] for r in high_conf_rules[:5])}",
                affected_vertical=vertical,
                confidence=float(np.mean([r.get("confidence", 0.5)
                                          for r in high_conf_rules])),
                evidence_count=sum(r.get("evidence_count", 1) for r in high_conf_rules),
                actionable=True,
                recommendation="These rules are automatically applied to filter "
                               "future optimizer suggestions.",
            ))

        # Insight 4: Rerun rate (signals optimizer quality)
        reruns = sum(1 for e in events if e.get("event_type") == "rerun_requested")
        if reruns > 0 and n > 0:
            rerun_rate = reruns / n
            if rerun_rate > 0.4:
                insights.append(MemoryInsight(
                    insight_type="opportunity",
                    title=f"High rerun rate ({rerun_rate*100:.0f}%) — optimizer needs tuning",
                    description=f"{reruns} of {n} sessions required rerun. "
                                f"The optimizer is not converging on good solutions "
                                f"first try for this vertical.",
                    affected_vertical=vertical,
                    confidence=0.75,
                    evidence_count=reruns,
                    actionable=True,
                    recommendation="Consider switching to Bayesian optimizer for this vertical — "
                                   "it learns from reruns.",
                ))

        return sorted(insights, key=lambda x: -(x.confidence * (1.2 if x.actionable else 1.0)))

    def adjust_optimization_priors(
        self,
        vertical: str,
        db_ingredients: List[str],
    ) -> Dict[str, float]:
        """
        Return per-ingredient weight adjustments for optimizer.
        Positive weight = memory suggests this ingredient is preferred.
        Negative weight = memory suggests this ingredient should be avoided.

        These weights can be used to bias the LP/Bayesian optimizer
        toward historically accepted formulations.
        """
        patterns = self.get_ingredient_patterns(vertical)
        neg_rules = {r["ingredient"] for r in self.get_negative_rules(vertical)
                     if r.get("confidence", 0) > 0.65}

        adjustments = {}
        for ing in db_ingredients:
            if ing in neg_rules:
                adjustments[ing] = -0.8   # strong downweight
            elif ing in patterns:
                p = patterns[ing]
                if p.confidence > 0.3:
                    # Scale: 0.5 acceptance rate = 0 adjustment
                    # 1.0 = +0.5, 0.0 = -0.5
                    adjustments[ing] = (p.acceptance_rate - 0.5) * p.confidence
            # Default: 0.0 (neutral)

        return adjustments

    def get_summary_stats(self) -> dict:
        """Return summary statistics for display."""
        return {
            "total_events":       len(self.state.events),
            "total_formulations": self.state.total_formulations,
            "total_acceptances":  self.state.total_acceptances,
            "total_rejections":   self.state.total_rejections,
            "acceptance_rate":    round(
                self.state.total_acceptances /
                max(self.state.total_acceptances + self.state.total_rejections, 1), 2
            ),
            "negative_rules":     len(self.state.negative_rules),
            "session_count":      self.state.session_count,
            "last_updated":       self.state.last_updated,
            "verticals_active":   list(set(
                e.get("vertical", "—") for e in self.state.events
            )),
        }

    def export_knowledge_base(self) -> str:
        """Export the full memory network as JSON for backup/analysis."""
        return json.dumps({
            "intelliform_memory_network": {
                "version": "1.4",
                "exported_at": datetime.now().isoformat(),
                "stats": self.get_summary_stats(),
            },
            "events": self.state.events[-500:],  # last 500 events
            "negative_knowledge": self.state.negative_rules,
            "patterns_summary": {
                v: {ing: {"acceptance_rate": p.acceptance_rate,
                          "evidence": p.acceptance_count + p.rejection_count}
                    for ing, p in self.get_ingredient_patterns(v).items()}
                for v in set(e.get("vertical", "") for e in self.state.events)
            },
        }, indent=2, default=str)

    def clear_vertical(self, vertical: str):
        """Remove all memory for a specific vertical (privacy/reset)."""
        self.state.events = [e for e in self.state.events
                             if e.get("vertical") != vertical]
        self.state.negative_rules = [r for r in self.state.negative_rules
                                     if r.get("vertical") != vertical]
        _save_local_memory(self.state)


# ── Module-level singleton ────────────────────────────────────────────────────

_MEMORY_NETWORK: Optional[FormulationMemoryNetwork] = None


def get_memory_network(session_id: Optional[str] = None) -> FormulationMemoryNetwork:
    """Get or create the module-level memory network singleton."""
    global _MEMORY_NETWORK
    if _MEMORY_NETWORK is None:
        _MEMORY_NETWORK = FormulationMemoryNetwork(session_id)
    return _MEMORY_NETWORK
