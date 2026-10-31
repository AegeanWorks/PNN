"""Tiering and ranking for the PNN Phase 2 property screen."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

from . import rules
from .compounds import Compound
from .rules import CriterionResult

TIER_ORDER = {"A": 0, "B": 1, "C": 2, "Reject": 3}
TIER_A_SCORE = 0.80
TIER_B_SCORE = 0.60
MAX_CORE_FAILURES_FOR_TRIAGE = 3
MAX_UNKNOWNS_FOR_TIER_A = 2


@dataclass(frozen=True)
class ScreenResult:
    """A compound, its verdicts, and the resulting triage tier."""

    compound: Compound
    tier: str
    score: float
    criteria: Tuple[CriterionResult, ...]

    def failed_criteria(self) -> List[str]:
        return [item.name for item in self.criteria if item.status == rules.FAIL]

    def unknown_criteria(self) -> List[str]:
        return [item.name for item in self.criteria if item.status == rules.UNKNOWN]


def _weighted_score(criteria: Sequence[CriterionResult]) -> float:
    scored = [item for item in criteria if item.status != rules.UNKNOWN]
    total_weight = sum(rules.weight(item.name) for item in scored)
    if total_weight == 0:
        return 0.0
    return sum(item.score * rules.weight(item.name) for item in scored) / total_weight


def _tier(criteria: Sequence[CriterionResult], score: float) -> str:
    hard_failures = [
        item for item in criteria if item.name in rules.HARD_REJECT_CRITERIA and item.status == rules.FAIL
    ]
    core_failures = [
        item for item in criteria if item.name in rules.CORE_CRITERIA and item.status == rules.FAIL
    ]
    unknowns = [item for item in criteria if item.status == rules.UNKNOWN]
    if hard_failures or len(core_failures) >= MAX_CORE_FAILURES_FOR_TRIAGE:
        return "Reject"
    if not core_failures and score >= TIER_A_SCORE and len(unknowns) <= MAX_UNKNOWNS_FOR_TIER_A:
        return "A"
    if not core_failures and score >= TIER_B_SCORE:
        return "B"
    return "C"


def evaluate_compound(compound: Compound) -> ScreenResult:
    """Evaluate every rule for one compound and assign its tier."""
    criteria = tuple(rules.evaluate(compound))
    score = _weighted_score(criteria)
    return ScreenResult(compound=compound, tier=_tier(criteria, score), score=score, criteria=criteria)


def screen_compounds(compounds: Iterable[Compound]) -> List[ScreenResult]:
    """Evaluate and rank compounds by tier, then score (highest first), then ID."""
    results = [evaluate_compound(compound) for compound in compounds]
    return sorted(results, key=lambda item: (TIER_ORDER[item.tier], -item.score, item.compound.compound_id))


def summarize(results: Sequence[ScreenResult]) -> Dict[str, object]:
    """Count tiers and report the scored average, including rejected compounds."""
    counts = {tier: 0 for tier in TIER_ORDER}
    for result in results:
        counts[result.tier] += 1
    average = sum(result.score for result in results) / len(results) if results else 0.0
    return {"total": len(results), "by_tier": counts, "average_score": round(average, 4)}
