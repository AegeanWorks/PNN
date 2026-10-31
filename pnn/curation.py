"""Merge, deduplicate and diversify compound libraries for PNN Phase 2."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .compounds import Compound
from .screen import TIER_ORDER, ScreenResult, evaluate_compound

_TIER_CEILING = {"A": 0, "B": 1, "C": 2}
_PROPERTY_RANGES = (
    (200.0, 450.0),
    (1.5, 3.5),
    (0.0, 70.0),
    (0.0, 3.0),
    (0.0, 7.0),
    (0.0, 7.0),
)


class CurationError(ValueError):
    """Raised when libraries cannot be merged or curated safely."""


@dataclass(frozen=True)
class Candidate:
    """A compound with its source file and screening verdict."""

    compound: Compound
    source: str
    result: ScreenResult


def _preference(candidate: Candidate) -> tuple:
    return (
        TIER_ORDER[candidate.result.tier],
        -candidate.result.score,
        not bool(candidate.compound.smiles),
        candidate.compound.compound_id,
    )


def _merge(left: Candidate, right: Candidate) -> Candidate:
    """Keep the better record, filling a blank SMILES/name from the other source."""
    winner = min(left, right, key=_preference)
    loser = right if winner is left else left
    compound = winner.compound
    if not compound.smiles and loser.compound.smiles:
        compound = replace(compound, smiles=loser.compound.smiles)
    if not compound.name and loser.compound.name:
        compound = replace(compound, name=loser.compound.name)
    return winner if compound is winner.compound else replace(winner, compound=compound)


def prepare_candidates(inputs: Iterable[Tuple[str, Sequence[Compound]]]) -> List[Candidate]:
    """Screen every input compound while keeping its source provenance."""
    candidates: List[Candidate] = []
    for source, compounds in inputs:
        for compound in compounds:
            candidates.append(Candidate(compound=compound, source=source, result=evaluate_compound(compound)))
    if not candidates:
        raise CurationError("no compounds were supplied")
    return candidates


def deduplicate(candidates: Sequence[Candidate]) -> Tuple[List[Candidate], Dict[str, int]]:
    """Collapse duplicate IDs first, then exact non-empty SMILES strings."""
    by_id: Dict[str, Candidate] = {}
    duplicate_ids = 0
    for candidate in candidates:
        compound_id = candidate.compound.compound_id
        existing = by_id.get(compound_id)
        if existing is None:
            by_id[compound_id] = candidate
            continue
        left, right = existing.compound.smiles, candidate.compound.smiles
        if left and right and left != right:
            raise CurationError(
                f"compound_id '{compound_id}' has two different SMILES "
                f"('{left}' from {existing.source}, '{right}' from {candidate.source})"
            )
        by_id[compound_id] = _merge(existing, candidate)
        duplicate_ids += 1

    by_smiles: Dict[str, Candidate] = {}
    without_smiles: List[Candidate] = []
    duplicate_smiles = 0
    for candidate in by_id.values():
        smiles = candidate.compound.smiles.strip()
        if not smiles:
            without_smiles.append(candidate)
            continue
        existing = by_smiles.get(smiles)
        if existing is None:
            by_smiles[smiles] = candidate
            continue
        by_smiles[smiles] = _merge(existing, candidate)
        duplicate_smiles += 1

    return (
        without_smiles + list(by_smiles.values()),
        {"duplicate_ids": duplicate_ids, "duplicate_smiles": duplicate_smiles},
    )


def property_bin(compound: Compound, width: float) -> Tuple[int, ...]:
    """Bin the six core properties on a normalized 0..1 grid."""
    bins = max(1, round(1.0 / width))
    cells = []
    for value, (low, high) in zip(
        (compound.mw, compound.logp, compound.tpsa, compound.hbd, compound.hba, compound.rotatable_bonds),
        _PROPERTY_RANGES,
    ):
        normalized = min(1.0, max(0.0, (value - low) / (high - low)))
        cells.append(min(int(normalized * bins), bins - 1))
    return tuple(cells)


def _validate(min_tier: str, per_bin: int, bin_width: float, max_compounds: Optional[int]) -> None:
    if min_tier not in _TIER_CEILING:
        raise CurationError(f"unknown screening tier: {min_tier}")
    if per_bin < 1:
        raise CurationError("per-bin cap must be at least 1")
    if not 0.0 < bin_width <= 1.0:
        raise CurationError("bin width must be greater than 0 and at most 1")
    if max_compounds is not None and max_compounds < 1:
        raise CurationError("max compounds must be at least 1")


def curate(
    candidates: Sequence[Candidate],
    min_tier: str = "B",
    per_bin: int = 25,
    bin_width: float = 0.1,
    max_compounds: Optional[int] = None,
) -> Tuple[List[Candidate], Dict[str, object]]:
    """Filter by tier, cap each property bin, then apply the global size limit."""
    _validate(min_tier, per_bin, bin_width, max_compounds)
    ceiling = _TIER_CEILING[min_tier]
    eligible = [
        candidate
        for candidate in candidates
        if TIER_ORDER[candidate.result.tier] <= ceiling
    ]
    skipped_by_tier = len(candidates) - len(eligible)

    groups: Dict[Tuple[int, ...], List[Candidate]] = {}
    for candidate in eligible:
        groups.setdefault(property_bin(candidate.compound, bin_width), []).append(candidate)

    kept: List[Candidate] = []
    for group in groups.values():
        kept.extend(sorted(group, key=_preference)[:per_bin])
    kept.sort(key=_preference)
    truncated = len(kept) > max_compounds if max_compounds is not None else False
    if max_compounds is not None:
        kept = kept[:max_compounds]

    summary = {
        "eligible": len(eligible),
        "skipped_by_tier": skipped_by_tier,
        "bins": len(groups),
        "per_bin": per_bin,
        "bin_width": bin_width,
        "max_compounds": max_compounds,
        "selected": len(kept),
        "truncated_by_max": truncated,
    }
    return kept, summary
