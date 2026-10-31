"""Docking-score ingestion and consensus ranking for PNN Phase 2."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

SCORE_FIELDS = ("vina_affinity", "glide_score", "mmgbsa_dg")
WEIGHTS = {"vina_affinity": 0.5, "glide_score": 0.3, "mmgbsa_dg": 0.2}
TIER_ORDER = {"A": 0, "B": 1, "C": 2, "Reject": 3}

_FIELD_ALIASES = {
    "compoundid": "compound_id",
    "id": "compound_id",
    "compound": "compound_id",
    "target": "target",
    "targetname": "target",
    "receptor": "target",
    "vina": "vina_affinity",
    "vinaaffinity": "vina_affinity",
    "affinity": "vina_affinity",
    "glide": "glide_score",
    "glidescore": "glide_score",
    "mmgbsa": "mmgbsa_dg",
    "mmgbsadg": "mmgbsa_dg",
    "deltag": "mmgbsa_dg",
    "deltagbind": "mmgbsa_dg",
    "rmsd": "rmsd",
    "rmsdligand": "rmsd",
    "hbondcount": "hbond_count",
    "hbonds": "hbond_count",
    "keyinteractions": "key_interactions",
    "interactions": "key_interactions",
    "tier": "tier",
    "score": "score",
    "name": "name",
    "compoundname": "name",
    "smiles": "smiles",
    "canonicalsmiles": "smiles",
}
_REQUIRED_FIELDS = ("compound_id", "target")
_NUMERIC_FIELDS = SCORE_FIELDS + ("rmsd", "hbond_count")
_NON_NEGATIVE_FIELDS = frozenset({"rmsd", "hbond_count"})


class DockingError(ValueError):
    """Raised when docking input cannot be interpreted."""


@dataclass(frozen=True)
class DockingScore:
    """One docking result for one compound against one target."""

    compound_id: str
    target: str
    vina_affinity: Optional[float] = None
    glide_score: Optional[float] = None
    mmgbsa_dg: Optional[float] = None
    rmsd: Optional[float] = None
    hbond_count: Optional[float] = None
    key_interactions: str = ""


@dataclass(frozen=True)
class DockingRank:
    """A docking score with its normalized consensus quality."""

    score: DockingScore
    consensus: float
    coverage: float


@dataclass(frozen=True)
class CompoundRank:
    """The best target row for one compound and how many targets it reached."""

    compound_id: str
    best: DockingRank
    targets_docked: int


def _token(name: str) -> str:
    base = name.split("(")[0].split("[")[0]
    return "".join(ch for ch in base.strip().lower() if ch.isalnum())


def _field(name: str) -> Optional[str]:
    token = _token(name)
    if token in _FIELD_ALIASES:
        return _FIELD_ALIASES[token]
    return token if token in _FIELD_ALIASES.values() else None


def _number(raw: str, field: str, where: str) -> float:
    text = raw.strip().replace(",", "")
    try:
        value = float(text)
    except ValueError as exc:
        raise DockingError(f"{where}: field '{field}' is not a number: {raw!r}") from exc
    if not math.isfinite(value):
        raise DockingError(f"{where}: field '{field}' must be finite, got {raw!r}")
    if field in _NON_NEGATIVE_FIELDS and value < 0:
        raise DockingError(f"{where}: field '{field}' must not be negative, got {raw!r}")
    return value


def _optional(values: Mapping[str, str], field: str, where: str) -> Optional[float]:
    raw = values.get(field, "").strip()
    return None if not raw else _number(raw, field, where)


def _build(values: Mapping[str, str], where: str) -> DockingScore:
    missing = [field for field in _REQUIRED_FIELDS if not values.get(field, "").strip()]
    if missing:
        raise DockingError(f"{where}: missing required field(s): {', '.join(missing)}")
    scores = {field: _optional(values, field, where) for field in SCORE_FIELDS}
    if all(value is None for value in scores.values()):
        raise DockingError(
            f"{where}: needs at least one of {', '.join(SCORE_FIELDS)}"
        )
    return DockingScore(
        compound_id=values["compound_id"].strip(),
        target=values["target"].strip(),
        rmsd=_optional(values, "rmsd", where),
        hbond_count=_optional(values, "hbond_count", where),
        key_interactions=values.get("key_interactions", "").strip(),
        **scores,
    )


def parse_docking_scores(rows: Iterable[Mapping[str, str]], source: str = "input") -> List[DockingScore]:
    """Turn raw CSV rows into docking records; refuse duplicates and unusable rows."""
    scores: List[DockingScore] = []
    seen: Dict[Tuple[str, str], int] = {}
    for row_number, raw in enumerate(rows, start=2):
        values: Dict[str, str] = {}
        for key, value in raw.items():
            if key is None:
                continue
            field = _field(key)
            if field is not None:
                values[field] = "" if value is None else str(value).strip()
        score = _build(values, f"{source} row {row_number}")
        key = (score.compound_id, score.target)
        if key in seen:
            raise DockingError(
                f"{source} row {row_number}: duplicate docking result for "
                f"'{score.compound_id}' / '{score.target}' (first seen on row {seen[key]})"
            )
        seen[key] = row_number
        scores.append(score)
    if not scores:
        raise DockingError(f"{source}: no docking rows found")
    return scores


def load_docking_scores(path: Path) -> List[DockingScore]:
    """Read a UTF-8 docking CSV."""
    path = Path(path)
    if not path.is_file():
        raise DockingError(f"docking file not found: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise DockingError(f"{path}: file has no header row")
            rows = list(reader)
    except UnicodeDecodeError as exc:
        raise DockingError(f"{path}: file is not valid UTF-8") from exc
    except csv.Error as exc:
        raise DockingError(f"{path}: malformed CSV: {exc}") from exc
    except OSError as exc:
        raise DockingError(f"{path}: cannot read file: {exc}") from exc
    return parse_docking_scores(rows, source=str(path))


def load_screen_table(path: Path) -> Dict[str, Dict[str, str]]:
    """Read the screening ranked CSV into a compound_id keyed table."""
    path = Path(path)
    if not path.is_file():
        raise DockingError(f"screening file not found: {path}")
    table: Dict[str, Dict[str, str]] = {}
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise DockingError(f"{path}: file has no header row")
            for row_number, raw in enumerate(reader, start=2):
                values = {_field(key): str(value or "").strip() for key, value in raw.items() if key and _field(key)}
                compound_id = values.get("compound_id", "")
                if not compound_id:
                    raise DockingError(f"{path} row {row_number}: field 'compound_id' must not be empty")
                if compound_id in table:
                    raise DockingError(f"{path} row {row_number}: duplicate compound_id '{compound_id}'")
                tier = values.get("tier", "")
                if tier and tier not in TIER_ORDER:
                    raise DockingError(f"{path} row {row_number}: unknown tier '{tier}'")
                table[compound_id] = values
    except UnicodeDecodeError as exc:
        raise DockingError(f"{path}: file is not valid UTF-8") from exc
    except csv.Error as exc:
        raise DockingError(f"{path}: malformed CSV: {exc}") from exc
    except OSError as exc:
        raise DockingError(f"{path}: cannot read file: {exc}") from exc
    if not table:
        raise DockingError(f"{path}: no screening rows found")
    return table


def filter_scores(
    scores: Iterable[DockingScore], max_rmsd: Optional[float] = None, min_vina: Optional[float] = None
) -> Tuple[List[DockingScore], int]:
    """Drop rows without Vina or outside the optional RMSD/Vina limits."""
    kept: List[DockingScore] = []
    dropped = 0
    for score in scores:
        if score.vina_affinity is None:
            dropped += 1
            continue
        if max_rmsd is not None and score.rmsd is not None and score.rmsd > max_rmsd:
            dropped += 1
            continue
        if min_vina is not None and score.vina_affinity > min_vina:
            dropped += 1
            continue
        kept.append(score)
    return kept, dropped


def normalize_scores(scores: Iterable[DockingScore]) -> List[DockingRank]:
    """Scale each scoring method to 0..1 and combine available methods by weight."""
    rows = list(scores)
    bounds: Dict[str, Tuple[float, float]] = {}
    for field in SCORE_FIELDS:
        values = [getattr(score, field) for score in rows if getattr(score, field) is not None]
        if values:
            bounds[field] = (min(values), max(values))
    ranks: List[DockingRank] = []
    for score in rows:
        available = [field for field in SCORE_FIELDS if getattr(score, field) is not None]
        if not available:
            continue
        quality = 0.0
        total_weight = sum(WEIGHTS[field] for field in available)
        for field in available:
            value = getattr(score, field)
            low, high = bounds[field]
            normalized = 1.0 if high == low else (high - value) / (high - low)
            quality += normalized * WEIGHTS[field]
        ranks.append(DockingRank(score=score, consensus=quality / total_weight, coverage=total_weight))
    return ranks


def aggregate_best(ranks: Iterable[DockingRank]) -> Dict[str, CompoundRank]:
    """Keep each compound's highest-consensus target and count its docked targets."""
    grouped: Dict[str, List[DockingRank]] = {}
    for rank in ranks:
        if rank.score.vina_affinity is None:
            continue
        grouped.setdefault(rank.score.compound_id, []).append(rank)
    best: Dict[str, CompoundRank] = {}
    for compound_id, candidates in grouped.items():
        winner = max(
            candidates,
            key=lambda item: (item.consensus, item.coverage, -(item.score.vina_affinity or 0.0)),
        )
        best[compound_id] = CompoundRank(
            compound_id=compound_id, best=winner, targets_docked=len(candidates)
        )
    return best
