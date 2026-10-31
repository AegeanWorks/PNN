"""Command-line consensus ranking of docking results for PNN Phase 2."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .docking import (
    TIER_ORDER,
    WEIGHTS,
    CompoundRank,
    DockingError,
    aggregate_best,
    filter_scores,
    load_docking_scores,
    load_screen_table,
    normalize_scores,
)
from .reports import write_csv, write_json

RANK_FIELDS = (
    "rank",
    "compound_id",
    "name",
    "smiles",
    "tier",
    "property_score",
    "target",
    "consensus_score",
    "coverage",
    "vina_affinity",
    "glide_score",
    "mmgbsa_dg",
    "rmsd",
    "hbond_count",
    "key_interactions",
    "targets_docked",
)
_TIER_CEILING = {"A": 0, "B": 1, "C": 2}


@dataclass(frozen=True)
class ShortlistRow:
    """One ranked compound ready for CSV output."""

    rank: int
    compound_id: str
    name: str
    smiles: str
    tier: str
    property_score: str
    target: str
    consensus: float
    coverage: float
    vina_affinity: float
    glide_score: Optional[float]
    mmgbsa_dg: Optional[float]
    rmsd: Optional[float]
    hbond_count: Optional[float]
    key_interactions: str
    targets_docked: int

    def as_row(self) -> Dict[str, object]:
        return {
            "rank": self.rank,
            "compound_id": self.compound_id,
            "name": self.name,
            "smiles": self.smiles,
            "tier": self.tier,
            "property_score": self.property_score,
            "target": self.target,
            "consensus_score": f"{self.consensus:.4f}",
            "coverage": f"{self.coverage:.2f}",
            "vina_affinity": self.vina_affinity,
            "glide_score": "" if self.glide_score is None else self.glide_score,
            "mmgbsa_dg": "" if self.mmgbsa_dg is None else self.mmgbsa_dg,
            "rmsd": "" if self.rmsd is None else self.rmsd,
            "hbond_count": "" if self.hbond_count is None else self.hbond_count,
            "key_interactions": self.key_interactions,
            "targets_docked": self.targets_docked,
        }


def _property_score(info: Optional[Dict[str, str]]) -> float:
    raw = (info or {}).get("score", "")
    try:
        return float(raw)
    except ValueError:
        return -1.0


def _selected(
    best: Dict[str, CompoundRank],
    screen: Optional[Dict[str, Dict[str, str]]],
    min_tier: str,
) -> tuple:
    """Filter aggregated compounds by screening evidence and sort them by consensus."""
    selected: List[tuple] = []
    skipped_tier = 0
    skipped_unlisted = 0
    for compound_id, rank in best.items():
        info = screen.get(compound_id) if screen is not None else None
        if screen is not None:
            if info is None:
                skipped_unlisted += 1
                continue
            tier = info.get("tier", "")
            if tier not in _TIER_CEILING or _TIER_CEILING[tier] > _TIER_CEILING[min_tier]:
                skipped_tier += 1
                continue
        selected.append((rank, info))
    selected.sort(
        key=lambda item: (
            -item[0].best.consensus,
            -item[0].best.coverage,
            -_property_score(item[1]),
            item[0].compound_id,
        )
    )
    return selected, skipped_tier, skipped_unlisted


def build_shortlist(
    scores: Sequence,
    screen: Optional[Dict[str, Dict[str, str]]] = None,
    min_tier: str = "B",
    max_rmsd: Optional[float] = None,
    min_vina: Optional[float] = None,
) -> tuple:
    """Return shortlist rows and a summary of what each filter removed."""
    if screen is not None and min_tier not in _TIER_CEILING:
        raise DockingError(f"unknown screening tier: {min_tier}")
    kept, dropped = filter_scores(scores, max_rmsd=max_rmsd, min_vina=min_vina)
    best = aggregate_best(normalize_scores(kept))
    selected, skipped_tier, skipped_unlisted = _selected(best, screen, min_tier)
    rows = [
        ShortlistRow(
            rank=index,
            compound_id=rank.compound_id,
            name=(info or {}).get("name", ""),
            smiles=(info or {}).get("smiles", ""),
            tier=(info or {}).get("tier", ""),
            property_score=(info or {}).get("score", ""),
            target=rank.best.score.target,
            consensus=rank.best.consensus,
            coverage=rank.best.coverage,
            vina_affinity=rank.best.score.vina_affinity or 0.0,
            glide_score=rank.best.score.glide_score,
            mmgbsa_dg=rank.best.score.mmgbsa_dg,
            rmsd=rank.best.score.rmsd,
            hbond_count=rank.best.score.hbond_count,
            key_interactions=rank.best.score.key_interactions,
            targets_docked=rank.targets_docked,
        )
        for index, (rank, info) in enumerate(selected, start=1)
    ]
    summary = {
        "docking_rows": len(scores),
        "rows_after_filters": len(kept),
        "skipped_by_filters": dropped,
        "compounds_ranked": len(best),
        "skipped_by_tier": skipped_tier,
        "skipped_not_in_screen": skipped_unlisted,
        "emitted": len(rows),
        "weights": dict(WEIGHTS),
        "filters": {"min_tier": min_tier, "max_rmsd": max_rmsd, "min_vina": min_vina},
    }
    return rows, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.rank",
        description=(
            "Rank docking results with normalized Vina/Glide/MM-GBSA consensus and an optional "
            "screening-tier filter."
        ),
    )
    parser.add_argument("docking", help="CSV with compound_id, target and at least one docking score")
    parser.add_argument("--screen", help="optional ranked CSV from `python -m pnn` for tier/name/SMILES")
    parser.add_argument("--out", default="shortlist.csv", help="shortlist CSV path")
    parser.add_argument("--json", dest="json_path", help="optional JSON summary path")
    parser.add_argument(
        "--min-tier",
        choices=("A", "B", "C"),
        default="B",
        help="keep this screening tier and better (default: B; ignored without --screen)",
    )
    parser.add_argument("--max-rmsd", type=float, help="drop target rows with ligand RMSD above this value")
    parser.add_argument("--min-vina", type=float, help="drop target rows with Vina affinity above this cutoff")
    return parser


def _failure(what: str, why: str, fix: str) -> int:
    print(f"What: {what}", file=sys.stderr)
    print(f"Why: {why}", file=sys.stderr)
    print(f"Fix: {fix}", file=sys.stderr)
    return 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    out = Path(args.out)
    try:
        scores = load_docking_scores(Path(args.docking))
        screen = load_screen_table(Path(args.screen)) if args.screen else None
        rows, summary = build_shortlist(
            scores,
            screen=screen,
            min_tier=args.min_tier,
            max_rmsd=args.max_rmsd,
            min_vina=args.min_vina,
        )
        if not rows:
            return _failure(
                "no compounds survived the docking and screening filters.",
                "every docking row was missing Vina, outside the RMSD/Vina limits, or below the tier floor.",
                "relax --max-rmsd, --min-vina or --min-tier, or supply docking rows with Vina scores.",
            )
        write_csv(out, RANK_FIELDS, [row.as_row() for row in rows])
        if args.json_path:
            write_json(Path(args.json_path), {"output": str(out), **summary})
    except DockingError as exc:
        return _failure(
            str(exc),
            "the ranker needs UTF-8 docking (and optional screening) CSVs matching the documented contract.",
            "check the contracts in computational/docking_ranking.md, then rerun with valid files.",
        )
    except OSError as exc:
        return _failure(
            str(exc),
            "the shortlist CSV or JSON summary could not be written.",
            "choose a writable --out/--json path and rerun.",
        )
    print(
        f"Aggregated {summary['compounds_ranked']} compound(s) from {summary['docking_rows']} docking row(s): "
        f"emitted {summary['emitted']} -> {out}"
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())