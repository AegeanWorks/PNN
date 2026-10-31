"""Command-line end-to-end Phase 2 lead selection for the PNN project."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .compounds import CompoundError
from .curation import CurationError
from .docking import DockingError
from .selection import PipelineError, run_pipeline


def _positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return value


def _bin_width(text: str) -> float:
    value = float(text)
    if not 0.0 < value <= 1.0:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.pipeline",
        description=(
            "Execute the end-to-end PNN Phase 2 lead selection run: "
            "curate compound libraries, screen CNS/BBB properties, join docking consensus scores, "
            "and write hashed manifest and report artifacts."
        ),
    )
    parser.add_argument("inputs", nargs="+", help="one or more compound library CSV paths")
    parser.add_argument("--docking", required=True, help="path to multi-target docking scores CSV")
    parser.add_argument(
        "--out-dir",
        default="phase2_run",
        help="directory to write curated library, screened set, shortlist, report and manifest",
    )
    parser.add_argument(
        "--min-tier",
        choices=("A", "B", "C"),
        default="B",
        help="minimum screening tier floor for curation and shortlist (default: B)",
    )
    parser.add_argument(
        "--per-bin",
        type=_positive_int,
        default=25,
        help="maximum compounds per property-space bin (default: 25)",
    )
    parser.add_argument(
        "--bin-width",
        type=_bin_width,
        default=0.1,
        help="normalized property bin width (default: 0.1)",
    )
    parser.add_argument(
        "--max-compounds",
        type=_positive_int,
        help="optional global cap on curated compounds",
    )
    parser.add_argument(
        "--max-rmsd",
        type=float,
        help="maximum ligand pose RMSD in Angstroms",
    )
    parser.add_argument(
        "--min-vina",
        type=float,
        help="minimum (worst) Vina affinity in kcal/mol",
    )
    parser.add_argument(
        "--top-n-report",
        type=_positive_int,
        default=25,
        help="maximum lead rows to list in lead_report.md (default: 25)",
    )
    return parser


def _failure(what: str, why: str, fix: str) -> int:
    print(f"What: {what}", file=sys.stderr)
    print(f"Why: {why}", file=sys.stderr)
    print(f"Fix: {fix}", file=sys.stderr)
    return 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        run = run_pipeline(
            inputs=args.inputs,
            docking=args.docking,
            out_dir=args.out_dir,
            min_tier=args.min_tier,
            per_bin=args.per_bin,
            bin_width=args.bin_width,
            max_compounds=args.max_compounds,
            max_rmsd=args.max_rmsd,
            min_vina=args.min_vina,
            top_n_report=args.top_n_report,
        )
    except (CompoundError, CurationError, DockingError, PipelineError) as exc:
        return _failure(
            str(exc),
            "lead selection requires valid compound CSVs, compatible docking rows with Vina scores, "
            "and surviving compounds after curation, screening and docking consensus filters.",
            "check inputs against data/compounds/screening_template.csv and "
            "data/docking/scores_template.csv, or adjust --min-tier, --max-rmsd, or --min-vina.",
        )
    except OSError as exc:
        return _failure(
            str(exc),
            "the pipeline artifacts could not be written to disk.",
            "ensure the --out-dir directory path is writable and rerun.",
        )

    shortlist_count = len(run.selection.shortlist)
    print(
        f"Phase 2 lead selection complete: {shortlist_count} lead(s) shortlisted "
        f"from {run.selection.curation.get('input_rows', 0)} input compound row(s) -> {args.out_dir}"
    )
    print("Artifacts generated:")
    for name, path in run.outputs.items():
        print(f"  • {name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
