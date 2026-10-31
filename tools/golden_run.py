#!/usr/bin/env python3
"""Execute the PNN Phase 2 golden reference run and verify artifact integrity."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence

# Anchor to PNN root
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pnn.assay import (
    load_cytotoxicity_data,
    load_degradation_data,
    load_permeability_data,
    qualify_leads,
)
from pnn.selection import SCHEMA, run_pipeline
from pnn.target import list_targets


def _check_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def execute_golden_run(out_dir: Optional[Path] = None) -> int:
    """Run Phase 2 pipeline & Phase 3 assay checks against reference benchmarks."""
    vendor_csv = PROJECT_ROOT / "data" / "compounds" / "vendor_library.csv"
    repurposing_csv = PROJECT_ROOT / "data" / "compounds" / "repurposing_library.csv"
    docking_csv = PROJECT_ROOT / "data" / "docking" / "multi_target_scores.csv"

    deg_csv = PROJECT_ROOT / "data" / "assays" / "degradation_template.csv"
    perm_csv = PROJECT_ROOT / "data" / "assays" / "permeability_template.csv"
    cyto_csv = PROJECT_ROOT / "data" / "assays" / "cytotoxicity_template.csv"

    for required in (vendor_csv, repurposing_csv, docking_csv, deg_csv, perm_csv, cyto_csv):
        if not required.is_file():
            print(f"Error: missing reference dataset: {required}", file=sys.stderr)
            return 1

    temp_dir_obj = None
    if out_dir is None:
        temp_dir_obj = tempfile.TemporaryDirectory()
        dest_dir = Path(temp_dir_obj.name)
    else:
        dest_dir = Path(out_dir)

    try:
        # 1. Target Profile verification
        targets = list_targets()
        if len(targets) < 5:
            print("Error: insufficient registered targets.", file=sys.stderr)
            return 1

        # 2. Phase 2 Lead Selection Pipeline
        run = run_pipeline(
            inputs=[str(vendor_csv), str(repurposing_csv)],
            docking=str(docking_csv),
            out_dir=str(dest_dir),
            min_tier="B",
            per_bin=25,
            max_rmsd=2.0,
        )

        manifest_path = dest_dir / "manifest.json"
        if not manifest_path.is_file():
            print("Error: manifest.json was not generated.", file=sys.stderr)
            return 1

        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("schema") != SCHEMA:
            print(f"Error: manifest schema mismatch: {manifest.get('schema')}", file=sys.stderr)
            return 1

        for name, entry in manifest.get("outputs", {}).items():
            artifact_path = Path(entry["path"])
            if not artifact_path.is_file():
                print(f"Error: missing artifact {name} at {artifact_path}", file=sys.stderr)
                return 1
            calculated_sha = _check_digest(artifact_path)
            if calculated_sha != entry["sha256"]:
                print(f"Error: SHA-256 mismatch for {name}: {calculated_sha} != {entry['sha256']}", file=sys.stderr)
                return 1

        # 3. Phase 3 In Vitro Assay Qualification
        deg = load_degradation_data(deg_csv)
        perm = load_permeability_data(perm_csv)
        cyto = load_cytotoxicity_data(cyto_csv)
        hits = qualify_leads(deg, perm, cyto)
        confirmed_hits = [h for h in hits if h.status == "CONFIRMED_HIT"]

        print("PNN Golden Run Verification PASSED:")
        print(f"  • Registered Targets:  {len(targets)}")
        print(f"  • Curated Candidates: {len(run.selection.curated)}")
        print(f"  • Screened Compounds: {len(run.selection.screened)}")
        print(f"  • Shortlisted Leads:  {len(run.selection.shortlist)}")
        print(f"  • Verified Artifacts: {len(manifest.get('outputs', {}))} files (SHA-256 intact)")
        print(f"  • Confirmed In Vitro: {len(confirmed_hits)} hit(s)")
        print("\nTop Shortlisted Leads:")
        for row in run.selection.shortlist[:5]:
            print(f"  [{row.rank}] {row.compound_id:<12} {row.name:<24} Target: {row.target:<10} Consensus: {row.consensus:.4f} (Vina: {row.vina_affinity:.1f})")
        print("\nTop Confirmed In Vitro Hits:")
        for h in confirmed_hits[:3]:
            print(f"  • {h.compound_id:<12} EC50: {h.ec50_nm:.1f} nM (Max: {h.max_degradation_pct:.1f}%) | {h.bbb_status} | SI: {h.selectivity_index:.1f}")
    finally:
        if temp_dir_obj is not None:
            temp_dir_obj.cleanup()

    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Execute and verify the PNN Phase 2/3 golden reference run.")
    parser.add_argument("--out-dir", help="optional output directory (defaults to ephemeral temp dir)")
    args = parser.parse_args(argv)
    return execute_golden_run(out_dir=Path(args.out_dir) if args.out_dir else None)


if __name__ == "__main__":
    raise SystemExit(main())
