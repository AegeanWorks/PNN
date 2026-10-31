"""Unit and CLI integration tests for the end-to-end Phase 2 lead selection pipeline."""

from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from pnn.pipeline import main
from pnn.selection import (
    LIMITS,
    SCHEMA,
    PipelineError,
    render_lead_report,
    run_pipeline,
    select_leads,
)

HEADER = "compound_id,mw,logp,tpsa,hbd,hba,rotatable_bonds,pka,charge,logbb,brain_plasma_ratio,oral_bioavailability,cyp3a4_inhibition,herg_liability,pnn_ec50_nm\n"
ROW_A = "CMP-A,350,2.5,55,2,5,4,8.0,neutral,-0.2,0.5,0.35,10,5,5\n"
ROW_B = "CMP-B,250,3.0,60,3,6,6,7.5,neutral,-0.7,0.2,0.18,35,15,400\n"
ROW_C = "CMP-C,210,3.4,68,3,7,7,7.2,neutral,-0.9,0.15,0.12,45,25,800\n"
ROW_REJECT = "CMP-REJ,550,5.5,120,6,10,12,12.0,acidic,-2.0,0.01,0.02,80,90,5000\n"

DOCKING_HEADER = "compound_id,target,vina_affinity,glide_score,mmgbsa_dg,rmsd\n"
DOCKING_ROWS = (
    "CMP-A,PTPsigma,-10.5,-8.2,-52.0,1.1\n"
    "CMP-A,LAR,-6.5,-5.0,-22.0,1.4\n"
    "CMP-B,PTPsigma,-8.8,-4.5,-32.0,1.0\n"
    "CMP-C,PTPsigma,-6.0,-3.0,-12.0,1.2\n"
)


class PipelineTest(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

        self.input1 = self.root / "vendor1.csv"
        self.input1.write_text(HEADER + ROW_A + ROW_B, encoding="utf-8")

        self.input2 = self.root / "vendor2.csv"
        self.input2.write_text(HEADER + ROW_B + ROW_C + ROW_REJECT, encoding="utf-8")

        self.docking = self.root / "docking.csv"
        self.docking.write_text(DOCKING_HEADER + DOCKING_ROWS, encoding="utf-8")

    def read_rows(self, path: Path):
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    def run_cli(self, *argv: str):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            code = main(list(argv))
        return code, stdout.getvalue(), stderr.getvalue()

    def test_select_leads_end_to_end(self) -> None:
        selection = select_leads(
            inputs=[str(self.input1), str(self.input2)],
            docking=str(self.docking),
            min_tier="B",
        )
        self.assertEqual(len(selection.shortlist), 3)
        lead_ids = [row.compound_id for row in selection.shortlist]
        self.assertEqual(lead_ids, ["CMP-A", "CMP-B", "CMP-C"])

        # Check provenance digests
        self.assertEqual(len(selection.inputs), 2)
        for item in selection.inputs:
            self.assertIn("sha256", item)
            self.assertEqual(len(item["sha256"]), 64)

        # Curation checks
        self.assertEqual(selection.curation["duplicate_ids"], 1)
        self.assertEqual(selection.curation["skipped_by_tier"], 1)  # CMP-REJ
        self.assertEqual(selection.curation["selected"], 3)  # CMP-A, CMP-B, CMP-C

    def test_select_leads_raises_on_empty_input(self) -> None:
        with self.assertRaises(PipelineError) as ctx:
            select_leads(inputs=[], docking=str(self.docking))
        self.assertIn("no input", str(ctx.exception))

    def test_select_leads_raises_when_curation_empties_set(self) -> None:
        only_reject = self.root / "rejects.csv"
        only_reject.write_text(HEADER + ROW_REJECT, encoding="utf-8")
        with self.assertRaises(PipelineError) as ctx:
            select_leads(inputs=[str(only_reject)], docking=str(self.docking), min_tier="A")
        self.assertIn("no compound survived curation", str(ctx.exception))

    def test_select_leads_raises_when_no_shortlist_survives(self) -> None:
        docking_unrelated = self.root / "docking_unrelated.csv"
        docking_unrelated.write_text(
            DOCKING_HEADER + "OTHER-1,PTPsigma,-9.0,-5.0,-40.0,1.0\n", encoding="utf-8"
        )
        with self.assertRaises(PipelineError) as ctx:
            select_leads(inputs=[str(self.input1)], docking=str(docking_unrelated))
        self.assertIn("no compound survived docking", str(ctx.exception))

    def test_render_lead_report_guidelines_and_evidence(self) -> None:
        selection = select_leads(
            inputs=[str(self.input1), str(self.input2)],
            docking=str(self.docking),
            min_tier="B",
        )
        options = {"min_tier": "B", "per_bin": 25}
        report = render_lead_report(selection, options)

        # Project guidelines G-001, G-002, G-003 compliance
        self.assertIn("G-001", report)
        self.assertIn("G-002", report)
        self.assertIn("G-003", report)
        self.assertIn("Screening tiers (A, B, C) are triage signals only", report)
        self.assertIn("Docking consensus scores are relative to the docked candidate set", report)
        self.assertIn("exact-SMILES deduplication is not chemical canonicalization", report)

        # Input provenance
        self.assertIn(str(self.input1), report)
        self.assertIn(str(self.docking), report)

        # Table rows
        self.assertIn("CMP-A", report)
        self.assertIn("CMP-B", report)

        # Evidence gaps
        self.assertIn("PAMPA-BBB", report)
        self.assertIn("MDCK-MDR1", report)
        self.assertIn("SPR/ITC", report)
        self.assertIn("WFA-positive", report)

    def test_run_pipeline_produces_all_artifacts_and_hashed_manifest(self) -> None:
        out_dir = self.root / "run_out"
        run = run_pipeline(
            inputs=[str(self.input1), str(self.input2)],
            docking=str(self.docking),
            out_dir=str(out_dir),
            min_tier="B",
        )

        # Output files exist
        self.assertTrue(run.outputs["curated"].is_file())
        self.assertTrue(run.outputs["ranked"].is_file())
        self.assertTrue(run.outputs["shortlist"].is_file())
        self.assertTrue(run.outputs["report"].is_file())
        self.assertTrue(run.outputs["manifest"].is_file())

        # Verify manifest schema and hashes
        manifest = json.loads(run.outputs["manifest"].read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema"], SCHEMA)
        self.assertEqual(manifest["limits"], list(LIMITS))

        for name, item in manifest["outputs"].items():
            path = Path(item["path"])
            self.assertTrue(path.is_file(), f"{name} file exists")
            content = path.read_bytes()
            expected_hash = hashlib.sha256(content).hexdigest()
            self.assertEqual(item["sha256"], expected_hash, f"Hash matches for {name}")

    def test_cli_execution_success(self) -> None:
        out_dir = self.root / "cli_run"
        code, stdout, stderr = self.run_cli(
            str(self.input1),
            str(self.input2),
            "--docking",
            str(self.docking),
            "--out-dir",
            str(out_dir),
            "--min-tier",
            "B",
        )
        self.assertEqual(code, 0)
        self.assertIn("Phase 2 lead selection complete: 3 lead(s) shortlisted", stdout)
        self.assertTrue((out_dir / "shortlist.csv").is_file())
        self.assertTrue((out_dir / "lead_report.md").is_file())
        self.assertTrue((out_dir / "manifest.json").is_file())

    def test_cli_execution_failure_triad(self) -> None:
        broken_input = self.root / "broken.csv"
        broken_input.write_text("invalid,header\n1,2\n", encoding="utf-8")
        code, stdout, stderr = self.run_cli(
            str(broken_input),
            "--docking",
            str(self.docking),
        )
        self.assertEqual(code, 1)
        self.assertIn("What:", stderr)
        self.assertIn("Why:", stderr)
        self.assertIn("Fix:", stderr)

    def test_cli_help(self) -> None:
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), self.assertRaises(SystemExit) as ctx:
            main(["--help"])
        self.assertEqual(ctx.exception.code, 0)
        self.assertIn("Execute the end-to-end PNN Phase 2 lead selection run", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
