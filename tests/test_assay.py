"""Unit and CLI tests for Phase 3 in vitro assay ingestion and hit qualification."""

from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from pnn.assay import (
    AssayError,
    CytotoxicityResult,
    DegradationResult,
    PermeabilityResult,
    load_cytotoxicity_data,
    load_degradation_data,
    load_permeability_data,
    main,
    qualify_leads,
    render_assay_report,
)

DEG_HEADER = "compound_id,ec50_nm,max_degradation_pct,hill_slope,test_concentration_um,substrate\n"
DEG_ROWS = (
    "CMP-1,25.0,90.0,1.1,10.0,DQ-CS\n"
    "CMP-2,150.0,75.0,1.0,10.0,DQ-CS\n"
    "CMP-3,800.0,40.0,0.9,10.0,DQ-CS\n"
    "CMP-4,65000.0,5.0,0.5,10.0,DQ-CS\n"
)

PERM_HEADER = "compound_id,assay_type,pe_1e6_cm_s,logbb,papp_ab_1e6_cm_s,efflux_ratio,pgp_substrate\n"
PERM_ROWS = (
    "CMP-1,PAMPA,5.5,0.1,,,false\n"
    "CMP-2,MDCK,,,8.0,1.5,false\n"
    "CMP-3,MDCK,,,2.0,5.0,true\n"
    "CMP-4,PAMPA,0.2,-2.0,,,true\n"
)

CYTO_HEADER = "compound_id,cc50_um,cell_line,exposure_hours\n"
CYTO_ROWS = (
    "CMP-1,80.0,mHypoA-60/3,24\n"
    "CMP-2,60.0,mHypoA-60/3,24\n"
    "CMP-3,15.0,mHypoA-60/3,24\n"
    "CMP-4,5.0,mHypoA-60/3,24\n"
)


class AssayTest(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

        self.deg_file = self.root / "deg.csv"
        self.deg_file.write_text(DEG_HEADER + DEG_ROWS, encoding="utf-8")

        self.perm_file = self.root / "perm.csv"
        self.perm_file.write_text(PERM_HEADER + PERM_ROWS, encoding="utf-8")

        self.cyto_file = self.root / "cyto.csv"
        self.cyto_file.write_text(CYTO_HEADER + CYTO_ROWS, encoding="utf-8")

    def run_cli(self, *argv: str):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            code = main(list(argv))
        return code, stdout.getvalue(), stderr.getvalue()

    def test_load_degradation_data(self) -> None:
        data = load_degradation_data(self.deg_file)
        self.assertEqual(len(data), 4)
        self.assertEqual(data["CMP-1"].ec50_nm, 25.0)
        self.assertEqual(data["CMP-1"].max_degradation_pct, 90.0)
        self.assertEqual(data["CMP-1"].substrate, "DQ-CS")

    def test_load_degradation_invalid_numeric_raises(self) -> None:
        broken = self.root / "broken_deg.csv"
        broken.write_text(DEG_HEADER + "CMP-ERR,not_a_number,50.0,1.0,10.0,DQ-CS\n", encoding="utf-8")
        with self.assertRaises(AssayError) as ctx:
            load_degradation_data(broken)
        self.assertIn("invalid numeric value", str(ctx.exception))

        broken_inf = self.root / "broken_inf.csv"
        broken_inf.write_text(DEG_HEADER + "CMP-INF,inf,50.0,1.0,10.0,DQ-CS\n", encoding="utf-8")
        with self.assertRaises(AssayError) as ctx:
            load_degradation_data(broken_inf)
        self.assertIn("non-finite number", str(ctx.exception))

    def test_load_permeability_and_cytotoxicity(self) -> None:
        perm = load_permeability_data(self.perm_file)
        self.assertEqual(len(perm), 4)
        self.assertEqual(perm["CMP-1"].pe_1e6_cm_s, 5.5)
        self.assertEqual(perm["CMP-2"].efflux_ratio, 1.5)

        cyto = load_cytotoxicity_data(self.cyto_file)
        self.assertEqual(len(cyto), 4)
        self.assertEqual(cyto["CMP-1"].cc50_um, 80.0)

    def test_qualify_leads_classifications(self) -> None:
        deg = load_degradation_data(self.deg_file)
        perm = load_permeability_data(self.perm_file)
        cyto = load_cytotoxicity_data(self.cyto_file)

        hits = qualify_leads(deg, perm, cyto)
        status_map = {h.compound_id: h.status for h in hits}

        # CMP-1: potent (25nM), PAMPA 5.5, CC50 80uM -> CONFIRMED_HIT
        self.assertEqual(status_map["CMP-1"], "CONFIRMED_HIT")

        # CMP-2: potent (150nM), MDCK efflux 1.5, CC50 60uM -> CONFIRMED_HIT
        self.assertEqual(status_map["CMP-2"], "CONFIRMED_HIT")

        # CMP-3: CC50 is 15uM (<20uM) -> CYTOTOXIC
        self.assertEqual(status_map["CMP-3"], "CYTOTOXIC")

        # CMP-4: EC50 65uM, CC50 5uM -> CYTOTOXIC
        self.assertEqual(status_map["CMP-4"], "CYTOTOXIC")

    def test_qualify_leads_untested_safety_yields_qualified_lead(self) -> None:
        deg = {"LEAD-A": DegradationResult("LEAD-A", ec50_nm=50.0, max_degradation_pct=85.0)}
        perm = {"LEAD-A": PermeabilityResult("LEAD-A", assay_type="PAMPA", pe_1e6_cm_s=6.0)}
        hits = qualify_leads(deg, perm, cytotoxicity=None)
        self.assertEqual(hits[0].status, "QUALIFIED_LEAD")
        self.assertIn("toxicity status: UNTESTED", hits[0].findings)

    def test_render_assay_report_contains_guideline_and_tables(self) -> None:
        deg = load_degradation_data(self.deg_file)
        perm = load_permeability_data(self.perm_file)
        cyto = load_cytotoxicity_data(self.cyto_file)
        hits = qualify_leads(deg, perm, cyto)

        report = render_assay_report(hits, {"ec50_ceiling_nm": 10000.0, "pampa_pe_floor": 3.0})
        self.assertIn("G-001", report)
        self.assertIn("In vitro assay qualification confirms target degradation activity", report)
        self.assertIn("CMP-1", report)
        self.assertIn("CONFIRMED_HIT", report)
        self.assertIn("Phase 4 In Vivo Proof-of-Concept", report)

    def test_cli_execution_success(self) -> None:
        out_csv = self.root / "hits.csv"
        json_out = self.root / "summary.json"
        report_out = self.root / "report.md"

        code, stdout, stderr = self.run_cli(
            "--degradation",
            str(self.deg_file),
            "--permeability",
            str(self.perm_file),
            "--cytotoxicity",
            str(self.cyto_file),
            "--out",
            str(out_csv),
            "--json",
            str(json_out),
            "--report",
            str(report_out),
        )
        self.assertEqual(code, 0)
        self.assertIn("In vitro hit qualification complete: 2 confirmed hit(s)", stdout)
        self.assertTrue(out_csv.is_file())
        self.assertTrue(json_out.is_file())
        self.assertTrue(report_out.is_file())

        with out_csv.open("r", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]["compound_id"], "CMP-1")
        self.assertEqual(rows[0]["status"], "CONFIRMED_HIT")

    def test_cli_failure_triad_on_missing_file(self) -> None:
        code, stdout, stderr = self.run_cli(
            "--degradation",
            str(self.root / "nonexistent.csv"),
        )
        self.assertEqual(code, 1)
        self.assertIn("What:", stderr)
        self.assertIn("Why:", stderr)
        self.assertIn("Fix:", stderr)


if __name__ == "__main__":
    unittest.main()
