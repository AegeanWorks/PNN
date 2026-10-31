"""Unit tests for the IND-enabling safety pharmacology and toxicology package."""

import json
import tempfile
import unittest
from pathlib import Path

from pnn.safety import (
    SAFETY_FIELDS,
    SafetyError,
    SafetyPharm,
    ToxStudy,
    assess_safety,
    load_safety_pharmacology,
    load_tox_studies,
    main,
    most_sensitive_study,
    noael_hed_mg_kg,
    render_safety_report,
)

TOX_HEADER = "compound_id,species,study_type,duration_days,noael_mg_kg,loael_mg_kg,target_organ,findings\n"
PHARM_HEADER = (
    "compound_id,herg_ic50_um,free_cmax_um,cyp3a4_inhibition_pct,cyp2d6_inhibition_pct,"
    "ames_positive,micronucleus_positive,respiratory_noael_mg_kg,cns_noael_mg_kg\n"
)


class TestSafetyIngestion(unittest.TestCase):
    def test_missing_tox_file_raises(self) -> None:
        with self.assertRaises(SafetyError):
            load_tox_studies(Path("does_not_exist_tox.csv"))

    def test_missing_tox_columns_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "tox.csv"
            path.write_text("compound_id,species\nC1,rat\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_tox_studies(path)

    def test_unknown_species_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "tox.csv"
            path.write_text(TOX_HEADER + "C1,zebrafish,REPEAT_DOSE,28,10.0,,none,\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_tox_studies(path)

    def test_non_positive_noael_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "tox.csv"
            path.write_text(TOX_HEADER + "C1,rat,REPEAT_DOSE,28,0,,none,\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_tox_studies(path)

    def test_loael_below_noael_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "tox.csv"
            path.write_text(TOX_HEADER + "C1,rat,REPEAT_DOSE,28,50.0,25.0,liver,\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_tox_studies(path)

    def test_invalid_duration_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "tox.csv"
            path.write_text(TOX_HEADER + "C1,rat,REPEAT_DOSE,1.5,50.0,,none,\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_tox_studies(path)

    def test_load_valid_tox_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "tox.csv"
            path.write_text(
                TOX_HEADER
                + "C1,rat,REPEAT_DOSE_28D,28,50.0,100.0,liver,Reversible\n"
                + "C1,dog,REPEAT_DOSE_28D,28,30.0,60.0,none,None at NOAEL\n",
                encoding="utf-8",
            )
            studies = load_tox_studies(path)
            self.assertEqual(len(studies), 2)
            self.assertEqual(studies[0].species, "rat")
            self.assertEqual(studies[0].noael_mg_kg, 50.0)

    def test_pharm_duplicate_compound_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "pharm.csv"
            path.write_text(
                PHARM_HEADER
                + "C1,30,0.5,10,5,false,false,50,40\n"
                + "C1,30,0.5,10,5,false,false,50,40\n",
                encoding="utf-8",
            )
            with self.assertRaises(SafetyError):
                load_safety_pharmacology(path)

    def test_pharm_invalid_boolean_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "pharm.csv"
            path.write_text(PHARM_HEADER + "C1,30,0.5,10,5,maybe,false,50,40\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_safety_pharmacology(path)

    def test_pharm_percent_out_of_range_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "pharm.csv"
            path.write_text(PHARM_HEADER + "C1,30,0.5,101,5,false,false,50,40\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_safety_pharmacology(path)

    def test_pharm_missing_compound_id_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "pharm.csv"
            path.write_text("herg_ic50_um,free_cmax_um\n30,0.5\n", encoding="utf-8")
            with self.assertRaises(SafetyError):
                load_safety_pharmacology(path)


class TestSafetyAssessment(unittest.TestCase):
    def _safe_pharm(self, compound_id: str = "C1") -> SafetyPharm:
        return SafetyPharm(
            compound_id=compound_id,
            herg_ic50_um=30.0,
            free_cmax_um=0.5,
            cyp3a4_inhibition_pct=10.0,
            cyp2d6_inhibition_pct=5.0,
            ames_positive=False,
            micronucleus_positive=False,
            respiratory_noael_mg_kg=50.0,
            cns_noael_mg_kg=40.0,
        )

    def _tox(self, compound_id: str = "C1") -> ToxStudy:
        return ToxStudy(
            compound_id=compound_id,
            species="rat",
            study_type="REPEAT_DOSE_28D",
            duration_days=28,
            noael_mg_kg=50.0,
            loael_mg_kg=100.0,
            target_organ="liver",
            findings="Reversible at LOAEL only",
        )

    def test_safe_compound_projection_and_status(self) -> None:
        assessment = assess_safety([self._tox()], {"C1": self._safe_pharm()})[0]
        self.assertEqual(assessment.status, "SAFE_TO_PROGRESS")
        self.assertAlmostEqual(assessment.noael_hed_mg_kg or 0.0, 50.0 * 6.0 / 37.0, places=6)
        self.assertAlmostEqual(assessment.phase1_mrsd_mg or 0.0, (50.0 * 6.0 / 37.0 / 10.0) * 60.0, places=6)
        self.assertAlmostEqual(assessment.herg_margin or 0.0, 60.0, places=6)
        self.assertEqual(assessment.genotoxicity, "NEGATIVE")

    def test_most_sensitive_species_is_selected_by_hed_not_raw_noael(self) -> None:
        # Rat 50 mg/kg -> HED 8.108; dog 20 mg/kg -> HED 10.811. Raw mg/kg would wrongly pick dog.
        dog = ToxStudy("C1", "dog", "REPEAT_DOSE_28D", 28, 20.0, None, "none", "")
        assessment = assess_safety([self._tox(), dog], {"C1": self._safe_pharm()})[0]
        self.assertEqual(assessment.most_sensitive_species, "rat")
        self.assertAlmostEqual(assessment.noael_mg_kg or 0.0, 50.0, places=6)
        self.assertAlmostEqual(assessment.noael_hed_mg_kg or 0.0, 50.0 * 6.0 / 37.0, places=6)
        self.assertAlmostEqual(assessment.phase1_mrsd_mg or 0.0, (50.0 * 6.0 / 37.0 / 10.0) * 60.0, places=6)

    def test_dog_wins_when_its_hed_is_lower(self) -> None:
        rat = ToxStudy("C1", "rat", "REPEAT_DOSE_28D", 28, 100.0, None, "none", "")
        dog = ToxStudy("C1", "dog", "REPEAT_DOSE_28D", 28, 20.0, None, "none", "")
        assessment = assess_safety([rat, dog], {"C1": self._safe_pharm()})[0]
        self.assertEqual(assessment.most_sensitive_species, "dog")
        self.assertAlmostEqual(assessment.noael_hed_mg_kg or 0.0, 20.0 * 20.0 / 37.0, places=6)

    def test_most_sensitive_study_helper_ranks_on_hed(self) -> None:
        rat = ToxStudy("C1", "rat", "REPEAT_DOSE_28D", 28, 50.0, None, "", "")
        dog = ToxStudy("C1", "dog", "REPEAT_DOSE_28D", 28, 30.0, None, "", "")
        self.assertEqual(most_sensitive_study([rat, dog]), rat)
        self.assertIsNone(most_sensitive_study([]))
        self.assertAlmostEqual(noael_hed_mg_kg(50.0, "rat"), 50.0 * 6.0 / 37.0, places=9)
        with self.assertRaises(SafetyError):
            noael_hed_mg_kg(50.0, "zebrafish")

    def test_missing_tox_study_is_caution(self) -> None:
        assessment = assess_safety([], {"C1": self._safe_pharm()})[0]
        self.assertEqual(assessment.status, "CAUTION")
        self.assertIn("No repeat-dose toxicology study", assessment.findings)

    def test_missing_safety_pharmacology_is_caution(self) -> None:
        assessment = assess_safety([self._tox()], {})[0]
        self.assertEqual(assessment.status, "CAUTION")
        self.assertIn("Genotoxicity data incomplete", assessment.findings)
        self.assertIn("hERG margin not computable", assessment.findings)

    def test_positive_genotoxicity_fails(self) -> None:
        pharm = SafetyPharm(
            compound_id="C1",
            herg_ic50_um=30.0,
            free_cmax_um=0.5,
            cyp3a4_inhibition_pct=10.0,
            ames_positive=True,
            micronucleus_positive=False,
        )
        assessment = assess_safety([self._tox()], {"C1": pharm})[0]
        self.assertEqual(assessment.status, "FAIL_SAFETY")
        self.assertEqual(assessment.genotoxicity, "POSITIVE")

    def test_narrow_herg_margin_fails(self) -> None:
        pharm = SafetyPharm(
            compound_id="C1",
            herg_ic50_um=2.0,
            free_cmax_um=1.0,
            cyp3a4_inhibition_pct=10.0,
            ames_positive=False,
            micronucleus_positive=False,
        )
        assessment = assess_safety([self._tox()], {"C1": pharm})[0]
        self.assertEqual(assessment.status, "FAIL_SAFETY")
        self.assertIn("Narrow hERG margin", assessment.findings)

    def test_cyp3a4_ceiling_fails(self) -> None:
        pharm = SafetyPharm(
            compound_id="C1",
            herg_ic50_um=30.0,
            free_cmax_um=0.5,
            cyp3a4_inhibition_pct=60.0,
            ames_positive=False,
            micronucleus_positive=False,
        )
        assessment = assess_safety([self._tox()], {"C1": pharm})[0]
        self.assertEqual(assessment.status, "FAIL_SAFETY")

    def test_marginal_herg_and_cyp_are_caution(self) -> None:
        pharm = SafetyPharm(
            compound_id="C1",
            herg_ic50_um=10.0,
            free_cmax_um=0.5,
            cyp3a4_inhibition_pct=35.0,
            ames_positive=False,
            micronucleus_positive=False,
            respiratory_noael_mg_kg=50.0,
            cns_noael_mg_kg=40.0,
        )
        assessment = assess_safety([self._tox()], {"C1": pharm})[0]
        self.assertEqual(assessment.status, "CAUTION")
        self.assertIn("hERG margin below target", assessment.findings)
        self.assertIn("approaching", assessment.findings)

    def test_invalid_thresholds_raise(self) -> None:
        with self.assertRaises(SafetyError):
            assess_safety([self._tox()], {}, safety_factor=0.0)
        with self.assertRaises(SafetyError):
            assess_safety([self._tox()], {}, herg_margin_min=5.0, herg_margin_fail=10.0)
        with self.assertRaises(SafetyError):
            assess_safety([self._tox()], {}, cyp_inhibition_caution_pct=60.0, cyp_inhibition_max_pct=50.0)

    def test_report_includes_matrix_and_limits(self) -> None:
        assessments = assess_safety([self._tox()], {"C1": self._safe_pharm()})
        report = render_safety_report(assessments, {"safety_factor": 10.0})
        self.assertIn("IND-Enabling Safety Pharmacology", report)
        self.assertIn("G-001", report)
        self.assertIn("Candidate Safety Matrix", report)
        self.assertIn("Package Limits", report)


class TestSafetyCli(unittest.TestCase):
    def test_cli_requires_an_input(self) -> None:
        self.assertEqual(main([]), 1)

    def test_cli_missing_file_returns_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out = Path(tmp_dir) / "safety.csv"
            self.assertEqual(main(["--tox", str(Path(tmp_dir) / "missing.csv"), "--out", str(out)]), 1)

    def test_cli_writes_csv_json_and_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            tox = root / "tox.csv"
            pharm = root / "pharm.csv"
            tox.write_text(TOX_HEADER + "C1,rat,REPEAT_DOSE_28D,28,50.0,100.0,liver,\n", encoding="utf-8")
            pharm.write_text(
                PHARM_HEADER + "C1,30,0.5,10,5,false,false,50,40\n",
                encoding="utf-8",
            )
            out_csv = root / "safety.csv"
            out_json = root / "safety.json"
            out_md = root / "safety.md"

            rc = main(
                [
                    "--tox", str(tox),
                    "--pharm", str(pharm),
                    "--out", str(out_csv),
                    "--json", str(out_json),
                    "--report", str(out_md),
                ]
            )
            self.assertEqual(rc, 0)
            self.assertTrue(out_csv.is_file())
            self.assertTrue(out_md.is_file())
            header = out_csv.read_text(encoding="utf-8").splitlines()[0]
            self.assertEqual(header, ",".join(SAFETY_FIELDS))
            summary = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(summary["total_assessed"], 1)
            self.assertEqual(summary["by_status"]["SAFE_TO_PROGRESS"], 1)


if __name__ == "__main__":
    unittest.main()
