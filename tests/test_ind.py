"""Unit tests for human allometric dose projection and IND regulatory dossier generator."""

import json
import tempfile
import unittest
from pathlib import Path

from pnn.ind import (
    IndError,
    calculate_allometric_scaling,
    compile_ind_dossier,
    load_candidate_records,
    main,
    project_clinical_candidates,
    render_ind_dossier,
)

COMPLETE_PK = {
    "dose_mg_kg": 10.0,
    "cl_l_h_kg": 1.8,
    "vd_l_kg": 5.2,
    "brain_kp": 0.45,
    "oral_f": 35.0,
    "noael_mg_kg": 20.0,
    "noael_km": 6.0,
}


class TestIndProjection(unittest.TestCase):
    def test_calculate_allometric_scaling(self) -> None:
        (
            hed,
            daily_mg,
            cl_h,
            vd_h,
            t12_h,
            noael,
            noael_hed,
            mrsd,
        ) = calculate_allometric_scaling(
            animal_dose_mg_kg=10.0,
            mouse_cl_l_h_kg=1.8,
            mouse_vd_l_kg=5.2,
            noael_mg_kg=50.0,
            noael_km=3.0,
        )
        # HED = 10 * (3 / 37) = 0.8108...
        self.assertAlmostEqual(hed, 0.811, places=2)
        # Daily mg = 0.8108 * 60 = 48.65...
        self.assertAlmostEqual(daily_mg, 48.6, places=1)
        self.assertGreater(cl_h, 0.0)
        self.assertGreater(vd_h, 0.0)
        self.assertGreater(t12_h, 0.0)
        self.assertEqual(noael, 50.0)
        self.assertAlmostEqual(mrsd, (noael_hed / 10.0) * 60.0, places=3)

    def test_allometric_scaling_rejects_missing_measurements(self) -> None:
        for kwargs in (
            {"animal_dose_mg_kg": 0.0},
            {"mouse_cl_l_h_kg": 0.0},
            {"mouse_vd_l_kg": -1.0},
            {"noael_mg_kg": 0.0},
        ):
            values = {
                "animal_dose_mg_kg": 10.0,
                "mouse_cl_l_h_kg": 1.8,
                "mouse_vd_l_kg": 5.2,
                "noael_mg_kg": 50.0,
            }
            values.update(kwargs)
            with self.assertRaises(IndError):
                calculate_allometric_scaling(**values)

    def test_project_clinical_candidates_with_complete_evidence(self) -> None:
        compounds = [
            {"compound_id": "PNN-101", "name": "Compound 101"},
            {"compound_id": "PNN-102", "name": "Compound 102"},
        ]
        pk_data = {
            "PNN-101": dict(COMPLETE_PK),
            "PNN-102": {**COMPLETE_PK, "dose_mg_kg": 15.0, "brain_kp": 0.32, "oral_f": 22.0},
        }
        projections = project_clinical_candidates(compounds, pk_data, lead_id="PNN-101")
        self.assertEqual(len(projections), 2)
        self.assertEqual(projections[0].compound_id, "PNN-101")
        self.assertEqual(projections[0].evidence_status, "COMPLETE")
        self.assertIsNotNone(projections[0].phase1_mrsd_mg)
        self.assertIsNotNone(projections[0].therapeutic_index)
        row = projections[0].as_row()
        self.assertEqual(row["evidence_status"], "COMPLETE")
        self.assertIn("human_total_daily_dose_mg", row)

    def test_missing_evidence_is_never_invented(self) -> None:
        projections = project_clinical_candidates([{"compound_id": "PNN-101", "name": "X"}], {})
        projection = projections[0]
        self.assertEqual(projection.evidence_status, "INSUFFICIENT_EVIDENCE")
        for field in (
            "animal_dose_mg_kg",
            "human_equivalent_dose_mg_kg",
            "human_total_daily_dose_mg",
            "projected_human_cl_l_h",
            "projected_human_vd_l",
            "projected_human_half_life_h",
            "noael_mg_kg",
            "noael_hed_mg_kg",
            "phase1_mrsd_mg",
            "therapeutic_index",
            "brain_kp",
            "oral_f_pct",
        ):
            self.assertIsNone(getattr(projection, field), field)
        row = projection.as_row()
        self.assertEqual(row["brain_kp"], "")
        self.assertEqual(row["phase1_mrsd_mg"], "")
        self.assertIn("Brain partitioning unmeasured", projection.findings)
        self.assertIn("NOAEL unmeasured", projection.findings)

    def test_lead_is_nominated_first(self) -> None:
        compounds = [
            {"compound_id": "A"},
            {"compound_id": "B"},
            {"compound_id": "C"},
        ]
        pk_data = {
            "A": dict(COMPLETE_PK),
            "B": {**COMPLETE_PK, "brain_kp": 0.9},
            "C": {**COMPLETE_PK, "brain_kp": 0.7},
        }
        projections = project_clinical_candidates(compounds, pk_data, lead_id="C")
        self.assertEqual(projections[0].compound_id, "C")

    def test_dossier_marks_unmeasured_and_omits_claims(self) -> None:
        projections = project_clinical_candidates([{"compound_id": "PNN-101", "name": "X"}], {})
        dossier = render_ind_dossier(projections, lead_id="PNN-101")
        self.assertIn("UNMEASURED", dossier)
        self.assertIn("INSUFFICIENT_EVIDENCE", dossier)
        self.assertNotIn("Robust central brain exposure", dossier)
        self.assertNotIn("Broad therapeutic window", dossier)
        self.assertIn("Starting Dose**: pending", dossier)

    def test_load_candidate_records_ignores_pk_timecourse_without_measurements(self) -> None:
        vendor_csv = Path("data/compounds/vendor_library.csv")
        pk_csv = Path("data/preclinical/pk_timecourse_template.csv")
        if not vendor_csv.is_file() or not pk_csv.is_file():
            return
        _, pk_data = load_candidate_records(vendor_csv, pk_csv)
        for cid in ("PNN-101", "PNN-102"):
            self.assertNotIn("brain_kp", pk_data.get(cid, {}))
            self.assertNotIn("oral_f", pk_data.get(cid, {}))
            self.assertNotIn("cl_l_h_kg", pk_data.get(cid, {}))

    def test_tox_path_supplies_most_sensitive_noael(self) -> None:
        vendor_csv = Path("data/compounds/vendor_library.csv")
        tox_csv = Path("data/toxicology/repeat_dose_tox_template.csv")
        if not vendor_csv.is_file() or not tox_csv.is_file():
            return
        _, pk_data = load_candidate_records(vendor_csv, tox_path=tox_csv)
        # PNN-101 has rat 50 mg/kg (Km 6 -> HED 8.108) and dog 30 mg/kg (Km 20 -> HED 16.216).
        # The rat study is more sensitive on the HED scale even though 50 > 30 mg/kg.
        self.assertAlmostEqual(pk_data["PNN-101"]["noael_mg_kg"], 50.0)
        self.assertAlmostEqual(pk_data["PNN-101"]["noael_km"], 6.0)

    def test_compile_ind_dossier_uses_measured_candidates(self) -> None:
        vendor_csv = Path("data/compounds/vendor_library.csv")
        tox_csv = Path("data/toxicology/repeat_dose_tox_template.csv")
        if not vendor_csv.is_file() or not tox_csv.is_file():
            return
        extra_pk = {"PNN-101": {"dose_mg_kg": 10.0, "cl_l_h_kg": 1.8, "vd_l_kg": 5.2}}
        projections, dossier_md = compile_ind_dossier(
            "PNN-101", vendor_csv, tox_path=tox_csv, extra_pk=extra_pk
        )
        lead = next(p for p in projections if p.compound_id == "PNN-101")
        self.assertEqual(lead.evidence_status, "PARTIAL")  # Kp and F remain unmeasured
        self.assertAlmostEqual(lead.noael_mg_kg or 0.0, 50.0)
        self.assertAlmostEqual(lead.noael_hed_mg_kg or 0.0, 50.0 * 6.0 / 37.0, places=6)
        self.assertIsNotNone(lead.phase1_mrsd_mg)
        self.assertIn("Investigational New Drug (IND)", dossier_md)
        self.assertIn("ICH S7A/S7B", dossier_md)
        self.assertIn("Phase 1 Clinical Protocol Architecture", dossier_md)


class TestIndCli(unittest.TestCase):
    def test_ind_cli(self) -> None:
        vendor_csv = Path("data/compounds/vendor_library.csv")
        pk_csv = Path("data/preclinical/pk_timecourse_template.csv")
        tox_csv = Path("data/toxicology/repeat_dose_tox_template.csv")
        if not vendor_csv.is_file():
            return

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_csv = Path(tmp_dir) / "doses.csv"
            out_md = Path(tmp_dir) / "dossier.md"
            out_json = Path(tmp_dir) / "summary.json"

            argv = [
                "--compounds", str(vendor_csv),
                "--out", str(out_csv),
                "--dossier", str(out_md),
                "--json", str(out_json),
            ]
            if pk_csv.is_file():
                argv.extend(["--candidates", str(pk_csv)])
            if tox_csv.is_file():
                argv.extend(["--tox", str(tox_csv)])

            exit_code = main(argv)
            self.assertEqual(exit_code, 0)
            self.assertTrue(out_csv.is_file())
            self.assertTrue(out_md.is_file())
            self.assertTrue(out_json.is_file())
            summary = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertIn("by_evidence_status", summary)


if __name__ == "__main__":
    unittest.main()
