import io
import math
import tempfile
import unittest
from pathlib import Path

from pnn.preclinical import (
    ConcentrationPoint,
    PdOutcome,
    PkProfile,
    PreclinicalCandidate,
    PreclinicalError,
    calculate_nca,
    evaluate_candidates,
    load_pd_data,
    load_pk_data,
    main,
    render_preclinical_report,
)


class TestPreclinicalNCA(unittest.TestCase):
    def test_nca_empty_points_raises(self) -> None:
        with self.assertRaises(PreclinicalError):
            calculate_nca([], dose_mg_kg=10.0, compound_id="PNN-TEST", route="PO", tissue="plasma")

    def test_nca_trapezoidal_auc_exact(self) -> None:
        # Simple trapezoid: (0, 0), (1, 100), (2, 50)
        # Segment 1: (0 + 100)/2 * 1 = 50
        # Segment 2: (100 + 50)/2 * 1 = 75
        # Total AUC_0_t = 125
        pts = [
            ConcentrationPoint(time_h=0.0, concentration=0.0),
            ConcentrationPoint(time_h=1.0, concentration=100.0),
            ConcentrationPoint(time_h=2.0, concentration=50.0),
        ]
        profile = calculate_nca(pts, dose_mg_kg=10.0, compound_id="PNN-TEST", route="PO", tissue="plasma")
        self.assertEqual(profile.cmax, 100.0)
        self.assertEqual(profile.tmax_h, 1.0)
        self.assertAlmostEqual(profile.auc_0_t, 125.0, places=4)

    def test_nca_terminal_elimination_and_clearance(self) -> None:
        # Exponential decay: C(t) = 1000 * exp(-0.693147 * t) -> half_life = 1.0 h, lambda_z = 0.693147
        pts = [
            ConcentrationPoint(time_h=0.0, concentration=1000.0),
            ConcentrationPoint(time_h=1.0, concentration=500.0),
            ConcentrationPoint(time_h=2.0, concentration=250.0),
            ConcentrationPoint(time_h=3.0, concentration=125.0),
        ]
        profile = calculate_nca(pts, dose_mg_kg=2.0, compound_id="PNN-IV", route="IV", tissue="plasma")
        self.assertIsNotNone(profile.lambda_z)
        self.assertAlmostEqual(profile.lambda_z or 0.0, 0.693147, places=2)
        self.assertIsNotNone(profile.half_life_h)
        self.assertAlmostEqual(profile.half_life_h or 0.0, 1.0, places=2)
        self.assertIsNotNone(profile.auc_0_inf)
        self.assertIsNotNone(profile.cl_l_h_kg)
        self.assertIsNotNone(profile.vd_l_kg)


class TestPreclinicalDataIngestion(unittest.TestCase):
    def test_load_pk_data_missing_file_raises(self) -> None:
        with self.assertRaises(PreclinicalError):
            load_pk_data(Path("non_existent_pk.csv"))

    def test_load_pk_data_missing_columns_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            f = Path(tmp_dir) / "invalid_pk.csv"
            f.write_text("compound_id,dose_mg_kg\nPNN-101,10.0\n", encoding="utf-8")
            with self.assertRaises(PreclinicalError):
                load_pk_data(f)

    def test_load_pk_data_invalid_numeric_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            f = Path(tmp_dir) / "corrupt_pk.csv"
            f.write_text(
                "compound_id,route,tissue,dose_mg_kg,time_h,concentration_ng_ml\n"
                "PNN-101,PO,plasma,10.0,1.0,NOT_A_NUM\n",
                encoding="utf-8",
            )
            with self.assertRaises(PreclinicalError):
                load_pk_data(f)

    def test_load_pd_data_valid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            f = Path(tmp_dir) / "pd.csv"
            f.write_text(
                "compound_id,dose_mg_kg,wfa_pnn_reduction_pct,itt_glucose_auc_reduction_pct\n"
                "PNN-101,10.0,35.0,28.0\n",
                encoding="utf-8",
            )
            outcomes = load_pd_data(f)
            self.assertIn("PNN-101", outcomes)
            self.assertEqual(outcomes["PNN-101"].wfa_pnn_reduction_pct, 35.0)


class TestPreclinicalCandidateEvaluation(unittest.TestCase):
    def setUp(self) -> None:
        # Create mock PK profiles
        self.iv_plasma = PkProfile(
            compound_id="PNN-101",
            route="IV",
            tissue="plasma",
            dose_mg_kg=2.0,
            cmax=1000.0,
            tmax_h=0.08,
            auc_0_t=1000.0,
            auc_0_inf=1100.0,
            half_life_h=2.0,
            lambda_z=0.35,
            cl_l_h_kg=1.8,
            vd_l_kg=5.1,
        )
        self.po_plasma = PkProfile(
            compound_id="PNN-101",
            route="PO",
            tissue="plasma",
            dose_mg_kg=10.0,
            cmax=400.0,
            tmax_h=1.0,
            auc_0_t=2500.0,  # F = (2500/10) / (1000/2) * 100 = 250 / 500 * 100 = 50%
            auc_0_inf=2600.0,
            half_life_h=2.2,
            lambda_z=0.32,
            cl_l_h_kg=None,
            vd_l_kg=None,
        )
        self.po_brain = PkProfile(
            compound_id="PNN-101",
            route="PO",
            tissue="brain",
            dose_mg_kg=10.0,
            cmax=250.0,
            tmax_h=1.5,
            auc_0_t=1500.0,  # Kp,brain = 1500 / 2500 = 0.60
            auc_0_inf=1600.0,
            half_life_h=2.4,
            lambda_z=0.29,
            cl_l_h_kg=None,
            vd_l_kg=None,
        )

    def test_evaluate_clinical_candidate(self) -> None:
        pk_map = {"PNN-101": {"IV_plasma": self.iv_plasma, "PO_plasma": self.po_plasma, "PO_brain": self.po_brain}}
        pd_map = {
            "PNN-101": PdOutcome(
                compound_id="PNN-101",
                dose_mg_kg=10.0,
                wfa_pnn_reduction_pct=35.0,
                itt_glucose_auc_reduction_pct=30.0,
            )
        }
        candidates = evaluate_candidates(pk_map, pd_map, min_kp_brain=0.30, min_oral_f_pct=20.0)
        self.assertEqual(len(candidates), 1)
        c = candidates[0]
        self.assertEqual(c.status, "CLINICAL_CANDIDATE")
        self.assertAlmostEqual(c.oral_bioavailability_pct or 0.0, 50.0, places=1)
        self.assertAlmostEqual(c.brain_plasma_ratio or 0.0, 0.60, places=2)

    def test_evaluate_insufficient_brain_penetration(self) -> None:
        # Low brain AUC -> Kp = 100 / 2500 = 0.04 < 0.20
        po_brain_low = PkProfile(
            compound_id="PNN-101",
            route="PO",
            tissue="brain",
            dose_mg_kg=10.0,
            cmax=20.0,
            tmax_h=1.5,
            auc_0_t=100.0,
            auc_0_inf=110.0,
            half_life_h=2.4,
            lambda_z=0.29,
            cl_l_h_kg=None,
            vd_l_kg=None,
        )
        pk_map = {"PNN-101": {"IV_plasma": self.iv_plasma, "PO_plasma": self.po_plasma, "PO_brain": po_brain_low}}
        candidates = evaluate_candidates(pk_map, min_kp_brain=0.30, min_oral_f_pct=20.0)
        self.assertEqual(candidates[0].status, "INSUFFICIENT_BRAIN_PENETRATION")

    def test_render_preclinical_report(self) -> None:
        candidate = PreclinicalCandidate(
            compound_id="PNN-101",
            status="CLINICAL_CANDIDATE",
            oral_bioavailability_pct=45.0,
            brain_plasma_ratio=0.65,
            plasma_cmax_ng_ml=500.0,
            plasma_auc_h_ng_ml=2400.0,
            brain_auc_h_ng_ml=1560.0,
            half_life_h=3.0,
            wfa_pnn_reduction_pct=38.0,
            itt_glucose_auc_reduction_pct=32.0,
            findings="High brain exposure; Robust oral F; Potent metabolic reversal",
        )
        report = render_preclinical_report([candidate], {"min_kp_brain": 0.30})
        self.assertIn("# PNN Phase 4 In Vivo Preclinical Candidate Selection Report", report)
        self.assertIn("G-001", report)
        self.assertIn("CLINICAL_CANDIDATE", report)
        self.assertIn("PNN-101", report)


class TestPreclinicalCLI(unittest.TestCase):
    def test_cli_execution_success(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            pk_path = Path("data/preclinical/pk_timecourse_template.csv")
            pd_path = Path("data/preclinical/pd_metabolic_template.csv")
            if not pk_path.is_file():
                return
            out_csv = Path(tmp_dir) / "cand.csv"
            out_json = Path(tmp_dir) / "cand.json"
            out_rep = Path(tmp_dir) / "cand.md"
            rc = main([
                "--pk", str(pk_path),
                "--pd", str(pd_path),
                "--out", str(out_csv),
                "--json", str(out_json),
                "--report", str(out_rep),
            ])
            self.assertEqual(rc, 0)
            self.assertTrue(out_csv.is_file())
            self.assertTrue(out_json.is_file())
            self.assertTrue(out_rep.is_file())


if __name__ == "__main__":
    unittest.main()
