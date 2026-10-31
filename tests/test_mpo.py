import tempfile
import unittest
from pathlib import Path

from pnn.mpo import (
    MpoError,
    MpoProfile,
    calculate_cns_mpo,
    main,
    profile_candidate,
    render_mpo_report,
    run_mpo_analysis,
)


class TestMpoCalculations(unittest.TestCase):
    def test_cns_mpo_optimal_properties(self) -> None:
        # Optimal CNS property values: MW <= 360, logP <= 3.0, logD <= 2.0, TPSA 40-90, HBD <= 0.5, pKa <= 8.0
        score = calculate_cns_mpo(mw=310.0, logp=2.2, tpsa=55.0, hbd=0, pka=7.6, logd=1.7)
        self.assertAlmostEqual(score, 6.0, places=2)

    def test_cns_mpo_penalized_properties(self) -> None:
        # Exceeding all bounds: MW > 500, logP > 5.0, logD > 4.0, TPSA > 120, HBD > 3.5, pKa > 10.0
        score = calculate_cns_mpo(mw=550.0, logp=5.8, tpsa=135.0, hbd=4, pka=10.5, logd=4.5)
        self.assertLess(score, 1.0)

    def test_profile_development_candidate(self) -> None:
        profile = profile_candidate(
            compound_id="PNN-101",
            name="HypoSprint-1",
            mw=349.46,
            logp=2.65,
            tpsa=48.4,
            hbd=2,
            pka=8.2,
            ec50_nm=24.0,
            brain_kp=0.65,
            oral_f_pct=42.0,
            selectivity_index=3540.0,
        )
        self.assertEqual(profile.mpo_tier, "DEVELOPMENT_CANDIDATE")
        self.assertGreaterEqual(profile.composite_desirability, 0.70)
        self.assertGreaterEqual(profile.cns_mpo_score, 4.5)

    def test_profile_deprioritized_liability_detection(self) -> None:
        profile = profile_candidate(
            compound_id="PNN-REJ-1",
            name="HeavyCluster-X",
            mw=515.72,
            logp=6.80,
            tpsa=95.2,
            hbd=2,
            pka=4.2,
            ec50_nm=8500.0,
            brain_kp=0.02,
            oral_f_pct=5.0,
            selectivity_index=1.2,
        )
        self.assertEqual(profile.mpo_tier, "DEPRIORITIZED")
        self.assertLess(profile.composite_desirability, 0.45)
        self.assertIn("High molecular weight", profile.primary_liability)


class TestMpoAnalysisAndCLI(unittest.TestCase):
    def test_missing_compounds_raises(self) -> None:
        with self.assertRaises(MpoError):
            run_mpo_analysis(Path("non_existent_compounds.csv"))

    def test_run_mpo_analysis_valid(self) -> None:
        comp_path = Path("data/compounds/vendor_library.csv")
        cand_path = Path("data/preclinical/pk_timecourse_template.csv")
        if not comp_path.is_file():
            return
        profiles = run_mpo_analysis(compounds_path=comp_path, candidates_path=cand_path)
        self.assertTrue(len(profiles) > 0)
        pnn101 = next((p for p in profiles if p.compound_id == "PNN-101"), None)
        self.assertIsNotNone(pnn101)
        self.assertEqual(pnn101.mpo_tier, "DEVELOPMENT_CANDIDATE")  # type: ignore[union-attr]

    def test_render_mpo_report(self) -> None:
        profile = profile_candidate(
            compound_id="PNN-101",
            name="HypoSprint-1",
            mw=349.46,
            logp=2.65,
            tpsa=48.4,
            hbd=2,
            pka=8.2,
            ec50_nm=24.0,
            brain_kp=0.65,
            oral_f_pct=42.0,
            selectivity_index=3540.0,
        )
        report = render_mpo_report([profile])
        self.assertIn("# PNN Phase 5 Lead Optimization", report)
        self.assertIn("DEVELOPMENT_CANDIDATE", report)
        self.assertIn("PNN-101", report)

    def test_cli_execution_success(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            comp_path = Path("data/compounds/vendor_library.csv")
            if not comp_path.is_file():
                return
            out_csv = Path(tmp_dir) / "mpo.csv"
            out_rep = Path(tmp_dir) / "mpo.md"
            out_json = Path(tmp_dir) / "mpo.json"
            rc = main([
                "--compounds", str(comp_path),
                "--out", str(out_csv),
                "--report", str(out_rep),
                "--json", str(out_json),
            ])
            self.assertEqual(rc, 0)
            self.assertTrue(out_csv.is_file())
            self.assertTrue(out_rep.is_file())
            self.assertTrue(out_json.is_file())


if __name__ == "__main__":
    unittest.main()
