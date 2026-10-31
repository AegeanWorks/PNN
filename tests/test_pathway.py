"""Unit tests for hypothalamic circuit dynamics and dosing simulator."""

import tempfile
import unittest
from pathlib import Path

from pnn.pathway import (
    PathwayError,
    main,
    render_pathway_report,
    simulate_circuit_dynamics,
)


class TestPathwaySimulator(unittest.TestCase):
    def test_simulate_circuit_dynamics_pulsatile(self) -> None:
        points, summary = simulate_circuit_dynamics(
            compound_id="PNN-101",
            dose_mg_kg=10.0,
            regimen="INTERMITTENT_PULSATILE",
            duration_days=7,
            dt_h=0.5,
        )
        self.assertGreater(len(points), 10)
        self.assertEqual(summary.compound_id, "PNN-101")
        self.assertEqual(summary.regimen_name, "INTERMITTENT_PULSATILE")
        self.assertGreater(summary.mean_pnn_reduction_pct, 0.0)
        self.assertGreater(summary.peak_pnn_reduction_pct, summary.mean_pnn_reduction_pct)
        self.assertGreater(summary.plasticity_safety_score, 80.0)
        self.assertIn("RECOMMENDED", summary.recommendation)

        row = points[0].as_row()
        self.assertIn("time_h", row)
        self.assertIn("pnn_density_pct", row)
        self.assertIn("akt_signaling_pct", row)

    def test_simulate_circuit_dynamics_continuous(self) -> None:
        points, summary = simulate_circuit_dynamics(
            compound_id="PNN-102",
            dose_mg_kg=15.0,
            regimen="CONTINUOUS_DAILY",
            duration_days=5,
            dt_h=0.5,
        )
        self.assertGreater(len(points), 5)
        self.assertEqual(summary.regimen_name, "CONTINUOUS_DAILY")
        self.assertGreaterEqual(summary.peak_pnn_reduction_pct, 40.0)

    def test_regimen_variations(self) -> None:
        for reg in ("TWICE_WEEKLY", "ONCE_WEEKLY", "SINGLE_DOSE"):
            pts, sm = simulate_circuit_dynamics(
                compound_id="PNN-101",
                dose_mg_kg=10.0,
                regimen=reg,
                duration_days=7,
                dt_h=1.0,
            )
            self.assertGreater(len(pts), 0)
            self.assertEqual(sm.regimen_name, reg)

    def test_invalid_parameters_raise(self) -> None:
        with self.assertRaises(PathwayError):
            simulate_circuit_dynamics(dose_mg_kg=-5.0)

        with self.assertRaises(PathwayError):
            simulate_circuit_dynamics(duration_days=0)

        with self.assertRaises(PathwayError):
            simulate_circuit_dynamics(duration_days=100)

        with self.assertRaises(PathwayError):
            simulate_circuit_dynamics(dt_h=0.0)

        with self.assertRaises(PathwayError):
            simulate_circuit_dynamics(regimen="UNKNOWN_REGIMEN")

        with self.assertRaises(PathwayError):
            simulate_circuit_dynamics(half_life_h=0.0)

    def test_untreated_baseline_has_no_restoration_or_weight_loss(self) -> None:
        _, summary = simulate_circuit_dynamics(
            compound_id="PNN-101",
            dose_mg_kg=0.0001,
            regimen="SINGLE_DOSE",
            duration_days=14,
            dt_h=0.5,
        )
        self.assertAlmostEqual(summary.mean_akt_restoration_pct, 0.0, places=1)
        self.assertAlmostEqual(summary.projected_weight_loss_pct, 0.0, places=2)
        self.assertAlmostEqual(summary.mean_pnn_reduction_pct, 0.0, places=1)
        self.assertIn("SUB-THERAPEUTIC", summary.recommendation)

    def test_continuous_regimen_is_flagged_for_plasticity_risk(self) -> None:
        _, summary = simulate_circuit_dynamics(
            compound_id="PNN-101",
            dose_mg_kg=10.0,
            regimen="CONTINUOUS_DAILY",
            duration_days=14,
            dt_h=0.5,
        )
        self.assertLess(summary.plasticity_safety_score, 80.0)
        self.assertIn("CAUTION", summary.recommendation)

    def test_render_pathway_report(self) -> None:
        points, summary = simulate_circuit_dynamics(
            compound_id="PNN-101",
            dose_mg_kg=10.0,
            duration_days=3,
            dt_h=1.0,
        )
        report = render_pathway_report(summary, points)
        self.assertIn("# PNN Hypothalamic Circuit & Dosing Dynamics Report", report)
        self.assertIn("PNN-101", report)
        self.assertIn("Plasticity Risk Avoidance", report)
        self.assertIn("## Interpretation", report)
        self.assertNotIn("optimal therapeutic index", report)
        self.assertNotIn("Runge-Kutta", report)

    def test_pathway_cli(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_csv = Path(tmp_dir) / "sim.csv"
            out_rep = Path(tmp_dir) / "rep.md"
            out_json = Path(tmp_dir) / "sum.json"

            exit_code = main([
                "--compound", "PNN-101",
                "--dose", "10.0",
                "--regimen", "INTERMITTENT_PULSATILE",
                "--days", "3",
                "--out", str(out_csv),
                "--report", str(out_rep),
                "--json", str(out_json),
            ])
            self.assertEqual(exit_code, 0)
            self.assertTrue(out_csv.is_file())
            self.assertTrue(out_rep.is_file())
            self.assertTrue(out_json.is_file())


if __name__ == "__main__":
    unittest.main()
