"""Unit tests for MD trajectory and contact stability analyzer."""

import tempfile
import unittest
from pathlib import Path

from pnn.dynamics import (
    DynamicsError,
    TrajectoryPoint,
    analyze_trajectory,
    generate_trajectory_svg,
    load_trajectory_csv,
    main,
    render_trajectory_report,
    run_dynamics_pipeline,
)


class TestDynamics(unittest.TestCase):
    def test_load_trajectory_csv(self) -> None:
        template_csv = Path("data/dynamics/ptpsigma_md_trajectory_template.csv")
        if not template_csv.is_file():
            return

        points = load_trajectory_csv(template_csv)
        self.assertGreaterEqual(len(points), 10)
        self.assertEqual(points[0].time_ns, 0.0)
        self.assertGreater(points[-1].time_ns, 50.0)

    def test_load_trajectory_missing_raises(self) -> None:
        with self.assertRaises(DynamicsError):
            load_trajectory_csv(Path("non_existent_trajectory.csv"))

    def test_analyze_trajectory_stable(self) -> None:
        points = [
            TrajectoryPoint(time_ns=float(i), ligand_rmsd=1.2 + 0.1 * (i % 3), protein_rmsd=1.1, contacts_count=4, mmgbsa_dg=-34.0, residue_contacts={"K201": True, "R204": True})
            for i in range(20)
        ]
        analysis = analyze_trajectory(points, compound_id="PNN-101", target_id="PTPsigma")
        self.assertEqual(analysis.stability_class, "STABLE_POSE")
        self.assertGreaterEqual(analysis.stability_score, 85.0)
        self.assertEqual(analysis.compound_id, "PNN-101")
        self.assertAlmostEqual(analysis.residue_persistence["K201"], 100.0)

        row = analysis.as_row()
        self.assertIn("mean_ligand_rmsd", row)
        self.assertIn("stability_class", row)

    def test_analyze_trajectory_unstable_drift(self) -> None:
        # Simulate severe drift unbinding event (RMSD creeping up to 4.5 A)
        points = [
            TrajectoryPoint(time_ns=float(i), ligand_rmsd=1.2 + (i * 0.2), protein_rmsd=1.5, contacts_count=max(0, 4 - i // 4), mmgbsa_dg=-10.0)
            for i in range(20)
        ]
        analysis = analyze_trajectory(points, compound_id="DRIFT-01", target_id="PTPsigma")
        self.assertEqual(analysis.stability_class, "UNSTABLE_DRIFT")
        self.assertLess(analysis.stability_score, 60.0)

    def test_generate_trajectory_svg(self) -> None:
        points = [
            TrajectoryPoint(time_ns=float(i), ligand_rmsd=1.3, protein_rmsd=1.2, contacts_count=3, mmgbsa_dg=-32.0)
            for i in range(10)
        ]
        analysis = analyze_trajectory(points, "PNN-101", "PTPsigma")
        svg = generate_trajectory_svg(analysis, points)
        self.assertIn("<svg", svg)
        self.assertIn("Stability Gate (2.0", svg)
        self.assertIn("PTPsigma", svg)

    def test_render_trajectory_report(self) -> None:
        points = [
            TrajectoryPoint(time_ns=float(i), ligand_rmsd=1.4, protein_rmsd=1.2, contacts_count=4, mmgbsa_dg=-35.0, residue_contacts={"K201": True})
            for i in range(10)
        ]
        analysis = analyze_trajectory(points, "PNN-101", "PTPsigma")
        report = render_trajectory_report(analysis)
        self.assertIn("# Molecular Dynamics Trajectory Stability Report", report)
        self.assertIn("PNN-101", report)
        self.assertIn("PTPsigma", report)
        self.assertIn("Key Metric Summary", report)

    def test_run_dynamics_pipeline_and_cli(self) -> None:
        template_csv = Path("data/dynamics/ptpsigma_md_trajectory_template.csv")
        if not template_csv.is_file():
            return

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "dyn_out"
            results = run_dynamics_pipeline(
                trajectory_path=template_csv,
                compound_id="PNN-101",
                target_id="PTPsigma",
                out_dir=out_dir,
            )
            self.assertIn("analysis", results)
            self.assertTrue((out_dir / "PNN-101_trajectory_metrics.csv").is_file())
            self.assertTrue((out_dir / "PNN-101_dynamics_report.md").is_file())
            self.assertTrue((out_dir / "PNN-101_rmsd_trajectory.svg").is_file())

            # Test CLI
            exit_code = main([
                "--trajectory", str(template_csv),
                "--compound", "PNN-101",
                "--out-dir", str(out_dir),
            ])
            self.assertEqual(exit_code, 0)

    def test_missing_optional_columns_produce_insufficient_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            trajectory = Path(tmp_dir) / "ligand_only.csv"
            trajectory.write_text(
                "time_ns,ligand_rmsd\n0,1.0\n5,1.1\n10,1.2\n15,1.1\n20,1.0\n",
                encoding="utf-8",
            )
            points = load_trajectory_csv(trajectory)
            analysis = analyze_trajectory(points, "PNN-101", "PTPsigma")
            self.assertEqual(analysis.stability_class, "INSUFFICIENT_EVIDENCE")
            self.assertEqual(analysis.evidence_status, "INSUFFICIENT_EVIDENCE")
            self.assertIsNone(analysis.stability_score)
            self.assertIsNone(analysis.mean_contacts)
            self.assertIsNone(analysis.min_contacts)
            self.assertIsNone(analysis.mean_protein_rmsd)
            self.assertIsNone(analysis.mean_mmgbsa_dg)
            row = analysis.as_row()
            self.assertEqual(row["stability_score"], "")
            self.assertEqual(row["mean_contacts"], "")
            self.assertEqual(row["mean_mmgbsa_dg"], "")
            report = render_trajectory_report(analysis)
            self.assertIn("UNMEASURED", report)
            self.assertNotIn("STABLE_POSE", report)

    def test_partial_evidence_keeps_verdict_but_flags_gaps(self) -> None:
        points = [
            TrajectoryPoint(time_ns=float(i), ligand_rmsd=1.2, contacts_count=4)
            for i in range(10)
        ]
        analysis = analyze_trajectory(points, "PNN-101", "PTPsigma")
        self.assertEqual(analysis.stability_class, "STABLE_POSE")
        self.assertEqual(analysis.evidence_status, "PARTIAL")
        self.assertIsNone(analysis.mean_protein_rmsd)
        self.assertIsNone(analysis.mean_mmgbsa_dg)


if __name__ == "__main__":
    unittest.main()
