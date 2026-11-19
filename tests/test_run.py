import contextlib
import csv
import hashlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pnn.preclinical import PreclinicalCandidate
from pnn.runner import (
    RunSummary,
    RunnerError,
    _candidate_exposures,
    execute_turnkey_pipeline,
    is_daily_timepoint,
    main,
    print_run_banner,
)


class TestTurnkeyRunner(unittest.TestCase):
    def test_missing_vendor_library_raises(self) -> None:
        with self.assertRaises(RunnerError):
            execute_turnkey_pipeline(
                vendor_csv=Path("non_existent_vendor.csv"),
                repurposing_csv=None,
                docking_csv=Path("data/docking/multi_target_scores.csv"),
                degradation_csv=None,
                permeability_csv=None,
                cytotoxicity_csv=None,
                pk_csv=None,
                pd_csv=None,
                out_dir=Path("tmp_out"),
            )

    def test_execute_turnkey_pipeline_success(self) -> None:
        vendor_csv = Path("data/compounds/vendor_library.csv")
        repurposing_csv = Path("data/compounds/repurposing_library.csv")
        docking_csv = Path("data/docking/multi_target_scores.csv")
        degradation_csv = Path("data/assays/degradation_template.csv")
        permeability_csv = Path("data/assays/permeability_template.csv")
        cytotoxicity_csv = Path("data/assays/cytotoxicity_template.csv")
        pk_csv = Path("data/preclinical/pk_timecourse_template.csv")
        pd_csv = Path("data/preclinical/pd_metabolic_template.csv")
        tox_csv = Path("data/toxicology/repeat_dose_tox_template.csv")
        safety_pharm_csv = Path("data/toxicology/safety_pharmacology_template.csv")

        if not vendor_csv.is_file() or not docking_csv.is_file():
            return

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "output"
            summary = execute_turnkey_pipeline(
                vendor_csv=vendor_csv,
                repurposing_csv=repurposing_csv,
                docking_csv=docking_csv,
                degradation_csv=degradation_csv,
                permeability_csv=permeability_csv,
                cytotoxicity_csv=cytotoxicity_csv,
                pk_csv=pk_csv,
                pd_csv=pd_csv,
                out_dir=out_dir,
                tox_csv=tox_csv,
                safety_pharm_csv=safety_pharm_csv,
            )
            self.assertGreater(summary.curated_count, 0)
            self.assertGreater(summary.shortlisted_count, 0)
            self.assertGreater(summary.scaffolds_count, 0)
            self.assertGreater(summary.safety_assessed_count, 0)
            self.assertTrue((out_dir / "shortlist.csv").is_file())
            self.assertTrue((out_dir / "lead_report.md").is_file())
            self.assertTrue((out_dir / "qualified_hits.csv").is_file())
            self.assertTrue((out_dir / "preclinical_candidates.csv").is_file())
            self.assertTrue((out_dir / "sar_clusters.csv").is_file())
            self.assertTrue((out_dir / "studio.html").is_file())
            studio_html = (out_dir / "studio.html").read_text(encoding="utf-8")
            self.assertIn(
                "Circuit dynamics data: PNN density and pAKT signaling over time", studio_html
            )
            self.assertTrue((out_dir / "safety_assessments.csv").is_file())
            self.assertTrue((out_dir / "safety_package.md").is_file())

            # The dossier must nominate the computed clinical candidate, not a hardcoded default.
            dossier = (out_dir / "ind_dossier.md").read_text(encoding="utf-8")
            self.assertIn("**Nominated Lead Compound**: `PNN-REP-01`", dossier)
            self.assertIn("**Evidence Status**: **COMPLETE**", dossier)
            with (out_dir / "human_dose_projections.csv").open("r", encoding="utf-8", newline="") as handle:
                projection_rows = list(csv.DictReader(handle))
            self.assertEqual(projection_rows[0]["compound_id"], "PNN-REP-01")
            self.assertEqual(projection_rows[0]["evidence_status"], "COMPLETE")
            self.assertEqual(projection_rows[0]["brain_kp"], "0.73")

            # The IND MRSD and the safety NOAEL-based dose must reconcile for the same lead.
            with (out_dir / "safety_assessments.csv").open("r", encoding="utf-8", newline="") as handle:
                safety_rows = {row["compound_id"]: row for row in csv.DictReader(handle)}
            self.assertAlmostEqual(
                float(safety_rows["PNN-REP-01"]["phase1_mrsd_mg"]),
                float(projection_rows[0]["phase1_mrsd_mg"]),
                places=2,
            )

            # Verify banner printing without UnicodeEncodeError
            print_run_banner(summary)

    def test_runner_cli_execution_no_browser(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "run_out"
            rc = main(["--no-browser", "--out-dir", str(out_dir)])
            self.assertEqual(rc, 0)
            self.assertTrue((out_dir / "studio.html").is_file())
            self.assertTrue((out_dir / "safety_assessments.csv").is_file())
            self.assertTrue((out_dir / "safety_package.md").is_file())

    def test_sar_reports_vina_affinity_not_consensus_score(self) -> None:
        # G-002: a relative 0..1 consensus score must not be reported as a
        # Vina binding affinity (negative kcal/mol) in sar_clusters.csv.
        if not Path("data/compounds/vendor_library.csv").is_file():
            return
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "output"
            execute_turnkey_pipeline(
                vendor_csv=Path("data/compounds/vendor_library.csv"),
                repurposing_csv=Path("data/compounds/repurposing_library.csv"),
                docking_csv=Path("data/docking/multi_target_scores.csv"),
                degradation_csv=None,
                permeability_csv=None,
                cytotoxicity_csv=None,
                pk_csv=None,
                pd_csv=None,
                out_dir=out_dir,
            )
            with (out_dir / "shortlist.csv").open("r", encoding="utf-8", newline="") as handle:
                vina_values = [round(float(r["vina_affinity"]), 2) for r in csv.DictReader(handle)]
            with (out_dir / "sar_clusters.csv").open("r", encoding="utf-8", newline="") as handle:
                clusters = list(csv.DictReader(handle))
            reported = [float(c["best_vina_docking"]) for c in clusters if c["best_vina_docking"]]
            self.assertTrue(reported)
            self.assertTrue(
                all(v < 0.0 for v in reported),
                f"SAR 'best_vina_docking' holds non-kcal/mol values: {reported}",
            )
            self.assertTrue(all(round(v, 2) in vina_values for v in reported))

    def test_turnkey_qualification_uses_reported_ec50_ceiling(self) -> None:
        # assay_report.md states EC50 <= 100 nM; classification must use that
        # same ceiling, not qualify_leads' permissive 10000 nM library default.
        degradation_csv = Path("data/assays/degradation_template.csv")
        if not degradation_csv.is_file():
            return
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "output"
            execute_turnkey_pipeline(
                vendor_csv=Path("data/compounds/vendor_library.csv"),
                repurposing_csv=Path("data/compounds/repurposing_library.csv"),
                docking_csv=Path("data/docking/multi_target_scores.csv"),
                degradation_csv=degradation_csv,
                permeability_csv=Path("data/assays/permeability_template.csv"),
                cytotoxicity_csv=Path("data/assays/cytotoxicity_template.csv"),
                pk_csv=Path("data/preclinical/pk_timecourse_template.csv"),
                pd_csv=Path("data/preclinical/pd_metabolic_template.csv"),
                out_dir=out_dir,
            )
            with (out_dir / "qualified_hits.csv").open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            confirmed = [r for r in rows if r["status"] == "CONFIRMED_HIT"]
            self.assertTrue(confirmed)
            self.assertTrue(
                all(float(r["ec50_nm"]) <= 100.0 for r in confirmed),
                f"CONFIRMED_HIT with EC50 above reported 100 nM ceiling: {confirmed}",
            )
            self.assertNotIn("PNN-104", [r["compound_id"] for r in confirmed])

    def test_confirmed_hits_count_excludes_qualified_leads(self) -> None:
        # QUALIFIED_LEAD means "softening confirmed; pending measured
        # permeability" — explicitly NOT confirmed. The banner label "Confirmed
        # In Vitro Hits" must count only CONFIRMED_HIT, or it overclaims evidence.
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "output"
            deg = Path(tmp_dir) / "deg.csv"
            deg.write_text(
                "compound_id,ec50_nm,max_degradation_pct\n"
                "PNN-101,50.0,80.0\n"
                "PNN-102,60.0,75.0\n",
                encoding="utf-8",
            )
            summary = execute_turnkey_pipeline(
                vendor_csv=Path("data/compounds/vendor_library.csv"),
                repurposing_csv=None,
                docking_csv=Path("data/docking/multi_target_scores.csv"),
                degradation_csv=deg,
                permeability_csv=None,  # no permeability -> QUALIFIED_LEAD
                cytotoxicity_csv=None,
                pk_csv=None,
                pd_csv=None,
                out_dir=out_dir,
            )
            with (out_dir / "qualified_hits.csv").open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            strict_confirmed = sum(1 for r in rows if r["status"] == "CONFIRMED_HIT")
            qualified_leads = [r for r in rows if r["status"] == "QUALIFIED_LEAD"]
            self.assertTrue(qualified_leads, "expected QUALIFIED_LEAD rows without permeability data")
            self.assertEqual(
                summary.qualified_hits_count,
                strict_confirmed,
                f"banner counted QUALIFIED_LEAD as confirmed: {summary.qualified_hits_count} vs {strict_confirmed}",
            )

    def test_banner_accounts_for_safety_caution_results(self) -> None:
        # 5 assessed = 3 safe + 1 caution + 1 failed. The banner must not
        # silently drop the caution and leave 1 of 5 unaccounted for.
        summary = RunSummary(
            curated_count=11,
            screened_count=11,
            shortlisted_count=11,
            qualified_hits_count=5,
            preclinical_candidates_count=1,
            scaffolds_count=14,
            elapsed_seconds=0.2,
            output_dir=Path("output"),
            studio_html_path=Path("output/studio.html"),
            safety_assessed_count=5,
            safety_safe_count=3,
            safety_caution_count=1,
            safety_failed_count=1,
        )
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            print_run_banner(summary)
        banner = out.getvalue()
        self.assertIn("3/5", banner)
        self.assertIn("1 failed", banner)
        self.assertIn("1 caution", banner)

    def test_daily_timepoint_tolerance(self) -> None:
        self.assertTrue(is_daily_timepoint(0.0))
        self.assertTrue(is_daily_timepoint(24.0))
        self.assertTrue(is_daily_timepoint(24.0000000005))
        self.assertTrue(is_daily_timepoint(47.9999999995))
        self.assertFalse(is_daily_timepoint(12.0))
        self.assertFalse(is_daily_timepoint(36.5))

    @patch("webbrowser.open_new_tab")
    def test_runner_cli_execution_auto_opens_browser(self, mock_open) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir) / "run_out"
            rc = main(["--out-dir", str(out_dir)])
            self.assertEqual(rc, 0)
            self.assertTrue(mock_open.called)
            call_arg = mock_open.call_args[0][0]
            self.assertIn("studio.html", call_arg)


class TurnkeyDeterminismTest(unittest.TestCase):
    def test_identical_inputs_produce_identical_artifacts(self) -> None:
        # Repeated turnkey runs with the same inputs must be byte-identical;
        # hidden mutable state or unordered iteration would silently make
        # reports, rankings, or the generated dashboard drift between runs.
        with tempfile.TemporaryDirectory() as tmp:
            def run_once(out_dir: Path) -> dict:
                execute_turnkey_pipeline(
                    vendor_csv=Path("data/compounds/vendor_library.csv"),
                    repurposing_csv=Path("data/compounds/repurposing_library.csv"),
                    docking_csv=Path("data/docking/multi_target_scores.csv"),
                    degradation_csv=Path("data/assays/degradation_template.csv"),
                    permeability_csv=Path("data/assays/permeability_template.csv"),
                    cytotoxicity_csv=Path("data/assays/cytotoxicity_template.csv"),
                    pk_csv=Path("data/preclinical/pk_timecourse_template.csv"),
                    pd_csv=Path("data/preclinical/pd_metabolic_template.csv"),
                    out_dir=out_dir,
                )
                # manifest.json embeds the output path, which differs between dirs.
                return {
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in out_dir.iterdir()
                    if p.is_file() and p.name != "manifest.json"
                }

            first = run_once(Path(tmp) / "run_a")
            second = run_once(Path(tmp) / "run_b")
            self.assertEqual(set(first), set(second))
            differing = [name for name in first if first[name] != second[name]]
            self.assertEqual(differing, [], f"non-deterministic artifacts: {differing}")


class CandidateExposureBridgeTest(unittest.TestCase):
    def _candidate(self, cid: str, kp, f_pct) -> PreclinicalCandidate:
        return PreclinicalCandidate(
            compound_id=cid,
            status="ADVANCED_LEAD",
            oral_bioavailability_pct=f_pct,
            brain_plasma_ratio=kp,
            plasma_cmax_ng_ml=None,
            plasma_auc_h_ng_ml=None,
            brain_auc_h_ng_ml=None,
            half_life_h=None,
            wfa_pnn_reduction_pct=None,
            itt_glucose_auc_reduction_pct=None,
            findings="",
        )

    def test_unmeasured_exposures_are_not_bridged_as_zero(self) -> None:
        # An unmeasured brain Kp (None) must stay absent so SAR reports it as
        # UNMEASURED. `or 0.0` fabricated a measured "0.00 no brain uptake",
        # violating G-001/ADR-010 (never invent evidence).
        candidates = [self._candidate("CMP-X", None, None)]
        self.assertEqual(_candidate_exposures(candidates), {})

    def test_measured_zero_exposures_are_preserved(self) -> None:
        # A measured zero is real evidence and must be kept, not dropped.
        candidates = [self._candidate("CMP-Z", 0.0, 0.0)]
        self.assertEqual(
            _candidate_exposures(candidates),
            {"CMP-Z": {"brain_kp": 0.0, "oral_f": 0.0}},
        )


if __name__ == "__main__":
    unittest.main()
