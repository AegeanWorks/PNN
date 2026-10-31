import tempfile
import unittest
from pathlib import Path

from pnn.studio import (
    StudioCompound,
    StudioError,
    _get_field,
    _get_val,
    generate_circuit_dynamics_svg,
    generate_docking_heatmap_svg,
    generate_pk_curve_svg,
    generate_radar_svg,
    load_studio_data,
    main,
    render_studio_html,
)
from pnn.target import list_targets


class TestStudioRendering(unittest.TestCase):
    def setUp(self) -> None:
        self.compound = StudioCompound(
            compound_id="PNN-101",
            name="HypoSprint-1",
            mw=349.46,
            logp=2.65,
            tpsa=48.4,
            hbd=2,
            hba=4,
            rotb=4,
            ec50_nm=24.0,
            brain_kp=0.65,
            oral_f_pct=42.0,
            docking_scores={"PTPsigma": -8.5, "LAR": -7.2, "NgR1": -7.5},
            pk_points={"PO_plasma": [(0.0, 0.0), (1.0, 490.0), (4.0, 210.0), (24.0, 3.2)]},
        )

    def test_generate_radar_svg(self) -> None:
        svg = generate_radar_svg([self.compound])
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("<polygon", svg)
        self.assertIn("MW", svg)
        self.assertIn("logP", svg)
        self.assertIn("TPSA", svg)

    def test_generate_docking_heatmap_svg(self) -> None:
        svg = generate_docking_heatmap_svg([self.compound], targets=["PTPsigma", "LAR", "NgR1"])
        self.assertIn("<svg", svg)
        self.assertIn("PTPsigma", svg)
        self.assertIn("-8.5", svg)
        self.assertIn("<rect", svg)

    def test_generate_pk_curve_svg(self) -> None:
        svg = generate_pk_curve_svg([self.compound])
        self.assertIn("<svg", svg)
        self.assertIn("<polyline", svg)
        self.assertIn("Conc (ng/mL)", svg)
        self.assertIn("Time post-dose", svg)

    def test_render_studio_html(self) -> None:
        r_svg = generate_radar_svg([self.compound])
        h_svg = generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"])
        p_svg = generate_pk_curve_svg([self.compound])
        html_doc = render_studio_html([self.compound], r_svg, h_svg, p_svg)
        self.assertIn("<!DOCTYPE html>", html_doc)
        self.assertIn("PNN Discovery Studio", html_doc)
        self.assertIn("PNN-101", html_doc)
        self.assertIn("24.0 nM", html_doc)

    def test_generate_circuit_dynamics_svg(self) -> None:
        sim_points = [
            {"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0},
            {"time_h": 24.0, "pnn_density_pct": 55.0, "akt_signaling_pct": 78.0},
            {"time_h": 72.0, "pnn_density_pct": 42.0, "akt_signaling_pct": 92.0},
            {"time_h": 336.0, "pnn_density_pct": 60.0, "akt_signaling_pct": 85.0},
        ]
        svg = generate_circuit_dynamics_svg(sim_points)
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("PNN Density", svg)
        self.assertIn("pAKT", svg)
        self.assertIn("<polyline", svg)

    def test_render_studio_html_with_translational_panels(self) -> None:
        r_svg = generate_radar_svg([self.compound])
        h_svg = generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"])
        p_svg = generate_pk_curve_svg([self.compound])
        c_svg = generate_circuit_dynamics_svg([
            {"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0},
            {"time_h": 168.0, "pnn_density_pct": 45.0, "akt_signaling_pct": 88.0},
        ])
        dossier_data = [
            {
                "compound_id": "PNN-REP-01",
                "chemotype": "Trodusquemine-Analog",
                "human_dose_mg_day": "48.6",
                "human_t12_h": "17.0",
                "mrsd_mg": "58.38",
                "safety_margin": "6.0x",
                "evidence_status": "COMPLETE",
            }
        ]
        safety_data = [
            {
                "compound_id": "PNN-REP-01",
                "status": "SAFE_TO_PROGRESS",
                "herg_margin": "125.0x",
                "cyp3a4_inhibition_pct": "8.0",
                "genotoxicity": "NEGATIVE",
                "most_sensitive_species": "rat",
                "most_sensitive_noael_mg_kg": "60.0",
            }
        ]
        html_doc = render_studio_html(
            [self.compound],
            radar_svg=r_svg,
            heatmap_svg=h_svg,
            pk_svg=p_svg,
            circuit_svg=c_svg,
            dossier_data=dossier_data,
            safety_data=safety_data,
            targets_data=list_targets(),
        )
        self.assertIn("Investigational New Drug (IND)", html_doc)
        self.assertIn("Safety Pharmacology", html_doc)
        self.assertIn("Hypothalamic Circuit Dosing Dynamics", html_doc)
        self.assertIn("ADAMTS4", html_doc)
        self.assertIn("MMP9", html_doc)
        self.assertIn("Trodusquemine-Analog", html_doc)


class TestStudioDataLoader(unittest.TestCase):
    def test_missing_compounds_file_raises(self) -> None:
        with self.assertRaises(StudioError):
            load_studio_data(Path("non_existent.csv"))

    def test_load_studio_data_valid(self) -> None:
        project_root = Path(__file__).resolve().parent.parent
        comp_path = project_root / "data/compounds/vendor_library.csv"
        dock_path = project_root / "data/docking/multi_target_scores.csv"
        pk_path = project_root / "data/preclinical/pk_timecourse_template.csv"
        if not comp_path.is_file():
            return
        compounds, targets = load_studio_data(
            compounds_path=comp_path,
            docking_path=dock_path if dock_path.is_file() else None,
            pk_path=pk_path if pk_path.is_file() else None,
        )
        self.assertTrue(len(compounds) > 0)
        pnn101 = next((c for c in compounds if c.compound_id == "PNN-101"), None)
        self.assertIsNotNone(pnn101)


class TestStudioCLI(unittest.TestCase):
    def test_cli_execution_success(self) -> None:
        project_root = Path(__file__).resolve().parent.parent
        with tempfile.TemporaryDirectory() as tmp_dir:
            comp_path = project_root / "data/compounds/vendor_library.csv"
            if not comp_path.is_file():
                return
            out_html = Path(tmp_dir) / "studio.html"
            rc = main([
                "--compounds", str(comp_path),
                "--out", str(out_html),
            ])
            self.assertEqual(rc, 0)
            self.assertTrue(out_html.is_file())
            content = out_html.read_text(encoding="utf-8")
            self.assertIn("<!DOCTYPE html>", content)
            self.assertIn("Receptor & Aggrecanase Target Invariants", content)

    def test_cli_execution_with_translational_flags(self) -> None:
        project_root = Path(__file__).resolve().parent.parent
        comp_path = project_root / "data/compounds/vendor_library.csv"
        if not comp_path.is_file():
            return
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_html = Path(tmp_dir) / "studio_trans.html"
            sim_csv = Path(tmp_dir) / "sim.csv"
            sim_csv.write_text("time_h,pnn_density_pct,akt_signaling_pct\n0.0,100.0,45.0\n48.0,40.0,90.0\n", encoding="utf-8")
            dos_csv = Path(tmp_dir) / "dossier.csv"
            dos_csv.write_text("compound_id,name,animal_dose_mg_kg,human_total_daily_dose_mg,projected_human_half_life_h,phase1_mrsd_mg,therapeutic_index,evidence_status\nPNN-101,HypoSprint-1,10.0,48.6,17.0,58.4,6,COMPLETE\n", encoding="utf-8")
            saf_csv = Path(tmp_dir) / "safety.csv"
            saf_csv.write_text("compound_id,status,most_sensitive_species,noael_mg_kg,herg_margin,cyp3a4_inhibition_pct,genotoxicity\nPNN-101,SAFE_TO_PROGRESS,rat,50.0,51.4,12.5,NEGATIVE\n", encoding="utf-8")

            rc = main([
                "--compounds", str(comp_path),
                "--simulation", str(sim_csv),
                "--dossier", str(dos_csv),
                "--safety", str(saf_csv),
                "--out", str(out_html),
            ])
            self.assertEqual(rc, 0)
            self.assertTrue(out_html.is_file())
            content = out_html.read_text(encoding="utf-8")
            self.assertIn("Hypothalamic Circuit Dosing Dynamics", content)
            self.assertIn("Investigational New Drug (IND)", content)
            self.assertIn("Safety Pharmacology", content)
            self.assertIn("48.6 mg", content)
            self.assertIn("50.0 mg/kg", content)
            self.assertIn("51.4x", content)

    def test_cli_missing_translational_files_error(self) -> None:
        project_root = Path(__file__).resolve().parent.parent
        comp_path = project_root / "data/compounds/vendor_library.csv"
        if not comp_path.is_file():
            return
        rc = main([
            "--compounds", str(comp_path),
            "--simulation", "non_existent_sim.csv",
        ])
        self.assertEqual(rc, 1)


class TestFieldExtraction(unittest.TestCase):
    def test_get_val_with_mapping(self) -> None:
        d = {"key": "val", "empty": "", "none": None}
        self.assertEqual(_get_val(d, "key"), "val")
        self.assertEqual(_get_val(d, "empty"), "-")
        self.assertEqual(_get_val(d, "none"), "-")
        self.assertEqual(_get_val(d, "missing", "def"), "def")

    def test_get_val_with_object(self) -> None:
        class Dummy:
            key = "val"
            none = None
        obj = Dummy()
        self.assertEqual(_get_val(obj, "key"), "val")
        self.assertEqual(_get_val(obj, "none"), "-")
        self.assertEqual(_get_val(obj, "missing", "fallback"), "fallback")

    def test_get_field_aliases(self) -> None:
        d = {"b": "found"}
        self.assertEqual(_get_field(d, ["a", "b", "c"]), "found")
        self.assertEqual(_get_field(d, ["x", "y"], "default"), "default")
        self.assertEqual(_get_field(None, ["a", "b"], "default"), "default")


if __name__ == "__main__":
    unittest.main()
