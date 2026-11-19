import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
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

    def test_heatmap_svg_escapes_user_names_and_stays_valid_xml(self) -> None:
        # Compound IDs and target names come from user CSVs; SVG is XML, so
        # raw '&'/'<' would produce malformed output. They must be escaped.
        special = StudioCompound(
            compound_id="A&B<C>",
            name="x",
            mw=300.0,
            logp=2.5,
            tpsa=50.0,
            hbd=2,
            hba=4,
            rotb=3,
            ec50_nm=10.0,
            brain_kp=None,
            oral_f_pct=None,
            docking_scores={"T&X": -9.0},
            pk_points={},
        )
        svg = generate_docking_heatmap_svg([special], targets=["T&X"])
        # Well-formed XML is the real contract; raw '&' would raise ParseError.
        ET.fromstring(svg)
        self.assertIn("A&amp;B&lt;C&gt;", svg)
        self.assertIn("T&amp;X", svg)

    def test_generate_pk_curve_svg(self) -> None:
        svg = generate_pk_curve_svg([self.compound])
        self.assertIn("<svg", svg)
        self.assertIn("<polyline", svg)
        self.assertIn("Conc (ng/mL)", svg)
        self.assertIn("Time post-dose", svg)

    def test_pk_curve_keeps_long_trajectories_on_canvas(self) -> None:
        # A hardcoded 24h time axis drew PK samples beyond 24h off-canvas,
        # where they are invisible. The axis must expand to fit the data.
        pts = {
            "IV_plasma": [
                (0.0, 500.0), (12.0, 200.0), (24.0, 80.0), (36.0, 30.0), (48.0, 10.0)
            ]
        }
        c = StudioCompound("C", "n", 300, 2.5, 50, 2, 4, 3, 10.0, None, None, {}, pts)
        svg = generate_pk_curve_svg([c])
        polys = re.findall(r'<polyline[^>]*points="([^"]+)"', svg)
        xs = [float(pt.split(",")[0]) for p in polys for pt in p.split() if "," in pt]
        self.assertTrue(xs)
        self.assertLessEqual(max(xs), 540 - 30 + 0.5)  # width - padding_right

    def test_render_studio_html(self) -> None:
        r_svg = generate_radar_svg([self.compound])
        h_svg = generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"])
        p_svg = generate_pk_curve_svg([self.compound])
        html_doc = render_studio_html([self.compound], r_svg, h_svg, p_svg)
        self.assertIn("<!DOCTYPE html>", html_doc)
        self.assertIn("PNN Discovery Studio", html_doc)
        self.assertIn("PNN-101", html_doc)
        self.assertIn("24.0 nM", html_doc)

    def test_pk_heading_names_the_compound_actually_plotted(self) -> None:
        other = StudioCompound(
            compound_id="TEST-PK-9",
            name="Alternative-Lead",
            mw=320.0,
            logp=2.4,
            tpsa=55.0,
            hbd=1,
            hba=4,
            rotb=4,
            ec50_nm=30.0,
            brain_kp=0.5,
            oral_f_pct=35.0,
            docking_scores={},
            pk_points={"IV_plasma": [(0.0, 0.0), (1.0, 400.0)]},
        )
        html_doc = render_studio_html(
            [other],
            generate_radar_svg([other]),
            generate_docking_heatmap_svg([other], targets=[]),
            generate_pk_curve_svg([other]),
        )
        self.assertIn("PK Profile: TEST-PK-9", html_doc)
        self.assertNotIn("PK Profile: PNN-101", html_doc)

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

    def test_circuit_axis_uses_hours_when_ticks_are_not_day_aligned(self) -> None:
        # Last point 84h -> ticks at 0, 12, 24, 36, 48, 60, 72, 84h.
        # Rounding t/24 to whole days mislabels a 12h tick as "0d" and a 36h tick as "2d".
        sim_points = [
            {"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0},
            {"time_h": 84.0, "pnn_density_pct": 60.0, "akt_signaling_pct": 85.0},
        ]
        svg = generate_circuit_dynamics_svg(sim_points)
        self.assertIn(">12h</text>", svg)
        self.assertIn(">36h</text>", svg)
        self.assertNotIn(">1d</text>", svg)
        self.assertIn("Simulation Duration (Hours post-initiation)", svg)

    def test_circuit_axis_uses_days_when_ticks_are_day_aligned(self) -> None:
        sim_points = [
            {"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0},
            {"time_h": 168.0, "pnn_density_pct": 60.0, "akt_signaling_pct": 85.0},
        ]
        svg = generate_circuit_dynamics_svg(sim_points)
        self.assertIn(">1d</text>", svg)
        self.assertIn(">7d</text>", svg)
        self.assertNotIn(">12h</text>", svg)
        self.assertIn("Simulation Duration (Days post-initiation)", svg)

    def test_unknown_genotoxicity_is_not_a_failure_badge(self) -> None:
        # Reconciles reforge finding review-unknown-genotoxicity-fail: an
        # unmeasured genotoxicity result must be a caution, never a fail.
        html_doc = render_studio_html(
            [self.compound],
            radar_svg=generate_radar_svg([self.compound]),
            heatmap_svg=generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"]),
            pk_svg=generate_pk_curve_svg([self.compound]),
            safety_data=[{"compound_id": "PNN-101", "status": "CAUTION", "genotoxicity": "UNKNOWN"}],
        )
        self.assertIn('<span class="badge-caution">UNKNOWN</span>', html_doc)
        self.assertNotIn('<span class="badge-fail">UNKNOWN</span>', html_doc)

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


class TestStudioAccessibility(unittest.TestCase):
    """The dashboard is read by assistive technology, not only by eye."""

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
            docking_scores={"PTPsigma": -8.5},
            pk_points={
                "PO_plasma": [(0.0, 0.0), (1.0, 490.0)],
                "PO_brain": [(0.0, 0.0), (1.0, 120.0)],
            },
        )

    def assert_labelled(self, svg: str, stem: str) -> None:
        self.assertIn('role="img"', svg)
        self.assertIn(f'aria-labelledby="{stem}-title {stem}-desc"', svg)
        self.assertIn(f'<title id="{stem}-title">', svg)
        self.assertIn(f'<desc id="{stem}-desc">', svg)

    def test_svg_charts_expose_accessible_names(self) -> None:
        self.assert_labelled(generate_radar_svg([self.compound]), "radar")
        self.assert_labelled(
            generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"]), "heatmap"
        )
        self.assert_labelled(generate_pk_curve_svg([self.compound]), "pk")
        self.assert_labelled(
            generate_circuit_dynamics_svg(
                [{"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0}]
            ),
            "circuit",
        )

    def test_empty_circuit_charts_keep_accessible_names(self) -> None:
        self.assert_labelled(generate_circuit_dynamics_svg([]), "circuit")
        self.assert_labelled(generate_circuit_dynamics_svg([{"time_h": "bad"}]), "circuit")

    def test_tables_have_captions_and_column_scopes(self) -> None:
        html_doc = render_studio_html(
            [self.compound],
            radar_svg=generate_radar_svg([self.compound]),
            heatmap_svg=generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"]),
            pk_svg=generate_pk_curve_svg([self.compound]),
            dossier_data=[{"compound_id": "PNN-REP-01", "chemotype": "Analog"}],
            safety_data=[{"compound_id": "PNN-REP-01", "status": "SAFE_TO_PROGRESS"}],
            targets_data=list_targets(),
        )
        self.assertEqual(html_doc.count('<caption class="sr-only">'), 3)
        self.assertNotIn("<th>", html_doc)
        self.assertGreaterEqual(html_doc.count('<th scope="col">'), 20)
        self.assertIn("IND dose projection", html_doc)
        self.assertIn("Safety pharmacology and repeat-dose toxicology", html_doc)
        self.assertIn("Registered receptor and aggrecanase target", html_doc)

    def test_charts_expose_data_tables(self) -> None:
        html_doc = render_studio_html(
            [self.compound],
            radar_svg=generate_radar_svg([self.compound]),
            heatmap_svg=generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"]),
            pk_svg=generate_pk_curve_svg([self.compound]),
        )
        self.assertIn("Radar chart data: CNS physicochemical properties", html_doc)
        self.assertIn("Docking heatmap data: predicted binding affinity", html_doc)
        self.assertIn("Pharmacokinetic chart data: concentration over time", html_doc)
        self.assertIn("<td>349.5</td>", html_doc)
        self.assertIn("<td>-8.5</td>", html_doc)
        self.assertIn("<td>PO plasma</td>", html_doc)
        self.assertIn('<th scope="row">PNN-101</th>', html_doc)

    def test_circuit_chart_exposes_data_table(self) -> None:
        points = [
            {"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0},
            {"time_h": 168.0, "pnn_density_pct": 45.0, "akt_signaling_pct": 88.0},
        ]
        html_doc = render_studio_html(
            [self.compound],
            radar_svg=generate_radar_svg([self.compound]),
            heatmap_svg=generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"]),
            pk_svg=generate_pk_curve_svg([self.compound]),
            circuit_svg=generate_circuit_dynamics_svg(points),
            circuit_points=points,
        )
        self.assertIn(
            "Circuit dynamics data: PNN density and pAKT signaling over time", html_doc
        )
        self.assertIn('<th scope="row">0</th>', html_doc)
        self.assertIn('<th scope="row">168</th>', html_doc)

    def test_chart_labels_use_white_background_contrast(self) -> None:
        """The SVGs sit on a white wrapper; every label clears 4.5:1 there."""
        radar = generate_radar_svg([self.compound])
        heatmap = generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"])
        pk = generate_pk_curve_svg([self.compound])
        circuit = generate_circuit_dynamics_svg(
            [{"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0}]
        )
        self.assertIn('fill="#475569"', radar)
        self.assertIn('fill="#334155"', heatmap)
        self.assertIn('fill="#1e293b"', heatmap)
        self.assertIn('fill="#475569"', pk)
        self.assertIn('stroke="#047857"', pk)
        self.assertIn('stroke="#b45309"', pk)
        self.assertIn('fill="#b91c1c"', circuit)

    def test_chart_series_use_dash_patterns_and_radar_legend(self) -> None:
        html_doc = render_studio_html(
            [self.compound],
            radar_svg=generate_radar_svg([self.compound]),
            heatmap_svg=generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"]),
            pk_svg=generate_pk_curve_svg([self.compound]),
            circuit_svg=generate_circuit_dynamics_svg(
                [{"time_h": 0.0, "pnn_density_pct": 100.0, "akt_signaling_pct": 45.0}]
            ),
        )
        self.assertIn('class="chart-legend"', html_doc)
        self.assertIn('class="legend-swatch"', html_doc)
        self.assertIn("Green dashed: PO Plasma", html_doc)
        self.assertIn("Cyan dashed: pAKT Signaling", html_doc)
        self.assertIn('stroke-dasharray="6,3"', html_doc)
        self.assertIn('stroke-dasharray="2,3"', html_doc)
        self.assertIn("Lead candidates", html_doc)

    def test_page_has_main_landmark_and_skip_link(self) -> None:
        html_doc = render_studio_html(
            [self.compound],
            radar_svg=generate_radar_svg([self.compound]),
            heatmap_svg=generate_docking_heatmap_svg([self.compound], targets=["PTPsigma"]),
            pk_svg=generate_pk_curve_svg([self.compound]),
        )
        self.assertIn('<a class="skip-link" href="#main-content">', html_doc)
        self.assertIn('<main id="main-content">', html_doc)
        self.assertLess(html_doc.find("<main"), html_doc.find("<footer>"))
        self.assertLess(html_doc.find("</main>"), html_doc.find("<footer>"))


class TestStudioDataLoader(unittest.TestCase):
    def test_missing_compounds_file_raises(self) -> None:
        with self.assertRaises(StudioError):
            load_studio_data(Path("non_existent.csv"))

    def test_missing_docking_stays_missing(self) -> None:
        """Absent docking evidence must not be synthesized (ADR-010/ADR-013)."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp_dir:
            comp_path = Path(tmp_dir) / "compounds.csv"
            comp_path.write_text(
                "compound_id,name,mw,logp,tpsa,hbd,hba,rotatable_bonds\n"
                "PNN-999,SynthCheck,350,2.5,55,1,4,3\n",
                encoding="utf-8",
            )
            compounds, _ = load_studio_data(compounds_path=comp_path)
            self.assertEqual(len(compounds), 1)
            self.assertEqual(compounds[0].docking_scores, {})

    def test_malformed_numeric_raises_studio_error(self) -> None:
        """Empty required physicochemical cells raise StudioError, not ValueError."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp_dir:
            comp_path = Path(tmp_dir) / "bad.csv"
            comp_path.write_text(
                "compound_id,name,mw,logp,tpsa,hbd,hba,rotatable_bonds\n"
                "PNN-1,Bad,,2.0,50,1,3,2\n",
                encoding="utf-8",
            )
            with self.assertRaises(StudioError):
                load_studio_data(compounds_path=comp_path)

    def test_short_row_none_cells_raise_studio_error(self) -> None:
        """Short CSV rows (None cells) raise StudioError instead of AttributeError."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp_dir:
            comp_path = Path(tmp_dir) / "short.csv"
            comp_path.write_text(
                "compound_id,name,mw,logp,tpsa,hbd,hba,rotatable_bonds\n"
                "PNN-2,Short,350\n",
                encoding="utf-8",
            )
            with self.assertRaises(StudioError):
                load_studio_data(compounds_path=comp_path)

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
        # The docking CSV uses the canonical 'vina_affinity' header; the studio
        # loader must recognize it or the heatmap shows N/A for every compound.
        if dock_path.is_file():
            self.assertTrue(
                any(c.docking_scores for c in compounds),
                "studio loaded zero docking scores from the reference docking CSV",
            )
            self.assertTrue(pnn101.docking_scores)  # type: ignore[union-attr]


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
            # The standalone CLI dashboard must include the same Wave Ledger
            # panel as the turnkey dashboard; it was silently missing.
            self.assertIn("Wave Ledger", content)

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
