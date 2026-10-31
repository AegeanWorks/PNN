import tempfile
import unittest
from pathlib import Path

from pnn.sar import (
    CompoundSar,
    SarError,
    ScaffoldCluster,
    analyze_sar,
    extract_scaffold,
    load_compounds,
    main,
    parse_smiles_graph,
    render_sar_report,
)


class TestSarSmilesScaffolding(unittest.TestCase):
    def test_parse_empty_smiles_raises(self) -> None:
        with self.assertRaises(SarError):
            parse_smiles_graph("")

    def test_parse_unmatched_parens_raises(self) -> None:
        with self.assertRaises(SarError):
            parse_smiles_graph("CC(CC")

    def test_acyclic_molecule_scaffold(self) -> None:
        core, name, r_groups = extract_scaffold("CCCCCC(=O)O")
        self.assertEqual(name, "Acyclic-Scaffold")
        self.assertEqual(r_groups, "None")

    def test_monocyclic_benzene_scaffold(self) -> None:
        core, name, r_groups = extract_scaffold("Cc1ccccc1")
        self.assertIn("Phenyl", name)
        self.assertEqual(r_groups, "C")

    def test_bemis_murcko_extraction_multiring(self) -> None:
        # PNN-101: CC1=C(C=CC(=C1)NC(=O)C2=CC=C(C=C2)CN3CCNCC3)C
        smi = "CC1=C(C=CC(=C1)NC(=O)C2=CC=C(C=C2)CN3CCNCC3)C"
        core, name, r_groups = extract_scaffold(smi)
        self.assertIn("Piperazine", name)
        self.assertIn("Carboxamide", name)
        self.assertTrue(len(core) > 0)
        # Substituents should contain the two methyls
        self.assertIn("C", r_groups)


class TestSarAnalysisAndClustering(unittest.TestCase):
    def test_load_compounds_missing_file_raises(self) -> None:
        with self.assertRaises(SarError):
            load_compounds(Path("non_existent_compounds.csv"))

    def test_analyze_sar_scaffold_clustering(self) -> None:
        records = [
            {"compound_id": "PNN-01", "name": "Analog-A", "smiles": "CC1=CC=CC=C1NC(=O)C2=CC=CC=C2", "mw": "211.3", "logp": "2.1", "tpsa": "29.1", "pnn_ec50_nm": "35.0", "brain_plasma_ratio": "0.60"},
            {"compound_id": "PNN-02", "name": "Analog-B", "smiles": "FC1=CC=CC=C1NC(=O)C2=CC=CC=C2", "mw": "215.2", "logp": "2.3", "tpsa": "29.1", "pnn_ec50_nm": "45.0", "brain_plasma_ratio": "0.55"},
            {"compound_id": "PNN-03", "name": "Peripheral-1", "smiles": "c1ccccc1S(=O)(=O)Nc2ccccc2", "mw": "249.3", "logp": "1.8", "tpsa": "54.5", "pnn_ec50_nm": "850.0", "brain_plasma_ratio": "0.10"},
        ]
        dock_scores = {"PNN-01": -8.4, "PNN-02": -8.2, "PNN-03": -6.5}
        assays = {"PNN-01": 35.0, "PNN-02": 45.0, "PNN-03": 850.0}

        compounds, clusters = analyze_sar(records, docking_scores=dock_scores, assay_potencies=assays)
        self.assertEqual(len(compounds), 3)
        # PNN-01 and PNN-02 share the same biphenyl carboxamide core
        self.assertEqual(len(clusters), 2)

        # The priority lead cluster should have 2 compounds and best EC50 = 35.0
        lead_cluster = next(c for c in clusters if c.compound_count == 2)
        self.assertEqual(lead_cluster.best_compound, "PNN-01")
        self.assertEqual(lead_cluster.best_ec50_nm, 35.0)
        self.assertEqual(lead_cluster.sar_category, "PRIORITY_LEAD_SERIES")

    def test_render_sar_report(self) -> None:
        compounds = [
            CompoundSar(
                compound_id="PNN-101",
                name="HypoSprint-1",
                smiles="CC1=CC=CC=C1",
                scaffold_id="SCF-01",
                scaffold_name="Phenyl",
                core_smiles="c1ccccc1",
                r_groups="C",
                mw=92.1,
                logp=2.1,
                tpsa=0.0,
                vina_score=-8.5,
                ec50_nm=24.0,
                brain_kp=0.65,
                oral_f_pct=42.0,
                activity_tier="HIGH_POTENCY (<50 nM)",
            )
        ]
        clusters = [
            ScaffoldCluster(
                scaffold_id="SCF-01",
                scaffold_name="Phenyl",
                core_smiles="c1ccccc1",
                compound_count=1,
                compound_ids="PNN-101",
                best_compound="PNN-101",
                best_ec50_nm=24.0,
                mean_ec50_nm=24.0,
                best_vina_docking=-8.5,
                mean_brain_kp=0.65,
                sar_category="PRIORITY_LEAD_SERIES",
                key_sar_insights="High potency CNS penetrant series",
            )
        ]
        report = render_sar_report(compounds, clusters)
        self.assertIn("# PNN Structure-Activity Relationship (SAR)", report)
        self.assertIn("G-002", report)
        self.assertIn("SCF-01", report)
        self.assertIn("PRIORITY_LEAD_SERIES", report)


class TestSarCLI(unittest.TestCase):
    def test_cli_execution_success(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            comp_path = Path("data/compounds/vendor_library.csv")
            if not comp_path.is_file():
                return
            out_csv = Path(tmp_dir) / "clusters.csv"
            out_rep = Path(tmp_dir) / "report.md"
            out_json = Path(tmp_dir) / "summary.json"
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
