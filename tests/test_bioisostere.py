"""Unit tests for bioisosteric replacement and metabolic soft-spot optimizer."""

import tempfile
import unittest
from pathlib import Path

from pnn.bioisostere import (
    BIOISOSTERE_RULES,
    BioisostereError,
    _ring_label_occurrences,
    _ring_labels_balanced,
    generate_bioisosteres,
    main,
    render_bioisostere_report,
    run_bioisostere_pipeline,
)


class TestBioisostere(unittest.TestCase):
    def test_generate_bioisosteres(self) -> None:
        smiles = "c1ccccc1C(=O)NC1CCOCC1"
        analogs = generate_bioisosteres(
            parent_id="TEST-01",
            smiles=smiles,
            mw=250.0,
            logp=2.1,
            tpsa=45.0,
            hbd=1,
            hba=3,
            rotb=3,
            pka=7.5,
        )
        self.assertGreater(len(analogs), 0)
        self.assertEqual(analogs[0].parent_id, "TEST-01")
        for analog in analogs:
            self.assertTrue(analog.analog_id.startswith("TEST-01-BIO"))
            self.assertTrue(_ring_labels_balanced(analog.analog_smiles), analog.analog_smiles)
            self.assertTrue(any(rule.rule_id == analog.rule_id for rule in BIOISOSTERE_RULES))

        row = analogs[0].as_row()
        self.assertIn("analog_smiles", row)
        self.assertIn("cns_mpo_score", row)
        self.assertIn("rule_id", row)

    def test_ring_label_parser_ignores_isotope_brackets(self) -> None:
        counts = _ring_label_occurrences("C([2H])([2H])N1CCNCC1")
        self.assertEqual(counts, {"1": 2})

    def test_replacement_does_not_split_ring_context(self) -> None:
        # The amide carbon in PNN-101 is the aromatic ring atom C2; replacing it would corrupt the ring.
        smiles = "CC1=C(C=CC(=C1)NC(=O)C2=CC=C(C=C2)CN3CCNCC3)C"
        analogs = generate_bioisosteres(
            parent_id="PNN-101",
            smiles=smiles,
            mw=349.46,
            logp=2.65,
            tpsa=48.4,
            hbd=2,
            hba=4,
            rotb=4,
            pka=8.2,
        )
        rule_ids = {a.rule_id for a in analogs}
        self.assertNotIn("BIO-03", rule_ids)
        for analog in analogs:
            self.assertNotIn("C1CC12", analog.analog_smiles)
            self.assertTrue(_ring_labels_balanced(analog.analog_smiles), analog.analog_smiles)

    def test_property_deltas_match_rule_table(self) -> None:
        analogs = generate_bioisosteres(
            parent_id="P",
            smiles="c1ccccc1C(=O)OC",
            mw=200.0,
            logp=2.0,
            tpsa=40.0,
            hbd=0,
            hba=3,
            rotb=4,
        )
        by_rule = {a.rule_id: a for a in analogs}
        self.assertIn("BIO-06", by_rule)
        ester = by_rule["BIO-06"]
        rule = next(r for r in BIOISOSTERE_RULES if r.rule_id == "BIO-06")
        self.assertAlmostEqual(ester.mw, 200.0 + rule.delta_mw)
        self.assertAlmostEqual(ester.tpsa, 40.0 + rule.delta_tpsa)

    def test_cyclopropyl_delta_is_not_a_rotatable_increase(self) -> None:
        rule = next(r for r in BIOISOSTERE_RULES if r.rule_id == "BIO-03")
        self.assertLess(rule.delta_rotb, 0)
        self.assertEqual(rule.delta_hba, 0)
        self.assertAlmostEqual(rule.delta_mw, 40.0)

    def test_run_bioisostere_pipeline_missing_file_raises(self) -> None:
        with self.assertRaises(BioisostereError):
            run_bioisostere_pipeline(Path("non_existent_library.csv"))

    def test_run_bioisostere_pipeline_vendor_lib(self) -> None:
        vendor_csv = Path("data/compounds/vendor_library.csv")
        if not vendor_csv.is_file():
            return

        analogs = run_bioisostere_pipeline(vendor_csv, target_parent_id="PNN-101")
        self.assertGreater(len(analogs), 0)
        for a in analogs:
            self.assertEqual(a.parent_id, "PNN-101")
            self.assertGreater(a.mw, 0.0)
            self.assertTrue(_ring_labels_balanced(a.analog_smiles), a.analog_smiles)

    def test_render_bioisostere_report(self) -> None:
        analogs = generate_bioisosteres(
            parent_id="PNN-101",
            smiles="c1ccccc1C(=O)NC1CCOCC1",
            mw=320.0,
            logp=2.4,
            tpsa=50.0,
            hbd=1,
            hba=3,
            rotb=4,
        )
        report = render_bioisostere_report(analogs, "PNN-101")
        self.assertIn("Bioisosteric Replacement", report)
        self.assertIn("PNN-101", report)
        self.assertIn("Method Limits", report)
        self.assertIn("Synthetic Priority Recommendations", report)

    def test_bioisostere_cli(self) -> None:
        vendor_csv = Path("data/compounds/vendor_library.csv")
        if not vendor_csv.is_file():
            return

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_csv = Path(tmp_dir) / "analogs.csv"
            out_md = Path(tmp_dir) / "report.md"
            out_json = Path(tmp_dir) / "summary.json"

            exit_code = main([
                "--compounds", str(vendor_csv),
                "--parent", "PNN-101",
                "--out", str(out_csv),
                "--report", str(out_md),
                "--json", str(out_json),
            ])
            self.assertEqual(exit_code, 0)
            self.assertTrue(out_csv.is_file())
            self.assertTrue(out_md.is_file())
            self.assertTrue(out_json.is_file())


if __name__ == "__main__":
    unittest.main()
