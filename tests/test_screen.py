"""Unit tests for tiering, ranking and CSV ingestion."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from pnn.compounds import Compound, CompoundError, load_compounds, parse_compounds
from pnn.screen import evaluate_compound, screen_compounds, summarize


def ideal(compound_id: str = "CMP-001", **overrides) -> Compound:
    values = {
        "compound_id": compound_id,
        "mw": 350.0,
        "logp": 2.5,
        "tpsa": 55.0,
        "hbd": 2.0,
        "hba": 5.0,
        "rotatable_bonds": 4.0,
        "charge": "neutral",
        "pka": 8.0,
        "logbb": -0.2,
        "brain_plasma_ratio": 0.5,
        "oral_bioavailability": 0.35,
        "cyp3a4_inhibition": 10.0,
        "herg_liability": 5.0,
        "pnn_ec50_nm": 5.0,
    }
    values.update(overrides)
    return Compound(**values)


def acceptable(compound_id: str = "CMP-002", **overrides) -> Compound:
    values = {
        "compound_id": compound_id,
        "mw": 220.0,
        "logp": 3.3,
        "tpsa": 65.0,
        "hbd": 3.0,
        "hba": 7.0,
        "rotatable_bonds": 7.0,
        "charge": "neutral",
        "pka": 7.5,
        "logbb": -0.8,
        "brain_plasma_ratio": 0.2,
        "oral_bioavailability": 0.16,
        "cyp3a4_inhibition": 40.0,
        "herg_liability": 20.0,
        "pnn_ec50_nm": 500.0,
    }
    values.update(overrides)
    return Compound(**values)


class TierTest(unittest.TestCase):
    def test_complete_ideal_compound_is_tier_a(self) -> None:
        self.assertEqual(evaluate_compound(ideal()).tier, "A")

    def test_property_only_compound_cannot_be_tier_a(self) -> None:
        bare = Compound("CMP-003", 350.0, 2.5, 55.0, 2.0, 5.0, 4.0)
        result = evaluate_compound(bare)
        self.assertEqual(result.tier, "B")
        self.assertEqual(len(result.unknown_criteria()), 8)

    def test_acceptable_compound_is_tier_b(self) -> None:
        self.assertEqual(evaluate_compound(acceptable()).tier, "B")

    def test_two_core_failures_are_tier_c(self) -> None:
        self.assertEqual(evaluate_compound(ideal(mw=500.0, logp=6.0)).tier, "C")

    def test_three_core_failures_are_rejected(self) -> None:
        self.assertEqual(evaluate_compound(ideal(mw=500.0, logp=6.0, tpsa=120.0)).tier, "Reject")

    def test_hard_liability_rejects_even_an_optimal_profile(self) -> None:
        self.assertEqual(evaluate_compound(ideal(herg_liability=50.0)).tier, "Reject")


class RankingTest(unittest.TestCase):
    def test_ranking_orders_by_tier_then_score(self) -> None:
        results = screen_compounds(
            [
                ideal("CMP-C", mw=500.0, logp=6.0),
                ideal("CMP-B", herg_liability=50.0),
                acceptable("CMP-B2"),
                ideal("CMP-A"),
            ]
        )
        self.assertEqual(
            [result.compound.compound_id for result in results],
            ["CMP-A", "CMP-B2", "CMP-C", "CMP-B"],
        )
        summary = summarize(results)
        self.assertEqual(summary["total"], 4)
        self.assertEqual(summary["by_tier"], {"A": 1, "B": 1, "C": 1, "Reject": 1})


class CsvIngestionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)

    def write_csv(self, text: str) -> Path:
        path = Path(self.directory.name) / "compounds.csv"
        path.write_text(text, encoding="utf-8")
        return path

    def test_aliases_and_blank_optionals_are_supported(self) -> None:
        path = self.write_csv(
            "Compound ID,Molecular Weight,CLogP,TPSA,H Bond Donors,H Bond Acceptors,Rotors,Brain to Plasma\n"
            "CMP-1,350,2.5,55,2,5,4,0.5\n"
        )
        compound = load_compounds(path)[0]
        self.assertEqual(compound.compound_id, "CMP-1")
        self.assertEqual(compound.brain_plasma_ratio, 0.5)
        self.assertIsNone(compound.logbb)

    def test_duplicate_ids_are_refused(self) -> None:
        rows = [
            {"compound_id": "X", "mw": "350", "logp": "2", "tpsa": "50", "hbd": "2", "hba": "5", "rotatable_bonds": "4"},
            {"compound_id": "X", "mw": "350", "logp": "2", "tpsa": "50", "hbd": "2", "hba": "5", "rotatable_bonds": "4"},
        ]
        with self.assertRaisesRegex(CompoundError, "duplicate compound_id"):
            parse_compounds(rows)

    def test_bad_number_names_the_row(self) -> None:
        rows = [{"compound_id": "X", "mw": "heavy", "logp": "2", "tpsa": "50", "hbd": "2", "hba": "5", "rotatable_bonds": "4"}]
        with self.assertRaisesRegex(CompoundError, "row 2"):
            parse_compounds(rows)

    def test_missing_field_is_refused(self) -> None:
        rows = [{"compound_id": "X", "logp": "2", "tpsa": "50", "hbd": "2", "hba": "5", "rotatable_bonds": "4"}]
        with self.assertRaisesRegex(CompoundError, "missing required field"):
            parse_compounds(rows)

    def test_non_finite_number_is_refused(self) -> None:
        rows = [{"compound_id": "X", "mw": "nan", "logp": "2", "tpsa": "50", "hbd": "2", "hba": "5", "rotatable_bonds": "4"}]
        with self.assertRaisesRegex(CompoundError, "must be finite"):
            parse_compounds(rows)

    def test_negative_property_is_refused(self) -> None:
        rows = [{"compound_id": "X", "mw": "350", "logp": "2", "tpsa": "-50", "hbd": "2", "hba": "5", "rotatable_bonds": "4"}]
        with self.assertRaisesRegex(CompoundError, "must not be negative"):
            parse_compounds(rows)

    def test_fractional_count_is_refused(self) -> None:
        rows = [{"compound_id": "X", "mw": "350", "logp": "2", "tpsa": "50", "hbd": "2.5", "hba": "5", "rotatable_bonds": "4"}]
        with self.assertRaisesRegex(CompoundError, "whole number"):
            parse_compounds(rows)

    def test_missing_file_is_refused(self) -> None:
        with self.assertRaisesRegex(CompoundError, "not found"):
            load_compounds(Path(self.directory.name) / "absent.csv")

    def test_header_only_file_is_refused(self) -> None:
        path = self.write_csv("compound_id,mw,logp,tpsa,hbd,hba,rotatable_bonds\n")
        with self.assertRaisesRegex(CompoundError, "no compound rows"):
            load_compounds(path)

    def test_non_utf8_file_is_refused(self) -> None:
        path = Path(self.directory.name) / "binary.csv"
        path.write_bytes(b"\xff\xfe\x00not text\x00")
        with self.assertRaisesRegex(CompoundError, "valid UTF-8"):
            load_compounds(path)


if __name__ == "__main__":
    unittest.main()
