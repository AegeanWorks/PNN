"""Unit tests for docking ingestion and consensus normalization."""

import unittest

from pnn.docking import (
    DockingError,
    DockingScore,
    aggregate_best,
    filter_scores,
    normalize_scores,
    parse_docking_scores,
)


class DockingIngestionTest(unittest.TestCase):
    def test_aliases_units_and_optional_columns(self) -> None:
        scores = parse_docking_scores(
            [
                {
                    "Compound ID": "CMP-1",
                    "Target": "PTPsigma",
                    "Vina affinity (kcal/mol)": "-8.1",
                    "Glide Score": "-7.2",
                    "MM-GBSA dG": "-40.0",
                    "RMSD (A)": "1.2",
                    "H-bond count": "3",
                }
            ]
        )
        self.assertEqual(scores[0].compound_id, "CMP-1")
        self.assertEqual(scores[0].target, "PTPsigma")
        self.assertEqual(scores[0].vina_affinity, -8.1)
        self.assertEqual(scores[0].glide_score, -7.2)
        self.assertEqual(scores[0].mmgbsa_dg, -40.0)
        self.assertEqual(scores[0].rmsd, 1.2)
        self.assertEqual(scores[0].hbond_count, 3.0)

    def test_row_needs_a_score(self) -> None:
        with self.assertRaisesRegex(DockingError, "at least one"):
            parse_docking_scores([{"compound_id": "X", "target": "T"}])

    def test_row_needs_target(self) -> None:
        with self.assertRaisesRegex(DockingError, "missing required field"):
            parse_docking_scores([{"compound_id": "X", "vina_affinity": "-8"}])

    def test_duplicate_compound_target_is_refused(self) -> None:
        rows = [
            {"compound_id": "X", "target": "T", "vina_affinity": "-8"},
            {"compound_id": "X", "target": "T", "vina_affinity": "-7"},
        ]
        with self.assertRaisesRegex(DockingError, "duplicate docking result"):
            parse_docking_scores(rows)

    def test_non_finite_score_is_refused(self) -> None:
        with self.assertRaisesRegex(DockingError, "must be finite"):
            parse_docking_scores([{"compound_id": "X", "target": "T", "vina_affinity": "inf"}])

    def test_negative_rmsd_is_refused(self) -> None:
        with self.assertRaisesRegex(DockingError, "must not be negative"):
            parse_docking_scores(
                [{"compound_id": "X", "target": "T", "vina_affinity": "-8", "rmsd": "-1"}]
            )


class ConsensusTest(unittest.TestCase):
    def test_vina_only_is_normalized_to_zero_and_one(self) -> None:
        ranks = normalize_scores(
            [
                DockingScore("A", "T", vina_affinity=-10.0),
                DockingScore("B", "T", vina_affinity=-5.0),
            ]
        )
        self.assertEqual([rank.consensus for rank in ranks], [1.0, 0.0])
        self.assertEqual([rank.coverage for rank in ranks], [0.5, 0.5])

    def test_multiple_methods_are_weighted_and_coverage_is_reported(self) -> None:
        ranks = normalize_scores(
            [
                DockingScore("A", "T", vina_affinity=-10.0, glide_score=-8.0, mmgbsa_dg=-50.0),
                DockingScore("B", "T", vina_affinity=-5.0, glide_score=-4.0, mmgbsa_dg=-20.0),
            ]
        )
        self.assertAlmostEqual(ranks[0].consensus, 1.0)
        self.assertAlmostEqual(ranks[1].consensus, 0.0)
        self.assertAlmostEqual(ranks[0].coverage, 1.0)

    def test_best_target_is_kept_and_targets_are_counted(self) -> None:
        ranks = normalize_scores(
            [
                DockingScore("A", "PTPsigma", vina_affinity=-10.0),
                DockingScore("A", "LAR", vina_affinity=-6.0),
                DockingScore("B", "PTPsigma", vina_affinity=-8.0),
            ]
        )
        best = aggregate_best(ranks)
        self.assertEqual(best["A"].best.score.target, "PTPsigma")
        self.assertEqual(best["A"].targets_docked, 2)
        self.assertEqual(best["B"].targets_docked, 1)

    def test_filters_drop_missing_vina_and_limit_breakers(self) -> None:
        scores = [
            DockingScore("A", "T", vina_affinity=-9.0, rmsd=1.0),
            DockingScore("B", "T", vina_affinity=-4.0),
            DockingScore("C", "T", glide_score=-6.0),
            DockingScore("D", "T", vina_affinity=-8.0, rmsd=3.0),
        ]
        kept, dropped = filter_scores(scores, max_rmsd=2.0, min_vina=-7.0)
        self.assertEqual([score.compound_id for score in kept], ["A"])
        self.assertEqual(dropped, 3)


if __name__ == "__main__":
    unittest.main()
