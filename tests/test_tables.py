"""Tests for the shared CSV row contracts used by every report writer."""

from __future__ import annotations

import unittest

from pnn.compounds import Compound
from pnn.curation import prepare_candidates
from pnn.screen import evaluate_compound
from pnn.tables import CURATED_FIELDS, SCREEN_FIELDS, _optional, curated_row, screen_row


class TestTableContracts(unittest.TestCase):
    def test_row_keys_match_declared_field_tuples(self) -> None:
        # A drift between a *_row() dict and its FIELDS tuple fails at runtime
        # inside csv.DictWriter (or silently writes empty columns).
        compound = Compound("X", 350.0, 2.5, 55.0, 2.0, 5.0, 4.0, smiles="CCO", name="Xname")
        screen = evaluate_compound(compound)
        candidate = prepare_candidates([("src", [compound])])[0]
        for label, fields, row in (
            ("curated", CURATED_FIELDS, curated_row(candidate, 0.1)),
            ("screen", SCREEN_FIELDS, screen_row(screen)),
        ):
            self.assertEqual(set(fields) - set(row), set(), f"{label} row missing fields")
            self.assertEqual(set(row) - set(fields), set(), f"{label} row has extra fields")
            self.assertEqual(list(row), list(fields), f"{label} row key order differs")

    def test_optional_maps_none_to_empty_string(self) -> None:
        self.assertEqual(_optional(None), "")
        self.assertEqual(_optional(0.0), 0.0)
        self.assertEqual(_optional(42.5), 42.5)


if __name__ == "__main__":
    unittest.main()
