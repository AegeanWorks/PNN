"""Unit and CLI tests for compound-library curation."""

import contextlib
import csv
import io
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from pnn.compounds import Compound
from pnn.curate import main
from pnn.curation import (
    Candidate,
    CurationError,
    curate,
    deduplicate,
    prepare_candidates,
    property_bin,
)
from pnn.screen import evaluate_compound


def compound(compound_id: str, smiles: str = "", **overrides) -> Compound:
    values = {
        "compound_id": compound_id,
        "mw": 350.0,
        "logp": 2.5,
        "tpsa": 55.0,
        "hbd": 2.0,
        "hba": 5.0,
        "rotatable_bonds": 4.0,
    }
    values.update(overrides)
    return Compound(smiles=smiles, **values)


def candidate(compound_id: str, source: str = "lib.csv", smiles: str = "", **overrides) -> Candidate:
    item = compound(compound_id, smiles=smiles, **overrides)
    return Candidate(compound=item, source=source, result=evaluate_compound(item))


class DeduplicationTest(unittest.TestCase):
    def test_duplicate_id_keeps_the_better_candidate(self) -> None:
        unique, counts = deduplicate(
            [candidate("A", smiles="CC"), candidate("A", smiles="CC", mw=500.0, logp=6.0)]
        )
        self.assertEqual(len(unique), 1)
        self.assertEqual(unique[0].result.tier, "B")
        self.assertEqual(counts["duplicate_ids"], 1)
        self.assertEqual(counts["duplicate_smiles"], 0)

    def test_better_record_keeps_smiles_from_the_other_source(self) -> None:
        complete = candidate(
            "A",
            source="a.csv",
            pka=8.0,
            charge="neutral",
            logbb=-0.2,
            brain_plasma_ratio=0.5,
            oral_bioavailability=0.35,
            cyp3a4_inhibition=10.0,
            herg_liability=5.0,
            pnn_ec50_nm=5.0,
        )
        unique, _ = deduplicate([complete, candidate("A", source="b.csv", smiles="CC")])
        self.assertEqual(unique[0].result.tier, "A")
        self.assertEqual(unique[0].compound.smiles, "CC")

    def test_same_id_with_different_smiles_is_refused(self) -> None:
        with self.assertRaisesRegex(CurationError, "two different SMILES"):
            deduplicate([candidate("A", smiles="CC"), candidate("A", smiles="CCC")])

    def test_same_smiles_with_different_ids_is_collapsed(self) -> None:
        unique, counts = deduplicate([candidate("A", smiles="CC"), candidate("B", smiles="CC")])
        self.assertEqual(len(unique), 1)
        self.assertEqual(counts["duplicate_smiles"], 1)

    def test_blank_smiles_do_not_collapse_distinct_ids(self) -> None:
        unique, counts = deduplicate([candidate("A"), candidate("B")])
        self.assertEqual(len(unique), 2)
        self.assertEqual(counts["duplicate_smiles"], 0)

    def test_prepare_candidates_keeps_sources(self) -> None:
        prepared = prepare_candidates([("a.csv", [compound("A")]), ("b.csv", [compound("B")])])
        self.assertEqual([item.source for item in prepared], ["a.csv", "b.csv"])


class CurationTest(unittest.TestCase):
    def test_tier_filter_removes_rejected_candidates(self) -> None:
        selected, summary = curate([candidate("A"), candidate("R", herg_liability=60.0)])
        self.assertEqual([item.compound.compound_id for item in selected], ["A"])
        self.assertEqual(summary["skipped_by_tier"], 1)

    def test_per_bin_cap_keeps_one_best_per_property_cell(self) -> None:
        selected, summary = curate(
            [candidate("A"), candidate("B", mw=352.0), candidate("C", mw=210.0)],
            per_bin=1,
        )
        self.assertEqual(len(selected), 2)
        self.assertEqual(summary["bins"], 2)

    def test_max_compounds_truncates_and_reports(self) -> None:
        selected, summary = curate([candidate("A"), candidate("B", mw=210.0)], max_compounds=1)
        self.assertEqual(len(selected), 1)
        self.assertTrue(summary["truncated_by_max"])

    def test_property_bin_is_deterministic(self) -> None:
        self.assertEqual(property_bin(compound("A"), 0.1), property_bin(compound("A"), 0.1))

    def test_invalid_bin_width_is_refused(self) -> None:
        with self.assertRaisesRegex(CurationError, "bin width"):
            curate([candidate("A")], bin_width=0.0)


class CurateCliTest(unittest.TestCase):
    HEADER = "compound_id,name,smiles,mw,logp,tpsa,hbd,hba,rotatable_bonds,herg_liability\n"

    def setUp(self) -> None:
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.first = self.root / "first.csv"
        self.first.write_text(
            self.HEADER
            + "CMP-A,Alpha,CC,350,2.5,55,2,5,4,\n"
            + "CMP-A,Alpha duplicate,CC,350,2.5,55,2,5,4,\n",
            encoding="utf-8",
        )
        self.second = self.root / "second.csv"
        self.second.write_text(
            self.HEADER
            + "CMP-B,Beta,CC,350,2.5,55,2,5,4,\n"
            + "CMP-R,Rejected,CCC,350,2.5,55,2,5,4,60\n",
            encoding="utf-8",
        )

    def run_main(self, *argv: str):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(io.StringIO()):
            code = main(list(argv))
        return code, stderr.getvalue()

    def read_rows(self, path: Path):
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    def test_merges_deduplicates_and_filters(self) -> None:
        out = self.root / "curated.csv"
        report = self.root / "summary.json"
        code, _ = self.run_main(
            str(self.first), str(self.second), "--out", str(out), "--json", str(report)
        )
        self.assertEqual(code, 0)
        rows = self.read_rows(out)
        self.assertEqual([row["compound_id"] for row in rows], ["CMP-A"])
        data = json.loads(report.read_text(encoding="utf-8"))
        self.assertEqual(data["input_rows"], 4)
        self.assertEqual(data["duplicate_ids"], 1)
        self.assertEqual(data["duplicate_smiles"], 1)
        self.assertEqual(data["skipped_by_tier"], 1)

    def test_conflicting_id_returns_one(self) -> None:
        broken = self.root / "broken.csv"
        broken.write_text(
            self.HEADER + "CMP-A,Other,CCC,350,2.5,55,2,5,4,\n", encoding="utf-8"
        )
        code, stderr = self.run_main(str(self.first), str(broken), "--out", str(self.root / "out.csv"))
        self.assertEqual(code, 1)
        self.assertIn("What:", stderr)

    def test_no_eligible_compounds_returns_one(self) -> None:
        code, stderr = self.run_main(
            str(self.second), "--out", str(self.root / "out.csv"), "--min-tier", "A"
        )
        self.assertEqual(code, 1)
        self.assertIn("no compounds survived", stderr)


if __name__ == "__main__":
    unittest.main()
