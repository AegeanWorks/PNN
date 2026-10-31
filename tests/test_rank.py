"""Unit and CLI tests for the docking consensus shortlist."""

import contextlib
import csv
import io
import json
import runpy
import sys
import unittest
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from pnn.docking import DockingScore, load_screen_table
from pnn.rank import build_shortlist, main

DOCKING_HEADER = "compound_id,target,vina_affinity,glide_score,mmgbsa_dg,rmsd\n"
DOCKING_ROWS = (
    "CMP-A,PTPsigma,-10,-8,-50,1.0\n"
    "CMP-A,LAR,-6,-5,-20,1.5\n"
    "CMP-B,PTPsigma,-8,-4,-30,1.0\n"
    "CMP-C,PTPsigma,-6,-3,-10,3.5\n"
)
SCREEN_HEADER = "compound_id,name,smiles,tier,score\n"
SCREEN_ROWS = (
    "CMP-A,Alpha,CC, A,0.95\n"
    "CMP-B,Beta,CCC,B,0.75\n"
    "CMP-C,Gamma,CCCC,C,0.55\n"
)


def screen_table():
    return {
        "CMP-A": {"compound_id": "CMP-A", "tier": "A", "score": "0.95", "name": "Alpha", "smiles": "CC"},
        "CMP-B": {"compound_id": "CMP-B", "tier": "B", "score": "0.75", "name": "Beta", "smiles": "CCC"},
        "CMP-C": {"compound_id": "CMP-C", "tier": "C", "score": "0.55", "name": "Gamma", "smiles": "CCCC"},
    }


def scores():
    return [
        DockingScore("CMP-A", "PTPsigma", vina_affinity=-10.0, glide_score=-8.0, mmgbsa_dg=-50.0, rmsd=1.0),
        DockingScore("CMP-A", "LAR", vina_affinity=-6.0, glide_score=-5.0, mmgbsa_dg=-20.0, rmsd=1.5),
        DockingScore("CMP-B", "PTPsigma", vina_affinity=-8.0, glide_score=-4.0, mmgbsa_dg=-30.0, rmsd=1.0),
        DockingScore("CMP-C", "PTPsigma", vina_affinity=-6.0, glide_score=-3.0, mmgbsa_dg=-10.0, rmsd=3.5),
    ]


class ShortlistTest(unittest.TestCase):
    def test_tier_filter_and_consensus_order(self) -> None:
        rows, summary = build_shortlist(scores(), screen=screen_table(), min_tier="B")
        self.assertEqual([row.compound_id for row in rows], ["CMP-A", "CMP-B"])
        self.assertEqual([row.tier for row in rows], ["A", "B"])
        self.assertEqual(summary["emitted"], 2)
        self.assertEqual(summary["skipped_by_tier"], 1)

    def test_rmsd_filter_removes_a_target_row(self) -> None:
        rows, summary = build_shortlist(
            scores(), screen=screen_table(), min_tier="C", max_rmsd=2.0
        )
        self.assertEqual([row.compound_id for row in rows], ["CMP-A", "CMP-B"])
        self.assertEqual(summary["skipped_by_filters"], 1)

    def test_vina_filter_removes_weak_rows(self) -> None:
        rows, _ = build_shortlist(scores(), screen=screen_table(), min_tier="B", min_vina=-9.0)
        self.assertEqual([row.compound_id for row in rows], ["CMP-A"])

    def test_compound_missing_from_screen_is_skipped(self) -> None:
        screen = screen_table()
        screen.pop("CMP-B")
        rows, summary = build_shortlist(scores(), screen=screen, min_tier="B")
        self.assertEqual([row.compound_id for row in rows], ["CMP-A"])
        self.assertEqual(summary["skipped_not_in_screen"], 1)

    def test_docking_only_mode_has_no_tier(self) -> None:
        rows, summary = build_shortlist(scores())
        self.assertEqual([row.compound_id for row in rows], ["CMP-A", "CMP-B", "CMP-C"])
        self.assertEqual(rows[0].tier, "")
        self.assertEqual(summary["skipped_by_tier"], 0)


class RankCliTest(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.docking = self.root / "docking.csv"
        self.docking.write_text(DOCKING_HEADER + DOCKING_ROWS, encoding="utf-8")
        self.screen = self.root / "screen.csv"
        self.screen.write_text(SCREEN_HEADER + SCREEN_ROWS, encoding="utf-8")

    def run_main(self, *argv: str):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(io.StringIO()):
            code = main(list(argv))
        return code, stderr.getvalue()

    def read_rows(self, path: Path):
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    def test_cli_writes_shortlist_and_json(self) -> None:
        out = self.root / "reports" / "shortlist.csv"
        report = self.root / "reports" / "summary.json"
        code, _ = self.run_main(
            str(self.docking), "--screen", str(self.screen), "--out", str(out), "--json", str(report),
            "--max-rmsd", "2.0",
        )
        self.assertEqual(code, 0)
        rows = self.read_rows(out)
        self.assertEqual([row["compound_id"] for row in rows], ["CMP-A", "CMP-B"])
        self.assertEqual(rows[0]["tier"], "A")
        data = json.loads(report.read_text(encoding="utf-8"))
        self.assertEqual(data["emitted"], 2)
        self.assertEqual(data["weights"], {"vina_affinity": 0.5, "glide_score": 0.3, "mmgbsa_dg": 0.2})

    def test_cli_returns_one_when_filters_remove_everything(self) -> None:
        code, stderr = self.run_main(
            str(self.docking), "--screen", str(self.screen), "--out", str(self.root / "out.csv"),
            "--min-tier", "A", "--min-vina", "-20",
        )
        self.assertEqual(code, 1)
        self.assertIn("What:", stderr)

    def test_cli_returns_one_for_invalid_docking(self) -> None:
        broken = self.root / "broken.csv"
        broken.write_text("compound_id,target\nCMP-X,PTPsigma\n", encoding="utf-8")
        code, stderr = self.run_main(str(broken), "--out", str(self.root / "out.csv"))
        self.assertEqual(code, 1)
        self.assertIn("Fix:", stderr)

    def test_screen_table_loader_reads_tier_and_score(self) -> None:
        table = load_screen_table(self.screen)
        self.assertEqual(table["CMP-A"]["tier"], "A")
        self.assertEqual(table["CMP-B"]["score"], "0.75")

    def test_module_entry_point_runs_main(self) -> None:
        out = self.root / "module_entry.csv"
        argv = [
            "pnn.rank",
            str(self.docking),
            "--screen",
            str(self.screen),
            "--out",
            str(out),
        ]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                with self.assertRaises(SystemExit) as exit_code:
                    runpy.run_module("pnn.rank", run_name="__main__")
        self.assertEqual(exit_code.exception.code, 0)
        self.assertTrue(out.is_file())


if __name__ == "__main__":
    unittest.main()
