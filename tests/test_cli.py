"""Integration tests for the `python -m pnn` CLI contract."""

import contextlib
import csv
import io
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from pnn.cli import main

HEADER = "compound_id,mw,logp,tpsa,hbd,hba,rotatable_bonds,pka,charge,logbb,brain_plasma_ratio,oral_bioavailability,cyp3a4_inhibition,herg_liability,pnn_ec50_nm\n"
IDEAL = "CMP-A,350,2.5,55,2,5,4,8,neutral,-0.2,0.5,0.35,10,5,5\n"
ACCEPTABLE = "CMP-B,220,3.3,65,3,7,7,7.5,neutral,-0.8,0.2,0.16,40,20,500\n"
REJECTED = "CMP-C,350,2.5,55,2,5,4,8,neutral,-0.2,0.5,0.35,10,50,5\n"


class CliTest(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.input = self.root / "compounds.csv"
        self.input.write_text(HEADER + IDEAL + ACCEPTABLE + REJECTED, encoding="utf-8")

    def run_main(self, *argv: str):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(io.StringIO()):
            code = main(list(argv))
        return code, stderr.getvalue()

    def read_rows(self, path: Path):
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    def test_writes_ranked_csv_and_json_summary(self) -> None:
        out = self.root / "reports" / "ranked.csv"
        report = self.root / "reports" / "summary.json"
        code, _ = self.run_main(str(self.input), "--out", str(out), "--json", str(report))
        self.assertEqual(code, 0)
        rows = self.read_rows(out)
        self.assertEqual([row["compound_id"] for row in rows], ["CMP-A", "CMP-B", "CMP-C"])
        self.assertEqual([row["tier"] for row in rows], ["A", "B", "Reject"])
        self.assertEqual(rows[0]["charge"], "neutral")
        self.assertEqual(rows[0]["pka"], "8.0")
        self.assertIn("smiles", rows[0])
        data = json.loads(report.read_text(encoding="utf-8"))
        self.assertEqual(data["summary"]["by_tier"], {"A": 1, "B": 1, "C": 0, "Reject": 1})
        self.assertEqual(data["emitted"], 3)

    def test_min_tier_filters_lower_tiers_and_rejects(self) -> None:
        out = self.root / "ranked.csv"
        code, _ = self.run_main(str(self.input), "--out", str(out), "--min-tier", "A")
        self.assertEqual(code, 0)
        rows = self.read_rows(out)
        self.assertEqual([row["compound_id"] for row in rows], ["CMP-A"])

    def test_invalid_input_returns_one_with_a_remediation_triad(self) -> None:
        broken = self.root / "broken.csv"
        broken.write_text("compound_id,mw\nCMP-X,350\n", encoding="utf-8")
        code, stderr = self.run_main(str(broken), "--out", str(self.root / "out.csv"))
        self.assertEqual(code, 1)
        self.assertIn("What:", stderr)
        self.assertIn("Why:", stderr)
        self.assertIn("Fix:", stderr)

    def test_missing_input_returns_one(self) -> None:
        code, stderr = self.run_main(str(self.root / "absent.csv"), "--out", str(self.root / "out.csv"))
        self.assertEqual(code, 1)
        self.assertIn("not found", stderr)

    def test_help_exits_zero_and_unknown_flag_exits_two(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as help_exit:
                main(["--help"])
            with self.assertRaises(SystemExit) as flag_exit:
                main([str(self.input), "--unknown"])
        self.assertEqual(help_exit.exception.code, 0)
        self.assertEqual(flag_exit.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
