"""Integration test for the PNN Phase 2 golden reference run."""

from __future__ import annotations

import contextlib
import io
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from tools.golden_run import execute_golden_run


class GoldenRunTest(unittest.TestCase):
    def test_golden_run_execution(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                code = execute_golden_run(out_dir=Path(tmp_dir))
            self.assertEqual(code, 0)
            output = stdout.getvalue()
            self.assertIn("PNN Golden Run Verification PASSED", output)
            self.assertIn("Shortlisted Leads:", output)
            self.assertTrue((Path(tmp_dir) / "manifest.json").is_file())
            self.assertTrue((Path(tmp_dir) / "shortlist.csv").is_file())
            self.assertTrue((Path(tmp_dir) / "lead_report.md").is_file())


if __name__ == "__main__":
    unittest.main()
