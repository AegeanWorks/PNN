"""Tests for the shared atomic CSV/JSON/text report writers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pnn.reports import write_csv, write_json, write_text


class TestAtomicWriters(unittest.TestCase):
    def test_writers_create_parent_and_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            write_text(root / "a" / "b" / "note.txt", "hello")
            write_json(root / "a" / "b" / "data.json", {"x": 1})
            write_csv(root / "a" / "b" / "rows.csv", ["id", "v"], [{"id": "r1", "v": 2}])
            self.assertEqual((root / "a/b/note.txt").read_text(encoding="utf-8"), "hello")
            self.assertIn('"x": 1', (root / "a/b/data.json").read_text(encoding="utf-8"))
            self.assertIn("r1,2", (root / "a/b/rows.csv").read_text(encoding="utf-8"))
            self.assertFalse([p for p in (root / "a/b").iterdir() if p.name.startswith(".")])

    def test_fsync_failure_leaves_no_temp_or_partial_destination(self) -> None:
        # _finish must clean the temp file even when flush/fsync/close fails,
        # not only when os.replace fails.
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            dest = root / "out.txt"
            with patch("pnn.reports.os.fsync", side_effect=OSError("simulated fsync failure")):
                with self.assertRaises(OSError):
                    write_text(dest, "hello")
            self.assertFalse(dest.exists())
            leftovers = [p.name for p in root.iterdir() if p.name.startswith(".out.txt.")]
            self.assertEqual(leftovers, [], f"temp file leaked: {leftovers}")


if __name__ == "__main__":
    unittest.main()
