"""Unit and CLI tests for target receptor profiles and Vina grid generation."""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from pnn.target import (
    TARGET_PROFILES,
    generate_vina_config,
    get_target,
    list_targets,
    main,
)


class TargetProfileTest(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def run_cli(self, *argv: str):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            code = main(list(argv))
        return code, stdout.getvalue(), stderr.getvalue()

    def test_list_targets_contains_primary_receptors(self) -> None:
        targets = list_targets()
        self.assertGreaterEqual(len(targets), 5)
        target_ids = {t.target_id for t in targets}
        self.assertIn("PTPsigma", target_ids)
        self.assertIn("LAR", target_ids)
        self.assertIn("NgR1", target_ids)
        self.assertIn("MMP2", target_ids)
        self.assertIn("ChABC", target_ids)

    def test_get_target_case_insensitive(self) -> None:
        ptp = get_target("ptpsigma")
        self.assertEqual(ptp.target_id, "PTPsigma")
        self.assertEqual(ptp.primary_pdb, "5K4P")
        self.assertEqual(ptp.gene_symbol, "PTPRS")
        self.assertEqual(ptp.grid_center, (12.5, 18.3, 22.1))
        self.assertEqual(ptp.grid_size, (20.0, 20.0, 20.0))

    def test_get_target_unknown_raises_key_error(self) -> None:
        with self.assertRaises(KeyError) as ctx:
            get_target("non_existent_target")
        self.assertIn("Unknown target ID", str(ctx.exception))
        self.assertIn("PTPsigma", str(ctx.exception))

    def test_generate_vina_config_format(self) -> None:
        target = get_target("LAR")
        config = generate_vina_config(
            profile=target,
            receptor_pdbqt="lar.pdbqt",
            out_pdbqt="lar_docked.pdbqt",
            exhaustiveness=16,
        )
        self.assertIn("receptor = lar.pdbqt", config)
        self.assertIn("out = lar_docked.pdbqt", config)
        self.assertIn("center_x = 18.100", config)
        self.assertIn("center_y = 10.400", config)
        self.assertIn("center_z = 15.800", config)
        self.assertIn("size_x = 20.0", config)
        self.assertIn("exhaustiveness = 16", config)

    def test_cli_list(self) -> None:
        code, stdout, _ = self.run_cli("list")
        self.assertEqual(code, 0)
        self.assertIn("PTPsigma", stdout)
        self.assertIn("5K4P", stdout)
        self.assertIn("LAR", stdout)

    def test_cli_show_success_and_failure(self) -> None:
        code, stdout, _ = self.run_cli("show", "NgR1")
        self.assertEqual(code, 0)
        self.assertIn("Nogo-66 Receptor 1", stdout)
        self.assertIn("2OZ7", stdout)

        code_fail, _, stderr = self.run_cli("show", "invalid_receptor")
        self.assertEqual(code_fail, 1)
        self.assertIn("What:", stderr)
        self.assertIn("Why:", stderr)
        self.assertIn("Fix:", stderr)

    def test_cli_config_export_to_file(self) -> None:
        dest = self.root / "vina_grid.conf"
        code, stdout, _ = self.run_cli("config", "PTPsigma", "--save", str(dest))
        self.assertEqual(code, 0)
        self.assertTrue(dest.is_file())
        text = dest.read_text(encoding="utf-8")
        self.assertIn("center_x = 12.500", text)
        self.assertIn("Saved Vina grid config", stdout)

    def test_aggrecanase_targets_present(self) -> None:
        targets = list_targets()
        target_ids = {t.target_id for t in targets}
        self.assertIn("ADAMTS4", target_ids)
        self.assertIn("MMP9", target_ids)

    def test_adamts4_target_profile(self) -> None:
        adamts = get_target("ADAMTS4")
        self.assertEqual(adamts.target_id, "ADAMTS4")
        self.assertEqual(adamts.gene_symbol, "ADAMTS4")
        self.assertEqual(adamts.primary_pdb, "2RJP")
        self.assertEqual(adamts.uniprot_id, "O75173")
        self.assertEqual(adamts.pocket_type, "catalytic_zinc")
        self.assertIn("His202", adamts.key_residues)
        self.assertIn("Glu203", adamts.key_residues)

    def test_mmp9_target_profile(self) -> None:
        mmp9 = get_target("MMP9")
        self.assertEqual(mmp9.target_id, "MMP9")
        self.assertEqual(mmp9.gene_symbol, "MMP9")
        self.assertEqual(mmp9.primary_pdb, "1GKC")
        self.assertEqual(mmp9.uniprot_id, "P14780")
        self.assertEqual(mmp9.pocket_type, "catalytic_zinc")
        self.assertIn("His401", mmp9.key_residues)
        self.assertIn("Glu402", mmp9.key_residues)


if __name__ == "__main__":
    unittest.main()
