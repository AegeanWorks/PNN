import tempfile
import unittest
from pathlib import Path

from pnn.waves import (
    WaveLedgerError,
    _classify_proof_row,
    _count_wave_compounds,
    collect_wave_ledger,
    render_wave_ledger_markdown,
    render_wave_section,
)


WAVE_CSV = """compound_id,name,smiles,mw,logp,tpsa,hbd,hba,rotatable_bonds,pka,charge,logbb,brain_plasma_ratio,oral_bioavailability,cyp3a4_inhibition,herg_liability,pnn_ec50_nm
W999-K01,Kill-a,,300.0,5.0,30.0,0,2,2,,,,,,,
W999-P01,Park-a,,300.0,2.0,30.0,1,3,2,,,,,,,
W999-S01,Short-a,,300.0,2.0,30.0,1,3,2,,,,,,,
"""

PROOFS_CSV = """compound_id,name,smiles,mw,logp,tpsa
W999-S01,,,FDA-approved SHORTLIST test entry with direct datum
W999-P01,,,FDA-approved PARKED hook gap test entry
W999-P01,MCE,HY-00000,,second vendor proof line for same park
W999-KILL-A,,,2-kills-mapped Kill-a logP plus one more kill
"""


class TestWaveLedger(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.waves_dir = Path(self.tmp.name)
        (self.waves_dir / "wave999_test.csv").write_text(WAVE_CSV, encoding="utf-8")
        self.proofs = self.waves_dir / "obtainability_proofs.csv"
        self.proofs.write_text(PROOFS_CSV, encoding="utf-8")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_collect_counts(self) -> None:
        rows, highlights = collect_wave_ledger(self.waves_dir, self.proofs)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.wave, "999")
        self.assertEqual(row.compounds, 3)
        self.assertEqual(row.shortlist, 1)
        self.assertEqual(row.parks, 1)
        self.assertEqual(row.kills, 2)
        kinds = sorted(h.kind for h in highlights)
        self.assertEqual(kinds, ["park", "shortlist"])

    def test_markdown_report(self) -> None:
        rows, highlights = collect_wave_ledger(self.waves_dir, self.proofs)
        report = render_wave_ledger_markdown(rows, highlights)
        self.assertIn("Wave Ledger", report)
        self.assertIn("999", report)
        self.assertIn("W999-S01", report)

    def test_html_section(self) -> None:
        rows, highlights = collect_wave_ledger(self.waves_dir, self.proofs)
        section = render_wave_section(
            [row.as_row() for row in rows], [h.as_row() for h in highlights]
        )
        self.assertIn("Wave Ledger", section)
        self.assertIn("wave999_test.csv", section)
        self.assertIn("W999-S01", section)

    def test_missing_inputs_yield_empty_ledger(self) -> None:
        rows, highlights = collect_wave_ledger(
            Path(self.tmp.name) / "nope", Path(self.tmp.name) / "missing.csv"
        )
        self.assertEqual(rows, [])
        self.assertEqual(highlights, [])

    def test_unreadable_wave_file_raises_actionable_error(self) -> None:
        with self.assertRaises(WaveLedgerError):
            _count_wave_compounds(Path(self.tmp.name))

    def test_stray_byte_does_not_break_ledger(self) -> None:
        proofs = self.waves_dir / "obtainability_proofs.csv"
        raw = proofs.read_bytes().replace(b"FDA-approved", b"FDA\x97approved")
        proofs.write_bytes(raw)
        rows, highlights = collect_wave_ledger(self.waves_dir, proofs)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].kills, 2)

    def test_sched_and_box_pass_kill_counts_expand(self) -> None:
        proofs = self.waves_dir / "obtainability_proofs.csv"
        proofs.write_text(
            "compound_id,supplier,catalog_no,approval_ref,evidence,verified_date\n"
            "W999-KILL-SCHED,,,7-sched-kills,Seven scheduled kills,2026-11-06\n"
            "W999-KILL-BOX,,,3-box-pass-kills,Three box-pass kills,2026-11-06\n",
            encoding="utf-8",
        )
        rows, _ = collect_wave_ledger(self.waves_dir, proofs)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].kills, 10)

    def test_route_kills_expand(self) -> None:
        proofs = self.waves_dir / "obtainability_proofs.csv"
        proofs.write_text(
            "compound_id,supplier,catalog_no,approval_ref,evidence,verified_date\n"
            "W999-KILL-ROUTE,,,8-route-kills,Eight injectable-only kills,2026-11-06\n",
            encoding="utf-8",
        )
        rows, _ = collect_wave_ledger(self.waves_dir, proofs)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].kills, 8)

    def test_ledger_explains_administrative_only_waves(self) -> None:
        # Proof rows such as buy-lists, MPO triages and reagent/synonym
        # confirmations carry no shortlist/park/kill outcome, so their wave
        # numbers are intentionally absent. Both reports must say so, or a
        # reader comparing against W347 proof ids assumes data was lost.
        rows, _ = collect_wave_ledger(self.waves_dir, self.proofs)
        markdown = render_wave_ledger_markdown(rows, [])
        html = render_wave_section([r.as_row() for r in rows], [])
        self.assertIn("administrative", markdown.lower())
        self.assertIn("administrative", html.lower())

    def test_kill_id_wins_over_status_words_in_note(self) -> None:
        # "parked" and "holds" appear in real kill-row evidence as prose / a
        # pipeline tally ("B=1(parked)", "wall holds"), not as classifications.
        # The authoritative KILL id must win over those note words.
        wave, kind, serial, kills = _classify_proof_row(
            "W322-KILL-COU",
            "3-kills-mapped,Acenocoumarol ... pipeline B=1(parked) C=2 Reject=1",
        )
        self.assertEqual((wave, kind, kills), ("322", "kill", 3))

        wave, kind, serial, kills = _classify_proof_row(
            "W216-KILL-GLYCO",
            "7-kills-mapped,... glyco/lipopeptide wall holds; pipeline Reject 7/7",
        )
        self.assertEqual((wave, kind, kills), ("216", "kill", 7))

        proofs = self.waves_dir / "obtainability_proofs.csv"
        proofs.write_text(
            "compound_id,supplier,catalog_no,approval_ref,evidence,verified_date\n"
            "W322-KILL-COU,,,3-kills-mapped,... pipeline B=1(parked) C=2 Reject=1,2026-11-06\n"
            "W216-KILL-GLYCO,,,7-kills-mapped,... wall holds; pipeline Reject 7/7,2026-11-06\n",
            encoding="utf-8",
        )
        rows, _ = collect_wave_ledger(self.waves_dir, proofs)
        self.assertEqual(sum(r.kills for r in rows), 10)
        self.assertEqual(sum(r.parks for r in rows), 0)


if __name__ == "__main__":
    unittest.main()
