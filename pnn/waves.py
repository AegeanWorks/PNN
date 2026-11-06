"""Wave-ledger aggregation: search-wave CSVs plus obtainability proofs into report artifacts.

Every search wave drops a ``wave<NNN>_*.csv`` screening set into
``data/compounds/waves/`` and appends proof lines to
``obtainability_proofs.csv``. Nothing in the turnkey pipeline read those
files, so wave findings (shortlist entries, parks, kills) never reached
``output/`` or ``studio.html``. This module closes that gap: it counts
per-wave compounds, parses shortlist/park/kill tags out of the proofs,
and renders the ledger as CSV, Markdown, and a studio HTML panel.
Standard library only.
"""

from __future__ import annotations

import csv
import html
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

WAVE_FILE_RE = re.compile(r"wave(\d+)_.*\.csv$", re.IGNORECASE)
WAVE_TAG_RE = re.compile(r"^W(\d+)-")
KILLS_MAPPED_RE = re.compile(r"(\d+)\s*-\s*kills-mapped", re.IGNORECASE)
ROUTE_KILLS_RE = re.compile(r"(\d+)\s*-\s*route-kills?", re.IGNORECASE)

WAVE_LEDGER_FIELDS = (
    "wave",
    "source_file",
    "compounds",
    "shortlist",
    "parks",
    "kills",
)


class WaveLedgerError(ValueError):
    """Raised when wave-ledger inputs cannot be read or parsed."""

SHORTLIST_KINDS = ("SHORTLIST",)
PARK_KINDS = ("PARKED",)


@dataclass(frozen=True)
class WaveRow:
    wave: str
    source_file: str
    compounds: int
    shortlist: int
    parks: int
    kills: int

    def as_row(self) -> Dict[str, object]:
        return {
            "wave": self.wave,
            "source_file": self.source_file,
            "compounds": self.compounds,
            "shortlist": self.shortlist,
            "parks": self.parks,
            "kills": self.kills,
        }


@dataclass(frozen=True)
class WaveHighlight:
    tag: str
    wave: str
    kind: str
    note: str

    def as_row(self) -> Dict[str, object]:
        return {"tag": self.tag, "wave": self.wave, "kind": self.kind, "note": self.note}


def _wave_number_from_filename(filename: str) -> Optional[str]:
    match = WAVE_FILE_RE.match(filename)
    return match.group(1) if match else None


def _count_wave_compounds(path: Path) -> int:
    """Count data rows in a wave screening CSV (header excluded)."""
    try:
        # errors="replace": one stray byte must never break turnkey reporting.
        # Raw line counting (not DictReader) keeps the full-ledger scan cheap.
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            lines = 0
            non_empty = 0
            for line in handle:
                lines += 1
                if line.strip():
                    non_empty += 1
            return max(0, non_empty - (1 if lines else 0))
    except OSError as exc:
        raise WaveLedgerError(
            f"What: cannot read wave file {path.name}. "
            f"Why: {exc}. Fix: restore the file or remove it before re-running."
        ) from exc


def _classify_proof_row(entry_id: str, note: str) -> Tuple[Optional[str], str, str, int]:
    """Return (wave, kind, serial, kill_count) for one obtainability_proofs row id.

    kind is one of "shortlist", "park", "kill", or "" when unrecognized.
    serial identifies the distinct candidate (e.g. "S01", "SL-01", "P01") so
    repeated vendor-proof lines for one candidate count once. Kill lines may
    map several kills at once ("13-kills-mapped").
    """
    tag_match = WAVE_TAG_RE.match(entry_id.strip())
    if not tag_match:
        return None, "", "", 0
    wave = tag_match.group(1)
    serial = entry_id.strip().split("-", 1)[1].upper() if "-" in entry_id else ""
    upper_note = note.upper()
    if re.match(r"^S\d+$", serial) or re.match(r"^SL-\d+$", serial):
        return wave, "shortlist", serial, 0
    if re.match(r"^P\d+$", serial):
        return wave, "park", serial, 0
    if "HOLD" in upper_note:
        return wave, "park", serial, 0
    if "SHORTLIST" in upper_note:
        return wave, "shortlist", serial, 0
    if "PARKED" in upper_note:
        return wave, "park", serial, 0
    if "KILL" in serial or "DEAD" in serial:
        mapped = KILLS_MAPPED_RE.search(note)
        if mapped:
            return wave, "kill", serial, int(mapped.group(1))
        routed = ROUTE_KILLS_RE.search(note)
        if routed:
            return wave, "kill", serial, int(routed.group(1))
        if "ROUTE-KILL" in upper_note or "ROUTE KILL" in upper_note:
            return wave, "kill", serial, 1
        return wave, "kill", serial, 1
    return wave, "", "", 0


def collect_wave_ledger(
    waves_dir: Path, proofs_path: Path
) -> Tuple[List[WaveRow], List[WaveHighlight]]:
    """Aggregate per-wave compound counts plus shortlist/park/kill tags."""
    compound_counts: Dict[str, int] = {}
    source_files: Dict[str, str] = {}
    if waves_dir.is_dir():
        for path in sorted(waves_dir.glob("wave*.csv")):
            if path.name.lower() == "obtainability_proofs.csv":
                continue
            wave = _wave_number_from_filename(path.name)
            if wave is None:
                continue
            compound_counts[wave] = compound_counts.get(wave, 0) + _count_wave_compounds(path)
            source_files.setdefault(wave, path.name)

    shortlist_serials: Dict[str, set] = {}
    park_serials: Dict[str, set] = {}
    kill_counts: Dict[str, int] = {}
    highlight_seen: set = set()
    highlights: List[WaveHighlight] = []
    if proofs_path.is_file():
        try:
            # errors="replace": one stray byte must never break turnkey reporting.
            with proofs_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
                proof_rows = list(csv.reader(handle))
        except (OSError, csv.Error) as exc:
            raise WaveLedgerError(
                f"What: cannot read obtainability proofs {proofs_path.name}. "
                f"Why: {exc}. Fix: restore the file before re-running."
            ) from exc
        for row in proof_rows:
            if not row:
                continue
            entry_id = (row[0] or "").strip()
            note = " ".join(cell.strip() for cell in row[3:] if cell and cell.strip())
            wave, kind, serial, kill_count = _classify_proof_row(entry_id, note)
            if wave is None or not kind:
                continue
            if kind == "shortlist":
                key = (wave, serial or entry_id)
                shortlist_serials.setdefault(wave, set()).add(key)
                if key not in highlight_seen:
                    highlight_seen.add(key)
                    highlights.append(
                        WaveHighlight(
                            tag=entry_id, wave=wave, kind="shortlist", note=note[:400]
                        )
                    )
            elif kind == "park":
                key = (wave, serial or entry_id)
                park_serials.setdefault(wave, set()).add(key)
                if key not in highlight_seen:
                    highlight_seen.add(key)
                    highlights.append(
                        WaveHighlight(
                            tag=entry_id, wave=wave, kind="park", note=note[:400]
                        )
                    )
            elif kind == "kill":
                kill_counts[wave] = kill_counts.get(wave, 0) + max(kill_count, 1)

    shortlist_counts = {wave: len(keys) for wave, keys in shortlist_serials.items()}
    park_counts = {wave: len(keys) for wave, keys in park_serials.items()}

    waves = sorted(
        set(compound_counts) | set(shortlist_counts) | set(park_counts) | set(kill_counts),
        key=lambda w: (len(w), w),
    )
    rows = [
        WaveRow(
            wave=wave,
            source_file=source_files.get(wave, "-"),
            compounds=compound_counts.get(wave, 0),
            shortlist=shortlist_counts.get(wave, 0),
            parks=park_counts.get(wave, 0),
            kills=kill_counts.get(wave, 0),
        )
        for wave in waves
    ]
    return rows, highlights


def render_wave_ledger_markdown(
    rows: Sequence[WaveRow], highlights: Sequence[WaveHighlight]
) -> str:
    """Render the wave ledger as a Markdown report section."""
    total_compounds = sum(r.compounds for r in rows)
    total_shortlist = sum(r.shortlist for r in rows)
    total_parks = sum(r.parks for r in rows)
    total_kills = sum(r.kills for r in rows)
    lines = [
        "# Wave Ledger — Search-Wave Findings",
        "",
        f"Waves recorded: {len(rows)} | Compounds screened: {total_compounds} | "
        f"Shortlist: {total_shortlist} | Parked: {total_parks} | Kills: {total_kills}",
        "",
        "| Wave | Source File | Compounds | Shortlist | Parks | Kills |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row.wave} | {row.source_file} | {row.compounds} | "
            f"{row.shortlist} | {row.parks} | {row.kills} |"
        )
    if highlights:
        lines += ["", "## Highlights (shortlist + parks)", ""]
        for highlight in highlights:
            lines.append(f"- `{highlight.tag}` (wave {highlight.wave}, {highlight.kind}): {highlight.note}")
    lines.append("")
    return "\n".join(lines)


def render_wave_section(
    rows: Sequence[Mapping[str, object]], highlights: Sequence[Mapping[str, object]]
) -> str:
    """Render the wave ledger as a standalone studio HTML panel section."""
    total_compounds = sum(int(r.get("compounds", 0) or 0) for r in rows)
    total_shortlist = sum(int(r.get("shortlist", 0) or 0) for r in rows)
    total_parks = sum(int(r.get("parks", 0) or 0) for r in rows)
    total_kills = sum(int(r.get("kills", 0) or 0) for r in rows)
    body = [
        '<div class="panel" style="margin-bottom: 24px;">',
        f"<h2>Wave Ledger — Search-Wave Findings "
        f'<span style="font-size:0.8rem; font-weight:400; color:var(--text-muted);">'
        f"{len(rows)} waves | {total_compounds} compounds | "
        f"{total_shortlist} shortlist | {total_parks} parked | {total_kills} kills</span></h2>",
        '<div class="table-wrapper"><table class="data-table">',
        "<thead><tr><th>Wave</th><th>Source File</th><th>Compounds</th>"
        "<th>Shortlist</th><th>Parks</th><th>Kills</th></tr></thead><tbody>",
    ]
    for row in rows:
        body.append(
            "<tr><td>{wave}</td><td>{source}</td><td>{compounds}</td>"
            "<td>{shortlist}</td><td>{parks}</td><td>{kills}</td></tr>".format(
                wave=html.escape(str(row.get("wave", "-"))),
                source=html.escape(str(row.get("source_file", "-"))),
                compounds=html.escape(str(row.get("compounds", 0))),
                shortlist=html.escape(str(row.get("shortlist", 0))),
                parks=html.escape(str(row.get("parks", 0))),
                kills=html.escape(str(row.get("kills", 0))),
            )
        )
    body.append("</tbody></table></div>")
    if highlights:
        body.append("<h3>Highlights (shortlist + parks)</h3>")
        body.append('<div class="table-wrapper"><table class="data-table">')
        body.append(
            "<thead><tr><th>Tag</th><th>Wave</th><th>Kind</th><th>Evidence</th></tr></thead><tbody>"
        )
        for highlight in highlights:
            body.append(
                "<tr><td>{tag}</td><td>{wave}</td><td>{kind}</td><td>{note}</td></tr>".format(
                    tag=html.escape(str(highlight.get("tag", "-"))),
                    wave=html.escape(str(highlight.get("wave", "-"))),
                    kind=html.escape(str(highlight.get("kind", "-"))),
                    note=html.escape(str(highlight.get("note", "-"))),
                )
            )
        body.append("</tbody></table></div>")
    body.append("</div>")
    return "\n".join(body)
