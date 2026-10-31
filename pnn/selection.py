"""End-to-end Phase 2 lead selection for the PNN project."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple

from .compounds import load_compounds
from .curation import Candidate, curate, deduplicate, prepare_candidates
from .docking import load_docking_scores
from .rank import RANK_FIELDS, ShortlistRow, build_shortlist
from .reports import write_csv, write_json, write_text
from .screen import TIER_ORDER, ScreenResult, summarize
from .tables import CURATED_FIELDS, SCREEN_FIELDS, curated_row, screen_row

SCHEMA = "pnn.phase2.pipeline/1"
LIMITS = (
    "Tier labels and consensus scores are prioritisation signals, not measured "
    "permeability, binding or PNN softening.",
    "Consensus normalization depends on the docked set; compare runs only with a fixed set.",
    "Exact-SMILES curation is not chemical canonicalization.",
)

_CHUNK = 1 << 20


class PipelineError(ValueError):
    """Raised when a Phase 2 run cannot produce a defensible lead set."""


@dataclass(frozen=True)
class LeadSelection:
    """The selected compounds and the audit counts for one Phase 2 run."""

    curated: Tuple[Candidate, ...]
    screened: Tuple[ScreenResult, ...]
    shortlist: Tuple[ShortlistRow, ...]
    inputs: Tuple[Mapping[str, object], ...]
    docking: Mapping[str, object]
    curation: Mapping[str, object]
    screening: Mapping[str, object]
    ranking: Mapping[str, object]


@dataclass(frozen=True)
class LeadRun:
    """A completed run: selected compounds, written artifacts and the manifest."""

    selection: LeadSelection
    outputs: Mapping[str, Path]
    manifest: Mapping[str, object]


def _digest(path: Path) -> str:
    """Return the SHA-256 hex digest of a file without loading it whole."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_inputs(paths: Sequence[str]) -> Tuple[Tuple[str, list], ...]:
    if not paths:
        raise PipelineError("no input compound CSV was supplied")
    return tuple((str(Path(path)), load_compounds(Path(path), allow_duplicates=True)) for path in paths)


def _screen_table(results: Sequence[ScreenResult]) -> Dict[str, Dict[str, str]]:
    return {
        result.compound.compound_id: {
            "tier": result.tier,
            "score": f"{result.score:.4f}",
            "name": result.compound.name,
            "smiles": result.compound.smiles,
        }
        for result in results
    }


def select_leads(
    inputs: Sequence[str],
    docking: str,
    min_tier: str = "B",
    per_bin: int = 25,
    bin_width: float = 0.1,
    max_compounds: Optional[int] = None,
    max_rmsd: Optional[float] = None,
    min_vina: Optional[float] = None,
) -> LeadSelection:
    """Curate, screen and rank one lead set from compound and docking CSVs."""
    loaded = _load_inputs(inputs)
    candidates = prepare_candidates(loaded)
    unique, duplicate_counts = deduplicate(candidates)
    selected, curation_summary = curate(
        unique,
        min_tier=min_tier,
        per_bin=per_bin,
        bin_width=bin_width,
        max_compounds=max_compounds,
    )
    if not selected:
        raise PipelineError(
            "no compound survived curation; every input row was rejected or the input set was empty"
        )

    screened = tuple(
        sorted(
            (candidate.result for candidate in selected),
            key=lambda result: (TIER_ORDER[result.tier], -result.score, result.compound.compound_id),
        )
    )
    scores = load_docking_scores(Path(docking))
    shortlist, ranking_summary = build_shortlist(
        scores,
        screen=_screen_table(screened),
        min_tier=min_tier,
        max_rmsd=max_rmsd,
        min_vina=min_vina,
    )
    if not shortlist:
        raise PipelineError(
            "no compound survived docking and screening filters; "
            "supply Vina-scored docking rows for curated compounds or relax the filters"
        )

    return LeadSelection(
        curated=tuple(selected),
        screened=screened,
        shortlist=tuple(shortlist),
        inputs=tuple(
            {"path": source, "sha256": _digest(Path(source)), "rows": len(compounds)}
            for source, compounds in loaded
        ),
        docking={"path": str(Path(docking)), "sha256": _digest(Path(docking)), "rows": len(scores)},
        curation={**duplicate_counts, **curation_summary, "input_rows": len(candidates)},
        screening={**summarize(screened), "emitted": len(screened)},
        ranking=ranking_summary,
    )


def _options(
    min_tier: str,
    per_bin: int,
    bin_width: float,
    max_compounds: Optional[int],
    max_rmsd: Optional[float],
    min_vina: Optional[float],
) -> Dict[str, object]:
    return {
        "min_tier": min_tier,
        "per_bin": per_bin,
        "bin_width": bin_width,
        "max_compounds": max_compounds,
        "max_rmsd": max_rmsd,
        "min_vina": min_vina,
    }


def _manifest(
    selection: LeadSelection, outputs: Mapping[str, Path], options: Mapping[str, object]
) -> Dict[str, object]:
    return {
        "schema": SCHEMA,
        "inputs": [dict(item) for item in selection.inputs],
        "docking": dict(selection.docking),
        "options": dict(options),
        "stages": {
            "curation": dict(selection.curation),
            "screening": dict(selection.screening),
            "ranking": dict(selection.ranking),
        },
        "outputs": {
            name: {"path": str(path), "sha256": _digest(path)} for name, path in outputs.items()
        },
        "limits": list(LIMITS),
    }


def render_lead_report(
    selection: LeadSelection,
    options: Mapping[str, object],
    top_n: int = 25,
) -> str:
    """Render a comprehensive Markdown lead report honoring G-001, G-002, and G-003."""
    lines = [
        "# PNN Phase 2 Lead Selection Report",
        "",
        "## Scientific Context & Triage Limits",
        "",
        "> [!IMPORTANT]",
        "> **Project Guidelines & Evidence Limits (G-001, G-002, G-003)**:",
        "> 1. **G-001 (ADR-001)**: Screening tiers (A, B, C) are triage signals only; never report a "
        "compound as BBB-penetrant or efficacious without measured permeability and PNN-softening data.",
        "> 2. **G-002 (ADR-002)**: Docking consensus scores are relative to the docked candidate set and "
        "never substitute for measured binding, PNN-softening or permeability data.",
        "> 3. **G-003 (ADR-003)**: Curated libraries retain input source provenance; exact-SMILES deduplication "
        "is not chemical canonicalization.",
        "",
        "## Input Provenance",
        "",
        "### Compound Sources",
        "",
        "| Source Path | Rows | SHA-256 Digest |",
        "|---|---|---|",
    ]
    for inp in selection.inputs:
        lines.append(f"| `{inp['path']}` | {inp['rows']} | `{inp['sha256']}` |")

    lines.extend([
        "",
        "### Docking Source",
        "",
        f"- **Path**: `{selection.docking['path']}`",
        f"- **Docking Rows**: {selection.docking['rows']}",
        f"- **SHA-256**: `{selection.docking['sha256']}`",
        "",
        "## Stage Summaries",
        "",
        "### 1. Library Curation & Diversity",
        "",
        f"- **Input Rows**: {selection.curation.get('input_rows', 0)}",
        f"- **Duplicate IDs Resolved**: {selection.curation.get('duplicate_ids', 0)}",
        f"- **Duplicate SMILES Resolved**: {selection.curation.get('duplicate_smiles', 0)}",
        f"- **Retained by Tier Floor ({options.get('min_tier', 'B')} or better)**: "
        f"{selection.curation.get('eligible', 0)} (removed {selection.curation.get('skipped_by_tier', 0)})",
        f"- **Surviving Bin Caps ({options.get('per_bin', 25)}/bin)**: "
        f"{selection.curation.get('selected', 0)} (removed {selection.curation.get('eligible', 0) - selection.curation.get('selected', 0)})",
    ])
    if selection.curation.get("truncated_by_max"):
        lines.append(f"- **Truncated by Global Cap**: True (capped at {selection.curation.get('max_compounds')})")

    by_tier = selection.screening.get("by_tier", {})
    lines.extend([
        "",
        "### 2. CNS/BBB Property Screening",
        "",
        f"- **Eligible Candidates Screened**: {selection.screening.get('emitted', 0)}",
        f"- **Tier A**: {by_tier.get('A', 0)}",
        f"- **Tier B**: {by_tier.get('B', 0)}",
        f"- **Tier C**: {by_tier.get('C', 0)}",
        f"- **Reject**: {by_tier.get('Reject', 0)}",
        "",
        "### 3. Docking Consensus Ranking",
        "",
        f"- **Docked Rows Ingested**: {selection.docking.get('rows', 0)}",
        f"- **Shortlist Compounds Emitted**: {len(selection.shortlist)}",
        f"- **Full Engine Coverage (Vina + Glide + MM-GBSA)**: {selection.ranking.get('full_coverage', 0)}",
        f"- **Partial Engine Coverage**: {selection.ranking.get('partial_coverage', 0)}",
        f"- **Skipped by Tier Filter**: {selection.ranking.get('skipped_by_tier', 0)}",
        f"- **Filtered by RMSD Threshold**: {selection.ranking.get('filtered_by_rmsd', 0)}",
        f"- **Filtered by Minimum Vina Affinity**: {selection.ranking.get('filtered_by_vina', 0)}",
        "",
        f"## Top Shortlisted Leads (Top {min(len(selection.shortlist), top_n)})",
        "",
        "| Rank | Compound ID | Name | Tier | Property Score | Best Target | Consensus | Coverage | Vina (kcal/mol) | Glide | MM-GBSA (kcal/mol) | RMSD (A) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ])

    for row in selection.shortlist[:top_n]:
        glide_str = f"{row.glide_score:.2f}" if row.glide_score is not None else "-"
        mmgbsa_str = f"{row.mmgbsa_dg:.2f}" if row.mmgbsa_dg is not None else "-"
        rmsd_str = f"{row.rmsd:.2f}" if row.rmsd is not None else "-"
        lines.append(
            f"| {row.rank} | `{row.compound_id}` | {row.name or '-'} | {row.tier} | "
            f"{row.property_score} | {row.target} | {row.consensus:.4f} | {row.coverage:.2f} | "
            f"{row.vina_affinity:.2f} | {glide_str} | {mmgbsa_str} | {rmsd_str} |"
        )

    lines.extend([
        "",
        "## Evidence Gaps & Recommended Experimental Progression",
        "",
        "1. **Physicochemical vs. Measured Permeability**: All tier allocations stem from in silico "
        "physicochemical limits (MW 200-450, logP 1.5-3.5, TPSA < 70 A2). High-priority leads require "
        "in vitro PAMPA-BBB and MDCK-MDR1 validation.",
        "2. **Target Engagement**: Consensus docking prioritises poses against PTPsigma, LAR, or NgR1. "
        "Direct biophysical binding (SPR/ITC) must confirm dissociation constants (Kd).",
        "3. **Functional Efficacy**: Verified binders must proceed to primary neuronal/glial cultures "
        "for WFA-positive perineuronal net digestion quantification.",
        "",
    ])
    return "\n".join(lines)


def run_pipeline(
    inputs: Sequence[str],
    docking: str,
    out_dir: str,
    min_tier: str = "B",
    per_bin: int = 25,
    bin_width: float = 0.1,
    max_compounds: Optional[int] = None,
    max_rmsd: Optional[float] = None,
    min_vina: Optional[float] = None,
    top_n_report: int = 25,
) -> LeadRun:
    """Run the three Phase 2 stages and write their artifacts, report and manifest."""
    selection = select_leads(
        inputs,
        docking,
        min_tier=min_tier,
        per_bin=per_bin,
        bin_width=bin_width,
        max_compounds=max_compounds,
        max_rmsd=max_rmsd,
        min_vina=min_vina,
    )
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=True)
    outputs = {
        "curated": root / "curated_library.csv",
        "ranked": root / "ranked_compounds.csv",
        "shortlist": root / "shortlist.csv",
        "report": root / "lead_report.md",
    }
    write_csv(
        outputs["curated"],
        CURATED_FIELDS,
        [curated_row(candidate, bin_width) for candidate in selection.curated],
    )
    write_csv(outputs["ranked"], SCREEN_FIELDS, [screen_row(result) for result in selection.screened])
    write_csv(outputs["shortlist"], RANK_FIELDS, [row.as_row() for row in selection.shortlist])

    options = _options(min_tier, per_bin, bin_width, max_compounds, max_rmsd, min_vina)
    report_text = render_lead_report(selection, options, top_n=top_n_report)
    write_text(outputs["report"], report_text)

    manifest = _manifest(selection, outputs, options)
    manifest_path = root / "manifest.json"
    write_json(manifest_path, manifest)

    all_outputs = {**outputs, "manifest": manifest_path}
    return LeadRun(selection=selection, outputs=all_outputs, manifest=manifest)
