"""Molecular Dynamics (MD) trajectory and contact stability analyzer for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .reports import write_csv, write_json, write_text


class DynamicsError(ValueError):
    """Raised when MD trajectory data or parameters are malformed."""


def _fmt(value: Optional[float], spec: str) -> str:
    return format(value, spec) if value is not None else ""


def _display(value: Optional[float], spec: str) -> str:
    return format(value, spec) if value is not None else "UNMEASURED"


@dataclass(frozen=True)
class TrajectoryPoint:
    """A single snapshot / frame in an MD simulation trajectory."""

    time_ns: float
    ligand_rmsd: float
    protein_rmsd: Optional[float] = None
    contacts_count: Optional[int] = None
    mmgbsa_dg: Optional[float] = None
    residue_contacts: Dict[str, bool] = field(default_factory=dict)


@dataclass(frozen=True)
class TrajectoryAnalysis:
    """Comprehensive statistical and structural analysis of an MD trajectory."""

    compound_id: str
    target_id: str
    total_time_ns: float
    frames_count: int
    mean_ligand_rmsd: float
    std_ligand_rmsd: float
    max_ligand_rmsd: float
    drift_ligand_rmsd: float
    mean_protein_rmsd: Optional[float]
    std_protein_rmsd: Optional[float]
    mean_contacts: Optional[float]
    min_contacts: Optional[int]
    mean_mmgbsa_dg: Optional[float]
    std_mmgbsa_dg: Optional[float]
    residue_persistence: Dict[str, float]
    stability_class: str
    stability_score: Optional[float]
    evidence_status: str
    verdict_rationale: str

    def as_dict(self) -> Dict[str, object]:
        def rounded(value: Optional[float], digits: int) -> Optional[float]:
            return round(value, digits) if value is not None else None

        return {
            "compound_id": self.compound_id,
            "target_id": self.target_id,
            "total_time_ns": round(self.total_time_ns, 2),
            "frames_count": self.frames_count,
            "mean_ligand_rmsd": round(self.mean_ligand_rmsd, 3),
            "std_ligand_rmsd": round(self.std_ligand_rmsd, 3),
            "max_ligand_rmsd": round(self.max_ligand_rmsd, 3),
            "drift_ligand_rmsd": round(self.drift_ligand_rmsd, 3),
            "mean_protein_rmsd": rounded(self.mean_protein_rmsd, 3),
            "std_protein_rmsd": rounded(self.std_protein_rmsd, 3),
            "mean_contacts": rounded(self.mean_contacts, 2),
            "min_contacts": self.min_contacts,
            "mean_mmgbsa_dg": rounded(self.mean_mmgbsa_dg, 2),
            "std_mmgbsa_dg": rounded(self.std_mmgbsa_dg, 2),
            "residue_persistence": {k: round(v, 1) for k, v in self.residue_persistence.items()},
            "stability_class": self.stability_class,
            "stability_score": rounded(self.stability_score, 1),
            "evidence_status": self.evidence_status,
            "verdict_rationale": self.verdict_rationale,
        }

    def as_row(self) -> Dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "target_id": self.target_id,
            "total_time_ns": f"{self.total_time_ns:.1f}",
            "frames_count": self.frames_count,
            "mean_ligand_rmsd": f"{self.mean_ligand_rmsd:.2f}",
            "std_ligand_rmsd": f"{self.std_ligand_rmsd:.2f}",
            "max_ligand_rmsd": f"{self.max_ligand_rmsd:.2f}",
            "drift_ligand_rmsd": f"{self.drift_ligand_rmsd:.2f}",
            "mean_protein_rmsd": _fmt(self.mean_protein_rmsd, ".2f"),
            "mean_contacts": _fmt(self.mean_contacts, ".1f"),
            "min_contacts": self.min_contacts if self.min_contacts is not None else "",
            "mean_mmgbsa_dg": _fmt(self.mean_mmgbsa_dg, ".1f"),
            "stability_class": self.stability_class,
            "stability_score": _fmt(self.stability_score, ".1f"),
            "evidence_status": self.evidence_status,
        }


DYNAMICS_SUMMARY_FIELDS: Sequence[str] = (
    "compound_id",
    "target_id",
    "total_time_ns",
    "frames_count",
    "mean_ligand_rmsd",
    "std_ligand_rmsd",
    "max_ligand_rmsd",
    "drift_ligand_rmsd",
    "mean_protein_rmsd",
    "mean_contacts",
    "min_contacts",
    "mean_mmgbsa_dg",
    "stability_class",
    "stability_score",
    "evidence_status",
)


def _calc_stats(values: Sequence[float]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Calculate mean, standard deviation, and maximum of the supplied (non-None) values."""
    if not values:
        return None, None, None
    n = len(values)
    mean_val = sum(values) / n
    variance = sum((x - mean_val) ** 2 for x in values) / n
    return mean_val, math.sqrt(variance), max(values)


def _parse_optional_float(raw: Optional[str], label: str, row_idx: int) -> Optional[float]:
    if raw is None:
        return None
    clean = raw.strip()
    if not clean or clean.upper() in ("NA", "NAN", "NULL", "NONE", "-"):
        return None
    try:
        value = float(clean)
    except ValueError:
        raise DynamicsError(f"Error parsing trajectory row {row_idx}: invalid value '{raw}' for {label}")
    if not math.isfinite(value):
        raise DynamicsError(f"Error parsing trajectory row {row_idx}: non-finite value for {label}")
    return value


def load_trajectory_csv(path: Path) -> List[TrajectoryPoint]:
    """Load MD trajectory snapshots from a CSV file; absent optional columns stay unmeasured."""
    path = Path(path)
    if not path.is_file():
        raise DynamicsError(f"Trajectory file does not exist: {path}")

    points: List[TrajectoryPoint] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise DynamicsError(f"Trajectory CSV file is empty: {path}")

        headers = [h.strip() for h in reader.fieldnames]

        time_col = next((c for c in headers if c in ("time_ns", "time", "frame_ns")), None)
        lig_col = next(
            (c for c in headers if c in ("ligand_rmsd_angstrom", "ligand_rmsd", "rmsd_ligand", "rmsd_lig")), None
        )
        prot_col = next(
            (c for c in headers if c in ("protein_backbone_rmsd", "protein_rmsd", "rmsd_protein", "rmsd_rec")), None
        )
        contact_col = next(
            (c for c in headers if c in ("key_contacts_count", "contacts_count", "contacts", "num_contacts")), None
        )
        dg_col = next(
            (c for c in headers if c in ("mmgbsa_binding_dg", "mmgbsa_dg", "binding_dg", "dg_bind")), None
        )

        if not time_col or not lig_col:
            raise DynamicsError(
                f"Trajectory CSV must contain at least time_ns and ligand_rmsd columns. Found: {headers}"
            )

        contact_cols = [c for c in headers if c.startswith("contact_") or c.startswith("res_")]

        for row_idx, row in enumerate(reader, start=2):
            t = _parse_optional_float(row[time_col], time_col, row_idx)
            if t is None:
                raise DynamicsError(f"Error parsing trajectory row {row_idx}: {time_col} is required")
            lig_rmsd = _parse_optional_float(row[lig_col], lig_col, row_idx)
            if lig_rmsd is None:
                raise DynamicsError(f"Error parsing trajectory row {row_idx}: {lig_col} is required")

            prot_rmsd = _parse_optional_float(row.get(prot_col), prot_col or "", row_idx) if prot_col else None
            contacts: Optional[int] = None
            if contact_col:
                contact_value = _parse_optional_float(row.get(contact_col), contact_col, row_idx)
                if contact_value is not None:
                    if contact_value < 0:
                        raise DynamicsError(f"Error parsing trajectory row {row_idx}: {contact_col} must be >= 0")
                    contacts = int(contact_value)
            dg = _parse_optional_float(row.get(dg_col), dg_col or "", row_idx) if dg_col else None

            res_contacts: Dict[str, bool] = {}
            for c_col in contact_cols:
                val_str = (row.get(c_col) or "").strip().lower()
                res_contacts[c_col.replace("contact_", "").replace("res_", "").upper()] = val_str in (
                    "1",
                    "true",
                    "yes",
                    "t",
                    "y",
                )

            points.append(
                TrajectoryPoint(
                    time_ns=t,
                    ligand_rmsd=lig_rmsd,
                    protein_rmsd=prot_rmsd,
                    contacts_count=contacts,
                    mmgbsa_dg=dg,
                    residue_contacts=res_contacts,
                )
            )

    if len(points) < 2:
        raise DynamicsError(f"Trajectory must contain at least 2 frames for analysis. Found: {len(points)}")

    points.sort(key=lambda p: p.time_ns)
    return points


def analyze_trajectory(
    points: Sequence[TrajectoryPoint],
    compound_id: str = "PNN-Compound",
    target_id: str = "PTPsigma",
) -> TrajectoryAnalysis:
    """Analyze MD frames; metrics absent from the input stay unmeasured and block a stability verdict."""
    if len(points) < 2:
        raise DynamicsError("Need at least 2 trajectory points to analyze dynamics.")

    n_frames = len(points)
    total_time = points[-1].time_ns - points[0].time_ns
    if total_time <= 0:
        total_time = float(n_frames)

    lig_rmsds = [p.ligand_rmsd for p in points]
    prot_rmsds = [p.protein_rmsd for p in points if p.protein_rmsd is not None]
    contact_values = [float(p.contacts_count) for p in points if p.contacts_count is not None]
    dg_values = [p.mmgbsa_dg for p in points if p.mmgbsa_dg is not None]

    mean_lig, std_lig, max_lig = _calc_stats(lig_rmsds)
    assert mean_lig is not None and std_lig is not None and max_lig is not None
    mean_prot, std_prot, _ = _calc_stats(prot_rmsds)
    mean_cnt, _, _ = _calc_stats(contact_values)
    min_cnt = min(contact_values) if contact_values else None
    mean_dg, std_dg, _ = _calc_stats(dg_values)

    window_size = max(1, n_frames // 5)
    first_window = lig_rmsds[:window_size]
    last_window = lig_rmsds[-window_size:]
    drift = (sum(last_window) / len(last_window)) - (sum(first_window) / len(first_window))

    residue_names = sorted({name for p in points for name in p.residue_contacts})
    res_persistence: Dict[str, float] = {}
    for res_name in residue_names:
        active_count = sum(1 for p in points if p.residue_contacts.get(res_name, False))
        res_persistence[res_name] = (active_count / n_frames) * 100.0

    has_contacts = mean_cnt is not None
    has_protein = mean_prot is not None
    has_dg = mean_dg is not None

    if not has_contacts:
        stability_class = "INSUFFICIENT_EVIDENCE"
        rationale = (
            "Contact counts were not supplied, so ligand residence in the pocket cannot be assessed. "
            "This trajectory does not support a stability verdict."
        )
        score: Optional[float] = None
    else:
        if mean_lig < 2.0 and drift < 1.0 and mean_cnt >= 2.0 and max_lig < 3.2:
            stability_class = "STABLE_POSE"
            rationale = "Ligand maintains stable binding geometry in the pocket (RMSD < 2.0 Å, minimal drift, high contact persistence)."
        elif mean_lig <= 3.5 and drift < 2.0 and mean_cnt >= 1.0:
            stability_class = "MODERATE_MOBILITY"
            rationale = "Ligand exhibits moderate conformational mobility or pocket breathing (RMSD 2.0-3.5 Å) while preserving core contacts."
        else:
            stability_class = "UNSTABLE_DRIFT"
            rationale = "Ligand exhibits excessive conformational drift or dissociation from the binding pocket (RMSD > 3.5 Å or lost contacts)."

        score = 100.0
        if mean_lig > 1.5:
            score -= (mean_lig - 1.5) * 20.0
        if std_lig > 0.4:
            score -= (std_lig - 0.4) * 15.0
        if drift > 0.5:
            score -= (drift - 0.5) * 20.0
        if mean_cnt < 3.0:
            score -= (3.0 - mean_cnt) * 10.0
        if max_lig > 3.0:
            score -= (max_lig - 3.0) * 10.0
        score = max(0.0, min(100.0, score))

    if stability_class == "INSUFFICIENT_EVIDENCE":
        evidence_status = "INSUFFICIENT_EVIDENCE"
    elif has_protein and has_dg:
        evidence_status = "COMPLETE"
    else:
        evidence_status = "PARTIAL"

    return TrajectoryAnalysis(
        compound_id=compound_id,
        target_id=target_id,
        total_time_ns=total_time,
        frames_count=n_frames,
        mean_ligand_rmsd=mean_lig,
        std_ligand_rmsd=std_lig,
        max_ligand_rmsd=max_lig,
        drift_ligand_rmsd=drift,
        mean_protein_rmsd=mean_prot,
        std_protein_rmsd=std_prot,
        mean_contacts=mean_cnt,
        min_contacts=int(min_cnt) if min_cnt is not None else None,
        mean_mmgbsa_dg=mean_dg,
        std_mmgbsa_dg=std_dg,
        residue_persistence=res_persistence,
        stability_class=stability_class,
        stability_score=score,
        evidence_status=evidence_status,
        verdict_rationale=rationale,
    )


def generate_trajectory_svg(
    analysis: TrajectoryAnalysis,
    points: Sequence[TrajectoryPoint],
    width: int = 700,
    height: int = 320,
) -> str:
    """Generate a clean SVG chart displaying MD trajectory RMSD over simulation time."""
    if not points:
        return "<svg></svg>"

    padding_left = 60
    padding_right = 30
    padding_top = 40
    padding_bottom = 50

    plot_width = width - padding_left - padding_right
    plot_height = height - padding_top - padding_bottom

    times = [p.time_ns for p in points]
    min_t, max_t = min(times), max(times)
    t_span = max(1e-3, max_t - min_t)

    max_rmsd = max(5.0, max(p.ligand_rmsd for p in points) * 1.15)
    rmsd_span = max_rmsd

    def x_coord(t: float) -> float:
        return padding_left + ((t - min_t) / t_span) * plot_width

    def y_coord(r: float) -> float:
        clamped = max(0.0, min(max_rmsd, r))
        return padding_top + plot_height - (clamped / rmsd_span) * plot_height

    lig_points_str = " ".join(f"{x_coord(p.time_ns):.1f},{y_coord(p.ligand_rmsd):.1f}" for p in points)
    protein_points = [p for p in points if p.protein_rmsd is not None]
    prot_points_str = " ".join(
        f"{x_coord(p.time_ns):.1f},{y_coord(p.protein_rmsd or 0.0):.1f}" for p in protein_points
    )

    y_threshold = y_coord(2.0)

    if analysis.stability_class == "STABLE_POSE":
        status_color = "#10b981"
    elif analysis.stability_class == "MODERATE_MOBILITY":
        status_color = "#f59e0b"
    elif analysis.stability_class == "INSUFFICIENT_EVIDENCE":
        status_color = "#64748b"
    else:
        status_color = "#ef4444"

    grid_lines = []
    for tick in (1.0, 2.0, 3.0, 4.0):
        if tick < max_rmsd:
            yt = y_coord(tick)
            grid_lines.append(
                f'<line x1="{padding_left}" y1="{yt:.1f}" x2="{width - padding_right}" y2="{yt:.1f}" '
                'stroke="#e2e8f0" stroke-width="1" stroke-dasharray="3,3"/>'
                f'<text x="{padding_left - 8}" y="{yt + 4:.1f}" fill="#64748b" font-size="11" text-anchor="end">{tick:.1f}</text>'
            )

    time_ticks = []
    n_ticks = 4
    for i in range(n_ticks + 1):
        tick_t = min_t + (t_span / n_ticks) * i
        xt = x_coord(tick_t)
        time_ticks.append(
            f'<line x1="{xt:.1f}" y1="{height - padding_bottom}" x2="{xt:.1f}" y2="{height - padding_bottom + 5}" '
            'stroke="#64748b" stroke-width="1"/>'
            f'<text x="{xt:.1f}" y="{height - padding_bottom + 18}" fill="#64748b" font-size="11" text-anchor="middle">{tick_t:.0f}</text>'
        )

    svg_parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%" style="font-family: -apple-system, BlinkMacSystemFont, \'Segoe UI\', Roboto, sans-serif;">',
        f'<rect width="{width}" height="{height}" fill="#ffffff" rx="8" />',
        f'<text x="{padding_left}" y="24" fill="#0f172a" font-size="14" font-weight="600">{analysis.compound_id} on {analysis.target_id} - MD RMSD Trajectory</text>',
        f'<rect x="{width - padding_right - 150}" y="10" width="150" height="22" rx="11" fill="{status_color}22" stroke="{status_color}" stroke-width="1" />',
        f'<text x="{width - padding_right - 75}" y="25" fill="{status_color}" font-size="11" font-weight="600" text-anchor="middle">{analysis.stability_class}</text>',
        f'<line x1="{padding_left}" y1="{padding_top}" x2="{padding_left}" y2="{height - padding_bottom}" stroke="#94a3b8" stroke-width="1.5" />',
        f'<line x1="{padding_left}" y1="{height - padding_bottom}" x2="{width - padding_right}" y2="{height - padding_bottom}" stroke="#94a3b8" stroke-width="1.5" />',
        f'<text x="{width / 2}" y="{height - 12}" fill="#475569" font-size="12" font-weight="500" text-anchor="middle">Simulation Time (ns)</text>',
        f'<text transform="rotate(-90)" x="{-height / 2}" y="18" fill="#475569" font-size="12" font-weight="500" text-anchor="middle">RMSD (Å)</text>',
        "".join(grid_lines),
        "".join(time_ticks),
        f'<line x1="{padding_left}" y1="{y_threshold:.1f}" x2="{width - padding_right}" y2="{y_threshold:.1f}" stroke="#dc2626" stroke-width="1.5" stroke-dasharray="4,4" />',
        f'<text x="{width - padding_right - 4}" y="{y_threshold - 4:.1f}" fill="#dc2626" font-size="10" font-weight="500" text-anchor="end">Stability Gate (2.0 Å)</text>',
    ]
    if prot_points_str:
        svg_parts.append(f'<polyline fill="none" stroke="#8b5cf6" stroke-width="1.5" points="{prot_points_str}" opacity="0.8" />')
    svg_parts.extend(
        [
            f'<polyline fill="none" stroke="#0284c7" stroke-width="2.5" points="{lig_points_str}" />',
            f'<g transform="translate({padding_left + 10}, {padding_top + 10})">',
            '<rect x="0" y="0" width="12" height="12" fill="#0284c7" rx="2" />',
            f'<text x="18" y="10" fill="#334155" font-size="11">Ligand RMSD (Mean {analysis.mean_ligand_rmsd:.2f} Å)</text>',
        ]
    )
    if analysis.mean_protein_rmsd is not None:
        svg_parts.extend(
            [
                '<rect x="190" y="0" width="12" height="12" fill="#8b5cf6" rx="2" opacity="0.8" />',
                f'<text x="208" y="10" fill="#334155" font-size="11">Protein Backbone (Mean {analysis.mean_protein_rmsd:.2f} Å)</text>',
            ]
        )
    svg_parts.extend(["</g>", "</svg>"])
    return "\n".join(svg_parts)


def render_trajectory_report(analysis: TrajectoryAnalysis) -> str:
    """Render a Markdown report documenting the MD trajectory stability analysis."""
    persistence_rows = []
    if analysis.residue_persistence:
        for res, pct in analysis.residue_persistence.items():
            status = "Strong Contact" if pct >= 75.0 else ("Moderate Contact" if pct >= 50.0 else "Transient Contact")
            persistence_rows.append(f"| {res} | {pct:.1f}% | {status} |")
        persistence_table = (
            "| Binding Pocket Residue | Contact Persistence | Assessment |\n"
            "|:---|:---:|:---|\n" + "\n".join(persistence_rows)
        )
    else:
        persistence_table = "*No residue-specific contact columns were detected in trajectory input.*"

    def status(value: Optional[float], threshold: float, ok: str = "PASS", bad: str = "FAIL") -> str:
        if value is None:
            return "UNMEASURED"
        return ok if value < threshold else bad

    return f"""# Molecular Dynamics Trajectory Stability Report: {analysis.compound_id}

- **Target Receptor**: `{analysis.target_id}`
- **Simulation Duration**: `{analysis.total_time_ns:.1f} ns` (`{analysis.frames_count}` frames)
- **Evidence Status**: **`{analysis.evidence_status}`**
- **Pose Stability Class**: **`{analysis.stability_class}`**
- **Overall Stability Score**: **`{_display(analysis.stability_score, '.1f')}` / 100**

## 1. Key Metric Summary

| Trajectory Parameter | Observed Value | Recommended Threshold | Status |
|:---|:---:|:---:|:---|
| Ligand RMSD (Mean ± Std) | {analysis.mean_ligand_rmsd:.2f} ± {analysis.std_ligand_rmsd:.2f} Å | < 2.0 Å | {status(analysis.mean_ligand_rmsd, 2.0)} |
| Maximum Ligand RMSD | {analysis.max_ligand_rmsd:.2f} Å | < 3.2 Å | {status(analysis.max_ligand_rmsd, 3.2)} |
| RMSD Drift (End - Start) | {analysis.drift_ligand_rmsd:.2f} Å | < 1.0 Å | {status(analysis.drift_ligand_rmsd, 1.0)} |
| Protein Backbone RMSD | {_display(analysis.mean_protein_rmsd, '.2f')} Å | < 2.5 Å | {status(analysis.mean_protein_rmsd, 2.5)} |
| Intermolecular Contacts (Mean) | {_display(analysis.mean_contacts, '.1f')} contacts | ≥ 2.0 | {('UNMEASURED' if analysis.mean_contacts is None else ('PASS' if analysis.mean_contacts >= 2.0 else 'LOW'))} |
| Minimum Contact Count | {analysis.min_contacts if analysis.min_contacts is not None else 'UNMEASURED'} contacts | ≥ 1 | {('UNMEASURED' if analysis.min_contacts is None else ('PASS' if analysis.min_contacts >= 1 else 'UNBOUND'))} |
| MM-GBSA Binding Free Energy | {_display(analysis.mean_mmgbsa_dg, '.1f')} kcal/mol | < -25.0 kcal/mol | {('UNMEASURED' if analysis.mean_mmgbsa_dg is None else ('STRONG' if analysis.mean_mmgbsa_dg <= -25.0 else 'MODERATE'))} |

## 2. Pocket Residue Contact Persistence

{persistence_table}

## 3. Structural Stability Assessment

**Verdict Rationale**: {analysis.verdict_rationale}

> [!NOTE]
> In silico molecular dynamics simulations provide atomistic insights into pose stability and pocket residence time.
> Per guidelines G-001 and G-002, simulated kinetic metrics serve as mechanistic prioritization hypotheses
> and do not replace experimental biophysical validation (Surface Plasmon Resonance / ITC binding kinetics).
> UNMEASURED values were absent from the input trajectory and are never inferred.
"""


def run_dynamics_pipeline(
    trajectory_path: Path,
    compound_id: str = "PNN-101",
    target_id: str = "PTPsigma",
    out_dir: Optional[Path] = None,
) -> Dict[str, object]:
    """Execute the MD trajectory stability workflow, optionally writing reports and charts."""
    points = load_trajectory_csv(trajectory_path)
    analysis = analyze_trajectory(points, compound_id=compound_id, target_id=target_id)

    results: Dict[str, object] = {
        "analysis": analysis,
        "points_count": len(points),
    }

    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        csv_path = out_dir / f"{compound_id}_trajectory_metrics.csv"
        write_csv(csv_path, DYNAMICS_SUMMARY_FIELDS, [analysis.as_row()])
        results["csv_path"] = str(csv_path)

        json_path = out_dir / f"{compound_id}_dynamics_summary.json"
        write_json(json_path, analysis.as_dict())
        results["json_path"] = str(json_path)

        report_path = out_dir / f"{compound_id}_dynamics_report.md"
        write_text(report_path, render_trajectory_report(analysis))
        results["report_path"] = str(report_path)

        svg_path = out_dir / f"{compound_id}_rmsd_trajectory.svg"
        write_text(svg_path, generate_trajectory_svg(analysis, points))
        results["svg_path"] = str(svg_path)

    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.dynamics",
        description="Analyze molecular dynamics (MD) trajectories for ligand-target contact stability and RMSD drift.",
    )
    parser.add_argument("--trajectory", required=True, help="path to trajectory CSV file")
    parser.add_argument("--compound", default="PNN-101", help="compound identifier (default: PNN-101)")
    parser.add_argument("--target", default="PTPsigma", help="target receptor identifier (default: PTPsigma)")
    parser.add_argument("--out-dir", default="output/dynamics", help="directory to store analysis artifacts")
    parser.add_argument("--json", dest="json_path", help="optional additional destination path for the JSON summary")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        results = run_dynamics_pipeline(
            trajectory_path=Path(args.trajectory),
            compound_id=args.compound,
            target_id=args.target,
            out_dir=Path(args.out_dir),
        )

        analysis = results["analysis"]
        assert isinstance(analysis, TrajectoryAnalysis)

        if args.json_path:
            write_json(Path(args.json_path), analysis.as_dict())

    except DynamicsError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: MD trajectory analysis requires valid time and ligand RMSD columns; optional metrics stay unmeasured.", file=sys.stderr)
        print("Fix: verify that the trajectory CSV contains time_ns and ligand_rmsd headers.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: trajectory output files could not be written to disk.", file=sys.stderr)
        print("Fix: verify that output directory is writable.", file=sys.stderr)
        return 1

    print(
        f"[PNN] MD Dynamics analysis complete for {args.compound} on {args.target}: "
        f"class={analysis.stability_class}, evidence={analysis.evidence_status}, "
        f"mean_RMSD={analysis.mean_ligand_rmsd:.2f}A -> {args.out_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
