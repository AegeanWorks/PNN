"""Interactive Lead & Assay Visualization Studio and Standalone HTML Dashboard for PNN."""

from __future__ import annotations

import argparse
import csv
import html
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .reports import write_text
from .target import list_targets
from .waves import render_wave_section


def _get_val(obj: Any, key: str, default: Any = "-") -> Any:
    """Extract a field from either a dict/mapping or an object attribute with fallback."""
    if obj is None:
        return default
    if isinstance(obj, Mapping) or hasattr(obj, "get"):
        val = obj.get(key, default)
    else:
        val = getattr(obj, key, default)
    return default if val is None or val == "" else val


def _get_field(obj: Any, keys: Sequence[str], default: Any = "-") -> Any:
    """Extract first non-empty field matching one of several alias keys."""
    if obj is None:
        return default
    for k in keys:
        val = _get_val(obj, k, default=None)
        if val is not None and val != "" and val != "-":
            return val
    return default


class StudioError(ValueError):
    """Raised when studio inputs are malformed or missing required data."""


@dataclass(frozen=True)
class StudioCompound:
    compound_id: str
    name: str
    mw: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotb: int
    ec50_nm: Optional[float]
    brain_kp: Optional[float]
    oral_f_pct: Optional[float]
    docking_scores: Dict[str, float]
    pk_points: Dict[str, List[Tuple[float, float]]]  # route_tissue -> [(time, conc)]


RADAR_COLORS = ("#2563eb", "#047857", "#8b5cf6", "#b45309", "#ec4899")


def _normalize_cns_axis(val: float, optimal_center: float, half_width: float) -> float:
    diff = abs(val - optimal_center)
    score = max(0.1, 1.0 - (diff / half_width))
    return min(1.0, score)


def generate_radar_svg(compounds: Sequence[StudioCompound], width: int = 420, height: int = 400) -> str:
    """Generate SVG radar chart for CNS drug-likeness physicochemical properties."""
    cx, cy = width / 2.0, height / 2.0 + 10
    radius = min(cx, cy) - 55

    axes = [
        ("MW", 325.0, 150.0),      # optimal ~325 Da (range 175-475)
        ("logP", 2.2, 1.8),       # optimal ~2.2 (range 0.4-4.0)
        ("TPSA", 55.0, 35.0),     # optimal ~55 A^2 (range 20-90)
        ("HBD", 1.0, 2.0),        # optimal 0-2
        ("HBA", 4.0, 3.5),        # optimal 2-6
        ("RotB", 3.5, 4.0),       # optimal <=5
    ]
    n_axes = len(axes)
    angles = [(i * 2.0 * math.pi / n_axes) - (math.pi / 2.0) for i in range(n_axes)]

    svg_parts: List[str] = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="radar-title radar-desc">',
        '  <title id="radar-title">CNS physicochemical radar chart</title>',
        '  <desc id="radar-desc">Radar chart of MW, logP, TPSA, HBD, HBA and rotatable bonds for up to five '
        'compounds; a larger area is closer to the CNS-optimal window.</desc>',
        '  <defs>',
        '    <filter id="shadow" x="-10%" y="-10%" width="120%" height="120%">',
        '      <feDropShadow dx="0" dy="2" stdDeviation="3" flood-opacity="0.15"/>',
        '    </filter>',
        '  </defs>',
    ]

    # Background webs (concentric rings for 25%, 50%, 75%, 100%)
    for level in (0.25, 0.50, 0.75, 1.0):
        pts = [
            f"{cx + radius * level * math.cos(ang):.1f},{cy + radius * level * math.sin(ang):.1f}"
            for ang in angles
        ]
        svg_parts.append(f'  <polygon points="{" ".join(pts)}" fill="none" stroke="#e2e8f0" stroke-width="1.2"/>')

    # Spoke axes
    for i, (name, _, _) in enumerate(axes):
        ang = angles[i]
        x2 = cx + radius * math.cos(ang)
        y2 = cy + radius * math.sin(ang)
        svg_parts.append(f'  <line x1="{cx:.1f}" y1="{cy:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#cbd5e1" stroke-width="1.2"/>')
        lx = cx + (radius + 22) * math.cos(ang)
        ly = cy + (radius + 22) * math.sin(ang) + 4
        svg_parts.append(
            f'  <text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" font-size="11" font-weight="600" fill="#475569">{name}</text>'
        )

    # Plot compounds
    colors = RADAR_COLORS
    for idx, comp in enumerate(compounds[:5]):
        color = colors[idx % len(colors)]
        ratios = [
            _normalize_cns_axis(comp.mw, axes[0][1], axes[0][2]),
            _normalize_cns_axis(comp.logp, axes[1][1], axes[1][2]),
            _normalize_cns_axis(comp.tpsa, axes[2][1], axes[2][2]),
            _normalize_cns_axis(float(comp.hbd), axes[3][1], axes[3][2]),
            _normalize_cns_axis(float(comp.hba), axes[4][1], axes[4][2]),
            _normalize_cns_axis(float(comp.rotb), axes[5][1], axes[5][2]),
        ]
        poly_pts = [
            f"{cx + radius * ratios[i] * math.cos(angles[i]):.1f},{cy + radius * ratios[i] * math.sin(angles[i]):.1f}"
            for i in range(n_axes)
        ]
        pts_str = " ".join(poly_pts)
        svg_parts.append(
            f'  <polygon points="{pts_str}" fill="{color}" fill-opacity="0.22" stroke="{color}" stroke-width="2.2" filter="url(#shadow)"/>'
        )
        for pt in poly_pts:
            px, py = pt.split(",")
            svg_parts.append(f'  <circle cx="{px}" cy="{py}" r="3.5" fill="{color}"/>')

    svg_parts.append('</svg>')
    return "\n".join(svg_parts)


def generate_docking_heatmap_svg(
    compounds: Sequence[StudioCompound],
    targets: Sequence[str],
    cell_w: int = 75,
    cell_h: int = 36,
) -> str:
    """Generate multi-target docking affinity heatmap in SVG."""
    margin_left = 110
    margin_top = 40
    width = margin_left + len(targets) * cell_w + 30
    height = margin_top + len(compounds) * cell_h + 30

    svg_parts: List[str] = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="heatmap-title heatmap-desc">',
        '  <title id="heatmap-title">Multi-target docking affinity heatmap</title>',
        '  <desc id="heatmap-desc">Docking scores in kcal/mol for each compound and target; darker cells '
        'mark stronger predicted binding.</desc>',
    ]

    # Target column headers
    for j, tgt in enumerate(targets):
        x = margin_left + j * cell_w + cell_w / 2.0
        svg_parts.append(
            f'  <text x="{x:.1f}" y="{margin_top - 12}" text-anchor="middle" font-size="11" font-weight="600" fill="#334155">{tgt}</text>'
        )

    # Heatmap rows
    for i, comp in enumerate(compounds):
        y = margin_top + i * cell_h
        # Row label
        svg_parts.append(
            f'  <text x="{margin_left - 12}" y="{y + cell_h/2.0 + 4}" text-anchor="end" font-size="11" font-weight="600" fill="#1e293b">{comp.compound_id}</text>'
        )

        for j, tgt in enumerate(targets):
            x = margin_left + j * cell_w
            score = comp.docking_scores.get(tgt)
            if score is not None:
                # Normalizing -5.0 to -9.5 kcal/mol to [0, 1]
                ratio = max(0.0, min(1.0, (abs(score) - 5.0) / 4.5))
                # Hue from light cyan (190, 60%, 90%) to deep navy (215, 85%, 28%)
                lightness = int(92 - 62 * ratio)
                color = f"hsl(215, 80%, {lightness}%)"
                text_color = "#ffffff" if lightness < 52 else "#000000"
                val_text = f"{score:.1f}"
            else:
                color = "#f1f5f9"
                text_color = "#475569"
                val_text = "N/A"

            svg_parts.append(
                f'  <rect x="{x}" y="{y}" width="{cell_w - 3}" height="{cell_h - 3}" rx="4" fill="{color}"/>'
            )
            svg_parts.append(
                f'  <text x="{x + cell_w/2.0 - 1.5}" y="{y + cell_h/2.0 + 4}" text-anchor="middle" font-size="11" font-weight="500" fill="{text_color}">{val_text}</text>'
            )

    svg_parts.append('</svg>')
    return "\n".join(svg_parts)


def generate_pk_curve_svg(
    compounds: Sequence[StudioCompound],
    width: int = 540,
    height: int = 300,
) -> str:
    """Generate in vivo pharmacokinetic concentration-time profiles SVG."""
    ml, mr, mt, mb = 65, 30, 30, 45
    pw = width - ml - mr
    ph = height - mt - mb

    # Find max time and max concentration across compounds
    all_pts: List[Tuple[float, float]] = []
    for c in compounds:
        for pts in c.pk_points.values():
            all_pts.extend(pts)

    max_t = 24.0
    max_c = max([p[1] for p in all_pts], default=1000.0)
    max_c = max(max_c, 500.0)

    def tx(t: float) -> float:
        return ml + (t / max_t) * pw

    def ty(c: float) -> float:
        # Linear or semi-log mapping
        c_clamped = max(0.0, min(c, max_c))
        return mt + ph - (c_clamped / max_c) * ph

    svg_parts: List[str] = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="pk-title pk-desc">',
        '  <title id="pk-title">In vivo pharmacokinetic concentration-time curves</title>',
        '  <desc id="pk-desc">Plasma and brain concentration over 24 hours for the lead compound after '
        'intravenous and oral dosing.</desc>',
    ]

    # Gridlines
    for step in range(5):
        val = (max_c / 4.0) * step
        y_pos = ty(val)
        svg_parts.append(f'  <line x1="{ml}" y1="{y_pos:.1f}" x2="{ml + pw}" y2="{y_pos:.1f}" stroke="#e2e8f0" stroke-width="1"/>')
        svg_parts.append(f'  <text x="{ml - 8}" y="{y_pos + 4:.1f}" text-anchor="end" font-size="10" fill="#475569">{val:.0f}</text>')

    for t_step in (0, 4, 8, 12, 16, 20, 24):
        x_pos = tx(float(t_step))
        svg_parts.append(f'  <line x1="{x_pos:.1f}" y1="{mt}" x2="{x_pos:.1f}" y2="{mt + ph}" stroke="#e2e8f0" stroke-width="1"/>')
        svg_parts.append(f'  <text x="{x_pos:.1f}" y="{mt + ph + 18}" text-anchor="middle" font-size="10" fill="#475569">{t_step}h</text>')

    # Axes
    svg_parts.append(f'  <line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt + ph}" stroke="#94a3b8" stroke-width="1.5"/>')
    svg_parts.append(f'  <line x1="{ml}" y1="{mt + ph}" x2="{ml + pw}" y2="{mt + ph}" stroke="#94a3b8" stroke-width="1.5"/>')
    svg_parts.append(f'  <text x="{ml - 14}" y="{mt - 10}" text-anchor="start" font-size="10" font-weight="600" fill="#475569">Conc (ng/mL)</text>')
    svg_parts.append(f'  <text x="{ml + pw/2}" y="{height - 6}" text-anchor="middle" font-size="10" font-weight="600" fill="#475569">Time post-dose (hours)</text>')

    # Plot curves for top candidate
    styles = {
        "IV_plasma": ("#2563eb", "none", "2.5"),        # Blue solid
        "PO_plasma": ("#047857", "6,3", "2.5"),         # Green dashed
        "PO_brain": ("#b45309", "2,3", "2.5"),          # Amber dotted
    }

    if compounds:
        lead = compounds[0]
        for route_tissue, pts in sorted(lead.pk_points.items()):
            color, dash, stroke_w = styles.get(route_tissue, ("#475569", "none", "1.5"))
            dash_attr = f' stroke-dasharray="{dash}"' if dash != "none" else ""
            sorted_pts = sorted(pts, key=lambda p: p[0])
            poly_coords = [f"{tx(p[0]):.1f},{ty(p[1]):.1f}" for p in sorted_pts]
            pts_str = " ".join(poly_coords)
            svg_parts.append(f'  <polyline points="{pts_str}" fill="none" stroke="{color}" stroke-width="{stroke_w}"{dash_attr}/>')
            for p in sorted_pts:
                svg_parts.append(f'  <circle cx="{tx(p[0]):.1f}" cy="{ty(p[1]):.1f}" r="3.2" fill="{color}"/>')

    svg_parts.append('</svg>')
    return "\n".join(svg_parts)


def _safe_cell(value: object) -> str:
    """Return a stripped string for a CSV cell that may be None."""
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    return str(value).strip()


def _parse_required_float(raw: str, field: str, compound_id: str) -> float:
    """Parse a required numeric studio field, raising StudioError on missing data."""
    parsed = _parse_opt_float(raw)
    if parsed is None:
        raise StudioError(f"compound {compound_id} has missing or non-numeric {field}: {raw!r}")
    return parsed


def _parse_required_int(raw: str, field: str, compound_id: str) -> int:
    """Parse a required integer studio field, raising StudioError on missing data."""
    parsed = _parse_opt_float(raw)
    if parsed is None:
        raise StudioError(f"compound {compound_id} has missing or non-numeric {field}: {raw!r}")
    return int(parsed)


def _parse_circuit_points(
    simulation_points: Sequence[Mapping[str, object]],
) -> List[Tuple[float, float, float]]:
    """(time_h, pnn_pct, akt_pct) for every complete, numeric simulation row, time-sorted."""
    parsed: List[Tuple[float, float, float]] = []
    for pt in simulation_points:
        try:
            t_raw = _get_field(pt, ["time_h", "time", "t_h"], None)
            pnn_raw = _get_field(pt, ["pnn_density_pct", "pnn_density", "pnn"], None)
            akt_raw = _get_field(pt, ["akt_signaling_pct", "akt_signaling", "akt"], None)
            if t_raw is None or pnn_raw is None or akt_raw is None:
                continue
            parsed.append((float(t_raw), float(pnn_raw), float(akt_raw)))
        except (ValueError, TypeError):
            continue
    parsed.sort(key=lambda point: point[0])
    return parsed


def generate_circuit_dynamics_svg(
    simulation_points: Sequence[Mapping[str, object]],
    width: int = 540,
    height: int = 340,
) -> str:
    """Generate SVG time-course chart for hypothalamic circuit dynamics (PNN turnover vs pAKT restoration)."""
    ml, mr, mt, mb = 60, 25, 30, 45
    pw = width - ml - mr
    ph = height - mt - mb

    if not simulation_points:
        return (
            f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
            f'xmlns="http://www.w3.org/2000/svg" role="img" '
            f'aria-labelledby="circuit-title circuit-desc">'
            f'<title id="circuit-title">Hypothalamic circuit dosing dynamics</title>'
            f'<desc id="circuit-desc">No simulation data was provided.</desc>'
            f'<text x="{width/2:.1f}" y="{height/2:.1f}" text-anchor="middle" fill="#94a3b8">No simulation data</text></svg>'
        )

    parsed_pts = _parse_circuit_points(simulation_points)

    if not parsed_pts:
        return (
            f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
            f'xmlns="http://www.w3.org/2000/svg" role="img" '
            f'aria-labelledby="circuit-title circuit-desc">'
            f'<title id="circuit-title">Hypothalamic circuit dosing dynamics</title>'
            f'<desc id="circuit-desc">No valid simulation points were provided.</desc>'
            f'<text x="{width/2:.1f}" y="{height/2:.1f}" text-anchor="middle" fill="#94a3b8">No valid simulation points</text></svg>'
        )

    max_t = max(parsed_pts[-1][0], 24.0)

    def tx(t: float) -> float:
        return ml + (t / max_t) * pw

    def ty(pct: float) -> float:
        clamped = max(0.0, min(pct, 120.0))
        return mt + ph - (clamped / 120.0) * ph

    svg_parts: List[str] = [
        f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="circuit-title circuit-desc">',
        '  <title id="circuit-title">Hypothalamic circuit dosing dynamics</title>',
        '  <desc id="circuit-desc">PNN density and pAKT signaling over the simulation window; the red '
        'line marks the 30% safety floor.</desc>',
    ]

    for y_val in (0, 25, 50, 75, 100):
        y_pos = ty(float(y_val))
        svg_parts.append(
            f'  <line x1="{ml}" y1="{y_pos:.1f}" x2="{ml + pw}" y2="{y_pos:.1f}" stroke="#e2e8f0" stroke-width="1"/>'
        )
        svg_parts.append(
            f'  <text x="{ml - 8}" y="{y_pos + 4:.1f}" text-anchor="end" font-size="10" fill="#475569">{y_val}%</text>'
        )

    n_intervals = 7
    for step in range(n_intervals + 1):
        t_val = (max_t / n_intervals) * step
        x_pos = tx(t_val)
        svg_parts.append(
            f'  <line x1="{x_pos:.1f}" y1="{mt}" x2="{x_pos:.1f}" y2="{mt + ph}" stroke="#e2e8f0" stroke-width="1"/>'
        )
        label = f"{t_val/24:.0f}d" if max_t >= 48 else f"{t_val:.0f}h"
        svg_parts.append(
            f'  <text x="{x_pos:.1f}" y="{mt + ph + 18}" text-anchor="middle" font-size="10" fill="#475569">{label}</text>'
        )

    safety_y = ty(30.0)
    svg_parts.append(
        f'  <line x1="{ml}" y1="{safety_y:.1f}" x2="{ml + pw}" y2="{safety_y:.1f}" stroke="#ef4444" stroke-width="1.5" stroke-dasharray="4,4"/>'
    )
    svg_parts.append(
        f'  <text x="{ml + pw - 6}" y="{safety_y - 6:.1f}" text-anchor="end" font-size="9" fill="#b91c1c" font-weight="600">30% PNN Safety Floor</text>'
    )

    svg_parts.append(f'  <line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt + ph}" stroke="#94a3b8" stroke-width="1.5"/>')
    svg_parts.append(f'  <line x1="{ml}" y1="{mt + ph}" x2="{ml + pw}" y2="{mt + ph}" stroke="#94a3b8" stroke-width="1.5"/>')
    svg_parts.append(
        f'  <text x="{ml - 14}" y="{mt - 10}" text-anchor="start" font-size="10" font-weight="600" fill="#475569">% of Physiological Baseline</text>'
    )
    svg_parts.append(
        f'  <text x="{ml + pw/2}" y="{height - 6}" text-anchor="middle" font-size="10" font-weight="600" fill="#475569">Simulation Duration (Days post-initiation)</text>'
    )

    pnn_coords = [f"{tx(t):.1f},{ty(pnn):.1f}" for t, pnn, _ in parsed_pts]
    svg_parts.append(
        f'  <polyline points="{" ".join(pnn_coords)}" fill="none" stroke="#8b5cf6" stroke-width="2.5"/>'
    )

    akt_coords = [f"{tx(t):.1f},{ty(akt):.1f}" for t, _, akt in parsed_pts]
    svg_parts.append(
        f'  <polyline points="{" ".join(akt_coords)}" fill="none" stroke="#0e7490" stroke-width="2.5" stroke-dasharray="6,3"/>'
    )

    svg_parts.append(
        f'  <line x1="{ml + 10}" y1="{mt + 10}" x2="{ml + 35}" y2="{mt + 10}" stroke="#8b5cf6" stroke-width="2.5"/>'
    )
    svg_parts.append(
        f'  <text x="{ml + 40}" y="{mt + 14}" font-size="10" font-weight="600" fill="#475569">PNN Density (%)</text>'
    )
    svg_parts.append(
        f'  <line x1="{ml + 140}" y1="{mt + 10}" x2="{ml + 165}" y2="{mt + 10}" stroke="#0e7490" stroke-width="2.5" stroke-dasharray="6,3"/>'
    )
    svg_parts.append(
        f'  <text x="{ml + 170}" y="{mt + 14}" font-size="10" font-weight="600" fill="#475569">pAKT Signaling (%)</text>'
    )

    svg_parts.append("</svg>")
    return "\n".join(svg_parts)


def render_studio_html(
    compounds: Sequence[StudioCompound],
    radar_svg: str,
    heatmap_svg: str,
    pk_svg: str,
    circuit_svg: Optional[str] = None,
    circuit_points: Optional[Sequence[Mapping[str, object]]] = None,
    dossier_data: Optional[Sequence[Mapping[str, object]]] = None,
      safety_data: Optional[Sequence[Mapping[str, object]]] = None,
      targets_data: Optional[Sequence[object]] = None,
      wave_ledger: Optional[Sequence[Mapping[str, object]]] = None,
      wave_highlights: Optional[Sequence[Mapping[str, object]]] = None,
  ) -> str:
    """Render full standalone HTML5 studio dashboard."""
    cards_html = []
    colors = ["#3b82f6", "#10b981", "#8b5cf6", "#f59e0b", "#ec4899"]
    for idx, c in enumerate(compounds[:5]):
        color = colors[idx % len(colors)]
        ec50_str = f"{c.ec50_nm:.1f} nM" if c.ec50_nm is not None else "N/A"
        kp_str = f"{c.brain_kp:.2f}" if c.brain_kp is not None else "N/A"
        f_str = f"{c.oral_f_pct:.1f}%" if c.oral_f_pct is not None else "N/A"
        badge_cls = "badge-candidate" if (c.brain_kp and c.brain_kp >= 0.30 and c.oral_f_pct and c.oral_f_pct >= 20) else "badge-lead"
        badge_label = "Clinical Candidate" if badge_cls == "badge-candidate" else "Advanced Lead"

        cards_html.append(f"""
        <div class="lead-card" style="border-top: 4px solid {color};">
          <div class="lead-header">
            <h3>{html.escape(c.compound_id)} <span class="lead-name">({html.escape(c.name)})</span></h3>
            <span class="badge {badge_cls}">{badge_label}</span>
          </div>
          <div class="lead-grid">
            <div class="metric"><span class="m-lbl">DQ-CS EC50</span><span class="m-val">{ec50_str}</span></div>
            <div class="metric"><span class="m-lbl">Brain Kp</span><span class="m-val">{kp_str}</span></div>
            <div class="metric"><span class="m-lbl">Oral F</span><span class="m-val">{f_str}</span></div>
            <div class="metric"><span class="m-lbl">MW / logP</span><span class="m-val">{c.mw:.0f} / {c.logp:.1f}</span></div>
          </div>
        </div>
        """)

    dossier_section = ""
    if dossier_data:
        dossier_rows = []
        for row in dossier_data:
            cid = html.escape(str(_get_field(row, ["compound_id", "id"], "-")))
            chem = html.escape(str(_get_field(row, ["chemotype", "name"], "-")))
            evid = html.escape(str(_get_field(row, ["evidence_status", "evidence"], "-")))
            evid_cls = "badge-safe" if evid == "COMPLETE" else "badge-caution"

            dose_val = _get_field(row, ["animal_dose_mg_kg", "animal_dose"], "-")
            dose = f"{html.escape(str(dose_val))} mg/kg" if dose_val != "-" else "-"

            hdose_val = _get_field(row, ["human_total_daily_dose_mg", "human_dose_mg_day", "human_dose"], "-")
            hdose = f"{html.escape(str(hdose_val))} mg" if hdose_val != "-" else "-"

            ht12_val = _get_field(row, ["projected_human_half_life_h", "human_t12_h", "human_t12"], "-")
            ht12 = f"{html.escape(str(ht12_val))} h" if ht12_val != "-" else "-"

            mrsd_val = _get_field(row, ["phase1_mrsd_mg", "mrsd_mg", "mrsd"], "-")
            mrsd = f"{html.escape(str(mrsd_val))} mg" if mrsd_val != "-" else "-"

            margin_raw = str(_get_field(row, ["therapeutic_index", "safety_margin"], "-"))
            if margin_raw != "-":
                margin = html.escape(margin_raw if margin_raw.endswith("x") else f"{margin_raw}x")
            else:
                margin = "-"

            dossier_rows.append(f"""
            <tr>
              <td><strong>{cid}</strong></td>
              <td>{chem}</td>
              <td><span class="{evid_cls}">{evid}</span></td>
              <td>{dose}</td>
              <td>{hdose}</td>
              <td>{ht12}</td>
              <td>{mrsd}</td>
              <td><strong>{margin}</strong></td>
            </tr>""")
        dossier_section = f"""
    <div class="panel" style="margin-bottom: 24px;">
      <h2>Investigational New Drug (IND) & Translational Dossier <span class="badge badge-candidate">ICH M3(R2) / S7A / S7B</span></h2>
      <div class="table-wrapper">
        <table class="data-table">
          <caption class="sr-only">IND dose projection, half-life and safety margins for each candidate compound</caption>
          <thead>
            <tr>
              <th scope="col">Compound ID</th>
              <th scope="col">Chemotype</th>
              <th scope="col">Evidence</th>
              <th scope="col">Animal Dose</th>
              <th scope="col">Projected Human Dose</th>
              <th scope="col">Human t1/2</th>
              <th scope="col">Phase 1 MRSD</th>
              <th scope="col">Safety Margin</th>
            </tr>
          </thead>
          <tbody>
            {"".join(dossier_rows)}
          </tbody>
        </table>
      </div>
    </div>
        """

    safety_section = ""
    if safety_data:
        safety_rows = []
        for row in safety_data:
            cid = html.escape(str(_get_field(row, ["compound_id", "id"], "-")))
            stat = html.escape(str(_get_field(row, ["status"], "-")))
            stat_cls = "badge-safe" if stat == "SAFE_TO_PROGRESS" else ("badge-fail" if stat == "FAIL_SAFETY" else "badge-caution")
            spec = html.escape(str(_get_field(row, ["most_sensitive_species", "species"], "-")))

            noael_val = _get_field(row, ["noael_mg_kg", "most_sensitive_noael_mg_kg", "noael"], "-")
            noael = f"{html.escape(str(noael_val))} mg/kg" if noael_val != "-" else "-"

            herg_raw = str(_get_field(row, ["herg_margin", "herg"], "-"))
            if herg_raw != "-":
                herg = html.escape(herg_raw if herg_raw.endswith("x") else f"{herg_raw}x")
            else:
                herg = "-"

            cyp_raw = str(_get_field(row, ["cyp3a4_inhibition_pct", "cyp3a4_inhibition", "cyp"], "-"))
            if cyp_raw != "-":
                cyp = html.escape(cyp_raw if cyp_raw.endswith("%") else f"{cyp_raw}%")
            else:
                cyp = "-"

            geno = html.escape(str(_get_field(row, ["genotoxicity"], "-")))
            geno_cls = "badge-safe" if geno == "NEGATIVE" else ("badge-fail" if geno == "POSITIVE" else "badge-caution")
            safety_rows.append(f"""
            <tr>
              <td><strong>{cid}</strong></td>
              <td><span class="{stat_cls}">{stat}</span></td>
              <td>{spec}</td>
              <td>{noael}</td>
              <td>{herg}</td>
              <td>{cyp}</td>
              <td><span class="{geno_cls}">{geno}</span></td>
            </tr>""")
        safety_section = f"""
    <div class="panel" style="margin-bottom: 24px;">
      <h2>IND-Enabling Safety Pharmacology & Repeat-Dose Toxicology <span class="badge badge-lead">ICH S7A/S7B Panel</span></h2>
      <div class="table-wrapper">
        <table class="data-table">
          <caption class="sr-only">Safety pharmacology and repeat-dose toxicology outcomes for each compound</caption>
          <thead>
            <tr>
              <th scope="col">Compound ID</th>
              <th scope="col">Status</th>
              <th scope="col">Sensitive Species</th>
              <th scope="col">NOAEL</th>
              <th scope="col">hERG Margin</th>
              <th scope="col">CYP3A4 Inh.</th>
              <th scope="col">Genotoxicity</th>
            </tr>
          </thead>
          <tbody>
            {"".join(safety_rows)}
          </tbody>
        </table>
      </div>
    </div>
        """

    wave_section = ""
    if wave_ledger:
        wave_section = render_wave_section(wave_ledger, wave_highlights or [])

    targets_section = ""
    if targets_data:
        target_rows = []
        for t in targets_data:
            tid = html.escape(str(_get_field(t, ["target_id", "id"], "-")))
            gene = html.escape(str(_get_field(t, ["gene_symbol", "gene"], "-")))
            pdb = html.escape(str(_get_field(t, ["primary_pdb", "pdb"], "-")))
            res_val = _get_field(t, ["resolution_angstrom", "resolution"], "-")
            res = f"{html.escape(str(res_val))} &Aring;" if res_val != "-" else "-"
            ptype = html.escape(str(_get_field(t, ["pocket_type", "pocket"], "-")))
            drug = html.escape(str(_get_field(t, ["druggability"], "-")))
            bio = html.escape(str(_get_field(t, ["biological_function", "function"], "-")))
            target_rows.append(f"""
            <tr>
              <td><strong>{tid}</strong></td>
              <td><code>{gene}</code></td>
              <td><strong>{pdb}</strong></td>
              <td>{res}</td>
              <td>{ptype}</td>
              <td>{drug}</td>
              <td style="font-size:0.8rem; color:#94a3b8;">{bio}</td>
            </tr>""")
        targets_section = f"""
    <div class="panel" style="margin-bottom: 24px;">
      <h2>Receptor & Aggrecanase Target Invariants <span class="badge badge-candidate">{len(targets_data)} Registered Targets</span></h2>
      <div class="table-wrapper">
        <table class="data-table">
          <caption class="sr-only">Registered receptor and aggrecanase target structural invariants</caption>
          <thead>
            <tr>
              <th scope="col">Target ID</th>
              <th scope="col">Gene Symbol</th>
              <th scope="col">Primary PDB</th>
              <th scope="col">Resolution</th>
              <th scope="col">Pocket Type</th>
              <th scope="col">Druggability</th>
              <th scope="col">Biological Role</th>
            </tr>
          </thead>
          <tbody>
            {"".join(target_rows)}
          </tbody>
        </table>
      </div>
    </div>
        """

    circuit_table = ""
    if circuit_svg and circuit_points:
        circuit_rows = "".join(
            f'<tr><th scope="row">{time_h:g}</th><td>{pnn:g}</td><td>{akt:g}</td></tr>'
            for time_h, pnn, akt in _parse_circuit_points(circuit_points)
        )
        if circuit_rows:
            circuit_table = f"""
        <table class="sr-only">
          <caption>Circuit dynamics data: PNN density and pAKT signaling over time</caption>
          <thead><tr><th scope="col">Time (h)</th><th scope="col">PNN density (%)</th><th scope="col">pAKT signaling (%)</th></tr></thead>
          <tbody>{circuit_rows}</tbody>
        </table>"""

    circuit_section = ""
    if circuit_svg:
        circuit_section = f"""
    <div class="panel" style="margin-bottom: 24px;">
      <h2>Hypothalamic Circuit Dosing Dynamics <span style="font-size:0.8rem; font-weight:400; color:var(--text-muted);">Purple solid: PNN Density (%) | Cyan dashed: pAKT Signaling (%) | Red dashed: 30% Safety Floor</span></h2>
      <div class="svg-wrapper">
        {circuit_svg}
      </div>
      {circuit_table}
    </div>
        """

    radar_rows = []
    for compound in compounds[:5]:
        radar_rows.append(
            f'<tr><th scope="row">{html.escape(compound.compound_id)}</th>'
            f'<td>{compound.mw:.1f}</td><td>{compound.logp:.2f}</td>'
            f'<td>{compound.tpsa:.1f}</td><td>{compound.hbd}</td>'
            f'<td>{compound.hba}</td><td>{compound.rotb}</td></tr>'
        )
    radar_table = ""
    if radar_rows:
        radar_table = f"""
        <table class="sr-only">
          <caption>Radar chart data: CNS physicochemical properties</caption>
          <thead><tr><th scope="col">Compound</th><th scope="col">MW (Da)</th><th scope="col">logP</th><th scope="col">TPSA (A^2)</th><th scope="col">HBD</th><th scope="col">HBA</th><th scope="col">RotB</th></tr></thead>
          <tbody>{"".join(radar_rows)}</tbody>
        </table>"""

    radar_legend = ""
    if compounds:
        legend_items = "".join(
            f'<li><span class="legend-swatch" style="background:{RADAR_COLORS[index % len(RADAR_COLORS)]}"></span>'
            f'{html.escape(compound.compound_id)}</li>'
            for index, compound in enumerate(compounds[:5])
        )
        radar_legend = f'<ul class="chart-legend">{legend_items}</ul>'

    docking_targets = sorted({target for c in compounds for target in c.docking_scores})
    heatmap_rows = []
    for compound in compounds:
        cells = "".join(
            f'<td>{compound.docking_scores[target]:.1f}</td>'
            if compound.docking_scores.get(target) is not None else '<td>N/A</td>'
            for target in docking_targets
        )
        heatmap_rows.append(
            f'<tr><th scope="row">{html.escape(compound.compound_id)}</th>{cells}</tr>'
        )
    heatmap_table = ""
    if docking_targets:
        target_headers = "".join(
            f'<th scope="col">{html.escape(target)}</th>' for target in docking_targets
        )
        heatmap_table = f"""
        <table class="sr-only">
          <caption>Docking heatmap data: predicted binding affinity (kcal/mol)</caption>
          <thead><tr><th scope="col">Compound</th>{target_headers}</tr></thead>
          <tbody>{"".join(heatmap_rows)}</tbody>
        </table>"""

    pk_rows = []
    if compounds:
        lead = compounds[0]
        for route, points in sorted(lead.pk_points.items()):
            route_label = route.replace("_", " ")
            for time_h, concentration in sorted(points):
                pk_rows.append(
                    f'<tr><th scope="row">{html.escape(lead.compound_id)}</th>'
                    f'<td>{html.escape(route_label)}</td>'
                    f'<td>{time_h:g}</td><td>{concentration:g}</td></tr>'
                )
    pk_table = ""
    if pk_rows:
        pk_table = f"""
        <table class="sr-only">
          <caption>Pharmacokinetic chart data: concentration over time</caption>
          <thead><tr><th scope="col">Compound</th><th scope="col">Route</th><th scope="col">Time (h)</th><th scope="col">Concentration (ng/mL)</th></tr></thead>
          <tbody>{"".join(pk_rows)}</tbody>
        </table>"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>PNN Preclinical & Lead Discovery Studio</title>
  <style>
    :root {{
      --bg: #0f172a;
      --surface: #1e293b;
      --surface-border: #334155;
      --text: #f8fafc;
      --text-muted: #94a3b8;
      --primary: #38bdf8;
      --success: #34d399;
      --warning: #fbbf24;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
      padding: 24px;
      line-height: 1.5;
    }}
    .container {{
      max-width: 1300px;
      margin: 0 auto;
    }}
    header {{
      margin-bottom: 24px;
      padding-bottom: 16px;
      border-bottom: 1px solid var(--surface-border);
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 12px;
    }}
    h1 {{ font-size: 1.6rem; font-weight: 700; color: #ffffff; }}
    .subtitle {{ font-size: 0.9rem; color: var(--text-muted); }}
    .funnel-container {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 14px;
      margin-bottom: 24px;
    }}
    .funnel-card {{
      background: var(--surface);
      border: 1px solid var(--surface-border);
      border-radius: 8px;
      padding: 14px;
      text-align: center;
    }}
    .funnel-step {{ font-size: 0.75rem; text-transform: uppercase; color: var(--primary); font-weight: 700; }}
    .funnel-title {{ font-size: 0.95rem; font-weight: 600; margin: 4px 0; }}
    .funnel-num {{ font-size: 1.6rem; font-weight: 700; color: #ffffff; }}
    .grid-2col {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(min(480px, 100%), 1fr));
      gap: 20px;
      margin-bottom: 24px;
    }}
    .panel {{
      background: var(--surface);
      border: 1px solid var(--surface-border);
      border-radius: 8px;
      padding: 18px;
    }}
    .panel h2 {{
      font-size: 1.1rem;
      font-weight: 600;
      margin-bottom: 12px;
      color: #f1f5f9;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }}
    .svg-wrapper {{
      display: flex;
      flex-direction: column;
      justify-content: center;
      align-items: center;
      background: #ffffff;
      border-radius: 6px;
      padding: 12px;
      overflow-x: auto;
    }}
    .svg-wrapper svg {{
      max-width: 100%;
      height: auto;
    }}
    .chart-legend {{
      display: flex;
      flex-wrap: wrap;
      justify-content: center;
      gap: 12px;
      margin: 8px 0 0;
      padding: 0;
      list-style: none;
      font-size: 0.8rem;
      color: #334155;
    }}
    .legend-swatch {{
      display: inline-block;
      width: 12px;
      height: 12px;
      border-radius: 3px;
      margin-right: 5px;
      vertical-align: middle;
    }}
    .cards-row {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(min(280px, 100%), 1fr));
      gap: 16px;
      margin-bottom: 24px;
    }}
    .lead-card {{
      background: var(--surface);
      border: 1px solid var(--surface-border);
      border-radius: 8px;
      padding: 14px;
    }}
    .lead-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 10px;
    }}
    .lead-header h3 {{ font-size: 1.05rem; font-weight: 600; }}
    .lead-name {{ font-size: 0.8rem; color: var(--text-muted); font-weight: 400; }}
    .badge {{
      font-size: 0.72rem;
      font-weight: 700;
      padding: 3px 8px;
      border-radius: 9999px;
      text-transform: uppercase;
    }}
    .badge-candidate {{ background: #065f46; color: #6ee7b7; }}
    .badge-lead {{ background: #1e3a8a; color: #93c5fd; }}
    .badge-safe {{ background: #065f46; color: #6ee7b7; font-weight: 700; font-size: 0.72rem; padding: 2px 7px; border-radius: 9999px; text-transform: uppercase; }}
    .badge-caution {{ background: #78350f; color: #fde68a; font-weight: 700; font-size: 0.72rem; padding: 2px 7px; border-radius: 9999px; text-transform: uppercase; }}
    .badge-fail {{ background: #7f1d1d; color: #fca5a5; font-weight: 700; font-size: 0.72rem; padding: 2px 7px; border-radius: 9999px; text-transform: uppercase; }}
    .lead-grid {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 8px;
    }}
    .metric {{
      display: flex;
      flex-direction: column;
      background: #0f172a;
      padding: 6px 10px;
      border-radius: 4px;
    }}
    .m-lbl {{ font-size: 0.7rem; color: var(--text-muted); }}
    .m-val {{ font-size: 0.88rem; font-weight: 600; color: #f8fafc; }}
    .table-wrapper {{
      overflow-x: auto;
      margin-top: 10px;
    }}
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.85rem;
      text-align: left;
    }}
    .data-table th {{
      background: #0f172a;
      color: var(--primary);
      padding: 10px 12px;
      font-weight: 600;
      border-bottom: 2px solid var(--surface-border);
      white-space: nowrap;
    }}
    .data-table td {{
      padding: 9px 12px;
      border-bottom: 1px solid var(--surface-border);
      color: #cbd5e1;
    }}
    .data-table tr:hover {{
      background: rgba(255, 255, 255, 0.03);
    }}
    .sr-only {{
      position: absolute;
      width: 1px;
      height: 1px;
      padding: 0;
      margin: -1px;
      overflow: hidden;
      clip: rect(0, 0, 0, 0);
      white-space: nowrap;
      border: 0;
    }}
    footer {{
      margin-top: 24px;
      padding-top: 16px;
      border-top: 1px solid var(--surface-border);
      font-size: 0.8rem;
      color: var(--text-muted);
      text-align: center;
    }}
    .skip-link {{
      position: absolute;
      left: -9999px;
      top: 0;
      z-index: 100;
      background: #ffffff;
      color: #0f172a;
      padding: 10px 14px;
      border-radius: 0 0 6px 0;
      font-weight: 600;
    }}
    .skip-link:focus {{
      left: 0;
    }}
  </style>
</head>
<body>
  <a class="skip-link" href="#main-content">Skip to main content</a>
  <div class="container">
    <header>
      <div>
        <h1>PNN Discovery Studio & Preclinical Dashboard</h1>
        <div class="subtitle">Phase 4-5 Perineuronal Net Degradation Candidate Triage & Translational Profile</div>
      </div>
      <div class="badge badge-candidate">Turnkey Interactive Report</div>
    </header>

    <main id="main-content">
    <div style="margin-bottom: 24px; padding: 12px 16px; background: rgba(56, 189, 248, 0.08); border-left: 4px solid var(--primary); border-radius: 4px;">
      <strong style="color: var(--primary);">Evidence Invariant (Guidelines G-001 & G-002):</strong>
      In silico property screening, docking consensus scores, circuit ODE simulations, and allometric dose estimates represent prioritization hypotheses. Measured in vitro binding, CSPG degradation, animal PK/PD, and GLP toxicology are strictly required before clinical administration.
    </div>

    <div class="funnel-container">
      <div class="funnel-card">
        <div class="funnel-step">Stage 1</div>
        <div class="funnel-title">In Silico Triaged</div>
        <div class="funnel-num">{len(compounds)}</div>
      </div>
      <div class="funnel-card">
        <div class="funnel-step">Stage 2</div>
        <div class="funnel-title">Multi-Target Docked</div>
        <div class="funnel-num">{sum(1 for c in compounds if c.docking_scores)}</div>
      </div>
      <div class="funnel-card">
        <div class="funnel-step">Stage 3</div>
        <div class="funnel-title">In Vitro Qualified</div>
        <div class="funnel-num">{sum(1 for c in compounds if c.ec50_nm and c.ec50_nm <= 100.0)}</div>
      </div>
      <div class="funnel-card">
        <div class="funnel-step">Stage 4</div>
        <div class="funnel-title">Clinical Nominees</div>
        <div class="funnel-num">{sum(1 for c in compounds if c.brain_kp and c.brain_kp >= 0.30 and c.oral_f_pct and c.oral_f_pct >= 20.0)}</div>
      </div>
    </div>

    <h2 class="sr-only">Lead candidates</h2>
    <div class="cards-row">
      {"".join(cards_html)}
    </div>

    <div class="grid-2col">
      <div class="panel">
        <h2>CNS Physicochemical Radar (Spider Chart)</h2>
        <div class="svg-wrapper">
          {radar_svg}
          {radar_legend}
        </div>
        {radar_table}
      </div>
      <div class="panel">
        <h2>Multi-Target Docking Affinity (kcal/mol)</h2>
        <div class="svg-wrapper">
          {heatmap_svg}
        </div>
        {heatmap_table}
      </div>
    </div>

    <div class="panel" style="margin-bottom: 24px;">
      <h2>In Vivo Pharmacokinetics (PK Profile: PNN-101) <span style="font-size:0.8rem; font-weight:400; color:var(--text-muted);">Blue solid: IV Plasma | Green dashed: PO Plasma | Amber dotted: PO Brain</span></h2>
      <div class="svg-wrapper">
        {pk_svg}
      </div>
      {pk_table}
    </div>

    {circuit_section}

    {dossier_section}

    {safety_section}

    {wave_section}

    {targets_section}
    </main>

    <footer>
      PNN Project Preclinical Pipeline &bull; Standard Library Implementation &bull; Principles G-001, G-002, G-003 Compliant
    </footer>
  </div>
</body>
</html>
"""


def _parse_opt_float(val: Optional[str]) -> Optional[float]:
    """Parse optional numeric float string, returning None for empty or non-numeric cells."""
    if val is None:
        return None
    clean = val.strip()
    if not clean or clean.upper() in ("NA", "NAN", "NULL", "NONE", "-"):
        return None
    try:
        f = float(clean)
        return f if math.isfinite(f) else None
    except ValueError:  # hygiene: allow non-numeric cell in optional CSV column
        return None


def load_studio_data(
    compounds_path: Path,
    docking_path: Optional[Path] = None,
    assays_path: Optional[Path] = None,
    pk_path: Optional[Path] = None,
    preclinical_candidates_path: Optional[Path] = None,
) -> Tuple[List[StudioCompound], List[str]]:
    """Assemble all project datasets into StudioCompound instances."""
    if not compounds_path.is_file():
        raise StudioError(f"Compound CSV not found: {compounds_path}")

    # Load base compounds
    raw_compounds: Dict[str, Dict[str, str]] = {}
    with open(compounds_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise StudioError(f"Empty compounds CSV: {compounds_path}")
        norm = {c.strip().lower(): c for c in reader.fieldnames}
        cid_col = norm.get("compound_id")
        if not cid_col:
            raise StudioError(f"Compounds CSV missing compound_id: {compounds_path}")
        for row in reader:
            cid = _safe_cell(row.get(cid_col))
            if cid:
                raw_compounds[cid] = {k: _safe_cell(row.get(v)) for k, v in norm.items()}

    # Load docking scores
    dock_map: Dict[str, Dict[str, float]] = {}
    target_names: List[str] = ["PTPsigma", "LAR", "NgR1", "MMP2", "ChABC"]
    if docking_path and docking_path.is_file():
        with open(docking_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                norm = {c.strip().lower(): c for c in reader.fieldnames}
                cid_col = norm.get("compound_id")
                tgt_col = norm.get("target") or norm.get("target_id")
                score_col = norm.get("vina_affinity_kcal_mol") or norm.get("score")
                if cid_col and tgt_col and score_col:
                    for row in reader:
                        cid = _safe_cell(row.get(cid_col))
                        tgt = _safe_cell(row.get(tgt_col))
                        val = _parse_opt_float(row.get(score_col))
                        if val is not None:
                            dock_map.setdefault(cid, {})[tgt] = val

    # Load assays
    ec50_map: Dict[str, float] = {}
    if assays_path and assays_path.is_file():
        with open(assays_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                norm = {c.strip().lower(): c for c in reader.fieldnames}
                cid_col = norm.get("compound_id")
                ec_col = norm.get("ec50_nm") or norm.get("pnn_ec50_nm")
                if cid_col and ec_col:
                    for row in reader:
                        cid = _safe_cell(row.get(cid_col))
                        ec_val = _parse_opt_float(row.get(ec_col))
                        if ec_val is not None:
                            ec50_map[cid] = ec_val

    # Load PK points
    pk_map: Dict[str, Dict[str, List[Tuple[float, float]]]] = {}
    if pk_path and pk_path.is_file():
        with open(pk_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                norm = {c.strip().lower(): c for c in reader.fieldnames}
                cid_col = norm.get("compound_id")
                route_col = norm.get("route")
                tissue_col = norm.get("tissue")
                time_col = norm.get("time_h")
                conc_col = norm.get("concentration_ng_ml")
                if cid_col and route_col and tissue_col and time_col and conc_col:
                    for row in reader:
                        cid = _safe_cell(row.get(cid_col))
                        route = _safe_cell(row.get(route_col)).upper()
                        tissue = _safe_cell(row.get(tissue_col)).lower()
                        t = _parse_opt_float(row.get(time_col))
                        c = _parse_opt_float(row.get(conc_col))
                        if t is not None and c is not None:
                            pk_map.setdefault(cid, {}).setdefault(f"{route}_{tissue}", []).append((t, c))

    # Load preclinical candidates summary if available
    kp_map: Dict[str, float] = {}
    f_map: Dict[str, float] = {}
    if preclinical_candidates_path and preclinical_candidates_path.is_file():
        with open(preclinical_candidates_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                norm = {c.strip().lower(): c for c in reader.fieldnames}
                cid_col = norm.get("compound_id")
                kp_col = norm.get("brain_plasma_ratio")
                f_col = norm.get("oral_bioavailability_pct")
                if cid_col:
                    for row in reader:
                        cid = _safe_cell(row.get(cid_col))
                        if kp_col:
                            kp_val = _parse_opt_float(row.get(kp_col))
                            if kp_val is not None:
                                kp_map[cid] = kp_val
                        if f_col:
                            f_val = _parse_opt_float(row.get(f_col))
                            if f_val is not None:
                                f_map[cid] = f_val

    compounds: List[StudioCompound] = []
    for cid, row in raw_compounds.items():
        name = row.get("name") or cid
        mw = _parse_required_float(row.get("mw", ""), "mw", cid)
        logp = _parse_required_float(row.get("logp", ""), "logp", cid)
        tpsa = _parse_required_float(row.get("tpsa", ""), "tpsa", cid)
        hbd = _parse_required_int(row.get("hbd", ""), "hbd", cid)
        hba = _parse_required_int(row.get("hba", ""), "hba", cid)
        rotb = _parse_required_int(row.get("rotatable_bonds", ""), "rotatable_bonds", cid)

        ec50 = ec50_map.get(cid)
        if ec50 is None and "pnn_ec50_nm" in row:
            ec50 = _parse_opt_float(row["pnn_ec50_nm"])

        kp = kp_map.get(cid)
        if kp is None and "brain_plasma_ratio" in row:
            kp = _parse_opt_float(row["brain_plasma_ratio"])

        oral_f = f_map.get(cid)
        if oral_f is None and "oral_bioavailability" in row:
            frac = _parse_opt_float(row["oral_bioavailability"])
            if frac is not None:
                oral_f = frac * 100.0

        # Missing docking evidence stays missing: the heatmap renders N/A.
        # See ADR-010/ADR-013: never synthesize binding affinities for display.
        scores = dict(dock_map.get(cid, {}))

        c_pk = pk_map.get(cid, {})

        compounds.append(
            StudioCompound(
                compound_id=cid,
                name=name,
                mw=mw,
                logp=logp,
                tpsa=tpsa,
                hbd=hbd,
                hba=hba,
                rotb=rotb,
                ec50_nm=ec50,
                brain_kp=kp,
                oral_f_pct=oral_f,
                docking_scores=scores,
                pk_points=c_pk,
            )
        )

    # Sort compounds: best EC50 first, then highest Brain Kp
    compounds.sort(key=lambda c: (c.ec50_nm if c.ec50_nm is not None else 1e9, -(c.brain_kp or 0.0)))
    return compounds, target_names


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.studio",
        description="Generate an interactive, standalone HTML/SVG visualization studio and dashboard for PNN leads.",
    )
    parser.add_argument("--compounds", required=True, help="path to compound library CSV")
    parser.add_argument("--docking", help="optional path to multi-target docking CSV")
    parser.add_argument("--assays", help="optional path to in vitro assays CSV")
    parser.add_argument("--pk", help="optional path to in vivo PK time-course CSV")
    parser.add_argument("--candidates", help="optional path to preclinical candidates CSV")
    parser.add_argument("--dossier", help="optional path to human dose projections CSV")
    parser.add_argument("--safety", help="optional path to safety assessments CSV")
    parser.add_argument("--simulation", help="optional path to pathway simulation CSV")
    parser.add_argument("--out", default="studio.html", help="output destination HTML file (default: studio.html)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        compounds, targets = load_studio_data(
            compounds_path=Path(args.compounds),
            docking_path=Path(args.docking) if args.docking else None,
            assays_path=Path(args.assays) if args.assays else None,
            pk_path=Path(args.pk) if args.pk else None,
            preclinical_candidates_path=Path(args.candidates) if args.candidates else None,
        )

        radar_svg = generate_radar_svg(compounds)
        heatmap_svg = generate_docking_heatmap_svg(compounds[:7], targets)
        pk_svg = generate_pk_curve_svg(compounds)

        circuit_svg = None
        if args.simulation:
            sim_path = Path(args.simulation)
            if not sim_path.is_file():
                raise StudioError(f"simulation file not found: {sim_path}")
            with open(sim_path, "r", encoding="utf-8") as handle:
                circuit_svg = generate_circuit_dynamics_svg(list(csv.DictReader(handle)))

        dossier_data = None
        if args.dossier:
            dos_path = Path(args.dossier)
            if not dos_path.is_file():
                raise StudioError(f"dossier file not found: {dos_path}")
            with open(dos_path, "r", encoding="utf-8") as handle:
                dossier_data = list(csv.DictReader(handle))

        safety_data = None
        if args.safety:
            saf_path = Path(args.safety)
            if not saf_path.is_file():
                raise StudioError(f"safety assessments file not found: {saf_path}")
            with open(saf_path, "r", encoding="utf-8") as handle:
                safety_data = list(csv.DictReader(handle))

        html_out = render_studio_html(
            compounds,
            radar_svg=radar_svg,
            heatmap_svg=heatmap_svg,
            pk_svg=pk_svg,
            circuit_svg=circuit_svg,
            dossier_data=dossier_data,
            safety_data=safety_data,
            targets_data=list_targets(),
        )
        out_path = Path(args.out)
        write_text(out_path, html_out)

    except StudioError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: Studio could not load one of the requested input datasets.", file=sys.stderr)
        print("Fix: verify each supplied path and its CSV columns, then rerun.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: Studio dashboard HTML could not be written to disk.", file=sys.stderr)
        print("Fix: check write permissions for the destination path.", file=sys.stderr)
        return 1

    print(f"Studio dashboard rendered successfully: {len(compounds)} compounds -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
