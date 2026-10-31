"""Hypothalamic circuit, perineuronal net turnover, and dosing dynamics simulator for PNN."""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .reports import write_csv, write_json, write_text


class PathwayError(ValueError):
    """Raised when simulation parameters are physically invalid or malformed."""


@dataclass(frozen=True)
class SimulationPoint:
    time_h: float
    brain_conc_ng_ml: float
    pnn_density_pct: float
    ir_accessibility_pct: float
    akt_signaling_pct: float
    food_intake_reduction_pct: float
    cumulative_weight_delta_pct: float
    plasticity_risk_tier: str

    def as_row(self) -> Dict[str, object]:
        return {
            "time_h": f"{self.time_h:.1f}",
            "brain_conc_ng_ml": f"{self.brain_conc_ng_ml:.1f}",
            "pnn_density_pct": f"{self.pnn_density_pct:.1f}",
            "ir_accessibility_pct": f"{self.ir_accessibility_pct:.1f}",
            "akt_signaling_pct": f"{self.akt_signaling_pct:.1f}",
            "food_intake_reduction_pct": f"{self.food_intake_reduction_pct:.1f}",
            "cumulative_weight_delta_pct": f"{self.cumulative_weight_delta_pct:.2f}",
            "plasticity_risk_tier": self.plasticity_risk_tier,
        }


@dataclass(frozen=True)
class RegimenSummary:
    regimen_name: str
    compound_id: str
    dose_mg_kg: float
    duration_days: int
    mean_pnn_reduction_pct: float
    peak_pnn_reduction_pct: float
    mean_akt_restoration_pct: float
    projected_weight_loss_pct: float
    plasticity_safety_score: float  # 0 to 100
    recommendation: str

    def as_row(self) -> Dict[str, object]:
        return {
            "regimen_name": self.regimen_name,
            "compound_id": self.compound_id,
            "dose_mg_kg": f"{self.dose_mg_kg:.1f}",
            "duration_days": self.duration_days,
            "mean_pnn_reduction_pct": f"{self.mean_pnn_reduction_pct:.1f}",
            "peak_pnn_reduction_pct": f"{self.peak_pnn_reduction_pct:.1f}",
            "mean_akt_restoration_pct": f"{self.mean_akt_restoration_pct:.1f}",
            "projected_weight_loss_pct": f"{self.projected_weight_loss_pct:.1f}",
            "plasticity_safety_score": f"{self.plasticity_safety_score:.1f}",
            "recommendation": self.recommendation,
        }


SIMULATION_FIELDS = (
    "time_h",
    "brain_conc_ng_ml",
    "pnn_density_pct",
    "ir_accessibility_pct",
    "akt_signaling_pct",
    "food_intake_reduction_pct",
    "cumulative_weight_delta_pct",
    "plasticity_risk_tier",
)

REGIMEN_FIELDS = (
    "regimen_name",
    "compound_id",
    "dose_mg_kg",
    "duration_days",
    "mean_pnn_reduction_pct",
    "peak_pnn_reduction_pct",
    "mean_akt_restoration_pct",
    "projected_weight_loss_pct",
    "plasticity_safety_score",
    "recommendation",
)


def _get_dosing_times(regimen: str, duration_h: float) -> List[float]:
    """Return schedule of dose administration times in hours."""
    clean = regimen.strip().upper()
    times: List[float] = []
    t = 0.0
    if clean == "CONTINUOUS_DAILY":
        while t < duration_h:
            times.append(t)
            t += 24.0
    elif clean in ("INTERMITTENT_PULSATILE", "ALTERNATE_DAY"):
        while t < duration_h:
            times.append(t)
            t += 48.0
    elif clean == "TWICE_WEEKLY":
        week = 0.0
        while week < duration_h:
            times.append(week)
            if week + 72.0 < duration_h:
                times.append(week + 72.0)
            week += 168.0
    elif clean == "ONCE_WEEKLY":
        while t < duration_h:
            times.append(t)
            t += 168.0
    elif clean == "SINGLE_DOSE":
        times.append(0.0)
    else:
        raise PathwayError(f"Unknown dosing regimen '{regimen}'. Supported: CONTINUOUS_DAILY, INTERMITTENT_PULSATILE, TWICE_WEEKLY, ONCE_WEEKLY, SINGLE_DOSE")
    return times


def simulate_circuit_dynamics(
    compound_id: str = "PNN-101",
    dose_mg_kg: float = 10.0,
    regimen: str = "INTERMITTENT_PULSATILE",
    duration_days: int = 14,
    dt_h: float = 0.25,
    half_life_h: float = 3.2,
    cmax_per_dose: float = 32.5,  # (ng/mL brain) per (mg/kg dose)
    ec50_nm: float = 24.0,
    mw: float = 349.46,
) -> Tuple[List[SimulationPoint], RegimenSummary]:
    """
    Run an explicit forward-Euler integration of hypothalamic circuit dynamics.
    Models drug absorption/elimination, PNN turnover, insulin receptor accessibility,
    downstream AKT phosphorylation, and body weight progression.
    """
    if dose_mg_kg <= 0:
        raise PathwayError(f"Dose must be > 0 mg/kg, got {dose_mg_kg}")
    if duration_days <= 0 or duration_days > 90:
        raise PathwayError(f"Duration must be between 1 and 90 days, got {duration_days}")
    if dt_h <= 0 or dt_h > 2.0:
        raise PathwayError(f"Step dt must be between 0.01 and 2.0 hours, got {dt_h}")
    if half_life_h <= 0:
        raise PathwayError(f"Half-life must be > 0 h, got {half_life_h}")
    if cmax_per_dose <= 0 or ec50_nm <= 0 or mw <= 0:
        raise PathwayError("cmax_per_dose, ec50_nm and mw must all be > 0")

    duration_h = duration_days * 24.0
    dose_times = set(_get_dosing_times(regimen, duration_h))

    # Biological constants
    # ke elimination rate constant
    ke = math.log(2.0) / half_life_h
    ka = 1.2  # absorption rate constant (h^-1)
    # Convert ec50 from nM to ng/mL: ng/mL = nM * MW / 1000
    ec50_ng_ml = max(1.0, ec50_nm * mw / 1000.0)

    # PNN synthesis and degradation rates
    k_synth = 0.012  # endogenous turnover: recovery half-life ~58h
    k_deg_max = 0.085  # maximum enzymatic / competitive degradation rate (h^-1)

    # Signaling & weight constants. AKT is homeostatic around the DIO baseline, so an untreated
    # animal reports zero restoration rather than drifting to the model's own equilibrium.
    k_akt = 0.15
    k_akt_decay = 0.12
    akt_baseline = 0.25
    ir_baseline = max(0.1, 1.0 - (0.75 * 1.0))

    # State variables
    # c_depot (unabsorbed oral dose in GI), c_brain, pnn (1.0 = 100% baseline), akt (relative to max)
    c_depot = 0.0
    c_brain = 0.0
    pnn = 1.0  # normalized density: 1.0 = baseline hypertrophic PNN
    akt = akt_baseline
    weight_delta = 0.0

    points: List[SimulationPoint] = []
    severe_plasticity_hours = 0.0
    pnn_values: List[float] = []
    akt_values: List[float] = []

    steps = int(duration_h / dt_h)
    for step in range(steps + 1):
        t = step * dt_h

        # Administer dose if at scheduled time
        # Find if any dose time falls within [t - dt_h/2, t + dt_h/2)
        for dt_time in dose_times:
            if abs(t - dt_time) < (dt_h / 2.0):
                c_depot += dose_mg_kg * cmax_per_dose * 1.5

        # Pharmacokinetics ODEs:
        # d(depot)/dt = -ka * depot
        # d(brain)/dt = ka * depot - ke * brain
        d_depot = -ka * c_depot
        d_brain = ka * c_depot - ke * c_brain

        # Pharmacodynamics ODEs:
        # PNN net reduction:
        # drug effect = c_brain / (c_brain + ec50)
        drug_effect = c_brain / (c_brain + ec50_ng_ml)
        d_pnn = k_synth * (1.0 - pnn) - (k_deg_max * drug_effect * pnn)

        # IR accessibility is inverse of PNN density: A = 1.0 - 0.7 * pnn
        ir_access = max(0.1, 1.0 - (0.75 * pnn))

        # AKT activation is driven by the gain in receptor accessibility over the DIO baseline.
        d_akt = (k_akt * (ir_access - ir_baseline)) - (k_akt_decay * (akt - akt_baseline))

        # Integration step (Euler-Heun predictor)
        c_depot = max(0.0, c_depot + d_depot * dt_h)
        c_brain = max(0.0, c_brain + d_brain * dt_h)
        pnn = max(0.05, min(1.0, pnn + d_pnn * dt_h))
        akt = max(0.1, min(1.0, akt + d_akt * dt_h))

        # Metrics
        pnn_pct = pnn * 100.0
        ir_access_pct = ir_access * 100.0
        akt_pct = akt * 100.0
        # Food intake reduction proportional to AKT signaling restoration (0 to 35%)
        food_red_pct = max(0.0, min(40.0, (akt - akt_baseline) * 50.0))
        # Daily weight delta integration (1% weight loss per ~3 days of 25% food intake reduction)
        d_weight_day = -(food_red_pct / 25.0) * (1.0 / 3.0)  # % weight delta per day
        weight_delta += d_weight_day * (dt_h / 24.0)

        # Plasticity Risk Assessment
        if pnn_pct < 20.0:
            severe_plasticity_hours += dt_h
            risk = "HIGH_PLASTICITY_RISK"
        elif pnn_pct < 40.0:
            risk = "MODERATE_PLASTICITY_WINDOW"
        else:
            risk = "SAFE_PHYSIOLOGICAL_RANGE"

        pnn_values.append(pnn_pct)
        akt_values.append(akt_pct)

        # Record simulation point (sample every 1 hour or at final step)
        if step % max(1, int(1.0 / dt_h)) == 0 or step == steps:
            points.append(
                SimulationPoint(
                    time_h=t,
                    brain_conc_ng_ml=c_brain,
                    pnn_density_pct=pnn_pct,
                    ir_accessibility_pct=ir_access_pct,
                    akt_signaling_pct=akt_pct,
                    food_intake_reduction_pct=food_red_pct,
                    cumulative_weight_delta_pct=weight_delta,
                    plasticity_risk_tier=risk,
                )
            )

    mean_pnn = sum(pnn_values) / len(pnn_values)
    mean_pnn_red = max(0.0, 100.0 - mean_pnn)
    peak_pnn_red = max(0.0, 100.0 - min(pnn_values))
    mean_akt_level = sum(akt_values) / len(akt_values)
    # Restoration is relative to the DIO baseline, not the absolute pAKT level.
    mean_akt_restoration = max(0.0, min(100.0, ((mean_akt_level / 100.0) - akt_baseline) / (1.0 - akt_baseline) * 100.0))

    # Plasticity Safety Score (100 = optimal therapeutic window, penalizes prolonged <20% PNN)
    plasticity_penalty = min(80.0, (severe_plasticity_hours / 24.0) * 15.0)
    safety_score = max(10.0, 100.0 - plasticity_penalty)

    # Recommendation verdict is derived entirely from the computed metrics.
    if mean_pnn_red < 15.0:
        verdict = "SUB-THERAPEUTIC: Insufficient PNN digestion; consider increasing dose or frequency."
    elif safety_score < 80.0 or mean_pnn_red < 20.0:
        verdict = "CAUTION: Prolonged deep PNN depletion; reduce exposure or extend the dosing interval."
    elif mean_akt_restoration < 40.0:
        verdict = "CAUTION: PNN reduction outweighs insulin signaling recovery; verify the downstream coupling."
    else:
        verdict = "RECOMMENDED: Balanced metabolic reversal with physiological plasticity preservation."

    summary = RegimenSummary(
        regimen_name=regimen,
        compound_id=compound_id,
        dose_mg_kg=dose_mg_kg,
        duration_days=duration_days,
        mean_pnn_reduction_pct=mean_pnn_red,
        peak_pnn_reduction_pct=peak_pnn_red,
        mean_akt_restoration_pct=mean_akt_restoration,
        projected_weight_loss_pct=abs(weight_delta),
        plasticity_safety_score=safety_score,
        recommendation=verdict,
    )

    return points, summary


def render_pathway_report(summary: RegimenSummary, sample_points: Sequence[SimulationPoint]) -> str:
    """Render comprehensive markdown simulation report for dosing regimen."""
    lines = [
        f"# PNN Hypothalamic Circuit & Dosing Dynamics Report: {summary.compound_id}",
        "",
        "## Executive Summary",
        "",
        f"- **Compound ID**: `{summary.compound_id}`",
        f"- **Simulated Regimen**: **{summary.regimen_name}** ({summary.dose_mg_kg} mg/kg PO for {summary.duration_days} days)",
        f"- **Mean Hypothalamic PNN Reduction**: **{summary.mean_pnn_reduction_pct:.1f}%** (Peak: {summary.peak_pnn_reduction_pct:.1f}%)",
        f"- **Insulin Signaling (pAKT) Restoration**: **{summary.mean_akt_restoration_pct:.1f}%** of the recoverable gap above the DIO baseline",
        f"- **Projected Cumulative Weight Loss**: **{summary.projected_weight_loss_pct:.1f}%**",
        f"- **Plasticity Safety Score**: **{summary.plasticity_safety_score:.1f} / 100**",
        f"- **Clinical Recommendation**: *{summary.recommendation}*",
        "",
        "## Scientific Rationale & Plasticity Risk Avoidance",
        "",
        "> [!IMPORTANT]",
        "> **Plasticity Risk Mitigation (ADR-008)**: Hypothalamic perineuronal nets maintain structural synaptic stability. "
        "Continuous unremitting degradation (<20% baseline density) risks anomalous synaptic remodeling. Intermittent "
        "pulsatile dosing allows periodic net stabilization while sustaining POMC insulin receptor accessibility.",
        "> This report models a single regimen. The recommendation is derived from that regimen's computed metrics, "
        "not from a cross-regimen comparison.",
        "",
        "## Simulated Circuit Trajectory (Sampled Windows)",
        "",
        "| Time (h) | Brain Conc (ng/mL) | PNN Density (%) | IR Access (%) | pAKT Signal (%) | Weight Delta (%) | Plasticity Risk |",
        "|---|---|---|---|---|---|---|",
    ]

    # Sample points: every 24h
    for pt in sample_points:
        if abs(pt.time_h % 24.0) < 0.3 or pt.time_h == sample_points[-1].time_h:
            lines.append(
                f"| {pt.time_h:.0f}h | {pt.brain_conc_ng_ml:.1f} | {pt.pnn_density_pct:.1f}% | "
                f"{pt.ir_accessibility_pct:.1f}% | {pt.akt_signaling_pct:.1f}% | {pt.cumulative_weight_delta_pct:.2f}% | "
                f"*{pt.plasticity_risk_tier}* |"
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- **PNN effect**: mean reduction of {summary.mean_pnn_reduction_pct:.1f}% with a peak of "
            f"{summary.peak_pnn_reduction_pct:.1f}%.",
            f"- **Signaling effect**: {summary.mean_akt_restoration_pct:.1f}% of the recoverable pAKT gap was restored "
            "above the untreated DIO baseline.",
            f"- **Plasticity margin**: {summary.plasticity_safety_score:.1f}/100, penalised for hours spent below 20% "
            "PNN density.",
            "- **Model limits**: this is a deterministic ODE abstraction with literature-derived constants. It ranks "
            "dosing hypotheses; it does not predict efficacy, safety or human exposure, and it must not be reported "
            "as measured pharmacology (G-001).",
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.pathway",
        description="Simulate hypothalamic circuit ODE dynamics, PNN turnover, insulin signaling, and dosing regimens.",
    )
    parser.add_argument("--compound", default="PNN-101", help="compound identifier (default: PNN-101)")
    parser.add_argument("--dose", type=float, default=10.0, help="dose in mg/kg PO (default: 10.0)")
    parser.add_argument(
        "--regimen",
        choices=("CONTINUOUS_DAILY", "INTERMITTENT_PULSATILE", "TWICE_WEEKLY", "ONCE_WEEKLY", "SINGLE_DOSE"),
        default="INTERMITTENT_PULSATILE",
        help="dosing schedule (default: INTERMITTENT_PULSATILE)",
    )
    parser.add_argument("--days", type=int, default=14, help="simulation duration in days (default: 14)")
    parser.add_argument("--dt", type=float, default=0.25, help="ODE integration step in hours (default: 0.25)")
    parser.add_argument("--out", default="dosing_simulation.csv", help="output destination CSV for time series")
    parser.add_argument("--report", dest="report_path", help="optional destination path for Markdown report")
    parser.add_argument("--json", dest="json_path", help="optional destination path for JSON summary")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        points, summary = simulate_circuit_dynamics(
            compound_id=args.compound,
            dose_mg_kg=args.dose,
            regimen=args.regimen,
            duration_days=args.days,
            dt_h=args.dt,
        )

        out_path = Path(args.out)
        write_csv(out_path, SIMULATION_FIELDS, [pt.as_row() for pt in points])

        if args.report_path:
            report_text = render_pathway_report(summary, points)
            write_text(Path(args.report_path), report_text)

        if args.json_path:
            json_data = {
                "summary": summary.as_row(),
                "time_points_count": len(points),
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), json_data)

    except PathwayError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: simulation requires physically plausible doses (>0 mg/kg), valid regimens, and durations <=90 days.", file=sys.stderr)
        print("Fix: verify CLI arguments (--dose, --regimen, --days).", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: simulation output files could not be written to disk.", file=sys.stderr)
        print("Fix: specify a writable destination for --out, --report, or --json.", file=sys.stderr)
        return 1

    print(
        f"Circuit simulation complete: {summary.compound_id} ({summary.regimen_name}, {summary.dose_mg_kg} mg/kg) -> "
        f"Mean PNN reduction: {summary.mean_pnn_reduction_pct:.1f}%, Weight loss: {summary.projected_weight_loss_pct:.1f}%, "
        f"Safety score: {summary.plasticity_safety_score:.1f}/100 -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
