"""Phase 4 in vivo pharmacokinetics (PK), pharmacodynamics (PD), and preclinical candidate selection for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .reports import write_csv, write_json, write_text


class PreclinicalError(ValueError):
    """Raised when in vivo preclinical data is malformed or violates contract."""


@dataclass(frozen=True)
class ConcentrationPoint:
    time_h: float
    concentration: float


@dataclass(frozen=True)
class PkProfile:
    compound_id: str
    route: str  # PO or IV
    tissue: str  # plasma or brain
    dose_mg_kg: float
    cmax: float
    tmax_h: float
    auc_0_t: float
    auc_0_inf: Optional[float]
    half_life_h: Optional[float]
    lambda_z: Optional[float]
    cl_l_h_kg: Optional[float]
    vd_l_kg: Optional[float]


@dataclass(frozen=True)
class PdOutcome:
    compound_id: str
    dose_mg_kg: float
    wfa_pnn_reduction_pct: float
    itt_glucose_auc_reduction_pct: float
    food_intake_reduction_pct: Optional[float] = None
    body_weight_delta_pct: Optional[float] = None


@dataclass(frozen=True)
class PreclinicalCandidate:
    compound_id: str
    status: str  # CLINICAL_CANDIDATE, ADVANCED_LEAD, INSUFFICIENT_BRAIN_PENETRATION, POOR_BIOAVAILABILITY, INSUFFICIENT_EFFICACY
    oral_bioavailability_pct: Optional[float]  # F (%)
    brain_plasma_ratio: Optional[float]  # Kp,brain
    plasma_cmax_ng_ml: Optional[float]
    plasma_auc_h_ng_ml: Optional[float]
    brain_auc_h_ng_ml: Optional[float]
    half_life_h: Optional[float]
    wfa_pnn_reduction_pct: Optional[float]
    itt_glucose_auc_reduction_pct: Optional[float]
    findings: str

    def as_row(self) -> Dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "status": self.status,
            "oral_bioavailability_pct": f"{self.oral_bioavailability_pct:.1f}" if self.oral_bioavailability_pct is not None else "",
            "brain_plasma_ratio": f"{self.brain_plasma_ratio:.2f}" if self.brain_plasma_ratio is not None else "",
            "plasma_cmax_ng_ml": f"{self.plasma_cmax_ng_ml:.1f}" if self.plasma_cmax_ng_ml is not None else "",
            "plasma_auc_h_ng_ml": f"{self.plasma_auc_h_ng_ml:.1f}" if self.plasma_auc_h_ng_ml is not None else "",
            "brain_auc_h_ng_ml": f"{self.brain_auc_h_ng_ml:.1f}" if self.brain_auc_h_ng_ml is not None else "",
            "half_life_h": f"{self.half_life_h:.1f}" if self.half_life_h is not None else "",
            "wfa_pnn_reduction_pct": f"{self.wfa_pnn_reduction_pct:.1f}" if self.wfa_pnn_reduction_pct is not None else "",
            "itt_glucose_auc_reduction_pct": f"{self.itt_glucose_auc_reduction_pct:.1f}" if self.itt_glucose_auc_reduction_pct is not None else "",
            "findings": self.findings,
        }


PRECLINICAL_FIELDS = (
    "compound_id",
    "status",
    "oral_bioavailability_pct",
    "brain_plasma_ratio",
    "plasma_cmax_ng_ml",
    "plasma_auc_h_ng_ml",
    "brain_auc_h_ng_ml",
    "half_life_h",
    "wfa_pnn_reduction_pct",
    "itt_glucose_auc_reduction_pct",
    "findings",
)


def _parse_num(val: str, field_name: str, row_id: str) -> Optional[float]:
    clean = val.strip()
    if not clean or clean.upper() in ("NA", "NAN", "NULL", "NONE"):
        return None
    try:
        f = float(clean)
    except ValueError:
        raise PreclinicalError(f"Row '{row_id}': invalid numeric value '{val}' for {field_name}")
    if not math.isfinite(f):
        raise PreclinicalError(f"Row '{row_id}': non-finite number '{val}' for {field_name}")
    return f


def calculate_nca(
    points: Sequence[ConcentrationPoint], dose_mg_kg: float, compound_id: str, route: str, tissue: str
) -> PkProfile:
    """Perform non-compartmental pharmacokinetic analysis (NCA) on concentration-time points."""
    sorted_pts = sorted(points, key=lambda p: p.time_h)
    if not sorted_pts:
        raise PreclinicalError(f"Compound '{compound_id}' ({route}/{tissue}): no time points provided")

    cmax_pt = max(sorted_pts, key=lambda p: p.concentration)
    cmax = cmax_pt.concentration
    tmax = cmax_pt.time_h

    # Linear trapezoidal rule for AUC_0_t
    auc_0_t = 0.0
    for i in range(len(sorted_pts) - 1):
        dt = sorted_pts[i + 1].time_h - sorted_pts[i].time_h
        if dt < 0:
            continue
        c_mean = (sorted_pts[i].concentration + sorted_pts[i + 1].concentration) / 2.0
        auc_0_t += c_mean * dt

    # Terminal rate constant estimation (linear regression on log(C) of terminal phase)
    # Require at least 3 positive concentration points after Tmax
    terminal_candidates = [p for p in sorted_pts if p.time_h >= tmax and p.concentration > 0]
    lambda_z = None
    half_life = None
    auc_0_inf = None
    cl = None
    vd = None

    if len(terminal_candidates) >= 3:
        n = len(terminal_candidates)
        sum_t = sum(p.time_h for p in terminal_candidates)
        sum_ln_c = sum(math.log(p.concentration) for p in terminal_candidates)
        sum_t2 = sum(p.time_h**2 for p in terminal_candidates)
        sum_t_ln_c = sum(p.time_h * math.log(p.concentration) for p in terminal_candidates)

        denom = (n * sum_t2) - (sum_t**2)
        if denom != 0:
            slope = ((n * sum_t_ln_c) - (sum_t * sum_ln_c)) / denom
            if slope < 0:
                lambda_z = -slope
                half_life = math.log(2.0) / lambda_z

    elif len(terminal_candidates) == 2:
        p1, p2 = terminal_candidates[0], terminal_candidates[1]
        dt = p2.time_h - p1.time_h
        if dt > 0 and p1.concentration > p2.concentration:
            lambda_z = (math.log(p1.concentration) - math.log(p2.concentration)) / dt
            half_life = math.log(2.0) / lambda_z

    last_pt = sorted_pts[-1]
    if lambda_z is not None and lambda_z > 0 and last_pt.concentration > 0:
        auc_extrap = last_pt.concentration / lambda_z
        auc_0_inf = auc_0_t + auc_extrap
        # Clearance = Dose / AUC (converting dose mg/kg and conc ng/mL -> L/h/kg)
        # dose (mg/kg) * 1e6 (ng/mg) / (auc (ng*h/mL) * 1e3 (mL/L)) = dose * 1000 / auc
        if auc_0_inf > 0:
            cl = (dose_mg_kg * 1000.0) / auc_0_inf
            vd = cl / lambda_z

    return PkProfile(
        compound_id=compound_id,
        route=route.upper(),
        tissue=tissue.lower(),
        dose_mg_kg=dose_mg_kg,
        cmax=cmax,
        tmax_h=tmax,
        auc_0_t=auc_0_t,
        auc_0_inf=auc_0_inf,
        half_life_h=half_life,
        lambda_z=lambda_z,
        cl_l_h_kg=cl,
        vd_l_kg=vd,
    )


def load_pk_data(path: Path) -> Dict[str, Dict[str, PkProfile]]:
    """
    Ingest concentration-time series CSV and calculate NCA PK profiles.
    Returns nested dictionary: results[compound_id][route_tissue] -> PkProfile.
    Keys e.g.: 'PO_plasma', 'PO_brain', 'IV_plasma'.
    """
    if not path.is_file():
        raise PreclinicalError(f"PK data file not found: {path}")

    series_map: Dict[Tuple[str, str, str, float], List[ConcentrationPoint]] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise PreclinicalError(f"Empty PK data CSV: {path}")

        norm = {col.strip().lower(): col for col in reader.fieldnames}
        req = {"compound_id", "route", "tissue", "dose_mg_kg", "time_h", "concentration_ng_ml"}
        missing = req - set(norm.keys())
        if missing:
            raise PreclinicalError(f"PK CSV missing required column(s): {', '.join(missing)}")

        for row in reader:
            cid = row[norm["compound_id"]].strip()
            route = row[norm["route"]].strip().upper()
            tissue = row[norm["tissue"]].strip().lower()
            if not cid or not route or not tissue:
                continue

            dose = _parse_num(row[norm["dose_mg_kg"]], "dose_mg_kg", cid)
            if dose is None or dose <= 0:
                raise PreclinicalError(f"Row '{cid}': dose_mg_kg must be > 0")

            t = _parse_num(row[norm["time_h"]], "time_h", cid)
            if t is None or t < 0:
                raise PreclinicalError(f"Row '{cid}': time_h must be >= 0")

            conc = _parse_num(row[norm["concentration_ng_ml"]], "concentration_ng_ml", cid)
            if conc is None or conc < 0:
                raise PreclinicalError(f"Row '{cid}': concentration_ng_ml must be >= 0")

            key = (cid, route, tissue, dose)
            series_map.setdefault(key, []).append(ConcentrationPoint(time_h=t, concentration=conc))

    if not series_map:
        raise PreclinicalError(f"No valid PK concentration rows found in {path}")

    profiles: Dict[str, Dict[str, PkProfile]] = {}
    for (cid, route, tissue, dose), pts in series_map.items():
        profile = calculate_nca(pts, dose_mg_kg=dose, compound_id=cid, route=route, tissue=tissue)
        sub_key = f"{route}_{tissue}"
        profiles.setdefault(cid, {})[sub_key] = profile

    return profiles


def load_pd_data(path: Path) -> Dict[str, PdOutcome]:
    """Ingest in vivo pharmacodynamic (PD) and metabolic outcomes CSV."""
    if not path.is_file():
        raise PreclinicalError(f"PD outcomes file not found: {path}")

    outcomes: Dict[str, PdOutcome] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise PreclinicalError(f"Empty PD outcomes CSV: {path}")

        norm = {col.strip().lower(): col for col in reader.fieldnames}
        req = {"compound_id", "wfa_pnn_reduction_pct", "itt_glucose_auc_reduction_pct"}
        missing = req - set(norm.keys())
        if missing:
            raise PreclinicalError(f"PD CSV missing required column(s): {', '.join(missing)}")

        for row in reader:
            cid = row[norm["compound_id"]].strip()
            if not cid:
                continue

            dose = 10.0
            if "dose_mg_kg" in norm:
                d = _parse_num(row[norm["dose_mg_kg"]], "dose_mg_kg", cid)
                if d is not None and d > 0:
                    dose = d

            wfa = _parse_num(row[norm["wfa_pnn_reduction_pct"]], "wfa_pnn_reduction_pct", cid)
            if wfa is None or wfa < 0:
                raise PreclinicalError(f"Row '{cid}': wfa_pnn_reduction_pct must be >= 0")

            itt = _parse_num(row[norm["itt_glucose_auc_reduction_pct"]], "itt_glucose_auc_reduction_pct", cid)
            if itt is None:
                raise PreclinicalError(f"Row '{cid}': itt_glucose_auc_reduction_pct required")

            food = None
            if "food_intake_reduction_pct" in norm:
                food = _parse_num(row[norm["food_intake_reduction_pct"]], "food_intake_reduction_pct", cid)

            wt = None
            if "body_weight_delta_pct" in norm:
                wt = _parse_num(row[norm["body_weight_delta_pct"]], "body_weight_delta_pct", cid)

            outcomes[cid] = PdOutcome(
                compound_id=cid,
                dose_mg_kg=dose,
                wfa_pnn_reduction_pct=wfa,
                itt_glucose_auc_reduction_pct=itt,
                food_intake_reduction_pct=food,
                body_weight_delta_pct=wt,
            )

    return outcomes


def evaluate_candidates(
    pk_profiles: Dict[str, Dict[str, PkProfile]],
    pd_outcomes: Optional[Dict[str, PdOutcome]] = None,
    min_kp_brain: float = 0.30,
    min_oral_f_pct: float = 20.0,
    min_wfa_reduction_pct: float = 15.0,
    min_itt_improvement_pct: float = 20.0,
) -> List[PreclinicalCandidate]:
    """Triage compounds against preclinical PK exposure and PD metabolic criteria."""
    pd_map = pd_outcomes or {}
    candidates: List[PreclinicalCandidate] = []

    all_compound_ids = sorted(set(list(pk_profiles.keys()) + list(pd_map.keys())))

    for cid in all_compound_ids:
        pk_sub = pk_profiles.get(cid, {})
        pd = pd_map.get(cid)

        po_plasma = pk_sub.get("PO_plasma")
        iv_plasma = pk_sub.get("IV_plasma")
        po_brain = pk_sub.get("PO_brain")

        # Oral Bioavailability F (%)
        oral_f = None
        if po_plasma and iv_plasma and iv_plasma.auc_0_t > 0 and iv_plasma.dose_mg_kg > 0:
            oral_f = (po_plasma.auc_0_t / po_plasma.dose_mg_kg) / (iv_plasma.auc_0_t / iv_plasma.dose_mg_kg) * 100.0

        # Brain / Plasma ratio Kp,brain
        kp_brain = None
        if po_brain and po_plasma and po_plasma.auc_0_t > 0:
            kp_brain = po_brain.auc_0_t / po_plasma.auc_0_t
        elif po_brain and iv_plasma and iv_plasma.auc_0_t > 0:
            kp_brain = po_brain.auc_0_t / iv_plasma.auc_0_t

        plasma_cmax = po_plasma.cmax if po_plasma else (iv_plasma.cmax if iv_plasma else None)
        plasma_auc = po_plasma.auc_0_t if po_plasma else (iv_plasma.auc_0_t if iv_plasma else None)
        brain_auc = po_brain.auc_0_t if po_brain else None
        half_life = po_plasma.half_life_h if (po_plasma and po_plasma.half_life_h) else (iv_plasma.half_life_h if iv_plasma else None)

        wfa_red = pd.wfa_pnn_reduction_pct if pd else None
        itt_imp = pd.itt_glucose_auc_reduction_pct if pd else None

        # Qualification decision
        findings: List[str] = []
        is_candidate = True

        if kp_brain is not None:
            if kp_brain < 0.20:
                is_candidate = False
                findings.append(f"Restricted brain exposure (Kp={kp_brain:.2f} < 0.20)")
            elif kp_brain >= min_kp_brain:
                findings.append(f"High brain exposure (Kp={kp_brain:.2f})")
            else:
                findings.append(f"Moderate brain exposure (Kp={kp_brain:.2f})")
        else:
            findings.append("Brain PK unmeasured")

        if oral_f is not None:
            if oral_f < 15.0:
                is_candidate = False
                findings.append(f"Low oral bioavailability (F={oral_f:.1f}% < 15%)")
            elif oral_f >= min_oral_f_pct:
                findings.append(f"Robust oral F ({oral_f:.1f}%)")
            else:
                findings.append(f"Marginal oral F ({oral_f:.1f}%)")

        if wfa_red is not None and itt_imp is not None:
            if wfa_red < 10.0 or itt_imp < 10.0:
                is_candidate = False
                findings.append(f"Sub-therapeutic PD efficacy (WFA={wfa_red:.1f}%, ITT={itt_imp:.1f}%)")
            elif wfa_red >= min_wfa_reduction_pct and itt_imp >= min_itt_improvement_pct:
                findings.append(f"Potent metabolic reversal (WFA={wfa_red:.1f}%, ITT={itt_imp:.1f}%)")
            else:
                findings.append(f"Moderate PD response (WFA={wfa_red:.1f}%, ITT={itt_imp:.1f}%)")
        else:
            findings.append("In vivo PD unmeasured")

        # Assign final category
        if kp_brain is not None and kp_brain < 0.20:
            status = "INSUFFICIENT_BRAIN_PENETRATION"
        elif oral_f is not None and oral_f < 15.0:
            status = "POOR_BIOAVAILABILITY"
        elif wfa_red is not None and (wfa_red < 10.0 or (itt_imp is not None and itt_imp < 10.0)):
            status = "INSUFFICIENT_EFFICACY"
        elif (
            kp_brain is not None
            and kp_brain >= min_kp_brain
            and oral_f is not None and oral_f >= min_oral_f_pct
            and wfa_red is not None and wfa_red >= min_wfa_reduction_pct
            and itt_imp is not None and itt_imp >= min_itt_improvement_pct
        ):
            status = "CLINICAL_CANDIDATE"
        else:
            status = "ADVANCED_LEAD"

        candidates.append(
            PreclinicalCandidate(
                compound_id=cid,
                status=status,
                oral_bioavailability_pct=oral_f,
                brain_plasma_ratio=kp_brain,
                plasma_cmax_ng_ml=plasma_cmax,
                plasma_auc_h_ng_ml=plasma_auc,
                brain_auc_h_ng_ml=brain_auc,
                half_life_h=half_life,
                wfa_pnn_reduction_pct=wfa_red,
                itt_glucose_auc_reduction_pct=itt_imp,
                findings="; ".join(findings),
            )
        )

    # Sort candidates: CLINICAL_CANDIDATE first, then ADVANCED_LEAD, then by Kp_brain descending
    priority = {
        "CLINICAL_CANDIDATE": 0,
        "ADVANCED_LEAD": 1,
        "INSUFFICIENT_BRAIN_PENETRATION": 2,
        "POOR_BIOAVAILABILITY": 3,
        "INSUFFICIENT_EFFICACY": 4,
    }
    candidates.sort(key=lambda c: (priority.get(c.status, 9), -(c.brain_plasma_ratio or 0.0), -(c.wfa_pnn_reduction_pct or 0.0)))
    return candidates


def render_preclinical_report(candidates: Sequence[PreclinicalCandidate], criteria: Mapping[str, object]) -> str:
    """Render a comprehensive Phase 4 in vivo preclinical qualification report."""
    clinical = [c for c in candidates if c.status == "CLINICAL_CANDIDATE"]
    advanced = [c for c in candidates if c.status == "ADVANCED_LEAD"]
    failed = [c for c in candidates if c.status not in ("CLINICAL_CANDIDATE", "ADVANCED_LEAD")]

    lines = [
        "# PNN Phase 4 In Vivo Preclinical Candidate Selection Report",
        "",
        "## Scientific Context & Translation Invariants",
        "",
        "> [!IMPORTANT]",
        "> **Preclinical Guidance (G-001)**: Preclinical efficacy demonstrated in mouse DIO models verifies "
        "target engagement, oral bioavailability, and metabolic phenotype rescue; human clinical translation "
        "requires formal GLP safety pharmacology (ICH S7A/S7B) and Phase 1 human PK dose-escalation trials.",
        "",
        "## Candidate Advancement Summary",
        "",
        f"- **Evaluated Candidates**: {len(candidates)}",
        f"- **Clinical Candidates (Kp >= {criteria.get('min_kp_brain', 0.30)}, F >= {criteria.get('min_oral_f_pct', 20.0)}%, WFA >= {criteria.get('min_wfa_reduction_pct', 15.0)}%)**: {len(clinical)}",
        f"- **Advanced Leads (Sub-optimal exposure or pending complete PD)**: {len(advanced)}",
        f"- **Failed Preclinical Thresholds (Low Kp / Poor F / Inactive)**: {len(failed)}",
        "",
        "## Preclinical Decision Thresholds",
        "",
        f"- **Brain Partitioning**: Kp,brain (AUC_brain / AUC_plasma) >= {criteria.get('min_kp_brain', 0.30)}",
        f"- **Oral Bioavailability**: F >= {criteria.get('min_oral_f_pct', 20.0)}%",
        f"- **PNN Remodeling Efficacy**: Hypothalamic WFA+ area reduction >= {criteria.get('min_wfa_reduction_pct', 15.0)}%",
        f"- **Insulin Sensitivity**: ITT glucose clearance AUC improvement >= {criteria.get('min_itt_improvement_pct', 20.0)}%",
        "",
        "## Preclinical Candidate Shortlist",
        "",
        "| Compound ID | Status | Brain Kp | Oral F (%) | Plasma Cmax (ng/mL) | Half-life (h) | WFA Red (%) | ITT Imp (%) | Findings |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for c in clinical + advanced:
        kp_str = f"{c.brain_plasma_ratio:.2f}" if c.brain_plasma_ratio is not None else "-"
        f_str = f"{c.oral_bioavailability_pct:.1f}%" if c.oral_bioavailability_pct is not None else "-"
        cmax_str = f"{c.plasma_cmax_ng_ml:.1f}" if c.plasma_cmax_ng_ml is not None else "-"
        t12_str = f"{c.half_life_h:.1f}" if c.half_life_h is not None else "-"
        wfa_str = f"{c.wfa_pnn_reduction_pct:.1f}%" if c.wfa_pnn_reduction_pct is not None else "-"
        itt_str = f"{c.itt_glucose_auc_reduction_pct:.1f}%" if c.itt_glucose_auc_reduction_pct is not None else "-"
        lines.append(
            f"| `{c.compound_id}` | **{c.status}** | {kp_str} | {f_str} | {cmax_str} | {t12_str} | {wfa_str} | {itt_str} | {c.findings} |"
        )

    lines.extend([
        "",
        "## Next Steps for Nominated Candidates",
        "",
        "1. **Formulation & Salt Selection**: Perform salt screening and crystal polymorph stability studies to optimize solid-state oral formulation.",
        "2. **Dose-Ranging Repeat Toxicology**: Complete 14-day and 28-day repeated-dose toxicity in rodents and non-rodents to establish NOAEL.",
        "3. **Regulatory IND-Enabling Package**: Prepare documentation for Investigational New Drug (IND) filing.",
        "",
    ])
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.preclinical",
        description="Ingest in vivo PK/PD time-series, perform non-compartmental analysis, and qualify preclinical candidates.",
    )
    parser.add_argument("--pk", required=True, help="path to in vivo PK concentration-time series CSV")
    parser.add_argument("--pd", help="optional path to in vivo PD and metabolic outcomes CSV")
    parser.add_argument("--out", default="preclinical_candidates.csv", help="output destination CSV")
    parser.add_argument("--json", dest="json_path", help="optional output path for JSON summary")
    parser.add_argument("--report", dest="report_path", help="optional output path for Markdown qualification report")
    parser.add_argument("--min-kp", type=float, default=0.30, help="minimum brain-to-plasma ratio (default: 0.30)")
    parser.add_argument("--min-f", type=float, default=20.0, help="minimum oral bioavailability pct (default: 20.0)")
    parser.add_argument("--min-wfa", type=float, default=15.0, help="minimum WFA net reduction pct (default: 15.0)")
    parser.add_argument("--min-itt", type=float, default=20.0, help="minimum ITT glucose AUC reduction pct (default: 20.0)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        pk_profiles = load_pk_data(Path(args.pk))
        pd_outcomes = load_pd_data(Path(args.pd)) if args.pd else None

        criteria = {
            "min_kp_brain": args.min_kp,
            "min_oral_f_pct": args.min_f,
            "min_wfa_reduction_pct": args.min_wfa,
            "min_itt_improvement_pct": args.min_itt,
        }
        candidates = evaluate_candidates(
            pk_profiles=pk_profiles,
            pd_outcomes=pd_outcomes,
            min_kp_brain=args.min_kp,
            min_oral_f_pct=args.min_f,
            min_wfa_reduction_pct=args.min_wfa,
            min_itt_improvement_pct=args.min_itt,
        )

        out_path = Path(args.out)
        write_csv(out_path, PRECLINICAL_FIELDS, [c.as_row() for c in candidates])

        if args.json_path:
            summary = {
                "total_evaluated": len(candidates),
                "by_status": {
                    status: sum(1 for c in candidates if c.status == status)
                    for status in ("CLINICAL_CANDIDATE", "ADVANCED_LEAD", "INSUFFICIENT_BRAIN_PENETRATION", "POOR_BIOAVAILABILITY", "INSUFFICIENT_EFFICACY")
                },
                "criteria": criteria,
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), summary)

        if args.report_path:
            report_text = render_preclinical_report(candidates, criteria)
            write_text(Path(args.report_path), report_text)

    except PreclinicalError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: preclinical candidate evaluation requires valid PK and PD CSV datasets.", file=sys.stderr)
        print("Fix: verify columns against templates in data/preclinical/ and ensure non-negative numeric metrics.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: preclinical output files could not be written to disk.", file=sys.stderr)
        print("Fix: provide a writable --out, --json, or --report directory path.", file=sys.stderr)
        return 1

    clinical_count = sum(1 for c in candidates if c.status == "CLINICAL_CANDIDATE")
    advanced_count = sum(1 for c in candidates if c.status == "ADVANCED_LEAD")
    print(
        f"Preclinical evaluation complete: {clinical_count} clinical candidate(s), "
        f"{advanced_count} advanced lead(s) from {len(candidates)} compound(s) -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
