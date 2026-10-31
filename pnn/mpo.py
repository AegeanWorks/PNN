"""Phase 5 Multiparametric Optimization (MPO) and candidate profiling engine for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .reports import write_csv, write_json, write_text


class MpoError(ValueError):
    """Raised when MPO input data is invalid or missing required metrics."""


@dataclass(frozen=True)
class MpoProfile:
    compound_id: str
    name: str
    cns_mpo_score: float  # 0.0 - 6.0
    composite_desirability: float  # 0.0 - 1.0
    mpo_tier: str  # DEVELOPMENT_CANDIDATE, OPTIMIZATION_LEAD, DEPRIORITIZED
    mw: float
    logp: float
    tpsa: float
    hbd: int
    pka: Optional[float]
    ec50_nm: Optional[float]
    brain_kp: Optional[float]
    oral_f_pct: Optional[float]
    selectivity_index: Optional[float]
    potency_desirability: float
    brain_desirability: float
    oral_desirability: float
    safety_desirability: float
    primary_liability: str
    optimization_recommendation: str

    def as_row(self) -> Dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "name": self.name,
            "cns_mpo_score": f"{self.cns_mpo_score:.2f}",
            "composite_desirability": f"{self.composite_desirability:.3f}",
            "mpo_tier": self.mpo_tier,
            "potency_desirability": f"{self.potency_desirability:.2f}",
            "brain_desirability": f"{self.brain_desirability:.2f}",
            "oral_desirability": f"{self.oral_desirability:.2f}",
            "safety_desirability": f"{self.safety_desirability:.2f}",
            "primary_liability": self.primary_liability,
            "optimization_recommendation": self.optimization_recommendation,
        }


MPO_FIELDS = (
    "compound_id",
    "name",
    "cns_mpo_score",
    "composite_desirability",
    "mpo_tier",
    "potency_desirability",
    "brain_desirability",
    "oral_desirability",
    "safety_desirability",
    "primary_liability",
    "optimization_recommendation",
)


def _ramp_desirability(val: float, low: float, high: float, increasing: bool = True) -> float:
    """Calculate linear ramp desirability between low and high bounds."""
    if high <= low:
        return 1.0 if val >= high else 0.0
    if val <= low:
        return 0.0 if increasing else 1.0
    if val >= high:
        return 1.0 if increasing else 0.0
    frac = (val - low) / (high - low)
    return frac if increasing else (1.0 - frac)


def calculate_cns_mpo(
    mw: float,
    logp: float,
    tpsa: float,
    hbd: int,
    pka: Optional[float] = None,
    logd: Optional[float] = None,
) -> float:
    """
    Calculate Pfizer CNS MPO score (0.0 to 6.0) from physicochemical properties.
    Uses continuous piecewise-linear desirability functions.
    """
    # 1. MW: <= 360 -> 1.0; 360-500 -> 1.0 to 0.0; > 500 -> 0.0
    d_mw = _ramp_desirability(mw, low=360.0, high=500.0, increasing=False)

    # 2. logP: <= 3.0 -> 1.0; 3.0-5.0 -> 1.0 to 0.0; > 5.0 -> 0.0
    d_logp = _ramp_desirability(logp, low=3.0, high=5.0, increasing=False)

    # 3. logD7.4: <= 2.0 -> 1.0; 2.0-4.0 -> 1.0 to 0.0; > 4.0 -> 0.0 (fallback to logp - 0.5)
    eff_logd = logd if logd is not None else (logp - 0.5)
    d_logd = _ramp_desirability(eff_logd, low=2.0, high=4.0, increasing=False)

    # 4. TPSA: 40-90 -> 1.0; < 40 -> 0.0 to 1.0 (from 20); 90-120 -> 1.0 to 0.0; > 120 -> 0.0
    if tpsa < 40.0:
        d_tpsa = _ramp_desirability(tpsa, low=20.0, high=40.0, increasing=True)
    elif tpsa <= 90.0:
        d_tpsa = 1.0
    else:
        d_tpsa = _ramp_desirability(tpsa, low=90.0, high=120.0, increasing=False)

    # 5. HBD: <= 0.5 -> 1.0; 0.5-3.5 -> 1.0 to 0.0; > 3.5 -> 0.0
    d_hbd = _ramp_desirability(float(hbd), low=0.5, high=3.5, increasing=False)

    # 6. pKa: <= 8.0 -> 1.0; 8.0-10.0 -> 1.0 to 0.0; > 10.0 -> 0.0 (fallback: neutral=7.4)
    eff_pka = pka if pka is not None else 7.4
    d_pka = _ramp_desirability(eff_pka, low=8.0, high=10.0, increasing=False)

    return d_mw + d_logp + d_logd + d_tpsa + d_hbd + d_pka


def profile_candidate(
    compound_id: str,
    name: str,
    mw: float,
    logp: float,
    tpsa: float,
    hbd: int,
    pka: Optional[float] = None,
    ec50_nm: Optional[float] = None,
    brain_kp: Optional[float] = None,
    oral_f_pct: Optional[float] = None,
    selectivity_index: Optional[float] = None,
) -> MpoProfile:
    """Calculate CNS MPO, multi-domain desirability functions, and identify optimization liabilities."""
    cns_mpo = calculate_cns_mpo(mw=mw, logp=logp, tpsa=tpsa, hbd=hbd, pka=pka)

    # Individual desirability components
    # Potency: EC50 <= 25 nM -> 1.0; >= 500 nM -> 0.0
    eff_ec50 = ec50_nm if ec50_nm is not None else 250.0
    d_potency = _ramp_desirability(eff_ec50, low=25.0, high=500.0, increasing=False)

    # Brain Partition Kp: >= 0.50 -> 1.0; <= 0.15 -> 0.0
    eff_kp = brain_kp if brain_kp is not None else 0.30
    d_brain = _ramp_desirability(eff_kp, low=0.15, high=0.50, increasing=True)

    # Oral Bioavailability F: >= 30% -> 1.0; <= 10% -> 0.0
    eff_f = oral_f_pct if oral_f_pct is not None else 20.0
    d_oral = _ramp_desirability(eff_f, low=10.0, high=30.0, increasing=True)

    # Safety Margin (Selectivity Index): >= 50 -> 1.0; <= 5 -> 0.0
    eff_si = selectivity_index if selectivity_index is not None else 25.0
    d_safety = _ramp_desirability(eff_si, low=5.0, high=50.0, increasing=True)

    # Derringer-Suich geometric composite desirability (avoid zero multiplication crash)
    vals = [max(0.01, v) for v in (d_potency, d_brain, d_oral, d_safety)]
    composite = math.pow(vals[0] * vals[1] * vals[2] * vals[3], 0.25)

    # Diagnostic liabilities & recommendations
    liabilities: List[str] = []
    recommendations: List[str] = []

    if cns_mpo < 4.5:
        if mw > 380.0:
            liabilities.append(f"High molecular weight ({mw:.1f} Da)")
            recommendations.append("Truncate peripheral side-chains to improve ligand efficiency.")
        if logp > 3.2:
            liabilities.append(f"Excessive lipophilicity (logP {logp:.2f})")
            recommendations.append("Incorporate polar heteroatoms or fluorinated motifs to control logD.")
        if tpsa > 75.0:
            liabilities.append(f"High polar surface area ({tpsa:.1f} Å²)")
            recommendations.append("Cap hydrogen bond donors to enhance blood-brain barrier passive transport.")

    if pka is not None and pka > 8.5:
        liabilities.append(f"High basic center pKa ({pka:.1f})")
        recommendations.append("Modulate amine pKa to 7.5-8.2 with adjacent electron-withdrawing groups.")

    if brain_kp is not None and brain_kp < 0.30:
        liabilities.append(f"Sub-optimal brain partition (Kp={brain_kp:.2f})")
        recommendations.append("Reduce rotatable bonds and optimize TPSA to boost brain uptake.")

    if oral_f_pct is not None and oral_f_pct < 20.0:
        liabilities.append(f"Marginal oral bioavailability (F={oral_f_pct:.1f}%)")
        recommendations.append("Screen oral salt forms (mesylate, besylate) or evaluate prodrug formulations.")

    primary_liability = liabilities[0] if liabilities else "Balanced physicochemical and exposure profile"
    rec_text = " ".join(recommendations) if recommendations else "Advance to dose-ranging toxicology and formulation development."

    # Tier assignment
    if composite >= 0.70 and cns_mpo >= 4.5:
        tier = "DEVELOPMENT_CANDIDATE"
    elif composite >= 0.45:
        tier = "OPTIMIZATION_LEAD"
    else:
        tier = "DEPRIORITIZED"

    return MpoProfile(
        compound_id=compound_id,
        name=name,
        cns_mpo_score=cns_mpo,
        composite_desirability=composite,
        mpo_tier=tier,
        mw=mw,
        logp=logp,
        tpsa=tpsa,
        hbd=hbd,
        pka=pka,
        ec50_nm=ec50_nm,
        brain_kp=brain_kp,
        oral_f_pct=oral_f_pct,
        selectivity_index=selectivity_index,
        potency_desirability=d_potency,
        brain_desirability=d_brain,
        oral_desirability=d_oral,
        safety_desirability=d_safety,
        primary_liability=primary_liability,
        optimization_recommendation=rec_text,
    )


def render_mpo_report(profiles: Sequence[MpoProfile]) -> str:
    """Render comprehensive Phase 5 Lead Optimization & MPO Profiling markdown report."""
    dev_candidates = [p for p in profiles if p.mpo_tier == "DEVELOPMENT_CANDIDATE"]
    opt_leads = [p for p in profiles if p.mpo_tier == "OPTIMIZATION_LEAD"]

    lines = [
        "# PNN Phase 5 Lead Optimization & Multiparametric Profiling (MPO) Report",
        "",
        "## Executive Summary",
        "",
        f"- **Total Compounds Profiled**: {len(profiles)}",
        f"- **Development Candidates (Composite D >= 0.70, CNS MPO >= 4.5)**: {len(dev_candidates)}",
        f"- **Optimization Leads (Composite D >= 0.45)**: {len(opt_leads)}",
        f"- **Deprioritized Compounds**: {len(profiles) - len(dev_candidates) - len(opt_leads)}",
        "",
        "## Multiparametric Optimization Taxonomy",
        "",
        "> [!NOTE]",
        "> **Pfizer CNS MPO & Derringer-Suich Desirability**: CNS MPO evaluates six central nervous system drug-likeness "
        "axes (MW, logP, logD, TPSA, HBD, pKa) on a 0.0 to 6.0 scale. Geometric composite desirability (D) weights target "
        "potency (DQ-CS EC50), brain partitioning (Kp,brain), oral bioavailability (F), and cytotoxicity safety margin (SI).",
        "",
        "## Profiled Candidates Shortlist",
        "",
        "| Compound ID | Chemotype Name | CNS MPO | Composite D | Status Tier | Primary Chemical Liability | Optimization Recommendation |",
        "|---|---|---|---|---|---|---|",
    ]

    for p in profiles:
        lines.append(
            f"| `{p.compound_id}` | **{p.name}** | {p.cns_mpo_score:.2f} / 6.0 | "
            f"**{p.composite_desirability:.3f}** | *{p.mpo_tier}* | {p.primary_liability} | {p.optimization_recommendation} |"
        )

    lines.extend([
        "",
        "## Key Medicinal Chemistry Action Items",
        "",
        "1. **Basic Nitrogen pKa Fine-Tuning**: Maintain pKa in the 7.5 - 8.2 window to maximize intracellular lysosomal penetration while avoiding systemic clearance.",
        "2. **Polar Surface Area Shielding**: Target TPSA between 45 and 65 Å² with minimal HBD count (<= 1) for optimal blood-brain barrier passive transport.",
        "3. **Solid-State Formulation & Polymorph Screen**: Advance nominated development candidates to crystalline salt screening for oral solid-dosage formulation stability.",
        "",
    ])
    return "\n".join(lines)


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


def run_mpo_analysis(
    compounds_path: Path,
    candidates_path: Optional[Path] = None,
    assays_path: Optional[Path] = None,
) -> List[MpoProfile]:
    """Ingest datasets, compute MPO and composite desirability, and generate profiles."""
    if not compounds_path.is_file():
        raise MpoError(f"Compounds CSV not found: {compounds_path}")

    # Load compounds
    raw_compounds: Dict[str, Dict[str, str]] = {}
    with open(compounds_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise MpoError(f"Empty compounds CSV: {compounds_path}")
        norm = {c.strip().lower(): c for c in reader.fieldnames}
        cid_col = norm.get("compound_id")
        if not cid_col:
            raise MpoError(f"Compounds CSV missing compound_id: {compounds_path}")
        for row in reader:
            cid = row[cid_col].strip()
            if cid:
                raw_compounds[cid] = {k: row[v].strip() for k, v in norm.items()}

    # Load candidates if available
    cand_data: Dict[str, Dict[str, float]] = {}
    if candidates_path and candidates_path.is_file():
        with open(candidates_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                norm = {c.strip().lower(): c for c in reader.fieldnames}
                cid_col = norm.get("compound_id")
                if cid_col:
                    for row in reader:
                        cid = row[cid_col].strip()
                        c_dict: Dict[str, float] = {}
                        if "brain_plasma_ratio" in norm:
                            val = _parse_opt_float(row[norm["brain_plasma_ratio"]])
                            if val is not None:
                                c_dict["brain_kp"] = val
                        if "oral_bioavailability_pct" in norm:
                            val = _parse_opt_float(row[norm["oral_bioavailability_pct"]])
                            if val is not None:
                                c_dict["oral_f"] = val
                        if c_dict:
                            cand_data[cid] = c_dict

    # Load assays if available
    assay_data: Dict[str, Dict[str, float]] = {}
    if assays_path and assays_path.is_file():
        with open(assays_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                norm = {c.strip().lower(): c for c in reader.fieldnames}
                cid_col = norm.get("compound_id")
                if cid_col:
                    for row in reader:
                        cid = row[cid_col].strip()
                        a_dict: Dict[str, float] = {}
                        if "ec50_nm" in norm or "pnn_ec50_nm" in norm:
                            col = norm.get("ec50_nm") or norm.get("pnn_ec50_nm")
                            val = _parse_opt_float(row[col])  # type: ignore[arg-type]
                            if val is not None:
                                a_dict["ec50_nm"] = val
                        if "selectivity_index" in norm:
                            val = _parse_opt_float(row[norm["selectivity_index"]])
                            if val is not None:
                                a_dict["si"] = val
                        if a_dict:
                            assay_data[cid] = a_dict

    profiles: List[MpoProfile] = []
    for cid, row in raw_compounds.items():
        name = row.get("name", cid)
        mw = _parse_opt_float(row.get("mw"))
        logp = _parse_opt_float(row.get("logp"))
        tpsa = _parse_opt_float(row.get("tpsa"))
        hbd = _parse_opt_float(row.get("hbd"))
        if None in (mw, logp, tpsa, hbd):
            raise MpoError(
                f"Compound {cid} lacks mw/logp/tpsa/hbd required for MPO profiling"
            )
        hbd = int(hbd)
        pka = _parse_opt_float(row.get("pka"))

        # In vitro & in vivo properties
        ec50 = assay_data.get(cid, {}).get("ec50_nm") or _parse_opt_float(row.get("pnn_ec50_nm"))
        brain_kp = cand_data.get(cid, {}).get("brain_kp") or _parse_opt_float(row.get("brain_plasma_ratio"))
        oral_f = cand_data.get(cid, {}).get("oral_f")
        if oral_f is None and "oral_bioavailability" in row:
            f_frac = _parse_opt_float(row.get("oral_bioavailability"))
            if f_frac is not None:
                oral_f = f_frac * 100.0

        si = assay_data.get(cid, {}).get("si")

        prof = profile_candidate(
            compound_id=cid,
            name=name,
            mw=mw,
            logp=logp,
            tpsa=tpsa,
            hbd=hbd,
            pka=pka,
            ec50_nm=ec50,
            brain_kp=brain_kp,
            oral_f_pct=oral_f,
            selectivity_index=si,
        )
        profiles.append(prof)

    # Sort profiles by composite desirability descending
    profiles.sort(key=lambda p: -p.composite_desirability)
    return profiles


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.mpo",
        description="Phase 5 Multiparametric Optimization (MPO) and candidate profiling engine.",
    )
    parser.add_argument("--compounds", required=True, help="path to compounds library CSV")
    parser.add_argument("--candidates", help="optional path to preclinical candidates CSV")
    parser.add_argument("--assays", help="optional path to in vitro assays or qualified hits CSV")
    parser.add_argument("--out", default="mpo_profiles.csv", help="output destination CSV")
    parser.add_argument("--report", dest="report_path", help="optional destination path for Markdown report")
    parser.add_argument("--json", dest="json_path", help="optional destination path for JSON summary")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        profiles = run_mpo_analysis(
            compounds_path=Path(args.compounds),
            candidates_path=Path(args.candidates) if args.candidates else None,
            assays_path=Path(args.assays) if args.assays else None,
        )

        out_path = Path(args.out)
        write_csv(out_path, MPO_FIELDS, [p.as_row() for p in profiles])

        if args.report_path:
            report_text = render_mpo_report(profiles)
            write_text(Path(args.report_path), report_text)

        if args.json_path:
            summary = {
                "total_profiled": len(profiles),
                "by_tier": {
                    tier: sum(1 for p in profiles if p.mpo_tier == tier)
                    for tier in ("DEVELOPMENT_CANDIDATE", "OPTIMIZATION_LEAD", "DEPRIORITIZED")
                },
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), summary)

    except MpoError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: MPO analysis requires valid compound and exposure parameters.", file=sys.stderr)
        print("Fix: verify that the compound dataset exists and contains required properties.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: MPO output files could not be written to disk.", file=sys.stderr)
        print("Fix: specify a writable destination for --out, --report, or --json.", file=sys.stderr)
        return 1

    dev_count = sum(1 for p in profiles if p.mpo_tier == "DEVELOPMENT_CANDIDATE")
    opt_count = sum(1 for p in profiles if p.mpo_tier == "OPTIMIZATION_LEAD")
    print(
        f"MPO profiling complete: {dev_count} development candidate(s), {opt_count} optimization lead(s) "
        f"across {len(profiles)} compound(s) -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
