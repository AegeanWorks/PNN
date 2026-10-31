"""Allometric human clinical dose projection and IND regulatory dossier compiler for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from .reports import write_csv, write_json, write_text
from .safety import KM_BY_SPECIES, load_tox_studies, most_sensitive_study

KM_MOUSE = 3.0
KM_HUMAN = 37.0
MOUSE_WEIGHT_KG = 0.020
HUMAN_WEIGHT_KG = 60.0
DEFAULT_SAFETY_FACTOR = 10.0


class IndError(ValueError):
    """Raised when clinical projection inputs violate regulatory or pharmacological contracts."""


def _fmt(value: Optional[float], spec: str) -> str:
    return format(value, spec) if value is not None else ""


def _display(value: Optional[float], spec: str) -> str:
    return format(value, spec) if value is not None else "UNMEASURED"


@dataclass(frozen=True)
class HumanDoseProjection:
    compound_id: str
    name: str
    animal_dose_mg_kg: Optional[float]
    human_equivalent_dose_mg_kg: Optional[float]
    human_total_daily_dose_mg: Optional[float]
    projected_human_cl_l_h: Optional[float]
    projected_human_vd_l: Optional[float]
    projected_human_half_life_h: Optional[float]
    noael_mg_kg: Optional[float]
    noael_hed_mg_kg: Optional[float]
    phase1_mrsd_mg: Optional[float]
    therapeutic_index: Optional[float]
    brain_kp: Optional[float]
    oral_f_pct: Optional[float]
    findings: str
    evidence_status: str  # COMPLETE, PARTIAL or INSUFFICIENT_EVIDENCE

    def as_row(self) -> Dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "name": self.name,
            "animal_dose_mg_kg": _fmt(self.animal_dose_mg_kg, ".1f"),
            "human_equivalent_dose_mg_kg": _fmt(self.human_equivalent_dose_mg_kg, ".3f"),
            "human_total_daily_dose_mg": _fmt(self.human_total_daily_dose_mg, ".1f"),
            "projected_human_cl_l_h": _fmt(self.projected_human_cl_l_h, ".2f"),
            "projected_human_vd_l": _fmt(self.projected_human_vd_l, ".1f"),
            "projected_human_half_life_h": _fmt(self.projected_human_half_life_h, ".1f"),
            "noael_mg_kg": _fmt(self.noael_mg_kg, ".1f"),
            "noael_hed_mg_kg": _fmt(self.noael_hed_mg_kg, ".3f"),
            "phase1_mrsd_mg": _fmt(self.phase1_mrsd_mg, ".3f"),
            "therapeutic_index": _fmt(self.therapeutic_index, ".0f"),
            "brain_kp": _fmt(self.brain_kp, ".2f"),
            "oral_f_pct": _fmt(self.oral_f_pct, ".1f"),
            "findings": self.findings,
            "evidence_status": self.evidence_status,
        }


# Backwards compatibility alias
ClinicalDoseProjection = HumanDoseProjection


DOSE_PROJECTION_FIELDS = (
    "compound_id",
    "name",
    "animal_dose_mg_kg",
    "human_equivalent_dose_mg_kg",
    "human_total_daily_dose_mg",
    "projected_human_cl_l_h",
    "projected_human_vd_l",
    "projected_human_half_life_h",
    "noael_mg_kg",
    "noael_hed_mg_kg",
    "phase1_mrsd_mg",
    "therapeutic_index",
    "brain_kp",
    "oral_f_pct",
    "findings",
    "evidence_status",
)


def calculate_allometric_scaling(
    animal_dose_mg_kg: float,
    mouse_cl_l_h_kg: float,
    mouse_vd_l_kg: float,
    noael_mg_kg: float,
    noael_km: float = KM_MOUSE,
    mouse_weight_kg: float = MOUSE_WEIGHT_KG,
    human_weight_kg: float = HUMAN_WEIGHT_KG,
    km_human: float = KM_HUMAN,
    safety_factor: float = DEFAULT_SAFETY_FACTOR,
) -> Tuple[float, float, float, float, float, float, float, float]:
    """
    Compute FDA / ICH allometric scaling parameters from measured mouse PK and NOAEL data.
    Returns: (hed_mg_kg, human_daily_mg, cl_human_l_h, vd_human_l, t12_human_h, noael_mg_kg, noael_hed, mrsd_mg)

    Every argument is measured input: no dose, clearance, volume or NOAEL is inferred.
    """
    if animal_dose_mg_kg <= 0.0:
        raise IndError("animal_dose_mg_kg must be > 0")
    if mouse_cl_l_h_kg <= 0.0:
        raise IndError("mouse_cl_l_h_kg must be > 0")
    if mouse_vd_l_kg <= 0.0:
        raise IndError("mouse_vd_l_kg must be > 0")
    if noael_mg_kg <= 0.0:
        raise IndError("noael_mg_kg must be > 0")
    if safety_factor <= 0.0:
        raise IndError("safety_factor must be > 0")
    if mouse_weight_kg <= 0.0 or human_weight_kg <= 0.0:
        raise IndError("body weights must be > 0")

    # 1. Human Equivalent Dose (HED) from the NOAEL species body-surface-area factor.
    noael_factor = noael_km / km_human
    hed_mg_kg = animal_dose_mg_kg * (KM_MOUSE / km_human)
    human_daily_mg = hed_mg_kg * human_weight_kg

    # 2. PK allometric scaling from measured mouse CL and Vd.
    weight_ratio = human_weight_kg / mouse_weight_kg
    cl_human = mouse_cl_l_h_kg * mouse_weight_kg * math.pow(weight_ratio, 0.75)
    vd_human = mouse_vd_l_kg * mouse_weight_kg * math.pow(weight_ratio, 1.0)
    t12_human = (math.log(2.0) * vd_human) / cl_human

    # 3. NOAEL HED and Phase 1 Maximum Recommended Starting Dose (MRSD).
    noael_hed = noael_mg_kg * noael_factor
    mrsd_mg = (noael_hed / safety_factor) * human_weight_kg

    return hed_mg_kg, human_daily_mg, cl_human, vd_human, t12_human, noael_mg_kg, noael_hed, mrsd_mg


def _project_one(
    cid: str,
    name: str,
    info: Mapping[str, float],
    safety_factor: float,
    human_weight_kg: float,
) -> HumanDoseProjection:
    dose = info.get("dose_mg_kg")
    cl_mouse = info.get("cl_l_h_kg")
    vd_mouse = info.get("vd_l_kg")
    noael = info.get("noael_mg_kg")
    noael_km = info.get("noael_km", KM_MOUSE)
    kp = info.get("brain_kp")
    oral_f = info.get("oral_f")

    findings: List[str] = []
    measured_fields = 0

    hed = daily_mg = None
    if dose is not None and dose > 0.0:
        hed = dose * (KM_MOUSE / KM_HUMAN)
        daily_mg = hed * human_weight_kg
        measured_fields += 1
    else:
        dose = None
        findings.append("Efficacious animal dose unmeasured")

    cl_h = vd_h = t12_h = None
    if cl_mouse is not None and vd_mouse is not None and cl_mouse > 0.0 and vd_mouse > 0.0:
        weight_ratio = human_weight_kg / MOUSE_WEIGHT_KG
        cl_h = cl_mouse * MOUSE_WEIGHT_KG * math.pow(weight_ratio, 0.75)
        vd_h = vd_mouse * MOUSE_WEIGHT_KG * weight_ratio
        t12_h = (math.log(2.0) * vd_h) / cl_h
        measured_fields += 1
    else:
        findings.append("Mouse CL/Vd unmeasured; human half-life not projected")

    noael_hed = mrsd_mg = ti = None
    if noael is not None and noael > 0.0:
        noael_hed = noael * (noael_km / KM_HUMAN)
        mrsd_mg = (noael_hed / safety_factor) * human_weight_kg
        measured_fields += 1
        if dose is not None:
            ti = noael / dose
    else:
        findings.append("NOAEL unmeasured; MRSD and therapeutic index not projected")

    if kp is not None:
        measured_fields += 1
        if kp >= 0.40:
            findings.append(f"Robust central brain exposure (Kp={kp:.2f})")
        elif kp >= 0.30:
            findings.append(f"Acceptable central brain exposure (Kp={kp:.2f})")
        else:
            findings.append(f"Restricted central brain exposure (Kp={kp:.2f})")
    else:
        findings.append("Brain partitioning unmeasured")

    if oral_f is not None:
        measured_fields += 1
        if oral_f >= 20.0:
            findings.append(f"Oral bioavailability above the progression floor (F={oral_f:.1f}%)")
        else:
            findings.append(f"Oral bioavailability below the progression floor (F={oral_f:.1f}%)")
    else:
        findings.append("Oral bioavailability unmeasured")

    if t12_h is not None:
        if t12_h >= 12.0:
            findings.append(f"Extended projected human t1/2 ({t12_h:.1f}h) supports once-daily dosing")
        else:
            findings.append(f"Moderate projected human t1/2 ({t12_h:.1f}h) suggests pulsatile dosing")

    if ti is not None and ti >= 5.0:
        findings.append(f"Broad therapeutic window (TI={ti:.1f}x)")

    required_fields = 5  # dose, mouse CL/Vd, NOAEL, brain Kp, oral F
    if measured_fields >= required_fields:
        status = "COMPLETE"
    elif measured_fields > 0:
        status = "PARTIAL"
    else:
        status = "INSUFFICIENT_EVIDENCE"

    return HumanDoseProjection(
        compound_id=cid,
        name=name,
        animal_dose_mg_kg=dose,
        human_equivalent_dose_mg_kg=hed,
        human_total_daily_dose_mg=daily_mg,
        projected_human_cl_l_h=cl_h,
        projected_human_vd_l=vd_h,
        projected_human_half_life_h=t12_h,
        noael_mg_kg=noael,
        noael_hed_mg_kg=noael_hed,
        phase1_mrsd_mg=mrsd_mg,
        therapeutic_index=ti,
        brain_kp=kp,
        oral_f_pct=oral_f,
        findings="; ".join(findings),
        evidence_status=status,
    )


def project_clinical_candidates(
    compounds: Sequence[Dict[str, str]],
    preclinical_data: Optional[Dict[str, Dict[str, float]]] = None,
    lead_id: Optional[str] = None,
    safety_factor: float = DEFAULT_SAFETY_FACTOR,
    human_weight_kg: float = HUMAN_WEIGHT_KG,
) -> List[HumanDoseProjection]:
    """Calculate clinical dose projections from measured preclinical data only."""
    pk_map = preclinical_data or {}
    projections: List[HumanDoseProjection] = []

    for row in compounds:
        cid = row.get("compound_id", "").strip()
        name = row.get("name", cid).strip() or cid
        if not cid:
            continue
        projections.append(
            _project_one(cid, name, pk_map.get(cid, {}), safety_factor, human_weight_kg)
        )

    evidence_rank = {"COMPLETE": 0, "PARTIAL": 1, "INSUFFICIENT_EVIDENCE": 2}
    projections.sort(
        key=lambda p: (
            0 if lead_id and p.compound_id == lead_id else 1,
            evidence_rank.get(p.evidence_status, 3),
            -p.brain_kp if p.brain_kp is not None else 0.0,
            p.compound_id,
        )
    )
    return projections


def render_ind_dossier(
    projections: Sequence[HumanDoseProjection],
    lead_id: Optional[str] = None,
    criteria: Optional[Mapping[str, object]] = None,
) -> str:
    """Render an IND nomination dossier; absent values are shown as UNMEASURED, never filled in."""
    criteria = criteria or {}
    lead = next((p for p in projections if lead_id and p.compound_id == lead_id), None)
    if lead is None:
        lead = projections[0] if projections else None

    lines = [
        "# Investigational New Drug (IND) Candidate Nomination Dossier",
        "",
        "## Module 2.4: Nonclinical Overview & Clinical Starting Dose Justification",
        "",
        "> [!IMPORTANT]",
        "> **Regulatory Guidance (ICH M3(R2) / S7A / S7B)**: Allometric scaling estimates Human Equivalent Dose (HED) "
        "and Maximum Recommended Starting Dose (MRSD) from measured animal dose, clearance, volume and NOAEL data. "
        "First-in-Human (FIH) Phase 1 initiation applies a 10-fold safety factor to the NOAEL HED. "
        "Values marked UNMEASURED were not supplied and must not be inferred.",
        "",
        "## Clinical Candidate Nomination Summary",
        "",
    ]

    if lead is None:
        lines.append("- **No candidate records were supplied.**")
    else:
        lines.extend(
            [
                f"- **Nominated Lead Compound**: `{lead.compound_id}` (**{lead.name}**)",
                f"- **Evidence Status**: **{lead.evidence_status}**",
                f"- **Rodent Efficacious Dose**: {_display(lead.animal_dose_mg_kg, '.1f')}"
                + (" mg/kg PO" if lead.animal_dose_mg_kg is not None else ""),
                f"- **Human Equivalent Dose (HED)**: {_display(lead.human_equivalent_dose_mg_kg, '.3f')}"
                + (f" mg/kg ({lead.human_total_daily_dose_mg:.1f} mg total daily for 60 kg adult)"
                   if lead.human_total_daily_dose_mg is not None else ""),
                f"- **Projected Human Half-Life ($t_{{1/2}}$)**: {_display(lead.projected_human_half_life_h, '.1f')}"
                + (f" hours (Projected CL: {lead.projected_human_cl_l_h:.2f} L/h, Vd: {lead.projected_human_vd_l:.1f} L)"
                   if lead.projected_human_half_life_h is not None else ""),
                f"- **No Observed Adverse Effect Level (NOAEL)**: {_display(lead.noael_mg_kg, '.1f')}"
                + (f" mg/kg (NOAEL HED: {lead.noael_hed_mg_kg:.3f} mg/kg)"
                   if lead.noael_hed_mg_kg is not None else ""),
                f"- **Phase 1 First-in-Human Starting Dose (MRSD)**: {_display(lead.phase1_mrsd_mg, '.2f')}"
                + (" mg (NOAEL HED / 10)" if lead.phase1_mrsd_mg is not None else ""),
                f"- **Therapeutic Window**: {_display(lead.therapeutic_index, '.1f')}"
                + ("-fold over the efficacious dose" if lead.therapeutic_index is not None else ""),
                f"- **Brain Partition ($K_{{p,\\text{{brain}}}}$)**: {_display(lead.brain_kp, '.2f')}",
                f"- **Oral Bioavailability ($F$)**: {_display(lead.oral_f_pct, '.1f')}"
                + ("%" if lead.oral_f_pct is not None else ""),
                "",
            ]
        )

    lines.extend(
        [
            "## Candidate Comparison & Scaling Matrix",
            "",
            "| Compound ID | Chemotype | Evidence | Animal Dose (mg/kg) | Human Dose (mg/day) | Human t1/2 (h) | NOAEL HED (mg/kg) | Phase 1 MRSD (mg) | Safety Margin | Key Findings |",
            "|---|---|---|---|---|---|---|---|---|---|",
        ]
    )

    for p in projections:
        lines.append(
            f"| `{p.compound_id}` | **{p.name}** | {p.evidence_status} | "
            f"{_display(p.animal_dose_mg_kg, '.1f')} | {_display(p.human_total_daily_dose_mg, '.1f')} | "
            f"{_display(p.projected_human_half_life_h, '.1f')} | {_display(p.noael_hed_mg_kg, '.3f')} | "
            f"{_display(p.phase1_mrsd_mg, '.2f')} | "
            f"{(_display(p.therapeutic_index, '.1f') + 'x') if p.therapeutic_index is not None else 'UNMEASURED'} | "
            f"{p.findings} |"
        )

    incomplete = [p for p in projections if p.evidence_status != "COMPLETE"]
    lines.extend(["", "## Evidence Gaps", ""])
    if incomplete:
        for p in incomplete:
            lines.append(f"- `{p.compound_id}` ({p.evidence_status}): {p.findings}")
        lines.append(
            "- Supply the missing measurements from `pnn.preclinical` (dose, CL, Vd, Kp, F) and "
            "`pnn.safety` (most sensitive NOAEL) before this dossier is used for filing."
        )
    else:
        lines.append(
            "- All projection inputs were supplied. Confirm study quality, GLP compliance and exposure "
            "margins with a qualified toxicologist before submission."
        )

    lines.extend(["", "## Safety Pharmacology Checklist (ICH S7A/S7B)", ""])
    lines.extend(
        [
            "> These items are filing requirements, not completed findings. `python -m pnn.safety` records the "
            "submitted results and progression thresholds; an unchecked item blocks IND filing until evidence exists.",
            "",
            "- [ ] **Cardiovascular Safety (hERG IC50)**: required, with the margin over free Cmax recorded by `pnn.safety`.",
            "- [ ] **Central Nervous System (Irwin / FOB Screen)**: required, with a CNS safety NOAEL on file.",
            "- [ ] **Respiratory Safety**: required, with a respiratory safety NOAEL on file.",
            "- [ ] **Genetic Toxicology**: required; bacterial reverse mutation (Ames) and in vitro micronucleus must be negative.",
            "",
            "## Recommended Phase 1 Clinical Protocol Architecture",
            "",
            "1. **Study Design**: Randomized, double-blind, placebo-controlled Single Ascending Dose (SAD) followed by "
            "14-day Multiple Ascending Dose (MAD) in healthy adult volunteers.",
        ]
    )
    if lead is not None and lead.phase1_mrsd_mg is not None:
        lines.append(
            f"2. **Starting Dose**: Initial cohort starts at **{lead.phase1_mrsd_mg:.2f} mg** oral capsule (fasted state)."
        )
    else:
        lines.append(
            "2. **Starting Dose**: pending. The NOAEL-based MRSD cannot be stated until the NOAEL is measured."
        )
    lines.extend(
        [
            "3. **Dose Escalation**: Modified Fibonacci escalation cohorts (e.g. 1x, 2x, 3.3x, 5x, 7.5x) with real-time "
            "safety monitoring and dense plasma PK sampling.",
            "4. **Biomarker Endpoints**: Lumbar puncture CSF-to-plasma drug partition ($K_{p,\\text{CSF}}$) and metabolic "
            "fasting glucose/insulin clamp.",
            "",
            "## Projection Limits",
            "",
            f"- HED uses the mouse K_m ({KM_MOUSE:g}) against the human K_m ({KM_HUMAN:g}); the NOAEL HED uses the "
            "NOAEL species K_m supplied with the toxicology data.",
            f"- The 10-fold safety factor is the default ({criteria.get('safety_factor', DEFAULT_SAFETY_FACTOR):g}). "
            "Reconcile the MRSD with the most sensitive species NOAEL in `pnn.safety` before filing.",
            "- Allometric scaling is a projection from animal data and never establishes human safety or efficacy.",
            "",
        ]
    )
    return "\n".join(lines)


def _parse_opt_float(val: Optional[str]) -> Optional[float]:
    """Parse an optional numeric float string, returning None for empty or non-numeric cells."""
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


def _as_paths(compounds_path: Union[Path, Sequence[Path]]) -> List[Path]:
    if isinstance(compounds_path, (str, Path)):
        return [Path(compounds_path)]
    return [Path(p) for p in compounds_path]


def load_candidate_records(
    compounds_path: Union[Path, Sequence[Path]],
    candidates_path: Optional[Path] = None,
    tox_path: Optional[Path] = None,
) -> Tuple[List[Dict[str, str]], Dict[str, Dict[str, float]]]:
    """Load compound records and measured preclinical exposures for allometric scaling."""
    records: List[Dict[str, str]] = []
    for path in _as_paths(compounds_path):
        if not path.is_file():
            raise IndError(f"Compounds dataset not found: {path}")
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise IndError(f"Empty compounds CSV: {path}")
            norm = {c.strip().lower(): c for c in reader.fieldnames}
            cid_col = norm.get("compound_id")
            if not cid_col:
                raise IndError(f"Compounds CSV missing compound_id: {path}")
            for row in reader:
                cid = row[cid_col].strip()
                if cid:
                    records.append({k: (row[v] or "").strip() for k, v in norm.items()})

    pk_data: Dict[str, Dict[str, float]] = {}

    if candidates_path is not None and candidates_path.is_file():
        candidate_fields = {
            "brain_plasma_ratio": "brain_kp",
            "oral_bioavailability_pct": "oral_f",
            "dose_mg_kg": "dose_mg_kg",
            "cl_l_h_kg": "cl_l_h_kg",
            "vd_l_kg": "vd_l_kg",
            "noael_mg_kg": "noael_mg_kg",
        }
        with open(candidates_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                norm = {c.strip().lower(): c for c in reader.fieldnames}
                cid_col = norm.get("compound_id")
                if cid_col:
                    for row in reader:
                        cid = row[cid_col].strip()
                        if not cid:
                            continue
                        c_dict: Dict[str, float] = {}
                        for column, key in candidate_fields.items():
                            if column in norm:
                                val = _parse_opt_float(row[norm[column]])
                                if val is not None:
                                    c_dict[key] = val
                        if c_dict:
                            pk_data[cid] = c_dict

    if tox_path is not None and tox_path.is_file():
        by_compound: Dict[str, List] = {}
        for study in load_tox_studies(tox_path):
            by_compound.setdefault(study.compound_id, []).append(study)
        for cid, studies in by_compound.items():
            sensitive = most_sensitive_study(studies)
            if sensitive is None:
                continue
            entry = pk_data.setdefault(cid, {})
            entry["noael_mg_kg"] = sensitive.noael_mg_kg
            entry["noael_km"] = KM_BY_SPECIES[sensitive.species]

    return records, pk_data


def compile_ind_dossier(
    lead_id: str,
    compounds_path: Union[Path, Sequence[Path]],
    candidates_path: Optional[Path] = None,
    tox_path: Optional[Path] = None,
    extra_pk: Optional[Mapping[str, Mapping[str, float]]] = None,
) -> Tuple[List[HumanDoseProjection], str]:
    """Compile human dose projections and the IND dossier from measured inputs only."""
    records, pk_data = load_candidate_records(
        compounds_path=compounds_path,
        candidates_path=candidates_path,
        tox_path=tox_path,
    )
    for cid, values in (extra_pk or {}).items():
        pk_data.setdefault(cid, {}).update(values)

    projections = project_clinical_candidates(records, pk_data, lead_id=lead_id)
    md_text = render_ind_dossier(projections, lead_id=lead_id)
    return projections, md_text


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.ind",
        description="Compute human allometric dose projections and compile an IND regulatory candidate dossier.",
    )
    parser.add_argument("--compounds", required=True, nargs="+", help="candidate compound library CSV path(s)")
    parser.add_argument("--candidates", help="optional preclinical candidates CSV with measured Kp/F/PK columns")
    parser.add_argument("--tox", help="optional repeat-dose toxicology CSV supplying measured NOAEL values")
    parser.add_argument("--lead", help="compound_id to nominate as the lead (default: best evidence)")
    parser.add_argument("--out", default="human_dose_projections.csv", help="output destination CSV for dose projections")
    parser.add_argument("--dossier", dest="dossier_path", default="ind_dossier.md", help="output destination for IND dossier report")
    parser.add_argument("--json", dest="json_path", help="optional output path for JSON summary")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        records, pk_data = load_candidate_records(
            compounds_path=[Path(p) for p in args.compounds],
            candidates_path=Path(args.candidates) if args.candidates else None,
            tox_path=Path(args.tox) if args.tox else None,
        )

        projections = project_clinical_candidates(records, pk_data)
        lead_id = args.lead or (projections[0].compound_id if projections else "")
        if args.lead:
            projections = project_clinical_candidates(records, pk_data, lead_id=lead_id)
        out_path = Path(args.out)
        write_csv(out_path, DOSE_PROJECTION_FIELDS, [p.as_row() for p in projections])

        if args.dossier_path:
            dossier_text = render_ind_dossier(projections, lead_id=lead_id)
            write_text(Path(args.dossier_path), dossier_text)

        if args.json_path:
            lead = next((p for p in projections if p.compound_id == lead_id), None)
            summary = {
                "total_projected": len(projections),
                "lead_candidate_id": lead_id or None,
                "lead_candidate": lead.as_row() if lead else None,
                "by_evidence_status": {
                    status: sum(1 for p in projections if p.evidence_status == status)
                    for status in ("COMPLETE", "PARTIAL", "INSUFFICIENT_EVIDENCE")
                },
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), summary)

    except IndError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: IND allometric projection requires measured dose, clearance, volume and NOAEL inputs.", file=sys.stderr)
        print("Fix: pass --candidates from pnn.preclinical and --tox from pnn.safety data; missing values stay UNMEASURED.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: IND output files could not be written to disk.", file=sys.stderr)
        print("Fix: specify a writable destination for --out, --dossier, or --json.", file=sys.stderr)
        return 1

    lead = next((p for p in projections if p.compound_id == lead_id), None)
    lead_label = f"{lead.compound_id} ({lead.name})" if lead else "None"
    mrsd = f"{lead.phase1_mrsd_mg:.2f} mg" if lead and lead.phase1_mrsd_mg is not None else "UNMEASURED"
    print(
        f"IND dose projection complete: Lead candidate {lead_label} -> "
        f"Evidence: {lead.evidence_status if lead else 'NONE'}, Phase 1 MRSD: {mrsd} -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
