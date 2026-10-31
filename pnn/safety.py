"""IND-enabling safety pharmacology and repeat-dose toxicology assessment for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from .reports import write_csv, write_json, write_text

KM_BY_SPECIES = {
    "mouse": 3.0,
    "rat": 6.0,
    "hamster": 5.0,
    "guinea pig": 7.0,
    "rabbit": 12.0,
    "dog": 20.0,
    "monkey": 12.0,
    "marmoset": 12.0,
    "human": 37.0,
}
HUMAN_KM = 37.0
HUMAN_WEIGHT_KG = 60.0

SAFETY_FIELDS = (
    "compound_id",
    "status",
    "most_sensitive_species",
    "noael_mg_kg",
    "noael_hed_mg_kg",
    "phase1_mrsd_mg",
    "herg_ic50_um",
    "free_cmax_um",
    "herg_margin",
    "cyp3a4_inhibition_pct",
    "genotoxicity",
    "target_organs",
    "findings",
)


class SafetyError(ValueError):
    """Raised when safety pharmacology or toxicology data is malformed or violates contract."""


@dataclass(frozen=True)
class ToxStudy:
    compound_id: str
    species: str
    study_type: str
    duration_days: int
    noael_mg_kg: float
    loael_mg_kg: Optional[float]
    target_organ: str
    findings: str


@dataclass(frozen=True)
class SafetyPharm:
    compound_id: str
    herg_ic50_um: Optional[float] = None
    free_cmax_um: Optional[float] = None
    cyp3a4_inhibition_pct: Optional[float] = None
    cyp2d6_inhibition_pct: Optional[float] = None
    ames_positive: Optional[bool] = None
    micronucleus_positive: Optional[bool] = None
    respiratory_noael_mg_kg: Optional[float] = None
    cns_noael_mg_kg: Optional[float] = None


@dataclass(frozen=True)
class SafetyAssessment:
    compound_id: str
    status: str
    most_sensitive_species: str
    noael_mg_kg: Optional[float]
    noael_hed_mg_kg: Optional[float]
    phase1_mrsd_mg: Optional[float]
    herg_ic50_um: Optional[float]
    free_cmax_um: Optional[float]
    herg_margin: Optional[float]
    cyp3a4_inhibition_pct: Optional[float]
    genotoxicity: str
    target_organs: str
    findings: str

    def as_row(self) -> Dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "status": self.status,
            "most_sensitive_species": self.most_sensitive_species,
            "noael_mg_kg": f"{self.noael_mg_kg:.1f}" if self.noael_mg_kg is not None else "",
            "noael_hed_mg_kg": f"{self.noael_hed_mg_kg:.3f}" if self.noael_hed_mg_kg is not None else "",
            "phase1_mrsd_mg": f"{self.phase1_mrsd_mg:.2f}" if self.phase1_mrsd_mg is not None else "",
            "herg_ic50_um": f"{self.herg_ic50_um:.2f}" if self.herg_ic50_um is not None else "",
            "free_cmax_um": f"{self.free_cmax_um:.3f}" if self.free_cmax_um is not None else "",
            "herg_margin": f"{self.herg_margin:.1f}" if self.herg_margin is not None else "",
            "cyp3a4_inhibition_pct": f"{self.cyp3a4_inhibition_pct:.1f}" if self.cyp3a4_inhibition_pct is not None else "",
            "genotoxicity": self.genotoxicity,
            "target_organs": self.target_organs,
            "findings": self.findings,
        }


def _normalize(name: str) -> str:
    return name.strip().lower()


def _species_key(species: str) -> str:
    return " ".join(species.strip().lower().replace("_", " ").split())


def _parse_required_float(value: str, field_name: str, row_id: str) -> float:
    clean = value.strip()
    if not clean or clean.upper() in ("NA", "NAN", "NULL", "NONE", "-"):
        raise SafetyError(f"Row '{row_id}': {field_name} is required")
    try:
        number = float(clean)
    except ValueError:
        raise SafetyError(f"Row '{row_id}': invalid numeric value '{value}' for {field_name}")
    if not math.isfinite(number):
        raise SafetyError(f"Row '{row_id}': non-finite number '{value}' for {field_name}")
    return number


def _parse_optional_float(value: Optional[str], field_name: str, row_id: str) -> Optional[float]:
    if value is None:
        return None
    clean = value.strip()
    if not clean or clean.upper() in ("NA", "NAN", "NULL", "NONE", "-"):
        return None
    try:
        number = float(clean)
    except ValueError:
        raise SafetyError(f"Row '{row_id}': invalid numeric value '{value}' for {field_name}")
    if not math.isfinite(number):
        raise SafetyError(f"Row '{row_id}': non-finite number '{value}' for {field_name}")
    return number


def _parse_optional_bool(value: Optional[str], field_name: str, row_id: str) -> Optional[bool]:
    if value is None:
        return None
    clean = value.strip().lower()
    if not clean or clean in ("na", "nan", "null", "none", "-"):
        return None
    if clean in ("true", "yes", "y", "1", "positive", "pos"):
        return True
    if clean in ("false", "no", "n", "0", "negative", "neg"):
        return False
    raise SafetyError(f"Row '{row_id}': invalid boolean value '{value}' for {field_name}")


def _check_percent(value: Optional[float], field_name: str, row_id: str) -> Optional[float]:
    if value is None:
        return None
    if value < 0.0 or value > 100.0:
        raise SafetyError(f"Row '{row_id}': {field_name} must be within 0-100 percent")
    return value


def load_tox_studies(path: Path) -> List[ToxStudy]:
    """Ingest repeat-dose toxicology studies with NOAEL/LOAEL and target-organ findings."""
    if not path.is_file():
        raise SafetyError(f"Repeat-dose toxicology file not found: {path}")

    studies: List[ToxStudy] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SafetyError(f"Empty toxicology CSV: {path}")
        columns = {_normalize(name): name for name in reader.fieldnames}
        missing = {"compound_id", "species", "noael_mg_kg"} - set(columns)
        if missing:
            raise SafetyError(f"Toxicology CSV missing required column(s): {', '.join(sorted(missing))}")

        for row in reader:
            cid = row[columns["compound_id"]].strip()
            if not cid:
                continue
            species = _species_key(row[columns["species"]])
            if species not in KM_BY_SPECIES:
                raise SafetyError(
                    f"Row '{cid}': unsupported species '{row[columns['species']].strip()}'; "
                    f"supported: {', '.join(sorted(KM_BY_SPECIES))}"
                )

            noael = _parse_required_float(row[columns["noael_mg_kg"]], "noael_mg_kg", cid)
            if noael <= 0.0:
                raise SafetyError(f"Row '{cid}': noael_mg_kg must be > 0")

            loael = None
            if "loael_mg_kg" in columns:
                loael = _parse_optional_float(row[columns["loael_mg_kg"]], "loael_mg_kg", cid)
                if loael is not None and loael < noael:
                    raise SafetyError(f"Row '{cid}': loael_mg_kg must be >= noael_mg_kg")

            duration = 28
            if "duration_days" in columns:
                raw_duration = _parse_optional_float(row[columns["duration_days"]], "duration_days", cid)
                if raw_duration is not None:
                    if raw_duration <= 0.0 or raw_duration != int(raw_duration):
                        raise SafetyError(f"Row '{cid}': duration_days must be a positive integer")
                    duration = int(raw_duration)

            study_type = "REPEAT_DOSE"
            if "study_type" in columns and row[columns["study_type"]].strip():
                study_type = row[columns["study_type"]].strip()

            target_organ = ""
            if "target_organ" in columns:
                target_organ = row[columns["target_organ"]].strip()

            findings = ""
            if "findings" in columns:
                findings = row[columns["findings"]].strip()

            studies.append(
                ToxStudy(
                    compound_id=cid,
                    species=species,
                    study_type=study_type,
                    duration_days=duration,
                    noael_mg_kg=noael,
                    loael_mg_kg=loael,
                    target_organ=target_organ,
                    findings=findings,
                )
            )

    if not studies:
        raise SafetyError(f"No valid repeat-dose toxicology rows found in {path}")
    return studies


def load_safety_pharmacology(path: Path) -> Dict[str, SafetyPharm]:
    """Ingest ICH S7A/S7B safety pharmacology and genotoxicity panel results."""
    if not path.is_file():
        raise SafetyError(f"Safety pharmacology file not found: {path}")

    records: Dict[str, SafetyPharm] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SafetyError(f"Empty safety pharmacology CSV: {path}")
        columns = {_normalize(name): name for name in reader.fieldnames}
        if "compound_id" not in columns:
            raise SafetyError("Safety pharmacology CSV missing required column: compound_id")

        for row in reader:
            cid = row[columns["compound_id"]].strip()
            if not cid:
                continue
            if cid in records:
                raise SafetyError(f"Duplicate safety pharmacology row for compound '{cid}'")

            def optional(field: str) -> Optional[str]:
                return row[columns[field]] if field in columns else None

            herg = _parse_optional_float(optional("herg_ic50_um"), "herg_ic50_um", cid)
            if herg is not None and herg <= 0.0:
                raise SafetyError(f"Row '{cid}': herg_ic50_um must be > 0")
            free_cmax = _parse_optional_float(optional("free_cmax_um"), "free_cmax_um", cid)
            if free_cmax is not None and free_cmax <= 0.0:
                raise SafetyError(f"Row '{cid}': free_cmax_um must be > 0")

            cyp3a4 = _check_percent(
                _parse_optional_float(optional("cyp3a4_inhibition_pct"), "cyp3a4_inhibition_pct", cid),
                "cyp3a4_inhibition_pct",
                cid,
            )
            cyp2d6 = _check_percent(
                _parse_optional_float(optional("cyp2d6_inhibition_pct"), "cyp2d6_inhibition_pct", cid),
                "cyp2d6_inhibition_pct",
                cid,
            )

            respiratory = _parse_optional_float(optional("respiratory_noael_mg_kg"), "respiratory_noael_mg_kg", cid)
            if respiratory is not None and respiratory < 0.0:
                raise SafetyError(f"Row '{cid}': respiratory_noael_mg_kg must be >= 0")
            cns = _parse_optional_float(optional("cns_noael_mg_kg"), "cns_noael_mg_kg", cid)
            if cns is not None and cns < 0.0:
                raise SafetyError(f"Row '{cid}': cns_noael_mg_kg must be >= 0")

            records[cid] = SafetyPharm(
                compound_id=cid,
                herg_ic50_um=herg,
                free_cmax_um=free_cmax,
                cyp3a4_inhibition_pct=cyp3a4,
                cyp2d6_inhibition_pct=cyp2d6,
                ames_positive=_parse_optional_bool(optional("ames_positive"), "ames_positive", cid),
                micronucleus_positive=_parse_optional_bool(
                    optional("micronucleus_positive"), "micronucleus_positive", cid
                ),
                respiratory_noael_mg_kg=respiratory,
                cns_noael_mg_kg=cns,
            )

    if not records:
        raise SafetyError(f"No valid safety pharmacology rows found in {path}")
    return records


def noael_hed_mg_kg(noael_mg_kg: float, species: str, human_km: float = HUMAN_KM) -> float:
    """Convert a species NOAEL to a human equivalent dose using its body-surface-area factor."""
    key = _species_key(species)
    if key not in KM_BY_SPECIES:
        raise SafetyError(f"Unsupported species '{species}'; supported: {', '.join(sorted(KM_BY_SPECIES))}")
    return noael_mg_kg * (KM_BY_SPECIES[key] / human_km)


def most_sensitive_study(studies: Sequence[ToxStudy]) -> Optional[ToxStudy]:
    """
    Return the study with the lowest NOAEL human equivalent dose.

    Sensitivity is ranked on HED, never on raw mg/kg: a 30 mg/kg dog NOAEL (HED 16.2) is less
    sensitive than a 50 mg/kg rat NOAEL (HED 8.1) because the K_m factors differ by more than 3x.
    Ties resolve by species name, then study type, for deterministic output.
    """
    if not studies:
        return None
    return min(
        studies,
        key=lambda s: (
            noael_hed_mg_kg(s.noael_mg_kg, s.species),
            s.species,
            s.study_type,
        ),
    )


def _genotoxicity(safety_pharm: Optional[SafetyPharm]) -> str:
    if safety_pharm is None:
        return "UNKNOWN"
    results = (safety_pharm.ames_positive, safety_pharm.micronucleus_positive)
    if any(result is True for result in results):
        return "POSITIVE"
    if all(result is False for result in results):
        return "NEGATIVE"
    return "UNKNOWN"


def assess_safety(
    tox_studies: Sequence[ToxStudy],
    safety_pharm: Optional[Mapping[str, SafetyPharm]] = None,
    safety_factor: float = 10.0,
    herg_margin_min: float = 30.0,
    herg_margin_fail: float = 10.0,
    cyp_inhibition_max_pct: float = 50.0,
    cyp_inhibition_caution_pct: float = 30.0,
    human_weight_kg: float = HUMAN_WEIGHT_KG,
) -> List[SafetyAssessment]:
    """Assess each compound against IND-enabling toxicology and safety pharmacology criteria."""
    if safety_factor <= 0.0:
        raise SafetyError("safety_factor must be > 0")
    if not 0.0 < herg_margin_fail <= herg_margin_min:
        raise SafetyError("herg_margin_fail must be > 0 and <= herg_margin_min")
    if not 0.0 <= cyp_inhibition_caution_pct < cyp_inhibition_max_pct <= 100.0:
        raise SafetyError("CYP thresholds must satisfy 0 <= caution < max <= 100")

    pharm_map = safety_pharm or {}
    by_compound: Dict[str, List[ToxStudy]] = {}
    for study in tox_studies:
        by_compound.setdefault(study.compound_id, []).append(study)

    assessments: List[SafetyAssessment] = []
    for cid in sorted(set(by_compound) | set(pharm_map)):
        studies = sorted(by_compound.get(cid, []), key=lambda s: (s.noael_mg_kg, s.species, s.study_type))
        sensitive_study = most_sensitive_study(studies)
        pharm = pharm_map.get(cid)
        failures: List[str] = []
        cautions: List[str] = []
        evidence: List[str] = []

        noael = sensitive_study.noael_mg_kg if sensitive_study else None
        sensitive_species = sensitive_study.species if sensitive_study else ""
        noael_hed = None
        mrsd_mg = None
        if sensitive_study is not None:
            noael_hed = noael_hed_mg_kg(sensitive_study.noael_mg_kg, sensitive_study.species)
            mrsd_mg = (noael_hed / safety_factor) * human_weight_kg
            evidence.append(
                f"Most sensitive NOAEL {sensitive_study.noael_mg_kg:.1f} mg/kg in {sensitive_species} "
                f"(HED {noael_hed:.3f} mg/kg)"
            )
        else:
            cautions.append("No repeat-dose toxicology study submitted")

        target_organs = sorted(
            {
                study.target_organ
                for study in studies
                if study.target_organ and study.target_organ.lower() not in ("none", "na", "n/a", "-")
            }
        )
        if target_organs:
            evidence.append(f"Target organs: {', '.join(target_organs)}")

        genotoxicity = _genotoxicity(pharm)
        if genotoxicity == "POSITIVE":
            failures.append("Positive genotoxicity finding (Ames or micronucleus)")
        elif genotoxicity == "UNKNOWN":
            cautions.append("Genotoxicity data incomplete")
        else:
            evidence.append("Ames and micronucleus negative")

        herg_ic50 = pharm.herg_ic50_um if pharm else None
        free_cmax = pharm.free_cmax_um if pharm else None
        herg_margin = None
        if herg_ic50 is not None and free_cmax is not None:
            herg_margin = herg_ic50 / free_cmax
            if herg_margin < herg_margin_fail:
                failures.append(f"Narrow hERG margin ({herg_margin:.1f}-fold < {herg_margin_fail:.0f}-fold)")
            elif herg_margin < herg_margin_min:
                cautions.append(f"hERG margin below target ({herg_margin:.1f}-fold < {herg_margin_min:.0f}-fold)")
            else:
                evidence.append(f"hERG margin {herg_margin:.1f}-fold")
        else:
            cautions.append("hERG margin not computable (missing IC50 or free Cmax)")

        cyp3a4 = pharm.cyp3a4_inhibition_pct if pharm else None
        if cyp3a4 is None:
            cautions.append("CYP3A4 inhibition not measured")
        elif cyp3a4 >= cyp_inhibition_max_pct:
            failures.append(f"CYP3A4 inhibition {cyp3a4:.1f}% >= {cyp_inhibition_max_pct:.0f}%")
        elif cyp3a4 >= cyp_inhibition_caution_pct:
            cautions.append(f"CYP3A4 inhibition {cyp3a4:.1f}% approaching the {cyp_inhibition_max_pct:.0f}% ceiling")
        else:
            evidence.append(f"CYP3A4 inhibition {cyp3a4:.1f}%")

        if pharm is None:
            cautions.append("Safety pharmacology panel not submitted")
        else:
            missing_safety_pharm = [
                label
                for label, value in (("respiratory", pharm.respiratory_noael_mg_kg), ("CNS", pharm.cns_noael_mg_kg))
                if value is None
            ]
            if missing_safety_pharm:
                cautions.append(f"Missing {'/'.join(missing_safety_pharm)} safety NOAEL")
            else:
                evidence.append("CNS and respiratory safety NOAELs recorded")

        if failures:
            status = "FAIL_SAFETY"
        elif cautions:
            status = "CAUTION"
        else:
            status = "SAFE_TO_PROGRESS"

        assessments.append(
            SafetyAssessment(
                compound_id=cid,
                status=status,
                most_sensitive_species=sensitive_species,
                noael_mg_kg=noael,
                noael_hed_mg_kg=noael_hed,
                phase1_mrsd_mg=mrsd_mg,
                herg_ic50_um=herg_ic50,
                free_cmax_um=free_cmax,
                herg_margin=herg_margin,
                cyp3a4_inhibition_pct=cyp3a4,
                genotoxicity=genotoxicity,
                target_organs="; ".join(target_organs),
                findings="; ".join(failures + cautions + evidence),
            )
        )

    priority = {"SAFE_TO_PROGRESS": 0, "CAUTION": 1, "FAIL_SAFETY": 2}
    assessments.sort(key=lambda a: (priority[a.status], a.compound_id))
    return assessments


def render_safety_report(assessments: Sequence[SafetyAssessment], criteria: Mapping[str, object]) -> str:
    """Render a Markdown IND-enabling safety package summary."""
    safe = [a for a in assessments if a.status == "SAFE_TO_PROGRESS"]
    caution = [a for a in assessments if a.status == "CAUTION"]
    failed = [a for a in assessments if a.status == "FAIL_SAFETY"]

    lines = [
        "# PNN IND-Enabling Safety Pharmacology & Toxicology Package",
        "",
        "> [!IMPORTANT]",
        "> **Evidence limit (G-001)**: This package summarizes submitted study results and applies "
        "deterministic progression thresholds. It is not a GLP audit, does not replace histopathology or "
        "pathologist review, and never establishes that a compound is safe for human administration. "
        "Missing values stay UNKNOWN and are never treated as passes.",
        "",
        "## Assessment Summary",
        "",
        f"- **Compounds assessed**: {len(assessments)}",
        f"- **Safe to progress**: {len(safe)}",
        f"- **Caution / incomplete evidence**: {len(caution)}",
        f"- **Failed safety thresholds**: {len(failed)}",
        "",
        "## Progression Thresholds",
        "",
        f"- **NOAEL to HED**: species K_m / human K_m = {criteria.get('human_km', HUMAN_KM)}; "
        f"safety factor {criteria.get('safety_factor', 10.0)}x",
        f"- **hERG margin**: >= {criteria.get('herg_margin_min', 30.0)}-fold; fail below "
        f"{criteria.get('herg_margin_fail', 10.0)}-fold",
        f"- **CYP3A4 inhibition**: caution at {criteria.get('cyp_caution_pct', 30.0)}%, "
        f"fail at {criteria.get('cyp_max_pct', 50.0)}%",
        "- **Genetic toxicology**: Ames and micronucleus must both be negative; any positive result fails",
        "",
        "## Candidate Safety Matrix",
        "",
        "| Compound ID | Status | Sensitive Species | NOAEL (mg/kg) | NOAEL HED (mg/kg) | Phase 1 MRSD (mg) | hERG Margin | CYP3A4 (%) | Genotoxicity | Target Organs |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]

    for a in assessments:
        noael = f"{a.noael_mg_kg:.1f}" if a.noael_mg_kg is not None else "-"
        hed = f"{a.noael_hed_mg_kg:.3f}" if a.noael_hed_mg_kg is not None else "-"
        mrsd = f"{a.phase1_mrsd_mg:.2f}" if a.phase1_mrsd_mg is not None else "-"
        margin = f"{a.herg_margin:.1f}x" if a.herg_margin is not None else "-"
        cyp = f"{a.cyp3a4_inhibition_pct:.1f}" if a.cyp3a4_inhibition_pct is not None else "-"
        organs = a.target_organs or "-"
        lines.append(
            f"| `{a.compound_id}` | **{a.status}** | {a.most_sensitive_species or '-'} | {noael} | {hed} | "
            f"{mrsd} | {margin} | {cyp} | {a.genotoxicity} | {organs} |"
        )

    lines.extend(["", "## Findings", ""])
    for a in assessments:
        lines.append(f"- `{a.compound_id}`: {a.findings}")

    blockers = [a for a in assessments if a.status != "SAFE_TO_PROGRESS"]
    lines.extend(["", "## Evidence Gaps & Required Follow-Up", ""])
    if blockers:
        for a in blockers:
            lines.append(
                f"- `{a.compound_id}` ({a.status}): resolve the flagged findings with a complete GLP "
                "repeat-dose study, ICH S7A/S7B battery and definitive genotoxicity package before IND filing."
            )
    else:
        lines.append(
            "- No gaps detected in the submitted datasets. Confirm the underlying studies are GLP-compliant, "
            "pathologist-reviewed and reflected in the final regulatory submission."
        )

    lines.extend(
        [
            "",
            "## Package Limits",
            "",
            "1. Study quality, GLP compliance and pathology peer review are outside this deterministic check.",
            "2. The MRSD here is a NOAEL-based projection; reconcile it with the allometric projection in "
            "`pnn.ind` before selecting a first-in-human starting dose.",
            "3. Species coverage, exposure margins and toxicokinetic data must be reviewed by a qualified "
            "toxicologist before any regulatory submission.",
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.safety",
        description="Ingest repeat-dose toxicology and safety pharmacology data and assemble an IND-enabling safety package.",
    )
    parser.add_argument("--tox", help="path to repeat-dose toxicology CSV (NOAEL/LOAEL studies)")
    parser.add_argument("--pharm", help="path to safety pharmacology and genotoxicity panel CSV")
    parser.add_argument("--out", default="safety_assessments.csv", help="output destination CSV")
    parser.add_argument("--json", dest="json_path", help="optional output path for JSON summary")
    parser.add_argument("--report", dest="report_path", help="optional output path for Markdown safety package")
    parser.add_argument("--safety-factor", type=float, default=10.0, help="NOAEL-to-MRSD safety divisor (default: 10)")
    parser.add_argument("--herg-margin", type=float, default=30.0, help="target hERG margin in fold (default: 30)")
    parser.add_argument("--cyp-max", type=float, default=50.0, help="CYP3A4 inhibition failure ceiling pct (default: 50)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if not args.tox and not args.pharm:
            raise SafetyError("at least one of --tox or --pharm is required")

        tox_studies = load_tox_studies(Path(args.tox)) if args.tox else []
        safety_pharm = load_safety_pharmacology(Path(args.pharm)) if args.pharm else {}

        criteria = {
            "safety_factor": args.safety_factor,
            "herg_margin_min": args.herg_margin,
            "herg_margin_fail": 10.0,
            "cyp_max_pct": args.cyp_max,
            "cyp_caution_pct": 30.0,
            "human_km": HUMAN_KM,
        }
        assessments = assess_safety(
            tox_studies=tox_studies,
            safety_pharm=safety_pharm,
            safety_factor=args.safety_factor,
            herg_margin_min=args.herg_margin,
            cyp_inhibition_max_pct=args.cyp_max,
        )

        out_path = Path(args.out)
        write_csv(out_path, SAFETY_FIELDS, [a.as_row() for a in assessments])

        if args.report_path:
            write_text(Path(args.report_path), render_safety_report(assessments, criteria))

        if args.json_path:
            summary = {
                "total_assessed": len(assessments),
                "by_status": {
                    status: sum(1 for a in assessments if a.status == status)
                    for status in ("SAFE_TO_PROGRESS", "CAUTION", "FAIL_SAFETY")
                },
                "criteria": criteria,
                "tox_study_count": len(tox_studies),
                "safety_pharmacology_count": len(safety_pharm),
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), summary)

    except SafetyError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: safety package assembly requires valid toxicology and pharmacology CSV datasets.", file=sys.stderr)
        print("Fix: verify columns against templates in data/toxicology/ and provide positive, finite measurements.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: safety package output files could not be written to disk.", file=sys.stderr)
        print("Fix: provide writable --out, --json, or --report directory paths.", file=sys.stderr)
        return 1

    safe_count = sum(1 for a in assessments if a.status == "SAFE_TO_PROGRESS")
    failed_count = sum(1 for a in assessments if a.status == "FAIL_SAFETY")
    print(
        f"Safety assessment complete: {safe_count} safe to progress, {failed_count} failed, "
        f"{len(assessments)} compound(s) assessed -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
