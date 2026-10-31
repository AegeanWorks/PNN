"""Phase 3 in vitro assay ingestion, hit qualification, and reporting for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .reports import write_csv, write_json, write_text


class AssayError(ValueError):
    """Raised when in vitro assay data is malformed or violates contract."""


@dataclass(frozen=True)
class DegradationResult:
    compound_id: str
    ec50_nm: float
    max_degradation_pct: float
    hill_slope: Optional[float] = None
    test_concentration_um: float = 10.0
    substrate: str = "DQ-CS"


@dataclass(frozen=True)
class PermeabilityResult:
    compound_id: str
    assay_type: str  # PAMPA or MDCK
    pe_1e6_cm_s: Optional[float] = None
    logbb: Optional[float] = None
    papp_ab_1e6_cm_s: Optional[float] = None
    efflux_ratio: Optional[float] = None
    pgp_substrate: Optional[bool] = None


@dataclass(frozen=True)
class CytotoxicityResult:
    compound_id: str
    cc50_um: float
    cell_line: str = "mHypoA-60/3"
    exposure_hours: int = 24


@dataclass(frozen=True)
class QualifiedHit:
    compound_id: str
    status: str  # CONFIRMED_HIT, QUALIFIED_LEAD, INCONCLUSIVE, CYTOTOXIC, INACTIVE
    ec50_nm: float
    max_degradation_pct: float
    pampa_pe: Optional[float]
    mdck_efflux: Optional[float]
    cc50_um: Optional[float]
    selectivity_index: Optional[float]
    bbb_status: str
    toxicity_status: str
    findings: str

    def as_row(self) -> Dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "status": self.status,
            "ec50_nm": f"{self.ec50_nm:.1f}",
            "max_degradation_pct": f"{self.max_degradation_pct:.1f}",
            "pampa_pe": f"{self.pampa_pe:.2f}" if self.pampa_pe is not None else "",
            "mdck_efflux": f"{self.mdck_efflux:.2f}" if self.mdck_efflux is not None else "",
            "cc50_um": f"{self.cc50_um:.1f}" if self.cc50_um is not None else "",
            "selectivity_index": f"{self.selectivity_index:.1f}" if self.selectivity_index is not None else "",
            "bbb_status": self.bbb_status,
            "toxicity_status": self.toxicity_status,
            "findings": self.findings,
        }


HIT_FIELDS = (
    "compound_id",
    "status",
    "ec50_nm",
    "max_degradation_pct",
    "pampa_pe",
    "mdck_efflux",
    "cc50_um",
    "selectivity_index",
    "bbb_status",
    "toxicity_status",
    "findings",
)


def _parse_float(val: str, field_name: str, row_id: str) -> Optional[float]:
    clean = val.strip()
    if not clean or clean.upper() in ("NA", "NAN", "NULL", "NONE"):
        return None
    try:
        f = float(clean)
    except ValueError:
        raise AssayError(f"Row '{row_id}': invalid numeric value '{val}' for {field_name}")
    if not math.isfinite(f):
        raise AssayError(f"Row '{row_id}': non-finite number '{val}' for {field_name}")
    return f


def load_degradation_data(path: Path) -> Dict[str, DegradationResult]:
    """Ingest PNN degradation assay CSV."""
    if not path.is_file():
        raise AssayError(f"Degradation assay file not found: {path}")

    results: Dict[str, DegradationResult] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise AssayError(f"Empty degradation assay CSV: {path}")

        req = {"compound_id", "ec50_nm", "max_degradation_pct"}
        normalized_cols = {col.strip().lower(): col for col in reader.fieldnames}
        missing = req - set(normalized_cols.keys())
        if missing:
            raise AssayError(f"Degradation CSV missing required column(s): {', '.join(missing)}")

        for row in reader:
            cid = row[normalized_cols["compound_id"]].strip()
            if not cid:
                continue

            ec50 = _parse_float(row[normalized_cols["ec50_nm"]], "ec50_nm", cid)
            if ec50 is None or ec50 <= 0:
                raise AssayError(f"Row '{cid}': ec50_nm must be a positive number")

            max_deg = _parse_float(row[normalized_cols["max_degradation_pct"]], "max_degradation_pct", cid)
            if max_deg is None or max_deg < 0:
                raise AssayError(f"Row '{cid}': max_degradation_pct must be >= 0")

            hill = None
            if "hill_slope" in normalized_cols:
                hill = _parse_float(row[normalized_cols["hill_slope"]], "hill_slope", cid)

            conc = 10.0
            if "test_concentration_um" in normalized_cols:
                parsed_c = _parse_float(row[normalized_cols["test_concentration_um"]], "test_concentration_um", cid)
                if parsed_c is not None and parsed_c > 0:
                    conc = parsed_c

            sub = "DQ-CS"
            if "substrate" in normalized_cols and row[normalized_cols["substrate"]].strip():
                sub = row[normalized_cols["substrate"]].strip()

            item = DegradationResult(
                compound_id=cid,
                ec50_nm=ec50,
                max_degradation_pct=max_deg,
                hill_slope=hill,
                test_concentration_um=conc,
                substrate=sub,
            )
            # If duplicate, keep the more potent (lower EC50)
            if cid not in results or results[cid].ec50_nm > item.ec50_nm:
                results[cid] = item

    if not results:
        raise AssayError(f"No valid degradation data rows found in {path}")
    return results


def load_permeability_data(path: Path) -> Dict[str, PermeabilityResult]:
    """Ingest BBB permeability assay CSV (PAMPA / MDCK-MDR1)."""
    if not path.is_file():
        raise AssayError(f"Permeability assay file not found: {path}")

    results: Dict[str, PermeabilityResult] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise AssayError(f"Empty permeability assay CSV: {path}")

        normalized_cols = {col.strip().lower(): col for col in reader.fieldnames}
        if "compound_id" not in normalized_cols:
            raise AssayError("Permeability CSV missing required 'compound_id' column")

        for row in reader:
            cid = row[normalized_cols["compound_id"]].strip()
            if not cid:
                continue

            assay_type = "PAMPA"
            if "assay_type" in normalized_cols and row[normalized_cols["assay_type"]].strip():
                assay_type = row[normalized_cols["assay_type"]].strip().upper()

            pe = None
            if "pe_1e6_cm_s" in normalized_cols:
                pe = _parse_float(row[normalized_cols["pe_1e6_cm_s"]], "pe_1e6_cm_s", cid)

            logbb = None
            if "logbb" in normalized_cols:
                logbb = _parse_float(row[normalized_cols["logbb"]], "logbb", cid)

            papp = None
            if "papp_ab_1e6_cm_s" in normalized_cols:
                papp = _parse_float(row[normalized_cols["papp_ab_1e6_cm_s"]], "papp_ab_1e6_cm_s", cid)

            efflux = None
            if "efflux_ratio" in normalized_cols:
                efflux = _parse_float(row[normalized_cols["efflux_ratio"]], "efflux_ratio", cid)

            pgp = None
            if "pgp_substrate" in normalized_cols:
                raw_pgp = row[normalized_cols["pgp_substrate"]].strip().lower()
                if raw_pgp in ("true", "yes", "1"):
                    pgp = True
                elif raw_pgp in ("false", "no", "0"):
                    pgp = False

            results[cid] = PermeabilityResult(
                compound_id=cid,
                assay_type=assay_type,
                pe_1e6_cm_s=pe,
                logbb=logbb,
                papp_ab_1e6_cm_s=papp,
                efflux_ratio=efflux,
                pgp_substrate=pgp,
            )
    return results


def load_cytotoxicity_data(path: Path) -> Dict[str, CytotoxicityResult]:
    """Ingest CellTiter-Glo cytotoxicity counter-screen CSV."""
    if not path.is_file():
        raise AssayError(f"Cytotoxicity assay file not found: {path}")

    results: Dict[str, CytotoxicityResult] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise AssayError(f"Empty cytotoxicity assay CSV: {path}")

        normalized_cols = {col.strip().lower(): col for col in reader.fieldnames}
        if "compound_id" not in normalized_cols or "cc50_um" not in normalized_cols:
            raise AssayError("Cytotoxicity CSV requires 'compound_id' and 'cc50_um' columns")

        for row in reader:
            cid = row[normalized_cols["compound_id"]].strip()
            if not cid:
                continue

            cc50 = _parse_float(row[normalized_cols["cc50_um"]], "cc50_um", cid)
            if cc50 is None or cc50 <= 0:
                raise AssayError(f"Row '{cid}': cc50_um must be a positive number")

            cell_line = "mHypoA-60/3"
            if "cell_line" in normalized_cols and row[normalized_cols["cell_line"]].strip():
                cell_line = row[normalized_cols["cell_line"]].strip()

            hours = 24
            if "exposure_hours" in normalized_cols:
                parsed_h = _parse_float(row[normalized_cols["exposure_hours"]], "exposure_hours", cid)
                if parsed_h is not None:
                    hours = int(parsed_h)

            results[cid] = CytotoxicityResult(
                compound_id=cid,
                cc50_um=cc50,
                cell_line=cell_line,
                exposure_hours=hours,
            )
    return results


def qualify_leads(
    degradation: Dict[str, DegradationResult],
    permeability: Optional[Dict[str, PermeabilityResult]] = None,
    cytotoxicity: Optional[Dict[str, CytotoxicityResult]] = None,
    ec50_ceiling_nm: float = 10000.0,
    max_deg_floor_pct: float = 20.0,
    pampa_pe_floor: float = 3.0,
    mdck_efflux_ceiling: float = 3.0,
    min_selectivity_index: float = 50.0,
) -> List[QualifiedHit]:
    """Qualify leads against Phase 3 in vitro degradation, permeability, and safety criteria."""
    perm_map = permeability or {}
    cyto_map = cytotoxicity or {}

    hits: List[QualifiedHit] = []
    for cid, deg in degradation.items():
        is_active = deg.ec50_nm <= ec50_ceiling_nm and deg.max_degradation_pct >= max_deg_floor_pct

        perm = perm_map.get(cid)
        pampa_pe = perm.pe_1e6_cm_s if perm else None
        mdck_efflux = perm.efflux_ratio if perm else None

        # Determine BBB Status
        bbb_pass = False
        bbb_status = "UNTESTED"
        if perm is not None:
            if pampa_pe is not None and pampa_pe >= pampa_pe_floor:
                bbb_pass = True
                bbb_status = "PERMEANT_PAMPA"
            elif mdck_efflux is not None and mdck_efflux < mdck_efflux_ceiling:
                bbb_pass = True
                bbb_status = "PERMEANT_MDCK"
            elif perm.logbb is not None and perm.logbb >= -1.0:
                bbb_pass = True
                bbb_status = "PERMEANT_LOGBB"
            elif mdck_efflux is not None and mdck_efflux >= mdck_efflux_ceiling:
                bbb_pass = False
                bbb_status = "HIGH_EFFLUX"
            else:
                bbb_pass = False
                bbb_status = "LOW_PERMEABILITY"

        # Determine Cytotoxicity and Selectivity
        cyto = cyto_map.get(cid)
        cc50 = cyto.cc50_um if cyto else None
        si = None
        tox_status = "UNTESTED"
        is_toxic = False

        if cc50 is not None:
            # SI = (CC50 * 1000) / EC50
            si = (cc50 * 1000.0) / deg.ec50_nm
            if cc50 < 20.0 or si < 10.0:
                is_toxic = True
                tox_status = "CYTOTOXIC"
            elif cc50 >= 50.0 or si >= min_selectivity_index:
                tox_status = "NON_TOXIC"
            else:
                tox_status = "MODERATE_SAFETY"

        # Determine Overall Hit Status
        findings: List[str] = []
        if is_toxic:
            status = "CYTOTOXIC"
            findings.append(f"Cell viability compromised (CC50={cc50:.1f} uM, SI={si:.1f})")
        elif not is_active:
            if deg.ec50_nm > 50000.0 or deg.max_degradation_pct < 10.0:
                status = "INACTIVE"
                findings.append(f"Insufficient PNN softening (EC50={deg.ec50_nm:.0f} nM, Max={deg.max_degradation_pct:.1f}%)")
            else:
                status = "INCONCLUSIVE"
                findings.append(f"Sub-threshold softening activity (EC50={deg.ec50_nm:.0f} nM)")
        else:
            findings.append(f"Confirmed PNN softening (EC50={deg.ec50_nm:.1f} nM, Max={deg.max_degradation_pct:.1f}%)")
            if bbb_pass and tox_status == "NON_TOXIC":
                status = "CONFIRMED_HIT"
                findings.append(f"Excellent permeability ({bbb_status}) and safety (SI={si:.1f})")
            elif bbb_pass:
                status = "QUALIFIED_LEAD"
                findings.append(f"Permeant ({bbb_status}); toxicity status: {tox_status}")
            elif bbb_status == "UNTESTED":
                status = "QUALIFIED_LEAD"
                findings.append("Softening confirmed; pending measured permeability")
            else:
                status = "INCONCLUSIVE"
                findings.append(f"Restricted BBB penetration ({bbb_status})")

        hits.append(
            QualifiedHit(
                compound_id=cid,
                status=status,
                ec50_nm=deg.ec50_nm,
                max_degradation_pct=deg.max_degradation_pct,
                pampa_pe=pampa_pe,
                mdck_efflux=mdck_efflux,
                cc50_um=cc50,
                selectivity_index=si,
                bbb_status=bbb_status,
                toxicity_status=tox_status,
                findings="; ".join(findings),
            )
        )

    # Sort: CONFIRMED_HIT first, then QUALIFIED_LEAD, then by EC50 ascending
    status_priority = {
        "CONFIRMED_HIT": 0,
        "QUALIFIED_LEAD": 1,
        "INCONCLUSIVE": 2,
        "INACTIVE": 3,
        "CYTOTOXIC": 4,
    }
    hits.sort(key=lambda h: (status_priority.get(h.status, 9), h.ec50_nm, -h.max_degradation_pct))
    return hits


def render_assay_report(hits: Sequence[QualifiedHit], criteria: Mapping[str, object]) -> str:
    """Render a comprehensive Phase 3 in vitro qualification Markdown report."""
    confirmed = [h for h in hits if h.status == "CONFIRMED_HIT"]
    qualified = [h for h in hits if h.status == "QUALIFIED_LEAD"]
    toxic = [h for h in hits if h.status == "CYTOTOXIC"]
    inactive = [h for h in hits if h.status == "INACTIVE"]

    lines = [
        "# PNN Phase 3 In Vitro Hit Qualification Report",
        "",
        "## Scientific Context & Evidence Invariants",
        "",
        "> [!IMPORTANT]",
        "> **Guideline G-001 (ADR-001)**: In vitro assay qualification confirms target degradation activity, "
        "membrane permeability, and safety windows in cell models. Never claim in vivo hypothalamic efficacy "
        "or therapeutic obesity reversal without pharmacokinetic brain-exposure and chronic DIO mouse study data.",
        "",
        "## Qualification Summary",
        "",
        f"- **Total Compounds Evaluated**: {len(hits)}",
        f"- **Confirmed Hits (Active + Permeant + Non-toxic)**: {len(confirmed)}",
        f"- **Qualified Leads (Active, pending complete safety/permeability)**: {len(qualified)}",
        f"- **Cytotoxic / Counter-Screen Failures**: {len(toxic)}",
        f"- **Inactive / Sub-threshold**: {len(inactive)}",
        "",
        "## Progression Thresholds",
        "",
        f"- **PNN Degradation**: EC50 <= {criteria.get('ec50_ceiling_nm', 10000)} nM and Max Cleavage >= {criteria.get('max_deg_floor_pct', 20)}%",
        f"- **BBB Permeability**: PAMPA Pe >= {criteria.get('pampa_pe_floor', 3.0)} x 10^-6 cm/s or MDCK Efflux Ratio < {criteria.get('mdck_efflux_ceiling', 3.0)}",
        f"- **Safety Window**: CC50 >= 50 uM or Selectivity Index (CC50/EC50) >= {criteria.get('min_selectivity_index', 50)}",
        "",
        "## Confirmed In Vitro Hits",
        "",
        "| Compound ID | EC50 (nM) | Max Deg (%) | PAMPA Pe | MDCK Efflux | CC50 (uM) | Selectivity Index | Status | Findings |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for h in confirmed + qualified:
        pampa_str = f"{h.pampa_pe:.2f}" if h.pampa_pe is not None else "-"
        mdck_str = f"{h.mdck_efflux:.2f}" if h.mdck_efflux is not None else "-"
        cc50_str = f"{h.cc50_um:.1f}" if h.cc50_um is not None else "-"
        si_str = f"{h.selectivity_index:.1f}" if h.selectivity_index is not None else "-"
        lines.append(
            f"| `{h.compound_id}` | {h.ec50_nm:.1f} | {h.max_degradation_pct:.1f}% | "
            f"{pampa_str} | {mdck_str} | {cc50_str} | {si_str} | **{h.status}** | {h.findings} |"
        )

    lines.extend([
        "",
        "## Recommended Next Steps (Phase 4 In Vivo Proof-of-Concept)",
        "",
        "1. **Primary ARC Neuron Target Engagement**: Evaluate top confirmed hits for pIR-Y1161 and pAKT-S473 restoration in primary mouse hypothalamic ARC neurons.",
        "2. **Single-Dose Pharmacokinetics**: Determine oral bioavailability (F), brain/plasma AUC ratio, and clearance in C57BL/6J mice.",
        "3. **Chronic DIO Proof-of-Concept**: Administer oral active leads for 14-28 days in diet-induced obese mice; measure body weight, food intake, and WFA-positive net degradation in the arcuate nucleus.",
        "",
    ])
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.assay",
        description="Ingest Phase 3 in vitro assay data, qualify candidate leads, and generate qualification reports.",
    )
    parser.add_argument("--degradation", required=True, help="path to PNN degradation assay CSV (EC50, max pct)")
    parser.add_argument("--permeability", help="optional path to BBB permeability assay CSV (PAMPA/MDCK)")
    parser.add_argument("--cytotoxicity", help="optional path to cytotoxicity counter-screen CSV (CC50)")
    parser.add_argument("--out", default="qualified_hits.csv", help="output destination for qualified hits CSV")
    parser.add_argument("--json", dest="json_path", help="optional output path for JSON summary")
    parser.add_argument("--report", dest="report_path", help="optional output path for Markdown qualification report")
    parser.add_argument("--ec50-ceiling", type=float, default=10000.0, help="maximum acceptable EC50 in nM (default: 10000)")
    parser.add_argument("--max-deg-floor", type=float, default=20.0, help="minimum acceptable degradation pct (default: 20)")
    parser.add_argument("--min-si", type=float, default=50.0, help="minimum selectivity index (CC50/EC50) (default: 50)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        deg = load_degradation_data(Path(args.degradation))
        perm = load_permeability_data(Path(args.permeability)) if args.permeability else None
        cyto = load_cytotoxicity_data(Path(args.cytotoxicity)) if args.cytotoxicity else None

        criteria = {
            "ec50_ceiling_nm": args.ec50_ceiling,
            "max_deg_floor_pct": args.max_deg_floor,
            "min_selectivity_index": args.min_si,
            "pampa_pe_floor": 3.0,
            "mdck_efflux_ceiling": 3.0,
        }
        hits = qualify_leads(
            degradation=deg,
            permeability=perm,
            cytotoxicity=cyto,
            ec50_ceiling_nm=args.ec50_ceiling,
            max_deg_floor_pct=args.max_deg_floor,
            min_selectivity_index=args.min_si,
        )

        out_path = Path(args.out)
        write_csv(out_path, HIT_FIELDS, [h.as_row() for h in hits])

        if args.json_path:
            summary = {
                "total_evaluated": len(hits),
                "by_status": {
                    status: sum(1 for h in hits if h.status == status)
                    for status in ("CONFIRMED_HIT", "QUALIFIED_LEAD", "INCONCLUSIVE", "CYTOTOXIC", "INACTIVE")
                },
                "criteria": criteria,
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), summary)

        if args.report_path:
            report_text = render_assay_report(hits, criteria)
            write_text(Path(args.report_path), report_text)

    except AssayError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: in vitro assay qualification requires valid CSV data with numeric assay parameters.", file=sys.stderr)
        print("Fix: verify columns against templates in data/assays/ and check numeric values.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: destination file could not be written to disk.", file=sys.stderr)
        print("Fix: specify a writable --out, --json, or --report file path.", file=sys.stderr)
        return 1

    confirmed_count = sum(1 for h in hits if h.status == "CONFIRMED_HIT")
    qualified_count = sum(1 for h in hits if h.status == "QUALIFIED_LEAD")
    print(
        f"In vitro hit qualification complete: {confirmed_count} confirmed hit(s), "
        f"{qualified_count} qualified lead(s) from {len(hits)} evaluated compound(s) -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
