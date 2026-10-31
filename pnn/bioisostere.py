"""Bioisosteric replacement and metabolic soft-spot optimization engine for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .mpo import calculate_cns_mpo
from .reports import write_csv, write_json, write_text


class BioisostereError(ValueError):
    """Raised when compound SMILES or bioisosteric transformation rules are invalid."""


@dataclass(frozen=True)
class BioisostereRule:
    rule_id: str
    name: str
    target_motif: str
    replacement_motif: str
    rationale: str
    delta_mw: float
    delta_logp: float
    delta_tpsa: float
    delta_hbd: int
    delta_hba: int
    delta_rotb: int


@dataclass(frozen=True)
class GeneratedAnalog:
    analog_id: str
    parent_id: str
    rule_id: str
    transformation_name: str
    parent_smiles: str
    analog_smiles: str
    mw: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotb: int
    cns_mpo_score: float
    delta_cns_mpo: float
    transformation_rationale: str
    optimization_verdict: str

    def as_row(self) -> Dict[str, object]:
        return {
            "analog_id": self.analog_id,
            "parent_id": self.parent_id,
            "rule_id": self.rule_id,
            "transformation_name": self.transformation_name,
            "parent_smiles": self.parent_smiles,
            "analog_smiles": self.analog_smiles,
            "mw": f"{self.mw:.1f}",
            "logp": f"{self.logp:.2f}",
            "tpsa": f"{self.tpsa:.1f}",
            "hbd": self.hbd,
            "hba": self.hba,
            "rotb": self.rotb,
            "cns_mpo_score": f"{self.cns_mpo_score:.2f}",
            "delta_cns_mpo": f"{self.delta_cns_mpo:+.2f}",
            "transformation_rationale": self.transformation_rationale,
            "optimization_verdict": self.optimization_verdict,
        }


ANALOG_FIELDS = (
    "analog_id",
    "parent_id",
    "rule_id",
    "transformation_name",
    "parent_smiles",
    "analog_smiles",
    "mw",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "rotb",
    "cns_mpo_score",
    "delta_cns_mpo",
    "transformation_rationale",
    "optimization_verdict",
)


BIOISOSTERE_RULES: Sequence[BioisostereRule] = (
    BioisostereRule(
        rule_id="BIO-01",
        name="Aromatic Fluorination (CYP Soft-Spot Block)",
        target_motif="c1ccccc1",
        replacement_motif="c1ccc(F)cc1",
        rationale="Fluorine blocks CYP3A4/2D6 aromatic para-hydroxylation, suppressing hepatic metabolic clearance.",
        delta_mw=18.0,
        delta_logp=0.14,
        delta_tpsa=0.0,
        delta_hbd=0,
        delta_hba=0,
        delta_rotb=0,
    ),
    BioisostereRule(
        rule_id="BIO-02",
        name="Trifluoromethyl Substitution (Lipophilic Shield)",
        target_motif="c1ccccc1",
        replacement_motif="c1ccc(C(F)(F)F)cc1",
        rationale="CF3 substitution confers metabolic inertness and enhances blood-brain barrier partition.",
        delta_mw=68.0,
        delta_logp=0.90,
        delta_tpsa=0.0,
        delta_hbd=0,
        delta_hba=0,
        delta_rotb=0,
    ),
    BioisostereRule(
        rule_id="BIO-03",
        name="Cyclopropyl Amide Rigidification",
        target_motif="NC(=O)C",
        replacement_motif="NC(=O)C1CC1",
        rationale="Cyclopropyl ring restricts torsional flexibility of the methyl group, improving conformational stability.",
        delta_mw=40.0,
        delta_logp=0.60,
        delta_tpsa=0.0,
        delta_hbd=0,
        delta_hba=0,
        delta_rotb=-1,
    ),
    BioisostereRule(
        rule_id="BIO-04",
        name="Piperazine-to-Morpholine Basicity Attenuation",
        target_motif="N1CCNCC1",
        replacement_motif="N1CCOCC1",
        rationale="Morpholine oxygen decreases the distal basic center pKa by ~1.5 units, mitigating lysosomal trapping and hERG liability.",
        delta_mw=1.0,
        delta_logp=-0.35,
        delta_tpsa=-2.8,
        delta_hbd=-1,
        delta_hba=0,
        delta_rotb=0,
    ),
    BioisostereRule(
        rule_id="BIO-05",
        name="Deuterium Switch (Benzylic Metabolic Shield)",
        target_motif="CN1CCNCC1",
        replacement_motif="C([2H])([2H])N1CCNCC1",
        rationale="Deuteration of the benzylic methylene exploits the primary kinetic isotope effect to attenuate N-dealkylation.",
        delta_mw=2.0,
        delta_logp=0.0,
        delta_tpsa=0.0,
        delta_hbd=0,
        delta_hba=0,
        delta_rotb=0,
    ),
    BioisostereRule(
        rule_id="BIO-06",
        name="Ester-to-Oxadiazole Replacement",
        target_motif="C(=O)OC",
        replacement_motif="c1nc(C)no1",
        rationale="1,2,4-Oxadiazole heterocycle replaces labile ester linkers, conferring resistance to plasma esterases.",
        delta_mw=24.0,
        delta_logp=-0.10,
        delta_tpsa=12.5,
        delta_hbd=0,
        delta_hba=1,
        delta_rotb=0,
    ),
)


def _ring_label_occurrences(smiles: str) -> Dict[str, int]:
    """Count single-digit and %nn ring labels in a SMILES string (label -> occurrences)."""
    counts: Dict[str, int] = {}
    i = 0
    while i < len(smiles):
        char = smiles[i]
        if char == "[":
            close = smiles.find("]", i)
            i = len(smiles) if close == -1 else close + 1
            continue
        if char == "%":
            label = smiles[i + 1 : i + 3]
            if len(label) == 2 and label.isdigit():
                counts[label] = counts.get(label, 0) + 1
            i += 3
            continue
        if char.isdigit():
            counts[char] = counts.get(char, 0) + 1
        i += 1
    return counts


def _ring_labels_balanced(smiles: str) -> bool:
    """Return True when every ring label is opened and closed (even occurrence count)."""
    return all(count % 2 == 0 for count in _ring_label_occurrences(smiles).values())


def _free_ring_label(parent_smiles: str) -> str:
    used = set(_ring_label_occurrences(parent_smiles))
    for label in "1234567890":
        if label not in used:
            return label
    raise BioisostereError("No free single-digit ring label remains in the parent SMILES")


def _replacement_token(parent_smiles: str, rule: BioisostereRule) -> str:
    """Return the replacement, renumbering its ring labels away from the parent's labels."""
    if not any(char.isdigit() or char == "%" for char in rule.replacement_motif):
        return rule.replacement_motif
    token = rule.replacement_motif
    used = sorted(set(_ring_label_occurrences(rule.replacement_motif)))
    for label in used:
        token = token.replace(label, _free_ring_label(parent_smiles))
    return token


def _target_variants(motif: str) -> List[str]:
    """Yield the motif with any ring labels relabelled, so a literal motif matches equivalent ring numbering."""
    labels = sorted({char for char in motif if char.isdigit()})
    if not labels:
        return [motif]
    variants: List[str] = []
    for candidate in "123456789":
        variant = motif
        for label in labels:
            variant = variant.replace(label, candidate)
        variants.append(variant)
    return variants


def _safe_match(smiles: str, rule: BioisostereRule) -> Optional[Tuple[int, int, str]]:
    """Find the first motif occurrence that can be replaced without splitting a ring context."""
    for target in _target_variants(rule.target_motif):
        start = smiles.find(target)
        while start != -1:
            end = start + len(target)
            following = smiles[end] if end < len(smiles) else ""
            # A ring label immediately after the motif belongs to the last matched atom: cutting it out
            # would leave the label attached to an unrelated replacement atom. Labels inside a balanced
            # motif are part of the replaced ring and are safe.
            touches_ring = following.isdigit() or following == "%"
            if not touches_ring:
                replacement = _replacement_token(smiles, rule)
                candidate = smiles[:start] + replacement + smiles[end:]
                if _ring_labels_balanced(candidate):
                    return start, end, replacement
            start = smiles.find(target, start + 1)
    return None


def generate_bioisosteres(
    parent_id: str,
    smiles: str,
    mw: float,
    logp: float,
    tpsa: float,
    hbd: int,
    hba: int,
    rotb: int,
    pka: Optional[float] = None,
) -> List[GeneratedAnalog]:
    """
    Apply rule-based bioisosteric transformations to generate virtual lead analogs.

    Matching is literal substring replacement, not chemical canonicalization. A rule is skipped when
    the motif would be cut out of a ring context or the result would leave unbalanced ring labels.
    Property deltas are rule estimates, not recomputed from the analog structure.
    """
    parent_cns_mpo = calculate_cns_mpo(mw=mw, logp=logp, tpsa=tpsa, hbd=hbd, pka=pka)
    analogs: List[GeneratedAnalog] = []

    for rule in BIOISOSTERE_RULES:
        match = _safe_match(smiles, rule)
        if match is None:
            continue
        start, end, replacement = match
        new_smiles = smiles[:start] + replacement + smiles[end:]

        new_mw = mw + rule.delta_mw
        new_logp = logp + rule.delta_logp
        new_tpsa = max(0.0, tpsa + rule.delta_tpsa)
        new_hbd = max(0, hbd + rule.delta_hbd)
        new_hba = max(0, hba + rule.delta_hba)
        new_rotb = max(0, rotb + rule.delta_rotb)
        new_pka = (pka - 1.5) if (pka is not None and rule.rule_id == "BIO-04") else pka
        new_cns_mpo = calculate_cns_mpo(mw=new_mw, logp=new_logp, tpsa=new_tpsa, hbd=new_hbd, pka=new_pka)
        delta_mpo = new_cns_mpo - parent_cns_mpo

        if new_cns_mpo >= 5.0 and delta_mpo >= 0:
            verdict = "FAVORABLE: Superior CNS drug-likeness with enhanced metabolic shielding."
        elif new_cns_mpo >= 4.5:
            verdict = "VIABLE: Maintained CNS MPO score; recommended for synthetic prioritization."
        else:
            verdict = "MARGINAL: Slight drop in CNS MPO due to increased molecular weight or TPSA."

        analogs.append(
            GeneratedAnalog(
                analog_id=f"{parent_id}-{rule.rule_id.replace('-', '')}",
                parent_id=parent_id,
                rule_id=rule.rule_id,
                transformation_name=rule.name,
                parent_smiles=smiles,
                analog_smiles=new_smiles,
                mw=new_mw,
                logp=new_logp,
                tpsa=new_tpsa,
                hbd=new_hbd,
                hba=new_hba,
                rotb=new_rotb,
                cns_mpo_score=new_cns_mpo,
                delta_cns_mpo=delta_mpo,
                transformation_rationale=rule.rationale,
                optimization_verdict=verdict,
            )
        )

    analogs.sort(key=lambda a: (-a.cns_mpo_score, a.analog_id))
    return analogs


def render_bioisostere_report(analogs: Sequence[GeneratedAnalog], parent_id: str) -> str:
    """Render a bioisostere generation and optimization Markdown report."""
    favorable = [a for a in analogs if "FAVORABLE" in a.optimization_verdict]

    lines = [
        f"# PNN Bioisosteric Replacement & Metabolic Soft-Spot Report: {parent_id}",
        "",
        "## Executive Summary",
        "",
        f"- **Parent Lead Compound**: `{parent_id}`",
        f"- **Generated Virtual Analogs**: {len(analogs)}",
        f"- **Favorable Optimization Candidates (CNS MPO >= 5.0)**: {len(favorable)}",
        "",
        "## Method Limits",
        "",
        "> [!IMPORTANT]",
        "> Transformations are literal SMILES substring replacements, not chemical canonicalization. A rule is "
        "skipped when the motif would be cut out of a ring context. Molecular weight, logP and TPSA are rule-estimated "
        "deltas rather than recomputed from the analog structure, so every analog requires cheminformatics "
        "validation and synthesis before any claim about potency, stability or exposure.",
        "",
        "## Generated Bioisosteric Analogs",
        "",
        "| Analog ID | Rule | Transformation | Analog SMILES | MW (Da) | logP | TPSA (Å²) | CNS MPO | Δ MPO | Optimization Verdict |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]

    for a in analogs:
        lines.append(
            f"| `{a.analog_id}` | {a.rule_id} | **{a.transformation_name}** | `{a.analog_smiles}` | "
            f"{a.mw:.1f} | {a.logp:.2f} | {a.tpsa:.1f} | **{a.cns_mpo_score:.2f}** | "
            f"{a.delta_cns_mpo:+.2f} | *{a.optimization_verdict.split(':')[0]}* |"
        )

    lines.extend(["", "## Detailed Transformation Rationale", ""])
    if not analogs:
        lines.extend(
            [
                "No rule in the table matched the parent structure, so no virtual analog was generated. "
                "This is a coverage limit of the literal-motif rule set, not evidence that the parent is "
                "unoptimizable.",
                "",
            ]
        )
    for a in analogs:
        lines.extend(
            [
                f"### {a.analog_id}: {a.transformation_name}",
                "",
                f"- **Parent Motif**: Replaced in `{a.parent_smiles}`",
                f"- **Optimized SMILES**: `{a.analog_smiles}`",
                f"- **Pharmacological Objective**: {a.transformation_rationale}",
                f"- **Evaluation**: {a.optimization_verdict}",
                "",
            ]
        )

    lines.extend(
        [
            "## Synthetic Priority Recommendations",
            "",
            "1. **Priority 1 (Fluorinated and CF3 Leads)**: Initiate parallel synthesis of para-fluorinated and "
            "para-CF3 analogs, then measure microsomal clearance (HLM/MLM).",
            "2. **Priority 2 (Morpholine Analogs)**: Synthesize basicity-attenuated morpholine analogs to test the "
            "predicted pKa reduction, hERG liability and brain free fraction.",
            "3. **Validation Gate**: Confirm every generated SMILES parses and the predicted property deltas match a "
            "cheminformatics calculation before any analog enters docking or synthesis.",
            "",
        ]
    )
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


def run_bioisostere_pipeline(
    compounds_path: Path,
    target_parent_id: Optional[str] = None,
) -> List[GeneratedAnalog]:
    """Load parent compound and generate bioisosteric analogs across all rules."""
    if not compounds_path.is_file():
        raise BioisostereError(f"Compounds library not found: {compounds_path}")

    all_analogs: List[GeneratedAnalog] = []
    with open(compounds_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise BioisostereError(f"Empty compounds CSV: {compounds_path}")
        norm = {c.strip().lower(): c for c in reader.fieldnames}
        cid_col = norm.get("compound_id")
        smi_col = norm.get("smiles")
        if not cid_col or not smi_col:
            raise BioisostereError(f"Compounds CSV missing compound_id or smiles: {compounds_path}")

        for row in reader:
            cid = row[cid_col].strip()
            smi = row[smi_col].strip()
            if not cid or not smi:
                continue
            if target_parent_id and cid != target_parent_id:
                continue

            mw = _parse_opt_float(row.get(norm.get("mw", "")))
            logp = _parse_opt_float(row.get(norm.get("logp", "")))
            tpsa = _parse_opt_float(row.get(norm.get("tpsa", "")))
            hbd = _parse_opt_float(row.get(norm.get("hbd", "")))
            hba = _parse_opt_float(row.get(norm.get("hba", "")))
            rotb = _parse_opt_float(row.get(norm.get("rotatable_bonds", "")))
            if None in (mw, logp, tpsa, hbd, hba, rotb):
                raise BioisostereError(
                    f"Compound {cid} lacks properties required for bioisostere deltas"
                )
            hbd = int(hbd)
            hba = int(hba)
            rotb = int(rotb)
            pka = _parse_opt_float(row.get(norm.get("pka", "")))

            all_analogs.extend(
                generate_bioisosteres(
                    parent_id=cid,
                    smiles=smi,
                    mw=mw,
                    logp=logp,
                    tpsa=tpsa,
                    hbd=hbd,
                    hba=hba,
                    rotb=rotb,
                    pka=pka,
                )
            )

    if not all_analogs and target_parent_id:
        raise BioisostereError(f"Parent compound '{target_parent_id}' not found or had no applicable bioisostere motifs.")

    return all_analogs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.bioisostere",
        description="Generate bioisosteric lead analogs with metabolic soft-spot shielding and basicity tuning.",
    )
    parser.add_argument("--compounds", required=True, help="path to compound library CSV")
    parser.add_argument("--parent", help="optional parent compound_id to generate analogs for (e.g. PNN-101)")
    parser.add_argument("--out", default="bioisosteric_analogs.csv", help="output destination CSV for generated analogs")
    parser.add_argument("--report", dest="report_path", help="optional destination path for Markdown report")
    parser.add_argument("--json", dest="json_path", help="optional output path for JSON summary")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        analogs = run_bioisostere_pipeline(
            compounds_path=Path(args.compounds),
            target_parent_id=args.parent,
        )

        out_path = Path(args.out)
        write_csv(out_path, ANALOG_FIELDS, [a.as_row() for a in analogs])

        if args.report_path:
            p_label = args.parent or "Candidate Library"
            write_text(Path(args.report_path), render_bioisostere_report(analogs, p_label))

        if args.json_path:
            summary = {
                "total_generated": len(analogs),
                "parent_id": args.parent,
                "analogs": [a.as_row() for a in analogs],
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), summary)

    except BioisostereError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: bioisostere generation requires valid SMILES strings and recognizable chemical motifs.", file=sys.stderr)
        print("Fix: verify that the compound library contains valid SMILES and compound_id headers.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: bioisostere output files could not be written to disk.", file=sys.stderr)
        print("Fix: specify a writable destination for --out, --report, or --json.", file=sys.stderr)
        return 1

    favorable_count = sum(1 for a in analogs if "FAVORABLE" in a.optimization_verdict)
    print(
        f"Bioisostere generation complete: {len(analogs)} analog(s) created ({favorable_count} favorable) -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
