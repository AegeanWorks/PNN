"""Target receptor profiles, binding-pocket invariants and Vina grid generation for PNN."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class TargetProfile:
    """Molecular target profile defining receptor invariants and docking grid coordinates."""

    target_id: str
    name: str
    gene_symbol: str
    uniprot_id: str
    primary_pdb: str
    resolution_angstrom: float
    domain: str
    pocket_type: str
    grid_center: Tuple[float, float, float]
    grid_size: Tuple[float, float, float]
    key_residues: Tuple[str, ...]
    druggability: str
    biological_function: str


TARGET_PROFILES: Dict[str, TargetProfile] = {
    "PTPsigma": TargetProfile(
        target_id="PTPsigma",
        name="Protein Tyrosine Phosphatase Sigma (PTPRS)",
        gene_symbol="PTPRS",
        uniprot_id="Q13332",
        primary_pdb="5K4P",
        resolution_angstrom=2.8,
        domain="D1 Ectodomain / CSPG Binding Site",
        pocket_type="cspg_binding",
        grid_center=(12.5, 18.3, 22.1),
        grid_size=(20.0, 20.0, 20.0),
        key_residues=("Lys201", "Arg204", "Lys207", "Arg210", "Tyr199", "Phe203"),
        druggability="Moderate",
        biological_function=(
            "Primary neuronal receptor for CSPGs; binding triggers receptor dimerization, "
            "inhibiting neurite outgrowth and dephosphorylating the insulin receptor / IRS-1."
        ),
    ),
    "PTPsigma_catalytic": TargetProfile(
        target_id="PTPsigma_catalytic",
        name="PTPsigma Catalytic Domain (PTPRS D1)",
        gene_symbol="PTPRS",
        uniprot_id="Q13332",
        primary_pdb="2PTS",
        resolution_angstrom=2.0,
        domain="Intracellular D1 Phosphatase Catalytic Pocket",
        pocket_type="catalytic",
        grid_center=(24.0, 15.2, -8.4),
        grid_size=(18.0, 18.0, 18.0),
        key_residues=("Cys267", "Arg268", "Asp326"),
        druggability="High",
        biological_function=(
            "Intracellular catalytic pocket hydrolyzing phosphotyrosine on IRS-1 and downstream adaptors."
        ),
    ),
    "LAR": TargetProfile(
        target_id="LAR",
        name="Leukocyte Common Antigen-Related Phosphatase (PTPRF)",
        gene_symbol="PTPRF",
        uniprot_id="P10586",
        primary_pdb="2LAR",
        resolution_angstrom=1.9,
        domain="D1 Catalytic / CSPG Co-Receptor Domain",
        pocket_type="catalytic",
        grid_center=(18.1, 10.4, 15.8),
        grid_size=(20.0, 20.0, 20.0),
        key_residues=("Cys1522", "Arg1528", "Asp1490"),
        druggability="High",
        biological_function=(
            "Co-receptor for CSPG inhibition; operates cooperatively with PTPsigma to suppress "
            "hypothalamic insulin signaling in diet-induced obesity."
        ),
    ),
    "NgR1": TargetProfile(
        target_id="NgR1",
        name="Nogo-66 Receptor 1 (RTN4R)",
        gene_symbol="RTN4R",
        uniprot_id="Q9BZF1",
        primary_pdb="2OZ7",
        resolution_angstrom=2.4,
        domain="Leucine-Rich Repeat (LRR) Ectodomain",
        pocket_type="allosteric_lrr",
        grid_center=(-4.5, 32.1, 12.0),
        grid_size=(24.0, 24.0, 24.0),
        key_residues=("Arg68", "Asn115", "Arg140"),
        druggability="Moderate",
        biological_function=(
            "High-affinity GPI-anchored co-receptor signaling growth cone collapse and restricting "
            "plasticity upon CSPG engagement."
        ),
    ),
    "MMP2": TargetProfile(
        target_id="MMP2",
        name="Matrix Metalloproteinase-2 (72 kDa Type IV Collagenase)",
        gene_symbol="MMP2",
        uniprot_id="P08253",
        primary_pdb="1CK7",
        resolution_angstrom=1.5,
        domain="Zinc-Dependent Catalytic Cleft",
        pocket_type="catalytic",
        grid_center=(14.2, 28.5, 34.0),
        grid_size=(18.0, 18.0, 18.0),
        key_residues=("His403", "His407", "His413", "Glu404"),
        druggability="High",
        biological_function=(
            "Endogenous enzyme responsible for physiological turnover and cleavage of CSPG core proteins."
        ),
    ),
    "ChABC": TargetProfile(
        target_id="ChABC",
        name="Chondroitinase ABC Lyase (Proteus vulgaris)",
        gene_symbol="PL8",
        uniprot_id="P59807",
        primary_pdb="1L5D",
        resolution_angstrom=1.8,
        domain="Polysaccharide Lyase Catalytic Domain",
        pocket_type="catalytic",
        grid_center=(30.5, 42.1, 15.6),
        grid_size=(22.0, 22.0, 22.0),
        key_residues=("His200", "Tyr247", "Arg287"),
        druggability="Moderate",
        biological_function=(
            "Bacterial polysaccharide lyase serving as the benchmark standard for complete PNN degradation; "
            "small-molecule mimics target its catalytic cleft mechanism to achieve non-enzymatic cleavage."
        ),
    ),
    "ADAMTS4": TargetProfile(
        target_id="ADAMTS4",
        name="ADAMTS-4 Aggrecanase-1",
        gene_symbol="ADAMTS4",
        uniprot_id="O75173",
        primary_pdb="2RJP",
        resolution_angstrom=1.5,
        domain="Zinc Metalloproteinase Catalytic & Disintegrin Domain",
        pocket_type="catalytic_zinc",
        grid_center=(28.4, -12.1, 45.6),
        grid_size=(20.0, 20.0, 20.0),
        key_residues=("His202", "His206", "His212", "Glu203", "Leu165", "Phe198"),
        druggability="High",
        biological_function=(
            "Endogenous CNS aggrecanase responsible for interglobular cleavage of CSPGs (aggrecan, brevican, "
            "neurocan); allosteric activation or selective engagement induces physiological PNN softening."
        ),
    ),
    "MMP9": TargetProfile(
        target_id="MMP9",
        name="Matrix Metalloproteinase-9 (92 kDa Gelatinase B)",
        gene_symbol="MMP9",
        uniprot_id="P14780",
        primary_pdb="1GKC",
        resolution_angstrom=2.3,
        domain="Zinc Catalytic Domain with Fibronectin Type II Repeats",
        pocket_type="catalytic_zinc",
        grid_center=(38.2, 16.5, 24.8),
        grid_size=(20.0, 20.0, 20.0),
        key_residues=("His401", "His405", "His411", "Glu402", "Leu188", "Met422"),
        druggability="High",
        biological_function=(
            "Activity-dependent synaptic matrix metalloproteinase in the hypothalamus; cleaves PNN CSPG components "
            "and matrix scaffolds to restore insulin signaling and facilitate synaptic plasticity."
        ),
    ),
}


def get_target(target_id: str) -> TargetProfile:
    """Retrieve a target profile by case-insensitive ID or alias."""
    normalized = target_id.strip()
    for key, profile in TARGET_PROFILES.items():
        if key.lower() == normalized.lower():
            return profile
    valid = ", ".join(TARGET_PROFILES.keys())
    raise KeyError(f"Unknown target ID '{target_id}'. Valid targets are: {valid}")


def list_targets() -> List[TargetProfile]:
    """Return all registered PNN target profiles."""
    return list(TARGET_PROFILES.values())


def generate_vina_config(
    profile: TargetProfile,
    receptor_pdbqt: str = "receptor.pdbqt",
    out_pdbqt: str = "docked_output.pdbqt",
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: int = 3,
) -> str:
    """Format AutoDock Vina grid parameter configuration text."""
    cx, cy, cz = profile.grid_center
    sx, sy, sz = profile.grid_size
    lines = [
        f"# AutoDock Vina configuration for PNN target: {profile.name} ({profile.target_id})",
        f"# Primary PDB structure: {profile.primary_pdb} ({profile.resolution_angstrom} A resolution)",
        f"# Binding pocket: {profile.pocket_type} [{profile.domain}]",
        f"# Key residues: {', '.join(profile.key_residues)}",
        "",
        f"receptor = {receptor_pdbqt}",
        f"out = {out_pdbqt}",
        "",
        "# Grid Center (Angstrom)",
        f"center_x = {cx:.3f}",
        f"center_y = {cy:.3f}",
        f"center_z = {cz:.3f}",
        "",
        "# Grid Dimensions (Angstrom)",
        f"size_x = {sx:.1f}",
        f"size_y = {sy:.1f}",
        f"size_z = {sz:.1f}",
        "",
        "# Search Parameters",
        f"exhaustiveness = {exhaustiveness}",
        f"num_modes = {num_modes}",
        f"energy_range = {energy_range}",
        "",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.target",
        description="Inspect molecular target profiles and export AutoDock Vina grid configurations.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: list
    subparsers.add_parser("list", help="list all registered PNN receptor targets")

    # Subcommand: show
    show_parser = subparsers.add_parser("show", help="display detailed receptor invariants for a target")
    show_parser.add_argument("target_id", help="target identifier (e.g. PTPsigma, LAR, NgR1)")

    # Subcommand: config
    config_parser = subparsers.add_parser("config", help="export AutoDock Vina grid config for a target")
    config_parser.add_argument("target_id", help="target identifier (e.g. PTPsigma, LAR, NgR1)")
    config_parser.add_argument("--receptor", default="receptor.pdbqt", help="receptor PDBQT file path")
    config_parser.add_argument("--out", default="docked.pdbqt", help="docking output PDBQT path")
    config_parser.add_argument("--save", help="optional destination path to write the config file")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "list":
        targets = list_targets()
        print(f"Registered PNN Target Profiles ({len(targets)}):")
        print("-" * 78)
        print(f"{'Target ID':<18} {'PDB':<6} {'Res(A)':<8} {'Druggability':<14} {'Domain'}")
        print("-" * 78)
        for t in targets:
            print(f"{t.target_id:<18} {t.primary_pdb:<6} {t.resolution_angstrom:<8.1f} {t.druggability:<14} {t.domain}")
        return 0

    if args.command == "show":
        try:
            profile = get_target(args.target_id)
        except KeyError as exc:
            print(f"What: {exc}", file=sys.stderr)
            print("Why: unknown target identifier.", file=sys.stderr)
            print("Fix: run 'python -m pnn.target list' to see available target IDs.", file=sys.stderr)
            return 1

        print(f"Target:       {profile.name} [{profile.target_id}]")
        print(f"Gene Symbol:  {profile.gene_symbol} (UniProt: {profile.uniprot_id})")
        print(f"Structure:    PDB {profile.primary_pdb} ({profile.resolution_angstrom} A resolution)")
        print(f"Domain:       {profile.domain}")
        print(f"Pocket Type:  {profile.pocket_type} (Druggability: {profile.druggability})")
        print(f"Grid Center:  X={profile.grid_center[0]}, Y={profile.grid_center[1]}, Z={profile.grid_center[2]}")
        print(f"Grid Size:    X={profile.grid_size[0]}, Y={profile.grid_size[1]}, Z={profile.grid_size[2]}")
        print(f"Key Residues: {', '.join(profile.key_residues)}")
        print(f"Function:     {profile.biological_function}")
        return 0

    if args.command == "config":
        try:
            profile = get_target(args.target_id)
        except KeyError as exc:
            print(f"What: {exc}", file=sys.stderr)
            print("Why: unknown target identifier.", file=sys.stderr)
            print("Fix: run 'python -m pnn.target list' to see available target IDs.", file=sys.stderr)
            return 1

        config_text = generate_vina_config(
            profile=profile,
            receptor_pdbqt=args.receptor,
            out_pdbqt=args.out,
        )
        if args.save:
            dest = Path(args.save)
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(config_text, encoding="utf-8")
            print(f"Saved Vina grid config for '{profile.target_id}' -> {dest}")
        else:
            print(config_text)
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
