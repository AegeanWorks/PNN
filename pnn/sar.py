"""Structure-Activity Relationship (SAR) and Bemis-Murcko scaffold decomposition for PNN."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from .reports import write_csv, write_json, write_text


class SarError(ValueError):
    """Raised when SAR data is malformed or violates scientific contract."""


@dataclass(frozen=True)
class CompoundSar:
    compound_id: str
    name: str
    smiles: str
    scaffold_id: str
    scaffold_name: str
    core_smiles: str
    r_groups: str
    mw: float
    logp: float
    tpsa: float
    vina_score: Optional[float]
    ec50_nm: Optional[float]
    brain_kp: Optional[float]
    oral_f_pct: Optional[float]
    activity_tier: str

    def as_row(self) -> Dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "name": self.name,
            "scaffold_id": self.scaffold_id,
            "scaffold_name": self.scaffold_name,
            "core_smiles": self.core_smiles,
            "r_groups": self.r_groups,
            "mw": f"{self.mw:.1f}",
            "logp": f"{self.logp:.2f}",
            "tpsa": f"{self.tpsa:.1f}",
            "vina_score": f"{self.vina_score:.2f}" if self.vina_score is not None else "",
            "ec50_nm": f"{self.ec50_nm:.1f}" if self.ec50_nm is not None else "",
            "brain_kp": f"{self.brain_kp:.2f}" if self.brain_kp is not None else "",
            "oral_f_pct": f"{self.oral_f_pct:.1f}" if self.oral_f_pct is not None else "",
            "activity_tier": self.activity_tier,
        }


@dataclass(frozen=True)
class ScaffoldCluster:
    scaffold_id: str
    scaffold_name: str
    core_smiles: str
    compound_count: int
    compound_ids: str
    best_compound: str
    best_ec50_nm: Optional[float]
    mean_ec50_nm: Optional[float]
    best_vina_docking: Optional[float]
    mean_brain_kp: Optional[float]
    sar_category: str
    key_sar_insights: str

    def as_row(self) -> Dict[str, object]:
        return {
            "scaffold_id": self.scaffold_id,
            "scaffold_name": self.scaffold_name,
            "core_smiles": self.core_smiles,
            "compound_count": self.compound_count,
            "compound_ids": self.compound_ids,
            "best_compound": self.best_compound,
            "best_ec50_nm": f"{self.best_ec50_nm:.1f}" if self.best_ec50_nm is not None else "",
            "mean_ec50_nm": f"{self.mean_ec50_nm:.1f}" if self.mean_ec50_nm is not None else "",
            "best_vina_docking": f"{self.best_vina_docking:.2f}" if self.best_vina_docking is not None else "",
            "mean_brain_kp": f"{self.mean_brain_kp:.2f}" if self.mean_brain_kp is not None else "",
            "sar_category": self.sar_category,
            "key_sar_insights": self.key_sar_insights,
        }


SCAFFOLD_CLUSTER_FIELDS = (
    "scaffold_id",
    "scaffold_name",
    "core_smiles",
    "compound_count",
    "compound_ids",
    "best_compound",
    "best_ec50_nm",
    "mean_ec50_nm",
    "best_vina_docking",
    "mean_brain_kp",
    "sar_category",
    "key_sar_insights",
)


_TOKEN_SPEC = [
    ("BRACKET", r"\[[^\]]+\]"),
    ("HALOGEN", r"Cl|Br"),
    ("ATOM", r"[BCNOPSbcnopsFI]"),
    ("RING_2", r"%\d{2}"),
    ("RING_1", r"\d"),
    ("BOND", r"[-=#:/\\.]"),
    ("BRANCH_OPEN", r"\("),
    ("BRANCH_CLOSE", r"\)"),
]
_TOK_REGEX = re.compile("|".join(f"(?P<{name}>{pattern})" for name, pattern in _TOKEN_SPEC))


def parse_smiles_graph(smiles: str) -> Tuple[List[str], Dict[int, Set[int]]]:
    """Parse SMILES into an atom list and adjacency graph."""
    clean = smiles.strip()
    if not clean:
        raise SarError("Empty SMILES string")

    atoms: List[str] = []
    adj: Dict[int, Set[int]] = {}
    ring_open: Dict[int, int] = {}
    branch_stack: List[int] = []
    prev_atom: Optional[int] = None
    pending_bond: Optional[str] = None

    for m in _TOK_REGEX.finditer(clean):
        kind = m.lastgroup
        val = m.group()
        if kind in ("BRACKET", "HALOGEN", "ATOM"):
            idx = len(atoms)
            atoms.append(val)
            adj[idx] = set()
            if prev_atom is not None and pending_bond != ".":
                adj[prev_atom].add(idx)
                adj[idx].add(prev_atom)
            prev_atom = idx
            pending_bond = None
        elif kind == "BOND":
            pending_bond = val
        elif kind in ("RING_1", "RING_2"):
            ring_id = int(val.replace("%", ""))
            if ring_id in ring_open:
                other_idx = ring_open.pop(ring_id)
                adj[other_idx].add(prev_atom)  # type: ignore[arg-type]
                adj[prev_atom].add(other_idx)  # type: ignore[index]
            else:
                ring_open[ring_id] = prev_atom  # type: ignore[assignment]
                pending_bond = None
        elif kind == "BRANCH_OPEN":
            if prev_atom is None:
                raise SarError(f"Branch open with no preceding atom in SMILES: {smiles}")
            branch_stack.append(prev_atom)
        elif kind == "BRANCH_CLOSE":
            if not branch_stack:
                raise SarError(f"Unmatched closing parenthesis in SMILES: {smiles}")
            prev_atom = branch_stack.pop()

    if not atoms:
        raise SarError(f"No valid atoms parsed from SMILES: {smiles}")
    if branch_stack:
        raise SarError(f"Unclosed branch parenthesis in SMILES: {smiles}")
    if ring_open:
        raise SarError(f"Unclosed ring identifier in SMILES: {smiles}")

    return atoms, adj


def find_biconnected_blocks(n: int, adj: Dict[int, Set[int]]) -> List[Set[int]]:
    """Hopcroft-Tarjan algorithm to find all biconnected components (cycle systems)."""
    visited = [False] * n
    depth = [0] * n
    low = [0] * n
    parent = [-1] * n
    stack: List[Tuple[int, int]] = []
    blocks: List[Set[int]] = []

    def dfs(u: int, d: int) -> None:
        visited[u] = True
        depth[u] = d
        low[u] = d
        for v in adj.get(u, set()):
            if v == parent[u]:
                continue
            if visited[v]:
                low[u] = min(low[u], depth[v])
                if depth[v] < depth[u]:
                    stack.append((u, v))
            else:
                parent[v] = u
                stack.append((u, v))
                dfs(v, d + 1)
                low[u] = min(low[u], low[v])
                if low[v] >= depth[u]:
                    block: Set[int] = set()
                    while stack:
                        edge = stack.pop()
                        block.add(edge[0])
                        block.add(edge[1])
                        if edge == (u, v):
                            break
                    if len(block) >= 3:
                        blocks.append(block)

    for i in range(n):
        if not visited[i]:
            dfs(i, 0)

    return blocks


def serialize_subgraph_smiles(atoms: Sequence[str], adj: Dict[int, Set[int]], nodes: Sequence[int]) -> str:
    """Serialize an induced subgraph of atoms into a SMILES string."""
    if not nodes:
        return ""
    node_set = set(nodes)
    sub_adj = {u: [v for v in adj[u] if v in node_set] for u in node_set}

    root = min(node_set, key=lambda u: (len(sub_adj[u]), u))
    visited: Set[int] = set()
    ring_closures: Dict[Tuple[int, int], int] = {}
    next_ring = 1

    tree_children: Dict[int, List[int]] = {u: [] for u in node_set}
    ring_edges: Dict[int, List[Tuple[int, int]]] = {u: [] for u in node_set}

    def dfs(u: int, p: Optional[int]) -> None:
        visited.add(u)
        for v in sub_adj[u]:
            if v == p:
                continue
            if v in visited:
                if (v, u) not in ring_closures:
                    nonlocal next_ring
                    rid = next_ring
                    next_ring += 1
                    ring_closures[(u, v)] = rid
                    ring_edges[u].append((v, rid))
                    ring_edges[v].append((u, rid))
            else:
                tree_children[u].append(v)
                dfs(v, u)

    dfs(root, None)

    def to_smiles(u: int) -> str:
        out = atoms[u]
        for _, rid in sorted(ring_edges[u], key=lambda x: x[1]):
            out += str(rid) if rid < 10 else f"%{rid}"
        children = tree_children[u]
        if not children:
            return out
        for c in children[:-1]:
            out += f"({to_smiles(c)})"
        out += to_smiles(children[-1])
        return out

    return to_smiles(root)


def _identify_scaffold_name(atoms: Sequence[str], cycle_blocks: Sequence[Set[int]], linkers: Set[int]) -> str:
    """Identify recognizable chemotype descriptor for Bemis-Murcko framework."""
    ring_names: List[str] = []
    for block in cycle_blocks:
        syms = [atoms[idx] for idx in block]
        has_n = "N" in syms or "n" in syms
        has_o = "O" in syms or "o" in syms
        has_s = "S" in syms or "s" in syms
        is_aromatic = any(s.islower() for s in syms)
        size = len(block)

        if size == 6:
            if has_n and has_o:
                ring_names.append("Morpholine")
            elif syms.count("N") >= 2 or syms.count("n") >= 2:
                ring_names.append("Piperazine" if not is_aromatic else "Pyrimidine")
            elif has_n:
                ring_names.append("Pyridine" if is_aromatic else "Piperidine")
            elif is_aromatic:
                ring_names.append("Phenyl")
            else:
                ring_names.append("Cyclohexyl")
        elif size == 5:
            if has_s and has_n:
                ring_names.append("Thiazole")
            elif has_o:
                ring_names.append("Furan")
            elif has_n:
                ring_names.append("Pyrrole" if is_aromatic else "Pyrrolidine")
            else:
                ring_names.append("Cyclopentyl")
        elif size == 3:
            ring_names.append("Cyclopropyl")
        elif size == 4:
            ring_names.append("Cyclobutyl")
        elif size > 6:
            ring_names.append(f"FusedPolycycle-{size}")
        else:
            ring_names.append(f"Ring-{size}")

    linker_syms = "".join(atoms[idx] for idx in sorted(linkers))
    linker_name = ""
    if "S" in linker_name or "S" in linker_syms:
        linker_name = "Sulfonamide"
    elif ("N" in linker_syms or "n" in linker_syms) and ("C" in linker_syms or "c" in linker_syms) and ("O" in linker_syms or "o" in linker_syms):
        linker_name = "Carboxamide"
    elif "O" in linker_syms:
        linker_name = "Ether"
    elif "N" in linker_syms:
        linker_name = "Amine"
    elif linker_syms:
        linker_name = "AlkylLinker"

    if ring_names:
        unique_rings = sorted(set(ring_names))
        base = "-".join(unique_rings)
        return f"{base}_{linker_name}" if linker_name else base
    return "Acyclic-Chain"


def extract_scaffold(smiles: str) -> Tuple[str, str, str]:
    """
    Decompose a SMILES into Bemis-Murcko framework SMILES, chemotype name, and R-group substituents.
    Returns: (core_smiles, scaffold_name, r_groups_str)
    """
    atoms, adj = parse_smiles_graph(smiles)
    blocks = find_biconnected_blocks(len(atoms), adj)

    if not blocks:
        # Acyclic molecule
        return smiles, "Acyclic-Scaffold", "None"

    all_cycle_atoms: Set[int] = set().union(*blocks)
    non_cycle = set(range(len(atoms))) - all_cycle_atoms

    visited_nc: Set[int] = set()
    linkers: Set[int] = set()
    substituent_strs: List[str] = []

    for node in non_cycle:
        if node not in visited_nc:
            comp: List[int] = []
            queue = [node]
            visited_nc.add(node)
            while queue:
                curr = queue.pop()
                comp.append(curr)
                for nbr in adj[curr]:
                    if nbr in non_cycle and nbr not in visited_nc:
                        visited_nc.add(nbr)
                        queue.append(nbr)

            attached_cycle_atoms: Set[int] = set()
            attached_blocks: Set[int] = set()
            for c_node in comp:
                for nbr in adj[c_node]:
                    if nbr in all_cycle_atoms:
                        attached_cycle_atoms.add(nbr)
                        for b_idx, block in enumerate(blocks):
                            if nbr in block:
                                attached_blocks.add(b_idx)

            comp_str = "".join(atoms[idx] for idx in comp)
            if len(attached_blocks) >= 2 or len(attached_cycle_atoms) >= 2:
                linkers.update(comp)
            else:
                substituent_strs.append(comp_str)

    framework_atoms = sorted(all_cycle_atoms.union(linkers))
    core_smiles = serialize_subgraph_smiles(atoms, adj, framework_atoms)
    scaffold_name = _identify_scaffold_name(atoms, blocks, linkers)
    r_groups = ", ".join(substituent_strs) if substituent_strs else "None"

    return core_smiles, scaffold_name, r_groups


def load_compounds(path: Path) -> List[Dict[str, str]]:
    """Load compound property records from CSV."""
    if not path.is_file():
        raise SarError(f"Compound dataset not found: {path}")
    records: List[Dict[str, str]] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SarError(f"Empty compound CSV: {path}")
        norm = {c.strip().lower(): c for c in reader.fieldnames}
        req = {"compound_id", "smiles"}
        if req - set(norm.keys()):
            raise SarError(f"Compound CSV missing compound_id or smiles: {path}")
        for row in reader:
            cid = row[norm["compound_id"]].strip()
            smi = row[norm["smiles"]].strip()
            if not cid or not smi:
                continue
            records.append({k: row[v].strip() for k, v in norm.items()})
    return records


def load_docking_scores(path: Path) -> Dict[str, float]:
    """Load multi-target or single-target docking scores (consensus or best affinity)."""
    if not path.is_file():
        return {}
    scores: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        norm = {c.strip().lower(): c for c in reader.fieldnames}
        cid_col = norm.get("compound_id")
        # Look for score or consensus column
        score_col = norm.get("consensus_score") or norm.get("vina_affinity_kcal_mol") or norm.get("score")
        if not cid_col or not score_col:
            return {}
        for row in reader:
            cid = row[cid_col].strip()
            val_str = row[score_col].strip()
            try:
                scores[cid] = float(val_str)
            except ValueError:
                continue
    return scores


def load_assay_potencies(path: Path) -> Dict[str, float]:
    """Load in vitro assay EC50 potencies."""
    if not path.is_file():
        return {}
    potencies: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        norm = {c.strip().lower(): c for c in reader.fieldnames}
        cid_col = norm.get("compound_id")
        ec_col = norm.get("ec50_nm") or norm.get("pnn_ec50_nm")
        if not cid_col or not ec_col:
            return {}
        for row in reader:
            cid = row[cid_col].strip()
            val_str = row[ec_col].strip()
            try:
                potencies[cid] = float(val_str)
            except ValueError:
                continue
    return potencies


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


def load_preclinical_exposures(path: Path) -> Dict[str, Dict[str, float]]:
    """Load in vivo preclinical exposures (brain_plasma_ratio, oral_bioavailability_pct)."""
    if not path.is_file():
        return {}
    exposures: Dict[str, Dict[str, float]] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        norm = {c.strip().lower(): c for c in reader.fieldnames}
        cid_col = norm.get("compound_id")
        if not cid_col:
            return {}
        for row in reader:
            cid = row[cid_col].strip()
            data: Dict[str, float] = {}
            if "brain_plasma_ratio" in norm:
                kp = _parse_opt_float(row[norm["brain_plasma_ratio"]])
                if kp is not None:
                    data["brain_kp"] = kp
            if "oral_bioavailability_pct" in norm:
                f_pct = _parse_opt_float(row[norm["oral_bioavailability_pct"]])
                if f_pct is not None:
                    data["oral_f"] = f_pct
            if data:
                exposures[cid] = data
    return exposures


def analyze_sar(
    compound_records: Sequence[Dict[str, str]],
    docking_scores: Optional[Dict[str, float]] = None,
    assay_potencies: Optional[Dict[str, float]] = None,
    preclinical_exposures: Optional[Dict[str, Dict[str, float]]] = None,
) -> Tuple[List[CompoundSar], List[ScaffoldCluster]]:
    """Decompose compounds into Bemis-Murcko scaffolds and perform SAR activity clustering."""
    dock_map = docking_scores or {}
    ec50_map = assay_potencies or {}
    pk_map = preclinical_exposures or {}

    compound_sars: List[CompoundSar] = []
    scaffold_map: Dict[str, List[CompoundSar]] = {}
    scaffold_meta: Dict[str, Tuple[str, str]] = {}  # scf_id -> (scaffold_name, core_smiles)
    scf_counter = 1
    scf_id_by_smiles: Dict[str, str] = {}

    for row in compound_records:
        cid = row["compound_id"]
        smi = row["smiles"]
        name = row.get("name", cid)
        mw = float(row.get("mw", 0.0))
        logp = float(row.get("logp", 0.0))
        tpsa = float(row.get("tpsa", 0.0))

        core_smi, scf_name, r_groups = extract_scaffold(smi)

        if core_smi not in scf_id_by_smiles:
            scf_id = f"SCF-{scf_counter:02d}"
            scf_counter += 1
            scf_id_by_smiles[core_smi] = scf_id
            scaffold_meta[scf_id] = (scf_name, core_smi)
        else:
            scf_id = scf_id_by_smiles[core_smi]

        vina = dock_map.get(cid)
        ec50 = ec50_map.get(cid)
        if ec50 is None and "pnn_ec50_nm" in row:
            ec50 = _parse_opt_float(row["pnn_ec50_nm"])

        brain_kp = None
        oral_f = None
        if cid in pk_map:
            brain_kp = pk_map[cid].get("brain_kp")
            oral_f = pk_map[cid].get("oral_f")
        if brain_kp is None and "brain_plasma_ratio" in row:
            brain_kp = _parse_opt_float(row["brain_plasma_ratio"])
        if oral_f is None and "oral_bioavailability" in row:
            f_frac = _parse_opt_float(row["oral_bioavailability"])
            if f_frac is not None:
                oral_f = f_frac * 100.0  # normalize fraction to %

        if ec50 is not None:
            if ec50 <= 50.0:
                tier = "HIGH_POTENCY (<50 nM)"
            elif ec50 <= 500.0:
                tier = "MODERATE_POTENCY (50-500 nM)"
            else:
                tier = "LOW_POTENCY (>500 nM)"
        else:
            tier = "UNTESTED"

        csar = CompoundSar(
            compound_id=cid,
            name=name,
            smiles=smi,
            scaffold_id=scf_id,
            scaffold_name=scf_name,
            core_smiles=core_smi,
            r_groups=r_groups,
            mw=mw,
            logp=logp,
            tpsa=tpsa,
            vina_score=vina,
            ec50_nm=ec50,
            brain_kp=brain_kp,
            oral_f_pct=oral_f,
            activity_tier=tier,
        )
        compound_sars.append(csar)
        scaffold_map.setdefault(scf_id, []).append(csar)

    # Build Scaffold Clusters
    clusters: List[ScaffoldCluster] = []
    for scf_id, analogs in scaffold_map.items():
        scf_name, core_smi = scaffold_meta[scf_id]
        cids = [a.compound_id for a in analogs]
        tested_ec50 = [a.ec50_nm for a in analogs if a.ec50_nm is not None]
        tested_vina = [a.vina_score for a in analogs if a.vina_score is not None]
        tested_kp = [a.brain_kp for a in analogs if a.brain_kp is not None]

        best_ec50 = min(tested_ec50) if tested_ec50 else None
        mean_ec50 = sum(tested_ec50) / len(tested_ec50) if tested_ec50 else None
        best_vina = min(tested_vina) if tested_vina else None
        mean_kp = sum(tested_kp) / len(tested_kp) if tested_kp else None

        best_analog = min(
            analogs,
            key=lambda a: (a.ec50_nm if a.ec50_nm is not None else 1e9, a.vina_score if a.vina_score is not None else 0.0),
        )

        # Categorize SAR
        insights: List[str] = []
        if best_ec50 is not None and best_ec50 <= 50.0 and mean_kp is not None and mean_kp >= 0.40:
            category = "PRIORITY_LEAD_SERIES"
            insights.append("Nanomolar target potency combined with high CNS/brain penetration.")
        elif best_ec50 is not None and best_ec50 <= 150.0:
            category = "VIABLE_OPTIMIZATION_SCAFFOLD"
            insights.append("Sub-micromolar potency; opportunities exist for R-group optimization.")
        elif mean_kp is not None and mean_kp < 0.20:
            category = "PERIPHERAL_RESTRICTED"
            insights.append("Restricted blood-brain barrier transport; requires lipophilic/polar surface balance.")
        else:
            category = "EXPLORATORY_SERIES"
            insights.append("Moderate affinity or uncharacterized brain pharmacokinetic exposure.")

        if "Piperazine" in scf_name:
            insights.append("Basic nitrogen center confers electrostatic interaction with catalytic cleft.")
        if "Sulfonamide" in scf_name:
            insights.append("Sulfonamide hinge provides planar docking alignment but lower brain permeability.")

        clusters.append(
            ScaffoldCluster(
                scaffold_id=scf_id,
                scaffold_name=scf_name,
                core_smiles=core_smi,
                compound_count=len(analogs),
                compound_ids="; ".join(cids),
                best_compound=best_analog.compound_id,
                best_ec50_nm=best_ec50,
                mean_ec50_nm=mean_ec50,
                best_vina_docking=best_vina,
                mean_brain_kp=mean_kp,
                sar_category=category,
                key_sar_insights=" ".join(insights),
            )
        )

    # Sort clusters by priority (best EC50 ascending, count descending)
    clusters.sort(key=lambda cl: (cl.best_ec50_nm if cl.best_ec50_nm is not None else 1e9, -cl.compound_count))
    return compound_sars, clusters


def render_sar_report(compounds: Sequence[CompoundSar], clusters: Sequence[ScaffoldCluster]) -> str:
    """Render a comprehensive Structure-Activity Relationship (SAR) markdown report."""
    lines = [
        "# PNN Structure-Activity Relationship (SAR) & Scaffold Analysis Report",
        "",
        "## Executive Summary",
        "",
        f"- **Total Compounds Analyzed**: {len(compounds)}",
        f"- **Unique Bemis-Murcko Scaffolds**: {len(clusters)}",
        f"- **Priority Lead Series**: {sum(1 for c in clusters if c.sar_category == 'PRIORITY_LEAD_SERIES')}",
        f"- **Viable Optimization Scaffolds**: {sum(1 for c in clusters if c.sar_category == 'VIABLE_OPTIMIZATION_SCAFFOLD')}",
        "",
        "## Scientific Context & Guidance Invariants",
        "",
        "> [!NOTE]",
        "> **Consensus & Exact Representation Invariants (G-002, G-003)**: Chemical scaffolds are decomposed into "
        "Bemis-Murcko frameworks (rings and linkers) while preserving original compound SMILES identifiers. Docking affinities "
        "represent rank-ordering indicators relative to the screened library rather than absolute thermodynamic free energies.",
        "",
        "## Scaffold Cluster Ranking",
        "",
        "| Scaffold ID | Chemotype Core | Count | Best Analog | Best EC50 (nM) | Mean EC50 (nM) | Best Vina | Brain Kp | SAR Category |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for cl in clusters:
        bec50 = f"{cl.best_ec50_nm:.1f}" if cl.best_ec50_nm is not None else "-"
        mec50 = f"{cl.mean_ec50_nm:.1f}" if cl.mean_ec50_nm is not None else "-"
        bvina = f"{cl.best_vina_docking:.2f}" if cl.best_vina_docking is not None else "-"
        mkp = f"{cl.mean_brain_kp:.2f}" if cl.mean_brain_kp is not None else "-"
        lines.append(
            f"| `{cl.scaffold_id}` | **{cl.scaffold_name}** | {cl.compound_count} | `{cl.best_compound}` | "
            f"{bec50} | {mec50} | {bvina} | {mkp} | *{cl.sar_category}* |"
        )

    lines.extend([
        "",
        "## Detailed Chemotype Insights",
        "",
    ])

    for cl in clusters:
        lines.extend([
            f"### {cl.scaffold_id}: {cl.scaffold_name}",
            "",
            f"- **Core Framework SMILES**: `{cl.core_smiles}`",
            f"- **Included Compounds**: `{cl.compound_ids}`",
            f"- **SAR Evaluation**: {cl.key_sar_insights}",
            "",
        ])

    lines.extend([
        "## Medicinal Chemistry Recommendations",
        "",
        "1. **Core Rigidification**: Introduction of piperazine/pyrimidine linkers stabilizes conformations favorable for PTPsigma catalytic pocket occupancy while preserving brain exposure.",
        "2. **Polar Surface Area Optimization**: Restrict polar surface area (TPSA < 75 Å²) to maintain brain-to-plasma ratios (Kp > 0.40).",
        "3. **Metabolic Soft-Spot Shielding**: Fluorination or methyl substitution on terminal aromatic rings attenuates CYP3A4 clearance while retaining nanomolar PNN degradation efficacy.",
        "",
    ])
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.sar",
        description="Decompose small molecules into Bemis-Murcko scaffolds and perform Structure-Activity Relationship (SAR) analysis.",
    )
    parser.add_argument("--compounds", required=True, help="path to compound library CSV")
    parser.add_argument("--docking", help="optional path to docking scores CSV")
    parser.add_argument("--assays", help="optional path to in vitro degradation assay CSV")
    parser.add_argument("--preclinical", help="optional path to in vivo preclinical candidates or PK CSV")
    parser.add_argument("--out", default="sar_clusters.csv", help="output destination CSV for scaffold clusters")
    parser.add_argument("--compounds-out", dest="compounds_out", help="optional output destination CSV for annotated compound SAR rows")
    parser.add_argument("--report", dest="report_path", help="optional output destination for Markdown SAR report")
    parser.add_argument("--json", dest="json_path", help="optional output path for JSON summary")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        compound_records = load_compounds(Path(args.compounds))
        docking_scores = load_docking_scores(Path(args.docking)) if args.docking else {}
        assay_potencies = load_assay_potencies(Path(args.assays)) if args.assays else {}
        preclinical_exposures = load_preclinical_exposures(Path(args.preclinical)) if args.preclinical else {}

        compounds_sar, clusters = analyze_sar(
            compound_records=compound_records,
            docking_scores=docking_scores,
            assay_potencies=assay_potencies,
            preclinical_exposures=preclinical_exposures,
        )

        out_path = Path(args.out)
        write_csv(out_path, SCAFFOLD_CLUSTER_FIELDS, [cl.as_row() for cl in clusters])

        if args.compounds_out:
            c_fields = list(compounds_sar[0].as_row().keys()) if compounds_sar else []
            write_csv(Path(args.compounds_out), c_fields, [c.as_row() for c in compounds_sar])

        if args.report_path:
            report_text = render_sar_report(compounds_sar, clusters)
            write_text(Path(args.report_path), report_text)

        if args.json_path:
            summary = {
                "total_compounds": len(compounds_sar),
                "total_scaffolds": len(clusters),
                "scaffolds": [cl.as_row() for cl in clusters],
                "output_csv": str(out_path),
            }
            write_json(Path(args.json_path), summary)

    except SarError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: SAR analysis requires valid SMILES strings and consistent CSV data headers.", file=sys.stderr)
        print("Fix: verify compound dataset format and ensure SMILES strings are chemically well-formed.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: SAR output files could not be written to disk.", file=sys.stderr)
        print("Fix: specify a writable destination for --out, --report, or --json.", file=sys.stderr)
        return 1

    priority_count = sum(1 for c in clusters if c.sar_category == "PRIORITY_LEAD_SERIES")
    print(
        f"SAR analysis complete: {len(clusters)} scaffold(s) identified ({priority_count} priority series) "
        f"across {len(compounds_sar)} compound(s) -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
