# Phase 2 Lead Selection Pipeline

## Purpose

`python -m pnn.pipeline` executes the end-to-end Phase 2 computational triage workflow in a single reproducible run. It combines library curation, CNS/BBB property screening, and docking consensus ranking, emitting standardized tabular datasets, a Markdown lead report, and a cryptographically hashed manifest.

This pipeline answers the core Phase 2 computational question: **Which candidate compounds exhibit the combined physicochemical profile and multi-target docking consensus required to justify progression to in vitro permeability and PNN-softening assays?**

## Pipeline Stages

```
Vendor / Repurposing Compound CSVs
                 │
                 ▼
       ┌───────────────────┐
       │ 1. Curation &     │   Deduplicate by ID & exact SMILES;
       │    Diversity      │   apply screening tier floor;
       │    (pnn.curate)   │   cap 6D property-space bins.
       └─────────┬─────────┘
                 │ -> curated_library.csv
                 ▼
       ┌───────────────────┐
       │ 2. CNS/BBB        │   Evaluate core & ADMET properties;
       │    Screening      │   assign Tiers A, B, C or Reject;
       │    (pnn.screen)   │   calculate composite property score.
       └─────────┬─────────┘
                 │ -> ranked_compounds.csv
                 ▼
       ┌───────────────────┐   Docking Scores CSV
       │ 3. Consensus      │   (Vina, Glide, MM-GBSA)
       │    Ranking        │         │
       │    (pnn.rank)     │◄────────┘
       └─────────┬─────────┘   Min-max normalize docking scales;
                 │             filter by pose RMSD & Vina affinity;
                 │             join screening tier.
                 │ -> shortlist.csv
                 │ -> lead_report.md
                 ▼ -> manifest.json
```

### Stage 1: Library Curation & Diversity
- Consumes one or more compound CSVs (e.g. vendor or repurposing libraries).
- Ingests compound structures and physicochemical properties.
- Deduplicates compounds by `compound_id` and exact non-empty SMILES string.
- Filters candidates by screening tier floor (`--min-tier`, default `B`).
- Discretizes the 6D core property space (`mw`, `logp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`) into normalized grid cells and caps each bin at `--per-bin` entries (default 25) to prevent property-space clustering.
- Emits `curated_library.csv`.

### Stage 2: CNS/BBB Property Screening
- Evaluates curated compounds against the CNS/BBB design windows documented in [screening_pipeline.md](screening_pipeline.md).
- Flags violations of core and optional ADMET criteria.
- Emits `ranked_compounds.csv` sorted by tier and composite property score.

### Stage 3: Docking Consensus Ranking & Shortlist
- Ingests multi-target docking scores (AutoDock Vina, Glide, MM-GBSA).
- Normalizes disparate scoring functions into a 0..1 scale relative to the run's docking set.
- Applies pose RMSD (`--max-rmsd`) and affinity (`--min-vina`) thresholds.
- Selects the optimal target interaction pose per compound, weighted across available methods.
- Emits `shortlist.csv`.

### Stage 4: Reporting & Cryptographic Manifest
- Renders `lead_report.md`, summarizing input provenance, curation metrics, tier distributions, top shortlisted candidates, and experimental evidence gaps.
- Compiles `manifest.json` under schema `pnn.phase2.pipeline/1`, recording SHA-256 digests for all input files, execution options, stage counts, and output artifacts.

## Usage

```bash
python -m pnn.pipeline vendor_library.csv repurposing_library.csv \
    --docking docking_results.csv \
    --out-dir phase2_lead_run \
    --min-tier B \
    --per-bin 25 \
    --max-rmsd 2.0
```

### Key Options

| Option | Default | Description |
|---|---|---|
| `inputs` | *(required)* | One or more compound library CSVs |
| `--docking` | *(required)* | Docking scores CSV with Vina, Glide, and/or MM-GBSA columns |
| `--out-dir` | `phase2_run` | Destination directory for generated artifacts |
| `--min-tier` | `B` | Minimum screening tier required for progression (`A`, `B`, or `C`) |
| `--per-bin` | `25` | Maximum compound capacity per 6D property-space bin |
| `--bin-width` | `0.1` | Normalized bin width (0..1) along each property axis |
| `--max-compounds` | `None` | Optional global ceiling on curated candidates |
| `--max-rmsd` | `None` | Optional pose RMSD cutoff in Angstroms |
| `--min-vina` | `None` | Minimum acceptable Vina affinity in kcal/mol |
| `--top-n-report` | `25` | Maximum top leads detailed in `lead_report.md` |

## Generated Artifacts

| Artifact | Format | Description |
|---|---|---|
| `curated_library.csv` | CSV | Deduplicated and property-bin-capped compound library with provenance |
| `ranked_compounds.csv` | CSV | Full screening results and property scores for all curated candidates |
| `shortlist.csv` | CSV | Prioritised leads with target consensus scores, coverage, and pose metrics |
| `lead_report.md` | Markdown | Comprehensive run report with G-001/G-002/G-003 triage disclaimers |
| `manifest.json` | JSON | Cryptographic SHA-256 run manifest under schema `pnn.phase2.pipeline/1` |

## Evidence Limits & Project Guidelines

This pipeline enforces the project's governance guidelines:
- **G-001 (ADR-001)**: Screening tiers are triage indicators only. Never represent a compound as BBB-penetrant or active in the absence of measured permeability (PAMPA-BBB / MDCK-MDR1) and PNN digestion assays.
- **G-002 (ADR-002)**: Docking consensus scores are relative to the docked candidate set and normalize across available scoring functions; they do not substitute for measured binding constants ($K_d$, $K_i$).
- **G-003 (ADR-003)**: Curated libraries preserve input source provenance. Exact-SMILES deduplication removes identical strings but is not chemical canonicalization.
