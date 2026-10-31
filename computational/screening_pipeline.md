# Phase 2 Property Screening Pipeline

## Purpose

`pnn` is a standard-library-only triage step that applies the project's CNS/BBB
physicochemical rules and optional ADMET limits to a CSV of candidate compounds.
It answers one narrow question: **which compounds deserve docking and assay
resources first?** It does not predict permeability, potency or efficacy.

## Input Contract

The CSV is UTF-8, has one header row and one row per compound. Column names are
case-insensitive; spaces, hyphens and underscores are ignored, so `Molecular Weight`,
`molecular_weight` and `MW` all map to `mw`.

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `compound_id` | yes | - | unique identifier; duplicates are refused |
| `mw` | yes | Da | molecular weight |
| `logp` | yes | - | computed logP (state the calculator in the source data) |
| `tpsa` | yes | A2 | topological polar surface area |
| `hbd` | yes | count | hydrogen-bond donors |
| `hba` | yes | count | hydrogen-bond acceptors |
| `rotatable_bonds` | yes | count | rotatable bonds |
| `name` | no | - | display name |
| `smiles` | no | - | structure, retained for traceability |
| `pka` | no | - | ionization constant |
| `charge` | no | - | `neutral`, `basic`, `acidic` or `zwitterion` |
| `logbb` | no | - | measured or modelled logBB |
| `brain_plasma_ratio` | no | - | brain/plasma concentration ratio |
| `oral_bioavailability` | no | fraction or % | values <= 1 are read as fractions |
| `cyp3a4_inhibition` | no | % at 10 uM | CYP3A4 inhibition |
| `herg_liability` | no | % at 10 uM | hERG inhibition |
| `pnn_ec50_nm` | no | nM | PNN-softening EC50, when measured |

A missing optional value is reported as `UNKNOWN` and is excluded from the score;
it is never treated as a pass. An unrecognized `charge` is also `UNKNOWN`. Negative
values for counts, areas, ratios, percentages or EC50, and non-finite numbers
(`nan`, `inf`), are refused so a malformed row cannot pass silently. Use
`data/compounds/screening_template.csv` as the starting header.

## Rules

Core design properties (weight 1.0):

| Rule | Acceptable | Optimal |
|---|---|---|
| `mw` | 200-450 Da | 300-400 Da |
| `logp` | 1.5-3.5 | 2-3 |
| `tpsa` | <= 70 A2 | <= 60 A2 |
| `hbd` | <= 3 | <= 3 |
| `hba` | <= 7 | <= 7 |
| `rotatable_bonds` | <= 7 | <= 7 |

Optional properties (weight 0.5): pKa windows (basic 7-9, acidic 3-5), neutral or
basic charge preferred, logBB >= -0.5 preferred (below -1.0 fails), brain/plasma
>= 0.3 preferred (below 0.1 fails), oral F >= 20% preferred (below 10% fails),
CYP3A4 inhibition < 50% and hERG liability < 30% (at or above either fails), and
PNN EC50 < 100 nM preferred.

A rule contributes 1.0 inside its optimal window, a partial score inside the
acceptable window, and 0.0 outside it.

## Tiers

| Tier | Meaning |
|---|---|
| `A` | every core rule passes, score >= 0.80, at most two optional values unknown, and no hard ADMET failure |
| `B` | every core rule passes and score >= 0.60 |
| `C` | one or two core rules fail, or the score is below 0.60, with no disqualifying liability |
| `Reject` | any hard failure (logBB, brain/plasma, oral F, CYP3A4, hERG) or three or more core failures |

Ranking is tier first, then score (descending), then `compound_id`. A `B` tier
therefore still means "properties are plausible"; only measured data can move a
compound beyond triage.

## Usage

```bash
python -m pnn data/compounds/screening_template.csv --out ranked_compounds.csv
python -m pnn compounds.csv --out ranked.csv --json summary.json --min-tier B
```

- `--out` is the ranked CSV (default `ranked_compounds.csv`); parent directories are created.
- `--json` writes counts, emitted rows and the active rule sets.
- `--min-tier A|B|C` keeps that tier and better and drops `Reject`.
- Exit codes: `0` success, `1` unreadable/invalid input or output failure, `2` misuse.
  Failures print `What:`, `Why:` and `Fix:` on stderr.

Output columns mirror the input plus `tier`, `score`, `failed_criteria` and
`unknown_criteria`, so a reviewer can see why a compound was placed.

## Limits and Validation

- The rules are heuristic design criteria from `agent_design/bbb_penetration_rules.md`
  and `research_framework.md`; they are not a validated QSAR or BBB model.
- Input properties must come from a named calculator or measurement. Record that
  provenance with the source data; do not mix calculators in one comparison.
- Every candidate that advances must be confirmed by the Phase 3 MDCK-MDR1/PAMPA
  and PNN-degradation assays. Tier labels are prioritisation only.
- Rule changes belong in `pnn/rules.py` with a failing-first test in
  `tests/test_rules.py`; changing the input contract requires an ADR.
