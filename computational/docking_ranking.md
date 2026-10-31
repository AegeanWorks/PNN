# Docking Consensus Ranking

## Purpose

Docking engines report scores on different scales. `python -m pnn.rank` turns
per-target Vina, Glide and MM-GBSA results into one reproducible shortlist without
pretending that raw scores from different engines are directly comparable. It is a
ranking tool, not a binding-free-energy calculation.

## Docking Input

UTF-8 CSV, one row per compound/target pair. Column names are case-insensitive and
parenthetical units are ignored.

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `compound_id` | yes | - | joins to the screening output |
| `target` | yes | - | e.g. `PTPsigma`, `LAR`, `NgR1` |
| `vina_affinity` | one score required | kcal/mol | AutoDock Vina affinity; lower is better |
| `glide_score` | one score required | - | Glide SP/XP score; lower is better |
| `mmgbsa_dg` | one score required | kcal/mol | MM-GBSA binding free energy; lower is better |
| `rmsd` | no | A | ligand RMSD for the chosen pose |
| `hbond_count` | no | count | hydrogen bonds in the chosen pose |
| `key_interactions` | no | - | residue interaction summary |

Aliases such as `Vina`, `Glide`, `MM-GBSA dG`, `RMSD (A)` and `H-bond count` are
accepted. Duplicate `compound_id`/`target` pairs and non-finite values are refused.
At least one score must be present, but **only rows with `vina_affinity` can enter
the consensus**; this keeps the primary engine consistent across the shortlist.

Start from `data/docking/scores_template.csv`.

## Normalization and Consensus

Within the filtered docking set, each method is min-max scaled to a quality score:

```text
quality = (max - value) / (max - min)      # 1.0 is the best value in this set
```

The consensus is the weighted average of the methods available on that row:

```text
consensus = 0.5 * vina + 0.3 * glide + 0.2 * mmgbsa
```

Weights are renormalized when a method is missing, so a Vina-only row is scored on
Vina alone. `coverage` reports the fraction of the full weight present: `0.5` for
Vina only, `0.8` for Vina + Glide, `1.0` for all three. Normalization is
set-dependent: a score is only meaningful relative to the compounds in that run.

For each compound, the target row with the highest consensus wins; ties prefer
higher coverage and then the stronger Vina affinity. `targets_docked` counts how
many targets the compound reached.

## Optional Screening Join

`--screen` joins the ranked property CSV from `python -m pnn`. The shortlist then
carries `name`, `smiles`, `tier` and `property_score`, and `--min-tier` (default
`B`) drops compounds whose screening tier is worse. Compounds absent from the
screening table are skipped. Without `--screen`, every docked compound is ranked
and the tier field stays empty.

## Usage

```bash
python -m pnn.rank data/docking/scores_template.csv \
    --screen ranked_compounds.csv \
    --out shortlist.csv --json shortlist_summary.json

python -m pnn.rank docking.csv --out shortlist.csv --max-rmsd 2.0 --min-vina -7.0
```

- `--max-rmsd` drops target rows above the pose-RMSD limit; rows without an RMSD are kept.
- `--min-vina` drops target rows with Vina affinity weaker than the cutoff.
- `--min-tier A|B|C` keeps that screening tier and better; `Reject` is never emitted.
- `--out` defaults to `shortlist.csv`; parent directories are created.
- Exit codes: `0` success, `1` invalid input/filters/output, `2` misuse.

## Output Columns

`rank`, `compound_id`, `name`, `smiles`, `tier`, `property_score`, `target`,
`consensus_score`, `coverage`, `vina_affinity`, `glide_score`, `mmgbsa_dg`, `rmsd`,
`hbond_count`, `key_interactions`, `targets_docked`.

## Limits

- Scores are normalized against the current set, so adding or removing compounds
  changes the ranking. Keep the input set fixed when comparing runs. A set with a
  single row for a method normalizes that method to `1.0`, which is not evidence
  of quality.
- A favorable consensus is prioritisation evidence only. It does not establish
  binding, PNN softening, BBB penetration or safety.
- Promote compounds to re-docking, MD and MM-GBSA, then to the Phase 3 assays; do
  not report a docking score as an experimental result.
- The ranker assumes lower is better for every score. Convert engines that report
  the opposite convention before export.
