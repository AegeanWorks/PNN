# Molecular Dynamics Trajectory & Contact Stability Analysis

## Overview

`pnn.dynamics` evaluates nanosecond-scale kinetic stability of docked ligand poses in perineuronal
net receptor pockets (for example the PTPσ CSPG binding pocket). Static docking scores do not
guarantee residence time; this module quantifies pose drift, pocket breathing and contact
persistence across trajectory frames.

The analyzer only reports what the input supplies. Optional columns that are absent stay
`UNMEASURED`; they are never replaced with defaults.

## Input Contract

Required columns (one alias each):

- `time_ns` (`time`, `frame_ns`)
- `ligand_rmsd_angstrom` (`ligand_rmsd`, `rmsd_ligand`, `rmsd_lig`)

Optional columns:

- `protein_backbone_rmsd` (`protein_rmsd`, `rmsd_protein`, `rmsd_rec`)
- `key_contacts_count` (`contacts_count`, `contacts`, `num_contacts`)
- `mmgbsa_binding_dg` (`mmgbsa_dg`, `binding_dg`, `dg_bind`)
- Residue contact flags `contact_<RES>` or `res_<RES>` (`1/true/yes` means active).

At least two frames are required.

## Key Trajectory Metrics

1. **Ligand RMSD**: mean, standard deviation and maximum; stable threshold mean < 2.0 Å and max
   < 3.2 Å.
2. **RMSD drift**: mean of the final 20% of frames minus the mean of the initial 20%.
3. **Protein backbone RMSD**: optional; confirms the receptor fold remains stable (< 2.5 Å).
4. **Pocket contact persistence**: optional per-residue percentage of frames with the contact active.
5. **MM-GBSA free energy**: optional; strong binding is ≤ −25.0 kcal/mol.

## Pose Stability Classification

Contact counts are required for a verdict. Without them the class is `INSUFFICIENT_EVIDENCE`, the
stability score is `UNMEASURED`, and the report says so.

- **`STABLE_POSE`**: mean ligand RMSD < 2.0 Å, drift < 1.0 Å, mean contacts ≥ 2.0, max RMSD < 3.2 Å.
- **`MODERATE_MOBILITY`**: mean ligand RMSD ≤ 3.5 Å, drift < 2.0 Å, mean contacts ≥ 1.0.
- **`UNSTABLE_DRIFT`**: mean ligand RMSD > 3.5 Å, drift ≥ 2.0 Å, or mean contacts < 1.0.

Evidence status is `COMPLETE` when protein RMSD, contacts and MM-GBSA are all present, `PARTIAL`
when the pose verdict is computable but one of those is absent, and `INSUFFICIENT_EVIDENCE` when the
verdict itself cannot be reached.

## Command-Line Usage

```bash
python -m pnn.dynamics \
  --trajectory data/dynamics/ptpsigma_md_trajectory_template.csv \
  --compound PNN-101 \
  --target PTPsigma \
  --out-dir output/dynamics
```

The command always writes:

- `{compound}_trajectory_metrics.csv`: summary statistics, stability class, score and evidence status.
- `{compound}_dynamics_summary.json`: JSON for programmatic consumption.
- `{compound}_dynamics_report.md`: Markdown report with `UNMEASURED` markers and residue persistence.
- `{compound}_rmsd_trajectory.svg`: RMSD time series (protein trace only when supplied).
