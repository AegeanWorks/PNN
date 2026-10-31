# Human Allometric Dose Projection & IND Regulatory Candidate Dossier

## Overview

`pnn.ind` translates measured preclinical dose, clearance, volume and NOAEL data into a human
equivalent dose (HED), a Phase 1 Maximum Recommended Starting Dose (MRSD), projected human PK and
an IND nomination dossier.

Every value is sourced from an input column. A missing measurement is reported as `UNMEASURED` and
the compound's `evidence_status` drops to `PARTIAL` or `INSUFFICIENT_EVIDENCE`. The module never
substitutes a default Kp, F, dose, clearance or NOAEL.

## Inputs

`--compounds` (required, one or more CSVs) supplies `compound_id` and `name`.

`--candidates` (optional) supplies measured preclinical values. Recognised columns:

| Column | Used as |
|---|---|
| `dose_mg_kg` | Efficacious animal dose |
| `cl_l_h_kg` | Mouse clearance for allometric scaling |
| `vd_l_kg` | Mouse volume of distribution |
| `brain_plasma_ratio` | Brain partition `Kp` |
| `oral_bioavailability_pct` | Oral bioavailability `F` |
| `noael_mg_kg` | NOAEL (when no `--tox` is supplied) |

`--tox` (optional) supplies repeat-dose toxicology via `pnn.safety`; the lowest NOAEL per compound
is selected and its species `K_m` is used for the NOAEL HED.

`--lead` nominates a compound by ID; otherwise the first compound with measured data is used. The
nominated lead is listed first and drives the dossier summary.

## Equations

1. **Human Equivalent Dose (HED)**:
   $$\text{HED}\ (\text{mg/kg}) = \text{Animal Dose}\ (\text{mg/kg}) \times \left( \frac{3}{37} \right)$$
2. **Total Human Daily Dose**:
   $$\text{Daily Dose}\ (\text{mg}) = \text{HED}\ (\text{mg/kg}) \times 60.0\text{ kg}$$
3. **Interspecies Scaling**:
   $$CL_{\text{human}} = CL_{\text{mouse}} \times W_{\text{mouse}} \times \left( \frac{60.0}{0.02} \right)^{0.75}$$
   $$V_{d,\text{human}} = V_{d,\text{mouse}} \times W_{\text{mouse}} \times \left( \frac{60.0}{0.02} \right)^{1.0}$$
4. **Projected Half-Life**:
   $$t_{1/2,\text{human}} = \frac{\ln(2) \times V_{d,\text{human}}}{CL_{\text{human}}}$$
5. **NOAEL HED and MRSD**:
   $$\text{NOAEL}_{\text{HED}} = \text{NOAEL} \times \frac{K_{m,\text{species}}}{37}$$
   $$\text{MRSD}\ (\text{mg}) = \frac{\text{NOAEL}_{\text{HED}}}{10} \times 60.0\text{ kg}$$

The 10-fold factor is the default ICH M3(R2) safety divisor and is recorded in the dossier limits.

## Evidence Status

- `COMPLETE`: dose, mouse CL/Vd, NOAEL, brain Kp and oral F were all supplied.
- `PARTIAL`: at least one measurement is present; missing values remain blank.
- `INSUFFICIENT_EVIDENCE`: nothing measurable was supplied.

The dossier lists every gap and its owning module (`pnn.preclinical` for dose/CL/Vd/Kp/F,
`pnn.safety` for the most sensitive NOAEL).

## Command-Line Usage

```bash
python -m pnn.ind \
  --compounds data/compounds/vendor_library.csv data/compounds/repurposing_library.csv \
  --candidates output/preclinical_candidates.csv \
  --tox data/toxicology/repeat_dose_tox_template.csv \
  --lead PNN-REP-01 \
  --out output/human_dose_projections.csv \
  --dossier output/ind_dossier.md \
  --json output/ind_summary.json
```

Outputs:

- `human_dose_projections.csv`: one row per compound, including `evidence_status`.
- `ind_dossier.md`: nomination summary, scaling matrix, evidence gaps, unchecked ICH S7A/S7B
  checklist and Phase 1 protocol skeleton.
- `ind_summary.json` (optional): counts by evidence status and the nominated lead row.

## Reconciliation with `pnn.safety`

`pnn.safety` selects the most sensitive species NOAEL and computes its own NOAEL-based MRSD. When
`--tox` is supplied, `pnn.ind` uses the same lowest NOAEL and species factor, so the two MRSD values
should match. Any difference means the two commands were given different toxicology inputs.

## Limits

- Allometric scaling is a projection from animal data and does not establish human safety or efficacy.
- Study quality, GLP compliance, exposure margins and pathology review belong to a qualified
  toxicologist; the dossier's safety checklist items remain unchecked until that evidence exists.
