# IND-Enabling Safety Pharmacology & Repeat-Dose Toxicology Package

## Overview

`pnn.safety` turns submitted repeat-dose toxicology and ICH S7A/S7B safety pharmacology results
into a deterministic, auditable safety assessment. It establishes the most sensitive species
NOAEL, projects a human equivalent dose (HED) and NOAEL-based starting dose, and classifies each
compound as `SAFE_TO_PROGRESS`, `CAUTION` or `FAIL_SAFETY` before the IND dossier is assembled.

The module is not a GLP audit. It applies explicit progression thresholds to the data supplied;
missing measurements remain UNKNOWN and are never treated as passes (G-001).

## Input Contracts

### Repeat-dose toxicology CSV

```text
compound_id,species,study_type,duration_days,noael_mg_kg,loael_mg_kg,target_organ,findings
PNN-101,rat,REPEAT_DOSE_28D,28,50.0,100.0,liver,Reversible ALT elevation at LOAEL only
```

- Required columns: `compound_id`, `species`, `noael_mg_kg`.
- Optional columns: `study_type` (default `REPEAT_DOSE`), `duration_days` (default `28`),
  `loael_mg_kg`, `target_organ`, `findings`.
- `noael_mg_kg` must be a finite value greater than zero; `loael_mg_kg`, when present, must be
  greater than or equal to the NOAEL.
- Supported species and body-surface-area conversion factors (`K_m`): mouse 3, rat 6, hamster 5,
  guinea pig 7, rabbit 12, dog 20, monkey 12, marmoset 12, human 37.
- Multiple rows per compound are expected (one per species and study). The most sensitive
  (lowest) NOAEL drives the projection; ties resolve by species name then study type.

### Safety pharmacology and genotoxicity CSV

```text
compound_id,herg_ic50_um,free_cmax_um,cyp3a4_inhibition_pct,cyp2d6_inhibition_pct,ames_positive,micronucleus_positive,respiratory_noael_mg_kg,cns_noael_mg_kg
PNN-101,18.0,0.35,12.5,8.0,false,false,75.0,60.0
```

- Required column: `compound_id`; one row per compound (duplicates are rejected).
- Optional numeric columns: `herg_ic50_um` (> 0), `free_cmax_um` (> 0),
  `cyp3a4_inhibition_pct` and `cyp2d6_inhibition_pct` (0-100), `respiratory_noael_mg_kg` and
  `cns_noael_mg_kg` (>= 0).
- Optional boolean columns: `ames_positive`, `micronucleus_positive`. Accepted values are
  `true/false`, `yes/no`, `1/0`, `positive/negative`; blank means UNKNOWN.

## Assessment Rules

1. **NOAEL HED**: `NOAEL_HED (mg/kg) = NOAEL (mg/kg) * K_m,species / 37`.
2. **NOAEL-based starting dose**: `MRSD (mg) = (NOAEL_HED / safety_factor) * 60 kg`, with a
   default `safety_factor` of 10.
3. **hERG margin**: `hERG IC50 / free Cmax`. Fail below 10-fold; caution below the 30-fold target.
4. **CYP3A4 inhibition**: caution at 30% or more; fail at 50% or more.
5. **Genetic toxicology**: any positive Ames or micronucleus result fails; both negative is
   required for a clean package; anything else is UNKNOWN.
6. **Missing evidence**: no toxicology study, absent NOAEL, incomplete genotoxicity, or an
   incomputable hERG margin produce `CAUTION`, never `SAFE_TO_PROGRESS`.
7. **Status precedence**: any failure produces `FAIL_SAFETY`; otherwise any caution produces
   `CAUTION`; otherwise `SAFE_TO_PROGRESS`.

## Command-Line Usage

```bash
python -m pnn.safety \
  --tox data/toxicology/repeat_dose_tox_template.csv \
  --pharm data/toxicology/safety_pharmacology_template.csv \
  --out output/safety_assessments.csv \
  --json output/safety_summary.json \
  --report output/safety_package.md
```

The command writes the assessment CSV, an optional JSON summary and an optional Markdown package
with thresholds, a candidate matrix, findings and required follow-up. At least one of `--tox` or
`--pharm` must be supplied.

## Evidence Limits

- Study quality, GLP compliance, toxicokinetic coverage and pathology peer review are outside this
  deterministic check.
- The NOAEL-based starting dose is a projection. Reconcile it with the allometric projection from
  `pnn.ind` before selecting a first-in-human dose.
- A `SAFE_TO_PROGRESS` result means the submitted data clear these thresholds. It does not
  establish clinical safety, and no regulatory submission should describe it as such.
