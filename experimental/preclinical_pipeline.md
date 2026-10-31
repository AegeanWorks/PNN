# Phase 4 In Vivo Pharmacokinetics, Pharmacodynamics & Preclinical Candidate Selection

## Purpose

`python -m pnn.preclinical` executes non-compartmental pharmacokinetic analysis (NCA), in vivo pharmacodynamic (PD) metabolic outcome evaluation, and preclinical candidate qualification for perineuronal net (PNN) disruption leads.

## Scientific Criteria & Decision Thresholds

| Preclinical Parameter | Clinical Candidate Threshold | Method / Biological Metric |
|---|---|---|
| **Brain Partition ($K_{p,\text{brain}}$)** | $\ge 0.30$ | Ratio of terminal $\text{AUC}_{\text{brain}} / \text{AUC}_{\text{plasma}}$ post-oral dosing |
| **Oral Bioavailability ($F$)** | $\ge 20.0\%$ | $(\text{AUC}_{\text{PO}} / \text{Dose}_{\text{PO}}) / (\text{AUC}_{\text{IV}} / \text{Dose}_{\text{IV}}) \times 100$ |
| **PNN Net Reduction** | $\ge 15.0\%$ | WFA+ area net reduction in mouse hypothalamic arcuate nucleus (ARC) |
| **Insulin Sensitivity** | $\ge 20.0\%$ | Insulin tolerance test (ITT) blood glucose clearance AUC improvement |

## Candidate Status Taxonomy

- **`CLINICAL_CANDIDATE`**: Meets or exceeds all brain partition ($K_p \ge 0.30$), oral bioavailability ($F \ge 20\%$), and in vivo metabolic efficacy thresholds ($\text{WFA} \ge 15\%$, $\text{ITT} \ge 20\%$).
- **`ADVANCED_LEAD`**: Sub-micromolar potency with robust brain exposure, pending full in vivo chronic profiling.
- **`INSUFFICIENT_BRAIN_PENETRATION`**: $K_{p,\text{brain}} < 0.20$ precludes central target engagement at safe systemic exposures.
- **`POOR_BIOAVAILABILITY`**: Oral $F < 15\%$ requires formulation enhancement or prodrug strategies.
- **`INSUFFICIENT_EFFICACY`**: Fails to demonstrate significant WFA+ net reduction ($<10\%$) in DIO mouse models.

## CLI Usage

```bash
python -m pnn.preclinical \
  --pk data/preclinical/pk_timecourse_template.csv \
  --pd data/preclinical/pd_metabolic_template.csv \
  --out preclinical_candidates.csv \
  --report preclinical_report.md
```
