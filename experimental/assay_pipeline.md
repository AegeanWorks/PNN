# Phase 3 In Vitro Assay Qualification Pipeline

## Purpose

`python -m pnn.assay` ingests experimental assay data, cross-evaluates PNN degradation potency against blood-brain barrier permeability and cytotoxicity counter-screens, and qualifies in silico leads into verified in vitro hits.

## Supported Assay Modalities

### 1. PNN Degradation Assay
- **Substrate**: DQ®-conjugated chondroitin sulfate (DQ-CS) or WFA lectin staining.
- **Metrics**: $EC_{50}$ (nM), maximum degradation percentage relative to ChABC positive control, Hill slope.
- **Qualification Threshold**: $EC_{50} \le 10,000\,\text{nM}$ ($10\,\mu\text{M}$) and Max Degradation $\ge 20\%$.

### 2. Blood-Brain Barrier (BBB) Permeability
- **PAMPA-BBB**: Artificial membrane permeability ($P_e \ge 3.0 \times 10^{-6}\text{ cm/s}$ or $\log BB \ge -1.0$).
- **MDCK-MDR1**: Cell monolayer transport ($P_{\text{app}} \ge 5.0 \times 10^{-6}\text{ cm/s}$ and Efflux Ratio $< 3.0$ with low P-gp liability).

### 3. Cytotoxicity Counter-Screen
- **CellTiter-Glo**: Viability in hypothalamic neuronal cells (`mHypoA-60/3`) or primary cultures.
- **Safety Window**: $CC_{50} \ge 50\,\mu\text{M}$ and Selectivity Index ($SI = CC_{50} / EC_{50} \ge 50$).

## Hit Qualification Statuses

| Status | Definition |
|---|---|
| **`CONFIRMED_HIT`** | Active PNN degradation + high BBB permeability + non-toxic ($SI \ge 50$ or $CC_{50} \ge 50\,\mu\text{M}$) |
| **`QUALIFIED_LEAD`** | Active PNN degradation + high BBB permeability, but cytotoxicity unmeasured or moderate |
| **`INCONCLUSIVE`** | Sub-threshold degradation, low permeability, or high MDCK efflux |
| **`CYTOTOXIC`** | Significant cell viability loss ($CC_{50} < 20\,\mu\text{M}$ or $SI < 10$) |
| **`INACTIVE`** | $EC_{50} > 50\,\mu\text{M}$ or maximum cleavage $< 10\%$ |

## CLI Usage

```bash
python -m pnn.assay \
    --degradation data/assays/degradation_template.csv \
    --permeability data/assays/permeability_template.csv \
    --cytotoxicity data/assays/cytotoxicity_template.csv \
    --out qualified_hits.csv \
    --json assay_summary.json \
    --report assay_report.md
```
