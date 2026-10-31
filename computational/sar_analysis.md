# Structure-Activity Relationship (SAR) & Scaffold Decomposition

## Purpose

`python -m pnn.sar` provides standard-library Bemis-Murcko framework extraction, cycle system decomposition, R-group substituent enumeration, and activity-property clustering for small-molecule PNN degraders and PTPsigma inhibitors.

## Algorithmic Architecture

1. **SMILES Graph Parsing**:
   - Parses SMILES into an atom list and adjacency graph supporting organic subset atoms, brackets, aromatic flags, branching parenthesis stacks, and cycle ring closures (`1`-`9`, `%10`+).
2. **Bemis-Murcko Framework Decomposition**:
   - Utilizes the Hopcroft-Tarjan biconnected component algorithm to extract all cycle systems (blocks $\ge 3$ vertices).
   - Identifies non-cycle connected components; components bridging $\ge 2$ cycle blocks or atoms are categorized as **linkers**; single-attachment components are categorized as **substituents / R-groups**.
   - Reconstructs canonical framework SMILES representing the core ring systems and linkers.
3. **SAR Activity-Property Correlation**:
   - Groups compounds by shared core framework (`scaffold_id`).
   - Calculates cluster-level metrics: compound count, best analog, best and mean $EC_{50}$ (nM), best Vina affinity (kcal/mol), and mean brain-to-plasma partition ($K_{p,\text{brain}}$).
   - Assigns medicinal chemistry categories: `PRIORITY_LEAD_SERIES`, `VIABLE_OPTIMIZATION_SCAFFOLD`, `PERIPHERAL_RESTRICTED`, or `EXPLORATORY_SERIES`.

## CLI Usage

```bash
# Decompose library and generate scaffold clusters and Markdown report
python -m pnn.sar \
  --compounds data/compounds/vendor_library.csv \
  --docking data/docking/multi_target_scores.csv \
  --assays data/assays/degradation_template.csv \
  --preclinical data/preclinical/pk_timecourse_template.csv \
  --out sar_clusters.csv \
  --report sar_report.md
```
