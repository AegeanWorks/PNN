# Phase 5 Lead Optimization: Multiparametric Optimization (MPO) & Candidate Profiling

## Purpose

`python -m pnn.mpo` provides continuous Pfizer CNS MPO scoring (0.0 to 6.0 scale), Derringer-Suich geometric composite desirability index ($D \in [0.0, 1.0]$), and automated chemical liability diagnosis for candidate optimization.

## Scoring Methodology

### 1. Pfizer CNS MPO Desirability Functions
Evaluates 6 physicochemical parameters using continuous piecewise-linear desirability functions ($d_i \in [0.0, 1.0]$):
- **Molecular Weight ($MW$)**: Optimal $\le 360$ Da ($d=1.0$), monotonic decay to 0 at 500 Da.
- **Lipophilicity ($\log P$)**: Optimal $\le 3.0$ ($d=1.0$), monotonic decay to 0 at 5.0.
- **Distribution Coefficient ($\log D_{7.4}$)**: Optimal $\le 2.0$ ($d=1.0$), decay to 0 at 4.0.
- **Topological Polar Surface Area ($TPSA$)**: Optimal $40-90$ Å² ($d=1.0$), penalty for $<40$ or $>90$ Å².
- **Hydrogen Bond Donors ($HBD$)**: Optimal $\le 0.5$ ($d=1.0$), decay to 0 at 3.5.
- **Most Basic Center ($pK_a$)**: Optimal $\le 8.0$ ($d=1.0$), decay to 0 at 10.0.

### 2. Derringer-Suich Composite Desirability ($D$)
Weights multi-domain criteria across the discovery cascade:
$$D = (d_{\text{potency}} \times d_{\text{brain}} \times d_{\text{oral}} \times d_{\text{safety}})^{1/4}$$
- $d_{\text{potency}}$: DQ-CS degradation $EC_{50} \le 25$ nM ($d=1.0$) to $\ge 500$ nM ($d=0.0$).
- $d_{\text{brain}}$: Brain-to-plasma partition $K_{p,\text{brain}} \ge 0.50$ ($d=1.0$) to $\le 0.15$ ($d=0.0$).
- $d_{\text{oral}}$: Oral bioavailability $F \ge 30\%$ ($d=1.0$) to $\le 10\%$ ($d=0.0$).
- $d_{\text{safety}}$: Cytotoxicity selectivity index $SI \ge 50$ ($d=1.0$) to $\le 5$ ($d=0.0$).

### 3. Candidate Tiering
- **`DEVELOPMENT_CANDIDATE`**: Composite $D \ge 0.70$ and $\text{CNS MPO} \ge 4.5$.
- **`OPTIMIZATION_LEAD`**: Composite $D \ge 0.45$.
- **`DEPRIORITIZED`**: Composite $D < 0.45$.

## CLI Usage

```bash
python -m pnn.mpo \
  --compounds data/compounds/vendor_library.csv \
  --candidates data/preclinical/pk_timecourse_template.csv \
  --assays data/assays/degradation_template.csv \
  --out mpo_profiles.csv \
  --report mpo_report.md
```
