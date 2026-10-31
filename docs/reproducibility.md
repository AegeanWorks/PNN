# Reproducibility Standards: PNN Project

## Purpose
This document defines the reproducibility standards for the PNN project, ensuring that all computational and experimental work can be replicated, verified, and built upon by any AI agent or human researcher, regardless of session or model.

---

## 1. Computational Reproducibility

### 1.1 Code Standards
| Requirement | Implementation |
|-------------|----------------|
| **Python version** | 3.9+ (specify in environment.yml) |
| **Package versions** | Pinned in requirements.txt or environment.yml |
| **Type hints** | All functions annotated |
| **Docstrings** | NumPy docstring format for all functions |
| **Unit tests** | pytest for critical functions (coverage ≥ 80%) |
| **Version control** | Git with descriptive commits |
| **Random seeds** | Set before any stochastic operation |

### 1.2 Required Environment
```yaml
# environment.yml (place in project root)
name: pnn-project
dependencies:
  - python=3.11
  - rdkit=2023.09.1
  - openmm=3.1.0
  - numpy=1.26.4
  - pandas=2.2.1
  - scipy=1.12.0
  - matplotlib=3.8.3
  - seaborn=0.13.2
  - jupyter=1.1.0
  - ipykernel=6.29.4
  - pip:
    - autoadock-vina=1.2.0
    - openmm-plumed=2.5.0  # if needed
    - mdanalysis=2.5.0
    - pandas-profiling=3.13.2  # optional
```

**Verification**: On cold start, run:
```bash
conda env create -f environment.yml
conda activate pnn-project
python -c "import rdkit; import openmm; print('Environment OK')"
```

### 1.3 Random Seed Protocol
**Set seed before any stochastic operation**:
```python
import numpy as np
import random

# At module level or script start
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

# For OpenMM (if used)
from openmm import app
# OpenMM uses its own random state management
```

**Document seed usage** in experiment reports:
```
Random seed: 42 (set at script start, used for: conformer generation, MD initialization, ML train/test split)
```

### 1.4 Computational Experiment Checklist

| Check | Required | Verification |
|-------|----------|--------------|
| **Environment file** | ✓ | environment.yml exists |
| **Package versions** | ✓ | `pip freeze > versions.txt` |
| **Random seed** | ✓ | Documented in report |
| **Input data** | ✓ | Paths absolute, not relative |
| **Output directory** | ✓ | Created if not exists, timestamped |
| **Intermediate files** | ✓ | Kept for audit trail |
| **Final plots** | ✓ | Saved as PNG/SVG with DPI, data in CSV |

### 1.5 Output Organization (Mandatory)
```
/output/
├── {experiment_name}/
│   ├── timestamp_YYYYMMDD_HHMMSS/  # Auto-generated
│   ├── input_data/                   # Copies of inputs used
│   ├── scripts/                      # Executed scripts
│   ├── results/                      # Processed results
│   │   ├── tables/
│   │   └── figures/
│   ├── logs/                         # Stdout/stderr captured
│   └── reproducibility/
│       ├── environment_snapshot.txt
│       ├── seed_log.txt
│       └── dependency_versions.txt
```

---

## 2. Experimental Reproducibility

### 2.1 Animal Studies

| Parameter | Standard Value | Documentation Requirement |
|-----------|---------------|--------------------------|
| **Species/Strain** | C57BL/6 or CD-1, male | Strain, source (Jackson, Charles River, etc.), age |
| **Age** | 8-12 weeks | Exact age at study start |
| **Weight range** | 20-35 g (mouse) | Pre-study weight recording |
| **Diet** | Chow (10% kcal fat) or HFD (60% kcal fat) | Diet type, manufacturer, duration |
| **Housing** | Single-housed or pair-housed | Cage type, enrichment, room conditions |
| **Light cycle** | 12:12 light-dark, lights on at 07:00 | ZT (zeitgeber time) specification |
| **Experimenter** | Named individual | Initials, training status |
| **Blinding** | Yes/no, who blinded | Group allocation, outcome assessment |
| **n per group** | Power-calculated (≥5) | Power analysis documented |
| **Ethics approval** | IACUC protocol # | Protocol number, approval date |

### 2.2 Cell Culture

| Parameter | Standard Value | Documentation Requirement |
|-----------|---------------|--------------------------|
| **Cell line** | mHypoA-60/3 or specified | RRID, passage number |
| **Passage number** | ≤ 20 from original | Documented at each experiment |
| **Culture medium** | DMEM/F-12 + 10% FBS + 1% pen-strep | Lot number if critical |
| **Coating** | Poly-D-lysine (10 μg/mL) | Coating concentration, duration |
| **Incubation** | 37°C, 5% CO2 | Exact conditions |
| **Serum batch** | Same lot for entire experiment | Lot number tracking |
| **Confluency at assay** | 70-80% | Microscopic confirmation |

### 2.3 Reagent Validation

| Reagent | Validation Requirement |
|---------|----------------------|
| **WFA-Biotin** | Lot-specific binding confirmed, store at -20°C |
| **ChABC** | Activity assay, units/mL confirmed |
| **Insulin (human)** | Batch potency, store at -20°C |
| **Antibodies** | Validated for species, application (IHC, WB, IF) |
| **Secondary antibodies** | Cross-validated, minimal cross-reactivity |

### 2.4 Assay Protocol Documentation

**Every assay must document**:
1. **Reagents** (catalog numbers, lot numbers, storage)
2. **Instrument settings** (plate reader parameters, microscope settings)
3. **Timing** (incubation durations, time points)
4. **Controls** (positive, negative, vehicle, blank)
5. **Replicates** (technical n, biological n)
6. **Analysis method** (software, version, parameters)
7. **Raw data** (preserved, not overwritten)
8. **Analysis script** (path, version, seed if applicable)

### 2.5 Data Integrity

| Requirement | Implementation |
|-------------|----------------|
| **Raw data preservation** | Never overwrite original files; `raw_data/` + `analysis/` |
| **File naming** | `{experiment}_{date}_{sample}_{replicate}.ext` |
| **Backup strategy** | Cloud + local, versioned |
| **Data dictionary** | Column headers documented in README of data folder |
| **Missing data** | Documented reason (equipment failure, outlier, etc.) |

---

## 3. Data Management Standards

### 3.1 Directory Structure (Mandatory)
```
PNN/
├── data/
│   ├── compounds/              # SDF/SMILES files
│   ├── assays/                 # Raw + analyzed assay data
│   ├── literature/             # Extracted references
│   └── raw_data/               # Original instrument outputs
├── docs/                       # All documentation
├── output/                     # Per-experiment outputs
├── results/                    # Aggregated results
└── sessions/                   # Session database (SQLite)
```

### 3.2 File Naming Convention
```
{project}_{date}_{phase}_{step}_{replicate}_{metric}.{ext}
```

**Examples**:
- `PNN_20261011_phase2_docking_run1_compound34_affinity.csv`
- `PNN_20261011_phase3_assayA_plate5_wellA12_raw.csv`
- `PNN_20261011_phase4_itt_mouse12_glucose.csv`
- `PNN_20261011_literature_review_brevican_2010.xlsx`

### 3.3 Data Dictionary (Per Dataset)
Every CSV/Excel file must have a companion `data_dictionary.md`:

```markdown
# Data Dictionary: PNN_20261011_phase3_assayA.csv

| Column | Type | Description | Units |
|--------|------|-------------|-------|
| compound_id | string | Unique compound identifier | - |
| concentration | float | Test compound concentration | μM |
| fluorescence | raw | Fluorescence intensity | AU |
| f_max | calculated | Max fluorescence (ChABC control) | AU |
| f_vehicle | raw | Vehicle fluorescence | AU |
| ec50 | calculated | Half-maximal effective concentration | μM |
| hbd | integer | H-bond donors count | - |
| tpsa | float | Topological polar surface area | Å² |
```

### 3.4 Session Database (SQLite)
**Use session_store_sql for cross-session state**:
```sql
-- Record experiment state
INSERT INTO session_state (key, value) VALUES 
  ('current_phase', '3'),
  ('last_compound_tested', 'Cmpd_042'),
  ('next_compound_index', '43'),
  ('pending_todos', '["Run assay on Cmpd_042-050"]');

-- Record checkpoint
INSERT INTO checkpoints (session_id, checkpoint_number, title, overview) VALUES
  ('current_session', 3, 'Phase 3 Complete', 'In vitro validation complete, 8 compounds advance to in vivo');
```

### 3.5 Backup Protocol
| Frequency | Scope |
|-----------|-------|
| **Per experiment** | output/, raw_data/, analysis scripts |
| **Per session** | /files/, session database, todos |
| **Per project** | Git repo (if initialized), /data/ archive |
| **Monthly** | Full project snapshot (compressed) |

---

## 4. Analysis Script Standards

### 4.1 Script Header (Required)
```python
#!/usr/bin/env python3
"""
Script: pnn_analyze_itt.py
Purpose: Analyze Insulin Tolerance Test data
Author: Yunus-Emre
Date: 2026-10-11
Version: 1.0

Dependencies:
  - Python 3.11
  - NumPy 1.26.4
  - Pandas 2.2.1
  - SciPy 1.12.0
  - Matplotlib 3.8.3

Seed: 42 (set at line 23)
"""

import numpy as np
import pandas as pd
from scipy import integrate
import matplotlib.pyplot as plt

# Set random seed for reproducibility
SEED = 42
np.random.seed(SEED)
```

### 4.2 Function Documentation
**All functions must have docstrings**:
```python
def calculate_auc(time: np.ndarray, conc: np.ndarray) -> float:
    """
    Calculate area under curve using trapezoidal rule.
    
    Parameters
    ----------
    time : np.ndarray
        Time points (minutes), monotonic increasing
    conc : np.ndarray
        Concentration values corresponding to time points
    
    Returns
    -------
    float
        AUC value (concentration × time units)
    
    Example
    -------
    >>> time = np.array([0, 15, 30, 60])
    >>> conc = np.array([180, 150, 120, 90])
    >>> calculate_auc(time, conc)
    5850.0
    """
    return integrate.trapz(conc, time)
```

### 4.3 Error Handling
```python
def safe_load_csv(path: str) -> pd.DataFrame:
    """Load CSV with error handling and logging."""
    try:
        df = pd.read_csv(path)
        return df
    except FileNotFoundError:
        print(f"ERROR: File not found: {path}")
        raise
    except pd.errors.EmptyDataError:
        print(f"ERROR: Empty file: {path}")
        raise
```

### 4.4 Output Directory Management
```python
import os
from datetime import datetime

def ensure_output_dir(experiment: str, phase: str) -> str:
    """Create timestamped output directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"output/{experiment}_{phase}_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)
    os.makedirs(f"{output_dir}/tables", exist_ok=True)
    return output_dir
```

---

## 5. Quality Control & Validation

### 5.1 QC Checklists by Phase

#### Phase 1: Literature Review
- [ ] All key papers extracted (n ≥ 20)
- [ ] Quantitative data in standardized table
- [ ] Human genetic associations documented
- [ ] Gaps identified and prioritized

#### Phase 2: Computational Screening
- [ ] Property filter applied (MW, logP, TPSA, etc.)
- [ ] Docking validation (redocking RMSD < 2.0 Å)
- [ ] ML model trained/tested (AUROC reported)
- [ ] Top 50 compounds selected with full property table

#### Phase 3: In Vitro Validation
- [ ] Assay positive controls validated
- [ ] Cytotoxicity CC50 determined for all compounds
- [ ] Triple screen (PNN + BBB + insulin) completed
- [ ] Top 5-10 ranked by MPO score

#### Phase 4: In Vivo Proof-of-Concept
- [ ] PK study complete with all parameters
- [ ] PD markers quantified (PNN, pIR, etc.)
- [ ] Metabolic outcomes (ITT, food intake) measured
- [ ] Go/No-Go decision documented

#### Phase 5: Lead Optimization
- [ ] SAR analysis complete
- [ ] Property optimization achieved
- [ ] Selectivity profile obtained
- [ ] IND-enabling package designed

### 5.2 Validation Metrics

| Metric | Acceptance Criterion |
|--------|---------------------|
| **Docking redocking RMSD** | < 2.0 Å (75% of test cases) |
| **ML model AUROC** | > 0.75 on external test set |
| **Assay Z'-factor** | > 0.5 (valid assay range) |
| **PK R² (terminal phase)** | > 0.95 |
| **ITT AUC n** | ≥ 6 mice per group |
| **WFA quantification** | Coefficient of variation < 20% |

---

## 6. Reporting Standards

### 6.1 Report Structure (Per Experiment)
```
1. Title & Date
2. Abstract (150-250 words)
3. Methods (detailed, reproducible)
4. Results (tables, figures)
5. Discussion (interpretation, limitations)
6. Conclusions
7. Data Availability
8. Code Availability
9. Reproducibility Statement
```

### 6.2 Reproducibility Statement Template
```
## Reproducibility

This experiment is fully reproducible as follows:

1. **Environment**: Conda environment `pnn-project` (environment.yml provided)
2. **Random seed**: 42 (set at script start, documented in analysis)
3. **Input data**: Copies preserved in `output/{experiment}_{timestamp}/input_data/`
4. **Code**: All analysis scripts in `scripts/` with version control
5. **Random number generation**: Seeded NumPy/Python random at module level
6. **Instrument settings**: Documented in methods (plate reader parameters, microscope settings)
7. **Controls**: Positive and negative controls included, results documented
8. **Data analysis**: Scripts available, deterministic given same inputs and seed
9. **Figures**: Generated from processed data, source data in `tables/` subdirectory
10. **Dependencies**: `pip freeze` output captured at experiment time
```

### 6.3 Required Figures/Tables
| Figure/Table | Minimum Content |
|-------------|-----------------|
| **Figure 1** | Study design overview |
| **Figure 2** | Primary results (main outcome measure) |
| **Figure 3** | Dose-response or time course |
| **Figure 4** | Controls and validation |
| **Figure 5** | Supplementary/secondary outcomes |
| **Table 1** | Participant/compound metadata |
| **Table 2** | Primary outcome statistics |
| **Table 3** | Adverse events / cytotoxicity |
| **Table 4** | Raw data summary |

---

## 7. Cold-Start Reproducibility Verification

### On New Session Startup:
```bash
# 1. Verify environment
conda env list | grep pnn
conda activate pnn-project

# 2. Check environment.yml exists
ls -la environment.yml

# 3. Verify key scripts exist
ls -la computational/ scripts/ docs/

# 4. Check session database state
python -c "from sql import get_session_state; print(get_session_state())"

# 5. Run reproducibility check
python -c "
import sys
errors = []

# Check environment
try:
    import rdkit
    import openmm
except ImportError as e:
    errors.append(f'Missing dependency: {e}')

# Check seed usage
try:
    with open('scripts/last_analysis.py') as f:
        content = f.read()
        if 'SEED' not in content:
            errors.append('No random seed found in analysis scripts')
except FileNotFoundError:
    errors.append('No scripts directory found')

if errors:
    print('REPRODUCIBILITY ISSUES:')
    for e in errors:
        print(f'  - {e}')
else:
    print('✓ Reproducibility check passed')
"
```

### 7.1 Common Reproducibility Failures & Fixes

| Failure | Cause | Fix |
|---------|-------|-----|
| **Different results each run** | Random seed not set | Add `np.random.seed(42)` at script start |
| **Import errors** | Environment not activated | `conda activate pnn-project` |
| **Relative path failures** | Working directory changed | Use absolute paths or `os.chdir()` at start |
| **Missing data** | Raw data not backed up | Copy raw_data/ before analysis, preserve originals |
| **Plot differs** | Matplotlib backend/version | Specify matplotlib version, save with explicit DPI |
| **AUC calculation varies** | Integration method | Use `scipy.integrate.trapz`, document method |

---

## 8. Git Repository (Optional but Recommended)

### 8.1 Initialization
```bash
cd /Users/Yunus-Emre/Desktop/PNN
git init
git remote add origin <repository-url>
```

### 8.2 Commit Convention
```
[type]: concise description

Types:
- feat: new feature or protocol
- fix: bug fix or error resolution
- docs: documentation update
- refactor: code restructuring without behavior change
- test: adding or fixing tests
- chore: housekeeping, configuration

Examples:
- feat: add PNN degradation assay protocol
- fix: correct logD7.4 calculation in filter script
- docs: update cold_start_guide.md with Phase 3 protocol
- refactor: reorganize computational/docking/ directory
- test: add unit tests for property_filter.py
- chore: update environment.yml with openmm 3.1.0
```

### 8.3 .gitignore (Essential)
```
# Dependencies
.env
environment.yml
**/*.pyc
__pycache__/
*.pkl
*.h5
*.hdf5

# Output
output/
**/timestamp_*/

# OS
.DS_Store
Thumbs.db

# Large data (use git-lfs or /files/)
*.sdf
*.csv
*.xlsx
*.mdcrd
*.mdcrd

# Session-specific
session_state.db
*.log
```

---

## 9. Version Control for Documents

### 9.1 Document Versioning
All markdown files follow: `v{major}.{minor}`

**Version bump triggers**:
- Major framework changes → v2.0.0
- New phase addition → v1.1.0
- Protocol overhaul → v1.2.0
- Typos/fixed → v1.0.1

**Version header** (add to top of each markdown file):
```markdown
# Document Version: 1.0
# Generated: 2026-10-11
# Compatible with: Recursive Language Model cold-start workflows
```

### 9.2 Checkpoint Versioning
```
checkpoint_phase1.md → v1.0 (initial)
checkpoint_phase1.md → v1.1 (if protocol modified)
checkpoint_phase2.md → v2.0 (new phase)
```

---

## 10. Ethical & Compliance Standards

### 10.1 Animal Research Compliance
- **IACUC protocol**: Must be approved before any in vivo work
- **Humane endpoints**: Define weight loss > 20%, morbidity, severe distress
- **Statistical rigor**: Power analysis before study, α = 0.05, β = 0.2 (80% power)
- **Data reporting**: Complete disclosure of all outcomes (positive and negative)

### 10.2 Computational Ethics
- **No human data** without IRB approval (if applicable)
- **Compound library sourcing**: Document origin (ChEMBL, ZINC, de novo)
- **ML model transparency**: Document training data, features, limitations
- **Open science**: Share protocols, code, data where possible

### 10.3 Data Privacy
- **No PHI** (protected health information) in session data
- **Anonymized identifiers** only for animal subjects (mouse_01, mouse_02, etc.)
- **Compound identifiers**: Use systematic naming, not personal data

---

## 11. Appendix: Standard Operating Procedures Quick Reference

### SOP: Property Filter (CNS + Oral)
```bash
# Run property filter on compound library
python scripts/property_filter.py \
    --input data/compounds/library.sdf \
    --output data/compounds/filtered.sdf \
    --log results/property_filter.log
```

### SOP: Docking Run (Vina)
```bash
# Run Vina docking on filtered library
for f in data/compounds/filtered_*.sdf; do
    vina --receptor computational/targets/P5K4P_prepared.pdb \
         --ligand "$f" \
         --out output/docking/$(basename $f .sdf)/output.pdbqt \
         --log output/docking/$(basename $f .sdf)/output.log \
         --exhaustiveness 8 \
         --center_x 12.5 --center_y 18.3 --center_z 22.1 \
         --size_x 20 --size_y 20 --size_z 20
done
```

### SOP: PK Parameter Calculation
```bash
# Run PK analysis Python script
python scripts/analyze_pk.py \
    --input data/assays/in_vivo/mouse_pk.csv \
    --output results/pk_parameters.xlsx \
    --seed 42
```

### SOP: WFA Quantification
```bash
# Run ImageJ macro analysis
python scripts/analyze_wfa.py \
    --input data/raw_data/wfa/*.tif \
    --output results/wfa_quantification.xlsx \
    --method puncta_density
```

---

*Document Version: 1.0*
*Generated: 2026-10-11*
*Compatible with Recursive Language Model cold-start workflows*
*Project: Oral BBB-Penetrant Agents for Hypothalamic PNN Softening*