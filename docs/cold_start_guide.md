# Cold-Start Guide: PNN Project Initialization for AI Agents

## Purpose
This guide ensures any AI agent (regardless of model, session, or cold start) can immediately begin productive work on the PNN project with full context and standardized procedures.

---

## 1. Project Identity & Mission

**Project**: PNN - Oral BBB-Penetrant Agents for Hypothalamic PNN Softening
**Ultimate Goal**: Fix insulin signaling in the hypothalamus to resolve obesity
**Approach**: Develop oral, blood-brain barrier (BBB) penetrant agents that soften perineuronal nets (PNNs) in the arcuate nucleus (ARC)

**Core Hypothesis**: Hypothalamic PNNs (primarily CSPGs: brevican, neurocan, versican) create a physical and signaling barrier that:
1. Impairs insulin receptor accessibility on POMC/NPY neurons
2. Disrupts insulin signaling cascades (IRS-PI3K-AKT pathway)
3. Contributes to leptin/insulin resistance in diet-induced obesity
4. Softening PNNs restores insulin sensitivity and normalizes energy homeostasis

---

## 2. Directory Structure (Must Exist)

```
PNN/
├── README.md                    # Project overview (READ FIRST)
├── literature_review.md         # Comprehensive literature findings
├── research_framework.md        # Methodology & experimental design
├── agent_design/                # Agent design specifications
│   ├── target_profiles.md       # Target protein profiles
│   ├── bbb_penetration_rules.md # BBB penetration guidelines
│   └── oral_bioavailability.md  # Oral drug design principles
├── computational/               # Computational approaches
│   ├── docking_protocols.md     # Molecular docking procedures
│   ├── md_simulations.md        # Molecular dynamics protocols
│   └── ml_models.md             # ML model specifications
├── experimental/                # Experimental validation
│   ├── in_vitro_assays.md       # In vitro assay protocols
│   ├── in_vivo_models.md        # Animal model specifications
│   └── biomarkers.md            # Biomarker panels
├── data/                        # Data storage
│   ├── compounds/               # Compound libraries (SDF/SMILES)
│   ├── assays/                  # Assay results
│   └── literature/              # Extracted literature data
└── docs/                        # Documentation
    ├── cold_start_guide.md      # THIS FILE
    └── reproducibility.md       # Reproducibility standards
```

**Verification**: On startup, run `ls -la` in project root to confirm structure exists.

---

## 3. Cold-Start Sequence (Execute in Order)

### Step 1: Load Core Context (3 files, ~5 min)
```bash
# Read in this exact order:
1. README.md                    # Project overview, structure, TPP
2. literature_review.md         # Scientific foundation, targets, gaps
3. research_framework.md        # Phases, milestones, protocols, go criteria
```

### Step 2: Initialize Session State
```bash
# 1. Set working directory
cd /Users/Yunus-Emre/Desktop/PNN  # or your project path

# 2. Check project structure
ls -la

# 3. Initialize session todos from framework
# (Use session SQL database if available)

# 4. Load agent design specs (if doing computational work)
cat agent_design/target_profiles.md
cat agent_design/bbb_penetration_rules.md
cat agent_design/oral_bioavailability.md
```

### Step 3: Determine Current Phase & Next Actions

**Phase Detection** (check latest checkpoint or todo status):
- **Phase 1** (Months 1-3): Target validation, literature extraction
- **Phase 2** (Months 3-6): In silico screening, compound library design
- **Phase 3** (Months 6-12): In vitro validation, assay development
- **Phase 4** (Months 12-18): In vivo proof-of-concept
- **Phase 5** (Months 18-24): Lead optimization, candidate selection

**Action Selection** based on phase:
| Phase | Primary Tasks |
|-------|---------------|
| **1** | Extract quantitative data from literature, validate PNN-insulin link |
| **2** | Build compound library, run docking, predict BBB penetration |
| **3** | Develop/validate assays, test top 50 compounds |
| **4** | Design PK/PD studies, execute in vivo experiments |
| **5** | SAR analysis, lead optimization, IND-enabling package |

---

## 4. Phase-Specific Quick-Start Protocols

### Phase 1: Literature & Target Validation
```bash
# 1. Read literature_review.md completely
# 2. Identify 20+ key papers for deep reading
# 3. Extract quantitative data tables:
#    - PNN density in lean vs DIO (WFA+ puncta/neuron)
#    - Insulin signaling metrics (pIR, pAKT fold-changes)
#    - ChABC efficacy data (dose, duration, outcomes)
#    - Human genetic associations (BCAN, NCAN loci)
# 4. Create data/literature/quantitative_data.xlsx
# 5. Update research_framework.md with validated numbers
```

### Phase 2: Computational Screening
```bash
# 1. Read computational/docking_protocols.md
# 2. Set up Python environment:
#    python -c "import rdkit; print('RDKit OK')"
#    vina --help  # Check AutoDock Vina
# 3. Download target structures:
#    wget https://files.rcsb.org/download/5K4P.pdb  # PTPσ + CSPG
#    wget https://files.rcsb.org/download/2PTS.pdb  # PTPσ catalytic
#    wget https://files.rcsb.org/download/2LAR.pdb  # LAR D1
# 4. Prepare receptors (see docking_protocols.md Section 1)
# 5. Build/load compound library (data/compounds/)
# 6. Run property filter (MW 300-450, logP 1.5-4, TPSA < 90)
# 7. Execute Vina docking (exhaustiveness 8)
# 8. Rank by affinity + BBB prediction + ML score
# 9. Select top 50 for MD simulation
```

### Phase 3: In Vitro Validation
```bash
# 1. Read experimental/experimental_protocols.md
# 2. Validate PNN degradation assay with ChABC positive control
# 3. Validate BBB permeability assay with known standards
# 4. Validate insulin signaling assay with insulin positive control
# 5. Test top 50 compounds in triple screen:
#    - PNN degradation (EC50)
#    - BBB permeability (Pe, ER)
#    - Insulin signaling rescue (pIR, pAKT EC50)
# 6. Cytotoxicity counter-screen (CC50)
# 7. Rank by MPO score (see experimental_protocols.md Section 7)
# 8. Select top 5-10 for in vivo
```

### Phase 4: In Vivo Proof-of-Concept
```bash
# 1. Read experimental/experimental_protocols.md (in vivo sections)
# 2. Design PK study: 4 doses (1, 3, 10, 30 mg/kg), 6 time points
# 3. Execute oral gavage PK in mice (n=3-4/time point)
# 4. Analyze brain/plasma ratios
# 5. If PK acceptable (F > 15%, brain/plasma > 0.3):
#    - Chronic dosing study (2-4 weeks)
#    - ITT, food intake, body weight
#    - Terminal WFA quantification
# 6. Go/No-Go decision (see research_framework.md Section 6)
```

### Phase 5: Lead Optimization
```bash
# 1. SAR analysis on active series
# 2. Property optimization (logD, solubility, metabolic stability)
# 3. Selectivity profiling (Cerep panel, hERG, CYP)
# 4. Scale-up synthesis of lead candidate
# 5. IND-enabling toxicology package design
# 6. Final candidate nomination
```

---

## 5. Standardized Output Formats

### 5.1 Results Documentation
All results saved as Markdown with tables:
```markdown
## Experiment: [Name]
**Date**: YYYY-MM-DD
**Phase**: X
**Compounds tested**: N
**Key findings**:
| Compound | Metric 1 | Metric 2 | Metric 3 | Pass/Fail |
|----------|----------|----------|----------|-----------|
| Cmpd_001 | 5.2 μM   | 12.5     | 8.7      | ✓         |
```

### 5.2 Data Files
| Data Type | Format | Location |
|-----------|--------|----------|
| Compound library | SDF / SMILES CSV | data/compounds/ |
| Docking results | CSV / Excel | computational/docking/results/ |
| MD trajectories | NetCDF / DCD | computational/md_simulations/ |
| Assay data | Excel (structured) | data/assays/ |
| PK/PD data | Excel / CSV | data/assays/in_vivo/ |
| Literature extractions | Excel / BibTeX | data/literature/ |

### 5.3 Code Standards
- Python 3.9+, type hints, docstrings
- Functions for repeated analyses
- Version control for analysis scripts
- Jupyter notebooks for exploration only (export to .py for production)

---

## 6. Session Artifacts & Persistence

### 6.1 Auto-Saved (Session Database)
- Todos and dependencies
- Checkpoint summaries
- Agent-to-agent handoff notes

### 6.2 Manual Save Required
**After each significant work session, save**:
- [ ] Updated compound library (new SMILES/SDF)
- [ ] Analysis scripts with results
- [ ] Protocol modifications
- [ ] Literature extraction tables
- [ ] Decision logs (why compound X was selected/dropped)
- [ ] Communication notes (if collaborating)

### 6.3 Checkpoint System
**Save after each phase completion**:
```
Phase 1 Complete → checkpoint_phase1.md
Phase 2 Complete → checkpoint_phase2.md  
Phase 3 Complete → checkpoint_phase3.md
Phase 4 Complete → checkpoint_phase4.md
```

**Checkpoint content**:
```markdown
# Checkpoint: Phase X Complete
**Date**: YYYY-MM-DD
**Completed milestones**: [List from research_framework.md]
**Key decisions**: [Why target Y chosen, why compound Z dropped]
**Data generated**: [Files created, locations]
**Next steps**: [Specific actions for Phase X+1]
**Open questions**: [Unresolved issues]
```

---

## 7. Inter-Session Continuity Protocol

### On New Session Startup:
1. **Read latest checkpoint** (highest phase number)
2. **Review pending todos** from session database
3. **Resume from last completed milestone**
4. **Update framework version** if protocols changed

### Version Bump Triggers:
- New target selection → framework vX.Y+1
- Protocol modification (>20% change) → framework vX.Y+1
- Major dataset addition (>50% increase) → literature vX.Y+1
- Software/tool upgrade → computational vX.Y+1

### Version Format: `v{major}.{minor}.{patch}`
- Major: Phase transition or fundamental approach change
- Minor: Protocol update, new target, significant data
- Patch: Bug fix, typo, minor clarification

---

## 8. Common Cold-Start Issues & Fixes

| Issue | Diagnosis | Fix |
|-------|-----------|-----|
| **Missing target structures** | PDB files not downloaded | Run download commands in Phase 2 protocol |
| **RDKit import error** | Environment not activated | `conda activate pnn_env` or install RDKit |
| **Vina not found** | Not in PATH | Add to PATH or use full path `/usr/local/bin/vina` |
| **No compound library** | data/compounds/ empty | Download from ChEMBL, ZINC, or generate |
| **Assay validation failed** | Positive controls not working | Check reagent expiry, protocol adherence |
| **PK below detection** | Compound metabolized too fast | Formulation optimization, prodrug strategy |

---

## 9. Decision Log Template

**Use for every major decision** (target selection, compound progression, assay choice):
```markdown
## Decision: [Title]
**Date**: YYYY-MM-DD
**Context**: [What prompted this decision]
**Options considered**: 
1. Option A - [pros/cons]
2. Option B - [pros/cons]
3. Option C - [pros/cons]
**Decision**: [Chosen option]
**Rationale**: [Scientific/strategic reasoning]
**Reversibility**: [Can this be changed later? How?]
**Dependencies**: [What depends on this decision]
```

---

## 10. Quick Reference: Key Parameters

### Target Product Profile (TPP)
| Parameter | Target Spec |
|-----------|-------------|
| **MW** | 300-400 Da |
| **logP** | 2.0-3.0 |
| **TPSA** | < 60 Å² |
| **HBD** | ≤ 2 |
| **HBA** | ≤ 5 |
| **pKa (basic)** | 7.5-8.5 |
| **Oral F** | > 20% |
| **Brain/plasma** | > 0.3 |
| **PNN EC50** | < 100 nM |
| **ITT improvement** | > 25% AUC reduction |
| **Food intake reduction** | > 15% |
| **hERG IC50** | > 30 μM |

### Primary Targets (Priority Order)
1. **PTPσ CSPG binding site** (Primary - allosteric, selective)
2. **PTPσ Catalytic site** (Backup - known inhibitors, selectivity challenge)
3. **LAR CSPG binding site** (Combination - dual PTPσ/LAR)
4. **ChABC mimic** (Long-term - enzymatic degradation)
5. **NgR1/3** (Alternative - peptide-based)

### Key Assays (Validation Order)
1. PNN degradation (fluorogenic DQ-CS)
2. BBB permeability (MDCK-MDR1 + PAMPA-BBB)
3. Insulin signaling (mHypoA cells, pIR/pAKT Western)
4. Cytotoxicity (CellTiter-Glo, CC50)
5. Hypothalamic slice (ex vivo WFA + pIR IHC)
6. In vivo PK (oral, brain/plasma)
7. In vivo PD (WFA, ITT, food intake)

---

## 11. Emergency Procedures

### If Completely Stuck:
1. Re-read README.md → literature_review.md → research_framework.md
2. Check latest checkpoint
3. Review open todos in session database
4. Pick ONE concrete next action (not "research more")
5. Execute and document

### If Data Loss Suspected:
1. Check `/files/` directory for session artifacts
2. Check session database for todos
3. Check git history (if repo initialized)
4. Reconstruct from checkpoints + literature

### If Results Contradict Hypothesis:
1. Document discrepancy precisely
2. Check positive controls worked
3. Verify compound identity (LC-MS)
4. Consider alternative mechanisms
5. Update hypothesis in literature_review.md
6. Adjust research_framework.md if needed

---

## 12. Contact & Escalation

**Project Owner**: Yunus-Emre
**Location**: c:\Users\Yunus-Emre\Desktop\PNN\
**Last Updated**: 2026-10-11
**Version**: 1.0

**For questions requiring human input**:
- Target selection changes
- Animal protocol modifications
- Budget/resource allocation
- Publication/patent decisions

---

*This guide ensures identical project initialization across all AI cold starts.*
*Compatible with Recursive Language Model workflows.*