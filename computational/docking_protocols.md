# Computational Docking Protocols: PNN Target Screening

## Overview
This document defines the docking workflows, target preparation, and screening pipelines for virtual screening against PNN-related targets.

---

## 1. Target Preparation Workflows

### 1.1 PTPσ D1 Domain (Primary Target - CSPG Binding Site)

#### 1.1.1 Structure Retrieval
| Action | Details |
|--------|---------|
| **PDB ID** | 5K4P (PTPσ + CSPG tetrasaccharide) |
| **Resolution** | 2.8 Å |
| **Chain** | A (PTPσ), B (CS chain) |
| **Download** | `https://files.rcsb.org/download/5K4P.pdb` |

#### 1.1.2 Preparation Steps
```python
import pdbtools
from openmm.app import PDBFile, ForceField
from openmm import XmlSerializer

# 1. Download and extract PDB
pdb = PDBFile('5K4P.pdb')

# 2. Remove water molecules, keep only protein + ligand
atoms_to_remove = [i for i, atom in enumerate(pdb.topology.atoms) 
                   if atom.resname == 'HOH']
# (remove via pdbtools or custom script)

# 3. Add hydrogens (critical for docking)
from openmm.app import Element
pdb.topology.addHydrogens()

# 4. Save prepared structure
PDBFile.savePdb(pdb, '5K4P_prepared.pdb')
```

#### 1.1.3 Binding Site Definition
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Center coordinates** | X: ~12.5, Y: ~18.3, Z: ~22.1 | CSPG binding pocket |
| **Box size** | 20 × 20 × 20 Å | Encompasses entire binding site |
| **Key residues** | Lys201, Arg204, Lys207, Arg210 | CS chain recognition |
| **Aromatic residues** | Tyr199, Phe203 | π-stacking with CS |
| **Excluded volumes** | Catalytic Cys267, active site | Avoid covalent docking |

#### 1.1.4 Receptor Grid Generation (AutoDock Vina)
```bash
# Prepare receptor
vina --receptor 5K4P_prepared.pdb \
     --out receptor.pdbqt \
     --center_x 12.5 --center_y 18.3 --center_z 22.1 \
     --size_x 20 --size_y 20 --size_z 20

# Generate grid
vina --genlig 5K4P_prepared.pdbqt \
     --out grid.xml
```

### 1.2 PTPσ Catalytic Domain (Backup Target)

#### 1.2.1 Structure Retrieval
| Action | Details |
|--------|---------|
| **PDB ID** | 2PTS (PTPσ D1 domain) |
| **Resolution** | 2.0 Å |
| **Chain** | A (catalytic domain) |
| **Download** | `https://files.rcsb.org/download/2PTS.pdb` |

#### 1.2.2 Preparation Steps
Same as Section 1.1.2, with additional:
- **Metal ion**: Remove catalytic Zn²⁺ or keep for zinc-docking
- **Cysteine modifications**: Carbamidomethylation (CAM) if needed
- **Loop regions**: Residues 250-260 often flexible, consider truncation

### 1.3 LAR D1 Domain (Secondary Target)

#### 1.3.1 Structure Retrieval
| Action | Details |
|--------|---------|
| **PDB ID** | 2LAR (LAR D1 domain) |
| **Resolution** | 1.9 Å |
| **Chain** | A (catalytic domain) |
| **Download** | `https://files.rcsb.org/download/2LAR.pdb` |

#### 1.3.2 Preparation Notes
- **LAR vs PTPσ**: 42% sequence identity in catalytic domain
- **Key differences**: Wedge domain conformation, substrate binding loop
- **Dual targeting**: Possible with conserved catalytic site inhibitors

### 1.4 ChABC (Enzyme Mimic) Target Preparation

#### 1.4.1 Structural Considerations
- **Not a traditional protein target**: Screening for small molecule mimetics
- **Alternative approach**: Dock against CSPG core protein domains
- **Available structures**: Brevican G1 domain, aggrecan G3 domain

#### 1.4.2 CSPG Core Protein Targets
| Target | PDB ID | Notes |
|--------|--------|-------|
| **Brevican G1** | 4ZIC | Hyaluronan-binding domain |
| **Aggrecan G3** | 1R6J | Link protein binding domain |
| **Link protein** | 2N8L | CSPG aggregate stabilization |

---

## 2. Ligand Preparation & Library Design

### 2.1 Compound Library Composition

| Library Type | Size | Source | Diversity Metric |
|--------------|------|--------|------------------|
| **Fragment library** | 1,000-5,000 | Commercially available (e.g., Enamine, Maybridge) | MAE < 1.0 Å |
| **Focused PNN library** | 10,000-50,000 | Designed around CSPG-binding motifs | Sphericity 0.7-0.9 |
| **Repurposing library** | 5,000-20,000 | FDA-approved drugs (DrugBank) | Known safety |
| **Natural products** | 2,000-10,000 | NPASS, TCMSP | Bioactivity precedent |
| **In silico generated** | 50,000-200,000 | De novo design, scaffold hopping | Property-tuned |

### 2.2 Ligand Preprocessing Workflow

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

def prepare_ligand(smiles):
    """Standardized ligand preparation for docking"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # 1. Remove salts/ionize
    cleaner = rdMolStandardize.ChargeParent()
    mol = cleaner.standardize(mol)
    
    # 2. Neutralize charges (for most docking)
    # (or keep charged for explicit receptor grids)
    
    # 3. Generate 3D conformers
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFF94OptimizeMolecule(mol)
    
    # 4. Generate multiple conformers
    AllChem.ConformerRange(mol, 0, 50)  # 50 conformers max
    
    # 5. Save as PDBQT
    # (use OpenBabel: obabel -i mol2 input.mol2 -o pdbqt output.pdbqt)
    
    return mol

# Batch processing example
compounds = ["CCOC(=O)c1ccccc1", "c1ccccc1C(=O)O"]  # example SMILES
prepared = [prepare_ligand(smiles) for smiles in compounds]
```

### 2.3 Property Filters (Pre-Docking)

| Filter | Threshold | Purpose |
|--------|-----------|---------|
| **MW** | 150-500 Da | Drug-likeness |
| **logP** | -2 to 8 | Permeability range |
| **TPSA** | < 140 Å² | General drug-likeness |
| **HBD** | ≤ 5 | Lipinski |
| **HBA** | ≤ 10 | Lipinski |
| **Rotatable bonds** | ≤ 10 | Flexibility |
| **PAINS** | Negative | Avoid false positives |
| **Brenk** | Negative | Substructure alerts |

### 2.4 Scaffold Hopping Strategies

**Goal**: Identify privileged scaffolds for PNN binding

**Approaches**:
1. **Template matching**: Known CSPG binders → scaffold replacement
2. **Fragment growing**: Start from fragments, add groups
3. **Linking**: Connect two fragments with optimal spacer
4. **Scaffold replacement**: Replace core while maintaining interactions

**Example**: If benzamide scaffold binds PTPσ CSPG site → replace with:
- Phthalimide (similar H-bond pattern)
- Isoquinolinone (extended π-system)
- Quinazolinone (additional H-bond acceptor)

---

## 3. Docking Protocols

### 3.1 AutoDock Vina Protocol

#### 3.1.1 Standard Protocol
```bash
# 1. Prepare receptor (already done in Section 1)
# 2. Prepare ligand (already done in Section 2)

# 3. Run docking
vina --receptor receptor.pdbqt \
     --ligand ligand.pdbqt \
     --out output.pdbqt \
     --log output.log \
     --center_x 12.5 --center_y 18.3 --center_z 22.1 \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 8 \
     --num_modes 20 \
     --energy_range 3

# 4. Parse results
# Top scoring modes in output.pdbqt
# Log contains: affinity, rmsd, and scoring terms
```

#### 3.1.2 Exhaustiveness Settings
| Use Case | Exhaustiveness | Time per Ligand | Recommendation |
|----------|---------------|-----------------|----------------|
| **High-throughput** | 4-6 | 1-3 min | Screen 1K+ compounds/day |
| **Standard screening** | 8-10 | 3-8 min | Screen 200-500 compounds |
| **Lead optimization** | 12-16 | 8-15 min | Top 50 compounds, detailed analysis |
| **Finalists** | 20-30 | 15-30 min | Top 5-10 compounds, full analysis |

#### 3.1.3 Scoring Terms & Interpretation
```
AutoDock Vina scoring function:
  E_total = E_binding + E_receptor + E_ligand

Reported in output.log:
  - affinity: Binding affinity (kcal/mol), more negative = better
  - rmsd_ligand: RMSD from original ligand pose
  - rmsd_receptor: RMSD from receptor starting structure
  - e_inter: Non-bonded interaction energy
  - e_rep: Pauli repulsion energy
  - e_att: Attractive van der Waals energy
  - e_hb: Hydrogen bonding energy
  - e_ele: Electrostatic energy
```

**Interpretation guidelines**:
- **affinity < -7 kcal/mol**: Potential binders (context-dependent)
- **affinity < -9 kcal/mol**: Strong candidates
- **affinity < -11 kcal/mol**: Very strong binders (rare, check for artifacts)
- **Hydrogen bonds**: ≥ 2-3 optimal in binding pose
- **Key interactions**: Salt bridges with Lys/Arg, π-stacking with aromatic residues

### 3.2 Glide SP/XP Protocol (Schrödinger, optional)

#### 3.2.1 Standard Precision (SP)
```bash
# Prepare receptor for Glide
# (Use Protein Preparation Wizard in Maestro)

# Run Glide SP
glide -receptor 5K4P_prepared.pdb -ligand ligand.mol2 \
      -gridfile grid.sp.grd -out output_sp.sdf \
      -precision SP

# Score results
# GlideScore: typically -0.5 to -11.0
# Higher (less negative) = better for SP
```

#### 3.2.2 Extra Precision (XP)
```bash
# Run Glide XP (same grid, different precision)
glide -receptor 5K4P_prepared.pdb -ligand ligand.mol2 \
      -gridfile grid.xp.grd -out output_xp.sdf \
      -precision XP

# XP scoring better than SP but more computationally expensive
# Recommended: SP first, XP top 20-50 hits
```

### 3.3 rDock (RNA/Protein Docking, alternative)

#### 3.3.1 Advantages for PNN Targets
- **Faster**: ~2-3× faster than AutoDock Vina
- **Better for charged ligands**: Explicit treatment of electrostatics
- **Open source**: No license required
- **Good for**: PTPσ catalytic site, LAR catalytic site

#### 3.3.2 rDock Protocol
```bash
# Prepare receptor (similar to AutoDock)
# rDock has its own preparation scripts

# Run rDock
rdock -r receptor.pdbqt -l ligand.pdbqt \
      -o output.sdf -s 10 -n 20

# -s: scoring function (1=GoldScore, 2=rDockScore)
# -n: number of poses to generate
```

### 3.4 Consensus Scoring

**Why consensus?** Single scoring functions have limitations:
- **AutoDock Vina**: Good for general screening, less accurate for charged species
- **Glide XP**: Better accuracy but proprietary
- **rDock**: Good alternative, different algorithm
- **MM-GBSA**: Post-docking refinement, more accurate but slower

**Consensus protocol**:
1. Dock with AutoDock Vina (all compounds)
2. Top 50-100: Dock with Glide SP
3. Top 20-30: Dock with rDock
4. Top 10: MM-GBSA refinement (Prime MM-GBSA or amber)
5. Rank by: 0.5 × Vina + 0.3 × Glide + 0.2 × MM-GBSA score

**Executable ranking**: export the per-target scores to the contract in
[docking_ranking.md](docking_ranking.md) and run `python -m pnn.rank`. The tool
min-max normalizes each engine within the run before applying these weights, so
raw values with different scales are not summed directly.

---

## 4. MD Simulation Protocols

### 4.1 System Setup for Protein-Ligand Complexes

#### 4.1.1 Force Field Selection
| Component | Force Field | Version | Notes |
|-----------|-------------|---------|-------|
| **Protein** | AMBER99SB-ILDN | 2.1 | Improved backbone dynamics |
| **Ligand** | GAFF2 | 2.1 | General AMBER force field |
| **Charges** | AM1-BCC | N/A | Derived from ESP fit |
| **Ions** | TIP3P water | N/A | Standard water model |

#### 4.1.2 System Building (AMBER mdtin)
```bash
# 1. Load prepared complex
parmed.loadpdb complex.pdb

# 2. Add missing residues, fix termini
parmed.edit "addmissing @1-@N"

# 3. Assign force fields
parmed.antechamber -i complex.in -fi pdb -o complex.prepi -fo prepi
parmed.tleap <<EOF
source leaprc.protein.ff19SB
source leaprc.gaff2
complex = loadpdb complex.prepi
# Solvate
complex = addSolventBox complex, type='TIP3P', size=[10, 10, 10] nm
# Add ions
complex = addIons complex, 'Na+', 0.15  # physiological concentration
complex = addIons complex, 'Cl-', 0.15
saveambercomplex complex, 'complex'
EOF

# 4. Minimization, equilibration, production
# (see full protocol below)
```

### 4.2 Minimization Protocol

| Stage | Atoms Restrained | Max Force | Steps |
|-------|-----------------|-----------|-------|
| **1. Heavy atom only** | All heavy atoms | 1000 kJ/mol/nm | 1000 |
| **2. Protein only** | Protein heavy atoms | 500 kJ/mol/nm | 500 |
| **3. No restraints** | None | 100 kJ/mol/nm | 1000 |
| **4. Optional** | Weak (10 kJ/mol/nm) | 10 kJ/mol/nm | 500 |

**Minimization script (AMBER)**:
```bash
# Minimization cycle 1: water only (no solute)
pmemd.cuda -O -i min1.in -p complex.prm -c complex.rst7 \
          -r min1.rst7 -inf min1.info

# Minimization cycle 2: solute with restraints
pmemd.cuda -O -i min2.in -p complex.prm -c min1.rst7 \
          -r min2.rst7 -inf min2.info -ref min1.rst7
```

### 4.3 Equilibration Protocol

| Phase | Duration | Restraints | Temperature | Pressure |
|-------|----------|------------|-------------|----------|
| **1. NVT** | 100 ps | Protein heavy atoms (5 kJ/mol/nm) | 0 → 310 K (Langevin) | - |
| **2. NPT** | 500 ps | Protein heavy atoms (1 kJ/mol/nm) | 310 K (Langevin) | 1 bar (Berendsen) |
| **3. NPT** | 1 ns | Protein Cα (1 kJ/mol/nm) | 310 K (Langevin) | 1 bar (Parrinello-Rahman) |
| **4. Production** | 10-100 ns | None | 310 K (NPT) | 1 bar |

**Equilibration script (AMBER)**:
```bash
# NVT equilibration
pmemd.cuda -O -i eq1.in -p complex.prm -c min2.rst7 \
           -r eq1.rst7 -x eq1.mdcrd -inf eq1.info \
           -cpi eq1.cpi -t 310.0

# NPT equilibration (gradual restraint removal)
pmemd.cuda -O -i eq2.in -p complex.prm -c eq1.rst7 \
           -r eq2.rst7 -x eq2.mdcrd -inf eq2.info \
           -cpi eq2.cpi -t 310.0
```

### 4.4 Production MD

| Parameter | Recommended Value |
|-----------|-------------------|
| **Duration** | 100 ns minimum per top pose |
| **Time step** | 2 fs (with SHAKE on bonds to H) |
| **Trajectory saving** | Every 10 ps |
| **Temperature** | 310 K (Langevin thermostat) |
| **Pressure** | 1 bar (Parrinello-Rahman barostat) |
| **Cutoffs** | 12 Å (short-range), PME for electrostatics |
| **Constraints** | SHAKE on all bonds to H |

**Production script (AMBER)**:
```bash
pmemd.cuda -O -i prod.in -p complex.prm -c eq2.rst7 \
          -r prod.rst7 -x prod.mdcrd -inf prod.info \
          -cpi prod.cpi -t 310.0 -p 1.0
```

### 4.5 MD Analysis Metrics

#### 4.5.1 Binding Stability
| Metric | Calculation | Success Criterion |
|--------|-------------|-------------------|
| **RMSD (ligand)** | Average over last 80 ns | < 2.0 Å (stable pose) |
| **RMSD (protein)** | Average over binding site | < 2.5 Å (fold maintained) |
| **RMSF (per residue)** | Fluctuation per residue | < 2.0 Å in binding site |
| **Distance metrics** | Key interaction distances | Stable over time |

#### 4.5.2 Interaction Analysis
| Interaction Type | Analysis Method | Longevity Criterion |
|-----------------|-----------------|---------------------|
| **Hydrogen bonds** | `hydrogen-bond` (MDAnalysis) | > 40% occupancy |
| **Salt bridges** | Distance < 4 Å between charged pairs | > 30% occupancy |
| **π-π stacking** | Angle + distance criteria | > 25% occupancy |
| **Hydrophobic contacts** | Carbon-carbon distances | > 30% occupancy |

#### 4.5.3 MM-GBSA Binding Free Energy
```bash
# Prime MM-GBSA (Schrödinger) or MMPBSA.py (AMBER)

# MMPBSA.py example
MMPBSA.py -O -i mmpbsa.in -sp complex.prm -cp complex.prm \
          -rp receptor.prm -lp ligand.prm -y prod.mdcrd

# Output: ΔG_binding = ΔE_MM + ΔG_solv - TΔS
# ΔE_MM: bonded + non-bonded energy difference
# ΔG_solv: GB + SA solvation free energy
# -TΔS: entropy (normal mode or interaction energy)
```

**Typical ΔG_binding values**:
- **Strong binders**: -40 to -60 kcal/mol
- **Moderate binders**: -20 to -40 kcal/mol
- **Weak binders**: -10 to -20 kcal/mol
- **Non-binders**: > -10 kcal/mol (or positive)

---

## 5. Machine Learning Models

### 5.1 Property Prediction Models

| Model | Input Features | Target Property | Accuracy (RMSE) |
|-------|---------------|-----------------|-----------------|
| **Yalkowsky logBB** | MW, TPSA, logP, HBD, HBA | logBB (brain/plasma) | ~0.34 |
| **adMET deep learning** | 2D/3D descriptors, SMILES | 12 ADMET properties | Varies by property |
| **CYP3A4 inhibition** | Structural fingerprints | Binary (inhibitor/non) | AUROC ~0.85 |
| **hERG liability** | 2D descriptors | Binary (liability/non) | AUROC ~0.82 |
| **Solubility (logS)** | Physicochemical descriptors | logS (mol/L) | ~0.5-0.7 |

### 5.2 QSAR for PNN Binding

**Goal**: Predict PNN binding probability from compound structure

**Feature Types**:
1. **2D fingerprints**: MACCS keys (166 bits), PubChem fingerprints (881 bits)
2. **Graph neural networks**: MPNN, GraphConv (SMILES → embedding)
3. **3D descriptors**: MOE, Dragon, custom interaction fingerprints
4. **Physicochemical**: logP, TPSA, MW, HBD/HBA, pKa

**Model Architectures**:
- **Random Forest**: Interpretable, good baseline
- **XGBoost**: High performance, handles imbalanced data
- **Graph Neural Network**: State-of-the-art for molecular property prediction
- **Transformer (MoleculeBERT, ChemBERTa)**: Large-scale pretrained models

**Training Protocol**:
1. **Dataset**: 500-2,000 actives + 5,000-20,000 inactives (PTPσ/CSPG binding)
2. **Positive data**: From docking scores, experimental binding, or CSPG competition assays
3. **Negative data**: Random compounds, decoys (DUD-E format)
4. **Split**: 80% train, 10% validation, 10% test (stratified by scaffold)
5. **Metrics**: AUROC, AUPRC, BEDROC (early enrichment focus)

### 5.3 De Novo Design

**Goal**: Generate novel compounds with predicted PNN binding + BBB penetration

**Tools**:
1. **MoleculeNet benchmarks**: Standardized evaluation
2. **Junction Tree VAE**: Graph generation with molecular constraints
3. **GraphGAN**: Adversarial graph generation
4. **Tensorial Transformer**: Large-scale molecular generation
5. **REINVENT**: Reinforcement learning for molecule optimization

**Design workflow**:
```
Step 1: Train on active compounds + property data
Step 2: Generate 1,000+ novel molecules
Step 3: Filter by property constraints:
  - MW: 300-450
  - logP: 1.5-4.0
  - TPSA: < 90
  - logBB: > -0.5 (predicted)
  - Synthetic accessibility: SA score < 3
Step 4: Dock top 100 against PTPσ/LAR targets
Step 5: Rescore with ML model
Step 6: Select 20-50 for synthesis
```

### 5.4 Virtual Screening Pipeline

```
1. Library preparation (50K-200K compounds)
2. Property filtering (Lipinski + CNS rules)
3. Receptor grid generation (AutoDock Vina)
4. High-throughput docking (exhaustiveness 8)
5. Score ranking (Vina affinity)
6. Top 1-5%: Re-dock with higher precision
7. Top 0.1-0.5%: MD simulation (50 ns)
8. Top 10-20: MM-GBSA refinement
9. Final ranking (consensus scoring)
10. Select 20-50 for in vitro testing
```

**Expected hit rates**:
- **Primary screen**: 0.1-1% actives (50-500 from 50K library)
- **Secondary screen**: 10-30% confirmation rate
- **Tertiary (MD)**: 20-50% stable binders
- **Experimental**: 30-60% active confirmation

---

## 6. Cold-Start Protocol for Computational Work

### 6.1 When Starting a New Session
1. **Load computational scripts** from `/computational/`
2. **Set up environment**: Python 3.9+, RDKit, OpenMM, AutoDock Vina
3. **Download target structures** (PDB IDs listed in Section 1)
4. **Initialize compound library** from `/data/compounds/`

### 6.2 Required Software & Versions
| Software | Version | Purpose |
|----------|---------|---------|
| **Python** | 3.9-3.11 | General scripting |
| **RDKit** | 2023.09.1 | Molecule preparation |
| **AutoDock Vina** | 1.2.0 | Docking (primary) |
| **OpenMM** | 8.0+ | MD simulations |
| **pymol** | 2.5+ | Visualization |
| **NumPy/Pandas** | 1.5+/2.0+ | Data analysis |
| **SciPy** | 1.10+ | Statistical analysis |

### 6.3 Initialization Script (Python)
```python
#!/usr/bin/env python3
"""Cold-start initialization for PNN computational project"""

import sys
import os
import platform

print("=" * 60)
print("PNN PROJECT - Computational Module Cold-Start")
print("=" * 60)

# 1. Check Python version
print(f"\nPython version: {sys.version}")
assert sys.version_info >= (3, 9), "Python >= 3.9 required"

# 2. Check RDKit availability
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    print("RDKit: ✓ Available")
except ImportError:
    print("RDKit: ✗ Not installed, installing...")
    os.system("pip install rdkit")
    from rdkit import Chem

# 3. Check AutoDock Vina
vina_path = "/usr/local/bin/vina"
if not os.path.exists(vina_path):
    print(f"Warning: AutoDock Vina not found at {vina_path}")
    print("Please install or update path")

# 4. Set project directories
project_dir = os.path.expanduser("~/Desktop/PNN")  # or cwd
os.chdir(project_dir)
print(f"\nWorking directory: {os.getcwd()}")

# 5. Load target structures
targets = {
    "PTPσ_CSPG": "5K4P_prepared.pdb",
    "PTPσ_catalytic": "2PTS_prepared.pdb",
    "LAR_D1": "2LAR_prepared.pdb",
}

print("\nTarget structures:")
for name, pdb in targets.items():
    if os.path.exists(pdb):
        print(f"  {name}: ✓ {pdb}")
    else:
        print(f"  {name}: ✗ Missing - download from RCSB")

# 6. List available compounds
compound_dir = os.path.join(project_dir, "data", "compounds")
if os.path.exists(compound_dir):
    sdf_files = [f for f in os.listdir(compound_dir) if f.endswith('.sdf')]
    print(f"\nCompound libraries: {len(sdf_files)} SDF files")
else:
    print(f"\nCompound directory: {compound_dir} (not found)")

print("\n" + "=" * 60)
print("Initialization complete. Ready for screening.")
print("=" * 60)
```

### 6.4 First-Run Checklist
- [ ] Python environment activated
- [ ] RDKit installed and working
- [ ] AutoDock Vina installed and executable
- [ ] Target PDB files downloaded (5K4P, 2PTS, 2LAR)
- [ ] Compound library in SDF/SMILES format
- [ ] Project directories created (README, data/, computational/)
- [ ] Session todos initialized from research_framework.md

### 6.5 Common Cold-Start Issues & Fixes

| Issue | Cause | Fix |
|-------|-------|-----|
| **RDKit import error** | Version mismatch, missing system deps | `conda install -c conda-forge rdkit`, install system libs |
| **Vina not found** | Not in PATH, wrong architecture | Add to PATH, install vina binary for your OS |
| **PDB download failed** | Network issue, PDB offline | Manual download from https://www.rcsb.org |
| **Memory error (MD)** | Too large system, insufficient RAM | Reduce system size, decrease simulation time |
| **Conformer generation fail** | Unusual scaffold, bad SMILES | Check SMILES validity, use alternative representation |

---

## 7. Project-Specific Screening Parameters

### 7.1 Primary Virtual Screening Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Library size** | 50,000-200,000 | Manageable, diverse |
| **Docking program** | AutoDock Vina (primary) | Open source, proven |
| **Exhaustiveness** | 8 | Standard screening balance |
| **Box center** | Target-specific | From Section 1 |
| **Box size** | 20 × 20 × 20 Å | Encompasses binding site |
| **Max modes per ligand** | 20 | Enough for pose selection |
| **Energy range** | 3 kcal/mol | Discard high-energy poses |

### 7.2 Property Filter Thresholds (Pre-Docking)

| Property | Lower | Upper | Rationale |
|----------|-------|-------|-----------|
| **MW** | 150 | 500 | Drug-likeness |
| **logP** | -2.0 | 8.0 | Permeability-solubility balance |
| **TPSA** | 0 | 140 | General drug-likeness |
| **HBD** | 0 | 5 | Lipinski |
| **HBA** | 0 | 10 | Lipinski |
| **Rotatable bonds** | 0 | 10 | Flexibility limit |
| **MW (CNS)** | 300 | 450 | BBB penetration |
| **logP (CNS)** | 1.5 | 4.0 | Oral + CNS optimization |
| **TPSA (CNS)** | 0 | 90 | CNS penetration rule |
| **logBB (predicted)** | -2.0 | 0.5 | Brain penetration |
| **YI score** | 1 | 100 | CNS multiplicity |

### 7.3 Post-Docking Filtering Cascade

```
Step 1: Property filter (all compounds) → ~80% pass
Step 2: Vina docking (all) → ~20% have poses
Step 3: Vina score cutoff (affinity < -7) → ~10% remain
Step 4: Interaction check (H-bonds + key residues) → ~5% remain
Step 5: MD stability (top 50, 50 ns) → ~20-30 stable
Step 6: MM-GBSA refinement (top 20) → ~10-15 ranked
Step 7: Final selection (top 10-20) → synthesize + test
```

### 7.4 Decoy Selection (For Enrichment Testing)

**Use DUD-E (Directory of Useful Decoys)**:
- **Format**: SDF files with known active/inactive pairs
- **Ratio**: 1 active : 9-20 decoys
- **Purpose**: Validate virtual screening pipeline enrichment

**Example**: If screening PTPσ CSPG site:
- 1 known actives from literature
- 50-100 DUD-E decoys matching physicochemical space
- Measure: ROC-AUC, BEDROC, EF1% (enrichment factor at 1%)

**Expected enrichment**: EF1% > 10 indicates good virtual screening performance.

---

## 8. Output Organization

### 8.1 Directory Structure
```
computational/
├── docking/
│   ├── vina_results/
│   │   ├── run_1/
│   │   │   ├── ligand_001.pdbqt
│   │   │   ├── output.pdbqt
│   │   │   └── output.log
│   │   ├── run_2/
│   │   │   └── ...
│   │   └── summary_vina.xlsx
│   ├── glide_results/ (optional)
│   └── r_dock_results/ (optional)
├── md_simulations/
│   ├── top_50/
│   │   ├── system_setup/
│   │   ├── min1/
│   │   ├── eq1/
│   │   ├── eq2/
│   │   └── prod/
│   │       └── traj.mdcrd
│   └── analysis/
│       ├── rmsd_protein.xvg
│       ├── rmsd_ligand.xvg
│       ├── interaction_occupancy.txt
│       └── mm_gbsa_results.csv
├── ml_models/
│   ├── qsar_rf.pkl (Random Forest model)
│   ├── qsar_xgb.pkl (XGBoost model)
│   ├── molecular_descriptors.csv
│   └── predictions_top50.csv
└── reports/
    ├── screening_summary.pdf
    ├── top_compounds_ranked.xlsx
    └── pipeline_metrics.txt
```

### 8.2 Key Output Files (Per Compound)

**Vina output.log** contains:
```
======================================
  Center x,y,z: 12.5 18.3 22.1
  Size x,y,z: 20 20 20
  Exhaustiveness: 8
  Number of modes: 20
  Energy range: 3

Affinity: -8.5 kcal/mol (best mode)
RMSD_ligam: 2.1 Å (from input)
RMSD_receptor: 0.5 Å

e_inter: -45.2 kcal/mol
e_rep: 12.1 kcal/mol
e_att: -32.4 kcal/mol
e_hb: -5.8 kcal/mol (3 H-bonds)
e_ele: -18.7 kcal/mol
```

**Summary Excel** columns:
| Column | Description |
|--------|-------------|
| compound_id | Unique identifier |
| smiles | SMILES string |
| vina_affinity | Docking score (kcal/mol) |
| vina_rmsd | Ligand RMSD (Å) |
| hbond_count | Number of H-bonds |
| key_interactions | Residue interaction summary |
| predicted_logBB | BBB penetration prediction |
| yi_score | CNS multiplicity index |
| ml_binding_prob | ML model binding probability |
| final_rank | Overall ranking position |

---

## 9. Quality Control & Validation

### 9.1 Docking Validation (Internal)

**Redocking test**: Remove ligand from complex, redock, compare to original

```bash
# Remove ligand, keep protein
# Run Vina redock
# Compare RMSD of redocked pose to original

Success criterion: RMSD < 2.0 Å for > 75% of cases
```

**Benchmark set**: 10-20 diverse protein-ligand complexes from PDBbind

### 9.2 ML Model Validation

**Training/test split**: 80/20 stratified

**Metrics to report**:
- **AUROC**: Area under ROC curve (overall discrimination)
- **AUPRC**: Area under PR curve (better for imbalanced data)
- **BEDROC**: Bias early enrichment (α = 20 for screening)
- **EF1%**: Enrichment factor at 1% of library
- **RMSE**: For regression properties (logBB, logS)

**Validation set requirements**:
- Minimum 50 test compounds
- Scaffold diversity (no near-duplicates in train/test)
- External test set if available (different source)

### 9.3 MD Simulation Quality Checks

| Check | Method | Acceptance Criterion |
|-------|--------|---------------------|
| **Energy conservation** | Plot TE vs time | Drift < 10% over production |
| **RMSD stability** | Plot RMSD vs time | Plateau after equilibration |
| **Temperature stability** | Time series analysis | Within ±5 K of 310 K |
| **Pressure stability** | Time series analysis | Within ±10% of 1 bar |
| **Water density** | Radial distribution function | Near TIP3P reference |
| **Ligand RMSD** | Over production interval | < 2.0 Å (stable pose) |

### 9.4 Reporting Standards

**Per screening campaign, report**:
1. **Library statistics**: Size, diversity metrics, property distribution
2. **Docking statistics**: % compounds docked, score distribution, validation RMSD
3. **ML model performance**: AUROC, AUPRC on validation set
4. **Top compound profiles**: All computed properties, predicted activities
5. **Discrepancies**: Unusual results, outliers, manual inspection notes

**Report format**: PDF summary + Excel data + JSON metadata

---

*Document Version: 1.0*
*Generated: 2026-10-11*
*Compatible with Recursive Language Model cold-start workflows*