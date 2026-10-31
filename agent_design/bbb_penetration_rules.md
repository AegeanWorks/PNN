# BBB Penetration Rules: Oral CNS Drug Design Guidelines

## Overview
This document defines the physicochemical rules, transporter strategies, and predictive models for achieving oral BBB penetration in small molecule drug design.

---

## 1. Physicochemical Rules for CNS Penetration

### 1.1 Modified Lipinski Rules (CNS Adaptation)
| Parameter | Standard Lipinski | CNS Adaptation | Rationale |
|-----------|------------------|----------------|-----------|
| **MW** | <500 Da | 300-450 Da | Balance between molecular recognition and passive diffusion |
| **logP (XlogP3)** | <5 | 1.5-3.5 | Optimal membrane partitioning without plasma protein binding |
| **TPSA** | <140 Å² | <70 Å² (optimal <60) | Reduced H-bonding for membrane crossing |
| **HBD** | <5 | ≤3 | Lower desolvation penalty at membrane interface |
| **HBA** | <10 | ≤7 | Balance permeability + solubility |
| **Rotatable bonds** | <10 | ≤7 | Limit entropy loss upon membrane partitioning |
| **pKa (basic)** | N/A | 7.5-8.5 | Brain trapping: protonated at pH 7.4, neutral for diffusion |
| **pKa (acidic)** | N/A | 3.5-5.0 | Neutral species at pH 7.4 for passive diffusion |
| **logD7.4** | N/A | 1.5-3.5 | Physiological pH partitioning coefficient |

### 1.2 Brain-to-Plasma Ratio (logBB) Predictions
| logBB Range | Classification | Typical Compound Properties |
|-------------|---------------|---------------------------|
| **logBB > 0.3** | Good penetration | MW < 400, logP 2-3, TPSA < 60 |
| **-0.3 < logBB < 0.3** | Moderate | MW 400-450, logP 1-4, TPSA 60-90 |
| **logBB < -0.3** | Poor penetration | MW > 450, logP < 1 or > 4, TPSA > 90 |
| **logBB < -1.0** | Very poor | Highly polar, charged, large molecules |

### 1.3 Yield Index (YI) for CNS Multiplicity
**Formula**: YI = 10^(logBB + 0.5 × (5 - TPSA/10) - 0.2 × |logP - 2.5|)

| YI Range | CNS Likelihood |
|----------|----------------|
| YI > 15  | High probability |
| 5 < YI < 15 | Moderate probability |
| YI < 5 | Low probability |

---

## 2. Transporter-Mediated Strategies

### 2.1 LAT1 (Large Neutral Amino Acid Transporter)
**Target Amino Acids**: L-leucine, L-phenylalanine, L-tyrosine, L-tryptophan

**Design Strategies**:
1. **Amino acid conjugates**: Replace polar groups with amino acid moieties
2. **D-amino acid mimetics**: Avoid peptidoglycan recognition
3. **Keto-acid analogs**: Pyruvate, α-ketoglutarate scaffolds
4. **Aromatic amino acids**: Higher affinity than aliphatic

**Affinity Requirements**: Km < 100 μM for effective competition
**Transporter Expression**: BBB endothelial cells, also in choroid plexus

**Design Examples**:
```
L-DOPA (levodopa):       F-COO-CH2-CH2-NH2   →   Brain uptake 5-10×
Amino acid esters:      R-COO-CH(NH2)-R'    →   Peptide transporter
Amino acid conjugates:  Drug-NH2 + amino acid → LAT1-mediated
```

### 2.2 MCT1 (Monocarboxylate Transporter-1)
**Target Substrates**: Pyruvate, lactate, acetoacetate, β-hydroxybutyrate

**Design Strategies**:
1. **Carboxylic acid prodrugs**: Attach acid to neutral scaffold
2. **Keto-analogs**: Pyruvate, acetoacetate scaffolds
3. **β-Hydroxybutyrate esters**: Ketone body mimics

**Affinity Requirements**: Km < 500 μM (lower affinity sufficient)
**Transporter Expression**: BBB, also in muscle, heart

**Design Examples**:
```
MCT1 substrates:        Pyruvate → Km ~0.1 mM
Prodrug strategy:       Drug-COOEt → hydrolyzed intracellularly
Keto-acid conjugates:   Drug-CO-COO⁻ → MCT1-mediated
```

### 2.3 GLUT1 (Glucose Transporter-1)
**Target Ligands**: Glucose, deoxyglucose, methylglucose

**Design Strategies**:
1. **Glucose conjugates**: Covalent attachment via linker
2. **2-DG analogs**: 2-deoxyglucose scaffolds
3. **SGLT1 substrates**: Sodium-glucose cotransporter

**Affinity Requirements**: Km < 10 mM (GLUT1 is high-capacity)
**Transporter Expression**: BBB (major), placenta, kidney

**Design Examples**:
```
Glucose conjugates:     Drug-O-Glc → GLUT1-mediated
Phenyl-β-D-glucoside:   Classic GLUT1 substrate
Deoxyglucose:           2-NBDG (fluorescent probe)
```

### 2.4 Transferrin Receptor (TfR) Mediated RMT
**Ligands**: Anti-TfR antibodies, transferrin, single-chain Fv (scFv), peptides

**Design Strategies**:
1. **Anti-TfR antibodies**: OX26 (rat), OX7 (humanized)
2. **Transferrin**: Iron-carrying glycoprotein
3. **Peptide mimetics**: 7-8 residue cyclic peptides
4. **scFv-Fc fusions**: Engineered for half-life

**Efficiency**: 2-5× brain uptake enhancement (anti-TfR)
**Challenge**: Large molecular weight, immunogenicity risk

**Design Examples**:
```
OX26 antibody:          Rat anti-rat TfR → 2-3× brain uptake
Transferrin conjugate:  Drug-Tf → receptor-mediated endocytosis
Cyclic peptide:         C(RADA)₃ → TfR binding, 3-5× enhancement
```

### 2.5 LRP1 (Low-Density Lipoprotein Receptor-Related Protein 1)
**Ligands**: Angiopep-2, ApoE, ApoJ, peptide mimetics

**Design Strategies**:
1. **Angiopep-2**: 18-amino acid peptide (TKPRWYSPYIWKLKDE)
2. **ApoE mimetics**: 10-12 residue amphipathic peptides
3. **Small molecule mimetics**: Being screened from libraries

**Efficiency**: 1.5-3× brain uptake enhancement
**Translocation mechanism**: Clathrin-mediated endocytosis

### 2.6 P-glycoprotein (P-gp) / BCRP Efflux Avoidance

**Efflux Ratio Thresholds**:
| Efflux Ratio (MDCK-MDR1) | Classification |
|-------------------------|----------------|
| < 2.0 | No significant efflux |
| 2.0 - 5.0 | Moderate efflux |
| 5.0 - 10.0 | High efflux |
| > 10.0 | Very high efflux |

**Structural Features to Avoid**:
- **Quaternary ammonium**: Permanent positive charge (P-gp substrate)
- **Multiple HBD**: >4 increases P-gp recognition
- **Planar aromatic systems**: Stack with P-gp binding pocket
- **Glucuronide conjugates**: P-gp efflux substrate
- **Fluoro-phenyl**: Common P-gp motif

**Structural Features to Emphasize**:
- **Aliphatic chains**: Reduce planar recognition
- **Saturated rings**: Less P-gp recognition than aromatic
- **Appropriate logP**: 2-3 optimal (too low or too high = P-gp)
- **Chiral centers**: Stereochemistry can avoid P-gp recognition

**P-gp Avoidance Rules**:
- MW < 400 + logP 2-3 + TPSA < 60 = typically low efflux
- Add methyl groups to block metabolic/P-gp sites
- Replace phenyl with cyclopropyl or other saturated rings

---

## 3. Predictive Modeling Tools

### 3.1 adMET (BAR AI Predictions)
| Predictor | Input | Output | Accuracy |
|-----------|-------|--------|----------|
| **Brain or AI** | SMILES, properties | logBB, brain/plasma ratio | ~85% (internal validation) |
| **adMET** | Structure | 12 ADMET properties including BBB | Broad coverage |
| **SPORCalc** | Physicochemical | logBB, BBB penetration | Historical benchmark |

### 3.2 QSAR Models
**Top-performing models** (from literature benchmarking):
1. **Yalkowsky logBB model**: logBB = 0.29 - 0.015 × TPSA - 0.047 × MW + 0.54 × logP
   - R² = 0.71, RMSE = 0.34
   
2. **SARpy model**: Machine learning (RF/XGBoost) on 2K+ compounds
   - Features: 2D descriptors, fingerprints, 3D conformers
   - Q² (cross-validated) = 0.68-0.79

3. ** adMET deep learning**: Neural network on PubChem + ChEMBL
   - Multi-task: BBB + solubility + metabolism
   - AUROC > 0.85 for binary classification

### 3.3 In Silico Screening Pipeline
```
1. Generate 3D conformers (OMEGA, ETKDG)
2. Calculate descriptors (RDKit, PaDEL, Dragon)
3. Predict logBB (adMET + Yalkowsky equation)
4. Predict efflux (P-gp, BCRP models)
5. Filter: logBB > -0.5 AND efflux ratio < 5
6. Rank by MPO score (multi-parameter optimization)
```

### 3.4 Property Prediction Commands (RDKit examples)
```python
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

mol = Chem.MolFromSmiles('CCOC(=O)c1ccccc1')  # example SMILES

# Basic properties
mw = Descriptors.MolWt(mol)
logp = Crippen.MolLogP(mol)
tpsa = Descriptors.TPSA(mol)
hbd = Lipinski.NumHDonors(mol)
hba = Lipinski.NumHAcceptors(mol)
rb = Descriptors.NumRotatableBonds(mol)

# logBB prediction (Yalkowsky)
logbb = 0.29 - 0.015*tpsa - 0.047*mw + 0.54*logp

print(f"MW: {mw:.1f}, logP: {logp:.2f}, TPSA: {tpsa:.1f}")
print(f"HBD: {hbd}, HBA: {hba}, RB: {rb}")
print(f"Predicted logBB: {logbb:.2f}")
```

---

## 4. Oral CNS Drug Design Database

### 4.1 FDA-Approved Oral CNS Drugs (n=35+)
| Drug | Class | logBB | MW | logP | TPSA | HBD | HBA |
|------|-------|-------|-----|------|------|-----|-----|
| **Loperamide** | Opioid | -0.9 | 477 | 5.08 | 46.5 | 1 | 5 |
| **Risperidone** | Antipsychotic | -0.3 | 410 | 5.89 | 57.3 | 2 | 7 |
| **Donepezil** | AChE inhibitor | -0.5 | 379 | 5.01 | 55.9 | 2 | 7 |
| **Venlafaxine** | SNRI | -0.4 | 278 | 3.01 | 56.7 | 2 | 6 |
| **Carbamazepine** | Anticonvulsant | -0.4 | 236 | 2.40 | 28.8 | 1 | 1 |
| **Valproic acid** | Anticonvulsant | -1.5 | 108 | -1.04 | 56.7 | 1 | 4 |
| **Gabapentin** | Anticonvulsant | -1.8 | 171 | -1.29 | 66.4 | 2 | 5 |
| **Pregabalin** | Anticonvulsant | -1.6 | 159 | 0.99 | 66.4 | 2 | 5 |
| **Lamotrigine** | Anticonvulsant | -0.8 | 256 | 2.90 | 58.8 | 0 | 5 |
| **Topiramate** | Anticonvulsant | -1.5 | 339 | 0.69 | 106.1 | 2 | 8 |

**Key Observations**:
- Most successful oral CNS drugs have MW < 450, logP 1-5, TPSA < 110
- Negative logBB compounds often have TPSA > 80 or ionizable groups
- Valproic acid (small, charged) has very poor logBB despite small MW

### 4.2 Successful Oral CNS Chemical Scaffolds
| Scaffold Type | Example Drugs | logBB Range | Notes |
|---------------|---------------|-------------|-------|
| **Benzamide** | Risperidone, amisulpride | -0.5 to -0.2 | Aromatic + basic side chain |
| **Phenylpiperidine** | Loperidine analogs | -0.8 to -0.2 | Moderate lipophilicity |
| **Benzocycloheptene** | Donepezil, tacrine | -0.6 to -0.3 | Fused ring system |
| **Indole** | Sumatriptan, vilazodone | -0.7 to -0.3 | Heteroaromatic |
| **Pyrimidine** | Various kinase inhibitors | -0.5 to 0.0 | Polar heterocycle |
| **Picolinamide** | Nicotine analogs | -0.3 to 0.2 | Pyridine + amide |

### 4.3 Failed CNS Orals (Lessons Learned)
| Drug | Issue | Lesson |
|------|-------|--------|
| **Selective serotonin reuptake inhibitors (early)** | Too polar (high TPSA) | TPSA < 90 recommended |
| **Many kinase inhibitors** | P-gp efflux + poor solubility | Balance logP, add solubilizing groups strategically |
| **Peptide hormones** | Proteolytic degradation + no passive diffusion | Use prodrugs or alternative delivery |
| **High-MW antibodies** | No BBB penetration | RMT strategies only, not passive diffusion |

---

## 5. Cold-Start Protocol for BBB Prediction

### 5.1 When Starting a New Session
1. **Load property prediction scripts** from `/docs/` or generate from scratch
2. **Set up compound library** in standardized format (SMILES)
3. **Run initial property filter**: MW 300-450, logP 1.5-4.0, TPSA < 90
4. **Apply YI score**: YI > 5 for further progression

### 5.2 Standardized SMILES Format
```
# All compounds stored as:
# - SMILES string (canonical)
# - Molecular weight
# - logP (Crippen)
# - TPSA
# - HBD, HBA, rotatable bonds
# - logBB (predicted)
# - YI score
```

### 5.3 Initial Filter Script (Python/RDKit)
```python
def cns_filter(smiles):
    """Apply CNS drug-likeness filters"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rb = Descriptors.NumRotatableBonds(mol)
    
    # CNS adaptive Lipinski
    violations = 0
    if not (300 <= mw <= 450):
        violations += 1
    if not (1.5 <= logp <= 4.0):
        violations += 1
    if tpsa > 90:
        violations += 1
    if hbd > 3:
        violations += 1
    if hba > 7:
        violations += 1
    if rb > 7:
        violations += 1
    
    return violations <= 2, f"MW={mw:.1f}, logP={logp:.2f}, TPSA={tpsa:.1f}, HBD={hbd}, HBA={hba}, RB={rb}"
```

### 5.4 YI Score Calculation
```python
def yi_score(logbb, tpsa, logp):
    """Yield Index for CNS multiplicity"""
    yi = 10**(logbb + 0.5*(5 - tpsa/10) - 0.2*abs(logp - 2.5))
    return yi

# Example usage
logbb = 0.29 - 0.015*tpsa - 0.047*mw + 0.54*logp
yi = yi_score(logbb, tpsa, logp)
```

---

## 6. Brain/Polar Surface Area Relationships

### 6.1 TPSA vs logBB Correlation
```
logBB = 0.29 - 0.015 × TPSA - 0.047 × MW + 0.54 × logP
         (R² = 0.71, n = 200+ compounds)
```

**Key thresholds**:
- **TPSA < 40 Å²**: Typically logBB > 0.5 (excellent penetration)
- **TPSA 40-60 Å²**: Typically logBB 0 to 0.3 (good penetration)
- **TPSA 60-80 Å²**: Typically logBB -0.3 to 0 (moderate)
- **TPSA 80-100 Å²**: Typically logBB -0.3 to -0.8 (poor)
- **TPSA > 100 Å²**: Typically logBB < -0.8 (very poor)

### 6.2 H-Bond Capacity vs Penetration
| HBD + HBA | Typical logBB | CNS Success Rate |
|-----------|---------------|------------------|
| ≤ 4 | > 0.3 | ~80% |
| 5-6 | 0 to 0.3 | ~50% |
| 7-8 | -0.3 to 0 | ~25% |
| ≥ 9 | < -0.3 | ~10% |

**Rule of thumb**: HBD + HBA ≤ 8 for reasonable CNS penetration

---

## 7. Metabolic Stability Considerations for CNS Drugs

### 7.1 CYP450 Interaction Risk
| Compound Feature | CYP3A4 Risk | CYP2D6 Risk | CYP2C9 Risk |
|------------------|-------------|-------------|-------------|
| **N-containing scaffolds** | Moderate | High (if basic) | Low |
| **Aryl heterocycles** | Moderate-High | Low-Moderate | Moderate |
| **Aliphatic esters** | High (esterase) | Low | Low |
| **Fluorinated aromatics** | Low-Moderate | Low | Low |
| **Sulfonic acids** | Very Low | Low | Low |

### 7.2 CNS-Specific Metabolic Challenges
1. **P-glycoprotein efflux** → increases apparent clearance
2. **Brain esterases** → prodrug activation in brain (desired) or periphery (undesired)
3. **Monoamine oxidase** → deamination of basic amines
4. **UGT/SULT conjugation** → glucuronide/sulfate formation, potential P-gp substrates

### 7.3 Stability Optimization Strategies
- **Fluorine substitution**: Blocks oxidative metabolism
- **Methyl blocking**: Prevents ω-oxidation of alkyl chains
- **Cyclization**: Reduces conformational flexibility + metabolic soft spots
- **Isosteres**: Replace labile groups (e.g., CF₃ for Cl, oxadiazole for amide)

---

## 8. Project-Specific BBB Penetration Criteria

### 7.1 Primary Criteria (Phase 1-2 Screening)
| Parameter | Minimum Requirement | Target |
|-----------|--------------------|--------|
| **MW** | 300-450 Da | 350 ± 50 Da |
| **logP** | 1.5-4.0 | 2.5 ± 0.5 |
| **TPSA** | < 90 Å² | < 70 Å² |
| **HBD** | ≤ 3 | ≤ 2 |
| **HBA** | ≤ 7 | ≤ 5 |
| **logBB (predicted)** | > -1.0 | > -0.5 |
| **YI score** | > 5 | > 15 |
| **Efflux ratio** | < 10 | < 5 |

### 7.2 Secondary Criteria (Phase 3+)
| Parameter | Minimum Requirement | Target |
|-----------|--------------------|--------|
| **Actual logBB (in vivo)** | > -1.0 | > -0.3 |
| **Brain/plasma ratio** | > 0.1 | > 0.3 |
| **Cmax/brain** | Detectable | > 10% of plasma Cmax |
| **Tmax (brain)** | < 4 h | < 2 h |

### 7.3 Failure Thresholds (Stop Criteria)
- **logBB < -1.5** → Reformulate or redesign
- **Efflux ratio > 15** → P-gp avoidance redesign
- **TPSA > 110** → Reduce polar surface
- **HBD > 4** → Mask or reduce H-bond donors
- **Oral F < 5%** → Formulation or prodrug strategy

---

*Document Version: 1.0*
*Generated: 2026-10-11*
*Compatible with Recursive Language Model cold-start workflows*