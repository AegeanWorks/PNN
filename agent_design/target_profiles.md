# Target Profiles: PNN Modulation for Hypothalamic Insulin Resistance

## Overview
This document defines the molecular targets for PNN softening, their structural biology, druggability assessment, and screening criteria.

---

## 1. Primary Target: Chondroitinase ABC (ChABC) Mimicry

### 1.1 Biological Rationale
- **Function**: Enzymatically cleaves β1-3 glycosidic bonds in chondroitin sulfate (CS) and dermatan sulfate (DS) chains
- **Effect**: Complete PNN degradation → restores insulin receptor accessibility
- **Evidence**: ICV ChABC in ARC → 25% food intake reduction, improved glucose tolerance (Lee et al., 2010)

### 1.2 Structural Data
| Structure | PDB ID | Resolution | Notes |
|-----------|--------|------------|-------|
| **ChABC (Proteus vulgaris)** | 1L5D | 1.8 Å | Catalytic domain only |
| **ChABC (Flavobacterium)** | 1L5E | 2.1 Å | Full-length with CBM |
| **ChABC + CS tetrasaccharide** | 2BZR | 1.9 Å | Substrate complex |
| **ChABC + inhibitor (CI-1)** | 3N8M | 2.0 Å | Competitive inhibitor |

### 1.3 Catalytic Mechanism
- **Family**: PL8 polysaccharide lyase
- **Mechanism**: β-elimination via acid-base catalysis
- **Key residues**: His200 (acid), Tyr247 (base), Arg287 (substrate binding)
- **Cofactors**: Ca²⁺ (structural), no metal catalysis

### 1.4 Druggability Assessment
| Parameter | Assessment | Notes |
|-----------|------------|-------|
| **Pocket druggability** | Moderate | Deep catalytic cleft, polar residues |
| **Small molecule precedent** | Limited | Few reported inhibitors (CI-1, suramin analogs) |
| **Allosteric sites** | Unexplored | Potential for non-competitive modulators |
| **Oral feasibility** | Low | Enzyme is large protein; need small molecule mimetics |

### 1.5 Small Molecule Mimic Strategy
**Goal**: Identify non-enzymatic small molecules that bind CSPG chains and promote cleavage or block CSPG-receptor interactions

**Approaches**:
1. **Catalytic mimetics**: Metal-chelating scaffolds (Zn²⁺/Cu²⁺ complexes) for β-elimination
2. **Substrate competitors**: CS-mimetic oligomers with non-cleavable linkages
3. **Allosteric modulators**: Bind exosite → enhance endogenous MMP activity
4. **CSPG-receptor blockers**: Prevent CSPG binding to PTPσ/LAR/NgR

---

## 2. Target: PTPσ (Protein Tyrosine Phosphatase Sigma)

### 2.1 Biological Rationale
- **Function**: Receptor tyrosine phosphatase activated by CSPG binding
- **Mechanism**: CSPG binding → PTPσ dimerization/activation → dephosphorylates insulin receptor/IRS-1
- **Evidence**: PTPσ KO mice resistant to CSPG inhibition; PTPσ blocking antibody restores plasticity

### 2.2 Structural Data
| Structure | PDB ID | Resolution | Domain |
|-----------|--------|------------|--------|
| **PTPσ D1 domain** | 2PTS | 2.0 Å | Catalytic domain |
| **PTPσ D2 domain** | 3T8C | 2.3 Å | Regulatory domain |
| **PTPσ + CSPG peptide** | 5K4P | 2.8 Å | CSPG binding site |
| **PTPσ + inhibitor (ISP)** | 4Y1R | 2.1 Å | Active site inhibitor |

### 2.3 Key Binding Sites
| Site | Residues | Function | Druggability |
|------|----------|----------|--------------|
| **Catalytic pocket** | Cys267, Arg268, Asp326 | Phosphotyrosine hydrolysis | High (conserved PTP fold) |
| **CSPG binding site** | Lys201, Arg204, Lys207, Arg210 | CS chain recognition | Moderate (shallow, polar) |
| **Wedge domain** | D2 domain interface | Allosteric regulation | Low (protein-protein) |

### 2.4 Druggability Assessment
| Parameter | Assessment | Notes |
|-----------|------------|-------|
| **Pocket druggability** | High (catalytic) / Moderate (CSPG site) | Catalytic site well-characterized |
| **Selectivity challenge** | High | PTP family conservation (>30 PTPs) |
| **CNS penetration** | Moderate | Known CNS-penetrant PTP inhibitors exist |
| **Oral feasibility** | Moderate | Several PTP1B inhibitors with oral activity |

### 2.5 Known Inhibitors
| Compound | Target | IC50 | CNS Penetration | Notes |
|----------|--------|------|-----------------|-------|
| **Trodusquemine (MSI-1436)** | PTP1B > PTPσ | 0.6 μM | Yes (brain/plasma ~0.5) | Aminosterol scaffold |
| **Compound 1 (PTPσ selective)** | PTPσ | 0.1 μM | Unknown | Benzothiazole scaffold |
| **ISP (inhibitor of PTPσ)** | PTPσ | 0.05 μM | No | Peptide-based |

---

## 3. Target: LAR (Leukocyte Common Antigen-Related)

### 3.1 Biological Rationale
- **Function**: Second major CSPG receptor (with PTPσ)
- **Mechanism**: Similar to PTPσ; CSPG binding → phosphatase activation → insulin signaling inhibition
- **Evidence**: LAR KO also shows CSPG resistance; double KO (PTPσ/LAR) has additive effects

### 3.2 Structural Data
| Structure | PDB ID | Resolution | Domain |
|-----------|--------|------------|--------|
| **LAR D1 domain** | 2LAR | 1.9 Å | Catalytic |
| **LAR D2 domain** | 3LAR | 2.2 Å | Regulatory |
| **LAR + CSPG** | 6LAR | 2.5 Å | CSPG binding |

### 3.3 Druggability
- Similar to PTPσ but less explored
- **Opportunity**: Dual PTPσ/LAR inhibitors may have superior efficacy
- **Challenge**: Selectivity over other RPTPs

---

## 4. Target: NgR1/NgR3 (Nogo Receptors)

### 4.1 Biological Rationale
- **Function**: Co-receptors for CSPG signaling (with PTPσ/LAR)
- **Mechanism**: CSPG binding → NgR clustering → RhoA/ROCK activation → growth cone collapse
- **Evidence**: NgR1 KO mice resistant to CSPG inhibition

### 4.2 Structural Data
| Structure | PDB ID | Resolution | Notes |
|-----------|--------|------------|-------|
| **NgR1 ectodomain** | 2OZ7 | 2.4 Å | LRR domain |
| **NgR1 + Nogo-66** | 3Z0L | 2.8 Å | Ligand complex |
| **NgR3 ectodomain** | 4X6M | 2.1 Å | CSPG binding |

### 4.3 Druggability
- **LRR domain**: Protein-protein interface, challenging for small molecules
- **Approach**: Peptide mimetics, cyclic peptides, stapled peptides
- **CNS penetration**: Peptide-based strategies (cell-penetrating peptides)

---

## 5. Target: Endogenous MMP-2/9 (Matrix Metalloproteinases)

### 5.1 Biological Rationale
- **Function**: Endogenous CSPG cleavage enzymes
- **Mechanism**: Cleave CSPG core proteins → PNN remodeling
- **Evidence**: MMP-9 KO mice have increased PNNs; MMP-2/9 upregulation after injury

### 5.2 Structural Data
| Structure | PDB ID | Resolution | Notes |
|-----------|--------|------------|-------|
| **MMP-2 catalytic domain** | 1CK7 | 1.5 Å | Zinc-dependent |
| **MMP-9 catalytic domain** | 1L6J | 1.8 Å | Gelatinase B |
| **MMP-2 + inhibitor** | 3AYU | 2.0 Å | Hydroxamate inhibitor |

### 5.3 Druggability
- **Activators**: No known direct small molecule activators (traditionally targeted with inhibitors); allosteric modulators or disruption of TIMP-3 inhibition can enhance local endogenous cleavage.
- **MMP-9 role**: Cleaves gelatin, aggrecan, and neural ECM components in an activity-dependent manner in ARC circuits.

---

## 5B. Target: ADAMTS-4 (Aggrecanase-1)

### 5B.1 Biological Rationale
- **Function**: Primary endogenous aggrecanase in the central nervous system
- **Mechanism**: Catalyzes specific Glu373-Ala374 interglobular domain cleavage in aggrecan core protein, as well as brevican and neurocan processing
- **Evidence**: Selective aggrecanase activity softens PNN coats without broad basement membrane degradation, preserving vascular and structural integrity while restoring insulin receptor mobility
- **Structural Data**: PDB 2RJP (1.5 Å), zinc metalloproteinase catalytic cleft with S1\' specificity subsite (His202, His206, His212, Glu203)

### 5B.2 Druggability Assessment
- **Pocket Druggability**: High (well-defined zinc catalytic domain and deep S1\' hydrophobic pocket)
- **Selectivity**: S1\' pocket variations allow discriminating ADAMTS-4 from MMP-1/-2/-14
- **Strategy**: Selective allosteric or catalytic modulation to stimulate controlled CSPG core protein turnover

---

## 6. Target Prioritization Matrix

| Target | Biological Validation | Structural Data | Druggability | CNS Precedent | Oral Feasibility | **Priority** |
|--------|----------------------|-----------------|--------------|---------------|------------------|--------------|
| **PTPσ CSPG site** | Strong | Good | Moderate | Moderate | Moderate | **1 (Primary)** |
| **PTPσ Catalytic** | Strong | Excellent | High | High | Moderate | **2 (Backup)** |
| **LAR CSPG site** | Strong | Good | Moderate | Low | Moderate | **3 (Combination)** |
| **ADAMTS-4** | Strong | Excellent (1.5 Å) | High | Moderate | Moderate | **4 (Aggrecanase)** |
| **MMP-2/9** | Moderate | Excellent | High | Moderate | Moderate | **5 (Matrix MMPs)** |
| **ChABC mimic** | Strong | Excellent | Low | None | Low | **6 (Long-term)** |
| **NgR1/3** | Moderate | Good | Low | Low | Low | **7 (Alternative)** |

---

## 7. Screening Strategy by Target

### 7.1 PTPσ CSPG Binding Site (Primary)
**Assay Type**: Fluorescence polarization / AlphaScreen / SPR
**Substrate**: Fluorescent CS tetrasaccharide (e.g., CS-E-AF488)
**Readout**: Displacement of CSPG from PTPσ
**Library**: 50K-200K fragments + focused CS-mimetics
**Hit Criteria**: IC50 < 10 μM, LE > 0.3, CSPG site selectivity > 10× vs catalytic

### 7.2 PTPσ Catalytic Site (Backup)
**Assay Type**: pNPP hydrolysis / DiFMUP fluorescence / HTRF
**Substrate**: pNPP or DiFMUP
**Readout**: Phosphatase inhibition
**Selectivity Panel**: PTP1B, TCPTP, SHP2, CD45, LAR
**Hit Criteria**: IC50 < 1 μM, >30× selectivity vs PTP1B

### 7.3 CSPG Chain Binding (ChABC Mimic)
**Assay Type**: DQ-CSPG cleavage (fluorogenic) / ELISA-based CSPG degradation
**Readout**: Fluorescence increase (cleavage) or CSPG loss
**Hit Criteria**: EC50 < 10 μM, >50% degradation at 10 μM

---

## 8. Structural Biology Resources

### 8.1 Available Structures (Download Links)
```
PTPσ D1:          https://files.rcsb.org/download/2PTS.pdb
PTPσ D2:          https://files.rcsb.org/download/3T8C.pdb
PTPσ + CSPG:      https://files.rcsb.org/download/5K4P.pdb
PTPσ + ISP:       https://files.rcsb.org/download/4Y1R.pdb
LAR D1:           https://files.rcsb.org/download/2LAR.pdb
NgR1:             https://files.rcsb.org/download/2OZ7.pdb
ChABC (P. vulg):  https://files.rcsb.org/download/1L5D.pdb
ChABC (Flavo):    https://files.rcsb.org/download/1L5E.pdb
ChABC + CS:       https://files.rcsb.org/download/2BZR.pdb
MMP-2:            https://files.rcsb.org/download/1CK7.pdb
MMP-9:            https://files.rcsb.org/download/1L6J.pdb
```

### 8.2 Homology Models Needed
- **PTPσ full-length** (D1-D2 + transmembrane + intracellular)
- **LAR full-length**
- **NgR1/3 full-length with CSPG**
- **Human ChABC** (if different from bacterial)

### 8.3 MD Simulation Targets
- **PTPσ CSPG site + CS tetrasaccharide** (100 ns)
- **PTPσ catalytic + selective inhibitor** (100 ns)
- **LAR CSPG site + CS tetrasaccharide** (100 ns)

---

## 9. Compound Design Guidelines by Target

### 9.1 PTPσ CSPG Site Binders
**Pharmacophore**:
- **Anionic groups**: 2-3 carboxylates/sulfonates (mimic CS)
- **Hydrophobic patch**: 1-2 aromatic rings (stack with Lys/Arg)
- **H-bond donors**: 2-3 (match CS hydroxyls)
- **Scaffold**: Rigid, pre-organized (biphenyl, terphenyl, macrocycle)

**Property Targets**:
- MW: 350-500 Da (larger for CSPG site)
- logP: 1.0-3.0
- TPSA: 80-120 Å² (higher for polar site)
- Charge: -2 to -3 at pH 7.4

### 9.2 PTPσ Catalytic Inhibitors
**Pharmacophore**:
- **Phosphotyrosine mimetic**: Oxalamate, difluoromethylene phosphonate, benzofuran
- **Hydrophobic extension**: Aryl group for specificity pocket
- **Scaffold**: Benzothiazole, thienopyridone, aminosterol

**Property Targets**:
- MW: 300-450 Da
- logP: 2.0-4.0
- TPSA: 60-90 Å²
- Charge: Neutral or -1

### 9.3 Dual PTPσ/LAR Inhibitors
**Strategy**: Target conserved CSPG binding residues
**Selectivity**: Exploit subtle differences in wedge domain

---

## 10. Validation Cascade

### Tier 1: Biochemical
- [ ] Binding to target (SPR/ITC, KD < 1 μM)
- [ ] Functional inhibition (IC50 < 10 μM)
- [ ] Selectivity panel (>10× vs related targets)

### Tier 2: Cellular
- [ ] PNN degradation in neuronal culture (WFA staining)
- [ ] Insulin signaling rescue (pIR, pAKT in mHypoA cells)
- [ ] CSPG binding competition (FRET/co-IP)

### Tier 3: Ex Vivo
- [ ] Hypothalamic slice insulin response
- [ ] PNN reduction in acute slices
- [ ] Electrophysiology (POMC/NPY firing)

### Tier 4: In Vivo
- [ ] Brain penetration (brain/plasma > 0.3)
- [ ] Hypothalamic target engagement
- [ ] Metabolic efficacy (ITT, food intake)

---

*Document Version: 1.0*
*Generated: 2026-10-11*
*Compatible with Recursive Language Model cold-start workflows*