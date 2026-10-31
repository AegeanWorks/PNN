# Oral Bioavailability Design Principles for CNS-Targeted Agents

## Overview
This document defines the physicochemical and formulation strategies for achieving oral bioavailability in BBB-penetrant small molecules.

---

## 1. Absorption & Bioavailability Fundamentals

### 1.1 Oral Absorption Mechanisms

| Mechanism | Description | Fraction of Drugs |
|-----------|-------------|-------------------|
| **Passive transcellular diffusion** | Lipid bilayer partitioning | ~60% of oral drugs |
| **P-glycoprotein efflux** | Active transport out of enterocytes | ~30% substrates |
| **Carrier-mediated transport** | LAT1, MCT1, GLUT1, PEPT1 | ~10% of drugs |
| **Paracellular diffusion** | Tight junction passage | <5% of drugs, limited by MW |

### 1.2 Key Absorption Parameters

| Parameter | Optimal Range | Impact on F (Bioavailability) |
|-----------|--------------|-------------------------------|
| **Solubility** | > 10 μM (BCS Class 1-2) | Primary rate limiter |
| **Permeability (Pe)** | > 10⁻⁶ cm/s (MDCK) | Primary rate limiter |
| **logD7.4** | 1-3 | Balances solubility + permeability |
| **pKa (basic)** | 7.5-8.5 | Affects dissolution + membrane partitioning |
| **Dissolution rate** | > 70% in 15 min (NBP) | Critical for BCS Class 2-4 |
| **First-pass metabolism** | < 50% loss | Key for high-extraction drugs |

### 1.3 Bioavailability (F) Equation
**F = f × Fa × Fg × Fh**

Where:
- **f**: Fraction of dose reaching absorption site (gastric emptying)
- **Fa**: Fraction absorbed (permeability × solubility limited)
- **Fg**: Fraction escaping gut wall metabolism (CYP3A4 predominant)
- **Fh**: Fraction escaping hepatic metabolism (CYP3A4, CYP2C9, CYP2C19)

**Typical F ranges**:
- **BCS Class 1** (high solubility, high permeability): 80-100%
- **BCS Class 2** (high solubility, low permeability): 20-50%
- **BCS Class 3** (low solubility, high permeability): 50-80%
- **BCS Class 4** (low solubility, low permeability): < 10%

---

## 2. Physicochemical Optimization for Oral Absorption

### 2.1 logD7.4 Target Ranges
| logD7.4 Range | Absorption Profile | CNS Penetration Profile |
|---------------|-------------------|------------------------|
| **< 1.0** | Solubility-limited | Poor (unless P-gp substrate) |
| **1.0 - 2.0** | Balanced | Good (optimal) |
| **2.0 - 3.0** | Permeability-limited | Good (optimal for CNS) |
| **3.0 - 4.0** | Solubility good, permeability variable | Moderate |
| **> 4.0** | Solubility-limited, formulation needed | Poor (P-gp substrate risk) |

**Rule of thumb**: logD7.4 = 2.0-3.0 optimizes both oral absorption AND CNS penetration for neutral/weakly basic compounds.

### 2.2 pKa Optimization for Basic Compounds

**Goal**: pKa 7.5-8.5 for optimal oral + CNS performance

**Rationale**:
- **Intestinal pH (6.5)**: ~10-50% ionized (depending on exact pKa)
- **Blood pH (7.4)**: ~50-80% ionized
- **Brain interstitial pH (7.2)**: ~30-60% ionized
- **Brain trapping**: Once ionized, limited efflux → accumulation

**Design Strategy**:
```
If pKa = 8.0 (optimal):
- At intestinal pH 6.5: 24% ionized → good dissolution
- At blood pH 7.4: 61% ionized → moderate solubility + some permeability
- At brain pH 7.2: 56% ionized → some brain retention

If pKa = 9.5 (too high):
- Intestinal ionization: 76% → poor dissolution
- Blood ionization: 86% → very poor permeability

If pKa = 7.0 (too low):
- Intestinal ionization: 9% → good dissolution but less absorption window
- Blood ionization: 44% → good permeability but less retention
```

### 2.3 pKa Optimization for Acidic Compounds

**Goal**: pKa 3.5-5.0 for optimal oral + CNS performance

**Rationale**:
- **Intestinal pH (6.5)**: Mostly unionized (favorable for absorption)
- **Blood pH (7.4)**: Mostly unionized (favorable for passive diffusion)
- **Brain pH (7.2)**: Mostly unionized (good penetration)
- **Advantage**: Neutral species at all relevant pHs

**Design Strategy**:
```
If pKa = 4.5 (optimal):
- At intestinal pH 6.5: 93% unionized → excellent absorption
- At blood pH 7.4: 97% unionized → excellent permeability
- At brain pH 7.2: 96% unionized → good brain penetration

If pKa = 3.0 (too low):
- Already mostly ionized at intestinal pH → solubility may help but absorption suffers

If pKa = 5.5 (too high):
- Some ionization at intestinal pH → reduced absorption window
```

### 2.4 Solubility-Permeability Balance

**The "Absorption Window" Concept**:
```
Optimal compound exists in the "absorption window" where:
- Sufficient solubility for dissolution in GI fluids
- Sufficient permeability for transcellular passage
- Neither property dominates to limit overall absorption
```

**Quantitative Model** (from Dressman & Fell):
```
Fa = 1 - exp(-Pe × t)
where Pe = permeability coefficient, t = transit time

Solubility requirement: 
Cmax,dissolved > Cmax,absorbed
or: Dose/Vd < Solubility × Fraction absorbed
```

**Practical Targets**:
- **BCS Class 1**: Solubility > 250 μg/mL (dose-adjusted) + Pe > 10⁻⁵ cm/s
- **Optimal for CNS**: Solubility > 50 μM + Pe > 5 × 10⁻⁶ cm/s

---

## 3. Formulation Strategies

### 3.1 pH-Dependent Solubilization

**Weak Base Strategies** (pKa 7.5-8.5):
- **Acidic co-formers**: Fumaric acid, maleic acid, succinic acid
- **Salt formation**: Hydrochloride, phosphate, mesylate
- **Co-solvents**: PEG 400, propylene glycol (≤20%)
- **Solid dispersions**: PVP K30, HPMCAS (pH-dependent release)

**Weak Acid Strategies** (pKa 3.5-5.0):
- **Basic co-formers**: Triethylamine, diethylamine
- **Cationic salts**: Meglumine, lysine, arginine
- **Co-crystals**: With nicotinamide, nicotinic acid

**Strong Acid/Bases**: 
- Limited use due to formulation stability
- Require enteric coating for GI protection

### 3.2 Solid Dispersion Technologies

| Technology | Carrier | Drug Loading | Release Profile | Maturity |
|------------|---------|--------------|-----------------|----------|
| **Melting (fusion)** | PVP K30, PEG 6000 | 20-50% | Rapid (5-15 min) | High |
| **Solvent evaporation** | PVP K30, HPMC | 10-30% | Moderate (15-30 min) | High |
| **Spray drying** | PVP K30, HPMCAS | 10-25% | Tailorable | High |
| **Melt electrospinning** | PVP, Eudragit | 10-20% | Sustained | Emerging |

**HPMCAS**: Hydroxypropyl methylcellulose acetate succinate
- **pH-dependent**: Soluble at pH < 5.5, gelates at pH > 6.5
- **Use case**: Protect acid labors, promote gastric dissolution of bases

### 3.3 Co-Crystal Formation

**Advantages**:
- Improves solubility 2-10× without adding ionizable groups
- Maintains molecular integrity (no covalent change)
- Can improve both solubility AND permeability

**Common Co-formers**:
- Nicotinamide (NA)
- Nicotinic acid (NACA)
- Pyrrolidone-2-one (P2O)
- Saccharin
- Theobromine

**Formation methods**:
- Grinding (neat mechanochemistry)
- Slurry anti-solvent
- Evaporative precipitation
- Melt co-crystallization

### 3.3 Nanoparticle Formulations

**Nanocrystal Approach**:
- **Sub-200 nm particles** of pure drug (no carrier)
- **Stabilizers**: PVP, Poloxamer 188, citrate
- **Dissolution enhancement**: Noyes-Whitney equation (surface area ↑ 100×)
- **Oral bioavailability**: 2-5× improvement typical

**Lipid-based nanoparticles**:
- **Solid lipid nanoparticles (SLN)**
- **Nanoemulsions** (oil + surfactant + drug)
- **Self-emulsifying drug delivery systems (SEDDS)**
- **Benefits**: Bypass P-gp, lymphatic transport, food effect mitigation

### 3.4 Enteric Coating

**When needed**:
- **Drug degradation** in stomach (pH < 4)
- **Gastric irritation** (NSAIDs, potassium chloride)
- **Targeting lower GI** (colon-specific delivery)

**Coating materials**:
- **Shellac** (pH 5.5 threshold)
- **Methacrylic acid copolymers** (Eudragit L100-55, S100 - pH 6)
- **Polyvinyl acetate phthalate** (Eudragit E - pH 5.5)

**Limitation**: Reduces absorption window (must survive stomach → dissolve in intestine)

---

## 4. First-Pass Metabolism Mitigation

### 4.1 CYP3A4 Interaction Risk

**High-risk structures**:
- **Grapfruit juice effect**: Furanocoumarins, certain alkaloids
- **Aryl heterocycles**: Pyridine, imidazole, triazole rings
- **N-containing scaffolds**: Secondary/tertiary amines

**Mitigation strategies**:
1. **Reduce basic amine basicity**: Lower pKa → less CYP3A4 binding
2. **Fluorine substitution**: Blocks oxidative metabolism at specific positions
3. **Cyclization**: Reduces flexibility + blocks metabolic soft spots
4. **Isosteric replacement**: Replace amide with oxadiazole, thiadiazole

**Low-risk structures**:
- **Perfluoroalkyl**: Highly resistant to metabolism
- **Saturated aliphatics**: Limited oxidation sites
- **Aryl fluorides**: Block adjacent metabolism

### 4.2 Gut Microbiome Metabolism

**Drugs affected by gut flora**:
- **Cardiac glycosides**: Digoxin (3β-hydroxysteroid oxidase)
- **Azo drugs**: Sulfasalazine → sulfapyridine + mesalamine
- **Isoflavones**: Equol production (varies by microbiome)
- **Bile acid conjugates**: Deconjugation by β-glucuronidases

**Mitigation**:
- **Antibiotic co-administration**: Not practical long-term
- **Probiotic co-administration**: May stabilize metabolism
- **Formulation protection**: Enteric coating, delayed release
- **Avoidance**: Don't design substrates for common gut enzymes

### 4.3 Hepatic Extraction Ratio (EH)

| EH Range | Hepatic Clearance | Fh (Escaping) | Design Priority |
|----------|------------------|---------------|-----------------|
| **EH < 0.3** | Low | > 70% | Standard oral dosing |
| **0.3 < EH < 0.7** | Moderate | 30-70% | Dose adjustment, modified release |
| **EH > 0.7** | High | < 30% | Alternative routes, prodrugs |

**Fh calculation**: Fh = e^(-EH) or Fh = 1/(1 + EH/Qh × fu)
Where Qh = hepatic blood flow (~90 L/h human), fu = unbound fraction

**Strategies for high EH drugs**:
1. **Reduce first-pass**: Sublingual, rectal, pulmonary routes
2. **Prodrug strategy**: Bypass gut/liver activation
3. **Inhibitor co-administration**: Low-dose CYP3A4 inhibitor (controversial)
4. **Formulation**: Lipid-based → lymphatic bypass → reduced first-pass

---

## 5. Food Effects on Absorption

### 5.1 Fed vs Fasted State

| Food Effect | Fed State | Fasted State | Design Implication |
|-------------|-----------|--------------|-------------------|
| **Gastric emptying** | Slower (2-3h) | Faster (1h) | Adjust Tmax timing |
| **Bile secretion** | Increased | Baseline | Improves solubility (lipophilic drugs) |
| **Peristalsis** | Increased | Baseline | Variable transit time |
| **pH changes** | Duodenal pH ↑ | More acidic | Affects ionizable drug solubility |
| **Plasma protein binding** | Variable | Baseline | Affects free drug concentration |

### 5.2 Food Effect Classification

| Food Effect Type | Example | Impact on F | Dosing Recommendation |
|-----------------|---------|---------------|----------------------|
| **No food effect** | Most BCS Class 1 | Minimal (±10%) | Take consistently, any time |
| **Decreased F (fed)** | Some BCS Class 2 | 30-50% reduction | Take fasted (1h before or 2h after meal) |
| **Increased F (fed)** | Lipophilic drugs, BCS Class 1 with fat | 20-50% increase | Take with meal (enhances absorption) |
| **Tmax delayed (fed)** | Many drugs | No F change | Adjust timing, consistent dosing |

### 5.3 Food Effect Mitigation Strategies

**For decreased F (fed)**:
1. **Fasted dosing**: 1h before or 2h after meal
2. **Formulation**: Take with moderate fat (may enhance some drugs)
3. **Coating**: Enteric if stomach degradation

**For increased F (fed)**:
1. **Consistent fed dosing**: Same meal composition each time
2. **Formulation**: Modify to reduce food sensitivity
3. **Timing**: Establish routine (always with breakfast, etc.)

**General recommendation**: 
- Define food effect in Phase 1
- Establish fasting/fed dosing in Phase 2-3
- Label consistently in IND/enabling studies

---

## 6. Project-Specific Oral Bioavailability Criteria

### 6.1 Phase 1-2 Screening Criteria

| Parameter | Minimum | Target | Measurement Method |
|-----------|---------|--------|-------------------|
| **Solubility** | > 10 μM | > 50 μM | Shake-flask, OECD 111 |
| **Permeability (Pe)** | > 10⁻⁶ cm/s | > 5 × 10⁻⁶ cm/s | MDCK, PAMPA-BBB |
| **logD7.4** | 1.0-4.0 | 2.0-3.0 | XlogP3-AA, ALOGPS |
| **pKa (basic)** | 7.0-9.0 | 7.5-8.5 | UV titration, software |
| **% Oral F (predicted)** | > 10% | > 20% | GastroPlus, Simcyp |
| **Food effect** | N/A | Documented | In vitro-in vivo correlation |

### 6.2 Formulation Development Priority

| Compound Profile | Formulation Strategy | Timeline |
|------------------|--------------------|----------|
| **BCS Class 1** (high sol + perm) | None required | Standard capsule |
| **BCS Class 2** (high sol, low perm) | Solid dispersion, nanoparticle | Month 3-6 |
| **BCS Class 3** (low sol, high perm) | Co-crystal, salt formation | Month 2-4 |
| **BCS Class 4** (low sol + low perm) | Nanoparticle + permeation enhancer | Month 6-12 |

### 6.3 Absorption Enhancement Strategies (If Needed)

| Strategy | Mechanism | Example Agents | Maturity |
|----------|-----------|----------------|----------|
| **P-gp inhibition** | Block efflux at intestine | Verapamil, tariquidar | Research |
| **Bile salt co-administration** | Micelle formation, lymphatic transport | Sodium taurocholate | Clinical |
| **Permeation enhancers** | Disrupt tight junctions | Sodium lauryl sulfate, CAPB | Cosmetics/OTC |
| **Mucoadhesive** | Prolong GI residence | Chitosan, Carbomer | Pharmaceutical |
| **Lymphatic targeting** | Bypass portal vein absorption | Medium-chain triglycerides, SLN | Clinical |

---

## 7. Cold-Start Protocol for Oral Bioavailability

### 7.1 When Starting a New Session
1. **Load solubility/permeability prediction scripts**
2. **Set up compound library** with standardized properties
3. **Run initial absorption filter**: logD 1-4, MW 300-450, TPSA < 90
4. **Calculate predicted Fa**: using MDCK/PAMPA data + solubility

### 7.2 Standardized Property Table
```
# All compounds tracked with:
# - SMILES string
# - MW, logP, TPSA, HBD, HBA, RB
# - pKa (acidic and basic)
# - logD7.4 (predicted)
# - Solubility (μM, pH 6.5)
# - Pe (MDCK, × 10⁻⁶ cm/s)
# - Predicted Fa (fraction absorbed)
# - Predicted F (oral bioavailability)
# - YI score (CNS multiplicity)
```

### 7.3 Initial Filter Script
```python
def oral_absorption_filter(smiles, pKa_base=None, pKa_acid=None):
    """Apply oral bioavailability filters"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rb = Descriptors.NumRotatableBonds(mol)
    
    # logD7.4 prediction (simplified)
    if pKa_base:
        # Henderson-Hasselbalch at pH 7.4
        fraction_ionized = 1 / (1 + 10**(pKa_base - 7.4))
        logd74 = logp - np.log10(1 + 10**(pKa_base - 7.4))
    else:
        logd74 = logp  # assume neutral
        fraction_ionized = 0
    
    # Violations
    violations = 0
    if not (300 <= mw <= 450):
        violations += 1
    if not (1.0 <= logd74 <= 4.0):
        violations += 1
    if tpsa > 90:
        violations += 1
    if hbd > 3:
        violations += 1
    if hba > 7:
        violations += 1
    if rb > 7:
        violations += 1
    
    # Fa prediction (empirical)
    fa = 0.8 if (logd74 >= 1.5 and tpsa <= 70 and hbd <= 2) else 0.5
    
    return violations <= 2, {
        "MW": f"{mw:.1f}", 
        "logP": f"{logp:.2f}", 
        "logD7.4": f"{logd74:.2f}",
        "TPSA": f"{tpsa:.1f}",
        "HBD": hbd, 
        "HBA": hba, 
        "RB": rb,
        "f_ionized": f"{fraction_ionized:.2f}",
        "predicted_Fa": fa
    }
```

---

## 8. Project-Specific TPP (Target Product Profile) for Oral F

| Parameter | Minimum Spec | Target Spec | Rationale |
|-----------|--------------|-------------|-----------|
| **MW** | 300-450 Da | 350 ± 50 Da | Balance absorption + molecular recognition |
| **logD7.4** | 1.0-4.0 | 2.0-3.0 | Optimizes solubility + permeability + CNS penetration |
| **TPSA** | < 90 Å² | < 70 Å² | Reduces desolvation penalty at membrane |
| **HBD** | ≤ 3 | ≤ 2 | Minimize permeability barrier |
| **HBA** | ≤ 7 | ≤ 5 | Balance permeability + solubility |
| **pKa (basic)** | 7.0-9.0 | 7.5-8.5 | Intestinal dissolution + brain retention |
| **pKa (acidic)** | 3.5-5.5 | 4.0-5.0 | Intestinal unionized fraction |
| **logBB (predicted)** | > -1.0 | > -0.5 | Basic BBB penetration criterion |
| **Predicted Fa** | > 10% | > 20% | Oral absorption estimate |
| **Predicted F** | > 5% | > 20% | Complete oral bioavailability |
| **Efflux ratio** | < 15 | < 5 | P-gp liability assessment |
| **Solubility (pH 6.5)** | > 10 μM | > 50 μM | Dissolution-limited absorption |
| **Dissolution (NBP, 15 min)** | > 70% | > 85% | Standard compendial test |

**Go/No-Go**: Advance compound if ≥70% of target specs met

---

*Document Version: 1.0*
*Generated: 2026-10-11*
*Compatible with Recursive Language Model cold-start workflows*