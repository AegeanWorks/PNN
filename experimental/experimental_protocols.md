# Experimental Validation Protocols: PNN Softening & Insulin Signaling

## Overview
This document defines the in vitro and in vivo experimental protocols for validating PNN-softening agents and their effects on hypothalamic insulin signaling.

---

## 1. In Vitro Assays

### 1.1 PNN Degradation Assay (Fluorogenic Substrate)

#### 1.1.1 Principle
- **Substrate**: DQ®-conjugated chondroitin sulfate (fluorescence quenched until cleavage)
- **Readout**: Fluorescence increase (excitation/emission 485/535 nm)
- **Principle**: Cleavage of CS chains by test compound → fluorescence dequenching

#### 1.1.2 Protocol

| Step | Details |
|------|---------|
| **1. Coating** | Add 50 μL DQ-CS (10 μg/mL) to black wall/clear bottom 96-well plate |
| **2. Incubation** | 37°C, 30 min (adsorption) |
| **3. Washing** | 3× with assay buffer (PBS + 0.1% BSA) |
| **4. Compound addition** | 100 μL compound dilutions (10 μM serial dilution) |
| **5. Positive control** | ChABC (1 U/mL) in separate wells |
| **6. Negative control** | Vehicle (0.1% DMSO) |
| **7. incubation** | 37°C, 2-24 h (time course) |
| **8. Reading** | Fluorescence plate reader (Ex/Em 485/535 nm) |

#### 1.1.3 Data Analysis
| Metric | Calculation |
|--------|-------------|
| **EC50** | Nonlinear regression, % degradation vs [compound] |
| **% Max degradation** | (F_sample - F_vehicle) / (F_ChABC - F_vehicle) × 100 |
| **Time to EC50** | Time point where 50% max degradation achieved |
| **Specificity** | Ratio: (compound / ChABC) at same % degradation |

#### 1.1.96 Expected Results
| Compound Type | Expected EC50 | Expected Max % |
|---------------|---------------|----------------|
| **ChABC** | ~0.01 U/mL | 100% |
| **PTPσ antagonist** | 1-100 μM | 20-50% (partial) |
| **LAR antagonist** | 1-50 μM | 15-40% (partial) |
| **Non-binders** | > 100 μM | < 10% |
| **Cytotoxic** | N/A | Variable (counter-screen) |

---

### 1.2 BBB Permeability Assays

#### 1.2.1 MDCK-MDR1 (Monkey Kidney Cells, MDR1-expressing)

##### 1.2.1.1 Principle
- **Monolayer integrity**: TEER (transendothelial electrical resistance)
- **Efflux ratio**: Pe_apical→basal / Pe_basolateral→apical
- **Classification**: ER < 2 = low efflux, ER > 5 = high efflux (P-gp substrate)

##### 1.2.1.2 Protocol

| Step | Details |
|------|---------|
| **1. Coating** | Collagen IV (10 μg/mL), 2 h at 37°C |
| **2. Seeding** | MDCK-MDR1 cells, 5 × 10⁴ cells/well, 24-well Transwell |
| **3. Differentiation** | 5-7 days, change medium daily |
| **4. TEER measurement** | EVOM2 voltohmmeter, ≥300 Ω·cm2 required |
| **5. Compound preparation** | Test compounds in transport medium (HBSS + 2% FBS) |
| **6. apical→basal** | Add 500 μL to apical chamber, 100 μM concentration |
| **7. basolateral→apical** | Add 500 μL to basolateral, 100 μM concentration |
| **8. Time points** | 30 min, 1 h, 2 h, 4 h |
| **9. Sampling** | 100 μL from receiving chamber, replace with fresh medium |
| **10. Analysis** | LC-MS/MS or UV-Vis (if chromophoric) |
| **11. Calculations** | Pe = (dQ/dt) / (A × Co), ER = Pe_AB / Pe_BA |

##### 1.2.1.3 Data Analysis
| Metric | Formula | Threshold |
|--------|---------|-----------|
| **Pe_apical→basal** | μmol/cm²/s | > 5 × 10⁻⁶ = good permeability |
| **Pe_basolateral→apical** | μmol/cm²/s | < 1 × 10⁻⁶ desirable |
| **Efflux Ratio (ER)** | Pe_BA / Pe_AB | < 2 = low efflux, > 5 = P-gp substrate |
| **Permeability classification** | ER ranges | See Section 2.6 |

#### 1.2.2 PAMPA-BBB (Parallel Artificial Membrane Permeability Assay)

##### 1.2.2.1 Principle
- **Artificial membrane**: Lipid-coated filter (POPC/POPE 3:1)
- **Readout**: Compound flux from donor to acceptor chamber
- **Advantage**: No living cells, no P-gp efflux (passive only)

##### 1.2.2.2 Protocol
| Step | Details |
|------|---------|
| **1. Membrane preparation** | Lipid film (POPC:POPE 3:1) in donor chamber |
| **2. Compound addition** | 100 μM in donor chamber (pH 7.4 buffer) |
| **3. Incubation** | 37°C, 2-4 h (linear range) |
| **3.5. Collection** | Sample acceptor chamber at intervals |
| **4. Analysis** | UV-Vis, fluorescence, or LC-MS |
| **5. Permeability calculation** | Pe = (dQ/dt) / (A × Co × Kp) |

##### 1.2.2.3 Data Analysis
| Metric | Formula | Target |
|--------|---------|--------|
| **logBB (predicted)** | From Pe and other properties | > -0.5 |
| **Pe (PAMPA)** | cm/s | > 1 × 10⁻⁶ = good |
| **Passive permeability** | No energy dependence | Confirms passive diffusion |

---

### 1.3 Insulin Signaling Assay (Neuronal Cells)

#### 1.3.1 Cell Type: mHypoA-60/3 (Mouse hypothalamic neuronal line)

##### 1.3.1.1 Culture Conditions
| Parameter | Details |
|-----------|---------|
| **Medium** | DMEM/F-12 + 10% FBS + 1% pen-strep |
| **Coating** | Poly-D-lysine (10 μg/mL), 2 h |
| **Passaging** | 1:4 ratio every 7-10 days |
| **Confluency** | 70-80% for assays |
| **Serum starvation** | 2 h before insulin treatment |

##### 1.3.1.2 Treatment Protocol

| Parameter | Details |
|-----------|---------|
| **Compound dilution** | 1, 10, 100, 1000 nM in serum-free medium |
| **Pre-treatment** | 1 h compound before insulin |
| **Insulin stimulation** | 100 nM human insulin, 10 min |
| **Controls** | - Vehicle (0.1% DMSO)<br>- Insulin only (no compound)<br>- ChABC positive control |
| **Replicates** | n = 3-4 per condition |
| **Lysis** | RIPA buffer + phosphatase inhibitors |

##### 1.3.1.3 Readouts (Western Blot)

| Target | Antibody | Expected Effect of PNN Softening |
|--------|----------|----------------------------------|
| **pIR-Y1161** | Anti-phospho-IR (Tyr1161) | ↑ 2-3 fold |
| **pIR-Y1162** | Anti-phospho-IR (Tyr1162) | ↑ 2-3 fold |
| **pAKT-S473** | Anti-phospho-AKT (Ser473) | ↑ 1.5-2 fold |
| **pAKT-T308** | Anti-phospho-AKT (Thr308) | ↑ 1.5-2 fold |
| **pFOXO1-T24** | Anti-phospho-FOXO1 (Thr24) | ↓ 40-50% (excluded from nucleus) |
| **Total IR** | Anti-IR (C-20) | No change (loading control) |
| **Total AKT** | Anti-AKT (N-19) | No change |
| **β-actin** | Anti-β-actin (AC-15) | Loading control |

##### 1.3.1.4 Quantification
| Step | Method |
|------|--------|
| **1. Band detection** | Li-COR Odyssey (fluorescent Western) |
| **2. Integrated density** | ImageJ "Analyze Gel" |
| **3. Normalization** | pIR/IR or pAKT/AKT, then /β-actin |
| **4. Fold-change** | Treated / Vehicle |
| **5. Statistics** | One-way ANOVA + Tukey post-hoc |

##### 1.3.1.5 Expected Dose-Response
```
Compound → PNN softening → Insulin receptor accessibility → pIR ↑
EC50 for pIR increase typically: 1-100 μM (depends on compound potency)
Max effect: 2-5 fold increase over vehicle (ChABC reference)
```

---

### 1.4 Cytotoxicity Assay (Counter-Screen)

#### 1.4.1 CellTiter-Glo 2.0 (ATP-based viability)

##### 1.4.1.1 Protocol
| Step | Details |
|------|---------|
| **1. Seeding** | 5,000 cells/well in 96-well white wall plate |
| **2. Treatment** | Test compounds, 24 h, 48 h, 72 h |
| **3. Incubation** | Add CellTiter-Glo reagent, 15 min room temperature |
| **4. Reading** | Luminescence (GLMax-Micro, integration 0.5 s) |
| **5. Standard curve** | Cell number standard (0-50,000 cells) |

##### 1.4.1.2 Data Analysis
| Metric | Calculation |
|--------|-------------|
| **% Viability** | (RLU_sample / RLU_vehicle) × 100 |
| **CC50** | Compound concentration causing 50% cell death |
| **Selectivity Index** | CC50 / EC50 (from PNN or signaling assay) |
| **Threshold for progression** | SI > 50 (preferable), > 10 (minimum) |

##### 1.4.1.2 Cytotoxicity Tiers
| Tier | CC50 Range | Action |
|------|------------|--------|
| **Safe** | > 100 μM | Proceed to in vivo |
| **Moderate** | 30-100 μM | Optimize structure, reduce toxicity |
| **Unsafe** | < 30 μM | Abandon or radically redesign |

---

## 2. Ex Vivo Protocols

### 2.1 Hypothalamic Slice Preparation

#### 2.1.1 Tissue Collection
| Parameter | Details |
|-----------|---------|
| **Animal** | Mouse (8-12 weeks), or rat |
| **Diet** | Chow or HFD (high-fat diet, 60% kcal from fat) |
| **Euthanasia** | CO2 followed by cervical dislocation |
| **Timing** | Same time of day (ZT4-ZT6, 4-6 h after lights on) |

#### 2.1.2 Slice Dissection
| Step | Details |
|------|---------|
| **1. Brain removal** | Rapid decapitation, immediate brain extraction |
| **2. Chilling** | Ice-cold slicing solution (in mM): 250 sucrose, 20 NaCl, 26 NaHCO3, 10 glucose, 1 CaCl2, 3 MgCl2, bubbled with 95% O2/5% CO2 |
| **3. ARC identification** | Anterior-posterior: -1.06 to -1.80 mm from Bregma |
| **4. Slicing** | Vibratome, 300-μm thickness, horizontal or coronal |
| **5. Recovery** | 37°C, 30 min in ACSF (in mM): 124 NaCl, 26 NaHCO3, 10 glucose, 2 CaCl2, 1.25 NaH2PO4, 2 MgCl2, 2.5 KCl |

#### 2.1.3 Treatment Protocol

| Parameter | Details |
|-----------|---------|
| **Compound** | Test compound, 10 μM (or dose-range) |
| **Vehicle** | 0.1% DMSO in ACSF |
| **ChABC** | Positive control, 0.1 U/mL |
| **Incubation** | 37°C, 2 h in ACSF |
| **Recovery** | Return to incubator, 30 min before assay |

#### 2.1.4 Readouts

| Readout | Method | Expected |
|---------|--------|----------|
| **WFA lectin staining** | Immunofluorescence, confocal | ↓ PNN area/volume |
| **pIR IHC** | Immunohistochemistry | ↑ IR phosphorylation |
| **pAKT IHC** | Immunohistochemistry | ↑ AKT phosphorylation |
| **Electrophysiology** | Whole-cell patch-clamp | ↑ POMC/NPY firing frequency |

---

### 2.2 WFA Lectin Quantification (Ex Vivo)

#### 2.2.1 Staining Protocol
| Step | Details |
|------|---------|
| **1. Fixation** | 4% PFA, 15 min RT |
| **2. Permeabilization** | 0.3% Triton X-100, 10 min |
| **3. Blocking** | 5% normal donkey serum, 1 h RT |
| **4. Primary** | WFA-Biotin (1:500), overnight 4°C |
| **5. Secondary** | Streptavidin-AlexaFluor 488 (1:1000), 2 h RT |
| **6. Counterstain** | DAPI (1:1000), 10 min |
| **7. Mounting** | Antifade medium, coverslip |

#### 2.2.2 Imaging & Analysis
| Parameter | Details |
|-----------|---------|
| **Microscope** | Confocal (63× oil objective) |
| **Settings** | Identical across all groups/conditions |
| **Region** | Arcuate nucleus (ARC) |
| **Fields** | 3-5 non-overlapping fields per slice |
| **Analysis** | ImageJ WFA puncta macro |

| Metric | Calculation |
|--------|-------------|
| **WFA+ area fraction** | (WFA+ area / total ARC area) × 100 |
| **Puncta per neuron** | (WFA+ puncta / NeuN+ cells) |
| **Mean puncta intensity** | Integrated density / puncta count |

#### 2.2.3 Statistical Power
| Parameter | Value |
|-----------|-------|
| **n (slices)** | 4-6 slices per mouse |
| **n (mice)** | 5-8 per group |
| **Effect size** | Based on ChABC: ~25% PNN reduction |
| **Power (1-β)** | 80% |
| **α** | 0.05 |
| **Test** | One-way ANOVA + Tukey |

---

## 3. In Vivo Protocols

### 3.1 Oral Gavage PK Study

#### 3.1.1 Experimental Design

| Parameter | Details |
|-----------|---------|
| **Animal** | Mouse (CD-1 or C57BL/6, male, 8-12 weeks) |
| **Dose** | 1, 3, 10, 30 mg/kg PO (oral gavage) |
| **Volume** | 10 μL/g body weight (max 200 μL) |
| **Fasting** | 4 h prior (water ad lib) |
| **Time points** | 0.5, 1, 2, 4, 8, 24 h post-dose |
| **Sample size** | n = 3-4 per time point per dose |
| **Anesthesia** | Avoid if possible (affects PK) |
| **Bleed method** | Retro-orbital (terminal) or tail vein (terminal) |

#### 3.1.2 Sample Collection

| Sample | Collection | Processing |
|--------|------------|------------|
| **Plasma** | Retro-orbital bleed | Heparinized tube, centrifuge 5 min @ 10,000 × g, plasma to -80°C |
| **Brain** | Terminal, separate cohort | Quick dissection (hypothalamus + rest), snap-freeze in liquid N2, store -80°C |
| **Brain regions** | Dissect ARC separately | May pool or analyze separately |

#### 3.1.3 PK Analysis

| Parameter | Formula/Method |
|-----------|----------------|
| **Cmax** | Maximum observed plasma concentration |
| **Tmax** | Time of Cmax |
| **AUC0-t** | Linear up-log down trapezoidal rule |
| **AUC0-∞** | AUC0-t + Ct/ke (Ct = last concentration, ke = terminal rate) |
| **λz** | Terminal rate constant (log-linear regression of last 3-4 points) |
| **t1/2** | ln2 / λz |
| **CL** | Dose / AUC0-∞ (clearance) |
| **Vd** | CL / ke (volume of distribution) |
| **F** | Oral bioavailability = (AUC_po × Dose_iv) / (AUC_iv × Dose_po) |
| **F_abs** | Fraction absorbed (relative to iv) |
| **Brain/plasma ratio** | C_brain / C_plasma (at Tmax or AUC ratio) |

#### 3.1.4 Data Processing (Python example)
```python
import pandas as pd
from scipy import integrate

# Load data
pk = pd.read_csv('mouse_pk.csv')

# Calculate AUC using trapezoidal rule
def trapezoidal_auc(time, conc):
    auc = 0
    for i in range(len(time)-1):
        auc += (time[i+1] - time[i]) * (conc[i] + conc[i+1]) / 2
    return auc

# Per animal AUC
pk['auc_trap'] = pk.groupby('animal').apply(
    lambda g: trapezoidal_auc(g['time'], g['conc']))

# Terminal elimination phase
def terminal_ke(time, conc):
    # Linear regression of last 3-4 points on log scale
    from numpy import polyfit, log
    n = min(4, len(time))
    t = time[-n:]
    c = conc[-n:]
    logc = log(c)
    ke = -polyfit(t, logc, 1)[0]  # slope = -ke
    return ke

pk['ke'] = pk.groupby('animal').apply(terminal_ke)
pk['t1_2'] = np.log(2) / pk['ke']
pk['auc_0_inf'] = pk['auc_trap'] + pk['conc'].iloc[-1] / pk['ke']
pk['cl'] = pk['dose'] / pk['auc_0_inf']
pk['vd'] = pk['cl'] / pk['ke']
pk['f'] = (pk['auc_po'] / pk['auc_iv']) * (pk['dose_iv'] / pk['dose_po'])
pk['brain_plasma'] = pk['brain_conc'] / pk['plasma_conc']
```

### 3.2 Insulin Tolerance Test (ITT)

#### 3.2.1 Protocol
| Step | Details |
|------|---------|
| **1. Baseline** | Fast mice 4 h, measure tail-bleed glucose (time 0) |
| **2. Insulin injection** | IP: 0.75 U/kg regular insulin |
| **3. Glucose measurement** | Tail-bleed at 0, 15, 30, 45, 60, 90, 120 min |
| **4. Glucometer** | Glucometer with animal strips |
| **5. Calculations** | % change from baseline, glucose AUC0-120 |

#### 3.2.2 Data Analysis
| Metric | Calculation |
|--------|-------------|
| **Glucose AUC** | Trapezoidal rule, 0-120 min |
| **% AUC reduction** | (AUC_vehicle - AUC_treated) / AUC_vehicle × 100 |
| **Time to normoglycemia** | Time for glucose < 150 mg/dL |
| **Overall effect** | Dose-dependent improvement |

#### 3.2.3 Success Criteria
| Metric | Threshold |
|--------|-----------|
| **% AUC reduction** | > 15% vs vehicle (dose-dependent) |
| **Statistical** | p < 0.05 (ANOVA + post-hoc) |
| **Consistency** | Reproducible across 2+ doses |

---

### 3.3 Food Intake & Body Weight Monitoring

#### 3.3.1 Protocol
| Parameter | Details |
|-----------|---------|
| **Housing** | Single-housed or pair-housed, automated feeding cages |
| **Baseline** | 7-day pre-treatment food intake measurement |
| **Treatment** | Daily oral gavage, 4 weeks minimum |
| **Food recording** | Gram consumed per mouse per day |
| **Body weight** | Digital scale, 3x/week (same time of day) |

#### 3.3.2 Data Analysis
| Metric | Calculation |
|--------|-------------|
| **% Food intake reduction** | (Baseline - Week 4) / Baseline × 100 |
| **Body weight change** | Week 4 - Baseline |
| **Weight loss %** | (Baseline - Week 4) / Baseline × 100 |
| **Food efficiency** | Weight change / Food intake |

#### 3.3.3 Expected Outcomes
| Outcome | Expected Range (optimistic) |
|---------|----------------------------|
| **Food intake reduction** | 15-25% at optimal dose |
| **Body weight reduction** | 10-20% over 4 weeks |
| **Dose-response** | Clear EC25/EC50 for metabolic effects |

---

### 3.4 Hypothalamic PNN Quantification (In Vivo)

#### 4.1 Tissue Collection
| Parameter | Details |
|-----------|---------|
| **Perfusion** | 4% PFA via cardiac catheter |
| **Post-fixation** | 4% PFA, 24 h at 4°C |
| **Cryoprotection** | 30% sucrose, 24 h at 4°C (or until sink) |
| **Sectioning** | Cryostat, 30-μm sections, ARC (bregma -1.06 to -1.80 mm) |
| **Storage** | Freezing buffer (30% sucrose, 30% glycerol), -20°C |

#### 4.2 WFA Staining (Free-Floating)
| Step | Details |
|------|---------|
| **1. Free-floating** | Sections in 24-well plate, PBS washes |
| **2. Permeabilization** | 0.3% Triton X-100, 15 min |
| **3. Blocking** | 5% normal serum, 1 h |
| **4. WFA primary** | WFA-Biotin, overnight 4°C, gentle rocking |
| **5. Secondary** | Streptavidin-AlexaFluor 488, 2 h RT |
| **6. DAPI** | 10 min, then wash |
| **7. Storage** | In PBS with azide, 4°C (short-term) |

#### 4.3 Image Analysis (ImageJ Macro)
```javascript
// WFA puncta counting macro
run("Threshold...", "method=Default white black white display=Standard 
    slice=1 pixel=Min max=255");
run("Set Measurements...", "mean standard redirect=None 
    histogram display=Histograms");
run("Convert to Mask");
run("Watershed");
run("Analyze Particles...", "size=50-Infinity 
    circularity=0.30-0.90 show=Outlines display=Overlay 
    clear");

// Measure results
selected = getResult("Area");
puncta_count = getResult("Count");
```

| Metric | Formula |
|--------|---------|
| **WFA+ area fraction** | Σ(WFA+ area) / Σ(total ARC area) × 100 |
| **Puncta density** | Σ(puncta) / Σ(NeuN+ cells) |
| **Mean puncta size** | Σ(area) / Σ(count) |

#### 4.4 Statistical Analysis
| Analysis | Test |
|----------|------|
| **Multiple groups** | One-way ANOVA + Tukey HSD |
| **Two groups** | Student's t-test (two-tailed) |
| **Repeated measures** | Two-way ANOVA (treatment × time) |
| **Correlation** | Pearson or Spearman r |
| **Effect size** | Cohen's d for t-test, η² for ANOVA |

---

## 4. Data Management & Analysis Pipeline

### 4.1 Data Structure
```
/data/
├── in_vitro/
│   ├── pnn_degradation/
│   │   ├── ec50_data.xlsx
│   │   └── time_course_images/
│   ├── bbB_permeability/
│   │   ├── mdck_results.xlsx
│   │   └── pampl_results.xlsx
│   └── insulin_signaling/
│       ├── western_blot_analysis/
│       └── plate_reads.xlsx
├── in_vivo/
│   ├── pk_studies/
│   │   ├── mouse_pk.xlsx
│   │   └── brain_plasma_ratios.xlsx
│   ├── pd_studies/
│   │   ├── wfa_quantification/
│   │   │   ├── imagej_analysis/
│   │   │   └── statistical_summary.xlsx
│   │   ├── itt_results.xlsx
│   │   └── food_intake_analysis.xlsx
│   └── raw_data/
│       ├── mouse_1/
│       │   ├── pk_0.5h_plasma.csv
│       │   ├── pk_0.5h_brain.csv
│       │   └── ...
│       └── mouse_2/
│           └── ...
└── analysis/
    ├── sAR_reports/
    ├── lead_profiles.xlsx
    └── candidate_ranking.xlsx
```

### 4.2 Analysis Scripts (Python)

#### 4.2.1 PK Analysis
```python
import pandas as pd
import numpy as np
from scipy.integrate import trapezoid

def calculate_pk_params(time, conc, dose):
    """Calculate comprehensive PK parameters"""
    
    # 1. AUC (trapezoidal)
    auc_trap = trapezoid(conc, time)
    
    # 2. Terminal phase rate constant
    # Use last 3 points, log-linear regression
    logc = np.log(conc[-3:])
    t = time[-3:]
    ke = -np.polyfit(t, logc, 1)[0]
    
    # 3. Half-life
    t12 = np.log(2) / ke
    
    # 4. Clearance
    cl = dose / auc_trap
    
    # 5. Volume of distribution
    vd = cl / ke
    
    # 6. Bioavailability (if iv data available)
    # f = (auc_po * dose_iv) / (auc_iv * dose_po)
    
    return {
        'auc_trap': auc_trap,
        'ke': ke,
        't12': t12,
        'cl': cl,
        'vd': vd
    }

# Example usage
time = [0.5, 1, 2, 4, 8, 24]
conc = [150, 120, 80, 40, 15, 3]  # ng/mL
dose = 10  # mg/kg

pk = calculate_pk_params(time, conc, dose)
print(pk)
```

#### 4.2.2 ITT Analysis
```python
def analyze_itt(baseline, glucose_data, time_points):
    """Analyze Insulin Tolerance Test data"""
    
    # Calculate % change from baseline
    pct_change = (glucose_data - baseline) / baseline * 100
    
    # Glucose AUC (trapezoidal)
    auc = trapezoid(glucose_data, time_points)
    
    # Time to normoglycemia (< 150 mg/dL)
    normo_time = None
    for i, t in enumerate(time_points):
        if glucose_data[i] <= 150:
            normo_time = t
            break
    
    # Overall stats
    auc_reduction = None  # compared to vehicle
    max_glucose = np.max(glucose_data)
    
    return {
        'pct_change': pct_change,
        'auc': auc,
        'normo_time_min': normo_time,
        'max_glucose': max_glucose
    }

# Example
baseline = 180  # mg/dL
glucose = [180, 165, 150, 140, 130, 145, 160]  # at 0, 15, 30, 45, 60, 90, 120 min
time = [0, 15, 30, 45, 60, 90, 120]

itt = analyze_itt(baseline, glucose, time)
print(itt)
```

#### 4.2.3 WFA Quantification
```python
def quantify_wfa(area_fractions, puncta_counts, neuron_counts):
    """Quantify WFA lectin staining data"""
    
    # Mean area fraction
    mean_area_frac = np.mean(area_fractions) if area_fractions else 0
    
    # Puncta per neuron (density)
    if puncta_counts and neuron_counts:
        puncta_density = np.mean([p/n for p, n in zip(puncta_counts, neuron_counts)])
    else:
        puncta_density = 0
    
    # Combined score (weighted)
    combined_score = 0.6 * mean_area_frac + 0.4 * puncta_density
    
    return {
        'mean_area_fraction': mean_area_frac,
        'puncta_density': puncta_density,
        'combined_score': combined_score
    }

# Example data from 3 mice per group
area_fracs = [0.05, 0.06, 0.04]  # 5-6% WFA+ area
punctas = [12, 15, 10]  # puncta per neuron
neurons = [100, 110, 95]  # NeuN+ cells counted

result = quantify_wfa(area_fracs, punctas, neurons)
print(result)
```

### 4.3 Statistical Reporting Standards

| Parameter | Recommended Practice |
|-----------|----------------------|
| **n per group** | Power-calculated (typically n ≥ 5 for in vivo, n ≥ 3 for in vitro) |
| **Outlier detection** | Grubbs' test (α = 0.05), report if removed |
| **Normality test** | Shapiro-Wilk, choose parametric vs non-parametric |
| **Multiple comparisons** | Correction: Tukey, Sidak, or Benjamini-Hochberg |
| **Effect size** | Report Cohen's d or η² alongside p-values |
| **Data exclusion** | Pre-specified criteria, document in methods |
| **Reproducibility** | Same batch/operator for entire experiment when possible |

### 4.4 Quality Control Checks

| Assay Type | QC Check | Acceptance Criterion |
|------------|----------|---------------------|
| **PNN degradation** | ChABC positive control | ≥ 80% degradation at EC50 |
| **BBB permeability** | Known P-gp substrate (verapamil) | ER 5-10 (high efflux) |
| **Insulin signaling** | Insulin-only positive control | ≥ 2-fold pIR increase |
| **Cytotoxicity** | Vehicle control only | ≥ 95% viability at test conc. |
| **WFA staining** | Fluorescence standardizer | Consistent across runs |
| **PK analysis** | Linear regression R² | > 0.95 for terminal phase |

---

## 5. Project-Specific In Vivo Success Criteria

### 5.1 Phase 4 Go Criteria (Proof-of-Concept)

| Metric | Minimum | Target | Measurement |
|--------|---------|--------|-------------|
| **Oral F** | > 10% | > 20% | PK study (mouse) |
| **Brain/plasma ratio** | > 0.1 | > 0.3 | PK study (mouse) |
| **Hypothalamic PNN reduction** | > 15% (WFA) | > 30% (WFA) | WFA quantification, ARC |
| **Insulin signaling rescue** | > 20% pIR increase | > 40% pIR increase | Western blot, ARC microdissection |
| **ITT AUC improvement** | > 15% | > 25% | ITT, 0.75 U/kg IP |
| **Food intake reduction** | > 10% | > 15% | Automated monitoring, 24 h |
| **Body weight change** | < 5% loss (or gain < 10%) | < 10% change | Digital scale, weekly |
| **Cytotoxicity (CC50)** | > 30 μM | > 100 μM | CellTiter-Glio, mHypoA cells |

**Go decision**: Advance if ≥ 6 of 8 metrics meet target criteria

### 5.2 Phase 5 Lead Optimization Criteria

| Metric | Lead Candidate | Optimized Candidate | Target |
|--------|---------------|---------------------|--------|
| **Oral F** | > 15% | > 30% | > 30% |
| **Brain/plasma** | > 0.3 | > 0.5 | > 0.5 |
| **PNN EC50** | < 1 μM | < 100 nM | < 100 nM |
| **ITT improvement** | > 20% | > 35% | > 35% |
| **Food intake reduction** | > 15% | > 25% | > 25% |
| **hERG IC50** | > 30 μM | > 100 μM | > 100 μM |
| **CYP3A4 inhibition** | < 50% | < 20% | < 20% |
| **Max tolerated dose** | > 30 mg/kg | > 100 mg/kg | > 100 mg/kg/kg |

**Decision**: Optimized candidate advances to IND-enabling package if all target specs met.

---

## 6. Cold-Start Protocol for Experimental Work

### 6.1 When Starting a New Session
1. **Load experimental protocols** from `/experimental/`
2. **Set up animal facility** (if doing in vivo)
3. **Prepare cell culture plates** (if doing in vitro)
4. **Initialize data tracking spreadsheets**

### 6.2 Required Materials & Equipment
| Category | Items |
|----------|-------|
| **Animals** | CD-1 or C57BL/6 mice, 8-12 weeks, chow or HFD |
| **Cell lines** | mHypoA-60/3, HEK293-IR (overexpressing) |
| **Reagents** | ChABC, WFA-Biotin, Streptavidin-AF488, Insulin |
| **Assay kits** | CellTiter-Glo 2.0, Western Blot kits |
| **Equipment** | Confocal microscope, plate reader, LC-MS/MS |
| **Software** | ImageJ, GraphPad Prism, Python (pandas, scipy) |

### 6.3 First-Day Checklist
- [ ] Animals acclimated (minimum 1 week)
- [ ] Cell culture plates seeded (24-48 h prior)
- [ ] All reagents thawed and quantified
- [ ] Data spreadsheets initialized (columns as per Section 4.1)
- [ ] Equipment calibrated (plate reader, microscope)
- [ ] Ethical approval confirmed (IACUC protocol)
- [ ] Randomization scheme prepared (for in vivo studies)

### 6.4 Common Issues & Troubleshooting

| Issue | Likely Cause | Remedy |
|-------|--------------|--------|
| **No PNN degradation in assay** | Compound insoluble, wrong pH | Check solubility, adjust buffer pH 7.4 |
| **High cytotoxicity** | Compound too concentrated, off-target | Reduce concentration, counter-screen |
| **No BBB penetration** | P-gp efflux, high TPSA | Reformulate, optimize physicochemical properties |
| **No insulin signaling rescue** | Target engagement failed | Verify compound brain levels, check dosing |
| **High variability in WFA** | Section quality, antibody batch | Standardize sectioning, validate antibody |
| **PK below detection** | Too rapid metabolism, poor absorption | Formulation optimization, prodrug strategy |

---

## 7. Inter-Assay Data Integration

### 7.1 Scoring Matrix

| Compound | PNN EC50 (μM) | BBB Pe (×10⁻⁶ cm/s) | Insulin EC50 (μM) | Cytotoxicity CC50 (μM) | SI = CC50/EC50 | Priority Score |
|----------|---------------|--------------------|-------------------|------------------------|----------------|----------------|
| **Compound A** | 5.2 | 12.5 | 8.7 | > 100 | > 11.5 | 0.82 |
| **Compound B** | 15.0 | 3.2 | 25.0 | 80 | 3.2 | 0.45 |
| **Compound C** | 2.1 | 8.7 | 3.5 | > 200 | > 57.1 | 0.91 |
| **Compound D** | > 50 | 1.5 | > 100 | 50 | < 1 | 0.12 |

**Priority Score formula**: Weighted combination of normalized metrics
- PNN potency (40%) + BBB penetration (30%) + Insulin signaling (30%)
- Higher score = better overall profile

### 7.2 Multi-Parameter Optimization (MPO)

| Metric | Weight | Normalization | Target |
|--------|--------|---------------|--------|
| **PNN EC50** | 25% | log(EC50) → 0-1 (best) | < 1 μM |
| **BBB Pe** | 25% | log(Pe) → 0-1 (best) | > 5 × 10⁻⁶ cm/s |
| **Insulin EC50** | 20% | log(EC50) → 0-1 (best) | < 10 μM |
| **Cytotoxicity SI** | 20% | log(SI) → 0-1 (best) | > 50 |
| **Oral F (predicted)** | 10% | 0-1 scale (best) | > 20% |

**MPO score**: Sum of (weight × normalized metric)
- **Excellent**: > 0.8
- **Good**: 0.6-0.8
- **Fair**: 0.4-0.6
- **Poor**: < 0.4

---

*Document Version: 1.0*
*Generated: 2026-10-11*
*Compatible with Recursive Language Model cold-start workflows*