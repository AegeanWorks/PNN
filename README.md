# PNN Project: Oral BBB-Penetrant Agents for Hypothalamic PNN Softening

<!-- aegean:media:start -->
<p align="center"><img src="assets/social-preview.png" alt="PNN social preview" width="100%"></p>
<p align="center"><a href="LICENSE"><img src="https://img.shields.io/badge/license-ASL--1.1-1c7ed6" alt="License"></a> <a href="https://github.com/AegeanWorks/AIKit"><img src="https://img.shields.io/badge/AIKit-integrated-6f42c1" alt="AIKit Integrated"></a></p>
<!-- aegean:media:end -->

Oral, blood-brain barrier-penetrant agents that soften hypothalamic perineuronal nets to restore insulin signaling in obesity.

## ⚡ 30-Second Cold-Start Quick Reference

- **Runtime Prerequisites**: Python 3.10+ (standard library only, for the screening pipeline and AIKit verification) and git 2.30+ (for version control).
- **One-Click Turnkey Start**:
  - Windows: Double-click start.bat or run python run.py
  - Linux/macOS: Run ./start.sh or python run.py
  - Executes end-to-end curation, screening, docking consensus, assay qualification, preclinical PK/PD, SAR scaffold analysis, and generates studio.html in <0.2s.
- **Canonical Regression Shield Gate**:
  ```bash
  python -m unittest discover -s tests -t .
  ```
  Run it from this project root; AIKit runs the same command as part of `python aikit.py check PNN`.

### 🛠️ AIKit Commands (run from the AIKit workspace root)
| Role | Command |
|---|---|
| Status, task and gates | `python aikit.py status PNN` |
| Delivery check (clean tree + map) | `python aikit.py check PNN --strict` |
| Repair memory and map, then check | `python aikit.py fix PNN` |
| Resume an interrupted task | `python aikit.py resume PNN` |

## Mission Statement
**Ultimate Goal**: Fix insulin signaling in the hypothalamus to resolve obesity by developing oral, blood-brain barrier (BBB)-penetrant agents that soften perineuronal nets (PNNs) in the hypothalamus.

## 📦 What's Inside

```text
PNN/
├── README.md                      <- This file - project overview
├── .gitignore                     <- AIKit transient-state baseline
├── ASSETS.md
├── LICENSE
├── LICENSES/
├── REUSE.toml
├── assets/
├── COMMERCIAL.md
├── literature_review.md           <- Comprehensive literature findings
├── research_framework.md          <- Methodology & experimental design
├── run.py                         <- One-click turnkey pipeline entrypoint
├── start.bat                      <- One-click Windows launcher script
├── start.sh                       <- One-click POSIX/macOS launcher script
├── .aikit/                        <- AIKit memory: guidelines, decisions, ledger, project record
├── agent_design/                  <- Agent design specifications
│   ├── target_profiles.md         <- Target protein profiles
│   ├── bbb_penetration_rules.md   <- BBB penetration guidelines
│   └── oral_bioavailability.md    <- Oral drug design principles
├── computational/                 <- Computational approaches
│   ├── bioisostere_optimization.md  <- Medicinal chemistry rules and metabolic soft spots
│   ├── docking_protocols.md       <- Molecular docking procedures
│   ├── docking_ranking.md         <- Docking consensus contract and normalization
│   ├── dynamics_analysis.md       <- MD trajectory stability and contact persistence analysis
│   ├── lead_selection.md          <- End-to-end lead selection workflow and manifest
│   ├── library_curation.md        <- Merge, dedupe and property-space diversity rules
│   ├── mpo_optimization.md        <- Phase 5 MPO and composite desirability
│   ├── sar_analysis.md            <- Bemis-Murcko framework extraction and SAR clustering
│   ├── screening_pipeline.md      <- Phase 2 property screening contract and rules
│   └── target_profiles.md         <- Receptor binding-pocket invariants and Vina grids
├── experimental/                  <- Experimental validation
│   ├── assay_pipeline.md          <- Phase 3 in vitro assay qualification and hit criteria
│   ├── experimental_protocols.md  <- In vitro/in vivo assay protocols
│   └── preclinical_pipeline.md    <- Phase 4 in vivo PK/PD and candidate qualification
├── translational/                 <- Clinical translation and regulatory filing
│   ├── ind_filing.md              <- Human allometric scaling and IND dossier design
│   └── safety_package.md          <- IND-enabling toxicology and safety pharmacology contract
├── data/                          <- Data storage
│   ├── compounds/                 <- Compound libraries
│   │   ├── repurposing_library.csv  <- Annotated CNS repurposing library
│   │   ├── screening_template.csv  <- Input header for the screening pipeline
│   │   └── vendor_library.csv      <- Diverse drug-like screening library
│   ├── docking/                   <- Docking results
│   │   ├── multi_target_scores.csv  <- Multi-target docking scoring dataset
│   │   └── scores_template.csv    <- Input header for the consensus ranker
│   ├── dynamics/                  <- Molecular dynamics trajectory time series
│   │   └── ptpsigma_md_trajectory_template.csv  <- PTPsigma MD trajectory template
│   ├── assays/                    <- Assay results and templates
│   │   ├── cytotoxicity_template.csv  <- CellTiter-Glo counter-screen template
│   │   ├── degradation_template.csv  <- Fluorogenic DQ-CS cleavage template
│   │   └── permeability_template.csv  <- PAMPA / MDCK-MDR1 transport template
│   ├── preclinical/               <- In vivo PK/PD templates and datasets
│   │   ├── pd_metabolic_template.csv  <- In vivo PD WFA/ITT metabolic outcome template
│   │   └── pk_timecourse_template.csv  <- In vivo PK concentration-time series template
│   ├── toxicology/                <- IND-enabling toxicology and safety pharmacology data
│   │   ├── repeat_dose_tox_template.csv  <- Repeat-dose NOAEL/LOAEL study template
│   │   └── safety_pharmacology_template.csv  <- hERG/CYP/genotoxicity panel template
│   └── literature/                <- Extracted literature data (empty until populated)
├── docs/                          <- Documentation
│   ├── cold_start_guide.md        <- AI cold-start instructions
│   └── reproducibility.md         <- Reproducibility standards
├── pnn/                           <- Phase 2 pipeline (standard library only)
│   ├── __init__.py                <- Public API
│   ├── __main__.py                <- `python -m pnn` entry point
│   ├── assay.py                   <- Phase 3 in vitro assay ingestion and hit qualification
│   ├── bioisostere.py             <- Bioisosteric replacement and metabolic soft-spot engine
│   ├── cli.py                     <- Screening command contract
│   ├── compounds.py               <- Compound record and CSV ingestion
│   ├── curate.py                  <- `python -m pnn.curate` library command
│   ├── curation.py                <- Merge, dedupe and property-space caps
│   ├── docking.py                 <- Docking ingestion and consensus normalization
│   ├── dynamics.py                <- MD trajectory and contact stability analyzer
│   ├── ind.py                     <- Allometric human dose projection and IND dossier
│   ├── mpo.py                     <- Phase 5 CNS MPO and candidate profiling
│   ├── pathway.py                 <- Hypothalamic circuit dynamics and dosing simulator
│   ├── pipeline.py                <- `python -m pnn.pipeline` lead selection CLI
│   ├── preclinical.py             <- Phase 4 in vivo PK/PD NCA analysis and candidate selection
│   ├── rank.py                    <- `python -m pnn.rank` shortlist command
│   ├── reports.py                 <- Atomic CSV/JSON writers
│   ├── rules.py                   <- CNS/BBB property rules
│   ├── runner.py                  <- One-click discovery pipeline runner
│   ├── safety.py                  <- IND-enabling toxicology and safety pharmacology assessment
│   ├── sar.py                     <- Bemis-Murcko scaffold extraction and SAR activity clustering
│   ├── screen.py                  <- Tiering, scoring and ranking
│   ├── selection.py               <- End-to-end pipeline execution and lead report
│   ├── studio.py                  <- Standalone HTML/SVG interactive visualization studio
│   ├── tables.py                  <- Shared tabular CSV row contracts
│   └── target.py                  <- Receptor profiles and AutoDock Vina grid generation
├── tools/                         <- Golden runs and automation tools
│   ├── golden_run.py              <- Golden reference run and manifest verification
│   └── update_readme_tree.py      <- Helper script to sync README directory tree
└── tests/                         <- Regression shield (docs and pipeline behaviour)
    ├── __init__.py                <- Makes the test package importable
    ├── test_assay.py              <- In vitro assay qualification and hit criteria
    ├── test_bioisostere.py        <- Bioisosteric replacement and soft-spot tests
    ├── test_cli.py                <- Screening command contract
    ├── test_curation.py           <- Merge, dedupe, bins and CLI contract
    ├── test_docking.py            <- Docking ingestion and normalization
    ├── test_docs.py               <- Required docs and README tree paths
    ├── test_dynamics.py           <- MD trajectory and contact stability tests
    ├── test_golden.py             <- Golden reference pipeline verification
    ├── test_ind.py                <- Human dose projection and IND dossier tests
    ├── test_mpo.py                <- Phase 5 MPO profiling tests
    ├── test_pathway.py            <- Hypothalamic circuit dynamics tests
    ├── test_pipeline.py           <- End-to-end pipeline and manifest tests
    ├── test_preclinical.py        <- In vivo PK/PD NCA analysis and candidate selection tests
    ├── test_rank.py               <- Consensus shortlist and CLI contract
    ├── test_rules.py              <- Rule boundaries
    ├── test_run.py                <- Turnkey runner and launcher tests
    ├── test_safety.py             <- IND-enabling toxicology and safety pharmacology tests
    ├── test_sar.py                <- Bemis-Murcko scaffold extraction and SAR tests
    ├── test_screen.py             <- Tiering, ranking and CSV ingestion
    ├── test_studio.py             <- Interactive HTML/SVG studio tests
    └── test_target.py             <- Receptor profiles and Vina config generation
```

## Phase 2: Property Screening Pipeline

`pnn` applies the project's CNS/BBB design windows to a compound CSV and ranks the
candidates into tiers `A`, `B`, `C` and `Reject`. It is a triage step, not a
permeability or efficacy model: tier labels decide which compounds reach docking
and assays, never whether a compound is claimed to be BBB-penetrant.

```bash
python -m pnn data/compounds/screening_template.csv --out ranked_compounds.csv
python -m pnn compounds.csv --out ranked.csv --json summary.json --min-tier B
```

The required columns are `compound_id`, `mw`, `logp`, `tpsa`, `hbd`, `hba` and
`rotatable_bonds`; optional ADMET and efficacy fields are listed with the full
rule table in [computational/screening_pipeline.md](computational/screening_pipeline.md).
Missing optional values become `UNKNOWN` and are never treated as passes.

## Phase 2: Docking Consensus Shortlist

`python -m pnn.rank` joins per-target docking results and ranks compounds by a
normalized Vina/Glide/MM-GBSA consensus. Raw engine scores are min-max scaled
within the run and renormalized for missing methods; only rows with a Vina score
enter the shortlist.

```bash
python -m pnn.rank data/docking/scores_template.csv \
    --screen ranked_compounds.csv \
    --out shortlist.csv --json shortlist_summary.json --max-rmsd 2.0
```

`--min-tier` (default `B`) keeps compounds that passed property triage, and
`coverage` shows whether Vina, Glide and MM-GBSA were all available. The contract,
weights and limits are in [computational/docking_ranking.md](computational/docking_ranking.md).

## Phase 2: Library Curation

`python -m pnn.curate` merges compound CSVs, resolves duplicate IDs and exact
SMILES, applies the screening-tier floor, and caps each six-property bin so the
screening set is diverse in property space.

```bash
python -m pnn.curate vendor_a.csv vendor_b.csv \
    --out curated_library.csv --json curation_summary.json \
    --min-tier B --per-bin 25 --max-compounds 200000
```

Every kept row retains its `source` path and `property_bin`; the JSON summary
reports duplicate, tier and cap removals. Exact-SMILES matching is not chemical
canonicalization; see [computational/library_curation.md](computational/library_curation.md).

## Phase 2: End-to-End Lead Selection Pipeline

`python -m pnn.pipeline` executes the full Phase 2 computational triage workflow
in one reproducible run: it curates multi-source compound libraries, screens CNS/BBB
physicochemical properties, joins docking consensus scores across targets, and
writes tabular datasets, a Markdown lead report, and a hashed manifest.

```bash
python -m pnn.pipeline vendor_library.csv repurposing_library.csv \
    --docking data/docking/scores_template.csv \
    --out-dir phase2_lead_run \
    --min-tier B --per-bin 25 --max-rmsd 2.0
```

The run produces `curated_library.csv`, `ranked_compounds.csv`, `shortlist.csv`,
`lead_report.md` (highlighting evidence gaps and G-001/G-002/G-003 limits), and
`manifest.json` under schema `pnn.phase2.pipeline/1` with SHA-256 digests for all
inputs and artifacts. Details are in [computational/lead_selection.md](computational/lead_selection.md).

## Target Profiles & Vina Grid Generation

`python -m pnn.target` codifies receptor structural invariants and exports AutoDock Vina grid configurations for primary targets (`PTPsigma`, `LAR`, `NgR1`, `MMP2`):

```bash
python -m pnn.target list
python -m pnn.target show PTPsigma
python -m pnn.target config PTPsigma --receptor receptor.pdbqt --out docked.pdbqt --save vina.conf
```

Details are in [computational/target_profiles.md](computational/target_profiles.md).

## Phase 3: In Vitro Assay Qualification

`python -m pnn.assay` ingests experimental PNN degradation data (DQ-CS/WFA EC50), BBB permeability (PAMPA Pe, MDCK-MDR1 efflux), and cytotoxicity counter-screens (CellTiter-Glo CC50, Selectivity Index) to qualify in silico leads into verified in vitro hits:

```bash
python -m pnn.assay \
    --degradation data/assays/degradation_template.csv \
    --permeability data/assays/permeability_template.csv \
    --cytotoxicity data/assays/cytotoxicity_template.csv \
    --out qualified_hits.csv \
    --json assay_summary.json \
    --report assay_report.md
```

Details are in [experimental/assay_pipeline.md](experimental/assay_pipeline.md).

## IND-Enabling Safety & Toxicology Package

`python -m pnn.safety` ingests repeat-dose NOAEL/LOAEL studies and the ICH S7A/S7B safety
pharmacology panel (hERG, CYP3A4, Ames, micronucleus, CNS/respiratory NOAEL), computes the
NOAEL-based HED and starting dose, and classifies candidates as `SAFE_TO_PROGRESS`, `CAUTION`
or `FAIL_SAFETY`. Missing values stay UNKNOWN and are never treated as passes.

```bash
python -m pnn.safety --tox data/toxicology/repeat_dose_tox_template.csv \
    --pharm data/toxicology/safety_pharmacology_template.csv \
    --out output/safety_assessments.csv --json output/safety_summary.json \
    --report output/safety_package.md
```

Contracts, thresholds and evidence limits are in [translational/safety_package.md](translational/safety_package.md).

## IND Dose Projection & Dossier

`python -m pnn.ind` projects HED, MRSD and human PK from measured preclinical inputs and compiles
an IND nomination dossier. Missing measurements stay `UNMEASURED` and lower the compound's
`evidence_status`; nothing is defaulted.

```bash
python -m pnn.ind --compounds data/compounds/vendor_library.csv \
    data/compounds/repurposing_library.csv \
    --candidates output/preclinical_candidates.csv \
    --tox data/toxicology/repeat_dose_tox_template.csv \
    --out output/human_dose_projections.csv --dossier output/ind_dossier.md
```

Details are in [translational/ind_filing.md](translational/ind_filing.md).

## Golden Run Verification

To verify the end-to-end pipeline against realistic benchmark datasets in `data/`:

```bash
python tools/golden_run.py
```

## Core Scientific Hypothesis
Hypothalamic PNNs (primarily composed of chondroitin sulfate proteoglycans - CSPGs) create a physical and signaling barrier that:
1. **Impairs insulin receptor accessibility** on POMC/NPY neurons
2. **Disrupts insulin signaling cascades** (IRS-PI3K-AKT pathway)
3. **Contributes to leptin/insulin resistance** in diet-induced obesity
4. **Softening PNNs restores insulin sensitivity** and normalizes energy homeostasis

## Key Molecular Targets

### Primary Targets (PNN Degradation/Softening)
| Target | Function | Modality | BBB Penetration Strategy |
|--------|----------|----------|-------------------------|
| **ChABC** (Chondroitinase ABC) | Degrades CSPG GAG chains | Enzyme | Angiopep-2/TfR fusion |
| **MMP-2/9** | Endogenous CSPG cleavage | Small-molecule activators | Lipophilic prodrugs |
| **PTPσ** (RPTPσ) | CSPG receptor | Antagonists/antibodies | Peptide mimetics |
| **LAR** | CSPG receptor | Antagonists | Small molecules |
| **NgR1/3** | CSPG co-receptors | Antagonists | Cyclic peptides |

### Secondary Targets (Insulin Signaling Restoration)
IR (allosteric modulators), IRS-1/2 and PTP1B (phosphatase inhibition), PI3K and AKT (activators),
mTORC1 (modulators), PTEN (inhibitors).

## BBB Penetration Strategies for Oral Agents

**Physicochemical windows (CNS-modified Lipinski)**: MW 200-450 Da (optimal 300-400), logP 1.5-3.5
(optimal 2-3), TPSA <70 Å² (optimal <60), HBD ≤3, HBA ≤7, rotatable bonds ≤7, pKa 7-9 (bases) or
3-5 (acids), neutral or weakly basic at pH 7.4.

**Transporter-mediated**: LAT1 amino-acid mimetics, MCT1 carboxylic-acid prodrugs, GLUT1 glucose
conjugates, nucleoside analogs.

**Receptor-mediated transcytosis**: anti-TfR antibodies/peptides, insulin conjugates, LRP1
(Angiopep-2, ApoE mimetics).

## Recursive Language Model Compatibility

Knowledge is markdown-structured, protocols are explicit, data is version-controlled, pipelines are
deterministic, and every document is self-contained for independent cold starts.

### Cold-Start Protocol
When an AI agent starts fresh on this project:
1. Read `README.md` → `literature_review.md` → `research_framework.md`
2. Load `agent_design/target_profiles.md` for target specifications
3. Follow `computational/docking_protocols.md` for virtual screening
4. Reference `docs/cold_start_guide.md` for session initialization
5. All outputs saved to structured directories for continuity

## Literature Search Strategy

### Primary Databases
- PubMed/MEDLINE
- Embase
- Web of Science
- Scopus
- ChEMBL (compound bioactivity)
- BindingDB (binding affinities)
- PDB (structural data)

### Key Search Terms
```
("perineuronal net" OR "PNN" OR "chondroitin sulfate proteoglycan" OR "CSPG")
AND ("hypothalamus" OR "arcuate nucleus" OR "POMC" OR "NPY" OR "AgRP")
AND ("insulin signaling" OR "insulin resistance" OR "obesity" OR "energy homeostasis")
AND ("blood-brain barrier" OR "BBB" OR "CNS penetration" OR "oral bioavailability")
AND ("chondroitinase" OR "MMP" OR "PTPσ" OR "LAR" OR "NgR" OR "small molecule")
```

### Inclusion Criteria
- Peer-reviewed publications (2010-present priority)
- In vivo rodent/primate studies
- Human genetic association studies
- Structural biology data (cryo-EM, X-ray)
- Quantitative BBB penetration data (logBB, Kp,uu)

## Success Metrics

### Computational
- [ ] Virtual screening library >1M compounds
- [ ] Docking validation: ROC-AUC >0.85 on known actives
- [ ] MD simulations: >100ns stability for top 50 hits
- [ ] BBB prediction: >80% accuracy on test set

### In Vitro
- [ ] PNN degradation assay: EC50 <100 nM
- [ ] BBB permeability (PAMPA/MDCK): Pe >5 × 10⁻⁶ cm/s
- [ ] Insulin signaling rescue: >2-fold pAKT increase
- [ ] Cytotoxicity: CC50 >50 μM

### In Vivo
- [ ] Oral bioavailability: F >20%
- [ ] Brain/plasma ratio: >0.3
- [ ] Hypothalamic PNN reduction: >50% (WFA staining)
- [ ] Insulin sensitivity: >30% improvement (hyperinsulinemic-euglycemic clamp)
- [ ] Body weight reduction: >15% in DIO mice

## Timeline & Milestones

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| **Phase 1: Literature & Target Validation** | Months 1-3 | Comprehensive review, target prioritization |
| **Phase 2: Computational Discovery** | Months 3-6 | 50 lead compounds with predicted BBB penetration |
| **Phase 3: In Vitro Validation** | Months 6-12 | 10 compounds with confirmed PNN softening + BBB penetration |
| **Phase 4: In Vivo Proof-of-Concept** | Months 12-18 | 3 candidates with oral efficacy in DIO mice |
| **Phase 5: Lead Optimization** | Months 18-24 | 1 clinical candidate with IND-enabling package |

## Risk Mitigation

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| BBB penetration failure | High | Critical | Multiple parallel strategies (passive + active transport) |
| PNN softening causes plasticity side effects | Medium | High | Titratable dosing, intermittent regimens |
| Compensatory mechanisms | Medium | Medium | Combination therapy (PNN + insulin sensitizer) |
| Species translation | High | High | Human iPSC-derived hypothalamic organoids |
| Oral bioavailability | High | Critical | Prodrug strategies, formulation optimization |

## References Management
References live in `data/literature/references.bib` (BibTeX), tagged `#PNN_structure`,
`#hypothalamus_insulin`, `#BBB_penetration`, `#oral_CNS_drugs`, `#chondroitinase`,
`#CSPG_receptors` and `#obesity_models`.

---

*Last Updated: 2026-10-24*
*Version: 1.5.0*
*Compatible with AI cold starts*

<!-- aegean:brand-assets:start -->
## 🎨 Brand Assets

> **AI-Generated Brand Artwork.** The banner and logo are **AI-generated brand assets**. They are illustrative and may not literally depict every feature or the full project scope.
<!-- aegean:brand-assets:end -->

## 📜 License

Free for non-commercial research and teaching under the Aegean Science License 1.1 (ASL-1.1); see [LICENSE](LICENSE). It allows publication, reproducibility sharing, text/data mining and non-commercial ML training. Commercial or industry-sponsored use needs the Aegean Commercial License 1.0 (ACL-1.0); see [COMMERCIAL.md](COMMERCIAL.md).

<!-- aegean:screenshots:start -->
## Screenshots

<a href="https://github.com/AegeanWorks/PNN/commit/70d7c48"><img src="assets/screenshots/2026-10-28-studio.png" alt="PNN interface at 2026-10-28"></a>

*commit `70d7c48` · 2026-10-28*
<!-- aegean:screenshots:end -->
