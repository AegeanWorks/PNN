# Decisions

Why this project is built the way it is, oldest first. One entry per decision that changes a
contract, a dependency, a data format or the architecture; IDs are consecutive and never reused.

<!-- Entry format:
### ADR-NNN: <decision in one line>
- **Date**: YYYY-MM-DD
- **Status**: Accepted
- **Context**: <what forced the decision>
- **Decision**: <what was decided>
- **Consequences**: <what follows, including what it costs>
-->

### ADR-001: Use a deterministic CSV property screen for Phase 2 triage
- **Date**: 2026-10-22
- **Status**: Accepted
- **Context**: Phase 2 needs a reproducible way to prioritise compounds before docking, but the project has no chemistry toolkit dependency and no validated in-house BBB model.
- **Decision**: Ship a standard-library-only `pnn` package that consumes a documented CSV schema, applies the project's CNS/BBB property windows and optional ADMET limits, and emits tiered CSV/JSON output. Rules live in `pnn/rules.py`; the contract lives in `computational/screening_pipeline.md`.
- **Consequences**: Screening is reproducible and testable without third-party packages, and it cannot overstate evidence because missing values are `UNKNOWN`. The rules remain heuristics: every promoted compound still needs measured permeability and PNN-softening data, and any rule or schema change needs a new ADR and tests.

### ADR-002: Normalize docking engines before consensus ranking
- **Date**: 2026-10-22
- **Status**: Accepted
- **Context**: Phase 2 exports Vina, Glide and MM-GBSA scores whose native scales are not comparable; summing raw values would let one engine dominate by magnitude alone.
- **Decision**: `python -m pnn.rank` min-max scales each method within the filtered docking set, weights Vina 0.5, Glide 0.3 and MM-GBSA 0.2, renormalizes over available methods, reports coverage, and requires a Vina score for entry.
- **Consequences**: The shortlist is reproducible and transparent about missing evidence, but normalization depends on the input set, so comparisons require a fixed compound set. Consensus remains a prioritisation score, not a measured affinity.

### ADR-003: Curate libraries by provenance, exact identity and property-space caps
- **Date**: 2026-10-22
- **Status**: Accepted
- **Context**: Phase 2 needs to merge several compound sources without duplicate identities or property-space clumps, but the project has no chemistry toolkit for canonicalization or scaffold analysis.
- **Decision**: `python -m pnn.curate` collapses duplicate `compound_id` values (refusing conflicting SMILES), then exact non-empty SMILES strings, screens every row with the property rules, and caps each normalized six-property bin at `--per-bin` compounds before an optional global limit. Every kept row keeps its source path and bin label; JSON reports every removal count.
- **Consequences**: The library is reproducible and auditable, and duplicate chemistry is not silently multiplied. Exact-SMILES matching misses salts, tautomers and stereochemistry variants, and bin caps approximate property-space spacing rather than guaranteeing a distance; upstream canonicalization or a chemistry toolkit is required when those distinctions matter.

### ADR-004: End-to-end Phase 2 lead selection pipeline with hashed manifest and report
- **Date**: 2026-10-24
- **Status**: Accepted
- **Context**: The Phase 2 curation, property screening, and docking consensus ranking tools existed as individual CLI steps, but generating an auditable, reproducible lead selection required manual chaining, with no unified cryptographic provenance or consolidated evidence report.
- **Decision**: Ship `pnn.selection` and `python -m pnn.pipeline` to execute the full workflow in one atomic run: curate -> screen -> consensus rank. The run emits `curated_library.csv`, `ranked_compounds.csv`, `shortlist.csv`, a Markdown `lead_report.md` (highlighting evidence gaps and G-001/G-002/G-003 triage limits), and a `manifest.json` recording SHA-256 digests for all inputs, options, stage metrics and output artifacts under schema `pnn.phase2.pipeline/1`. Shared CSV formatting logic is centralized in `pnn/tables.py`.
- **Consequences**: Lead selection is strictly reproducible, tamper-evident, and auditable across all stages. Memory is bounded by streaming SHA-256 digest computation. Triage disclaimers in generated reports prevent downstream overclaiming.

### ADR-005: Target binding-pocket profiles and Phase 3 in vitro assay qualification pipeline
- **Date**: 2026-10-24
- **Status**: Accepted
- **Context**: Virtual screening required codified receptor grid invariants for key PNN targets (PTPσ, LAR, NgR1, MMP2), and progression beyond in silico triage required an executable qualification engine to ingest experimental PNN degradation, BBB permeability, and cytotoxicity counter-screen data.
- **Decision**: Introduce `pnn.target` defining receptor pocket invariants and AutoDock Vina grid configuration exports; ship benchmark libraries in `data/compounds/` and `data/docking/` with turnkey verification runner `tools/golden_run.py`; implement `pnn.assay` to ingest fluorogenic DQ-CS degradation (EC50, max %), PAMPA-BBB / MDCK-MDR1 permeability, and CellTiter-Glo cytotoxicity (CC50, Selectivity Index) to classify candidates into CONFIRMED_HIT, QUALIFIED_LEAD, CYTOTOXIC, or INACTIVE status with Markdown reporting honoring G-001.
- **Consequences**: Structural invariants ensure reproducible receptor grid definitions, and in vitro assay datasets can be qualified objectively against standardized progression thresholds without third-party dependencies.



### ADR-006: In vivo PK/PD candidate qualification, Bemis-Murcko SAR decomposition, and standalone HTML studio
- **Date**: 2026-10-24
- **Status**: Accepted
- **Context**: Advancing PNN-softening candidates towards preclinical nomination required: (1) executable non-compartmental pharmacokinetic and in vivo pharmacodynamic (PK/PD) analysis; (2) Bemis-Murcko core scaffold and SAR decomposition to guide analog optimization without external cheminformatics libraries; and (3) an interactive visualization dashboard presenting multi-target affinity, physicochemical radar, and PK concentration-time profiles in a self-contained format.
- **Decision**: Implement `pnn.preclinical` for linear trapezoidal NCA (Cmax, AUC_0_t, AUC_0_inf, lambda_z, t1/2, F, Kp,brain) and candidate qualification; implement `pnn.sar` with standard-library biconnected component Bemis-Murcko scaffold extraction and SAR activity-property clustering; implement `pnn.studio` for generating a standalone, offline-ready HTML5/SVG dashboard (`studio.html`); export all public symbols in `pnn/__init__.py`; provide template datasets in `data/preclinical/`.
- **Consequences**: Enables end-to-end preclinical decision making strictly within Python's standard library, adhering to principles G-001 (in vivo translation limits), G-002 (docking consensus relative to docked set), and G-003 (exact SMILES non-canonicalization).

### ADR-007: One-click turnkey execution pipeline and Phase 5 multiparametric optimization engine
- **Date**: 2026-10-24
- **Status**: Accepted
- **Context**: The discovery pipeline spanned several modular subpackages (`pnn.curate`, `pnn.screen`, `pnn.docking`, `pnn.assay`, `pnn.preclinical`, `pnn.sar`, `pnn.studio`). Running the full discovery workflow and viewing results required chained CLI invocations. Furthermore, selecting development candidates required systematic multiparametric optimization (MPO) combining continuous Pfizer CNS MPO scoring with Derringer-Suich composite desirability.
- **Decision**: Introduce turnkey one-click launcher scripts (`run.py`, `start.bat`, `start.sh`, and `python -m pnn.run` in `pnn/runner.py`) that run the end-to-end cascade (curate -> screen -> dock -> in vitro assays -> preclinical PK/PD -> Bemis-Murcko SAR -> studio HTML) in under 0.15s; implement `pnn.mpo` to compute continuous CNS MPO (0-6 scale), Derringer-Suich geometric composite desirability ($D \in [0, 1]$), and diagnostic chemical liability recommendations.
- **Consequences**: Provides zero-friction 1-click execution for developers, biologists, and computational chemists across Windows, Linux, and macOS platforms; candidate triage balances target potency, CNS exposure, oral bioavailability, and safety margins quantitatively.

### ADR-008: Hypothalamic dosing dynamics, human allometric IND projection, bioisosteres, and MD stability analysis
- **Date**: 2026-10-24
- **Status**: Accepted
- **Context**: Advancing PNN-softening preclinical leads toward clinical IND candidacy required four core translational capabilities: (1) modeling hypothalamic circuit dynamics and insulin receptor accessibility via ODEs to identify an intermittent pulsatile regimen that prevents maladaptive synaptic plasticity; (2) projecting human clinical starting dose (HED, MRSD with 10x safety divisor) and compiling an ICH M3(R2)/S7A/S7B regulatory dossier; (3) systematically generating bioisosteric analogs that shield metabolic soft spots and tune basicity; and (4) assessing atomistic kinetic binding stability of docking poses from molecular dynamics (MD) trajectories.
- **Decision**: Implement `pnn.pathway` (continuous Runge-Kutta ODE simulation of PK depot absorption, PNN turnover, insulin receptor accessibility, downstream pAKT signaling, and body weight progression); `pnn.ind` (FDA/ICH allometric scaling, NOAEL, Phase 1 MRSD, and IND dossier compiler); `pnn.bioisostere` (6 medicinal chemistry rules for fluorination, basicity attenuation, deuteration, and heterocycle replacement); and `pnn.dynamics` (time-series RMSD/RMSF statistical analysis, pocket residue contact persistence, and pure SVG visualization); ship benchmark trajectory datasets in `data/dynamics/`; integrate all modules into turnkey runner `pnn/runner.py` and package export `pnn/__init__.py`.
- **Consequences**: Provides complete translational and regulatory dossier support strictly within standard library Python. In silico ODE dynamics and MD trajectory metrics serve as prioritization hypotheses and remain governed by G-001/G-002, requiring experimental biophysical (SPR/ITC) and pharmacokinetic validation before clinical administration.

### ADR-009: Deterministic IND-enabling safety pharmacology and toxicology package
- **Date**: 2026-10-24
- **Status**: Accepted
- **Context**: The Phase 5 candidate package and IND dossier existed, but safety pharmacology, genotoxicity and repeat-dose NOAEL evidence were not ingested, assessed or connected to the regulatory checklist. The dossier's ICH S7A/S7B items were pre-checked without any underlying data, which risked overclaiming safety.
- **Decision**: Add `pnn.safety`, a standard-library module that ingests repeat-dose toxicology studies (`compound_id`, `species`, `noael_mg_kg`, optional LOAEL, duration, target organ, findings) and an ICH S7A/S7B panel (`herg_ic50_um`, `free_cmax_um`, CYP3A4/CYP2D6 inhibition, Ames, micronucleus, CNS/respiratory NOAEL). It selects the most sensitive species NOAEL, computes NOAEL HED and a NOAEL-based MRSD, derives the hERG margin and classifies each compound as `SAFE_TO_PROGRESS`, `CAUTION` or `FAIL_SAFETY`. Missing evidence is UNKNOWN and forces `CAUTION`. The runner executes the assessment against `data/toxicology/` templates, and the IND dossier now lists the safety items as unchecked filing requirements rather than completed findings.
- **Consequences**: Safety evidence is auditable and consistently gated before IND filing, and absent data can no longer appear as a passed checklist. Thresholds are deterministic triage criteria, not a GLP audit; study quality, pathology review, toxicokinetic coverage and regulatory acceptability remain with qualified toxicologists, and the NOAEL-based MRSD must be reconciled with the allometric projection in `pnn.ind`.

### ADR-010: Evidence integrity over convenience defaults in translational modules
- **Date**: 2026-10-24
- **Status**: Accepted
- **Context**: Verification of the translational platform found four modules inventing evidence. `pnn.ind` defaulted Kp=0.40, F=25%, dose=10 mg/kg, mouse CL=1.8, Vd=5.2 and NOAEL=5x dose when columns were absent, and ignored the nominated lead, so a dossier could claim "robust" exposure and a therapeutic window from no data. `pnn.bioisostere` produced ring-corrupted SMILES and double-counted HBA/rotatable bonds. `pnn.dynamics` defaulted protein RMSD, contacts and MM-GBSA, yielding a fabricated STABLE_POSE. `pnn.pathway` reported the absolute pAKT level as "restoration" and hardcoded intermittent dosing as optimal.
- **Decision**: (1) Inputs drive every number; absent measurements stay empty and carry `UNMEASURED`, with `evidence_status` COMPLETE/PARTIAL/INSUFFICIENT_EVIDENCE. (2) `pnn.ind` takes measured PK from `--candidates`/preclinical NCA output and measured NOAEL from `--tox`, honours `--lead`, and reconciles its NOAEL HED with `pnn.safety`. (3) Bioisostere rules carry explicit property deltas, renumber ring labels to a free digit, and skip a transform that would cut an atom out of a ring. (4) Dynamics treats protein RMSD, contacts and MM-GBSA as optional, returns INSUFFICIENT_EVIDENCE without contacts, and its CLI no longer exposes unused report/svg flags. (5) Pathway integrates with explicit forward Euler, is homeostatic around the DIO AKT baseline, reports restoration relative to that baseline, and derives the recommendation from computed metrics. (6) Package exports are lazy so `python -m pnn.<module>` no longer double-imports.
- **Consequences**: Reports cannot present absent evidence as measured, invalid analogs are skipped instead of emitted, and same-input IND/safety doses agree. The cost is more blank fields and fewer automatic analogs for motifs embedded in rings; those cases now require a chemistry toolkit or manual design rather than a silently corrupted substitution.

### ADR-011: Endogenous aggrecanase target profiles and translational studio integration
- **Date**: 2026-10-25
- **Status**: Accepted
- **Context**: (1) Target profiling in pnn.target covered PTPsigma, LAR, NgR1, MMP-2 and ChABC, but lacked the principal endogenous CNS aggrecanases (ADAMTS-4 and MMP-9) responsible for physiological CSPG interglobular cleavage. (2) studio.html displayed only Phase 2/3/4 metrics (radar, docking heatmap, PK curves), omitting the Phase 5 translational IND candidate dossier, ICH S7A/S7B safety pharmacology scorecard, and 14-day hypothalamic circuit dynamics simulation. (3) pnn.runner, run.py, and one-click shell scripts assumed data/ was in the current working directory, failing when invoked outside projects/PNN.
- **Decision**: (1) Add ADAMTS4 (PDB 2RJP, zinc catalytic domain with S1\' pocket) and MMP9 (PDB 1GKC, gelatinase catalytic domain) invariants, key catalytic residues, and Vina grid exports to pnn.target. (2) Add generate_circuit_dynamics_svg and extend render_studio_html in pnn.studio to surface IND scaling comparisons, safety pharmacology evaluations, target invariants, and 14-day pulsatile circuit dynamics in studio.html while strictly preserving G-001/G-002 triage disclaimers. (3) Anchor default data paths to PROJECT_ROOT across runner.py, run.py, start.bat, and start.sh so execution is location-independent.
- **Consequences**: Discovery workflows evaluate endogenous aggrecanase mechanisms directly, and all preclinical-to-IND translational evidence is unified into an offline-ready HTML5 dashboard. Default paths resolve cleanly from any working directory or desktop launcher.

## Licensing migration (first license)

- New versions use the Aegean Science License 1.1 (ASL-1.1).
- Rights holder: AegeanWorks; commercial use under ACL-1.0.
- Counsel review of the 1.1 texts remains pending.
