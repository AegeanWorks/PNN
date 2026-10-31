# Codebase Map & Symbol Catalog 🗺️

> **AUTO-GENERATED & VERIFIED BY `scripts/update_knowledge_base.py`**
> This document provides an automated, up-to-date registry of all modules, classes, and exported functions across this repository.

---

## 📊 Codebase Metrics
- **Total Python Modules**: `51`
- **Total Classes**: `88`
- **Total Top-Level Functions**: `208`
- **Total Class Methods**: `239`
- **Total Source Lines**: `12695`

---

## 📦 Modules & Exported Symbols

### [`pnn/__init__.py`](pnn/__init__.py)
> *Phase 2-5 discovery, preclinical and translational pipeline for PNN.*

**Size**: `130` lines (`4319` bytes)

**Functions:**
- `def __getattr__(name)`
- `def __dir__()`

---

### [`pnn/__main__.py`](pnn/__main__.py)
> *Run the PNN screening pipeline as `python -m pnn`.*

**Size**: `6` lines (`137` bytes)

---

### [`pnn/assay.py`](pnn/assay.py)
> *Phase 3 in vitro assay ingestion, hit qualification, and reporting for PNN.*

**Size**: `525` lines (`21705` bytes)

**Classes:**
- `class AssayError` - *Raised when in vitro assay data is malformed or violates contract.*
- `class DegradationResult`
- `class PermeabilityResult`
- `class CytotoxicityResult`
- `class QualifiedHit`
  - Method `.as_row()`

**Functions:**
- `def _parse_float(val, field_name, row_id)`
- `def load_degradation_data(path)` - *Ingest PNN degradation assay CSV.*
- `def load_permeability_data(path)` - *Ingest BBB permeability assay CSV (PAMPA / MDCK-MDR1).*
- `def load_cytotoxicity_data(path)` - *Ingest CellTiter-Glo cytotoxicity counter-screen CSV.*
- `def qualify_leads(degradation, permeability, cytotoxicity, ec50_ceiling_nm, max_deg_floor_pct, pampa_pe_floor, mdck_efflux_ceiling, min_selectivity_index)` - *Qualify leads against Phase 3 in vitro degradation, permeability, and safety criteria.*
- `def render_assay_report(hits, criteria)` - *Render a comprehensive Phase 3 in vitro qualification Markdown report.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/bioisostere.py`](pnn/bioisostere.py)
> *Bioisosteric replacement and metabolic soft-spot optimization engine for PNN.*

**Size**: `533` lines (`19895` bytes)

**Classes:**
- `class BioisostereError` - *Raised when compound SMILES or bioisosteric transformation rules are invalid.*
- `class BioisostereRule`
- `class GeneratedAnalog`
  - Method `.as_row()`

**Functions:**
- `def _ring_label_occurrences(smiles)` - *Count single-digit and %nn ring labels in a SMILES string (label -> occurrences).*
- `def _ring_labels_balanced(smiles)` - *Return True when every ring label is opened and closed (even occurrence count).*
- `def _free_ring_label(parent_smiles)`
- `def _replacement_token(parent_smiles, rule)` - *Return the replacement, renumbering its ring labels away from the parent's labels.*
- `def _target_variants(motif)` - *Yield the motif with any ring labels relabelled, so a literal motif matches equivalent ring numbering.*
- `def _safe_match(smiles, rule)` - *Find the first motif occurrence that can be replaced without splitting a ring context.*
- `def generate_bioisosteres(parent_id, smiles, mw, logp, tpsa, hbd, hba, rotb, pka)` - *Apply rule-based bioisosteric transformations to generate virtual lead analogs.*
- `def render_bioisostere_report(analogs, parent_id)` - *Render a bioisostere generation and optimization Markdown report.*
- `def _parse_opt_float(val)` - *Parse optional numeric float string, returning None for empty or non-numeric cells.*
- `def run_bioisostere_pipeline(compounds_path, target_parent_id)` - *Load parent compound and generate bioisosteric analogs across all rules.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/cli.py`](pnn/cli.py)
> *Command-line entry point for the PNN Phase 2 screening pipeline.*

**Size**: `91` lines (`3406` bytes)

**Functions:**
- `def build_parser()`
- `def _filter(results, min_tier)`
- `def _report(results, kept, source, out)`
- `def _failure(what, why, fix)`
- `def main(argv)`

---

### [`pnn/compounds.py`](pnn/compounds.py)
> *Compound records and CSV ingestion for the PNN screening pipeline.*

**Size**: `234` lines (`8036` bytes)

**Classes:**
- `class CompoundError` - *Raised when input cannot be interpreted as a compound record.*
- `class Compound` - *One candidate with the properties the screening rules consume.*

**Functions:**
- `def normalize_field(name)` - *Map one CSV header to its canonical field name, or None when the column is not used.*
- `def _number(value, field, where)`
- `def _optional_number(values, field, where)`
- `def _charge(raw)`
- `def _build_compound(values, where)`
- `def parse_compounds(rows, source, allow_duplicates)` - *Turn raw CSV rows into compound records; raise CompoundError on the first bad row.*
- `def load_compounds(path, allow_duplicates)` - *Read a UTF-8 CSV file and return its compound records.*

---

### [`pnn/curate.py`](pnn/curate.py)
> *Command-line merge and diversity curation of PNN compound libraries.*

**Size**: `114` lines (`4217` bytes)

**Functions:**
- `def _positive_int(text)`
- `def _bin_width(text)`
- `def build_parser()`
- `def _failure(what, why, fix)`
- `def main(argv)`

---

### [`pnn/curation.py`](pnn/curation.py)
> *Merge, deduplicate and diversify compound libraries for PNN Phase 2.*

**Size**: `170` lines (`6038` bytes)

**Classes:**
- `class CurationError` - *Raised when libraries cannot be merged or curated safely.*
- `class Candidate` - *A compound with its source file and screening verdict.*

**Functions:**
- `def _preference(candidate)`
- `def _merge(left, right)` - *Keep the better record, filling a blank SMILES/name from the other source.*
- `def prepare_candidates(inputs)` - *Screen every input compound while keeping its source provenance.*
- `def deduplicate(candidates)` - *Collapse duplicate IDs first, then exact non-empty SMILES strings.*
- `def property_bin(compound, width)` - *Bin the six core properties on a normalized 0..1 grid.*
- `def _validate(min_tier, per_bin, bin_width, max_compounds)`
- `def curate(candidates, min_tier, per_bin, bin_width, max_compounds)` - *Filter by tier, cap each property bin, then apply the global size limit.*

---

### [`pnn/docking.py`](pnn/docking.py)
> *Docking-score ingestion and consensus ranking for PNN Phase 2.*

**Size**: `273` lines (`10305` bytes)

**Classes:**
- `class DockingError` - *Raised when docking input cannot be interpreted.*
- `class DockingScore` - *One docking result for one compound against one target.*
- `class DockingRank` - *A docking score with its normalized consensus quality.*
- `class CompoundRank` - *The best target row for one compound and how many targets it reached.*

**Functions:**
- `def _token(name)`
- `def _field(name)`
- `def _number(raw, field, where)`
- `def _optional(values, field, where)`
- `def _build(values, where)`
- `def parse_docking_scores(rows, source)` - *Turn raw CSV rows into docking records; refuse duplicates and unusable rows.*
- `def load_docking_scores(path)` - *Read a UTF-8 docking CSV.*
- `def load_screen_table(path)` - *Read the screening ranked CSV into a compound_id keyed table.*
- `def filter_scores(scores, max_rmsd, min_vina)` - *Drop rows without Vina or outside the optional RMSD/Vina limits.*
- `def normalize_scores(scores)` - *Scale each scoring method to 0..1 and combine available methods by weight.*
- `def aggregate_best(ranks)` - *Keep each compound's highest-consensus target and count its docked targets.*

---

### [`pnn/dynamics.py`](pnn/dynamics.py)
> *Molecular Dynamics (MD) trajectory and contact stability analyzer for PNN.*

**Size**: `586` lines (`25668` bytes)

**Classes:**
- `class DynamicsError` - *Raised when MD trajectory data or parameters are malformed.*
- `class TrajectoryPoint` - *A single snapshot / frame in an MD simulation trajectory.*
- `class TrajectoryAnalysis` - *Comprehensive statistical and structural analysis of an MD trajectory.*
  - Method `.as_dict()`
  - Method `.as_row()`

**Functions:**
- `def _fmt(value, spec)`
- `def _display(value, spec)`
- `def _calc_stats(values)` - *Calculate mean, standard deviation, and maximum of the supplied (non-None) values.*
- `def _parse_optional_float(raw, label, row_idx)`
- `def load_trajectory_csv(path)` - *Load MD trajectory snapshots from a CSV file; absent optional columns stay unmeasured.*
- `def analyze_trajectory(points, compound_id, target_id)` - *Analyze MD frames; metrics absent from the input stay unmeasured and block a stability verdict.*
- `def generate_trajectory_svg(analysis, points, width, height)` - *Generate a clean SVG chart displaying MD trajectory RMSD over simulation time.*
- `def render_trajectory_report(analysis)` - *Render a Markdown report documenting the MD trajectory stability analysis.*
- `def run_dynamics_pipeline(trajectory_path, compound_id, target_id, out_dir)` - *Execute the MD trajectory stability workflow, optionally writing reports and charts.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/ind.py`](pnn/ind.py)
> *Allometric human clinical dose projection and IND regulatory dossier compiler for PNN.*

**Size**: `597` lines (`25577` bytes)

**Classes:**
- `class IndError` - *Raised when clinical projection inputs violate regulatory or pharmacological contracts.*
- `class HumanDoseProjection`
  - Method `.as_row()`

**Functions:**
- `def _fmt(value, spec)`
- `def _display(value, spec)`
- `def calculate_allometric_scaling(animal_dose_mg_kg, mouse_cl_l_h_kg, mouse_vd_l_kg, noael_mg_kg, noael_km, mouse_weight_kg, human_weight_kg, km_human, safety_factor)` - *Compute FDA / ICH allometric scaling parameters from measured mouse PK and NOAEL data.*
- `def _project_one(cid, name, info, safety_factor, human_weight_kg)`
- `def project_clinical_candidates(compounds, preclinical_data, lead_id, safety_factor, human_weight_kg)` - *Calculate clinical dose projections from measured preclinical data only.*
- `def render_ind_dossier(projections, lead_id, criteria)` - *Render an IND nomination dossier; absent values are shown as UNMEASURED, never filled in.*
- `def _parse_opt_float(val)` - *Parse an optional numeric float string, returning None for empty or non-numeric cells.*
- `def _as_paths(compounds_path)`
- `def load_candidate_records(compounds_path, candidates_path, tox_path)` - *Load compound records and measured preclinical exposures for allometric scaling.*
- `def compile_ind_dossier(lead_id, compounds_path, candidates_path, tox_path, extra_pk)` - *Compile human dose projections and the IND dossier from measured inputs only.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/mpo.py`](pnn/mpo.py)
> *Phase 5 Multiparametric Optimization (MPO) and candidate profiling engine for PNN.*

**Size**: `463` lines (`19141` bytes)

**Classes:**
- `class MpoError` - *Raised when MPO input data is invalid or missing required metrics.*
- `class MpoProfile`
  - Method `.as_row()`

**Functions:**
- `def _ramp_desirability(val, low, high, increasing)` - *Calculate linear ramp desirability between low and high bounds.*
- `def calculate_cns_mpo(mw, logp, tpsa, hbd, pka, logd)` - *Calculate Pfizer CNS MPO score (0.0 to 6.0) from physicochemical properties.*
- `def profile_candidate(compound_id, name, mw, logp, tpsa, hbd, pka, ec50_nm, brain_kp, oral_f_pct, selectivity_index)` - *Calculate CNS MPO, multi-domain desirability functions, and identify optimization liabilities.*
- `def render_mpo_report(profiles)` - *Render comprehensive Phase 5 Lead Optimization & MPO Profiling markdown report.*
- `def _parse_opt_float(val)` - *Parse optional numeric float string, returning None for empty or non-numeric cells.*
- `def run_mpo_analysis(compounds_path, candidates_path, assays_path)` - *Ingest datasets, compute MPO and composite desirability, and generate profiles.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/pathway.py`](pnn/pathway.py)
> *Hypothalamic circuit, perineuronal net turnover, and dosing dynamics simulator for PNN.*

**Size**: `422` lines (`17610` bytes)

**Classes:**
- `class PathwayError` - *Raised when simulation parameters are physically invalid or malformed.*
- `class SimulationPoint`
  - Method `.as_row()`
- `class RegimenSummary`
  - Method `.as_row()`

**Functions:**
- `def _get_dosing_times(regimen, duration_h)` - *Return schedule of dose administration times in hours.*
- `def simulate_circuit_dynamics(compound_id, dose_mg_kg, regimen, duration_days, dt_h, half_life_h, cmax_per_dose, ec50_nm, mw)` - *Run an explicit forward-Euler integration of hypothalamic circuit dynamics.*
- `def render_pathway_report(summary, sample_points)` - *Render comprehensive markdown simulation report for dosing regimen.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/pipeline.py`](pnn/pipeline.py)
> *Command-line end-to-end Phase 2 lead selection for the PNN project.*

**Size**: `137` lines (`4489` bytes)

**Functions:**
- `def _positive_int(text)`
- `def _bin_width(text)`
- `def build_parser()`
- `def _failure(what, why, fix)`
- `def main(argv)`

---

### [`pnn/preclinical.py`](pnn/preclinical.py)
> *Phase 4 in vivo pharmacokinetics (PK), pharmacodynamics (PD), and preclinical candidate selection for PNN.*

**Size**: `558` lines (`24181` bytes)

**Classes:**
- `class PreclinicalError` - *Raised when in vivo preclinical data is malformed or violates contract.*
- `class ConcentrationPoint`
- `class PkProfile`
- `class PdOutcome`
- `class PreclinicalCandidate`
  - Method `.as_row()`

**Functions:**
- `def _parse_num(val, field_name, row_id)`
- `def calculate_nca(points, dose_mg_kg, compound_id, route, tissue)` - *Perform non-compartmental pharmacokinetic analysis (NCA) on concentration-time points.*
- `def load_pk_data(path)` - *Ingest concentration-time series CSV and calculate NCA PK profiles.*
- `def load_pd_data(path)` - *Ingest in vivo pharmacodynamic (PD) and metabolic outcomes CSV.*
- `def evaluate_candidates(pk_profiles, pd_outcomes, min_kp_brain, min_oral_f_pct, min_wfa_reduction_pct, min_itt_improvement_pct)` - *Triage compounds against preclinical PK exposure and PD metabolic criteria.*
- `def render_preclinical_report(candidates, criteria)` - *Render a comprehensive Phase 4 in vivo preclinical qualification report.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/rank.py`](pnn/rank.py)
> *Command-line consensus ranking of docking results for PNN Phase 2.*

**Size**: `245` lines (`8401` bytes)

**Classes:**
- `class ShortlistRow` - *One ranked compound ready for CSV output.*
  - Method `.as_row()`

**Functions:**
- `def _property_score(info)`
- `def _selected(best, screen, min_tier)` - *Filter aggregated compounds by screening evidence and sort them by consensus.*
- `def build_shortlist(scores, screen, min_tier, max_rmsd, min_vina)` - *Return shortlist rows and a summary of what each filter removed.*
- `def build_parser()`
- `def _failure(what, why, fix)`
- `def main(argv)`

---

### [`pnn/reports.py`](pnn/reports.py)
> *Atomic CSV, JSON and text report writers shared by the PNN command line tools.*

**Size**: `68` lines (`2038` bytes)

**Functions:**
- `def _open_temp(path)`
- `def _finish(path, handle)`
- `def write_csv(path, fieldnames, rows)` - *Write `rows` as CSV to a temporary file, then replace `path` atomically.*
- `def write_json(path, report)` - *Write `report` as JSON to a temporary file, then replace `path` atomically.*
- `def write_text(path, text)` - *Write `text` to a temporary file, then replace `path` atomically.*

---

### [`pnn/rules.py`](pnn/rules.py)
> *Property rules that turn the PNN BBB/CNS design criteria into verdicts.*

**Size**: `219` lines (`7784` bytes)

**Classes:**
- `class CriterionResult` - *One rule's verdict for one compound.*

**Functions:**
- `def _range(name, value, low, high, optimal_low, optimal_high, unit)`
- `def _ceiling(name, value, limit, optimal, unit)`
- `def _physicochemical(compound)`
- `def _pka(compound)`
- `def _charge(compound)`
- `def _logbb(compound)`
- `def _brain_plasma(compound)`
- `def _oral_bioavailability(compound)`
- `def _cyp3a4(compound)`
- `def _herg(compound)`
- `def _pnn_ec50(compound)`
- `def evaluate(compound)` - *Return every rule verdict for `compound` in a stable order.*
- `def weight(criterion_name)` - *Core design properties outweigh optional measured/ADMET properties.*

---

### [`pnn/runner.py`](pnn/runner.py)
> *One-click turnkey discovery pipeline runner for PNN.*

**Size**: `488` lines (`19183` bytes)

**Classes:**
- `class RunnerError` - *Raised when one-click pipeline runner cannot execute or find input data.*
- `class RunSummary`

**Functions:**
- `def execute_turnkey_pipeline(vendor_csv, repurposing_csv, docking_csv, degradation_csv, permeability_csv, cytotoxicity_csv, pk_csv, pd_csv, out_dir, min_tier, emit_root, tox_csv, safety_pharm_csv)` - *Execute complete end-to-end PNN discovery workflow from libraries to interactive studio.*
- `def print_run_banner(summary)` - *Print terminal summary banner for turnkey runner.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/safety.py`](pnn/safety.py)
> *IND-enabling safety pharmacology and repeat-dose toxicology assessment for PNN.*

**Size**: `644` lines (`26843` bytes)

**Classes:**
- `class SafetyError` - *Raised when safety pharmacology or toxicology data is malformed or violates contract.*
- `class ToxStudy`
- `class SafetyPharm`
- `class SafetyAssessment`
  - Method `.as_row()`

**Functions:**
- `def _normalize(name)`
- `def _species_key(species)`
- `def _parse_required_float(value, field_name, row_id)`
- `def _parse_optional_float(value, field_name, row_id)`
- `def _parse_optional_bool(value, field_name, row_id)`
- `def _check_percent(value, field_name, row_id)`
- `def load_tox_studies(path)` - *Ingest repeat-dose toxicology studies with NOAEL/LOAEL and target-organ findings.*
- `def load_safety_pharmacology(path)` - *Ingest ICH S7A/S7B safety pharmacology and genotoxicity panel results.*
- `def noael_hed_mg_kg(noael_mg_kg, species, human_km)` - *Convert a species NOAEL to a human equivalent dose using its body-surface-area factor.*
- `def most_sensitive_study(studies)` - *Return the study with the lowest NOAEL human equivalent dose.*
- `def _genotoxicity(safety_pharm)`
- `def assess_safety(tox_studies, safety_pharm, safety_factor, herg_margin_min, herg_margin_fail, cyp_inhibition_max_pct, cyp_inhibition_caution_pct, human_weight_kg)` - *Assess each compound against IND-enabling toxicology and safety pharmacology criteria.*
- `def render_safety_report(assessments, criteria)` - *Render a Markdown IND-enabling safety package summary.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/sar.py`](pnn/sar.py)
> *Structure-Activity Relationship (SAR) and Bemis-Murcko scaffold decomposition for PNN.*

**Size**: `772` lines (`29999` bytes)

**Classes:**
- `class SarError` - *Raised when SAR data is malformed or violates scientific contract.*
- `class CompoundSar`
  - Method `.as_row()`
- `class ScaffoldCluster`
  - Method `.as_row()`

**Functions:**
- `def parse_smiles_graph(smiles)` - *Parse SMILES into an atom list and adjacency graph.*
- `def find_biconnected_blocks(n, adj)` - *Hopcroft-Tarjan algorithm to find all biconnected components (cycle systems).*
- `def serialize_subgraph_smiles(atoms, adj, nodes)` - *Serialize an induced subgraph of atoms into a SMILES string.*
- `def _identify_scaffold_name(atoms, cycle_blocks, linkers)` - *Identify recognizable chemotype descriptor for Bemis-Murcko framework.*
- `def extract_scaffold(smiles)` - *Decompose a SMILES into Bemis-Murcko framework SMILES, chemotype name, and R-group substituents.*
- `def load_compounds(path)` - *Load compound property records from CSV.*
- `def load_docking_scores(path)` - *Load multi-target or single-target docking scores (consensus or best affinity).*
- `def load_assay_potencies(path)` - *Load in vitro assay EC50 potencies.*
- `def _parse_opt_float(val)` - *Parse optional numeric float string, returning None for empty or non-numeric cells.*
- `def load_preclinical_exposures(path)` - *Load in vivo preclinical exposures (brain_plasma_ratio, oral_bioavailability_pct).*
- `def analyze_sar(compound_records, docking_scores, assay_potencies, preclinical_exposures)` - *Decompose compounds into Bemis-Murcko scaffolds and perform SAR activity clustering.*
- `def render_sar_report(compounds, clusters)` - *Render a comprehensive Structure-Activity Relationship (SAR) markdown report.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/screen.py`](pnn/screen.py)
> *Tiering and ranking for the PNN Phase 2 property screen.*

**Size**: `79` lines (`3049` bytes)

**Classes:**
- `class ScreenResult` - *A compound, its verdicts, and the resulting triage tier.*
  - Method `.failed_criteria()`
  - Method `.unknown_criteria()`

**Functions:**
- `def _weighted_score(criteria)`
- `def _tier(criteria, score)`
- `def evaluate_compound(compound)` - *Evaluate every rule for one compound and assign its tier.*
- `def screen_compounds(compounds)` - *Evaluate and rank compounds by tier, then score (highest first), then ID.*
- `def summarize(results)` - *Count tiers and report the scored average, including rejected compounds.*

---

### [`pnn/selection.py`](pnn/selection.py)
> *End-to-end Phase 2 lead selection for the PNN project.*

**Size**: `336` lines (`12882` bytes)

**Classes:**
- `class PipelineError` - *Raised when a Phase 2 run cannot produce a defensible lead set.*
- `class LeadSelection` - *The selected compounds and the audit counts for one Phase 2 run.*
- `class LeadRun` - *A completed run: selected compounds, written artifacts and the manifest.*

**Functions:**
- `def _digest(path)` - *Return the SHA-256 hex digest of a file without loading it whole.*
- `def _load_inputs(paths)`
- `def _screen_table(results)`
- `def select_leads(inputs, docking, min_tier, per_bin, bin_width, max_compounds, max_rmsd, min_vina)` - *Curate, screen and rank one lead set from compound and docking CSVs.*
- `def _options(min_tier, per_bin, bin_width, max_compounds, max_rmsd, min_vina)`
- `def _manifest(selection, outputs, options)`
- `def render_lead_report(selection, options, top_n)` - *Render a comprehensive Markdown lead report honoring G-001, G-002, and G-003.*
- `def run_pipeline(inputs, docking, out_dir, min_tier, per_bin, bin_width, max_compounds, max_rmsd, min_vina, top_n_report)` - *Run the three Phase 2 stages and write their artifacts, report and manifest.*

---

### [`pnn/studio.py`](pnn/studio.py)
> *Interactive Lead & Assay Visualization Studio and Standalone HTML Dashboard for PNN.*

**Size**: `1117` lines (`44761` bytes)

**Classes:**
- `class StudioError` - *Raised when studio inputs are malformed or missing required data.*
- `class StudioCompound`

**Functions:**
- `def _get_val(obj, key, default)` - *Extract a field from either a dict/mapping or an object attribute with fallback.*
- `def _get_field(obj, keys, default)` - *Extract first non-empty field matching one of several alias keys.*
- `def _normalize_cns_axis(val, optimal_center, half_width)`
- `def generate_radar_svg(compounds, width, height)` - *Generate SVG radar chart for CNS drug-likeness physicochemical properties.*
- `def generate_docking_heatmap_svg(compounds, targets, cell_w, cell_h)` - *Generate multi-target docking affinity heatmap in SVG.*
- `def generate_pk_curve_svg(compounds, width, height)` - *Generate in vivo pharmacokinetic concentration-time profiles SVG.*
- `def generate_circuit_dynamics_svg(simulation_points, width, height)` - *Generate SVG time-course chart for hypothalamic circuit dynamics (PNN turnover vs pAKT restoration).*
- `def render_studio_html(compounds, radar_svg, heatmap_svg, pk_svg, circuit_svg, dossier_data, safety_data, targets_data)` - *Render full standalone HTML5 studio dashboard.*
- `def _parse_opt_float(val)` - *Parse optional numeric float string, returning None for empty or non-numeric cells.*
- `def load_studio_data(compounds_path, docking_path, assays_path, pk_path, preclinical_candidates_path)` - *Assemble all project datasets into StudioCompound instances.*
- `def build_parser()`
- `def main(argv)`

---

### [`pnn/tables.py`](pnn/tables.py)
> *CSV row contracts shared by the PNN Phase 2 command line tools.*

**Size**: `108` lines (`3053` bytes)

**Functions:**
- `def _optional(value)`
- `def curated_row(candidate, width)` - *Return one curated-library row, including source and property bin.*
- `def screen_row(result)` - *Return one ranked-screening row for a screened compound.*

---

### [`pnn/target.py`](pnn/target.py)
> *Target receptor profiles, binding-pocket invariants and Vina grid generation for PNN.*

**Size**: `313` lines (`12212` bytes)

**Classes:**
- `class TargetProfile` - *Molecular target profile defining receptor invariants and docking grid coordinates.*

**Functions:**
- `def get_target(target_id)` - *Retrieve a target profile by case-insensitive ID or alias.*
- `def list_targets()` - *Return all registered PNN target profiles.*
- `def generate_vina_config(profile, receptor_pdbqt, out_pdbqt, exhaustiveness, num_modes, energy_range)` - *Format AutoDock Vina grid parameter configuration text.*
- `def build_parser()`
- `def main(argv)`

---

### [`run.py`](run.py)
> *Turnkey one-click entrypoint for the PNN perineuronal net drug discovery project.*

**Size**: `19` lines (`476` bytes)

---

### [`tests/__init__.py`](tests/__init__.py)
**Size**: `0` lines (`0` bytes)

---

### [`tests/test_assay.py`](tests/test_assay.py)
> *Unit and CLI tests for Phase 3 in vitro assay ingestion and hit qualification.*

**Size**: `185` lines (`6942` bytes)

**Classes:**
- `class AssayTest`
  - Method `.setUp()`
  - Method `.run_cli(*argv)`
  - Method `.test_load_degradation_data()`
  - Method `.test_load_degradation_invalid_numeric_raises()`
  - Method `.test_load_permeability_and_cytotoxicity()`
  - Method `.test_qualify_leads_classifications()`
  - Method `.test_qualify_leads_untested_safety_yields_qualified_lead()`
  - Method `.test_render_assay_report_contains_guideline_and_tables()`
  - Method `.test_cli_execution_success()`
  - Method `.test_cli_failure_triad_on_missing_file()`

---

### [`tests/test_bioisostere.py`](tests/test_bioisostere.py)
> *Unit tests for bioisosteric replacement and metabolic soft-spot optimizer.*

**Size**: `150` lines (`5365` bytes)

**Classes:**
- `class TestBioisostere`
  - Method `.test_generate_bioisosteres()`
  - Method `.test_ring_label_parser_ignores_isotope_brackets()`
  - Method `.test_replacement_does_not_split_ring_context()`
  - Method `.test_property_deltas_match_rule_table()`
  - Method `.test_cyclopropyl_delta_is_not_a_rotatable_increase()`
  - Method `.test_run_bioisostere_pipeline_missing_file_raises()`
  - Method `.test_run_bioisostere_pipeline_vendor_lib()`
  - Method `.test_render_bioisostere_report()`
  - Method `.test_bioisostere_cli()`

---

### [`tests/test_cli.py`](tests/test_cli.py)
> *Integration tests for the `python -m pnn` CLI contract.*

**Size**: `84` lines (`3719` bytes)

**Classes:**
- `class CliTest`
  - Method `.setUp()`
  - Method `.run_main(*argv)`
  - Method `.read_rows(path)`
  - Method `.test_writes_ranked_csv_and_json_summary()`
  - Method `.test_min_tier_filters_lower_tiers_and_rejects()`
  - Method `.test_invalid_input_returns_one_with_a_remediation_triad()`
  - Method `.test_missing_input_returns_one()`
  - Method `.test_help_exits_zero_and_unknown_flag_exits_two()`

---

### [`tests/test_curation.py`](tests/test_curation.py)
> *Unit and CLI tests for compound-library curation.*

**Size**: `181` lines (`6910` bytes)

**Classes:**
- `class DeduplicationTest`
  - Method `.test_duplicate_id_keeps_the_better_candidate()`
  - Method `.test_better_record_keeps_smiles_from_the_other_source()`
  - Method `.test_same_id_with_different_smiles_is_refused()`
  - Method `.test_same_smiles_with_different_ids_is_collapsed()`
  - Method `.test_blank_smiles_do_not_collapse_distinct_ids()`
  - Method `.test_prepare_candidates_keeps_sources()`
- `class CurationTest`
  - Method `.test_tier_filter_removes_rejected_candidates()`
  - Method `.test_per_bin_cap_keeps_one_best_per_property_cell()`
  - Method `.test_max_compounds_truncates_and_reports()`
  - Method `.test_property_bin_is_deterministic()`
  - Method `.test_invalid_bin_width_is_refused()`
- `class CurateCliTest`
  - Method `.setUp()`
  - Method `.run_main(*argv)`
  - Method `.read_rows(path)`
  - Method `.test_merges_deduplicates_and_filters()`
  - Method `.test_conflicting_id_returns_one()`
  - Method `.test_no_eligible_compounds_returns_one()`

**Functions:**
- `def compound(compound_id, smiles, **overrides)`
- `def candidate(compound_id, source, smiles, **overrides)`

---

### [`tests/test_docking.py`](tests/test_docking.py)
> *Unit tests for docking ingestion and consensus normalization.*

**Size**: `113` lines (`4470` bytes)

**Classes:**
- `class DockingIngestionTest`
  - Method `.test_aliases_units_and_optional_columns()`
  - Method `.test_row_needs_a_score()`
  - Method `.test_row_needs_target()`
  - Method `.test_duplicate_compound_target_is_refused()`
  - Method `.test_non_finite_score_is_refused()`
  - Method `.test_negative_rmsd_is_refused()`
- `class ConsensusTest`
  - Method `.test_vina_only_is_normalized_to_zero_and_one()`
  - Method `.test_multiple_methods_are_weighted_and_coverage_is_reported()`
  - Method `.test_best_target_is_kept_and_targets_are_counted()`
  - Method `.test_filters_drop_missing_vina_and_limit_breakers()`

---

### [`tests/test_docs.py`](tests/test_docs.py)
> *Regression shield for the PNN research documents (PRINCIPLE-005, PRINCIPLE-014).*

**Size**: `50` lines (`1880` bytes)

**Classes:**
- `class ProjectDocsTest`
  - Method `.test_required_documents_exist()`
  - Method `.test_readme_tree_lists_only_existing_paths()`

---

### [`tests/test_dynamics.py`](tests/test_dynamics.py)
> *Unit tests for MD trajectory and contact stability analyzer.*

**Size**: `145` lines (`6268` bytes)

**Classes:**
- `class TestDynamics`
  - Method `.test_load_trajectory_csv()`
  - Method `.test_load_trajectory_missing_raises()`
  - Method `.test_analyze_trajectory_stable()`
  - Method `.test_analyze_trajectory_unstable_drift()`
  - Method `.test_generate_trajectory_svg()`
  - Method `.test_render_trajectory_report()`
  - Method `.test_run_dynamics_pipeline_and_cli()`
  - Method `.test_missing_optional_columns_produce_insufficient_evidence()`
  - Method `.test_partial_evidence_keeps_verdict_but_flags_gaps()`

---

### [`tests/test_golden.py`](tests/test_golden.py)
> *Integration test for the PNN Phase 2 golden reference run.*

**Size**: `30` lines (`1024` bytes)

**Classes:**
- `class GoldenRunTest`
  - Method `.test_golden_run_execution()`

---

### [`tests/test_ind.py`](tests/test_ind.py)
> *Unit tests for human allometric dose projection and IND regulatory dossier generator.*

**Size**: `216` lines (`8715` bytes)

**Classes:**
- `class TestIndProjection`
  - Method `.test_calculate_allometric_scaling()`
  - Method `.test_allometric_scaling_rejects_missing_measurements()`
  - Method `.test_project_clinical_candidates_with_complete_evidence()`
  - Method `.test_missing_evidence_is_never_invented()`
  - Method `.test_lead_is_nominated_first()`
  - Method `.test_dossier_marks_unmeasured_and_omits_claims()`
  - Method `.test_load_candidate_records_ignores_pk_timecourse_without_measurements()`
  - Method `.test_tox_path_supplies_most_sensitive_noael()`
  - Method `.test_compile_ind_dossier_uses_measured_candidates()`
- `class TestIndCli`
  - Method `.test_ind_cli()`

---

### [`tests/test_mpo.py`](tests/test_mpo.py)
**Size**: `120` lines (`4281` bytes)

**Classes:**
- `class TestMpoCalculations`
  - Method `.test_cns_mpo_optimal_properties()`
  - Method `.test_cns_mpo_penalized_properties()`
  - Method `.test_profile_development_candidate()`
  - Method `.test_profile_deprioritized_liability_detection()`
- `class TestMpoAnalysisAndCLI`
  - Method `.test_missing_compounds_raises()`
  - Method `.test_run_mpo_analysis_valid()`
  - Method `.test_render_mpo_report()`
  - Method `.test_cli_execution_success()`

---

### [`tests/test_pathway.py`](tests/test_pathway.py)
> *Unit tests for hypothalamic circuit dynamics and dosing simulator.*

**Size**: `141` lines (`5217` bytes)

**Classes:**
- `class TestPathwaySimulator`
  - Method `.test_simulate_circuit_dynamics_pulsatile()`
  - Method `.test_simulate_circuit_dynamics_continuous()`
  - Method `.test_regimen_variations()`
  - Method `.test_invalid_parameters_raise()`
  - Method `.test_untreated_baseline_has_no_restoration_or_weight_loss()`
  - Method `.test_continuous_regimen_is_flagged_for_plasticity_risk()`
  - Method `.test_render_pathway_report()`
  - Method `.test_pathway_cli()`

---

### [`tests/test_pipeline.py`](tests/test_pipeline.py)
> *Unit and CLI integration tests for the end-to-end Phase 2 lead selection pipeline.*

**Size**: `206` lines (`8185` bytes)

**Classes:**
- `class PipelineTest`
  - Method `.setUp()`
  - Method `.read_rows(path)`
  - Method `.run_cli(*argv)`
  - Method `.test_select_leads_end_to_end()`
  - Method `.test_select_leads_raises_on_empty_input()`
  - Method `.test_select_leads_raises_when_curation_empties_set()`
  - Method `.test_select_leads_raises_when_no_shortlist_survives()`
  - Method `.test_render_lead_report_guidelines_and_evidence()`
  - Method `.test_run_pipeline_produces_all_artifacts_and_hashed_manifest()`
  - Method `.test_cli_execution_success()`
  - Method `.test_cli_execution_failure_triad()`
  - Method `.test_cli_help()`

---

### [`tests/test_preclinical.py`](tests/test_preclinical.py)
**Size**: `224` lines (`8589` bytes)

**Classes:**
- `class TestPreclinicalNCA`
  - Method `.test_nca_empty_points_raises()`
  - Method `.test_nca_trapezoidal_auc_exact()`
  - Method `.test_nca_terminal_elimination_and_clearance()`
- `class TestPreclinicalDataIngestion`
  - Method `.test_load_pk_data_missing_file_raises()`
  - Method `.test_load_pk_data_missing_columns_raises()`
  - Method `.test_load_pk_data_invalid_numeric_raises()`
  - Method `.test_load_pd_data_valid()`
- `class TestPreclinicalCandidateEvaluation`
  - Method `.setUp()`
  - Method `.test_evaluate_clinical_candidate()`
  - Method `.test_evaluate_insufficient_brain_penetration()`
  - Method `.test_render_preclinical_report()`
- `class TestPreclinicalCLI`
  - Method `.test_cli_execution_success()`

---

### [`tests/test_rank.py`](tests/test_rank.py)
> *Unit and CLI tests for the docking consensus shortlist.*

**Size**: `158` lines (`6522` bytes)

**Classes:**
- `class ShortlistTest`
  - Method `.test_tier_filter_and_consensus_order()`
  - Method `.test_rmsd_filter_removes_a_target_row()`
  - Method `.test_vina_filter_removes_weak_rows()`
  - Method `.test_compound_missing_from_screen_is_skipped()`
  - Method `.test_docking_only_mode_has_no_tier()`
- `class RankCliTest`
  - Method `.setUp()`
  - Method `.run_main(*argv)`
  - Method `.read_rows(path)`
  - Method `.test_cli_writes_shortlist_and_json()`
  - Method `.test_cli_returns_one_when_filters_remove_everything()`
  - Method `.test_cli_returns_one_for_invalid_docking()`
  - Method `.test_screen_table_loader_reads_tier_and_score()`
  - Method `.test_module_entry_point_runs_main()`

**Functions:**
- `def screen_table()`
- `def scores()`

---

### [`tests/test_rules.py`](tests/test_rules.py)
> *Unit tests for the property rules (PRINCIPLE-014: failing-first boundaries).*

**Size**: `89` lines (`3378` bytes)

**Classes:**
- `class PropertyRuleTest`
  - Method `.test_ideal_compound_passes_every_rule()`
  - Method `.test_molecular_weight_boundaries()`
  - Method `.test_tpsa_and_rotor_ceilings()`
  - Method `.test_missing_optional_properties_are_unknown()`
  - Method `.test_basic_pka_outside_window_fails()`
  - Method `.test_unrecognized_charge_is_unknown()`
  - Method `.test_oral_bioavailability_accepts_fraction_and_percent()`
  - Method `.test_hard_liabilities_fail()`

**Functions:**
- `def ideal(**overrides)`
- `def by_name(compound)`

---

### [`tests/test_run.py`](tests/test_run.py)
**Size**: `116` lines (`5398` bytes)

**Classes:**
- `class TestTurnkeyRunner`
  - Method `.test_missing_vendor_library_raises()`
  - Method `.test_execute_turnkey_pipeline_success()`
  - Method `.test_runner_cli_execution_no_browser()`
  - Method `.test_runner_cli_execution_auto_opens_browser(mock_open)`

---

### [`tests/test_safety.py`](tests/test_safety.py)
> *Unit tests for the IND-enabling safety pharmacology and toxicology package.*

**Size**: `301` lines (`13158` bytes)

**Classes:**
- `class TestSafetyIngestion`
  - Method `.test_missing_tox_file_raises()`
  - Method `.test_missing_tox_columns_raises()`
  - Method `.test_unknown_species_raises()`
  - Method `.test_non_positive_noael_raises()`
  - Method `.test_loael_below_noael_raises()`
  - Method `.test_invalid_duration_raises()`
  - Method `.test_load_valid_tox_rows()`
  - Method `.test_pharm_duplicate_compound_raises()`
  - Method `.test_pharm_invalid_boolean_raises()`
  - Method `.test_pharm_percent_out_of_range_raises()`
  - Method `.test_pharm_missing_compound_id_raises()`
- `class TestSafetyAssessment`
  - Method `._safe_pharm(compound_id)`
  - Method `._tox(compound_id)`
  - Method `.test_safe_compound_projection_and_status()`
  - Method `.test_most_sensitive_species_is_selected_by_hed_not_raw_noael()`
  - Method `.test_dog_wins_when_its_hed_is_lower()`
  - Method `.test_most_sensitive_study_helper_ranks_on_hed()`
  - Method `.test_missing_tox_study_is_caution()`
  - Method `.test_missing_safety_pharmacology_is_caution()`
  - Method `.test_positive_genotoxicity_fails()`
  - Method `.test_narrow_herg_margin_fails()`
  - Method `.test_cyp3a4_ceiling_fails()`
  - Method `.test_marginal_herg_and_cyp_are_caution()`
  - Method `.test_invalid_thresholds_raise()`
  - Method `.test_report_includes_matrix_and_limits()`
- `class TestSafetyCli`
  - Method `.test_cli_requires_an_input()`
  - Method `.test_cli_missing_file_returns_error()`
  - Method `.test_cli_writes_csv_json_and_report()`

---

### [`tests/test_sar.py`](tests/test_sar.py)
**Size**: `138` lines (`5460` bytes)

**Classes:**
- `class TestSarSmilesScaffolding`
  - Method `.test_parse_empty_smiles_raises()`
  - Method `.test_parse_unmatched_parens_raises()`
  - Method `.test_acyclic_molecule_scaffold()`
  - Method `.test_monocyclic_benzene_scaffold()`
  - Method `.test_bemis_murcko_extraction_multiring()`
- `class TestSarAnalysisAndClustering`
  - Method `.test_load_compounds_missing_file_raises()`
  - Method `.test_analyze_sar_scaffold_clustering()`
  - Method `.test_render_sar_report()`
- `class TestSarCLI`
  - Method `.test_cli_execution_success()`

---

### [`tests/test_screen.py`](tests/test_screen.py)
> *Unit tests for tiering, ranking and CSV ingestion.*

**Size**: `167` lines (`6535` bytes)

**Classes:**
- `class TierTest`
  - Method `.test_complete_ideal_compound_is_tier_a()`
  - Method `.test_property_only_compound_cannot_be_tier_a()`
  - Method `.test_acceptable_compound_is_tier_b()`
  - Method `.test_two_core_failures_are_tier_c()`
  - Method `.test_three_core_failures_are_rejected()`
  - Method `.test_hard_liability_rejects_even_an_optimal_profile()`
- `class RankingTest`
  - Method `.test_ranking_orders_by_tier_then_score()`
- `class CsvIngestionTest`
  - Method `.setUp()`
  - Method `.write_csv(text)`
  - Method `.test_aliases_and_blank_optionals_are_supported()`
  - Method `.test_duplicate_ids_are_refused()`
  - Method `.test_bad_number_names_the_row()`
  - Method `.test_missing_field_is_refused()`
  - Method `.test_non_finite_number_is_refused()`
  - Method `.test_negative_property_is_refused()`
  - Method `.test_fractional_count_is_refused()`
  - Method `.test_missing_file_is_refused()`
  - Method `.test_header_only_file_is_refused()`
  - Method `.test_non_utf8_file_is_refused()`

**Functions:**
- `def ideal(compound_id, **overrides)`
- `def acceptable(compound_id, **overrides)`

---

### [`tests/test_studio.py`](tests/test_studio.py)
**Size**: `242` lines (`9955` bytes)

**Classes:**
- `class TestStudioRendering`
  - Method `.setUp()`
  - Method `.test_generate_radar_svg()`
  - Method `.test_generate_docking_heatmap_svg()`
  - Method `.test_generate_pk_curve_svg()`
  - Method `.test_render_studio_html()`
  - Method `.test_generate_circuit_dynamics_svg()`
  - Method `.test_render_studio_html_with_translational_panels()`
- `class TestStudioDataLoader`
  - Method `.test_missing_compounds_file_raises()`
  - Method `.test_load_studio_data_valid()`
- `class TestStudioCLI`
  - Method `.test_cli_execution_success()`
  - Method `.test_cli_execution_with_translational_flags()`
  - Method `.test_cli_missing_translational_files_error()`
- `class TestFieldExtraction`
  - Method `.test_get_val_with_mapping()`
  - Method `.test_get_val_with_object()`
  - Method `.test_get_field_aliases()`

---

### [`tests/test_target.py`](tests/test_target.py)
> *Unit and CLI tests for target receptor profiles and Vina grid generation.*

**Size**: `129` lines (`4857` bytes)

**Classes:**
- `class TargetProfileTest`
  - Method `.setUp()`
  - Method `.run_cli(*argv)`
  - Method `.test_list_targets_contains_primary_receptors()`
  - Method `.test_get_target_case_insensitive()`
  - Method `.test_get_target_unknown_raises_key_error()`
  - Method `.test_generate_vina_config_format()`
  - Method `.test_cli_list()`
  - Method `.test_cli_show_success_and_failure()`
  - Method `.test_cli_config_export_to_file()`
  - Method `.test_aggrecanase_targets_present()`
  - Method `.test_adamts4_target_profile()`
  - Method `.test_mmp9_target_profile()`

---

### [`tools/golden_run.py`](tools/golden_run.py)
> *Execute the PNN Phase 2 golden reference run and verify artifact integrity.*

**Size**: `131` lines (`5240` bytes)

**Functions:**
- `def _check_digest(path)`
- `def execute_golden_run(out_dir)` - *Run Phase 2 pipeline & Phase 3 assay checks against reference benchmarks.*
- `def main(argv)`

---

### [`tools/update_readme_tree.py`](tools/update_readme_tree.py)
> *Update README directory tree with translational modules and tests.*

**Size**: `132` lines (`9356` bytes)

**Functions:**
- `def update_tree()`

---
