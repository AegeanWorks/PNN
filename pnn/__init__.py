"""Phase 2-5 discovery, preclinical and translational pipeline for PNN.

Exports are resolved lazily so `python -m pnn.<module>` does not import the module it is about to
run (which would trigger runpy's double-import warning).
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS = {
    "AssayError": ".assay",
    "BioisostereError": ".bioisostere",
    "BioisostereRule": ".bioisostere",
    "Candidate": ".curation",
    "ClinicalDoseProjection": ".ind",
    "Compound": ".compounds",
    "CompoundError": ".compounds",
    "CompoundSar": ".sar",
    "CurationError": ".curation",
    "CytotoxicityResult": ".assay",
    "DOSE_PROJECTION_FIELDS": ".ind",
    "DegradationResult": ".assay",
    "DockingError": ".docking",
    "DockingScore": ".docking",
    "DynamicsError": ".dynamics",
    "GeneratedAnalog": ".bioisostere",
    "HumanDoseProjection": ".ind",
    "IndError": ".ind",
    "LeadRun": ".selection",
    "LeadSelection": ".selection",
    "MpoError": ".mpo",
    "MpoProfile": ".mpo",
    "PathwayError": ".pathway",
    "PdOutcome": ".preclinical",
    "PermeabilityResult": ".assay",
    "PipelineError": ".selection",
    "PkProfile": ".preclinical",
    "PreclinicalCandidate": ".preclinical",
    "PreclinicalError": ".preclinical",
    "QualifiedHit": ".assay",
    "REGIMEN_FIELDS": ".pathway",
    "RegimenSummary": ".pathway",
    "RunSummary": ".runner",
    "RunnerError": ".runner",
    "SAFETY_FIELDS": ".safety",
    "SIMULATION_FIELDS": ".pathway",
    "SafetyAssessment": ".safety",
    "SafetyError": ".safety",
    "SafetyPharm": ".safety",
    "SarError": ".sar",
    "ScaffoldCluster": ".sar",
    "ScreenResult": ".screen",
    "SimulationPoint": ".pathway",
    "StudioCompound": ".studio",
    "StudioError": ".studio",
    "TargetProfile": ".target",
    "ToxStudy": ".safety",
    "TrajectoryAnalysis": ".dynamics",
    "TrajectoryPoint": ".dynamics",
    "analyze_sar": ".sar",
    "analyze_trajectory": ".dynamics",
    "assess_safety": ".safety",
    "calculate_allometric_scaling": ".ind",
    "calculate_cns_mpo": ".mpo",
    "calculate_nca": ".preclinical",
    "compile_ind_dossier": ".ind",
    "curate": ".curation",
    "deduplicate": ".curation",
    "evaluate_candidates": ".preclinical",
    "execute_turnkey_pipeline": ".runner",
    "extract_scaffold": ".sar",
    "generate_bioisosteres": ".bioisostere",
    "generate_circuit_dynamics_svg": ".studio",
    "generate_trajectory_svg": ".dynamics",
    "generate_vina_config": ".target",
    "get_target": ".target",
    "list_targets": ".target",
    "load_candidate_records": ".ind",
    "load_compounds": ".compounds",
    "load_cytotoxicity_data": ".assay",
    "load_degradation_data": ".assay",
    "load_docking_scores": ".docking",
    "load_pd_data": ".preclinical",
    "load_permeability_data": ".assay",
    "load_pk_data": ".preclinical",
    "load_safety_pharmacology": ".safety",
    "load_studio_data": ".studio",
    "load_tox_studies": ".safety",
    "load_trajectory_csv": ".dynamics",
    "parse_compounds": ".compounds",
    "profile_candidate": ".mpo",
    "project_clinical_candidates": ".ind",
    "qualify_leads": ".assay",
    "render_assay_report": ".assay",
    "render_bioisostere_report": ".bioisostere",
    "render_ind_dossier": ".ind",
    "render_lead_report": ".selection",
    "render_mpo_report": ".mpo",
    "render_pathway_report": ".pathway",
    "render_preclinical_report": ".preclinical",
    "render_safety_report": ".safety",
    "render_sar_report": ".sar",
    "render_studio_html": ".studio",
    "render_trajectory_report": ".dynamics",
    "run_bioisostere_pipeline": ".bioisostere",
    "run_dynamics_pipeline": ".dynamics",
    "run_mpo_analysis": ".mpo",
    "run_pipeline": ".selection",
    "screen_compounds": ".screen",
    "select_leads": ".selection",
    "simulate_circuit_dynamics": ".pathway",
    "summarize": ".screen",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module_name, __name__), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return __all__
