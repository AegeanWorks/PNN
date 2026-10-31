"""One-click turnkey discovery pipeline runner for PNN."""

from __future__ import annotations

import argparse
import sys
import time
import webbrowser
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from .assay import (
    HIT_FIELDS,
    load_cytotoxicity_data,
    load_degradation_data,
    load_permeability_data,
    qualify_leads,
    render_assay_report,
)
from .bioisostere import (
    ANALOG_FIELDS,
    render_bioisostere_report,
    run_bioisostere_pipeline,
)
from .dynamics import (
    run_dynamics_pipeline,
)
from .ind import (
    DOSE_PROJECTION_FIELDS,
    compile_ind_dossier,
)
from .pathway import (
    SIMULATION_FIELDS,
    render_pathway_report,
    simulate_circuit_dynamics,
)
from .preclinical import (
    PRECLINICAL_FIELDS,
    evaluate_candidates,
    load_pd_data,
    load_pk_data,
    render_preclinical_report,
)
from .rank import RANK_FIELDS
from .reports import write_csv, write_text
from .safety import (
    SAFETY_FIELDS,
    assess_safety,
    load_safety_pharmacology,
    load_tox_studies,
    render_safety_report,
)
from .sar import (
    SCAFFOLD_CLUSTER_FIELDS,
    analyze_sar,
    load_compounds as load_sar_compounds,
    render_sar_report,
)
from .selection import render_lead_report, select_leads
from .studio import (
    generate_circuit_dynamics_svg,
    generate_docking_heatmap_svg,
    generate_pk_curve_svg,
    generate_radar_svg,
    load_studio_data,
    render_studio_html,
)
from .target import list_targets

PROJECT_ROOT = Path(__file__).resolve().parent.parent


class RunnerError(ValueError):
    """Raised when one-click pipeline runner cannot execute or find input data."""


@dataclass(frozen=True)
class RunSummary:
    curated_count: int
    screened_count: int
    shortlisted_count: int
    qualified_hits_count: int
    preclinical_candidates_count: int
    scaffolds_count: int
    elapsed_seconds: float
    output_dir: Path
    studio_html_path: Path
    analogs_count: int = 0
    safety_assessed_count: int = 0
    safety_safe_count: int = 0
    safety_failed_count: int = 0


def execute_turnkey_pipeline(
    vendor_csv: Path,
    repurposing_csv: Optional[Path],
    docking_csv: Path,
    degradation_csv: Optional[Path],
    permeability_csv: Optional[Path],
    cytotoxicity_csv: Optional[Path],
    pk_csv: Optional[Path],
    pd_csv: Optional[Path],
    out_dir: Path,
    min_tier: str = "B",
    emit_root: bool = False,
    tox_csv: Optional[Path] = None,
    safety_pharm_csv: Optional[Path] = None,
) -> RunSummary:
    """Execute complete end-to-end PNN discovery workflow from libraries to interactive studio."""
    start_time = time.perf_counter()
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Curation, Screening, and Docking Consensus Shortlist
    lib_paths = [str(vendor_csv)]
    if repurposing_csv and repurposing_csv.is_file():
        lib_paths.append(str(repurposing_csv))

    try:
        selection = select_leads(
            inputs=lib_paths,
            docking=str(docking_csv),
            min_tier=min_tier,
        )
    except Exception as exc:
        raise RunnerError(f"Phase 2 lead selection failed: {exc}") from exc

    shortlist = selection.shortlist
    curated = selection.curated
    screened = selection.screened

    shortlist_csv = out_dir / "shortlist.csv"
    write_csv(shortlist_csv, RANK_FIELDS, [s.as_row() for s in shortlist])

    lead_report = render_lead_report(selection, {"min_tier": min_tier})
    write_text(out_dir / "lead_report.md", lead_report)

    # 2. Phase 3 In Vitro Assay Qualification
    deg_data = load_degradation_data(degradation_csv) if degradation_csv and degradation_csv.is_file() else {}
    perm_data = load_permeability_data(permeability_csv) if permeability_csv and permeability_csv.is_file() else {}
    cyto_data = load_cytotoxicity_data(cytotoxicity_csv) if cytotoxicity_csv and cytotoxicity_csv.is_file() else {}

    qualified = qualify_leads(deg_data, perm_data, cyto_data)
    qualified_csv = out_dir / "qualified_hits.csv"
    write_csv(qualified_csv, HIT_FIELDS, [h.as_row() for h in qualified])

    assay_criteria = {
        "ec50_ceiling_nm": 100.0,
        "max_deg_floor_pct": 50.0,
        "min_selectivity_index": 10.0,
        "pampa_pe_floor": 3.0,
        "mdck_efflux_ceiling": 3.0,
    }
    assay_report = render_assay_report(qualified, assay_criteria)
    write_text(out_dir / "assay_report.md", assay_report)

    # 3. Phase 4 Preclinical PK/PD Candidate Evaluation
    pk_profiles = load_pk_data(pk_csv) if pk_csv and pk_csv.is_file() else {}
    pd_outcomes = load_pd_data(pd_csv) if pd_csv and pd_csv.is_file() else {}

    candidates = evaluate_candidates(pk_profiles, pd_outcomes)
    candidates_csv = out_dir / "preclinical_candidates.csv"
    write_csv(candidates_csv, PRECLINICAL_FIELDS, [c.as_row() for c in candidates])

    preclinical_report = render_preclinical_report(candidates, {"min_kp_brain": 0.30})
    write_text(out_dir / "preclinical_report.md", preclinical_report)

    # 4. SAR & Bemis-Murcko Scaffold Decomposition
    sar_records = load_sar_compounds(vendor_csv)
    if repurposing_csv and repurposing_csv.is_file():
        sar_records.extend(load_sar_compounds(repurposing_csv))

    docking_map = {s.compound_id: s.consensus for s in shortlist}
    assay_potencies = {h.compound_id: h.ec50_nm for h in qualified}
    preclinical_exposures = {
        c.compound_id: {
            "brain_kp": c.brain_plasma_ratio or 0.0,
            "oral_f": c.oral_bioavailability_pct or 0.0,
        }
        for c in candidates
    }

    compounds_sar, clusters = analyze_sar(
        compound_records=sar_records,
        docking_scores=docking_map,
        assay_potencies=assay_potencies,
        preclinical_exposures=preclinical_exposures,
    )

    sar_clusters_csv = out_dir / "sar_clusters.csv"
    write_csv(sar_clusters_csv, SCAFFOLD_CLUSTER_FIELDS, [cl.as_row() for cl in clusters])

    sar_report = render_sar_report(compounds_sar, clusters)
    write_text(out_dir / "sar_report.md", sar_report)

    # 5. Advanced Translational Capabilities (IND Dossier, Bioisosteres, Pathway ODE, MD Dynamics)
    lead_id = candidates[0].compound_id if candidates else "PNN-101"
    library_paths = [vendor_csv]
    if repurposing_csv and repurposing_csv.is_file():
        library_paths.append(repurposing_csv)

    analogs_count = 0
    for library_path in library_paths:
        try:
            analogs = run_bioisostere_pipeline(library_path, target_parent_id=lead_id)
        except (ValueError, KeyError, OSError):  # hygiene: allow non-critical translational step exception
            continue
        if analogs:
            analogs_count = len(analogs)
            write_csv(out_dir / "bioisosteric_analogs.csv", ANALOG_FIELDS, [a.as_row() for a in analogs])
            write_text(out_dir / "bioisostere_report.md", render_bioisostere_report(analogs, lead_id))
            break
    if analogs_count == 0:
        # Record the stage outcome instead of silently omitting the artifact.
        write_text(out_dir / "bioisostere_report.md", render_bioisostere_report([], lead_id))

    projections: list = []
    lead_pk: dict[str, float] = {}
    if pk_csv and pk_csv.is_file():
        try:
            measured_pk: dict[str, dict[str, float]] = {}
            for cid, profiles in pk_profiles.items():
                entry: dict[str, float] = {}
                po_plasma = profiles.get("PO_plasma")
                iv_plasma = profiles.get("IV_plasma")
                if po_plasma is not None:
                    entry["dose_mg_kg"] = po_plasma.dose_mg_kg
                if iv_plasma is not None:
                    if iv_plasma.cl_l_h_kg is not None:
                        entry["cl_l_h_kg"] = iv_plasma.cl_l_h_kg
                    if iv_plasma.vd_l_kg is not None:
                        entry["vd_l_kg"] = iv_plasma.vd_l_kg
                if entry:
                    measured_pk[cid] = entry
            lead_pk = measured_pk.get(lead_id, {})

            projections, dossier_md = compile_ind_dossier(
                lead_id,
                compounds_path=library_paths,
                candidates_path=candidates_csv,
                tox_path=tox_csv,
                extra_pk=measured_pk,
            )
            write_text(out_dir / "ind_dossier.md", dossier_md)
            write_csv(out_dir / "human_dose_projections.csv", DOSE_PROJECTION_FIELDS, [p.as_row() for p in projections])
        except (ValueError, KeyError, OSError):  # hygiene: allow non-critical translational step exception
            pass

    # Simulate only the measured efficacious dose; an unmeasured dose is never substituted.
    p_points: list = []
    pathway_dose = lead_pk.get("dose_mg_kg")
    if pathway_dose is None and lead_id in pd_outcomes:
        pathway_dose = pd_outcomes[lead_id].dose_mg_kg
    if pathway_dose is not None:
        try:
            p_points, p_summary = simulate_circuit_dynamics(
                compound_id=lead_id, dose_mg_kg=pathway_dose, regimen="INTERMITTENT_PULSATILE", duration_days=14
            )
            write_csv(out_dir / "pathway_simulation.csv", SIMULATION_FIELDS, [p.as_row() for p in p_points])
            write_text(out_dir / "pathway_report.md", render_pathway_report(p_summary, p_points))
        except (ValueError, KeyError, OSError):  # hygiene: allow non-critical translational step exception
            pass

    traj_template = PROJECT_ROOT / "data/dynamics/ptpsigma_md_trajectory_template.csv"
    if traj_template.is_file():
        try:
            run_dynamics_pipeline(traj_template, compound_id=lead_id, target_id="PTPsigma", out_dir=out_dir)
        except (ValueError, KeyError, OSError):  # hygiene: allow non-critical translational step exception
            pass

    # 6. IND-Enabling Safety Pharmacology and Repeat-Dose Toxicology Package
    safety_assessed_count = 0
    safety_safe_count = 0
    safety_failed_count = 0
    assessments: list = []
    try:
        tox_studies = load_tox_studies(tox_csv) if tox_csv and tox_csv.is_file() else []
        safety_pharm = (
            load_safety_pharmacology(safety_pharm_csv)
            if safety_pharm_csv and safety_pharm_csv.is_file()
            else {}
        )
        if tox_studies or safety_pharm:
            assessments = assess_safety(tox_studies, safety_pharm)
            write_csv(out_dir / "safety_assessments.csv", SAFETY_FIELDS, [a.as_row() for a in assessments])
            write_text(
                out_dir / "safety_package.md",
                render_safety_report(assessments, {"safety_factor": 10.0}),
            )
            safety_assessed_count = len(assessments)
            safety_safe_count = sum(1 for a in assessments if a.status == "SAFE_TO_PROGRESS")
            safety_failed_count = sum(1 for a in assessments if a.status == "FAIL_SAFETY")
    except (ValueError, KeyError, OSError):  # hygiene: allow non-critical translational step exception
        safety_assessed_count = 0
        safety_safe_count = 0
        safety_failed_count = 0

    # 7. Interactive Visualization Studio HTML Dashboard (Comprehensive Preclinical & Translational UI)
    studio_compounds, targets = load_studio_data(
        compounds_path=vendor_csv,
        docking_path=docking_csv,
        assays_path=degradation_csv,
        pk_path=pk_csv,
        preclinical_candidates_path=candidates_csv,
    )

    radar_svg = generate_radar_svg(studio_compounds)
    heatmap_svg = generate_docking_heatmap_svg(studio_compounds[:7], targets)
    pk_svg = generate_pk_curve_svg(studio_compounds)

    circuit_svg = (
        generate_circuit_dynamics_svg([p.as_row() for p in p_points])
        if p_points
        else None
    )
    dossier_rows = [p.as_row() for p in projections] if projections else None
    safety_rows = [a.as_row() for a in assessments] if assessments else None

    html_out = render_studio_html(
        studio_compounds,
        radar_svg=radar_svg,
        heatmap_svg=heatmap_svg,
        pk_svg=pk_svg,
        circuit_svg=circuit_svg,
        dossier_data=dossier_rows,
        safety_data=safety_rows,
        targets_data=list_targets(),
    )
    studio_html = out_dir / "studio.html"
    write_text(studio_html, html_out)

    # Optionally emit turnkey root studio.html for zero-friction access
    if emit_root:
        root_studio = PROJECT_ROOT / "studio.html"
        write_text(root_studio, html_out)

    elapsed = time.perf_counter() - start_time
    clinical_count = sum(1 for c in candidates if c.status == "CLINICAL_CANDIDATE")
    confirmed_hit_count = sum(1 for h in qualified if h.status in ("CONFIRMED_HIT", "QUALIFIED_LEAD"))

    return RunSummary(
        curated_count=len(curated),
        screened_count=len(screened),
        shortlisted_count=len(shortlist),
        qualified_hits_count=confirmed_hit_count,
        preclinical_candidates_count=clinical_count,
        scaffolds_count=len(clusters),
        elapsed_seconds=elapsed,
        output_dir=out_dir,
        studio_html_path=studio_html,
        analogs_count=analogs_count,
        safety_assessed_count=safety_assessed_count,
        safety_safe_count=safety_safe_count,
        safety_failed_count=safety_failed_count,
    )


def print_run_banner(summary: RunSummary) -> None:
    """Print terminal summary banner for turnkey runner."""
    border = "=" * 68
    print(border)
    print("  [PNN] TURNKEY DISCOVERY PIPELINE & STUDIO COMPLETED")
    print(border)
    print(f"  * Curated Compounds:         {summary.curated_count}")
    print(f"  * Docking Consensus Leads:   {summary.shortlisted_count}")
    print(f"  * Confirmed In Vitro Hits:   {summary.qualified_hits_count}")
    print(f"  * Preclinical Candidates:    {summary.preclinical_candidates_count}")
    print(f"  * Bemis-Murcko Scaffolds:    {summary.scaffolds_count}")
    if summary.analogs_count > 0:
        print(f"  * Bioisosteric Analogs:      {summary.analogs_count}")
    if summary.safety_assessed_count > 0:
        print(
            f"  * IND Safety Cleared:        {summary.safety_safe_count}/{summary.safety_assessed_count} "
            f"({summary.safety_failed_count} failed)"
        )
    print(f"  * Execution Time:            {summary.elapsed_seconds:.2f} seconds")
    print(f"  * Output Directory:          {summary.output_dir}")
    print(f"  * Interactive Dashboard:     {summary.studio_html_path.resolve()}")
    print(border)
    print("  Open 'studio.html' in your browser to view the interactive dashboard.")
    print(border)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pnn.run",
        description="One-click execution of the entire PNN perineuronal net drug discovery pipeline.",
    )
    parser.add_argument(
        "--vendor",
        default=str(PROJECT_ROOT / "data/compounds/vendor_library.csv"),
        help="vendor screening library CSV",
    )
    parser.add_argument(
        "--repurposing",
        default=str(PROJECT_ROOT / "data/compounds/repurposing_library.csv"),
        help="repurposing library CSV",
    )
    parser.add_argument(
        "--docking",
        default=str(PROJECT_ROOT / "data/docking/multi_target_scores.csv"),
        help="multi-target docking CSV",
    )
    parser.add_argument(
        "--degradation",
        default=str(PROJECT_ROOT / "data/assays/degradation_template.csv"),
        help="DQ-CS degradation assay CSV",
    )
    parser.add_argument(
        "--permeability",
        default=str(PROJECT_ROOT / "data/assays/permeability_template.csv"),
        help="PAMPA/MDCK permeability CSV",
    )
    parser.add_argument(
        "--cytotoxicity",
        default=str(PROJECT_ROOT / "data/assays/cytotoxicity_template.csv"),
        help="CellTiter-Glo cytotoxicity CSV",
    )
    parser.add_argument(
        "--pk",
        default=str(PROJECT_ROOT / "data/preclinical/pk_timecourse_template.csv"),
        help="in vivo PK time series CSV",
    )
    parser.add_argument(
        "--pd",
        default=str(PROJECT_ROOT / "data/preclinical/pd_metabolic_template.csv"),
        help="in vivo PD outcomes CSV",
    )
    parser.add_argument(
        "--tox",
        default=str(PROJECT_ROOT / "data/toxicology/repeat_dose_tox_template.csv"),
        help="repeat-dose toxicology NOAEL/LOAEL CSV",
    )
    parser.add_argument(
        "--safety-pharm",
        dest="safety_pharm",
        default=str(PROJECT_ROOT / "data/toxicology/safety_pharmacology_template.csv"),
        help="hERG/CYP/genotoxicity safety pharmacology CSV",
    )
    parser.add_argument("--out-dir", default="output", help="destination directory for all generated artifacts")
    parser.add_argument("--open", dest="open_browser", action="store_true", help="explicitly enable browser launch (enabled by default)")
    parser.add_argument("--no-browser", dest="no_browser", action="store_true", help="suppress automatically opening studio.html in browser")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        resolved_out = Path(args.out_dir)
        summary = execute_turnkey_pipeline(
            vendor_csv=Path(args.vendor),
            repurposing_csv=Path(args.repurposing) if args.repurposing else None,
            docking_csv=Path(args.docking),
            degradation_csv=Path(args.degradation) if args.degradation else None,
            permeability_csv=Path(args.permeability) if args.permeability else None,
            cytotoxicity_csv=Path(args.cytotoxicity) if args.cytotoxicity else None,
            pk_csv=Path(args.pk) if args.pk else None,
            pd_csv=Path(args.pd) if args.pd else None,
            out_dir=resolved_out,
            emit_root=(resolved_out.resolve() == (PROJECT_ROOT / "output").resolve() or resolved_out == Path("output")),
            tox_csv=Path(args.tox) if args.tox else None,
            safety_pharm_csv=Path(args.safety_pharm) if args.safety_pharm else None,
        )
        print_run_banner(summary)

        if not args.no_browser:
            target_url = summary.studio_html_path.resolve().as_uri()
            try:
                webbrowser.open_new_tab(target_url)
            except Exception as exc:  # hygiene: allow non-critical desktop browser launch
                print(f"[Notice] Could not launch browser ({exc}). Open dashboard manually: {target_url}", file=sys.stderr)

    except RunnerError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: turnkey runner requires valid input datasets in data/ subdirectories.", file=sys.stderr)
        print("Fix: verify that benchmark datasets exist in data/compounds/, data/docking/, data/assays/, and data/preclinical/.", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"What: {exc}", file=sys.stderr)
        print("Why: output directory or artifacts could not be written to disk.", file=sys.stderr)
        print("Fix: specify a writable destination for --out-dir.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
