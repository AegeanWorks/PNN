"""Update README directory tree with translational modules and tests."""

from pathlib import Path

README_PATH = Path(__file__).resolve().parents[1] / "README.md"

NEW_TREE = """```text
PNN/
├── README.md                      <- This file - project overview
├── .gitignore                     <- AIKit transient-state baseline
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
│   └── ind_filing.md              <- Human allometric scaling and IND dossier design
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
    ├── test_sar.py                <- Bemis-Murcko scaffold extraction and SAR tests
    ├── test_screen.py             <- Tiering, ranking and CSV ingestion
    ├── test_studio.py             <- Interactive HTML/SVG studio tests
    └── test_target.py             <- Receptor profiles and Vina config generation
```"""


def update_tree() -> None:
    content = README_PATH.read_text(encoding="utf-8")
    start_tag = "```text\n"
    end_tag = "\n```\n"

    start_idx = content.find(start_tag)
    if start_idx == -1:
        raise ValueError("Start of tree block not found")

    end_idx = content.find(end_tag, start_idx)
    if end_idx == -1:
        raise ValueError("End of tree block not found")

    new_content = content[:start_idx] + NEW_TREE + content[end_idx + len(end_tag) - 1 :]
    README_PATH.write_text(new_content, encoding="utf-8")
    print(f"Updated README tree. New length: {len(new_content)} characters.")


if __name__ == "__main__":
    update_tree()
