"""CSV row contracts shared by the PNN Phase 2 command line tools."""

from __future__ import annotations

from typing import Dict, Optional

from .curation import Candidate, property_bin
from .screen import ScreenResult

CURATED_FIELDS = (
    "compound_id",
    "name",
    "smiles",
    "source",
    "tier",
    "score",
    "mw",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "pka",
    "charge",
    "property_bin",
    "failed_criteria",
    "unknown_criteria",
)

SCREEN_FIELDS = (
    "compound_id",
    "name",
    "smiles",
    "tier",
    "score",
    "mw",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "pka",
    "charge",
    "logbb",
    "brain_plasma_ratio",
    "oral_bioavailability",
    "cyp3a4_inhibition",
    "herg_liability",
    "pnn_ec50_nm",
    "failed_criteria",
    "unknown_criteria",
)


def _optional(value: Optional[float]) -> object:
    return "" if value is None else value


def curated_row(candidate: Candidate, width: float) -> Dict[str, object]:
    """Return one curated-library row, including source and property bin."""
    compound = candidate.compound
    return {
        "compound_id": compound.compound_id,
        "name": compound.name,
        "smiles": compound.smiles,
        "source": candidate.source,
        "tier": candidate.result.tier,
        "score": f"{candidate.result.score:.4f}",
        "mw": compound.mw,
        "logp": compound.logp,
        "tpsa": compound.tpsa,
        "hbd": compound.hbd,
        "hba": compound.hba,
        "rotatable_bonds": compound.rotatable_bonds,
        "pka": _optional(compound.pka),
        "charge": compound.charge,
        "property_bin": "-".join(str(cell) for cell in property_bin(compound, width)),
        "failed_criteria": ";".join(candidate.result.failed_criteria()),
        "unknown_criteria": ";".join(candidate.result.unknown_criteria()),
    }


def screen_row(result: ScreenResult) -> Dict[str, object]:
    """Return one ranked-screening row for a screened compound."""
    compound = result.compound
    return {
        "compound_id": compound.compound_id,
        "name": compound.name,
        "smiles": compound.smiles,
        "tier": result.tier,
        "score": f"{result.score:.4f}",
        "mw": compound.mw,
        "logp": compound.logp,
        "tpsa": compound.tpsa,
        "hbd": compound.hbd,
        "hba": compound.hba,
        "rotatable_bonds": compound.rotatable_bonds,
        "pka": _optional(compound.pka),
        "charge": compound.charge,
        "logbb": _optional(compound.logbb),
        "brain_plasma_ratio": _optional(compound.brain_plasma_ratio),
        "oral_bioavailability": _optional(compound.oral_bioavailability),
        "cyp3a4_inhibition": _optional(compound.cyp3a4_inhibition),
        "herg_liability": _optional(compound.herg_liability),
        "pnn_ec50_nm": _optional(compound.pnn_ec50_nm),
        "failed_criteria": ";".join(result.failed_criteria()),
        "unknown_criteria": ";".join(result.unknown_criteria()),
    }
