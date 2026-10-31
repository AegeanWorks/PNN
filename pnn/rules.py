"""Property rules that turn the PNN BBB/CNS design criteria into verdicts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from .compounds import Compound

PASS = "PASS"
FAIL = "FAIL"
UNKNOWN = "UNKNOWN"

CORE_CRITERIA = ("mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds")
HARD_REJECT_CRITERIA = (
    "logbb",
    "brain_plasma_ratio",
    "oral_bioavailability",
    "cyp3a4_inhibition",
    "herg_liability",
)
CORE_WEIGHT = 1.0
OPTIONAL_WEIGHT = 0.5
KNOWN_CHARGES = frozenset({"neutral", "basic", "acidic", "zwitterion"})


@dataclass(frozen=True)
class CriterionResult:
    """One rule's verdict for one compound."""

    name: str
    status: str
    score: float
    value: Optional[float]
    detail: str


def _range(
    name: str,
    value: float,
    low: float,
    high: float,
    optimal_low: float,
    optimal_high: float,
    unit: str = "",
) -> CriterionResult:
    if value < low or value > high:
        return CriterionResult(
            name=name,
            status=FAIL,
            score=0.0,
            value=value,
            detail=f"{value:g}{unit} is outside {low:g}-{high:g}{unit}",
        )
    optimal = optimal_low <= value <= optimal_high
    target = f"{optimal_low:g}-{optimal_high:g}{unit}" if optimal_low != low or optimal_high != high else f"<= {high:g}{unit}"
    return CriterionResult(
        name=name,
        status=PASS,
        score=1.0 if optimal else 0.7,
        value=value,
        detail=f"{value:g}{unit} is acceptable; optimal is {target}",
    )


def _ceiling(name: str, value: float, limit: float, optimal: float, unit: str = "") -> CriterionResult:
    if value > limit:
        return CriterionResult(
            name=name,
            status=FAIL,
            score=0.0,
            value=value,
            detail=f"{value:g}{unit} exceeds the limit of {limit:g}{unit}",
        )
    return CriterionResult(
        name=name,
        status=PASS,
        score=1.0 if value <= optimal else 0.7,
        value=value,
        detail=f"{value:g}{unit} is within the limit of {limit:g}{unit} (optimal <= {optimal:g}{unit})",
    )


def _physicochemical(compound: Compound) -> List[CriterionResult]:
    return [
        _range("mw", compound.mw, 200.0, 450.0, 300.0, 400.0, " Da"),
        _range("logp", compound.logp, 1.5, 3.5, 2.0, 3.0),
        _ceiling("tpsa", compound.tpsa, 70.0, 60.0, " A2"),
        _ceiling("hbd", compound.hbd, 3.0, 3.0),
        _ceiling("hba", compound.hba, 7.0, 7.0),
        _ceiling("rotatable_bonds", compound.rotatable_bonds, 7.0, 7.0),
    ]


def _pka(compound: Compound) -> CriterionResult:
    if compound.pka is None:
        return CriterionResult("pka", UNKNOWN, 0.0, None, "no pKa supplied")
    if compound.charge == "basic" and not 7.0 <= compound.pka <= 9.0:
        return CriterionResult("pka", FAIL, 0.0, compound.pka, "basic compound outside pKa 7-9")
    if compound.charge == "acidic" and not 3.0 <= compound.pka <= 5.0:
        return CriterionResult("pka", FAIL, 0.0, compound.pka, "acidic compound outside pKa 3-5")
    if compound.charge not in ("basic", "acidic") and not (
        7.0 <= compound.pka <= 9.0 or 3.0 <= compound.pka <= 5.0
    ):
        return CriterionResult("pka", FAIL, 0.0, compound.pka, "pKa outside 3-5 and 7-9")
    return CriterionResult("pka", PASS, 1.0, compound.pka, "pKa fits the CNS window")


def _charge(compound: Compound) -> CriterionResult:
    if not compound.charge:
        return CriterionResult("charge", UNKNOWN, 0.0, None, "no ionization class supplied")
    if compound.charge not in KNOWN_CHARGES:
        return CriterionResult(
            "charge", UNKNOWN, 0.0, None, f"unrecognized ionization class: {compound.charge}"
        )
    if compound.charge in ("neutral", "basic"):
        return CriterionResult("charge", PASS, 1.0, None, f"{compound.charge} at pH 7.4")
    return CriterionResult("charge", PASS, 0.4, None, f"{compound.charge} is less favorable than neutral/basic")


def _logbb(compound: Compound) -> CriterionResult:
    value = compound.logbb
    if value is None:
        return CriterionResult("logbb", UNKNOWN, 0.0, None, "no logBB supplied")
    if value < -1.0:
        return CriterionResult("logbb", FAIL, 0.0, value, f"logBB {value:g} is below -1.0")
    return CriterionResult(
        "logbb",
        PASS,
        1.0 if value >= -0.5 else 0.7,
        value,
        f"logBB {value:g} (preferred >= -0.5)",
    )


def _brain_plasma(compound: Compound) -> CriterionResult:
    value = compound.brain_plasma_ratio
    if value is None:
        return CriterionResult("brain_plasma_ratio", UNKNOWN, 0.0, None, "no brain/plasma ratio supplied")
    if value < 0.1:
        return CriterionResult("brain_plasma_ratio", FAIL, 0.0, value, f"ratio {value:g} is below 0.1")
    return CriterionResult(
        "brain_plasma_ratio",
        PASS,
        1.0 if value >= 0.3 else 0.5,
        value,
        f"ratio {value:g} (progression target >= 0.3)",
    )


def _oral_bioavailability(compound: Compound) -> CriterionResult:
    value = compound.oral_bioavailability
    if value is None:
        return CriterionResult("oral_bioavailability", UNKNOWN, 0.0, None, "no oral F supplied")
    percent = value * 100.0 if value <= 1.0 else value
    if percent < 10.0:
        return CriterionResult("oral_bioavailability", FAIL, 0.0, percent, f"oral F {percent:g}% is below 10%")
    return CriterionResult(
        "oral_bioavailability",
        PASS,
        1.0 if percent >= 20.0 else 0.7 if percent >= 15.0 else 0.4,
        percent,
        f"oral F {percent:g}% (progression target >= 20%)",
    )


def _cyp3a4(compound: Compound) -> CriterionResult:
    value = compound.cyp3a4_inhibition
    if value is None:
        return CriterionResult("cyp3a4_inhibition", UNKNOWN, 0.0, None, "no CYP3A4 inhibition supplied")
    if value >= 50.0:
        return CriterionResult("cyp3a4_inhibition", FAIL, 0.0, value, f"{value:g}% inhibition at 10 uM")
    return CriterionResult(
        "cyp3a4_inhibition", PASS, 1.0 if value < 30.0 else 0.7, value, f"{value:g}% inhibition at 10 uM"
    )


def _herg(compound: Compound) -> CriterionResult:
    value = compound.herg_liability
    if value is None:
        return CriterionResult("herg_liability", UNKNOWN, 0.0, None, "no hERG liability supplied")
    if value >= 30.0:
        return CriterionResult("herg_liability", FAIL, 0.0, value, f"{value:g}% hERG inhibition at 10 uM")
    return CriterionResult("herg_liability", PASS, 1.0 if value < 10.0 else 0.7, value, f"{value:g}% hERG inhibition")


def _pnn_ec50(compound: Compound) -> CriterionResult:
    value = compound.pnn_ec50_nm
    if value is None:
        return CriterionResult("pnn_ec50_nm", UNKNOWN, 0.0, None, "no PNN-softening EC50 supplied")
    if value >= 1000.0:
        return CriterionResult("pnn_ec50_nm", FAIL, 0.0, value, f"EC50 {value:g} nM is above 1000 nM")
    return CriterionResult(
        "pnn_ec50_nm",
        PASS,
        1.0 if value < 10.0 else 0.7 if value < 100.0 else 0.4,
        value,
        f"EC50 {value:g} nM (progression target < 100 nM)",
    )


def evaluate(compound: Compound) -> List[CriterionResult]:
    """Return every rule verdict for `compound` in a stable order."""
    return [
        *_physicochemical(compound),
        _pka(compound),
        _charge(compound),
        _logbb(compound),
        _brain_plasma(compound),
        _oral_bioavailability(compound),
        _cyp3a4(compound),
        _herg(compound),
        _pnn_ec50(compound),
    ]


def weight(criterion_name: str) -> float:
    """Core design properties outweigh optional measured/ADMET properties."""
    return CORE_WEIGHT if criterion_name in CORE_CRITERIA else OPTIONAL_WEIGHT
