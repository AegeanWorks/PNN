"""Compound records and CSV ingestion for the PNN screening pipeline."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional

REQUIRED_FIELDS = ("compound_id", "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds")
OPTIONAL_NUMERIC_FIELDS = (
    "pka",
    "logbb",
    "brain_plasma_ratio",
    "oral_bioavailability",
    "cyp3a4_inhibition",
    "herg_liability",
    "pnn_ec50_nm",
)
WHOLE_NUMBER_FIELDS = frozenset({"hbd", "hba", "rotatable_bonds"})
NON_NEGATIVE_FIELDS = frozenset(
    {
        "mw",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "brain_plasma_ratio",
        "oral_bioavailability",
        "cyp3a4_inhibition",
        "herg_liability",
        "pnn_ec50_nm",
    }
)

FIELD_ALIASES = {
    "id": "compound_id",
    "compound": "compound_id",
    "compoundid": "compound_id",
    "molecule": "compound_id",
    "moleculeid": "compound_id",
    "molecularweight": "mw",
    "molweight": "mw",
    "mwda": "mw",
    "clogp": "logp",
    "logp": "logp",
    "xlogp": "logp",
    "polararea": "tpsa",
    "psa": "tpsa",
    "tpsaa2": "tpsa",
    "hbonddonors": "hbd",
    "hbdcount": "hbd",
    "donors": "hbd",
    "hbondacceptors": "hba",
    "hbacount": "hba",
    "acceptors": "hba",
    "rotors": "rotatable_bonds",
    "numrotatablebonds": "rotatable_bonds",
    "rotatablebonds": "rotatable_bonds",
    "logbb": "logbb",
    "bbratio": "brain_plasma_ratio",
    "braintoplasma": "brain_plasma_ratio",
    "brainplasmaratio": "brain_plasma_ratio",
    "brainplasma": "brain_plasma_ratio",
    "oralbioavailability": "oral_bioavailability",
    "oralf": "oral_bioavailability",
    "bioavailability": "oral_bioavailability",
    "fpercent": "oral_bioavailability",
    "cyp3a4": "cyp3a4_inhibition",
    "cyp3a4inhibition": "cyp3a4_inhibition",
    "cyp3a4inhibitionpercent": "cyp3a4_inhibition",
    "herg": "herg_liability",
    "herginhibition": "herg_liability",
    "hergliability": "herg_liability",
    "pnnec50": "pnn_ec50_nm",
    "ec50nm": "pnn_ec50_nm",
    "pnnsofteningec50nm": "pnn_ec50_nm",
}

_CANONICAL_FIELDS = REQUIRED_FIELDS + OPTIONAL_NUMERIC_FIELDS + ("name", "smiles", "charge")
_CANONICAL_BY_TOKEN = {
    "".join(ch for ch in field.lower() if ch.isalnum()): field for field in _CANONICAL_FIELDS
}

_CHARGE_ALIASES = {
    "neutral": "neutral",
    "neutralion": "neutral",
    "basic": "basic",
    "weaklybasic": "basic",
    "base": "basic",
    "acidic": "acidic",
    "acid": "acidic",
    "zwitterion": "zwitterion",
    "zwitterionic": "zwitterion",
}


class CompoundError(ValueError):
    """Raised when input cannot be interpreted as a compound record."""


@dataclass(frozen=True)
class Compound:
    """One candidate with the properties the screening rules consume."""

    compound_id: str
    mw: float
    logp: float
    tpsa: float
    hbd: float
    hba: float
    rotatable_bonds: float
    name: str = ""
    smiles: str = ""
    pka: Optional[float] = None
    charge: str = ""
    logbb: Optional[float] = None
    brain_plasma_ratio: Optional[float] = None
    oral_bioavailability: Optional[float] = None
    cyp3a4_inhibition: Optional[float] = None
    herg_liability: Optional[float] = None
    pnn_ec50_nm: Optional[float] = None


def normalize_field(name: str) -> Optional[str]:
    """Map one CSV header to its canonical field name, or None when the column is not used."""
    canonical = "".join(ch for ch in name.strip().lower() if ch.isalnum())
    if canonical in FIELD_ALIASES:
        return FIELD_ALIASES[canonical]
    if canonical in _CANONICAL_BY_TOKEN:
        return _CANONICAL_BY_TOKEN[canonical]
    if canonical == "compoundname":
        return "name"
    if canonical == "canonicalsmiles":
        return "smiles"
    if canonical in ("ionization", "ionclass"):
        return "charge"
    return None


def _number(value: str, field: str, where: str) -> float:
    text = value.strip().replace(",", "")
    try:
        number = float(text)
    except ValueError as exc:
        raise CompoundError(f"{where}: field '{field}' is not a number: {value!r}") from exc
    if not math.isfinite(number):
        raise CompoundError(f"{where}: field '{field}' must be finite, got {value!r}")
    if field in NON_NEGATIVE_FIELDS and number < 0:
        raise CompoundError(f"{where}: field '{field}' must not be negative, got {value!r}")
    if field in WHOLE_NUMBER_FIELDS and not number.is_integer():
        raise CompoundError(f"{where}: field '{field}' must be a whole number, got {value!r}")
    return number


def _optional_number(values: Mapping[str, str], field: str, where: str) -> Optional[float]:
    raw = values.get(field, "").strip()
    return None if not raw else _number(raw, field, where)


def _charge(raw: str) -> str:
    canonical = "".join(ch for ch in raw.strip().lower() if ch.isalnum())
    return _CHARGE_ALIASES.get(canonical, raw.strip().lower())


def _build_compound(values: Mapping[str, str], where: str) -> Compound:
    compound_id = values.get("compound_id", "").strip()
    if not compound_id:
        raise CompoundError(f"{where}: field 'compound_id' must not be empty")
    missing = [field for field in REQUIRED_FIELDS if not values.get(field, "").strip()]
    if missing:
        raise CompoundError(f"{where}: missing required field(s): {', '.join(missing)}")
    numeric = {field: _number(values[field], field, where) for field in REQUIRED_FIELDS[1:]}
    optional = {
        field: _optional_number(values, field, where) for field in OPTIONAL_NUMERIC_FIELDS
    }
    return Compound(
        compound_id=compound_id,
        name=values.get("name", "").strip(),
        smiles=values.get("smiles", "").strip(),
        charge=_charge(values.get("charge", "")),
        **numeric,
        **optional,
    )


def parse_compounds(
    rows: Iterable[Mapping[str, str]], source: str = "input", allow_duplicates: bool = False
) -> List[Compound]:
    """Turn raw CSV rows into compound records; raise CompoundError on the first bad row."""
    compounds: List[Compound] = []
    seen: Dict[str, int] = {}
    for row_number, raw in enumerate(rows, start=2):
        values: Dict[str, str] = {}
        for key, value in raw.items():
            if key is None:
                continue
            field = normalize_field(key)
            if field is not None:
                values[field] = "" if value is None else str(value).strip()
        where = f"{source} row {row_number}"
        compound = _build_compound(values, where)
        if not allow_duplicates:
            if compound.compound_id in seen:
                raise CompoundError(
                    f"{where}: duplicate compound_id '{compound.compound_id}' "
                    f"(first seen on row {seen[compound.compound_id]})"
                )
            seen[compound.compound_id] = row_number
        compounds.append(compound)
    if not compounds:
        raise CompoundError(f"{source}: no compound rows found")
    return compounds


def load_compounds(path: Path, allow_duplicates: bool = False) -> List[Compound]:
    """Read a UTF-8 CSV file and return its compound records."""
    path = Path(path)
    if not path.is_file():
        raise CompoundError(f"input file not found: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise CompoundError(f"{path}: file has no header row")
            rows = list(reader)
    except UnicodeDecodeError as exc:
        raise CompoundError(f"{path}: file is not valid UTF-8") from exc
    except csv.Error as exc:
        raise CompoundError(f"{path}: malformed CSV: {exc}") from exc
    except OSError as exc:
        raise CompoundError(f"{path}: cannot read file: {exc}") from exc
    return parse_compounds(rows, source=str(path), allow_duplicates=allow_duplicates)
