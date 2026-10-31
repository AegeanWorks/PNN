"""Atomic CSV, JSON and text report writers shared by the PNN command line tools."""

from __future__ import annotations

import csv
import json
import os
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence


def _open_temp(path: Path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path, tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", delete=False, dir=str(path.parent), prefix=f".{path.name}."
    )


def _finish(path: Path, handle) -> None:
    handle.flush()
    os.fsync(handle.fileno())
    handle.close()
    try:
        os.replace(handle.name, path)
    except OSError:
        Path(handle.name).unlink(missing_ok=True)
        raise


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    """Write `rows` as CSV to a temporary file, then replace `path` atomically."""
    destination, handle = _open_temp(path)
    try:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    except BaseException:
        handle.close()
        Path(handle.name).unlink(missing_ok=True)
        raise
    _finish(destination, handle)


def write_json(path: Path, report: Mapping[str, object]) -> None:
    """Write `report` as JSON to a temporary file, then replace `path` atomically."""
    destination, handle = _open_temp(path)
    try:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    except BaseException:
        handle.close()
        Path(handle.name).unlink(missing_ok=True)
        raise
    _finish(destination, handle)


def write_text(path: Path, text: str) -> None:
    """Write `text` to a temporary file, then replace `path` atomically."""
    destination, handle = _open_temp(path)
    try:
        handle.write(text)
    except BaseException:
        handle.close()
        Path(handle.name).unlink(missing_ok=True)
        raise
    _finish(destination, handle)
