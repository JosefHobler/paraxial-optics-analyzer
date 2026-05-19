"""YAML I/O for prescriptions. """
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import yaml

from paraxial_optics_analyzer.prescription import (
    ObjectSpec,
    Prescription,
    PrescriptionError,
    Surface,
)

_SURFACE_REQUIRED = {"radius", "thickness", "n", "semi_diameter"}
_TOP_REQUIRED = {"name", "wavelength_um", "units", "object", "surfaces"}

_INF_TOKENS = {".inf", "+.inf", "inf", "+inf", "infinity", "+infinity"}
_NEG_INF_TOKENS = {"-.inf", "-inf", "-infinity"}


def load_prescription(path: str | Path) -> Prescription:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as e:
        raise PrescriptionError(f"could not read {path}: {e}") from e
    try:
        raw = yaml.safe_load(text)
    except yaml.YAMLError as e:
        raise PrescriptionError(f"failed to parse YAML in {path}: {e}") from e
    if not isinstance(raw, dict):
        raise PrescriptionError(
            f"top-level YAML in {path} must be a mapping, got {type(raw).__name__}"
        )
    return prescription_from_dict(raw)


def prescription_from_dict(raw: dict[str, Any]) -> Prescription:
    _require_keys(raw, _TOP_REQUIRED, "top-level")

    obj_raw = raw["object"]
    if not isinstance(obj_raw, dict):
        raise PrescriptionError(f"'object' must be a mapping, got {type(obj_raw).__name__}")
    _require_keys(obj_raw, {"distance"}, "object")
    obj = ObjectSpec(
        distance=_to_float(obj_raw["distance"], "object.distance"),
        height=_to_float(obj_raw.get("height", 0.0), "object.height"),
    )

    surfs_raw = raw["surfaces"]
    if not isinstance(surfs_raw, list) or not surfs_raw:
        raise PrescriptionError("'surfaces' must be a non-empty list")

    surfaces: list[Surface] = []
    for i, s in enumerate(surfs_raw, start=1):
        if not isinstance(s, dict):
            raise PrescriptionError(f"surface {i} must be a mapping, got {type(s).__name__}")
        _require_keys(s, _SURFACE_REQUIRED, f"surface {i}")
        surfaces.append(Surface(
            radius=_to_float(s["radius"], f"surface {i}.radius"),
            thickness=_to_float(s["thickness"], f"surface {i}.thickness"),
            n=_to_float(s["n"], f"surface {i}.n"),
            semi_diameter=_to_float(s["semi_diameter"], f"surface {i}.semi_diameter"),
        ))

    return Prescription(
        name=str(raw["name"]),
        wavelength_um=_to_float(raw["wavelength_um"], "wavelength_um"),
        units=str(raw["units"]),
        obj=obj,
        surfaces=tuple(surfaces),
        stop=_to_int(raw.get("stop", 1), "stop"),
    )


def _require_keys(d: dict[str, Any], required: set[str], where: str) -> None:
    missing = required - d.keys()
    if missing:
        raise PrescriptionError(f"{where} is missing required keys: {sorted(missing)}")


def _to_float(v: Any, where: str) -> float:
    if isinstance(v, bool):
        raise PrescriptionError(f"{where} must be a number, got bool")
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        s = v.strip().lower()
        if s in _INF_TOKENS:
            return math.inf
        if s in _NEG_INF_TOKENS:
            return -math.inf
        raise PrescriptionError(f"{where} must be a number (or .inf), got string {v!r}")
    raise PrescriptionError(f"{where} must be a number, got {type(v).__name__}")


def _to_int(v: Any, where: str) -> int:
    if isinstance(v, bool):
        raise PrescriptionError(f"{where} must be an integer, got bool")
    if isinstance(v, int):
        return v
    raise PrescriptionError(f"{where} must be an integer, got {type(v).__name__}")
