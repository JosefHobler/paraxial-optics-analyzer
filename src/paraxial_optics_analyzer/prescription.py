"""Centered sequential lens prescription"""
from __future__ import annotations

import math
from dataclasses import dataclass


class PrescriptionError(ValueError):
    pass


_VALID_UNITS = frozenset({"mm", "cm", "m", "inch"})


@dataclass(frozen=True)
class Surface:
    """Refracting surface in a sequential system"""
    radius: float
    thickness: float
    n: float
    semi_diameter: float
 
    def __post_init__(self) -> None:
        r, t, n, sd = self.radius, self.thickness, self.n, self.semi_diameter
        if math.isnan(r):
            raise PrescriptionError("surface radius is NaN")
        if r == 0.0:
            raise PrescriptionError("surface radius cannot be 0 (use .inf for a plano surface)")
        if not math.isfinite(t):
            raise PrescriptionError(f"surface thickness must be finite, got {t}")
        if t < 0.0:
            raise PrescriptionError(f"surface thickness must be >= 0, got {t}")
        if not math.isfinite(n) or n < 1.0:
            raise PrescriptionError(f"refractive index must be finite and >= 1.0, got {n}")
        if not math.isfinite(sd) or sd <= 0.0:
            raise PrescriptionError(f"semi_diameter must be > 0, got {sd}")


@dataclass(frozen=True)
class ObjectSpec:
    distance: float
    height: float = 0.0

    def __post_init__(self) -> None:
        d, h = self.distance, self.height
        if math.isnan(d):
            raise PrescriptionError("object distance is NaN")
        if d <= 0.0:
            raise PrescriptionError(f"object distance must be > 0 (or .inf), got {d}")
        if math.isnan(h) or h < 0.0:
            raise PrescriptionError(f"object height must be >= 0, got {h}")


@dataclass(frozen=True)
class Prescription:
    name: str
    wavelength_um: float
    units: str
    obj: ObjectSpec
    surfaces: tuple[Surface, ...]
    stop: int

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise PrescriptionError("name must be a non-empty string")
        wl = self.wavelength_um
        if (not math.isfinite(wl)) or wl <= 0.0:
            raise PrescriptionError(f"wavelength_um must be finite and > 0, got {wl}")
        if self.units not in _VALID_UNITS:
            raise PrescriptionError(
                f"units must be one of {sorted(_VALID_UNITS)}, got {self.units!r}"
            )
        n = len(self.surfaces)
        if n < 1:
            raise PrescriptionError("prescription must have at least one optical surface")
        if not (1 <= self.stop <= n):
            raise PrescriptionError(
                f"stop must be a 1-indexed surface in [1, {n}], got {self.stop}"
            )

    @property
    def n_surfaces(self) -> int:
        return len(self.surfaces)
