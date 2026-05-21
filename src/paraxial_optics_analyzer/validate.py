"""Self-checks the user can run from the CLI. Each one is a tiny end-to-end
physics validation, not just a unit test."""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.paraxial import (
    efl_from_matrix,
    lensmaker_thick,
    trace_paraxial,
)
from paraxial_optics_analyzer.prescription import ObjectSpec, Prescription, Surface
from paraxial_optics_analyzer.raytrace import image_plane_z, trace_system

_BK7_R1 = 50.0
_BK7_R2 = math.inf
_BK7_N = 1.5168
_BK7_T = 5.0


@dataclass(frozen=True)
class CheckResult:
    name: str
    passed: bool
    metric: str   # short description, e.g. "relative error"
    value: float
    threshold: float

    def format_line(self) -> str:
        verdict = "PASS" if self.passed else "FAIL"
        v = "< 1e-15" if self.value == 0.0 else f"{self.value:.2g}"
        return f"{self.name}: {verdict}, {self.metric} {v}"


def _bk7_singlet(*, semi_d: float = 12.0) -> Prescription:
    return Prescription(
        name="bk7-singlet",
        wavelength_um=0.5876,
        units="mm",
        obj=ObjectSpec(distance=math.inf),
        surfaces=(
            Surface(radius=_BK7_R1, thickness=_BK7_T, n=_BK7_N, semi_diameter=semi_d),
            Surface(radius=_BK7_R2, thickness=100.0, n=1.0, semi_diameter=semi_d),
        ),
        stop=1,
    )


def check_lensmaker(threshold: float = 1e-10) -> CheckResult:
    """Paraxial trace EFL vs the thick-lens lensmaker equation."""
    pre = _bk7_singlet()
    f_trace = trace_paraxial(pre).efl
    f_lens = lensmaker_thick(_BK7_R1, _BK7_R2, _BK7_N, _BK7_T)
    rel = abs(f_trace - f_lens) / abs(f_lens)
    return CheckResult(
        name="Lensmaker validation",
        passed=rel <= threshold,
        metric="relative error",
        value=rel,
        threshold=threshold,
    )


def check_paraxial_limit(threshold: float = 1e-7) -> CheckResult:
    """Small-aperture real trace must land at the paraxial focus (spherical
    aberration ~ y³)."""
    pre = _bk7_singlet()
    para = trace_paraxial(pre)
    z_img = image_plane_z(pre)
    z_last = z_img - pre.surfaces[-1].thickness
    defocus = (z_last + para.bfl) - z_img

    worst = 0.0
    for y in (1e-3, 1e-4):
        r = trace_system(
            np.array([0.0, y, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            pre,
            image_plane_offset=defocus,
        )
        worst = max(worst, abs(float(r.image_point[1])))
    return CheckResult(
        name="Paraxial-limit validation",
        passed=worst <= threshold,
        metric="max deviation",
        value=worst,
        threshold=threshold,
    )


def check_cooke_triplet(threshold: float = 1e-10) -> CheckResult:
    """Cooke triplet: direct paraxial trace EFL vs ABCD-matrix EFL must agree.

    Two independent code paths over a multi-element system — a stronger
    cross-check than the singlet ones above.
    """
    path = Path(__file__).resolve().parents[2] / "examples" / "cooke_triplet.yaml"
    pre = load_prescription(path)
    f_trace = trace_paraxial(pre).efl
    f_matrix = efl_from_matrix(pre)
    rel = abs(f_trace - f_matrix) / abs(f_matrix)
    return CheckResult(
        name="Cooke triplet EFL check",
        passed=rel <= threshold,
        metric="relative error",
        value=rel,
        threshold=threshold,
    )


def run_all() -> list[CheckResult]:
    return [check_lensmaker(), check_paraxial_limit(), check_cooke_triplet()]
