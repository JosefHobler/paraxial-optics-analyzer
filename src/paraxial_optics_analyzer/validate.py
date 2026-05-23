from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from paraxial_optics_analyzer.paraxial import (
    efl_from_matrix,
    lensmaker_thick,
    trace_paraxial,
)
from paraxial_optics_analyzer.prescription import ObjectSpec, Prescription, Surface
from paraxial_optics_analyzer.raytrace impo rt image_plane_z, trace_system

_BK7_R1 = 50.0
_BK7_R2 = math.inf
_BK7_N = 1.5168
_BK7_T = 5.0


@dataclass(frozen=True)
class CheckResult:
    name: str
    passed: bool
    metric: str   
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


def _cooke_triplet() -> Prescription:
    return Prescription(
        name="cooke-triplet",
        wavelength_um=0.5876,
        units="mm",
        obj=ObjectSpec(distance=math.inf),
        surfaces=(
            Surface(radius=  21.25, thickness= 4.00, n=1.6116, semi_diameter=8.0),
            Surface(radius=-158.65, thickness= 7.76, n=1.0,    semi_diameter=8.0),
            Surface(radius= -20.25, thickness= 1.50, n=1.6053, semi_diameter=6.0),
            Surface(radius=  19.60, thickness= 4.90, n=1.0,    semi_diameter=6.0),
            Surface(radius= 141.25, thickness= 3.00, n=1.6116, semi_diameter=8.0),
            Surface(radius= -18.40, thickness=45.00, n=1.0,    semi_diameter=8.0),
        ),
        stop=3,
    )


def check_lensmaker(threshold: float = 1e-10) -> CheckResult:
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
    aberration = approx. y to the third power)."""
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
    pre = _cooke_triplet()
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
