"""Paraxial trace + ABCD matrix. Cross-checked against the lensmaker eq."""
from __future__ import annotations

import math
from dataclasses import dataclass

from paraxial_optics_analyzer.prescription import Prescription, Surface

N_OBJ = 1.0


@dataclass(frozen=True)
class ParaxialResult:
    efl: float
    bfl: float
    image_distance: float
    f_number: float
    last_medium_index: float


def trace_paraxial(pre: Prescription) -> ParaxialResult:
    n_last = pre.surfaces[-1].n

    ys, us = _trace_parallel_ray(pre)
    u_out = us[-1]
    y_out = ys[-1]

    if u_out == 0.0:
        efl = bfl = math.inf
    else:
        efl = -1.0 / u_out
        bfl = -y_out / u_out

    if math.isinf(pre.obj.distance):
        s_img = bfl
    else:
        s_img = _image_distance_finite_conjugate(pre)

    stop = pre.surfaces[pre.stop - 1]
    D_pupil = 2.0 * stop.semi_diameter
    fno = abs(efl) / D_pupil if (D_pupil > 0.0 and math.isfinite(efl)) else math.inf

    return ParaxialResult(efl, bfl, s_img, fno, n_last)


def system_matrix(pre: Prescription) -> tuple[float, float, float, float]:
    """ABCD in the reduced-angle [y, n*u] convention. Last-surface to image transfer is not applied."""
    A, B, C, D = 1.0, 0.0, 0.0, 1.0
    n_prev = N_OBJ

    surfs = pre.surfaces
    n_surf = len(surfs)
    for i, s in enumerate(surfs):
        P = 0.0 if math.isinf(s.radius) else (s.n - n_prev) / s.radius
        C, D = -P * A + C, -P * B + D

        if i < n_surf - 1:
            tn = s.thickness / s.n
            A, B = A + tn * C, B + tn * D

        n_prev = s.n

    return A, B, C, D


def efl_from_matrix(pre: Prescription) -> float:
    _, _, C, _ = system_matrix(pre)
    return math.inf if C == 0.0 else -pre.surfaces[-1].n / C


def bfl_from_matrix(pre: Prescription) -> float:
    A, _, C, _ = system_matrix(pre)
    return math.inf if C == 0.0 else -A * pre.surfaces[-1].n / C


def lensmaker_thick(R1: float, R2: float, n: float, thickness: float) -> float:
    if n <= 1.0:
        raise ValueError(f"refractive index must be > 1 for a lens in air, got {n}")
    iR1 = 0.0 if math.isinf(R1) else 1.0 / R1
    iR2 = 0.0 if math.isinf(R2) else 1.0 / R2
    if math.isinf(R1) or math.isinf(R2):
        thk = 0.0
    else:
        thk = (n - 1.0) * thickness / (n * R1 * R2)
    P = (n - 1.0) * (iR1 - iR2 + thk)
    return math.inf if P == 0.0 else 1.0 / P




def _refract(y: float, u: float, surf: Surface, n_before: float) -> tuple[float, float]:
    P = 0.0 if math.isinf(surf.radius) else (surf.n - n_before) / surf.radius
    return y, (n_before * u - y * P) / surf.n


def _trace_parallel_ray(pre: Prescription) -> tuple[list[float], list[float]]:
    ys: list[float] = []
    us: list[float] = []
    n_before = N_OBJ
    y, u = 1.0, 0.0
    last = len(pre.surfaces) - 1
    for i, surf in enumerate(pre.surfaces):
        ys.append(y)
        y, u = _refract(y, u, surf, n_before)
        us.append(u)
        if i < last:
            y += surf.thickness * u
        n_before = surf.n
    return ys, us


def _image_distance_finite_conjugate(pre: Prescription) -> float:
    n_before = N_OBJ
    y, u = 0.0, 1.0
    y += pre.obj.distance * u
    last = len(pre.surfaces) - 1
    for i, surf in enumerate(pre.surfaces):
        y, u = _refract(y, u, surf, n_before)
        if i < last:
            y += surf.thickness * u
        n_before = surf.n
    return math.inf if u == 0.0 else -y / u
