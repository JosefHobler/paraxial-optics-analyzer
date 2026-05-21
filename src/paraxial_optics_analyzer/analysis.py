"""Image-quality analyses: spot diagram, ray fan, best-focus search."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from paraxial_optics_analyzer.prescription import Prescription
from paraxial_optics_analyzer.raytrace import TraceError, trace_system
from paraxial_optics_analyzer.sampling import (
    hexapolar_pupil,
    launch_parallel,
    linear_pupil,
    pupil_radius_of,
)


@dataclass(frozen=True)
class SpotDiagram:
    points: np.ndarray
    centroid: np.ndarray
    rms: float
    field_angle_rad: float
    image_plane_offset: float
    n_failed: int


def spot_diagram(
    pre: Prescription,
    field_angle_rad: float = 0.0,
    *,
    n_rings: int = 6,
    image_plane_offset: float = 0.0,
) -> SpotDiagram:
    """Hex-sampled bundle from infinity -> image-plane spot statistics."""
    pupil = hexapolar_pupil(n_rings, pupil_radius_of(pre))
    positions, directions = launch_parallel(pupil, field_angle_rad, "y")

    pts: list[np.ndarray] = []
    n_failed = 0
    for P, d in zip(positions, directions, strict=True):
        try:
            r = trace_system(P, d, pre, image_plane_offset=image_plane_offset)
            pts.append(r.image_point[:2])
        except TraceError:
            n_failed += 1

    arr = np.array(pts) if pts else np.empty((0, 2))
    centroid = arr.mean(axis=0) if len(arr) else np.array([math.nan, math.nan])
    return SpotDiagram(
        points=arr,
        centroid=centroid,
        rms=_rms_radius(arr, centroid),
        field_angle_rad=field_angle_rad,
        image_plane_offset=image_plane_offset,
        n_failed=n_failed,
    )


def _rms_radius(pts: np.ndarray, centroid: np.ndarray) -> float:
    if len(pts) == 0:
        return math.nan
    delta = pts - centroid
    return float(math.sqrt((delta * delta).sum(axis=1).mean()))


@dataclass(frozen=True)
class RayFan:
    pupil_coords: np.ndarray
    transverse: np.ndarray
    chief_image_xy: np.ndarray
    axis: str


def ray_fan(
    pre: Prescription,
    field_angle_rad: float = 0.0,
    *,
    axis: str = "tangential",
    n_samples: int = 41,
    image_plane_offset: float = 0.0,
) -> RayFan:
    """Transverse aberration along one pupil meridian, chief-ray referenced."""
    if axis == "tangential":
        pupil = linear_pupil(n_samples, pupil_radius_of(pre), axis="y")
        comp = 1
    elif axis == "sagittal":
        pupil = linear_pupil(n_samples, pupil_radius_of(pre), axis="x")
        comp = 0
    else:
        raise ValueError(f"axis must be 'tangential' or 'sagittal', got {axis!r}")

    positions, directions = launch_parallel(pupil, field_angle_rad, "y")

    chief_pos, chief_dir = launch_parallel(np.array([[0.0, 0.0]]), field_angle_rad, "y")
    chief = trace_system(chief_pos[0], chief_dir[0], pre, image_plane_offset=image_plane_offset)
    chief_xy = chief.image_point[:2]

    tv = np.full(n_samples, math.nan)
    for i, (P, d) in enumerate(zip(positions, directions, strict=True)):
        try:
            r = trace_system(P, d, pre, image_plane_offset=image_plane_offset)
            tv[i] = r.image_point[comp] - chief_xy[comp]
        except TraceError:
            pass

    return RayFan(pupil_coords=pupil[:, comp], transverse=tv, chief_image_xy=chief_xy, axis=axis)


@dataclass(frozen=True)
class BestFocusResult:
    image_plane_offset: float
    rms_at_best: float
    rms_at_nominal: float


def find_best_focus(
    pre: Prescription,
    field_angle_rad: float = 0.0,
    *,
    n_rings: int = 6,
    search_range: tuple[float, float] | None = None,
    tol: float = 1e-7,
) -> BestFocusResult:
    """Golden-section search for the image-plane offset minimizing RMS spot.

    Default search range is ±10% of the paraxial EFL.
    """
    if search_range is None:
        from paraxial_optics_analyzer.paraxial import trace_paraxial  # avoid cycle
        para = trace_paraxial(pre)
        margin = max(1.0, 0.1 * abs(para.efl))
        search_range = (-margin, margin)

    def rms_at(off: float) -> float:
        return spot_diagram(pre, field_angle_rad, n_rings=n_rings, image_plane_offset=off).rms

    a, b = search_range
    if b <= a:
        raise ValueError(f"search_range must be increasing, got ({a}, {b})")

    phi = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - phi * (b - a)
    d = a + phi * (b - a)
    f_c, f_d = rms_at(c), rms_at(d)
    while (b - a) > tol:
        if f_c < f_d:
            b, d, f_d = d, c, f_c
            c = b - phi * (b - a)
            f_c = rms_at(c)
        else:
            a, c, f_c = c, d, f_d
            d = a + phi * (b - a)
            f_d = rms_at(d)

    return BestFocusResult(
        image_plane_offset=0.5 * (a + b),
        rms_at_best=min(f_c, f_d),
        rms_at_nominal=rms_at(0.0),
    )
