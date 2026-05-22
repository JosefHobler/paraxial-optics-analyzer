"""Image-quality analyses: spot diagram, ray fan, best-focus search."""
from __future__ import annotations

import math
from collections.abc import Callable
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
    n_rings: int = 16,
    image_plane_offset: float = 0.0,
) -> SpotDiagram:
    """Hex-sampled bundle from infinity to image-plane spot statistics."""
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
    n_rings: int = 16,
    search_range: tuple[float, float] | None = None,
    tol: float = 1e-7,
) -> BestFocusResult:
    """Search for the image-plane offset minimizing RMS spot.

    Default search range is +/-10% of the paraxial EFL around the paraxial
    focus, not around the nominal image plane. A coarse grid brackets the best
    local basin before golden-section polishing, so narrow minima near paraxial
    focus are not discarded by the unimodality assumption.
    """
    from paraxial_optics_analyzer.paraxial import trace_paraxial  # avoid cycle

    para = trace_paraxial(pre)
    paraxial_offset = _paraxial_focus_offset(pre, para.image_distance)

    if search_range is None:
        margin = max(1.0, 0.1 * abs(para.efl))
        search_range = (paraxial_offset - margin, paraxial_offset + margin)

    a, b = search_range
    if b <= a:
        raise ValueError(f"search_range must be increasing, got ({a}, {b})")

    rms_cache: dict[float, float] = {}

    def rms_at(off: float) -> float:
        if off not in rms_cache:
            rms_cache[off] = spot_diagram(
                pre, field_angle_rad, n_rings=n_rings, image_plane_offset=off,
            ).rms
        return rms_cache[off]

    grid = np.linspace(a, b, 81)
    grid_rms = np.array([rms_at(float(off)) for off in grid])
    finite = np.isfinite(grid_rms)
    if not finite.any():
        raise TraceError("no valid rays in best-focus search range")

    finite_indices = np.flatnonzero(finite)
    best_grid_index = int(finite_indices[np.argmin(grid_rms[finite])])

    candidates: list[tuple[float, float]] = [
        (float(grid[best_grid_index]), float(grid_rms[best_grid_index])),
        (paraxial_offset, rms_at(paraxial_offset)),
        (0.0, rms_at(0.0)),
    ]

    if 0 < best_grid_index < len(grid) - 1:
        candidates.append(
            _golden_section_minimize(
                rms_at,
                float(grid[best_grid_index - 1]),
                float(grid[best_grid_index + 1]),
                tol,
            )
        )

    best_off, best_rms = min(
        ((off, rms) for off, rms in candidates if math.isfinite(rms)),
        key=lambda item: item[1],
    )

    return BestFocusResult(
        image_plane_offset=best_off,
        rms_at_best=best_rms,
        rms_at_nominal=rms_at(0.0),
    )


def _paraxial_focus_offset(pre: Prescription, image_distance: float) -> float:
    from paraxial_optics_analyzer.raytrace import image_plane_z

    z_nominal = image_plane_z(pre)
    z_last = z_nominal - pre.surfaces[-1].thickness
    return (z_last + image_distance) - z_nominal


def _golden_section_minimize(
    func: Callable[[float], float],
    a: float,
    b: float,
    tol: float,
) -> tuple[float, float]:
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - phi * (b - a)
    d = a + phi * (b - a)
    f_c, f_d = func(c), func(d)
    while (b - a) > tol:
        if f_c < f_d:
            b, d, f_d = d, c, f_c
            c = b - phi * (b - a)
            f_c = func(c)
        else:
            a, c, f_c = c, d, f_d
            d = a + phi * (b - a)
            f_d = func(d)

    if f_c <= f_d:
        return c, f_c
    return d, f_d
