"""Pupil + field sampling -> ray launches"""
from __future__ import annotations

import math

import numpy as np

from paraxial_optics_analyzer.prescription import Prescription


def hexapolar_pupil(n_rings: int, pupil_radius: float) -> np.ndarray:
    """Hexapolar sampling: centre + 6·r points on ring r. Count = 1 + 3·N·(N+1)."""
    if n_rings < 0:
        raise ValueError("n_rings must be >= 0")
    if pupil_radius <= 0.0:
        raise ValueError("pupil_radius must be > 0")
    pts: list[tuple[float, float]] = [(0.0, 0.0)]
    for ring in range(1, n_rings + 1):
        r = ring * pupil_radius / n_rings
        n_pts = 6 * ring
        two_pi_over_n = 2.0 * math.pi / n_pts
        for k in range(n_pts):
            theta = k * two_pi_over_n
            pts.append((r * math.cos(theta), r * math.sin(theta)))
    return np.array(pts, dtype=float)


def linear_pupil(n_samples: int, pupil_radius: float, axis: str = "y") -> np.ndarray:
    """1D pupil sweep along one meridian. Returns (n_samples, 2) xy."""
    if n_samples < 2:
        raise ValueError("n_samples must be >= 2")
    if pupil_radius <= 0.0:
        raise ValueError("pupil_radius must be > 0")
    coords = np.linspace(-pupil_radius, pupil_radius, n_samples)
    z = np.zeros_like(coords)
    if axis == "y":
        return np.column_stack([z, coords])
    if axis == "x":
        return np.column_stack([coords, z])
    raise ValueError(f"axis must be 'x' or 'y', got {axis!r}")


def pupil_radius_of(pre: Prescription) -> float:
    """Stop semi-diameter as the entrance pupil radius. Exact only when stop
    is surface 1 — deeper stops need pupil imaging (not implemented)."""
    return pre.surfaces[pre.stop - 1].semi_diameter


def launch_parallel(
    pupil_xy: np.ndarray,
    field_angle_rad: float = 0.0,
    field_axis: str = "y",
    *,
    pupil_z: float = 0.0,
    z_start: float = -20.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Collimated bundle from infinity: rays cross z=pupil_z at the given (x_p, y_p).

    z_start is the launch plane (must be in front of the first surface).
    """
    if field_axis == "y":
        dx, dy = 0.0, math.sin(field_angle_rad)
    elif field_axis == "x":
        dx, dy = math.sin(field_angle_rad), 0.0
    else:
        raise ValueError(f"field_axis must be 'x' or 'y', got {field_axis!r}")
    dz2 = 1.0 - dx * dx - dy * dy
    if dz2 <= 0.0:
        raise ValueError("|sin(field_angle)| >= 1; field angle too large")
    dz = math.sqrt(dz2)

    pupil_xy = np.asarray(pupil_xy, dtype=float)
    if pupil_xy.ndim != 2 or pupil_xy.shape[1] != 2:
        raise ValueError(f"pupil_xy must be shape (N, 2), got {pupil_xy.shape}")
    n = pupil_xy.shape[0]

    # Back-translate the launch position so the ray crosses pupil_z at (x_p, y_p).
    dt = (z_start - pupil_z) / dz
    pos = np.empty((n, 3), dtype=float)
    pos[:, 0] = pupil_xy[:, 0] + dt * dx
    pos[:, 1] = pupil_xy[:, 1] + dt * dy
    pos[:, 2] = z_start

    dirs = np.broadcast_to(np.array([dx, dy, dz]), (n, 3)).copy()
    return pos, dirs
