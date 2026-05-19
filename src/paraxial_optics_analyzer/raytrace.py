"""non-paraxial sequential ray trace through centered spherical surfaces."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from paraxial_optics_analyzer.prescription import Prescription


N_OBJ = 1.0


class TraceError(RuntimeError):
    pass


@dataclass(frozen=True)
class TraceResult:
    image_point: np.ndarray
    image_direction: np.ndarray
    hits: list[tuple[float, np.ndarray]]


def intersect_sphere(P: np.ndarray, d: np.ndarray, vertex_z: float, R: float) -> float:
    if not math.isfinite(R):
        if d[2] == 0.0:
            raise TraceError("ray parallel to plane surface")
        return (vertex_z - P[2]) / d[2]

    C = np.array([0.0, 0.0, vertex_z + R])
    PC = P - C
    b = float(np.dot(PC, d))
    c = float(np.dot(PC, PC) - R * R)
    disc = b * b - c
    if disc < 0.0:
        raise TraceError("ray misses surface (no real intersection)")
    s = math.sqrt(disc)
    t1, t2 = -b - s, -b + s
    z1 = P[2] + t1 * d[2]
    z2 = P[2] + t2 * d[2]
    return t1 if abs(z1 - vertex_z) < abs(z2 - vertex_z) else t2


def surface_normal(P_on_surface: np.ndarray, vertex_z: float, R: float) -> np.ndarray:
    if not math.isfinite(R):
        return np.array([0.0, 0.0, -1.0])
    C = np.array([0.0, 0.0, vertex_z + R])
    n_vec = (P_on_surface - C) / R
    norm = float(np.linalg.norm(n_vec))
    if norm == 0.0:
        raise TraceError("degenerate surface normal")
    return n_vec / norm


def refract(d: np.ndarray, n_hat: np.ndarray, n_before: float, n_after: float) -> np.ndarray:
    """Vector Snell"""
    mu = n_before / n_after
    cos_i = -float(np.dot(n_hat, d))
    sin2_t = mu * mu * (1.0 - cos_i * cos_i)
    if sin2_t > 1.0:
        raise TraceError("total internal reflection")
    cos_t = math.sqrt(1.0 - sin2_t)
    out = mu * d + (mu * cos_i - cos_t) * n_hat
    return out / float(np.linalg.norm(out))


def surface_vertex_z(pre: Prescription) -> np.ndarray:
    z = np.zeros(pre.n_surfaces)
    for i in range(1, pre.n_surfaces):
        z[i] = z[i - 1] + pre.surfaces[i - 1].thickness
    return z


def image_plane_z(pre: Prescription) -> float:
    return float(sum(s.thickness for s in pre.surfaces))


def trace_system(
    position: np.ndarray,
    direction: np.ndarray,
    pre: Prescription,
    *,
    image_plane_offset: float = 0.0,
) -> TraceResult:
    """Trace one ray surface-by-surface, finishing on the image plane"""
    P = np.array(position, dtype=float).reshape(3)
    d = np.array(direction, dtype=float).reshape(3)
    dn = float(np.linalg.norm(d))
    if dn == 0.0:
        raise TraceError("zero-length direction vector")
    d = d / dn

    n_before = N_OBJ
    z = 0.0
    hits: list[tuple[float, np.ndarray]] = []

    for i, surf in enumerate(pre.surfaces):
        t = intersect_sphere(P, d, z, surf.radius)
        if t < 0.0:
            raise TraceError(f"ray must travel backwards to hit surface {i + 1}")
        P = P + t * d
        hits.append((z, P.copy()))
        n_hat = surface_normal(P, z, surf.radius)
        d = refract(d, n_hat, n_before, surf.n)
        n_before = surf.n
        z += surf.thickness

    z_img = z + image_plane_offset
    if d[2] == 0.0:
        raise TraceError("output ray parallel to image plane")
    t_img = (z_img - P[2]) / d[2]
    return TraceResult(image_point=P + t_img * d, image_direction=d, hits=hits)
