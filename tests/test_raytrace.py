"""Tests for vector ray tracing primitives and full-system trace

Validation strategy:
1. Hand-checked single-surface intersections and Snell refraction
2. Snell's law in scalar form (n1 sin θ1 = n2 sin θ2) holds
3. Reversibility: refraction through a surface and back recovers the original ray
4. In the small-aperture limit, the full trace converges to the paraxial
  prediction
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from paraxial_optics_analyzer.paraxial import trace_paraxial
from paraxial_optics_analyzer.prescription import (
    ObjectSpec,
    Prescription,
    Surface,
)
from paraxial_optics_analyzer.raytrace import (
    TraceError,
    intersect_sphere,
    refract,
    surface_normal,
    trace_system,
)


class TestIntersectSphere:
    def test_axial_ray_hits_vertex(self):
        # Ray on the optical axis hits a convex surface at its vertex
        P = np.array([0.0, 0.0, -10.0])
        d = np.array([0.0, 0.0, 1.0])
        t = intersect_sphere(P, d, vertex_z=0.0, R=50.0)
        assert t == pytest.approx(10.0, rel=1e-12)

    def test_axial_ray_concave_surface(self):
        P = np.array([0.0, 0.0, -5.0])
        d = np.array([0.0, 0.0, 1.0])
        t = intersect_sphere(P, d, vertex_z=0.0, R=-30.0)
        assert t == pytest.approx(5.0, rel=1e-12)

    def test_plane_intersection(self):
        P = np.array([1.0, 2.0, -5.0])
        d = np.array([0.0, 0.0, 1.0])
        t = intersect_sphere(P, d, vertex_z=3.0, R=math.inf)
        assert t == pytest.approx(8.0, rel=1e-12)

    def test_offaxis_sphere_intersection_satisfies_sphere_equation(self):
        P = np.array([0.0, 5.0, -20.0])
        d = np.array([0.0, 0.0, 1.0])
        R = 50.0
        vertex_z = 0.0
        t = intersect_sphere(P, d, vertex_z, R)
        Q = P + t * d
        C = np.array([0.0, 0.0, vertex_z + R])
        assert float(np.linalg.norm(Q - C)) == pytest.approx(abs(R), rel=1e-12)

    def test_missed_surface_raises(self):
        # Ray offset far from axis can't reach a small sphere
        P = np.array([0.0, 1000.0, -10.0])
        d = np.array([0.0, 0.0, 1.0])
        with pytest.raises(TraceError, match="misses surface"):
            intersect_sphere(P, d, vertex_z=0.0, R=10.0)


class TestSurfaceNormal:
    def test_vertex_normal_is_minus_z(self):
        # At the vertex of any sphere (or plane), the outward normal points in -z
        n = surface_normal(np.array([0.0, 0.0, 0.0]), vertex_z=0.0, R=50.0)
        assert n == pytest.approx(np.array([0.0, 0.0, -1.0]), abs=1e-12)

    def test_vertex_normal_concave(self):
        n = surface_normal(np.array([0.0, 0.0, 0.0]), vertex_z=0.0, R=-50.0)
        assert n == pytest.approx(np.array([0.0, 0.0, -1.0]), abs=1e-12)

    def test_plane_normal(self):
        n = surface_normal(np.array([3.0, 4.0, 5.0]), vertex_z=5.0, R=math.inf)
        assert n == pytest.approx(np.array([0.0, 0.0, -1.0]), abs=1e-12)


class TestRefract:
    def test_normal_incidence_passes_through(self):
        d_in = np.array([0.0, 0.0, 1.0])
        n_hat = np.array([0.0, 0.0, -1.0])
        d_out = refract(d_in, n_hat, n_before=1.0, n_after=1.5)
        assert d_out == pytest.approx(np.array([0.0, 0.0, 1.0]), abs=1e-12)

    def test_obeys_snells_law(self):
        # Ray at 30° from the optical axis enters a flat interface n1=1 -> n2=1.5
        theta_i = math.radians(30.0)
        d_in = np.array([math.sin(theta_i), 0.0, math.cos(theta_i)])
        n_hat = np.array([0.0, 0.0, -1.0])
        d_out = refract(d_in, n_hat, n_before=1.0, n_after=1.5)
        assert float(np.linalg.norm(d_out)) == pytest.approx(1.0, abs=1e-12)
        # angle measured from -n_hat
        cos_t = -float(np.dot(n_hat, d_out))
        sin_t = math.sqrt(1.0 - cos_t * cos_t)
        assert 1.0 * math.sin(theta_i) == pytest.approx(1.5 * sin_t, rel=1e-12)

    def test_time_reversal(self):
        theta_i = math.radians(25.0)
        d_in = np.array([math.sin(theta_i), 0.1, math.cos(theta_i)])
        d_in = d_in / float(np.linalg.norm(d_in))
        n_hat = np.array([0.0, 0.0, -1.0])
        d_after = refract(d_in, n_hat, 1.0, 1.6)
        d_back = refract(-d_after, -n_hat, 1.6, 1.0)
        assert d_back == pytest.approx(-d_in, abs=1e-12)

    def test_total_internal_reflection_raises(self):
        theta_i = math.radians(70.0)
        d_in = np.array([math.sin(theta_i), 0.0, math.cos(theta_i)])
        # Normal points into the (denser) incident medium.
        n_hat = np.array([0.0, 0.0, -1.0])
        with pytest.raises(TraceError, match="total internal reflection"):
            refract(d_in, n_hat, n_before=1.6, n_after=1.0)




def _singlet(R1: float, R2: float, n: float, thickness: float, image_distance: float,
*, semi_d: float = 12.0) -> Prescription:
    return Prescription(
        name="singlet",
        wavelength_um=0.5876,
        units="mm",
        obj=ObjectSpec(distance=math.inf),
        surfaces=(
            Surface(radius=R1, thickness=thickness, n=n, semi_diameter=semi_d),
            Surface(radius=R2, thickness=image_distance, n=1.0, semi_diameter=semi_d),
        ),
        stop=1,
    )


class TestFullTrace:
    def test_axial_ray_through_singlet(self):
        """A ray on the optical axis must come out on the optical axis."""
        pre = _singlet(50.0, math.inf, 1.5168, 5.0, image_distance=91.75)
        P0 = np.array([0.0, 0.0, 0.0])
        d0 = np.array([0.0, 0.0, 1.0])
        r = trace_system(P0, d0, pre)
        assert float(np.linalg.norm(r.image_point[:2])) < 1e-12

    def test_real_trace_converges_to_paraxial(self):
        pre = _singlet(50.0, math.inf, 1.5168, 5.0, image_distance=91.75)
        paraxial = trace_paraxial(pre)

        from paraxial_optics_analyzer.raytrace import image_plane_z
        nominal_image_z = image_plane_z(pre)
        last_surface_z = nominal_image_z - pre.surfaces[-1].thickness
        defocus = (last_surface_z + paraxial.bfl) - nominal_image_z

        y_in = 1e-3
        P0 = np.array([0.0, y_in, 0.0])
        d0 = np.array([0.0, 0.0, 1.0])
        r = trace_system(P0, d0, pre, image_plane_offset=defocus)
        # Transverse miss at the paraxial focus
        assert abs(r.image_point[1]) < 1e-9
        assert abs(r.image_point[0]) < 1e-12

    @pytest.mark.parametrize("y_in", [1e-2, 1e-3, 1e-4])
    def test_real_trace_cubic_convergence(self, y_in):
        """The miss at paraxial focus should scale roughly like y to the third power (spherical aberration)."""
        pre = _singlet(50.0, math.inf, 1.5168, 5.0, image_distance=91.75)
        paraxial = trace_paraxial(pre)
        from paraxial_optics_analyzer.raytrace import image_plane_z
        nominal_image_z = image_plane_z(pre)
        last_surface_z = nominal_image_z - pre.surfaces[-1].thickness
        defocus = (last_surface_z + paraxial.bfl) - nominal_image_z

        P0 = np.array([0.0, y_in, 0.0])
        d0 = np.array([0.0, 0.0, 1.0])
        r = trace_system(P0, d0, pre, image_plane_offset=defocus)
        # Miss must be small AND scale at worst as y_in to the third power times a constant.
        assert abs(r.image_point[1]) < 100.0 * y_in ** 3

    def test_marginal_ray_misses_paraxial_focus_by_more_than_axial(self):
        """A real marginal ray crosses the axis short of the paraxial focus"""
        pre = _singlet(50.0, math.inf, 1.5168, 5.0, image_distance=91.75)
        P0 = np.array([0.0, 10.0, 0.0])
        d0 = np.array([0.0, 0.0, 1.0])
        paraxial = trace_paraxial(pre)
        from paraxial_optics_analyzer.raytrace import image_plane_z
        nominal_image_z = image_plane_z(pre)
        last_surface_z = nominal_image_z - pre.surfaces[-1].thickness
        defocus = (last_surface_z + paraxial.bfl) - nominal_image_z
        r = trace_system(P0, d0, pre, image_plane_offset=defocus)
        # The marginal-ray transverse aberration must be much larger than the near-axis one.
        assert abs(r.image_point[1]) > 1e-3
