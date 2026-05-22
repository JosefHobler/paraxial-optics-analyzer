"""Tests for the paraxial trace and first-order properties

Validation strategy:
1. paraxial trace must agree with the thick-lens lensmaker equation to ~1e-12
2. direct refraction-transfer trace must agree with the system ABCD matrix to ~1e-12
3. image-distance for a finite conjugate must satisfy the Gaussian imaging
  equation 1/s + 1/s' = 1/f for a thin singlet
"""
from __future__ import annotations

import math

import pytest

from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.paraxial import (
    bfl_from_matrix,
    efl_from_matrix,
    lensmaker_thick,
    system_matrix,
    trace_paraxial,
)
from paraxial_optics_analyzer.prescription import (
    ObjectSpec,
    Prescription,
    Surface,
)


def _singlet(R1: float, R2: float, n: float, thickness: float, *, semi_d: float = 12.0,
             obj_distance: float = math.inf, stop: int = 1) -> Prescription:
    return Prescription(
        name="singlet",
        wavelength_um=0.5876,
        units="mm",
        obj=ObjectSpec(distance=obj_distance),
        surfaces=(
            Surface(radius=R1, thickness=thickness, n=n, semi_diameter=semi_d),
            Surface(radius=R2, thickness=100.0, n=1.0, semi_diameter=semi_d),
        ),
        stop=stop,
    )


class TestLensmakerAgreement:
    """Paraxial trace EFL must almost equal lensmaker EFL."""

    @pytest.mark.parametrize("R1, R2, n, t", [
        ( 50.0,  math.inf, 1.5168, 0.001),   # near-thin plano-convex
        ( 50.0,  math.inf, 1.5168, 5.0),     # thick plano-convex
        (100.0, -100.0,    1.5,    0.001),   # near-thin equi-convex
        (100.0, -100.0,    1.5,    8.0),     # thick equi-convex (non-trivial thickness term)
        ( 80.0,  200.0,    1.6,    3.0),     # meniscus
        (math.inf, -75.0,  1.7,    4.0),     # plano-convex with the curved side last
    ])
    def test_efl(self, R1, R2, n, t):
        f_paraxial = trace_paraxial(_singlet(R1, R2, n, t)).efl
        f_lensmaker = lensmaker_thick(R1, R2, n, t)
        assert f_paraxial == pytest.approx(f_lensmaker, rel=1e-12, abs=1e-12)


class TestTraceMatrixAgreement:
    """Direct trace and ABCD-matrix readouts must almost agree."""

    @pytest.mark.parametrize("R1, R2, n, t", [
        ( 50.0,  math.inf, 1.5168, 5.0),
        (100.0, -100.0,    1.5,    8.0),
        ( 80.0,  200.0,    1.6,    3.0),
    ])
    def test_efl_matches(self, R1, R2, n, t):
        pre = _singlet(R1, R2, n, t)
        efl_trace = trace_paraxial(pre).efl
        efl_matrix = efl_from_matrix(pre)
        assert efl_trace == pytest.approx(efl_matrix, rel=1e-12, abs=1e-12)

    @pytest.mark.parametrize("R1, R2, n, t", [
        ( 50.0,  math.inf, 1.5168, 5.0),
        (100.0, -100.0,    1.5,    8.0),
        ( 80.0,  200.0,    1.6,    3.0),
    ])
    def test_bfl_matches(self, R1, R2, n, t):
        pre = _singlet(R1, R2, n, t)
        bfl_trace = trace_paraxial(pre).bfl
        bfl_matrix = bfl_from_matrix(pre)
        assert bfl_trace == pytest.approx(bfl_matrix, rel=1e-12, abs=1e-12)


class TestImageDistance:
    def test_infinity_object_equals_bfl(self):
        r = trace_paraxial(_singlet(50.0, math.inf, 1.5168, 5.0))
        assert r.image_distance == pytest.approx(r.bfl, rel=1e-12, abs=1e-12)

    def test_gaussian_imaging_equation(self):
        """For finite object distance, 1/s + 1/s' = 1/f (Gaussian formula) at the principal planes.

        For a thin lens with object distance s in front, image distance s' behind,
        and focal length f, the Gaussian equation holds exactly. We approximate
        a thin lens by using a very small physical thickness.
        """
        thin = _singlet(50.0, math.inf, 1.5168, 1e-6, obj_distance=300.0)
        r = trace_paraxial(thin)
        s = 300.0
        sprime = r.image_distance
        f_implied = 1.0 / (1.0 / s + 1.0 / sprime)
        assert f_implied == pytest.approx(r.efl, rel=1e-6)


class TestFNumber:
    def test_f_number_for_stop_at_first_surface(self):
        # f = approximately 96.75 mm, semi-diameter 12 -> diameter 24 -> f/# = approximately 4.03
        pre = _singlet(50.0, math.inf, 1.5168, 5.0, semi_d=12.0, stop=1)
        r = trace_paraxial(pre)
        expected = r.efl / 24.0
        assert r.f_number == pytest.approx(expected, rel=1e-12)

    def test_f_number_for_internal_stop_uses_entrance_pupil(self):
        pre = load_prescription("examples/singlet_bk7.yaml")
        r = trace_paraxial(pre)
        entrance_pupil_radius = 8.856807495055554

        assert r.f_number == pytest.approx(r.efl / (2.0 * entrance_pupil_radius), rel=1e-12)


class TestSystemMatrix:
    def test_identity_for_zero_power_plate(self):
        """A parallel plate (R1=R2=inf) has refractive matrix = identity except for the slowed transfer."""
        pre = Prescription(
            name="plate",
            wavelength_um=0.5876,
            units="mm",
            obj=ObjectSpec(distance=math.inf),
            surfaces=(
                Surface(radius=math.inf, thickness=5.0,  n=1.5, semi_diameter=10.0),
                Surface(radius=math.inf, thickness=10.0, n=1.0, semi_diameter=10.0),
            ),
            stop=1,
        )
        A, B, C, D = system_matrix(pre)
        # For a parallel plate
        assert C == pytest.approx(0.0, abs=1e-14)
        assert A == pytest.approx(1.0, abs=1e-14)
        assert D == pytest.approx(1.0, abs=1e-14)
        assert B == pytest.approx(5.0 / 1.5, rel=1e-14)
