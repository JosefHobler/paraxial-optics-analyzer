"""Tests for pupil and field sampling utilities."""
from __future__ import annotations

import math

import numpy as np
import pytest

from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.sampling import (
    hexapolar_pupil,
    launch_parallel,
    linear_pupil,
    pupil_radius_of,
)


class TestHexapolar:
    def test_centre_only(self):
        pts = hexapolar_pupil(0, 10.0)
        assert pts.shape == (1, 2)
        assert pts[0, 0] == 0.0 and pts[0, 1] == 0.0

    def test_count_formula(self):
        # 1 + 3*n*(n+1)
        for n in range(1, 6):
            pts = hexapolar_pupil(n, 5.0)
            assert pts.shape == (1 + 3 * n * (n + 1), 2)

    def test_outermost_radius(self):
        pts = hexapolar_pupil(4, 5.0)
        # Outer ring radius should equal pupil radius
        r = np.linalg.norm(pts, axis=1)
        assert r.max() == pytest.approx(5.0, rel=1e-12)
        assert r.min() == 0.0

    def test_radius_must_be_positive(self):
        with pytest.raises(ValueError, match="pupil_radius"):
            hexapolar_pupil(3, 0.0)


class TestLinearPupil:
    def test_y_axis(self):
        pts = linear_pupil(5, 1.0, "y")
        assert pts.shape == (5, 2)
        assert np.allclose(pts[:, 0], 0)
        assert np.allclose(pts[:, 1], [-1.0, -0.5, 0.0, 0.5, 1.0])

    def test_x_axis(self):
        pts = linear_pupil(5, 1.0, "x")
        assert np.allclose(pts[:, 1], 0)
        assert np.allclose(pts[:, 0], [-1.0, -0.5, 0.0, 0.5, 1.0])

    def test_invalid_axis(self):
        with pytest.raises(ValueError, match="axis must be"):
            linear_pupil(5, 1.0, "z")


class TestLaunchParallel:
    def test_axial_field_directions_are_plus_z(self):
        pupil = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        positions, directions = launch_parallel(pupil, field_angle_rad=0.0)
        assert np.allclose(directions, np.array([0.0, 0.0, 1.0]))
        # All launched at z_start
        assert np.allclose(positions[:, 2], -20.0)
        assert np.allclose(positions[:, 0:2], pupil)

    def test_offaxis_field_directions(self):
        theta = math.radians(5.0)
        positions, directions = launch_parallel(np.array([[0.0, 0.0]]), field_angle_rad=theta)
        assert directions[0, 0] == pytest.approx(0.0, abs=1e-12)
        assert directions[0, 1] == pytest.approx(math.sin(theta), abs=1e-12)
        assert directions[0, 2] == pytest.approx(math.cos(theta), abs=1e-12)

    def test_rays_cross_pupil_plane_at_target_coords(self):
        """The whole point of the back-translation: at z=pupil_z, the ray is at (x_p, y_p)"""
        theta = math.radians(10.0)
        pupil = np.array([[2.0, -3.0], [0.0, 4.0], [1.5, 1.5]])
        positions, directions = launch_parallel(pupil, field_angle_rad=theta, pupil_z=0.0)
        # Propagate each ray to z=0
        for i in range(len(pupil)):
            t = (0.0 - positions[i, 2]) / directions[i, 2]
            xy_at_pupil = positions[i, 0:2] + t * directions[i, 0:2]
            assert xy_at_pupil == pytest.approx(pupil[i], abs=1e-12)

    def test_field_angle_too_large(self):
        with pytest.raises(ValueError, match="field angle too large"):
            launch_parallel(np.array([[0.0, 0.0]]), field_angle_rad=math.pi / 2)


class TestEntrancePupil:
    def test_internal_stop_is_imaged_by_front_group(self):
        pre = load_prescription("examples/cooke_triplet.yaml")

        assert pupil_radius_of(pre) == pytest.approx(8.856807495055554, rel=1e-12)
