"""Tests for spot diagram, ray fan, and best-focus search"""
from __future__ import annotations

import math

import numpy as np
import pytest

from paraxial_optics_analyzer.analysis import (
    find_best_focus,
    ray_fan,
    spot_diagram,
)
from paraxial_optics_analyzer.io import load_prescription
from paraxial_optics_analyzer.paraxial import trace_paraxial
from paraxial_optics_analyzer.prescription import (
    ObjectSpec,
    Prescription,
    Surface,
)
from paraxial_optics_analyzer.raytrace import image_plane_z


def _singlet(R1=50.0, R2=math.inf, n=1.5168, thickness=5.0,
             image_distance=91.75, *, semi_d=12.0) -> Prescription:
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


def _paraxial_focus_offset(prescription: Prescription) -> float:
    """Image-plane offset that places the focus exactly at the paraxial BFL."""
    paraxial = trace_paraxial(prescription)
    nominal_image_z = image_plane_z(prescription)
    last_surface_z = nominal_image_z - prescription.surfaces[-1].thickness
    return (last_surface_z + paraxial.bfl) - nominal_image_z


class TestSpotDiagram:
    def test_axial_spot_has_symmetric_centroid(self):
        """For an on-axis field, the centroid must lie on the optical axis."""
        pre = _singlet()
        sd = spot_diagram(pre)
        assert abs(sd.centroid[0]) < 1e-12
        assert abs(sd.centroid[1]) < 1e-12

    def test_spot_rms_smaller_at_paraxial_focus_than_far_away(self):
        pre = _singlet()
        offset = _paraxial_focus_offset(pre)
        sd_near = spot_diagram(pre, image_plane_offset=offset)
        sd_far = spot_diagram(pre, image_plane_offset=offset + 5.0)
        assert sd_near.rms < sd_far.rms

    def test_no_rays_failed_for_modest_aperture(self):
        sd = spot_diagram(_singlet(semi_d=8.0))
        assert sd.n_failed == 0


class TestRayFan:
    def test_axial_tangential_fan_is_odd(self):
        pre = _singlet()
        offset = _paraxial_focus_offset(pre)
        fan = ray_fan(pre, image_plane_offset=offset)
        # The vector should equal -reversed-self to high precision.
        tv = fan.transverse
        rev = tv[::-1]
        assert tv == pytest.approx(-rev, abs=1e-10)

    def test_axial_chief_ray_lands_on_axis(self):
        pre = _singlet()
        fan = ray_fan(pre)
        assert float(np.linalg.norm(fan.chief_image_xy)) < 1e-12

    def test_sagittal_and_tangential_match_on_axis(self):
        pre = _singlet()
        offset = _paraxial_focus_offset(pre)
        tan = ray_fan(pre, axis="tangential", image_plane_offset=offset)
        sag = ray_fan(pre, axis="sagittal", image_plane_offset=offset)
        # Absolute values must match (sign convention may differ by orientation)
        assert np.abs(tan.transverse) == pytest.approx(np.abs(sag.transverse), abs=1e-12)


class TestBestFocus:
    def test_finds_focus_better_than_nominal(self):
        pre = _singlet()
        result = find_best_focus(pre, n_rings=4)
        assert result.rms_at_best <= result.rms_at_nominal

    def test_best_focus_near_paraxial_for_modest_aperture(self):
        pre = _singlet(semi_d=2.0)  # very narrow pupil
        para_offset = _paraxial_focus_offset(pre)
        result = find_best_focus(pre, n_rings=4, search_range=(para_offset - 1.0, para_offset + 1.0))
        # Best focus should be within a small distance of paraxial focus for a thin pupil
        assert abs(result.image_plane_offset - para_offset) < 0.1

    def test_default_search_does_not_miss_cooke_paraxial_focus(self):
        pre = load_prescription("examples/cooke_triplet.yaml")
        para_offset = _paraxial_focus_offset(pre)
        paraxial = spot_diagram(pre, n_rings=4, image_plane_offset=para_offset)

        result = find_best_focus(pre, n_rings=4)

        assert result.rms_at_best <= paraxial.rms
        assert result.rms_at_best <= result.rms_at_nominal
        assert abs(result.image_plane_offset - para_offset) < 0.75


class TestSphericalAberrationConvergence:
    @staticmethod
    def _equiconvex_singlet(semi_d: float = 12.5) -> Prescription:
        n = 1.5168
        f_target = 100.0
        R = 2.0 * (n - 1.0) * f_target  # ≈ 103.36 mm
        return Prescription(
            name="equiconvex-f100",
            wavelength_um=0.5876,
            units="mm",
            obj=ObjectSpec(distance=math.inf),
            surfaces=(
                Surface(radius=R,  thickness=4.0,   n=n,   semi_diameter=semi_d),
                Surface(radius=-R, thickness=200.0, n=1.0, semi_diameter=semi_d),
            ),
            stop=1,
        )

    @staticmethod
    def _marginal_ta_at_paraxial(pre: Prescription) -> tuple[float, float, float]:
        from paraxial_optics_analyzer.paraxial import trace_paraxial
        from paraxial_optics_analyzer.raytrace import image_plane_z, trace_system
        para = trace_paraxial(pre)
        z_nominal = image_plane_z(pre)
        z_last = z_nominal - pre.surfaces[-1].thickness
        paraxial_offset = (z_last + para.bfl) - z_nominal
        h = pre.surfaces[0].semi_diameter
        r = trace_system(
            np.array([0.0, h, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            pre,
            image_plane_offset=paraxial_offset,
        )
        ta_max = float(r.image_point[1])
        u_marginal = h / abs(para.efl)
        return ta_max, paraxial_offset, u_marginal

    def test_paraxial_rms_matches_half_ta_max(self):
        pre = self._equiconvex_singlet()
        ta_max, paraxial_offset, _ = self._marginal_ta_at_paraxial(pre)
        sd = spot_diagram(pre, n_rings=24, image_plane_offset=paraxial_offset)
        expected = abs(ta_max) / 2.0
        # 7% tolerance: hexapolar sampling at N=24 + higher-order SA terms
        assert sd.rms == pytest.approx(expected, rel=0.07), (
            f"expected RMS ≈ |TA_max|/2 = {expected:.5f}, got {sd.rms:.5f}"
        )

    def _focused_search(self, pre, paraxial_offset, width=5.0):
        return find_best_focus(
            pre, n_rings=24, tol=1e-7,
            search_range=(paraxial_offset - width, paraxial_offset + width),
        )

    def test_best_focus_shift_matches_two_thirds_lsa(self):
        pre = self._equiconvex_singlet()
        ta_max, paraxial_offset, u_marginal = self._marginal_ta_at_paraxial(pre)
        # For a converging lens with under-corrected SA (marginal focus *before*
        # paraxial), best focus sits 2/3 of the way from paraxial toward marginal.
        lsa_magnitude = abs(ta_max) / u_marginal
        expected_shift = -(2.0 / 3.0) * lsa_magnitude
        fb = self._focused_search(pre, paraxial_offset)
        actual_shift = fb.image_plane_offset - paraxial_offset
        # 12% tolerance: hexapolar + higher-order SA + golden-section tol
        assert actual_shift == pytest.approx(expected_shift, rel=0.12), (
            f"expected best-focus shift from paraxial ~ {expected_shift:.4f}, "
            f"got {actual_shift:.4f}"
        )

    def test_best_focus_rms_matches_sixth_of_ta_max(self):
        pre = self._equiconvex_singlet()
        ta_max, paraxial_offset, _ = self._marginal_ta_at_paraxial(pre)
        fb = self._focused_search(pre, paraxial_offset)
        expected = abs(ta_max) / 6.0
        assert fb.rms_at_best == pytest.approx(expected, rel=0.10), (
            f"expected RMS at best focus ~ |TA_max|/6 = {expected:.5f}, "
            f"got {fb.rms_at_best:.5f}"
        )

    def test_low_ring_count_overestimates_rms(self):
        pre = self._equiconvex_singlet()
        _, paraxial_offset, _ = self._marginal_ta_at_paraxial(pre)
        rms_coarse = spot_diagram(pre, n_rings=4, image_plane_offset=paraxial_offset).rms
        rms_fine   = spot_diagram(pre, n_rings=24, image_plane_offset=paraxial_offset).rms
        assert rms_coarse > 1.25 * rms_fine
