"""Unit tests for the prescription data model and its validation."""
from __future__ import annotations

import dataclasses
import math

import pytest

from paraxial_optics_analyzer.prescription import (
    ObjectSpec,
    Prescription,
    PrescriptionError,
    Surface,
)


def _ok_surface(**overrides):
    base = dict(radius=50.0, thickness=5.0, n=1.5168, semi_diameter=15.0)
    base.update(overrides)
    return Surface(**base)


def _ok_prescription(**overrides):
    base = dict(
        name="test",
        wavelength_um=0.5876,
        units="mm",
        obj=ObjectSpec(distance=math.inf),
        surfaces=(_ok_surface(),),
        stop=1,
    )
    base.update(overrides)
    return Prescription(**base)


class TestSurface:
    def test_ok(self):
        s = _ok_surface()
        assert s.radius == 50.0
        assert s.n == pytest.approx(1.5168)

    def test_zero_radius_rejected(self):
        with pytest.raises(PrescriptionError, match="radius cannot be 0"):
            _ok_surface(radius=0.0)

    def test_nan_radius_rejected(self):
        with pytest.raises(PrescriptionError, match="NaN"):
            _ok_surface(radius=float("nan"))

    def test_inf_radius_ok(self):
        s = _ok_surface(radius=math.inf)
        assert math.isinf(s.radius)

    def test_negative_radius_ok(self):
        # Concave-toward-light surfaces have negative radius.
        s = _ok_surface(radius=-25.0)
        assert s.radius == -25.0

    def test_negative_thickness_rejected(self):
        with pytest.raises(PrescriptionError, match="thickness must be >= 0"):
            _ok_surface(thickness=-1.0)

    def test_zero_thickness_ok(self):
        # Cemented doublets need zero air-gap thickness between elements.
        _ok_surface(thickness=0.0)

    def test_infinite_thickness_rejected(self):
        with pytest.raises(PrescriptionError, match="thickness must be finite"):
            _ok_surface(thickness=math.inf)

    def test_index_below_one_rejected(self):
        with pytest.raises(PrescriptionError, match="refractive index"):
            _ok_surface(n=0.5)

    def test_index_air_ok(self):
        s = _ok_surface(n=1.0)
        assert s.n == 1.0

    def test_nonpositive_semi_diameter_rejected(self):
        with pytest.raises(PrescriptionError, match="semi_diameter"):
            _ok_surface(semi_diameter=0.0)
        with pytest.raises(PrescriptionError, match="semi_diameter"):
            _ok_surface(semi_diameter=-1.0)

    def test_frozen(self):
        s = _ok_surface()
        with pytest.raises(dataclasses.FrozenInstanceError):
            s.radius = 99.0  # type: ignore[misc]


class TestObjectSpec:
    def test_infinity_ok(self):
        obj = ObjectSpec(distance=math.inf)
        assert math.isinf(obj.distance)
        assert obj.height == 0.0

    def test_finite_with_height_ok(self):
        obj = ObjectSpec(distance=1000.0, height=10.0)
        assert obj.distance == 1000.0
        assert obj.height == 10.0

    def test_nonpositive_distance_rejected(self):
        with pytest.raises(PrescriptionError, match="object distance"):
            ObjectSpec(distance=0.0)
        with pytest.raises(PrescriptionError, match="object distance"):
            ObjectSpec(distance=-1.0)

    def test_negative_height_rejected(self):
        with pytest.raises(PrescriptionError, match="object height"):
            ObjectSpec(distance=1000.0, height=-1.0)


class TestPrescription:
    def test_ok(self):
        p = _ok_prescription()
        assert p.n_surfaces == 1
        assert p.stop == 1

    def test_empty_name_rejected(self):
        with pytest.raises(PrescriptionError, match="name"):
            _ok_prescription(name="")

    def test_invalid_units_rejected(self):
        with pytest.raises(PrescriptionError, match="units"):
            _ok_prescription(units="furlong")

    def test_nonpositive_wavelength_rejected(self):
        with pytest.raises(PrescriptionError, match="wavelength_um"):
            _ok_prescription(wavelength_um=0.0)

    def test_infinite_wavelength_rejected(self):
        with pytest.raises(PrescriptionError, match="wavelength_um"):
            _ok_prescription(wavelength_um=math.inf)

    def test_no_surfaces_rejected(self):
        with pytest.raises(PrescriptionError, match="at least one"):
            _ok_prescription(surfaces=())

    def test_stop_out_of_range_rejected(self):
        with pytest.raises(PrescriptionError, match="stop"):
            _ok_prescription(stop=2)  # only 1 surface in the fixture

    def test_stop_zero_rejected(self):
        with pytest.raises(PrescriptionError, match="stop"):
            _ok_prescription(stop=0)

    def test_multi_surface_stop_ok(self):
        s = _ok_surface()
        p = _ok_prescription(surfaces=(s, s, s), stop=2)
        assert p.n_surfaces == 3
        assert p.stop == 2
