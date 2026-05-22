"""Tests for the YAML loader."""
from __future__ import annotations

import math
import textwrap
from pathlib import Path

import pytest

from paraxial_optics_analyzer.io import load_prescription, prescription_from_dict
from paraxial_optics_analyzer.prescription import PrescriptionError

SINGLET_YAML = textwrap.dedent("""\
    name: Singlet BK7
    wavelength_um: 0.5876
    units: mm
    object:
      distance: .inf
    surfaces:
      - { radius:  50.0, thickness:  5.0,  n: 1.5168, semi_diameter: 15.0 }
      - { radius:  .inf, thickness: 91.75, n: 1.0,    semi_diameter: 15.0 }
    stop: 1
""")


def test_load_singlet_from_file(tmp_path: Path):
    p = tmp_path / "singlet.yaml"
    p.write_text(SINGLET_YAML, encoding="utf-8")

    pre = load_prescription(p)

    assert pre.name == "Singlet BK7"
    assert pre.wavelength_um == pytest.approx(0.5876)
    assert pre.units == "mm"
    assert math.isinf(pre.obj.distance)
    assert pre.obj.height == 0.0
    assert pre.n_surfaces == 2
    assert pre.surfaces[0].radius == 50.0
    assert pre.surfaces[0].n == pytest.approx(1.5168)
    assert math.isinf(pre.surfaces[1].radius)
    assert pre.stop == 1


def test_load_real_cooke_example():
    here = Path(__file__).resolve().parents[1] / "examples" / "cooke_triplet.yaml"
    pre = load_prescription(here)
    assert pre.n_surfaces == 6
    assert math.isinf(pre.obj.distance)
    assert pre.stop == 3


def test_missing_top_level_keys():
    with pytest.raises(PrescriptionError, match="missing required keys"):
        prescription_from_dict({"name": "x", "wavelength_um": 0.5, "units": "mm"})


def test_missing_surface_keys():
    raw = {
        "name": "x",
        "wavelength_um": 0.5876,
        "units": "mm",
        "object": {"distance": 100.0},
        "surfaces": [{"radius": 1.0, "thickness": 1.0}],
    }
    with pytest.raises(PrescriptionError, match="surface 1 is missing"):
        prescription_from_dict(raw)


def test_inf_as_string_accepted():
    raw = {
        "name": "x",
        "wavelength_um": 0.5876,
        "units": "mm",
        "object": {"distance": ".inf"},
        "surfaces": [
            {"radius": "inf", "thickness": 1.0, "n": 1.5, "semi_diameter": 5.0},
        ],
    }
    pre = prescription_from_dict(raw)
    assert math.isinf(pre.obj.distance)
    assert math.isinf(pre.surfaces[0].radius)


def test_bool_rejected():
    raw = {
        "name": "x",
        "wavelength_um": 0.5876,
        "units": "mm",
        "object": {"distance": True},
        "surfaces": [{"radius": 1.0, "thickness": 1.0, "n": 1.5, "semi_diameter": 5.0}],
    }
    with pytest.raises(PrescriptionError, match="got bool"):
        prescription_from_dict(raw)


def test_top_level_not_mapping(tmp_path: Path):
    p = tmp_path / "bad.yaml"
    p.write_text("- this is a list\n- not a mapping\n", encoding="utf-8")
    with pytest.raises(PrescriptionError, match="mapping"):
        load_prescription(p)


def test_invalid_yaml_syntax(tmp_path: Path):
    p = tmp_path / "bad.yaml"
    p.write_text("name: [unterminated\n", encoding="utf-8")
    with pytest.raises(PrescriptionError, match="failed to parse YAML"):
        load_prescription(p)


def test_surfaces_must_be_nonempty_list():
    raw = {
        "name": "x",
        "wavelength_um": 0.5876,
        "units": "mm",
        "object": {"distance": 100.0},
        "surfaces": [],
    }
    with pytest.raises(PrescriptionError, match="non-empty list"):
        prescription_from_dict(raw)


def test_surface_entry_must_be_mapping():
    raw = {
        "name": "x",
        "wavelength_um": 0.5876,
        "units": "mm",
        "object": {"distance": 100.0},
        "surfaces": ["not-a-dict"],
    }
    with pytest.raises(PrescriptionError, match="must be a mapping"):
        prescription_from_dict(raw)


def test_string_for_number_field_rejected():
    raw = {
        "name": "x",
        "wavelength_um": 0.5876,
        "units": "mm",
        "object": {"distance": 100.0},
        "surfaces": [{"radius": "fifty", "thickness": 1.0, "n": 1.5, "semi_diameter": 5.0}],
    }
    with pytest.raises(PrescriptionError, match="got string 'fifty'"):
        prescription_from_dict(raw)
