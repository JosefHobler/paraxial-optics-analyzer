# Paraxial Optics Analyzer

Scriptable sequential ray-tracing and image-quality analysis for centered spherical optical systems (a small open-source analogue of one Zemax-style workflow)

Small Python project that reads a YAML lens prescription, traces paraxial and real rays through centered spherical optics, and generates first-order image-quality numbers plus a visual report.

Built-in validation currently passes:

- lensmaker validation: relative error `< 1e-15`
- paraxial-limit validation: max deviation `1.2e-13`
- Cooke triplet EFL check: relative error `3.9e-16`

https://github.com/user-attachments/assets/995828a9-2453-4f28-967f-809c66763304

## Example Output

The input is a small lens prescription in YAML:

```yaml
name: Cooke triplet style example
wavelength_um: 0.5876
units: mm

object:
  distance: .inf
  height: 0.0

surfaces:
  - { radius: 26.0, thickness: 4.0, n: 1.6116, semi_diameter: 8.85 }
  - { radius: -48.0, thickness: 6.0, n: 1.0, semi_diameter: 8.25 }
  - { radius: -32.0, thickness: 2.0, n: 1.6040, semi_diameter: 6.0 }
  - { radius: 22.0, thickness: 8.0, n: 1.0, semi_diameter: 6.0 }
  - { radius: 45.0, thickness: 3.5, n: 1.5168, semi_diameter: 8.0 }
  - { radius: -65.0, thickness: 42.66, n: 1.0, semi_diameter: 8.0 }

stop: 3
```

Run the analysis:

```bash
analyze info examples/cooke_triplet.yaml
```

Current output:

```text
Cooke triplet style example
  EFL: 54.1813925245 mm
  BFL: 38.8823263606 mm
  image distance: 38.8823263606 mm
  f-number: f/3.05874

  nominal image plane z:  66.16 mm
  paraxial focus z:       60.0423 mm   (shift from nominal: -6.118 mm)
  best-focus z:           60.2614 mm   (shift from paraxial: +0.2191 mm)

  RMS spot at nominal:     0.705092 mm
  RMS spot at paraxial:    0.044284 mm
  RMS spot at best focus:  0.035734 mm
```

Generate a report:

```bash
analyze report examples/cooke_triplet.yaml -o cooke_report.png
```

The report contains a summary, spot diagram, tangential/sagittal ray fan, and best-focus RMS curve.

<img width="1100" height="850" alt="demo_report" src="https://github.com/user-attachments/assets/d0c3b1b9-d58e-4ee4-bcee-0db63e85317d" />

## Built-in Checks


The project is intentionally small enough that the numbers can be checked against independent paths.

Current built-in checks include:

- BK7 singlet paraxial EFL against the thick-lens lensmaker equation
- direct paraxial trace against an independent ABCD matrix calculation
- real vector-Snell ray trace converging to the paraxial result in the small-aperture limit
- Cooke-triplet internal-stop f-number regression, where the correct value depends on entrance-pupil diameter rather than raw stop diameter
- best-focus search regression for a case where the nominal image plane is far from the paraxial focus

Run:

```bash
analyze validate
```

Current output:

```text
Lensmaker validation: PASS, relative error < 1e-15
Paraxial-limit validation: PASS, max deviation 1.2e-13
Cooke triplet EFL check: PASS, relative error 3.9e-16
```

## Scope

This is not a replacement for Zemax or any professional optical-design tool. It deliberately excludes tolerancing, optimization, aspheres, tilts, decenters, diffraction, coating effects, polarization, and full glass catalogs.

Works:

- rotationally symmetric sequential systems
- spherical surfaces
- plane surfaces, written as `.inf` radius
- air object space, `n = 1.0`
- arbitrary refractive index per surface, supplied directly as `n`
- object at infinity
- finite object distance in the paraxial calculation
- axial and off-axis collimated bundles through `--field-angle-deg`
- paraxial EFL, BFL, image distance, and f-number
- ABCD matrix in reduced-angle convention
- thick-lens lensmaker helper for a singlet
- real non-paraxial trace through all surfaces
- vector Snell refraction
- TIR detection
- spot diagrams
- tangential and sagittal ray fans
- best-focus search by RMS spot radius
- PDF/PNG report using matplotlib
- YAML prescription files
- command line interface
- direct Python API
- CI with lint, tests, validation, wheel build, and wheel smoke test

Deliberate simplifications:

- no aspheres
- no tilts or decenters
- no mirrors
- no GRIN media
- no wavelength-dependent glass catalog; `wavelength_um` is metadata right now
- no dispersion/chromatic analysis
- no optimization
- no tolerancing
- no diffraction, MTF, PSF, OPD, or wavefront output
- no coatings, Fresnel losses, polarization, or ghost analysis
- no obscurations
- no vignetting model beyond rays failing when they miss geometry or hit TIR
- surface `semi_diameter` is validated, but the real tracer does not currently clip each surface by semi-diameter
- entrance pupil is paraxial and axial; good enough for current tests, not a full pupil solver for arbitrary fields
- finite object height exists in the schema, but the analysis tools mostly use field angle / collimated input

## Install

Python requirement: `>=3.10`.

Runtime dependencies:

- `numpy`
- `matplotlib`
- `pyyaml`

Dev dependencies:

- `pytest`
- `pytest-cov`
- `ruff`

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install -e ".[dev]"
```

Git Bash / Linux / macOS:

```bash
python -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e ".[dev]"
```

The installed console script is `analyze`.

```bash
analyze --version
```

If `analyze` is not on PATH, run the CLI through Python:

```bash
.venv/Scripts/python -m paraxial_optics_analyzer.cli --version
```

or on Unix-like shells:

```bash
.venv/bin/python -m paraxial_optics_analyzer.cli --version
```

## Quick Start

Analyze the bundled Cooke-triplet example:

```bash
analyze info examples/cooke_triplet.yaml
```

Generate a visual report:

```bash
analyze report examples/cooke_triplet.yaml -o cooke_report.png
```

Run built-in validation checks:

```bash
analyze validate
```

Run tests and lint locally:

```bash
pytest -q
ruff check src tests
```

## CLI

Top-level commands:

```bash
analyze info <prescription.yaml>
analyze report <prescription.yaml>
analyze validate
analyze --version
```

`info` prints first-order properties and spot statistics:

```bash
analyze info examples/cooke_triplet.yaml
analyze info examples/cooke_triplet.yaml --field-angle-deg 5
analyze info examples/cooke_triplet.yaml --rings 24
```

`report` writes a 4-panel PDF/PNG:

```bash
analyze report examples/cooke_triplet.yaml -o cooke.pdf
analyze report examples/cooke_triplet.yaml -o cooke.png
analyze report examples/cooke_triplet.yaml --field-angle-deg 5 --rings 12
```

If `-o/--output` is omitted, report output defaults to:

```text
<prescription_stem>_report.pdf
```

Shared analysis flags:

- `--field-angle-deg DEG`: field angle in degrees. Default `0.0`.
- `--rings N`: hexapolar pupil rings. CLI default `16`. Lower is faster; higher gives better RMS convergence.

Exit codes:

- `0`: ok
- `1`: validation command found failed checks
- `2`: bad prescription or bad CLI argument
- `3`: ray trace failure, usually TIR, total failed bundle, or impossible ray

User-facing errors are intended to be one-line `error: ...` messages, not Python tracebacks.

## YAML Prescription Format

Minimal singlet:

```yaml
name: Plano-convex singlet (BK7)
wavelength_um: 0.5876
units: mm

object:
  distance: .inf
  height: 0.0

surfaces:
  - { radius: 50.0, thickness: 5.0, n: 1.5168, semi_diameter: 15.0 }
  - { radius: .inf, thickness: 91.75, n: 1.0, semi_diameter: 15.0 }

stop: 1
```

Top-level required keys:

- `name`
- `wavelength_um`
- `units`
- `object`
- `surfaces`

Top-level optional key:

- `stop`, defaults to `1`

Accepted `units` values:

- `mm`
- `cm`
- `m`
- `inch`

`object` keys:

- `distance`: required. Must be positive or infinity.
- `height`: optional, defaults to `0.0`. Must be non-negative.

Each surface requires:

- `radius`: signed radius of curvature. Use `.inf` for a plane.
- `thickness`: axial distance to the next surface. Must be finite and `>= 0`.
- `n`: refractive index of the medium after the surface. Must be finite and `>= 1`.
- `semi_diameter`: clear-aperture half-width. Must be finite and `> 0`.

Bundled examples:

- `examples/singlet_bk7.yaml`
- `examples/cooke_triplet.yaml`

## Coordinate And Sign Conventions

- Optical axis is `+z`.
- Surface vertices are placed sequentially along `z`.
- First surface vertex is at `z = 0`.
- `surface_vertex_z(pre)` returns all surface vertex positions.
- `image_plane_z(pre)` is the sum of all surface thicknesses.
- Radius is signed.
- Positive radius means the center of curvature is downstream of the vertex.
- Plane surfaces are represented by infinite radius.
- BFL is measured from the last surface vertex to the rear focal point.
- EFL is positive for a converging system.
- Field angle is in the tangential `y-z` meridian by default.
- Sagittal fan samples along `x`.
- Object-space refractive index is fixed at `1.0`.
- ABCD matrices use the reduced-angle state `[y, n*u]`.

## Architecture

```text
lens.yaml
   |
   v
[parser / validation]
   |
   v
[Prescription dataclass]
   |
   +--> [paraxial core] ----> EFL / BFL / f-number / image distance
   |
   +--> [real-ray core] ----> spot diagram / ray fan / best focus
   |
   +--> [report writer] ---> PDF / PNG report
```

The computational core has no CLI or file-output dependency. YAML loading, command-line handling, and report writing are kept outside the lower-level physics routines.

## Tests

Run:

```bash
pytest -q
ruff check src tests
analyze validate
```

Current local count: 110 tests.

Test files:

- `test_prescription.py`: dataclass validation and frozen-ness
- `test_io.py`: YAML loader, required keys, type coercion, infinity parsing, bool rejection
- `test_paraxial.py`: lensmaker, ABCD/direct agreement, Gaussian imaging, f-number
- `test_raytrace.py`: sphere intersections, normals, vector Snell, TIR, full trace
- `test_sampling.py`: hexapolar counts, linear pupil, launch geometry, entrance pupil
- `test_analysis.py`: spots, fans, best focus, spherical aberration convergence
- `test_cli_report.py`: CLI info/report/error paths and report file writing
- `test_validate.py`: validation checks and validate CLI command

Some tests are physics tests, not just code-shape tests:

- BK7 singlet EFL against thick-lens lensmaker
- direct paraxial trace against ABCD matrix
- real trace convergence to paraxial focus in the small-aperture limit
- third-order spherical aberration relationships for an equiconvex singlet
- Cooke triplet EFL regression
- Cooke internal-stop entrance-pupil/f-number regression
- Cooke best-focus search regression

## CI

Workflow: `.github/workflows/ci.yml`.

Runs on:

- push to `main` or `master`
- pull request to `main` or `master`
- manual dispatch

Jobs:

- `ruff` on Python 3.12
- tests on Python 3.10, 3.11, 3.12
- `analyze validate` after tests
- build `sdist` and wheel
- install the wheel into a clean environment
- run `analyze --version`
- run `analyze validate`

The workflow cancels older in-progress runs on the same branch/PR.

## Repository Layout

```text
.
  .github/workflows/ci.yml
  docs/
    demo.tape
    demo_report.png
  examples/
    singlet_bk7.yaml
    cooke_triplet.yaml
  src/paraxial_optics_analyzer/
    __init__.py
    analysis.py
    cli.py
    io.py
    paraxial.py
    prescription.py
    raytrace.py
    report.py
    sampling.py
    validate.py
  tests/
    conftest.py
    test_analysis.py
    test_cli_report.py
    test_io.py
    test_paraxial.py
    test_prescription.py
    test_raytrace.py
    test_sampling.py
    test_validate.py
  Makefile
  pyproject.toml
  LICENSE
  README.md
```

Generated stuff such as `build/`, `dist/`, `*.egg-info`, `__pycache__/`, `.pytest_cache/`, and `.ruff_cache/` is not part of the source. The Makefile has a `clean` target for most of that.

## Bugs Caught During Validation

Internal stop f-number:

- Wrong behavior: `f/# = EFL / stop_diameter`
- Correct behavior: `f/# = EFL / entrance_pupil_diameter`
- This matters for `examples/cooke_triplet.yaml`, where the stop is at surface 3.

Best focus:

- Wrong behavior: search could miss paraxial focus when nominal image plane was far away.
- Correct behavior: default search is centered around paraxial focus, grid bracketed, then polished.

Test fixture drift:

- Several tests once loaded the singlet YAML while asserting Cooke behavior.
- If CI fails with an obvious title/surface-count mismatch, check fixture paths first.

