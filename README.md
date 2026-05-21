# Paraxial Optics Analyzer

Scriptable sequential ray tracing & image-quality analysis for centered spherical lens systems — a tiny open-source Zemax-style tool.

Two independent paraxial implementations (direct refraction-transfer trace + ABCD matrix) cross-validate against the thick-lens lensmaker equation. The real (non-paraxial) trace converges to the paraxial prediction in the small-aperture limit at ~1e-9.

## Install

```bash
pip install -e ".[dev]"
```

Or via the task runner:

```bash
make install
```

## Quick demo

```bash
make demo
```

Output:

```
>>> Running analysis on examples/singlet_bk7.yaml
Plano-convex singlet (BK7)
  EFL: 96.7492260062 mm
  BFL: 93.4528125041 mm
  image distance: 93.4528125041 mm
  f-number: f/3.22497
  nominal RMS spot: 0.101376806029 mm
  best-focus offset: -0.418181032053 mm
  best-focus RMS spot: 0.0866224230684 mm

>>> Built-in self-checks
Lensmaker validation: PASS, relative error < 1e-15
Paraxial-limit validation: PASS, max deviation 1.2e-13
Cooke triplet EFL check: PASS, relative error 3.9e-16
```

## Make targets

| target          | what it does                                                       |
|-----------------|--------------------------------------------------------------------|
| `make install`  | editable install with dev deps (pytest, ruff)                      |
| `make test`     | run the pytest suite (~3 s)                                        |
| `make demo`     | run the headline analysis + built-in self-checks                   |
| `make report`   | write a 4-panel PDF report (summary, spot, ray fan, focus search)  |
| `make lint`     | static checks with ruff                                            |
| `make validate` | only the physics self-checks                                       |
| `make clean`    | drop `__pycache__`, `*.egg-info`, `.pytest_cache`, `.ruff_cache`   |

Pick a different prescription with `EXAMPLE=`:

```bash
make report EXAMPLE=examples/cooke_triplet.yaml REPORT=cooke.pdf
```

## CLI

Three subcommands. Run `analyze --help` for the top-level overview, `analyze <command> --help` for command-specific options.

```bash
analyze info examples/singlet_bk7.yaml                   # numeric results
analyze info examples/cooke_triplet.yaml --field-angle-deg 5
analyze report examples/cooke_triplet.yaml -o cooke.pdf  # write PDF/PNG report
analyze validate                                         # built-in self-checks
analyze --version
```

Top-level help:

```
$ analyze --help
usage: analyze [-h] [--version] <command> ...

Paraxial optics analyzer — sequential ray tracing for centered lens systems.

options:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

commands:
  <command>
    info      print first-order properties and spot statistics
    report    write a 4-panel PDF/PNG report
    validate  run built-in physics self-checks

examples:
  analyze info examples/singlet_bk7.yaml
  analyze info examples/cooke_triplet.yaml --field-angle-deg 5
  analyze report examples/cooke_triplet.yaml -o cooke.pdf
  analyze validate
```

Exit codes: `0` ok, `2` bad input (prescription / args), `3` ray-trace failure (TIR, total vignetting, chief-ray failure). The CLI prints a one-line `error: ...` to stderr with a hint pointing at `--field-angle-deg` or aperture size — failures don't drop a traceback.

## What's in the box

- `paraxial.py` — direct trace, ABCD matrix, lensmaker equation
- `raytrace.py` — vector Snell's law, sequential trace through spheres + planes
- `analysis.py` — spot diagram, ray fan, golden-section best-focus search
- `sampling.py` — hexapolar pupil + collimated-bundle launches
- `prescription.py` / `io.py` — data model + YAML loader
- `report.py` — matplotlib PDF/PNG report
- `validate.py` — three end-to-end physics checks the user can run from the CLI

## License

MIT
