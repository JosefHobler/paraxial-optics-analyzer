# One-command task runner.
#   make test       run the test suite
#   make demo       run the headline analysis on the BK7 singlet
#   make report     write a PDF report for the BK7 singlet
#   make lint       static checks (ruff)
#   make validate   built-in physics self-checks
#   make install    editable install with dev deps
#   make clean      drop build / cache artefacts

PYTHON  ?= python
EXAMPLE ?= examples/singlet_bk7.yaml
REPORT  ?= singlet_report.pdf

.PHONY: help install test demo report lint validate clean

help:
	@echo "Targets: install test demo report lint validate clean"

install:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	$(PYTHON) -m pytest -q

demo:
	@echo ">>> Running analysis on $(EXAMPLE)"
	$(PYTHON) -m paraxial_optics_analyzer.cli $(EXAMPLE) --no-report --rings 4
	@echo ""
	@echo ">>> Built-in self-checks"
	$(PYTHON) -m paraxial_optics_analyzer.cli --validate

report:
	$(PYTHON) -m paraxial_optics_analyzer.cli $(EXAMPLE) -o $(REPORT)

lint:
	$(PYTHON) -m ruff check src tests

validate:
	$(PYTHON) -m paraxial_optics_analyzer.cli --validate

clean:
	$(PYTHON) -c "import pathlib, shutil; [shutil.rmtree(p, ignore_errors=True) for p in pathlib.Path('.').rglob('__pycache__')]; [shutil.rmtree(p, ignore_errors=True) for p in pathlib.Path('.').rglob('*.egg-info')]; [shutil.rmtree(p, ignore_errors=True) for p in ['.pytest_cache', '.ruff_cache', 'build', 'dist']]"
