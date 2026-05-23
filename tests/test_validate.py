from __future__ import annotations

import pytest

from paraxial_optics_analyzer.cli import main
from paraxial_optics_analyzer.validate import (
    check_cooke_triplet,
    check_lensmaker,
    check_paraxial_limit,
    run_all,
)


class TestChecks:
    def test_lensmaker_passes(self):
        r = check_lensmaker()
        assert r.passed
        assert r.value <= r.threshold

    def test_paraxial_limit_passes(self):
        r = check_paraxial_limit()
        assert r.passed
        assert r.value <= r.threshold

    def test_cooke_triplet_passes(self):
        r = check_cooke_triplet()
        assert r.passed
        assert r.value <= r.threshold

    def test_run_all_returns_three_checks(self):
        results = run_all()
        assert len(results) == 3
        assert all(r.passed for r in results)


class TestFormatLine:
    def test_pass_with_zero(self):
        r = check_lensmaker()
        # The plano singlet gives exact-zero relative error; format should
        # render that as "< 1e-15", not the bare "0".
        if r.value == 0.0:
            assert "< 1e-15" in r.format_line()
        assert "PASS" in r.format_line()

    def test_fail_renders_fail_verdict(self):
        from paraxial_optics_analyzer.validate import CheckResult
        r = CheckResult(name="x", passed=False, metric="err", value=1.0, threshold=0.1)
        assert "FAIL" in r.format_line()


class TestCli:
    def test_validate_subcommand(self, capsys):
        rc = main(["validate"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "Lensmaker validation: PASS" in out
        assert "Paraxial-limit validation: PASS" in out
        assert "Cooke triplet EFL check: PASS" in out

    def test_no_subcommand_errors(self):
        with pytest.raises(SystemExit):
            main([])

    def test_version_flag(self, capsys):
        with pytest.raises(SystemExit) as exc:
            main(["--version"])
        assert exc.value.code == 0
        out = capsys.readouterr().out
        assert "analyze" in out
