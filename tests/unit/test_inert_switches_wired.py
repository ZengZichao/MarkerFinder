"""Switches that existed but reached nothing.

Two CLI/config knobs were advertised ( /) yet had no
consumer in the code, which is the same defect class as baseline
("stored as an instance attribute and never read"). These tests give each
knob observable teeth plus a must-fail control, and pin the wiring so a future
refactor cannot silently disconnect it again.
"""

from __future__ import annotations

import argparse
import inspect
import logging
from types import SimpleNamespace

import pytest

from markerfinder.config import HGTConfig, ReportConfig
from markerfinder.modules.hgt_scan import run_threshold_scan
from markerfinder.models.marker import MarkerLevel
from markerfinder.utils.informative_sites import apply_pis_floor


def _ev(marker_id, risk, level=MarkerLevel.LEVEL_2):
    return SimpleNamespace(marker_id=marker_id, overall_risk=risk, level=level)


class TestScanStabilityMin:
    """--scan-stability-min must actually move the scanner's verdict."""

    def _bands(self):
        # Two bands whose selected marker sets overlap only halfway:
        # Jaccard = |{M1,M2} / {M1..M4}| = 2/4 = 0.5.
        return [
            _ev("M1", 0.10), _ev("M2", 0.10),
            _ev("M3", 0.50), _ev("M4", 0.50),
        ]

    def _warned(self, tmp_path, caplog, floor):
        with caplog.at_level(logging.WARNING):
            run_threshold_scan(
                self._bands(), str(tmp_path), "mf",
                bands=[(0.20, 0.60), (0.15, 0.40)], stability_min=floor,
            )
        return any("unstable across bands" in r.message for r in caplog.records)

    def test_low_floor_stays_silent(self, tmp_path, caplog):
        assert self._warned(tmp_path, caplog, 0.2) is False

    def test_high_floor_warns(self, tmp_path, caplog):
        # MUST-FAIL CONTROL: same data, stricter floor — the only difference is
        # The wired-through value, so a dead flag cannot pass this.
        assert self._warned(tmp_path, caplog, 0.9) is True

    def test_config_declares_the_key_with_the_documented_default(self):
        assert ReportConfig().scan_stability_min == 0.6

    def test_pipeline_forwards_the_configured_value(self):
        """Wiring lock: the scanner must be called with the config value."""
        import markerfinder.pipeline as pipeline_module

        src = inspect.getsource(pipeline_module)
        assert "stability_min=" in src, "scan floor not forwarded to scanner"
        assert "scan_stability_min" in src, "config key never read by the pipeline"


class TestPisFloor:
    """Too few parsimony-informative sites must not read as clean."""

    def _card_eval(self, marker_id="M1", level=MarkerLevel.LEVEL_1):
        return SimpleNamespace(
            marker_id=marker_id, level=level, notes="", decision_card={}
        )

    def test_floor_zero_is_a_noop(self):
        ev = self._card_eval()
        assert apply_pis_floor([ev], {"M1": 1}, 0) == []
        assert ev.notes == "" and ev.decision_card == {}

    def test_below_floor_marked_inconclusive(self):
        ev = self._card_eval()
        demoted = apply_pis_floor([ev], {"M1": 3}, 5)
        assert demoted == ["M1"]
        assert "inconclusive" in ev.notes
        assert ev.decision_card["pis_grade"] == "inconclusive"
        assert ev.decision_card["pis_floor"] == 5
        assert ev.decision_card["pis_measured"] == 3

    def test_at_or_above_floor_untouched(self):
        ev = self._card_eval()
        assert apply_pis_floor([ev], {"M1": 5}, 5) == []
        assert ev.notes == ""

    def test_unmeasured_pis_is_not_guessed(self):
        ev = self._card_eval()
        assert apply_pis_floor([ev], {}, 5) == []
        assert ev.decision_card == {}

    def test_unscreened_markers_are_not_regraded(self):
        ev = self._card_eval(level=MarkerLevel.UNKNOWN)
        assert apply_pis_floor([ev], {"M1": 0}, 5) == []

    def test_grading_is_never_altered(self):
        """The floor stamps inconclusive; it does not invent exclusion."""
        ev = self._card_eval()
        apply_pis_floor([ev], {"M1": 1}, 10)
        assert ev.level is MarkerLevel.LEVEL_1

    def test_hgt_config_default_keeps_the_screen_off(self):
        assert HGTConfig().min_informative_sites == 0

    def test_pipeline_consumes_the_configured_floor(self):
        import markerfinder.pipeline as pipeline_module

        src = inspect.getsource(pipeline_module)
        assert "apply_pis_floor" in src, "PIS floor never applied by the run"
        assert "min_informative_sites" in src, "floor config key unread"


class TestConsistencyModeWiring:
    """A labelled criterion must either run or refuse — never mislabel.

    ``screen`` had no caller while the gate let
    ``--hgt-mode consistency`` through, which would have emitted risk-mode
    grading under a consistency label. ``consistency`` is now actually computed
    in Phase 3; ``hybrid`` is the spec-fixed conjunction (both
    legs must pass) and is wired through the same screen block.
    """

    def test_consistency_criterion_is_called_by_the_run(self):
        import markerfinder.pipeline as pipeline_module

        src = inspect.getsource(pipeline_module)
        assert "consistency_results = screen(" in src, (
            "screen() is not called by the run: consistency mode would lie"
        )
        assert "run_consistency_stage(" in src, (
            "excluded-marker profile is never produced (the stage "
            "owns the write; the call must reach run_stage)"
        )
        assert "consistency_grades" in src, "grades never reach the products"

    def test_hybrid_runs_the_same_screen_and_records_the_conjunction(self):
        """Hybrid = risk AND consistency, both must pass.

        The old lock asserted hybrid *refuses*; the spec fixes the semantics
        ("两者都要过", applied 与 ``excluded`` 并列), so the lock flipped to
        assert the conjunction is real: the screen block admits both modes,
        the verdicts are annotated per marker, and no refusal branch remains.
        """
        import markerfinder.pipeline as pipeline_module

        src = inspect.getsource(pipeline_module)
        gate = src.index("guard(self.config.hgt_config.hgt_mode")
        block = src.index("if _hgt_mode in (\"consistency\", \"hybrid\"):", gate)
        assert src.index("_annotate_hybrid_verdicts(", block) > block, (
            "hybrid conjunction verdicts are never annotated onto the products"
        )
        assert "hybrid_passed" in src[block:], (
            "the consistency leg of the conjunction is not computed"
        )
        assert "UnsupportedCriterion" not in src[gate:block], (
            "a hybrid refusal branch survived the spec decision"
        )

    def test_risk_remains_the_untouched_default(self):
        # The default mode never auto-switches, so none of the above can
        # Change baseline behaviour.
        assert HGTConfig().hgt_mode == "risk"


class TestCogCategoryMap:
    """The loader exists, and its absence renders NA instead of a blank."""

    def test_loads_two_column_tsv_and_skips_noise(self, tmp_path):
        from markerfinder.modules.consistency_screen import load_cog_category_map

        f = tmp_path / "cog.tsv"
        f.write_text(
            "# comment line\n"
            "marker_id\tfunctional_category\n"
            "M1\tribosomal\n"
            "M2\ttranslation\n"
            "\n"
            "M3\t\n"
            "lonely_row\n",
            encoding="utf-8",
        )
        assert load_cog_category_map(str(f)) == {
            "M1": "ribosomal", "M2": "translation",
        }

    def test_missing_or_unset_path_yields_empty_mapping(self, tmp_path):
        from markerfinder.modules.consistency_screen import load_cog_category_map

        assert load_cog_category_map(None) == {}
        assert load_cog_category_map(str(tmp_path / "absent.tsv")) == {}

    def test_profile_product_carries_the_category(self, tmp_path):
        from markerfinder.modules.consistency_screen import (
            ConsistencyGrade, MarkerConsistency, write_excluded_profile,
        )

        rows = [MarkerConsistency(
            "M1", None, None, ConsistencyGrade.INCONSISTENT, 1,
            ["opposite_sides"],
        )]
        path = write_excluded_profile(
            rows, str(tmp_path), "mf",
            pis_map={"M1": 7}, category_map={"M1": "ribosomal"},
        )
        body = path.read_text(encoding="utf-8")
        assert path.name == "mf.excluded_profile.tsv"
        assert "M1" in body and "ribosomal" in body
        # A marker with no mapped category must show NA, not an empty cell.
        body2 = write_excluded_profile(
            rows, str(tmp_path), "mf2", pis_map={"M1": 7}, category_map={},
        ).read_text(encoding="utf-8")
        assert "NA" in body2


class TestMolTypeIsAnAssertionNotADisplay:
    """``--mol-type`` documents "force molecule type for sequence validation".

    It reached ``args`` and stopped there: the one place the sequence file is
    read (tree-sequence cross-validation) called the validator without it, so
    every input was auto-detected and the option could not fail a run that
    claimed to check its alphabet.
    """

    PROTEIN_FASTA = ">g1\nMKVLAAGTIGRA\n>g2\nMKVLAAGTIGRB\n"
    TREE = "((g1:0.1,g2:0.2):0.3);"

    def _files(self, tmp_path):
        seq = tmp_path / "tips.faa"
        seq.write_text(self.PROTEIN_FASTA, encoding="utf-8")
        tree = tmp_path / "tree.nwk"
        tree.write_text(self.TREE, encoding="utf-8")
        return str(tree), str(seq)

    def test_forced_dna_on_proteins_is_an_error(self, tmp_path):
        from markerfinder.exceptions import CrossValidationError
        from markerfinder.validation import cross_validate

        tree, seq = self._files(tmp_path)
        with pytest.raises(CrossValidationError, match="mol-type DNA"):
            cross_validate(tree, seq, strict=True, mol_type="DNA")

    def test_forced_protein_on_proteins_passes(self, tmp_path):
        from markerfinder.validation import cross_validate

        tree, seq = self._files(tmp_path)
        report = cross_validate(tree, seq, strict=True, mol_type="protein")
        assert report.is_valid, report.errors

    def test_unforced_input_is_auto_detected_and_only_warns(self, tmp_path):
        """Without the flag the old behaviour stands: no assertion, no refusal."""
        from markerfinder.validation import cross_validate

        tree, seq = self._files(tmp_path)
        report = cross_validate(tree, seq, strict=True)
        assert report.is_valid, report.errors

    def test_the_cli_forwards_the_flag(self, tmp_path, monkeypatch):
        """Wiring lock: the handler must pass what the user asked for."""
        from markerfinder.cli import validation as cli_validation

        tree, seq = self._files(tmp_path)
        seen = {}

        def fake_cross_validate(tree_path, seq_path, strict=True, **kwargs):
            seen.update(kwargs)
            seen["strict"] = strict

            class R:
                tree_tips = {"g1", "g2"}
            return R()

        monkeypatch.setattr("markerfinder.validation.cross_validate",
                            fake_cross_validate)
        args = argparse.Namespace(species_tree=tree, tree=None, sequences=seq,
                                 no_cross_check=False, mol_type="DNA",
                                 skip_length_check=True)
        cli_validation._handle_cross_validation(args)
        assert seen.get("mol_type") == "DNA", seen
        assert seen.get("skip_length_check") is True, seen
