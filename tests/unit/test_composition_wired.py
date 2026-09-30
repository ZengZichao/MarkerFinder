"""Composition evidence must reach the products.

The metrics in ``composition.py`` were implemented and unit-tested, but
nothing in the pipeline ever called them, and ``--hgt-steps`` still rejected
``composition`` — so 's product-level promise (outlier list, parallel
columns) did not exist. These tests exercise the real writer and pin the
wiring so the channel cannot go dead again.
"""

from __future__ import annotations

import inspect
from pathlib import Path

from markerfinder.config import ReportConfig
from markerfinder.modules.composition import (
    marker_composition_metrics,
    rcv,
    write_composition_tsv,
)

_HOMOGENEOUS = {
    "g1": "ACDEFGHIKLMNPQRSTVWY" * 8,
    "g2": "ACDEFGHIKLMNPQRSTVWY" * 8,
    "g3": "ACDEFGHIKLMNPQRSTVWY" * 8,
    "g4": "ACDEFGHIKLMNPQRSTVWY" * 8,
}


def _skewed():
    # Same 20 residues, but one genome is dominated by a handful of them:
    # A composition outlier, not a divergence outlier.
    return dict(
        _HOMOGENEOUS,
        g5="PPPPPPPPPPPPPPPPCCCCCCCCCCCCCCC" * 5,
    )


class TestMarkerCompositionMetrics:
    def test_outlier_marker_scores_higher_rcv(self):
        baseline = marker_composition_metrics(_HOMOGENEOUS)["rcv"]
        planted = marker_composition_metrics(_skewed())["rcv"]
        assert baseline is not None and planted is not None
        assert planted > baseline

    def test_protein_input_reports_gc_as_not_applicable(self):
        """GC on amino-acid data is NOT_APPLICABLE — never a fake 0.0."""
        metrics = marker_composition_metrics(_HOMOGENEOUS, protein=True)
        assert metrics["gc_bias"] is None

    def test_single_sequence_is_unmeasurable_not_zero(self):
        metrics = marker_composition_metrics({"only": "AAAA"})
        assert metrics["rcv"] is None and metrics["n_sequences"] == 1


class TestCompositionTsv:
    def test_worst_outlier_heads_the_list(self, tmp_path):
        rows = [
            {"marker_id": "clean_a", "rcv": 0.01, "gc_bias": None,
             "n_sequences": 4},
            {"marker_id": "outlier", "rcv": 0.92, "gc_bias": None,
             "n_sequences": 4},
            {"marker_id": "clean_b", "rcv": 0.02, "gc_bias": None,
             "n_sequences": 4},
        ]
        path = write_composition_tsv(
            rows, str(tmp_path), "mf", warn_threshold=0.15
        )
        assert path == tmp_path / "Phase5_reports" / "mf.composition.tsv"
        body = path.read_text(encoding="utf-8").strip().splitlines()
        assert body[0].split("\t") == [
            "marker_id", "rcv", "gc_bias", "n_sequences", "outlier_flag"
        ]
        assert body[1].split("\t")[0] == "outlier"
        assert body[1].split("\t")[-1] == "warn"
        assert body[-1].split("\t")[-1] == "NA"

    def test_missing_values_render_as_na(self, tmp_path):
        path = write_composition_tsv(
            [{"marker_id": "M1", "rcv": None, "gc_bias": None,
              "n_sequences": 0}],
            str(tmp_path), "mf", warn_threshold=0.15,
        )
        line = path.read_text(encoding="utf-8").strip().splitlines()[1]
        assert line.split("\t") == ["M1", "NA", "NA", "0", "NA"]

    def test_screen_defaults_to_off(self):
        """Opt-in only: shipping behaviour must stay at the baseline."""
        assert ReportConfig().composition_screen is False


class TestWiringLocks:
    """Source locks: the channel must stay connected to the run."""

    def test_hgt_steps_accepts_composition(self):
        import markerfinder.cli.config_build as config_build

        src = inspect.getsource(config_build)
        assert '"composition"' in src, "--hgt-steps no longer accepts composition"
        assert "composition_screen=" in src, "step name not mapped to config"

    def test_sidecar_carries_the_metrics(self):
        import markerfinder.utils.gtdb_tk_markers as gtdb

        src = inspect.getsource(gtdb)
        assert "marker_composition_metrics" in src, "metrics not computed"
        assert '"rcv"' in src, "rcv not stored in the sidecar"

    def test_pipeline_writes_the_product(self):
        import markerfinder.pipeline as pipeline_module

        src = inspect.getsource(pipeline_module)
        # The orchestrator no longer calls the writer itself -- it goes
        # Through the composition stage, which owns the artifact. The lock is
        # Unchanged in intent (the switch must reach a real product), only in
        # Which name it greps for; see tests/unit/test_stage_owns_its_products.py
        # For the byte-identical-output guarantee.
        assert "run_composition_stage(" in src, "composition file never written"
        assert "composition_screen" in src, "opt-in switch unread"

    def test_cli_help_advertises_the_step(self):
        import markerfinder.cli.parser as parser_module

        assert "composition" in inspect.getsource(parser_module)
