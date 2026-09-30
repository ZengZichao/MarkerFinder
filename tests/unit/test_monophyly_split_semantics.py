"""Monophyly screen semantics: root invariance, informative splits, requested rank.

Three claims, each pinned here at the numeric/structural layer while the
end-to-end consequence lives in ``validation/cases/test_v04_taxonomy_inputs.py``
and ``test_v06_hgt_screening.py``:

1. **Root-invariant verdict** (the manual states monophyly is root-invariant,
   and MAD rooting runs upstream of it): a taxon counts when the tree carries
   the split ``{taxon | rest}``, i.e. when the taxon *or* its complement is a
   clade. Grading only rooted-clade membership let the place where the rooting
   rule happened to put the root decide the HGT risk.
2. **An informative denominator**: a taxon is measured only when both sides of
   the split carry >= 2 representatives. With one tip on the other side every
   topology satisfies the claim, so counting it would encode "not measurable"
   as a proportion.
3. **``--monophyly-rank`` means what it says**: naming a rank measures at that
   rank; the scope-derived rank is what ``auto`` (the default) does. A
   documented knob the code overrides with its own automatic choice is the
   defect class this project audits.
"""

import pytest

from markerfinder.config import HGTConfig
from markerfinder.taxonomy import (
    AUTO_NO_SCOPE_RANK,
    AUTO_RANK,
    load_taxonomy_table,
    monophyly_proportion,
)
from markerfinder.modules.hgt_filter import PhylogeneticHGTDetector


TAX_2X2 = {
    "t1": {"domain": "Bacteria", "phylum": "P1", "class": "C1", "order": "O1",
           "family": "F1", "genus": "G1"},
    "t2": {"domain": "Bacteria", "phylum": "P1", "class": "C1", "order": "O1",
           "family": "F1", "genus": "G1"},
    "t3": {"domain": "Bacteria", "phylum": "P1", "class": "C1", "order": "O1",
           "family": "F2", "genus": "G2"},
    "t4": {"domain": "Bacteria", "phylum": "P1", "class": "C1", "order": "O1",
           "family": "F2", "genus": "G2"},
}


class TestRootInvariance:
    def test_the_same_unrooted_topology_gives_the_same_proportion(self):
        """Two Newick renderings of ONE unrooted tree, two different rootings.

        ``((t1,t2),(t3,t4))`` rooted on its internal edge, and the same quartet
        rooted on the branch leading to t1 — which makes ``(t3,t4)`` a clade and
        leaves ``{t1,t2}`` as the complement of a clade instead. The split set is
        identical, so the proportion must be identical — and it must be 1.0,
        because nothing about this topology conflicts with the taxonomy.
        """
        as_rooted_between = "((t1:0.1,t2:0.2):0.3,(t3:0.3,t4:0.4):0.5);"
        rooted_on_the_t1_branch = "(t1:0.1,(t2:0.2,(t3:0.3,t4:0.4):0.5):0.6);"

        p1, n1, m1 = monophyly_proportion(as_rooted_between, TAX_2X2, "family")
        p2, n2, m2 = monophyly_proportion(rooted_on_the_t1_branch, TAX_2X2,
                                        "family")

        assert (n1, m1) == (2, 2), (n1, m1)
        assert p1 == 1.0, p1
        # No node of the second rendering has leaf set F1; its complement F2 is
        # The clade. Rooted-clade membership alone reported 1/2 here, turning a
        # Concordant marker into an HGT call because of where the root sits.
        assert (n2, m2) == (2, 2), (n2, m2)
        assert p2 == 1.0, p2

    def test_genuine_discordance_still_scores_zero(self):
        """The complement rule must not rubber-stamp a conflicting tree."""
        nwk = "((t1:0.1,t3:0.2):0.3,(t2:0.3,t4:0.4):0.5);"
        prop, n_total, n_mono = monophyly_proportion(nwk, TAX_2X2, "family")
        assert (n_total, n_mono) == (2, 0), (n_total, n_mono)
        assert prop == 0.0


class TestInformativeDenominator:
    def test_a_single_tip_on_the_other_side_is_not_a_measurement(self):
        """7-vs-1 sampling cannot falsify anything about the 7.

        The group of all tips except one is a split of every unrooted topology,
        so it is neither a pass nor a failure: it is not measurable, and the
        caller reports it as such instead of a proportion.
        """
        tax = dict(TAX_2X2)
        tax["t5"] = {"domain": "Bacteria", "phylum": "P2", "class": "C2",
                     "order": "O2", "family": "F3", "genus": "G3"}
        nwk = "((t1:0.1,t2:0.2):0.3,(t3:0.3,t4:0.4):0.5,t5:0.9);"
        # Phylum: P1 (4 tips) vs P2 (1 tip) -> not informative either way.
        prop, n_total, n_mono = monophyly_proportion(nwk, tax, "phylum")
        assert prop is None, prop
        assert (n_total, n_mono) == (0, 0)

    def test_a_group_covering_every_tip_is_not_a_measurement(self):
        prop, n_total, _ = monophyly_proportion(
            "((t1:0.1,t2:0.2):0.3,(t3:0.3,t4:0.4):0.5);", TAX_2X2, "domain")
        assert prop is None
        assert n_total == 0


class TestRequestedRankIsAuthoritative:
    @staticmethod
    def _detector():
        return PhylogeneticHGTDetector(HGTConfig())

    def _measure(self, level):
        nwk = "((t1:0.1,t2:0.2):0.3,(t3:0.3,t4:0.4):0.5);"
        return self._detector().detect_monophyly(
            "M1", nwk, TAX_2X2, level, 0.5)

    def test_auto_derives_the_rank_from_the_scope(self):
        res = self._measure(AUTO_RANK)
        assert res is not None
        # All four tips share one order, so the scope is 'order' and the
        # Measured rank is one below it.
        assert res.test_level_used == "family", res.test_level_used
        assert res.cross_rank_comparison is False
        assert res.monophyly_proportion == 1.0

    def test_a_named_rank_is_measured_even_though_the_scope_says_otherwise(self):
        res = self._measure("genus")
        assert res is not None
        assert res.test_level_used == "genus", res.test_level_used
        assert res.cross_rank_comparison is False
        assert res.monophyly_proportion == 1.0

    def test_a_named_rank_that_cannot_be_measured_is_flagged_cross_rank(self):
        """One phylum covers the whole tree: no informative split at phylum.

        The walk may settle on a rank that can be measured, but that proportion
        is not comparable with a threshold set for the requested rank,
        so it is surfaced rather than graded.
        """
        res = self._measure("phylum")
        assert res is not None
        assert res.test_level_used != "phylum", res.test_level_used
        assert res.cross_rank_comparison is True


class TestRankVocabularyIsOneWord:
    """Config.py keeps the default as a literal; guard that it is the constant."""

    def test_config_default_is_the_auto_token(self):
        assert HGTConfig().monophyly_rank == AUTO_RANK == "auto"

    def test_auto_resolves_to_a_real_rank_without_a_scope(self):
        assert AUTO_NO_SCOPE_RANK in {
            "domain", "phylum", "class", "order", "family", "genus", "species"}

    def test_parser_offers_auto_and_every_rank(self):
        from markerfinder.cli.parser import _build_parser

        for action in _build_parser()._actions:      # Noqa: SLF001
            if action.dest == "monophyly_rank":
                assert action.default == AUTO_RANK
                assert AUTO_RANK in action.choices
                assert set(action.choices) >= {
                    "domain", "phylum", "class", "order", "family", "genus",
                    "species"}
                return
        pytest.fail("--monophyly-rank is no longer a parser option")


class TestCustomRanksAreParsed:
    """``--taxonomy-levels`` used to be recorded and never applied."""

    _ROW = ("{tip}\tk__Bacteria;d__Bacteria;p__P1;o__O1;f__F1;g__G1")
    TABLE = ("genome_id\ttaxonomy\n"
             + _ROW.format(tip="t1") + "\n"
             + _ROW.format(tip="t2") + "\n")

    def _write(self, path):
        path.write_text(self.TABLE, encoding="utf-8")
        return str(path)

    def test_custom_rank_reaches_the_parsed_mapping(self, tmp_path):
        table = self._write(tmp_path / "tax.tsv")
        got = load_taxonomy_table(table, custom_levels="kingdom:k__")
        assert got["t1"]["kingdom"] == "Bacteria", got["t1"]

    def test_without_the_option_the_rank_is_dropped_and_said(self, tmp_path,
                                                             caplog):
        table = self._write(tmp_path / "tax.tsv")
        with caplog.at_level("WARNING"):
            got = load_taxonomy_table(table)
        assert "kingdom" not in got["t1"], got["t1"]
        assert "Unknown level prefix" in caplog.text

    def test_the_extension_is_announced_once_per_table(self, tmp_path, caplog):
        table = self._write(tmp_path / "tax.tsv")
        with caplog.at_level("INFO"):
            load_taxonomy_table(table, custom_levels="kingdom:k__")
        announced = [r for r in caplog.records if "kingdom" in r.getMessage()]
        assert len(announced) == 1, [r.getMessage() for r in announced]
