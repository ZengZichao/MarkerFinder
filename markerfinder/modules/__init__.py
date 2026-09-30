"""The five analysis-stage modules wired into ``MarkerFinderPipeline``.

``pipeline.py`` imports exactly these five and no other stage:
``mag_optimization`` (Phase 0), ``marker_selection`` (Phase 1/1.5),
``hgt_filter`` (Phase 2), ``phylogenetic_inference`` (Phase 3),
``report_generator`` (Phase 4).

``OrthologResolver`` is imported here for backward compatibility only. It is a
**reserved, not wired** interface (see its module docstring,
``phases.PHASES["1.5"]`` and the ``--check`` disclosure), so it is deliberately
kept out of ``__all__``: ``from markerfinder.modules import *`` must not present
an unused capability as one of the pipeline stages.
"""

from markerfinder.modules.marker_selection import AdaptiveMarkerSelectionModule
from markerfinder.modules.hgt_filter import HGTFilterModule, HGTDecisionEngine
from markerfinder.modules.phylogenetic_inference import PhylogeneticInferenceModule
from markerfinder.modules.report_generator import ReportGeneratorModule
from markerfinder.modules.mag_optimization import MAGOptimizationModule
from markerfinder.modules.ortholog_resolver import OrthologResolver  # Reserved, not wired

__all__ = [
    "AdaptiveMarkerSelectionModule",
    "HGTFilterModule",
    "HGTDecisionEngine",
    "PhylogeneticInferenceModule",
    "ReportGeneratorModule",
    "MAGOptimizationModule",
]
