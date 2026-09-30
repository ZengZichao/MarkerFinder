"""Version and authorship information for MarkerFinder.

The version is a single constant: what ``markerfinder --version``, the startup
banner, the run metadata and the build metadata all report. Nothing here is
derived from a VCS checkout, so an unpacked source tree and an installed
distribution state the same thing.
"""

__version__ = "0.1.0"

__author__ = "Zengzichao"
__email__ = "zengzichao@sjtu.edu.cn"
__orcid__ = "0000-0001-6553-970X"


def get_version_string() -> str:
    """The version as it appears in the banner, the reports and the run state."""
    return __version__
