"""Banner and startup display for MarkerFinder."""

import re
import sys

from markerfinder._version import (
    __author__,
    __email__,
    __orcid__,
    get_version_string,
)

BANNER = r"""
 __  __       _ _     _ _____           _   _ _           _
|  \/  | __ _| (_) __| |  _ \ ___  ___| |_(_) |__  _   _| |_
| |\/| |/ _` | | |/ _` | |_) / _ \/ __| __| | '_ \| | | | __|
| |  | | (_| | | | (_| |  _ <  __/\__ \ |_| | |_) | |_| | |_
|_|  |_|\__,_|_|_|\__,_|_| \_\___||___/\__|_|_.__/ \__,_|\__|
  ___  _                    _       _
 | _ \| |__  __ _ _  _ __ _| |_  __| |___ _ _
 |  _/| '_ \/ _` | || / _` | ' \/ _` / -_) '_|
 |_|  |_.__/\__, |\_,_\__,_|_||_\__,_\___|_|
            |___/  Finder
"""

# SPDX-style identifiers for the libraries the distribution actually requires.
# A name that is not here prints as "see upstream" rather than being guessed.
LICENSE_BY_PACKAGE = {
    "biopython": "BSD-3-Clause",
    "pyyaml": "MIT",
    "ete3": "GPL-3.0",
    "tomli": "MIT",
}

_DISPLAY_NAMES = {"biopython": "Biopython", "pyyaml": "PyYAML",
                  "ete3": "ete3", "tomli": "tomli"}

# Only used when no distribution metadata can be read (a bare source tree).
# It must match pyproject.toml: validation/cases/test_v01_environment_and_selfcheck.py
# Compares what the banner prints against the installed distribution's
# Requirements, so a drift is caught by a run rather than by review.
_FALLBACK_REQUIREMENTS = ["biopython", "pyyaml", "ete3", "tomli"]


def runtime_dependencies() -> list:
    """Normalized names of the declared runtime dependencies.

    Read from the installed distribution metadata rather than typed, because a
    typed list is how this banner ended up attributing licences to pandas,
    numpy, scipy and rich -- libraries removed from the dependency set in
     and imported nowhere in the package. The startup screen was
    claiming them as third-party libraries of MarkerFinder on every run.
    """
    try:
        from importlib import metadata

        requires = metadata.distribution("markerfinder").requires or []
    except Exception:  # Noqa: BLE001 - no metadata: use the declared set
        return list(_FALLBACK_REQUIREMENTS)

    names: list = []
    for entry in requires:
        # "ete3>=3.1", "tomli>=2.0;python_version<'3.11'" -> the bare name.
        # Optional-dependency extras are excluded: the banner lists what a run
        # Needs, not what a developer installed alongside it.
        if "extra ==" in entry or "extra>=" in entry:
            continue
        match = re.match(r"[A-Za-z0-9_.\-]+", entry.strip())
        if not match:
            continue
        name = match.group(0).lower().replace("_", "-")
        if name not in names:
            names.append(name)
    return names


def third_party_licenses() -> list:
    """``(display name, license id)`` for each declared runtime dependency."""
    return [
        (_DISPLAY_NAMES.get(name, name),
         LICENSE_BY_PACKAGE.get(name, "see upstream"))
        for name in runtime_dependencies()
    ]


# Module attribute kept for callers that imported it, now derived not typed.
THIRD_PARTY_LICENSES = third_party_licenses()


def print_banner() -> None:
    sys.stdout.write(BANNER + "\n")
    sys.stdout.write("  MIT License\n")
    sys.stdout.write(f"  Version: {get_version_string()}\n")
    sys.stdout.write(f"  Author: {__author__} <{__email__}> (ORCID {__orcid__})\n")
    sys.stdout.write("  Runtime dependencies (declared):\n")
    for name, license_id in third_party_licenses():
        sys.stdout.write(f"    {name:20s} {license_id}\n")
    sys.stdout.write("\n")
    sys.stdout.flush()
