"""Python 3.13+ compatibility shim for the deprecated stdlib ``cgi`` module.

Python 3.13 removed :mod:`cgi` (PEP 594). ``ete3`` — a hard runtime dependency
of MarkerFinder — still imports it at package import time via
``ete3.webplugin.webapp``, so ``import ete3`` raises ``ModuleNotFoundError`` on
3.13+. That silently disabled the entire MAD/monophyly code path, which then
reported every affected marker as "unmeasured" instead of failing loudly.

Rather than pinning the project to an EOL interpreter, we install a minimal
``cgi`` stand-in **only when the module is genuinely absent**, before ``ete3``
is imported. On 3.12 and earlier this is a no-op because the real ``cgi``
already exists.

The shim covers the single symbol ``ete3`` touches (``cgi.FieldStorage``); it is
not a general-purpose reimplementation, and it is not exported for other code.
"""

from __future__ import annotations

import sys
import types


def _install_cgi_shim() -> bool:
    """Register a minimal ``cgi`` module when the stdlib one is missing.

    Returns ``True`` if the shim was installed, ``False`` if the real module is
    already importable (or could not be provided, in which case ``ete3`` will
    fail loudly rather than silently).
    """
    if "cgi" in sys.modules:
        return False
    try:  # Real module available (Python <= 3.12) — nothing to do.
        import cgi  # noqa: F401, PLC0415

        return False
    except ModuleNotFoundError:
        pass

    shim = types.ModuleType("cgi")

    class FieldStorage:  # noqa: D401 — ete3 only needs the name to exist.
        """Placeholder for :class:`cgi.FieldStorage`.

        ``ete3.webplugin.webapp`` subclasses/imports this symbol but the
        round-trip test and MarkerFinder's own code never exercise the web
        plugin's multipart parsing, so an inert stand-in is sufficient to let
        ``import ete3`` succeed.
        """

        def __init__(self, *args, **kwargs) -> None:
            raise NotImplementedError(
                "cgi.FieldStorage is not available on Python 3.13+. "
                "The ete3 web plugin is not used by MarkerFinder."
            )

    shim.FieldStorage = FieldStorage  # type: ignore[attr-defined]
    sys.modules["cgi"] = shim
    return True


_CGI_SHIM_ACTIVE = _install_cgi_shim()

__all__ = ["_CGI_SHIM_ACTIVE"]
