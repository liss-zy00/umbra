"""
Reporting layer — disk persistence, narrative generation, and matplotlib output.

The subpackage layout reflects three orthogonal concerns:

* :mod:`persistence` — tabular and JSON disk writes.
* :mod:`narrative`   — human-readable Scientific-Evidence narrative.
* :mod:`plotting`    — matplotlib figures and their automatic Markdown companions.
"""

from __future__ import annotations

from .narrative import NarrativeGenerator
from .persistence import DataPersister, ResultPersister
from .plotting import GalaxyPlotter

__all__ = ["DataPersister", "GalaxyPlotter", "NarrativeGenerator", "ResultPersister"]
