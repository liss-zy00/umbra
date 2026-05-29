"""
Tabular and JSON disk persistence.

* :mod:`exporter` — :class:`ResultPersister`: saves per-galaxy CSV/JSON, master
  inventory, validation cohort CSVs, scientific-evidence report, anchor sweep, and reliability partitions.
* :mod:`raw_data` — :class:`DataPersister`: thin writer for per-galaxy CSV/JSON and ad-hoc
  ``step_validation/`` text exports; population synthesis tables prefer :class:`ResultPersister`.
"""

from __future__ import annotations

from .exporter import ResultPersister
from .raw_data import DataPersister

__all__ = ["DataPersister", "ResultPersister"]
