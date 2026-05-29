"""
Executable entry point for the SPARC + Planck pipeline.

Adds the repository root to ``sys.path`` so ``python src/main.py`` resolves ``src.*`` imports.
"""

from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.pipeline.cli import parse_args
from src.pipeline.runner import run_pipeline

if __name__ == "__main__":
    run_pipeline(parse_args())
