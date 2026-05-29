"""Resolve galaxy lists for ``src.visualization`` and Phase 5 supreme-cohort exports."""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from src.config.thresholds import VIZ_DISPLAY
from src.pipeline.constants import (
    RESULTS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_RELIABLE_FILENAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)

logger = logging.getLogger(__name__)


def cohort_csv_path(project_root: Path, tier: str) -> Path:
    if tier == "evidence":
        name = VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME
    else:
        name = VALIDATION_RELIABLE_FILENAME
    return project_root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME / name


def read_galaxy_names_from_cohort(project_root: Path, *, tier: str, max_n: int) -> list[str]:
    path = cohort_csv_path(project_root, tier)
    if not path.is_file():
        raise FileNotFoundError(f"Missing cohort table for viz auto mode: {path}")
    df = pd.read_csv(path)
    if "galaxy_name" not in df.columns:
        raise KeyError(f"{path} has no galaxy_name column")
    names = [str(x) for x in df["galaxy_name"].tolist() if str(x).strip()]
    return names[: max(0, int(max_n))]


def expand_cli_galaxy_names(
    project_root: Path,
    tokens: Sequence[str],
    *,
    tier: str,
    max_n: int,
) -> list[str]:
    if not tokens:
        return []
    if len(tokens) == 1 and tokens[0].strip().lower() == "auto":
        return read_galaxy_names_from_cohort(project_root, tier=tier, max_n=max_n)
    return [str(t).strip() for t in tokens if str(t).strip()]


def write_viz_manifest(output_dir: Path, written_paths: Sequence[str]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": "1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "written": list(written_paths),
    }
    out = output_dir / "viz_3d_run_manifest.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return out


def run_supreme_museum_batch(*, project_root: Path, ui: Any | None) -> str:
    """
    Phase 5: supreme-cohort PNG + Museum HTML under ``results/graphs/viz_3d/``.

    Returns a short status line for the pipeline closing summary.
    """
    try:
        import pyvista as pv  # noqa: F401
    except ImportError:
        msg = "Phase 5 3D: skipped (pyvista is not installed in this environment)."
        if ui is not None:
            ui.emit(msg)
        else:
            print(msg, file=sys.stderr, flush=True)
        return msg

    from src.reporting.publication_bundle import cohort_csv_path, read_cohort_galaxy_names
    from src.visualization.render import _sanitize_filename, render_galaxy_png

    root = project_root.resolve()
    ev_path = cohort_csv_path(root, "scientific_evidence")
    if not ev_path.is_file():
        msg = f"Phase 5 3D: skipped (missing {ev_path})."
        if ui is not None:
            ui.emit(msg)
        else:
            print(msg, file=sys.stderr, flush=True)
        return msg

    names = read_cohort_galaxy_names(root, "scientific_evidence")
    if not names:
        msg = "Phase 5 3D: supreme cohort CSV has no galaxy_name entries."
        if ui is not None:
            ui.emit(msg)
        else:
            print(msg, file=sys.stderr, flush=True)
        return msg

    ev_df = pd.read_csv(ev_path)
    morph_by_name: dict[str, str | None] = {}
    if "galaxy_name" in ev_df.columns and "morphology_class" in ev_df.columns:
        for _, row in ev_df.iterrows():
            gn = str(row.get("galaxy_name") or "").strip()
            if not gn:
                continue
            mc = row.get("morphology_class")
            morph_by_name[gn] = str(mc).strip() if mc is not None and str(mc).strip() else None

    out = root / RESULTS_DIRNAME / "graphs" / "viz_3d"
    out.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for name in names:
        safe = _sanitize_filename(name)
        target = out / f"{safe}_viz_3d.png"
        try:
            render_galaxy_png(
                root,
                name,
                target,
                morphology_class_override=morph_by_name.get(name),
                gap_glyphs=None,
                wim_rings=None,
                streamlines=None,
                export_html=True,
            )
        except FileNotFoundError as exc:
            warn = f"Phase 5 3D: skipped {name} ({exc})"
            if ui is not None:
                ui.emit(warn)
            else:
                logger.warning("%s", warn)
            continue
        written.append(str(target.resolve()))
        html_p = target.with_suffix(".html")
        if html_p.is_file():
            written.append(str(html_p.resolve()))
        if ui is not None:
            ui.emit(f"Phase 5 3D: wrote {target} (+ HTML)")
        else:
            print(f"Phase 5 3D: wrote {target}", flush=True)

    if written:
        mp = write_viz_manifest(out, written)
        if ui is not None:
            ui.emit(f"Phase 5 3D manifest: {mp}")
    if not written:
        msg = "Phase 5 3D: no exports written (missing per-galaxy analysis folders?)."
        if ui is not None:
            ui.emit(msg)
        return msg
    try:
        rel_disp = str(out.relative_to(root))
    except ValueError:
        rel_disp = str(out)
    n_png = len([p for p in written if p.lower().endswith(".png")])
    return f"Phase 5 3D: {n_png} supreme galaxy export(s) under {rel_disp}."
