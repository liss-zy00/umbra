"""
Merge NED-hosted disk geometry and kinematic context into the SPARC-name-keyed PA overlay JSON.

Uses ``astroquery.ipac.ned.Ned`` (not the deprecated ``astroquery.ned`` import path). Intended
for offline curation only — import :func:`refresh_disk_pa_overlay_from_validation_csv`
from a REPL or notebook; the pipeline runtime never calls this module.

Manual PA overlay policy
------------------------
Any hand-edited rows in ``disk_position_angle_overlay.json`` (whether produced by this helper
or typed directly) must follow the same **literature-only** rule as ``sources/morphology/README.md``:
``position_angle_deg`` is accepted only when a reviewer can trace it to a verifiable primary
reference recorded in ``source`` (NED/LEDA bibcode, published table, DOI, or stable catalogue
field). Do not insert object PA values inferred from uncited secondary material.
"""

from __future__ import annotations

import json
import time
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.cosmology import Planck18
from astroquery.ipac.ned import Ned

from src.pipeline.constants import (
    MORPHOLOGY_DISK_PA_OVERLAY_JSON,
    RESULTS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_RELIABLE_FILENAME,
)
from src.providers.ned import _build_candidate_queries, _default_alias_map_path, _load_ned_alias_config

_DEFAULT_HTTP_PAUSE_SEC: float = 0.45
_MAX_RECENT_KINEMATIC_ROWS: int = 8


def _refcode_year(ref: Any) -> int:
    s = str(ref).strip()
    if len(s) >= 4 and s[:4].isdigit():
        y = int(s[:4])
        if 1800 <= y <= 2100:
            return y
    return 0


def _as_optional_float(val: Any) -> float | None:
    if val is None:
        return None
    s = str(val).strip()
    if s in ("", "--", "nan"):
        return None
    try:
        x = float(val)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(x):
        return None
    return float(x)


def _luminosity_distance_mpc(z: float | None) -> float | None:
    if z is None or z <= 0.0 or not np.isfinite(z):
        return None
    try:
        return float(Planck18.luminosity_distance(z).to_value(u.Mpc))
    except Exception:
        return None


def _first_resolved_ned_query(candidates: Iterable[tuple[str, str]]) -> tuple[str, Any] | None:
    for query, _tag in candidates:
        try:
            main = Ned.query_object(query)
        except Exception:
            continue
        if main is None or len(main) == 0:
            continue
        return query, main
    return None


def _pick_preferred_pa_deg(diam: Any) -> dict[str, Any] | None:
    if diam is None or len(diam) == 0:
        return None
    best: tuple[tuple[int, str], dict[str, Any]] | None = None
    for i in range(len(diam)):
        ref = diam["Refcode"][i]
        ned_pa = _as_optional_float(diam["NED Position Angle"][i]) if "NED Position Angle" in diam.colnames else None
        pub_pa = _as_optional_float(diam["Position Angle"][i]) if "Position Angle" in diam.colnames else None
        pa = ned_pa if ned_pa is not None else pub_pa
        if pa is None:
            continue
        unc_col = "NED Position Angle Uncertainty" if ned_pa is not None else "Position Angle Uncertainty"
        unc_raw = diam[unc_col][i] if unc_col in diam.colnames else None
        unc = _as_optional_float(unc_raw)
        row = {
            "position_angle_deg": pa,
            "position_angle_uncertainty_deg": unc,
            "refcode": str(ref).strip(),
            "ned_column": "NED Position Angle" if ned_pa is not None else "Position Angle",
        }
        key = (_refcode_year(ref), str(ref).strip())
        if best is None or key > best[0]:
            best = (key, row)
    return best[1] if best else None


def _kinematic_rows_recent(red: Any) -> list[dict[str, Any]]:
    if red is None or len(red) == 0:
        return []
    rows: list[tuple[tuple[int, str], dict[str, Any]]] = []
    for i in range(len(red)):
        ref = red["Refcode"][i]
        z = _as_optional_float(red["Published Redshift"][i]) if "Published Redshift" in red.colnames else None
        z_unc = _as_optional_float(red["Published Redshift Uncertainty"][i]) if "Published Redshift Uncertainty" in red.colnames else None
        v = _as_optional_float(red["Published Velocity"][i]) if "Published Velocity" in red.colnames else None
        v_unc = _as_optional_float(red["Published Velocity Uncertainty"][i]) if "Published Velocity Uncertainty" in red.colnames else None
        if z is None and v is None:
            continue
        dl = _luminosity_distance_mpc(z) if z is not None else None
        if dl is not None:
            dl = round(dl, 6)
        row = {
            "refcode": str(ref).strip(),
            "published_redshift": z,
            "published_redshift_uncertainty": z_unc,
            "published_velocity_kms": v,
            "published_velocity_uncertainty_kms": v_unc,
            "luminosity_distance_mpc_planck18": dl,
        }
        key = (_refcode_year(ref), str(ref).strip())
        rows.append((key, row))
    rows.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in rows[:_MAX_RECENT_KINEMATIC_ROWS]]


def _fetch_for_sparc_label(
    sparc_name: str,
    *,
    alias_config: dict[str, Any],
    pause_sec: float,
) -> dict[str, Any]:
    candidates = _build_candidate_queries(
        sparc_name,
        alias_map=alias_config["alias_map"],
        pattern_rules=alias_config["pattern_rules"],
        prefix_attempts=alias_config["ned_prefix_attempts"],
    )
    resolved = _first_resolved_ned_query(candidates)
    time.sleep(pause_sec)
    if resolved is None:
        return {"galaxy": sparc_name, "status": "unresolved", "ned_query_attempted": [q for q, _ in candidates]}
    query, main = resolved
    obj_name = str(main["Object Name"][0]).strip()
    out: dict[str, Any] = {
        "galaxy": sparc_name,
        "status": "ok",
        "ned_query_string": query,
        "ned_resolved_object_name": obj_name,
    }
    try:
        diam = Ned.get_table(query, table="diameters")
    except Exception as exc:
        diam = None
        out["diameters_error"] = str(exc)
    time.sleep(pause_sec)
    try:
        red = Ned.get_table(query, table="redshifts")
    except Exception as exc:
        red = None
        out["redshifts_error"] = str(exc)
    time.sleep(pause_sec)

    pa_pick = _pick_preferred_pa_deg(diam)
    if pa_pick:
        ref = pa_pick["refcode"]
        unc = pa_pick.get("position_angle_uncertainty_deg")
        out["position_angle_deg"] = float(pa_pick["position_angle_deg"])
        out["position_angle_uncertainty_deg"] = unc if unc is not None else None
        out["source"] = f"NED diameters ({pa_pick['ned_column']}); refcode={ref}"
    else:
        out["position_angle_deg"] = None
        out["position_angle_uncertainty_deg"] = None
        out["source"] = None
        out["position_angle_note"] = "No finite Position Angle / NED Position Angle in NED diameters table."

    out["ned_kinematic_measurements_by_recency"] = _kinematic_rows_recent(red)
    return out


def _merge_positions(
    base: dict[str, Any],
    fetch: dict[str, Any],
) -> dict[str, Any]:
    g = fetch["galaxy"]
    entry: dict[str, Any] = dict(base.get(g, {}))
    entry.pop("position_angle_source", None)
    entry["ned_fetch_status"] = fetch.get("status")
    entry["ned_query_string"] = fetch.get("ned_query_string")
    entry["ned_resolved_object_name"] = fetch.get("ned_resolved_object_name")
    entry["ned_kinematic_measurements_by_recency"] = fetch.get("ned_kinematic_measurements_by_recency") or []
    if fetch.get("status") != "ok":
        entry["ned_fetch_note"] = {k: fetch[k] for k in ("galaxy", "status", "ned_query_attempted") if k in fetch}
        return entry
    if fetch.get("position_angle_deg") is not None:
        entry["position_angle_deg"] = fetch["position_angle_deg"]
        u = fetch.get("position_angle_uncertainty_deg")
        if u is not None and np.isfinite(float(u)):
            entry["position_angle_uncertainty_deg"] = float(u)
        else:
            entry.pop("position_angle_uncertainty_deg", None)
        entry["source"] = str(fetch.get("source") or "NED")
    return entry


def refresh_disk_pa_overlay_from_validation_csv(
    *,
    project_root: Path,
    csv_path: Path | None = None,
    out_json: Path | None = None,
    ned_alias_path: Path | None = None,
    pause_sec: float = _DEFAULT_HTTP_PAUSE_SEC,
    dry_run: bool = False,
) -> list[dict[str, Any]]:
    """
    For each ``galaxy_name`` in the validation CSV, query NED and merge into the PA overlay JSON.

    ``project_root`` is the repository root (parent of ``sources/``). Returns a small summary
    list for logging. When ``dry_run`` is True, the summary is still computed but the JSON file
    is not written.
    """
    csv_path = csv_path or (
        project_root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME / VALIDATION_RELIABLE_FILENAME
    )
    out_path = out_json or (project_root / MORPHOLOGY_DISK_PA_OVERLAY_JSON)
    alias_path = ned_alias_path or _default_alias_map_path()
    alias_config = _load_ned_alias_config(Path(alias_path))

    df = pd.read_csv(csv_path)
    if "galaxy_name" not in df.columns:
        raise ValueError("CSV must contain a galaxy_name column.")
    names = [str(x).strip() for x in df["galaxy_name"].tolist() if str(x).strip()]
    if not names:
        raise ValueError("No galaxy names in CSV.")

    if out_path.is_file():
        payload = json.loads(out_path.read_text(encoding="utf-8"))
    else:
        payload = {}
    positions: dict[str, Any] = payload.get("positions") if isinstance(payload.get("positions"), dict) else {}

    summary: list[dict[str, Any]] = []
    for name in names:
        fetch = _fetch_for_sparc_label(name, alias_config=alias_config, pause_sec=float(pause_sec))
        merged = _merge_positions(positions, fetch)
        positions[name] = merged
        summary.append(
            {
                "galaxy": name,
                "status": fetch.get("status"),
                "pa_deg": fetch.get("position_angle_deg"),
                "kinematic_rows": len(fetch.get("ned_kinematic_measurements_by_recency") or []),
            }
        )

    today = date.today().isoformat()
    payload["positions"] = positions
    payload["_schema"] = {
        "<sparc_name>": {
            "position_angle_deg": "float (east of north, deg) — required for pipeline PA overlay",
            "position_angle_uncertainty_deg": "float or NaN — 1-sigma when published",
            "source": "string — catalogue / NED provenance (read by SparcMorphologyProvider as position_angle_source)",
            "ned_fetch_status": "ok | unresolved",
            "ned_query_string": "string — NED query that resolved",
            "ned_resolved_object_name": "string — NED Object Name",
            "ned_kinematic_measurements_by_recency": (
                "list[{refcode, published_redshift, published_redshift_uncertainty, "
                "published_velocity_kms, published_velocity_uncertainty_kms, "
                "luminosity_distance_mpc_planck18}] — audit only; SPARC D_Mpc remains canonical in the engine"
            ),
        }
    }
    prev_note = str(payload.get("_note") or "")
    ned_line = (
        f"NED bulk row ({today}, validation_reliable): disk PA from NED ``diameters`` (latest bibcode among finite "
        f"NED/ published PA); kinematic z and heliocentric velocity rows from NED ``redshifts`` (newest bibcodes "
        f"first, cap {_MAX_RECENT_KINEMATIC_ROWS}); luminosity_distance_mpc_planck18 is astropy.cosmology.Planck18 "
        "for each row's published_redshift (derived, not a direct NED distance modulus)."
    )
    if ned_line not in prev_note:
        payload["_note"] = (prev_note + " " + ned_line).strip()
    prev_prov = str(payload.get("_provenance") or "")
    prov_add = (
        f"NED/IPAC via astroquery.ipac.ned ({today}, UTC workflow); name candidates mirror "
        "sources/ned/ned_alias_map.json. Derived luminosity distances: Planck18 (astropy)."
    )
    if prov_add not in prev_prov:
        payload["_provenance"] = (prev_prov + " | " + prov_add).strip(" |")

    if not dry_run:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    return summary
