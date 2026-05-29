"""Extragalactic Distance Database (EDD) table I/O and luminosity-distance hierarchy.

This module is intentionally **catalogue-only**: it parses the HyperLEDA-formatted EDD CSV
drops and applies the tier order TRGB → CF4 individual → Cosmicflows-3 described in
``sources/extragalactic_distances/README.md``. Cosmology-specific inferences beyond the empty-space modulus
relation live in :mod:`src.config.thresholds`.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.config.thresholds import (
    luminosity_distance_mpc_from_distance_modulus,
    luminosity_distance_uncertainty_mpc_from_modulus,
)
from src.providers.sparc import SparcMorphologyProvider

logger = logging.getLogger(__name__)

_HYPERLEDA_DATA_SKIPROWS = 5


def _finite_positive(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(v) or v <= 0.0:
        return None
    return float(v)


def _finite_non_negative(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(v) or v < 0.0:
        return None
    return float(v)


def _read_trgb(path: Path) -> pd.DataFrame:
    cols = ["pgc", "name_cmd", "m_trgb", "dm_tip", "edm_lo", "edm_hi"]
    df = pd.read_csv(path, skiprows=_HYPERLEDA_DATA_SKIPROWS, header=None, names=cols)
    df["pgc"] = pd.to_numeric(df["pgc"], errors="coerce").astype("Int64")
    return df.dropna(subset=["pgc"]).drop_duplicates(subset=["pgc"], keep="first")


def _read_cf4(path: Path) -> pd.DataFrame:
    cols = ["pgc", "pgc1", "vcmb", "DM", "eDM", "DMtrgb", "eDMt", "DMcep", "eDMcep"]
    df = pd.read_csv(path, skiprows=_HYPERLEDA_DATA_SKIPROWS, header=None, names=cols)
    df["pgc"] = pd.to_numeric(df["pgc"], errors="coerce").astype("Int64")
    return df.dropna(subset=["pgc"]).drop_duplicates(subset=["pgc"], keep="first")


def _read_cf3(path: Path) -> pd.DataFrame:
    cols = ["pgc", "Dist", "DM", "eDM", "Ty", "Name"]
    df = pd.read_csv(path, skiprows=_HYPERLEDA_DATA_SKIPROWS, header=None, names=cols)
    df["pgc"] = pd.to_numeric(df["pgc"], errors="coerce").astype("Int64")
    return df.dropna(subset=["pgc"]).drop_duplicates(subset=["pgc"], keep="first")


@dataclass(frozen=True, slots=True)
class EddDistanceHit:
    """One resolved luminosity distance from the EDD hierarchy (tiers 1–3)."""

    tier: int
    distance_modulus_mag: float
    sigma_distance_modulus_mag: float
    distance_mpc: float
    sigma_distance_mpc: float
    source_table: str


class EddProvider:
    """In-memory EDD tables indexed by PGC."""

    def __init__(self, trgb: pd.DataFrame, cf4: pd.DataFrame, cf3: pd.DataFrame) -> None:
        self._trgb = trgb.set_index("pgc")
        self._cf4 = cf4.set_index("pgc")
        self._cf3 = cf3.set_index("pgc")

    @classmethod
    def from_csv_paths(cls, trgb_csv: Path, cf4_csv: Path, cf3_csv: Path) -> EddProvider:
        return cls(_read_trgb(trgb_csv), _read_cf4(cf4_csv), _read_cf3(cf3_csv))

    def lookup_distance_for_pgc(self, pgc: int) -> EddDistanceHit | None:
        """Return the highest-precedence finite distance for ``pgc``, or ``None``."""
        if pgc in self._trgb.index:
            row = self._trgb.loc[pgc]
            dm = _finite_positive(row.get("dm_tip"))
            if dm is not None:
                lo = _finite_non_negative(row.get("edm_lo"))
                hi = _finite_non_negative(row.get("edm_hi"))
                if lo is not None and hi is not None:
                    sigma_mu = float(max(abs(dm - lo), abs(hi - dm)))
                elif lo is not None:
                    sigma_mu = float(abs(dm - lo))
                elif hi is not None:
                    sigma_mu = float(abs(hi - dm))
                else:
                    sigma_mu = float("nan")
                d_mpc = luminosity_distance_mpc_from_distance_modulus(dm)
                sig_d = (
                    luminosity_distance_uncertainty_mpc_from_modulus(d_mpc, sigma_mu)
                    if np.isfinite(sigma_mu)
                    else float("nan")
                )
                return EddDistanceHit(
                    tier=1,
                    distance_modulus_mag=dm,
                    sigma_distance_modulus_mag=sigma_mu,
                    distance_mpc=d_mpc,
                    sigma_distance_mpc=sig_d,
                    source_table="edd_cmd_trgb",
                )

        if pgc in self._cf4.index:
            row = self._cf4.loc[pgc]
            dmcep = _finite_positive(row.get("DMcep"))
            ecep = _finite_non_negative(row.get("eDMcep"))
            if dmcep is not None:
                sigma_mu = ecep if ecep is not None else float("nan")
                d_mpc = luminosity_distance_mpc_from_distance_modulus(dmcep)
                sig_d = (
                    luminosity_distance_uncertainty_mpc_from_modulus(d_mpc, sigma_mu)
                    if sigma_mu is not None and np.isfinite(sigma_mu)
                    else float("nan")
                )
                return EddDistanceHit(
                    tier=2,
                    distance_modulus_mag=dmcep,
                    sigma_distance_modulus_mag=float(sigma_mu) if np.isfinite(sigma_mu) else float("nan"),
                    distance_mpc=d_mpc,
                    sigma_distance_mpc=sig_d,
                    source_table="edd_cf4_DMcep",
                )
            dmtrgb = _finite_positive(row.get("DMtrgb"))
            etrgb = _finite_non_negative(row.get("eDMt"))
            if dmtrgb is not None:
                sigma_mu = etrgb if etrgb is not None else float("nan")
                d_mpc = luminosity_distance_mpc_from_distance_modulus(dmtrgb)
                sig_d = (
                    luminosity_distance_uncertainty_mpc_from_modulus(d_mpc, sigma_mu)
                    if sigma_mu is not None and np.isfinite(sigma_mu)
                    else float("nan")
                )
                return EddDistanceHit(
                    tier=2,
                    distance_modulus_mag=dmtrgb,
                    sigma_distance_modulus_mag=float(sigma_mu) if np.isfinite(sigma_mu) else float("nan"),
                    distance_mpc=d_mpc,
                    sigma_distance_mpc=sig_d,
                    source_table="edd_cf4_DMtrgb",
                )
            dm = _finite_positive(row.get("DM"))
            edm = _finite_non_negative(row.get("eDM"))
            if dm is not None:
                sigma_mu = edm if edm is not None else float("nan")
                d_mpc = luminosity_distance_mpc_from_distance_modulus(dm)
                sig_d = (
                    luminosity_distance_uncertainty_mpc_from_modulus(d_mpc, sigma_mu)
                    if sigma_mu is not None and np.isfinite(sigma_mu)
                    else float("nan")
                )
                return EddDistanceHit(
                    tier=2,
                    distance_modulus_mag=dm,
                    sigma_distance_modulus_mag=float(sigma_mu) if np.isfinite(sigma_mu) else float("nan"),
                    distance_mpc=d_mpc,
                    sigma_distance_mpc=sig_d,
                    source_table="edd_cf4_DM",
                )

        if pgc in self._cf3.index:
            row = self._cf3.loc[pgc]
            dm = _finite_positive(row.get("DM"))
            edm = _finite_non_negative(row.get("eDM"))
            if dm is not None:
                sigma_mu = edm if edm is not None else float("nan")
                d_mpc = luminosity_distance_mpc_from_distance_modulus(dm)
                sig_d = (
                    luminosity_distance_uncertainty_mpc_from_modulus(d_mpc, sigma_mu)
                    if sigma_mu is not None and np.isfinite(sigma_mu)
                    else float("nan")
                )
                return EddDistanceHit(
                    tier=3,
                    distance_modulus_mag=dm,
                    sigma_distance_modulus_mag=float(sigma_mu) if np.isfinite(sigma_mu) else float("nan"),
                    distance_mpc=d_mpc,
                    sigma_distance_mpc=sig_d,
                    source_table="edd_cosmicflows3_DM",
                )
            dist_mpc = _finite_positive(row.get("Dist"))
            if dist_mpc is not None:
                dm_from_d = float(5.0 * np.log10(dist_mpc) + 25.0)
                sigma_mu = edm if edm is not None else float("nan")
                sig_d = (
                    luminosity_distance_uncertainty_mpc_from_modulus(dist_mpc, sigma_mu)
                    if sigma_mu is not None and np.isfinite(sigma_mu)
                    else float("nan")
                )
                return EddDistanceHit(
                    tier=3,
                    distance_modulus_mag=dm_from_d,
                    sigma_distance_modulus_mag=float(sigma_mu) if np.isfinite(sigma_mu) else float("nan"),
                    distance_mpc=dist_mpc,
                    sigma_distance_mpc=sig_d,
                    source_table="edd_cosmicflows3_Dist",
                )

        return None


def build_sparc_name_to_pgc_lookup(
    json_path: Path,
    cosmicflows_csv: Path,
    morphology_provider: SparcMorphologyProvider,
) -> dict[str, int]:
    """
    Build ``normalize_galaxy_name`` → ``pgc`` for EDD joins.

    JSON ``positions`` entries win over Cosmicflows-3 auto-ingest on key collision.
    """
    out: dict[str, int] = {}
    if json_path.is_file():
        try:
            raw = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("Could not read PGC override map %s: %s", json_path, exc)
        else:
            positions = raw.get("positions") if isinstance(raw, dict) else None
            if isinstance(positions, dict):
                for name, payload in positions.items():
                    if str(name).startswith("_") or not isinstance(payload, dict):
                        continue
                    try:
                        pgc_val = int(payload["pgc"])
                    except (KeyError, TypeError, ValueError):
                        continue
                    key = morphology_provider.normalize_galaxy_name(str(name))
                    out[key] = pgc_val

    if not cosmicflows_csv.is_file():
        return out

    try:
        cf3 = _read_cf3(cosmicflows_csv)
    except (OSError, ValueError) as exc:
        logger.warning("Could not ingest Cosmicflows-3 for PGC map %s: %s", cosmicflows_csv, exc)
        return out

    for _, row in cf3.iterrows():
        name_raw = row.get("Name")
        if name_raw is None or str(name_raw).strip() == "":
            continue
        key = morphology_provider.normalize_galaxy_name(str(name_raw).strip())
        if key in out:
            continue
        try:
            pgc_int = int(row["pgc"])
        except (TypeError, ValueError):
            continue
        out[key] = pgc_int
    return out


def apply_edd_distance_over_rotmod(
    deps: Any,
    *,
    morphology_key: str,
    galaxy_name: str,
    rotmod_distance_mpc: float,
    galaxy: Any,
) -> bool:
    """
    When EDD tiers 1–3 resolve a distance for the galaxy's ``pgc``, overwrite ``galaxy`` fields.

    Returns ``True`` when an EDD hit was applied (``prov_high_precision_dist`` metadata).
    """
    edd = getattr(deps, "edd", None)
    pmap: dict[str, int] = getattr(deps, "pgc_by_sparc_name", {}) or {}
    raw_key = galaxy_name.replace(" ", "").upper()
    pgc = pmap.get(morphology_key) or pmap.get(raw_key)
    galaxy.metadata["prov_high_precision_dist"] = False
    if edd is None or pgc is None:
        return False
    hit = edd.lookup_distance_for_pgc(int(pgc))
    if hit is None:
        return False
    galaxy.distance_mpc = float(hit.distance_mpc)
    if np.isfinite(hit.sigma_distance_mpc) and float(hit.sigma_distance_mpc) > 0.0:
        galaxy.distance_uncertainty_mpc = float(hit.sigma_distance_mpc)
    galaxy.metadata["edd_distance_tier"] = int(hit.tier)
    galaxy.metadata["edd_distance_source_table"] = str(hit.source_table)
    galaxy.metadata["prov_high_precision_dist"] = True
    _ = rotmod_distance_mpc  # SPARC header retained only when EDD misses
    return True


def build_edd_provenance_log_dict(
    deps: Any,
    *,
    morphology_key: str,
    galaxy_name: str,
    galaxy: Any,
) -> dict[str, Any]:
    """
    Snapshot for :mod:`src.reporting.galaxy_log` ``--logs`` rows (PGC join, tier, table, distances).

    Call after :func:`apply_edd_distance_over_rotmod` so ``galaxy.metadata`` and ``galaxy.distance_*``
    reflect any EDD overwrite.
    """
    pmap: dict[str, int] = getattr(deps, "pgc_by_sparc_name", {}) or {}
    raw_key = galaxy_name.replace(" ", "").upper()
    pgc = pmap.get(morphology_key) or pmap.get(raw_key)
    return {
        "edd_tables_loaded": getattr(deps, "edd", None) is not None,
        "pgc": int(pgc) if pgc is not None else None,
        "applied": bool(galaxy.metadata.get("prov_high_precision_dist")),
        "tier": galaxy.metadata.get("edd_distance_tier"),
        "source_table": galaxy.metadata.get("edd_distance_source_table"),
        "distance_mpc": galaxy.distance_mpc,
        "distance_uncertainty_mpc": galaxy.distance_uncertainty_mpc,
    }
