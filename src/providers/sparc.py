from __future__ import annotations

"""
SPARC rotation-curve ingestion (Lelli+2016 ``*_rotmod.dat`` and the morphology catalog).

The module exposes two providers with strict single-responsibility contracts:

* :class:`SparcProvider` parses ``*_rotmod.dat`` files and produces a per-radius
  :class:`pandas.DataFrame` with the canonical SPARC columns
  ``radius_kpc, observed_velocity_kms, gas_velocity_kms, disk_velocity_kms,
  bulge_velocity_kms`` plus the synthesized
  ``baryonic_velocity_kms = sqrt(V_gas² + V_disk² + V_bulge²)`` (the rotational analogue
  of the *baryonic* velocity used in the RAR/BTFR literature). All values are reported in
  km/s and kpc — SI conversions live in :mod:`src.physics.kinematics`.

* :class:`SparcMorphologyProvider` parses the JSON catalog produced from
  ``Lelli2016`` and exposes :py:meth:`load_morphology_lookup` mapping galaxy name → bin.
  This is also the source for the photometric-LSB pseudo-class flag used by the
  population-level morphology comparison.

All parsing is local-disk only (no network I/O); foreign-key lookups (NED, Planck) are
delegated to :mod:`src.providers.ned` and :mod:`src.providers.planck` respectively.

Inclination uncertainty propagation
-----------------------------------
SPARC distributes ``V_obs = V_LOS / sin i`` already deprojected, so a residual
inclination uncertainty ``δi`` becomes the dominant kinematic systematic. The provider
exposes :py:meth:`SparcProvider.attach_inclination_uncertainty` which writes a column
``observed_velocity_systematic_kms = V_obs · |cos i / sin i| · δi`` consumed by the
:class:`BaryonicGapCalculator` to populate ``centripetal_acceleration_systematic``.
"""

import json
import re
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

DISTANCE_PATTERN = re.compile(r"Distance\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*Mpc", re.IGNORECASE)
FLOAT_TOKEN_PATTERN = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def _parse_float_tokens(line: str, need_at_least: int = 6) -> list[float] | None:
    """Helper to extract raw float numbers from SPARC text lines."""
    toks = FLOAT_TOKEN_PATTERN.findall(line)
    if len(toks) < need_at_least:
        return None
    try:
        return [float(t) for t in toks[:need_at_least]]
    except ValueError:
        return None

class SparcProvider:
    """
    Provider responsible for loading SPARC rotation curves (*_rotmod.dat).
    SRP (Single Responsibility Principle): Handles exclusively the extraction of kinematics.
    """

    def load_rotation_curve(self, file_path: str | Path) -> tuple[pd.DataFrame, float]:
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"SPARC file not found: {path}")

        distance_mpc: float | None = None
        rows: list[list[float]] = []
        
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                
                # Extract distance from the header
                if line.startswith("#"):
                    if distance_mpc is None:
                        match = DISTANCE_PATTERN.search(line)
                        if match:
                            distance_mpc = float(match.group(1))
                    continue

                parts = _parse_float_tokens(line, 6)
                if parts is None:
                    continue

                radius_kpc, observed_velocity_kms, _errv, gas_velocity_kms, disk_velocity_kms, bulge_velocity_kms = parts
                
                # Filter out invalid rows (negative radius or NaNs)
                if not all(np.isfinite(p) for p in parts) or radius_kpc <= 0:
                    continue
                
                rows.append([radius_kpc, observed_velocity_kms, gas_velocity_kms, disk_velocity_kms, bulge_velocity_kms])

        if distance_mpc is None:
            raise ValueError(f"Distance not found in header: {path}")
        if not rows:
            raise ValueError(f"No valid rows in SPARC file: {path}")

        df = pd.DataFrame(
            rows,
            columns=["radius_kpc", "observed_velocity_kms", "gas_velocity_kms", "disk_velocity_kms", "bulge_velocity_kms"],
        )

        # PHYSICS EXPLANATION: Baryonic Velocity (Newtonian gravity linear addition)
        # Gravitational potentials (Phi) add up linearly: Phi_total = Phi_gas + Phi_disk + Phi_bulge.
        # Since Centripetal Force depends on V^2 / R, and F = -d(Phi)/dR,
        # velocities add up in quadrature: V_total^2 = V_gas^2 + V_disk^2 + V_bulge^2
        df["baryonic_velocity_kms"] = np.sqrt(
            df["gas_velocity_kms"] ** 2 + df["disk_velocity_kms"] ** 2 + df["bulge_velocity_kms"] ** 2
        )
        # Per-shell stellar fraction of the baryonic V². Used by GasDensityEstimator and the
        # reliability gate to flag shells where stars (not gas) dominate the rotational support
        # — the spherical-shell ``M ≈ V_gas² r / G`` density estimator becomes unreliable there.
        v_bar_sq = df["baryonic_velocity_kms"] ** 2
        v_stellar_sq = df["disk_velocity_kms"] ** 2 + df["bulge_velocity_kms"] ** 2
        df["stellar_velocity_fraction"] = np.where(
            v_bar_sq > 0,
            v_stellar_sq / v_bar_sq,
            0.0,
        )
        return df, distance_mpc

    @staticmethod
    def attach_inclination_uncertainty(
        df: pd.DataFrame,
        inclination_deg: float | None,
        inclination_uncertainty_deg: float | None,
    ) -> pd.DataFrame:
        """
        Attach per-shell systematic uncertainty bands on ``observed_velocity_kms`` propagated from
        the inclination uncertainty of the SPARC rotmod fit.

        SPARC rotmod files publish ``V_obs`` already deprojected with ``V_circ = V_los / sin(i)``.
        A small inclination error ``δi`` therefore induces a fractional velocity error of:

            δV_obs / V_obs = -cos(i) · δi / sin(i),

        which is the dominant kinematic systematic on outer-disk dynamics. This method writes the
        absolute systematic envelope ``observed_velocity_systematic_kms`` so downstream engines can
        report a band on ``acceleration_gap`` without having to re-implement the calculation.

        If inclination metadata is missing, the systematic column is zero — the function never
        crashes the pipeline.
        """
        out = df.copy()
        v_obs = out["observed_velocity_kms"].to_numpy(dtype=float)
        if (
            inclination_deg is None
            or inclination_uncertainty_deg is None
            or not np.isfinite(float(inclination_deg))
            or not np.isfinite(float(inclination_uncertainty_deg))
        ):
            out["observed_velocity_systematic_kms"] = np.zeros_like(v_obs)
            return out
        i_rad = np.radians(float(inclination_deg))
        di_rad = np.radians(float(inclination_uncertainty_deg))
        sin_i = np.sin(i_rad)
        # Avoid division by zero for unphysical i = 0; report inf so the reliability gate kicks in.
        if not np.isfinite(sin_i) or abs(sin_i) <= 0.0:
            out["observed_velocity_systematic_kms"] = np.full_like(v_obs, np.inf)
            return out
        relative_error = abs(np.cos(i_rad) / sin_i) * abs(di_rad)
        out["observed_velocity_systematic_kms"] = np.abs(v_obs) * float(relative_error)
        return out


def _load_pa_supplement(path: str | Path) -> dict[str, dict[str, Any]]:
    """
    Read the disk-PA overlay JSON (default ``MORPHOLOGY_DISK_PA_OVERLAY_JSON`` in :mod:`src.pipeline.constants`)
    and return a
    normalised ``{galaxy_name -> {position_angle_deg, position_angle_uncertainty_deg, source}}`` map.

    SPARC does not publish a disk position angle; values are curated externally (HyperLEDA,
    NED, etc.). Missing file or empty ``positions`` section returns ``{}`` so the engine
    falls back to the axisymmetric path and emits ``warn_no_pa`` per galaxy.
    """
    p = Path(path)
    if not p.exists():
        return {}
    try:
        payload = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    positions = payload.get("positions") if isinstance(payload, dict) else None
    if not isinstance(positions, dict):
        return {}

    out: dict[str, dict[str, Any]] = {}
    for raw_name, value in positions.items():
        if not isinstance(value, dict):
            continue
        try:
            position_angle_deg_parsed = float(value["position_angle_deg"])
        except (KeyError, TypeError, ValueError):
            continue
        if not np.isfinite(position_angle_deg_parsed):
            continue
        try:
            position_angle_uncertainty_deg_parsed = float(
                value.get("position_angle_uncertainty_deg", float("nan"))
            )
        except (TypeError, ValueError):
            position_angle_uncertainty_deg_parsed = float("nan")
        if not np.isfinite(position_angle_uncertainty_deg_parsed):
            position_angle_uncertainty_deg_parsed = float("nan")
        source = str(value.get("source") or "unknown")
        # ``SparcMorphologyProvider.normalize_galaxy_name`` is the canonical SPARC key style;
        # we apply it here so a single lookup table serves the morphology + PA join.
        norm = SparcMorphologyProvider.normalize_galaxy_name(str(raw_name))
        out[norm] = {
            "position_angle_deg": position_angle_deg_parsed,
            "position_angle_uncertainty_deg": position_angle_uncertainty_deg_parsed,
            "source": source,
        }
        raw_key = str(raw_name).replace(" ", "").upper()
        if raw_key != norm:
            out[raw_key] = out[norm]
    return out


class SparcMorphologyProvider:
    """
    Provider responsible for classifying galaxy types using the Lelli 2016 JSON catalog
    and merging the optional multi-archive disk-PA overlay (``MORPHOLOGY_DISK_PA_OVERLAY_JSON``).

    SRP: Parses the JSON to provide the morphology class mapping and attaches overlay PA.
    """

    def __init__(
        self,
        json_path: str | Path,
        *,
        pa_supplement_path: str | Path | None = None,
    ) -> None:
        self.json_path = Path(json_path).resolve()
        # Optional PA overlay (multi-catalog). Default: sibling ``disk_position_angle_overlay.json`` next
        # to the Lelli JSON (used by unit tests); production wiring passes ``MORPHOLOGY_DISK_PA_OVERLAY_JSON``.
        pa_path = (
            Path(pa_supplement_path).resolve()
            if pa_supplement_path is not None
            else self.json_path.parent / "disk_position_angle_overlay.json"
        )
        self._pa_supplement: dict[str, dict[str, Any]] = _load_pa_supplement(pa_path)

    @staticmethod
    def normalize_galaxy_name(name: str) -> str:
        """
        Strips spaces, makes it uppercase, and removes spurious leading zeros.
        Crucial for joining SPARC file names (e.g., 'NGC0024') with JSON keys ('NGC 24').
        """
        s = re.sub(r"\s+", "", name.strip().upper())
        m = re.match(r"^([A-Za-z]+)0*(\d+)(.*)$", s)
        if m:
            prefix, num, rest = m.groups()
            return f"{prefix}{int(num)}{rest}"
        return s

    @staticmethod
    def classify_morphology_from_t(hubble_type: int | float | None) -> str:
        """
        Maps the continuous SPARC Hubble stage T into standard taxonomic strings.
        This provides the final identity for the galaxy in reports.
        """
        if hubble_type is None:
            return "Unclassified"
        try:
            t = int(round(float(hubble_type)))
        except (TypeError, ValueError):
            return "Unclassified"
            
        if t <= 0: return "Lenticular/S0"
        if 1 <= t <= 2: return "Early-type Spirals (Sa, Sab)"
        if 3 <= t <= 4: return "Intermediate Spirals (Sb, Sbc)"
        if 5 <= t <= 6: return "Late-type Spirals (Sc, Scd)"
        if 7 <= t <= 9: return "Sm/Sd Spirals"
        if t == 10: return "Dwarf Irregulars (Im)"
        if t == 11: return "Blue Compact Dwarfs (BCD)"
        return "Unclassified"

    def load_morphology_lookup(self) -> dict[str, dict[str, Any]]:
        """
        Read the Lelli+2016 catalog JSON and return a dictionary keyed by normalized galaxy name.

        In addition to the morphology class, the entry now exposes:
          * ``inclination_deg`` and ``inclination_uncertainty_deg`` — used downstream to propagate
            a systematic uncertainty band on V_obs (SPARC rotmod files are already deprojected
            using this inclination, so its error is the dominant kinematic systematic).
          * ``vflat_kms`` and ``vflat_uncertainty_kms`` — fitted flat rotation velocity from the
            catalog, used by the BTFR / RAR diagnostics.
          * ``hi_mass_1e9_solar`` — atomic gas mass; useful for sanity-checking the rotmod-derived
            gas density against the catalog HI mass.
        """
        if not self.json_path.exists():
            return {}

        with self.json_path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)

        rows = payload.get("data", [])
        lookup: dict[str, dict[str, Any]] = {}

        for row in rows:
            galaxy_raw = row.get("Galaxy")
            if not galaxy_raw:
                continue

            galaxy_name = str(galaxy_raw).strip()
            hubble_type = row.get("T")

            # Map both normalized ("NGC24") and raw ("NGC0024") to prevent lookup misses.
            key_primary = self.normalize_galaxy_name(galaxy_name)
            key_raw = galaxy_name.replace(" ", "").upper()

            # Inject the PA supplement entry (HyperLEDA-curated) when present; absent keys
            # leave both PA fields ``None`` so the engine knows to use the axisymmetric path.
            pa_entry = (
                self._pa_supplement.get(key_primary)
                or self._pa_supplement.get(key_raw)
                or {}
            )

            entry = {
                "hubble_type": hubble_type,
                "morphology_class": self.classify_morphology_from_t(hubble_type),
                # Surface brightness (solar luminosities per pc²) — used by pipeline LSB proxy filters.
                "disk_surface_brightness_solar_per_pc2": row.get("SBdisk_solLum_pc2"),
                # SPARC luminosity distance (Mpc) and its 1-σ catalogue error. The error feeds
                # the Monte-Carlo band in :mod:`src.physics.diagnostics.uncertainty`.
                "distance_mpc": row.get("D_Mpc"),
                "distance_uncertainty_mpc": row.get("e_D_Mpc"),
                # Disk inclination (deg) and its error — needed for systematic V_obs propagation.
                "inclination_deg": row.get("Inc_deg"),
                "inclination_uncertainty_deg": row.get("e_Inc_deg"),
                # Catalogued flat rotation velocity (km/s) and its error — used by RAR / BTFR.
                "vflat_kms": row.get("Vflat_km_s"),
                "vflat_uncertainty_kms": row.get("e_Vflat_km_s"),
                # Atomic-gas mass (10⁹ M_⊙) for cross-checking the rotmod-derived density and
                # for the BTFR baryonic-mass reconstruction.
                "hi_mass_1e9_solar": row.get("MHI_1e9solMass"),
                # Spitzer 3.6 µm luminosity (10⁹ L_⊙) — feeds the BTFR stellar-mass branch.
                "luminosity_3p6um_1e9_solar": row.get("L36_1e9solLum"),
                # Disk major-axis position angle (deg) from the curated supplement. ``None``
                # ⇒ no PA available ⇒ engine falls back to axisymmetric gradient + warn_no_pa.
                "position_angle_deg": pa_entry.get("position_angle_deg"),
                "position_angle_uncertainty_deg": pa_entry.get("position_angle_uncertainty_deg"),
                "position_angle_source": pa_entry.get("source"),
            }

            lookup[key_primary] = entry
            if key_raw != key_primary:
                lookup[key_raw] = entry

        return lookup