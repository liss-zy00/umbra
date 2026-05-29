from __future__ import annotations

"""
Planck HEALPix polarization reader.

Opens FITS with ``memmap=True`` and samples only the central pixel + neighbors to limit RAM use.

Provenance contract
-------------------
``initialize()`` and ``extract_stokes_qu_mean()`` always include a ``reference_method`` field so
that downstream code (and the final reports) can audit whether the calibration used a real
dynamic median of the map (``CALIBRATION_REFERENCE_METHOD_DYNAMIC``) or fell back to the
emergency-only constant (``CALIBRATION_REFERENCE_METHOD_FALLBACK``). The pipeline turns this
into a ``calibration_metadata`` block in every per-galaxy summary.
"""

from pathlib import Path
from typing import Any
import contextlib
import warnings

import numpy as np
from astropy.io import fits

from src.config.thresholds import (
    CALIBRATION_REFERENCE_METHOD_DYNAMIC,
    CALIBRATION_REFERENCE_METHOD_FALLBACK,
    PLANCK_CALIBRATION,
)


@contextlib.contextmanager
def _fits_open_memmap(path: Path):
    """
    Open FITS with ``memmap=True`` while silencing the common Astropy warning on Windows /
    synced folders (readonly mmap falls back to denywrite — data remain read-only in practice).
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*Could not memory map array with mode='readonly'.*",
        )
        with fits.open(path, memmap=True) as hdul:
            yield hdul


def resolve_stokes_amp_reference(stokes_summary: dict[str, Any], fits_summary: dict[str, Any]) -> float:
    """
    Return the Stokes amplitude reference ``S_ref`` actually used for B-field scaling.

    Preference order:
        1. ``stokes_summary['stokes_amp_reference']`` — echoed by ``extract_stokes_qu_mean`` from the
           live ``PlanckProvider`` state (always consistent with the per-pixel extraction).
        2. ``fits_summary['stokes_amp_reference']`` — value from ``PlanckProvider.initialize()``.
        3. ``PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY`` — last resort.

    This prevents the scenario grid from silently using the emergency constant when the FITS
    dict omits the key but the provider already computed the dynamic map median.
    """
    for src in (stokes_summary, fits_summary):
        v = src.get("stokes_amp_reference")
        if isinstance(v, (int, float)) and np.isfinite(float(v)) and float(v) > 0:
            return float(v)
    return float(PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY)


class PlanckProvider:
    """
    Provider responsible for loading and extracting Stokes Q/U parameters from Planck FITS maps.
    SRP: Acts as the telescope interface, ensuring lazy loading (memmap) prevents RAM exhaustion.
    """

    def __init__(self, fits_candidates: list[Path], verbose: bool = True) -> None:
        self.fits_candidates = fits_candidates
        self.verbose = verbose
        self.fits_path: Path | None = None
        self.nside: int | None = None
        self.healpix_nested: bool = False
        self._healpy = None
        # Provenance state set during ``initialize()`` and reused by ``extract_stokes_qu_mean()``.
        self.stokes_amp_reference: float = PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY
        self.reference_method: str = CALIBRATION_REFERENCE_METHOD_FALLBACK

    def _ensure_healpy(self) -> bool:
        """Lazy load healpy to prevent crashes if the environment is missing it."""
        if self._healpy is not None:
            return True
        try:
            import healpy as hp
            self._healpy = hp
            return True
        except ImportError:
            return False

    def _compute_stokes_amplitude_reference_median(self) -> tuple[float, str]:
        """
        Median of :math:`\\sqrt{Q^2+U^2}` over a strided subsample of the map (RAM-bounded).

        Returns
        -------
        tuple[float, str]
            ``(reference_amplitude, method_label)`` where ``method_label`` is
            ``CALIBRATION_REFERENCE_METHOD_DYNAMIC`` on success or
            ``CALIBRATION_REFERENCE_METHOD_FALLBACK`` whenever the map cannot yield a usable
            positive median (no FITS, read failure, all NaN, all zero, etc.).
        """
        pc = PLANCK_CALIBRATION
        emergency = pc.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY
        if not self.fits_path:
            return emergency, CALIBRATION_REFERENCE_METHOD_FALLBACK
        try:
            with _fits_open_memmap(self.fits_path) as hdul:
                data = hdul[1].data
                q_col = data["Q_STOKES"]
                u_col = data["U_STOKES"]
                n = int(len(q_col))
                stride = max(1, int(np.ceil(n / float(pc.STOKES_MEDIAN_MAX_SAMPLES))))
                q = np.asarray(q_col[::stride], dtype=np.float64)
                u = np.asarray(u_col[::stride], dtype=np.float64)
                valid = np.isfinite(q) & np.isfinite(u)
                amp = np.hypot(q[valid], u[valid])
                if amp.size == 0:
                    return emergency, CALIBRATION_REFERENCE_METHOD_FALLBACK
                med = float(np.median(amp))
                if not np.isfinite(med) or med < pc.STOKES_MEDIAN_MIN_POSITIVE:
                    return emergency, CALIBRATION_REFERENCE_METHOD_FALLBACK
                return med, CALIBRATION_REFERENCE_METHOD_DYNAMIC
        except Exception:
            return emergency, CALIBRATION_REFERENCE_METHOD_FALLBACK

    def initialize(self) -> dict[str, Any]:
        """
        Find the FITS file, read the header, check HEALPix ordering and compute the dynamic
        Stokes amplitude reference (median of |Q+iU|).

        Raises
        ------
        RuntimeError
            If a FITS map is found but ``healpy`` cannot be imported (required for HEALPix
            pixel indexing).

        The returned dict always carries ``stokes_amp_reference`` plus ``reference_method``
        so the rest of the pipeline never has to guess whether the calibration was real or
        a fallback emergency value.
        """
        for candidate in self.fits_candidates:
            if candidate.exists():
                self.fits_path = candidate
                break

        if not self.fits_path:
            self.stokes_amp_reference = PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY
            self.reference_method = CALIBRATION_REFERENCE_METHOD_FALLBACK
            return {
                "fits_found": False,
                "stokes_amp_reference": self.stokes_amp_reference,
                "reference_method": self.reference_method,
            }

        with _fits_open_memmap(self.fits_path) as hdul:
            hdr = hdul[1].header
            self.nside = int(hdr.get("NSIDE", 0)) or None
            ord_raw = hdr.get("ORDERING") or hdul[0].header.get("ORDERING") or "RING"
            self.healpix_nested = str(ord_raw).upper().strip() == "NESTED"

        if not self._ensure_healpy():
            raise RuntimeError(
                "healpy is not available in the current Python environment but is required to "
                "sample Planck HEALPix maps. Activate the dark-matter conda env (or invoke "
                "`conda run --no-capture-output -n dark-matter python src/main.py ...`)."
            ) from None

        reference_amp, reference_method = self._compute_stokes_amplitude_reference_median()
        self.stokes_amp_reference = float(reference_amp)
        self.reference_method = reference_method

        return {
            "fits_found": True,
            "fits_filename": self.fits_path.name,
            "nside": self.nside,
            "healpy_ready": self._healpy is not None,
            "stokes_amp_reference": self.stokes_amp_reference,
            "reference_method": self.reference_method,
        }

    def extract_stokes_qu_mean(
        self, right_ascension_degrees: float, declination_degrees: float, galaxy_name: str | None = None
    ) -> dict[str, Any]:
        """
        Extract the magnetic polarization data (Stokes Q and U) for the exact galaxy coordinate.

        LOGIC:
            1. Find the target pixel using HEALPix grid.
            2. Find its immediate neighbors.
            3. Read ONLY those pixels from the FITS file using Advanced Array Indexing.
            4. Measure the **mean angular separation** (radians) between the center pixel and the
               valid neighbors so the engine can convert it to a **physical** Δx using the
               galaxy distance — replacing the previous arbitrary scenario knob with a real
               sky-to-disk gradient length.

        Provenance: the result always echoes ``reference_method`` and ``stokes_amp_reference``
        carried by the provider since ``initialize()`` so error / skipped paths still report the
        calibration status of the run.
        """
        provenance = {
            "stokes_amp_reference": self.stokes_amp_reference,
            "reference_method": self.reference_method,
        }
        empty = {"status": "missing inputs", **provenance}
        if not self.fits_path or not self._ensure_healpy():
            return empty

        hp = self._healpy
        nside = int(self.nside)
        nest = self.healpix_nested

        # If ``initialize()`` could not set a dynamic reference (missing FITS, read error, etc.),
        # retry the map-wide median here whenever a FITS path exists — keeps ``stokes_amp_reference``
        # in the returned dict aligned with what the provider would use after a successful init.
        if self.fits_path and self.reference_method == CALIBRATION_REFERENCE_METHOD_FALLBACK:
            ref_amp, ref_method = self._compute_stokes_amplitude_reference_median()
            self.stokes_amp_reference = float(ref_amp)
            self.reference_method = ref_method

        # Convert RA/Dec to Spherical Theta/Phi (Healpy convention: Theta is co-latitude [0, pi]).
        theta = np.radians(90.0 - declination_degrees)
        phi = np.radians(right_ascension_degrees)

        center_pix = int(hp.ang2pix(nside, theta, phi, nest=nest))
        neighbors = hp.get_all_neighbours(nside, theta, phi, nest=nest)
        valid_nbr = neighbors[neighbors >= 0].astype(np.int64)

        # Mean great-circle separation between center pixel and valid neighbors.
        # Healpy's ``pix2ang`` returns (theta, phi) for each pixel; the angular separation between
        # two unit-sphere points is acos(sin θ₁ sin θ₂ cos(φ₁-φ₂) + cos θ₁ cos θ₂).
        if valid_nbr.size:
            theta_n, phi_n = hp.pix2ang(nside, valid_nbr, nest=nest)
            cos_sep = np.clip(
                np.sin(theta) * np.sin(theta_n) * np.cos(phi - phi_n) + np.cos(theta) * np.cos(theta_n),
                -1.0,
                1.0,
            )
            neighbour_separations_rad = np.arccos(cos_sep)
            mean_neighbor_separation_rad = float(np.mean(neighbour_separations_rad))
        else:
            theta_n = np.array([], dtype=np.float64)
            phi_n = np.array([], dtype=np.float64)
            neighbour_separations_rad = np.array([], dtype=np.float64)
            mean_neighbor_separation_rad = float("nan")

        # Array of indices we actually need to read from disk.
        pix_union = np.unique(np.concatenate(([np.int64(center_pix)], valid_nbr)))

        stokes_by_pixel: dict[int, tuple[float, float]] = {}

        try:
            with _fits_open_memmap(self.fits_path) as hdul:
                data = hdul[1].data
                n_pix = int(len(data))

                idx = pix_union[(pix_union >= 0) & (pix_union < n_pix)]

                # ADVANCED INDEXING: pulls only the necessary rows from the memmap array
                # without loading the other ~50 million pixels.
                q_raw = np.asarray(data["Q_STOKES"][idx], dtype=np.float64).ravel()
                u_raw = np.asarray(data["U_STOKES"][idx], dtype=np.float64).ravel()

                for i in range(idx.size):
                    p = int(idx[i])
                    stokes_by_pixel[p] = (q_raw[i], u_raw[i])

        except Exception as exc:
            return {"status": f"fits read error: {exc}", **provenance}

        q_center, u_center = stokes_by_pixel.get(center_pix, (0.0, 0.0))
        stokes_amp_center = float(np.hypot(q_center, u_center))

        neighbor_amps = [
            float(np.hypot(stokes_by_pixel[p][0], stokes_by_pixel[p][1]))
            for p in valid_nbr if p in stokes_by_pixel
        ]
        stokes_amp_neighbors_mean = float(np.mean(neighbor_amps)) if neighbor_amps else stokes_amp_center

        # Per-pixel ring: every neighbour entry carries its raw Q/U, its HEALPix (theta, phi)
        # angles, and the great-circle separation from the central pixel. The engine consumes
        # this list to fit a disk-plane planar gradient (see
        # :meth:`src.physics.engine.PhysicsEngine.calculate_magnetic_pressure_acceleration`).
        # The mean-amplitude / mean-separation keys above are retained for back-compat with
        # the scenario grid; they are derived statistics of the same ring.
        neighbour_ring: list[dict[str, float]] = []
        for i, p in enumerate(valid_nbr.tolist()):
            qu = stokes_by_pixel.get(int(p))
            if qu is None:
                continue
            neighbour_ring.append(
                {
                    "q": float(qu[0]),
                    "u": float(qu[1]),
                    "theta_rad": float(theta_n[i]),
                    "phi_rad": float(phi_n[i]),
                    "separation_rad": float(neighbour_separations_rad[i]),
                }
            )

        return {
            "stokes_amp_center": stokes_amp_center,
            "stokes_amp_neighbors_mean": stokes_amp_neighbors_mean,
            "q_center": float(q_center),
            "u_center": float(u_center),
            "mean_neighbor_separation_rad": mean_neighbor_separation_rad,
            "neighbour_ring": neighbour_ring,
            # Central pixel ``(theta, phi)`` in HEALPix convention (radians); convenient
            # for downstream code that needs to recompute great-circle bearings against
            # the ring without re-doing the RA/Dec → spherical conversion.
            "center_theta_rad": float(theta),
            "center_phi_rad": float(phi),
            "status": "ok",
            **provenance,
        }
