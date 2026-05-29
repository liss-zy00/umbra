from __future__ import annotations

"""
NASA/IPAC Extragalactic Database (NED) metadata provider.

Resolves SPARC catalogue labels into ICRS coordinates, redshift and object type for the
rest of the pipeline (Planck HEALPix sampling needs RA/Dec; the engine needs the SPARC
distance, kept as a separate fallback). Name resolution is staged from the most
deterministic (static alias) to the most permissive (NED cone search around a SIMBAD
seed), so a noisy network never demotes the whole batch.

Resolution order (each step short-circuits on success)
------------------------------------------------------
1. **Exact alias map** (``ned_alias_map.json::alias_map``). One canonical query per SPARC
   label, hand-curated when NED's bare-name interpreter rejects an identifier.
2. **Pattern rules** (``ned_alias_map.json::pattern_rules``). First regex with
   ``re.fullmatch`` produces a templated query (``{name}`` substituted).
3. **Generic prefix permutation** (``ned_alias_map.json::ned_prefix_attempts``).
   ``f"{prefix}{name}"`` for each prefix; the empty prefix is always inserted at the
   front so the bare label is tried first.
4. **SIMBAD permutation** on the same identifier list.
5. **NED cone search** seeded by SIMBAD coordinates (when SIMBAD resolved RA/Dec) or by
   an explicit ``icrs_deg`` argument (e.g. from ``NED_ICRS_CONE_SEEDS_JSON``). The
   cone radius defaults to 2 arcminutes; the brightest catalogue match is selected.

Logging contract
----------------
Per-attempt failures (``RemoteServiceError``, network timeouts, empty tables) emit
``logger.debug`` only — the batch terminal stays clean. A single ``logger.warning`` is
issued once per galaxy when every avenue is exhausted (``raw_status == "unresolved"``).

The ``fetch_metadata`` contract never raises for typical query failures.
"""

import json
import logging
import re
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


# Default NED cone-search radius when seeding from a SIMBAD hit or an external ICRS
# supplement: 2 arcminutes is wide enough to catch SPARC-vs-NED astrometric offsets
# typical for low-surface-brightness disks while still narrow enough that the brightest
# match is, in practice, the target itself.
_DEFAULT_CONE_RADIUS_ARCMIN: float = 2.0

# Stable ``raw_status`` tokens from :func:`NedProvider.fetch_metadata`. The user-facing
# ``ned_status`` column (Discovery CSVs, ``validation_all.csv``) is produced by
# :func:`raw_status_to_ned_status` using **verbose** labels, e.g. ``OK_NED_VIA_STATIC_ALIAS_MAP``
# instead of a bare ``OK``, so CSV consumers see which resolution path succeeded.
STATUS_OK: str = "ok"
STATUS_OK_VIA_ALIAS_MAP: str = "ok_via_alias_map"
STATUS_OK_VIA_PATTERN: str = "ok_via_pattern_rule"
STATUS_OK_VIA_PREFIX: str = "ok_via_prefix_attempt"
STATUS_OK_VIA_SIMBAD: str = "ok_via_simbad"
STATUS_OK_VIA_NED_CONE: str = "ok_via_ned_cone"
STATUS_UNRESOLVED: str = "unresolved"
STATUS_ASTROQUERY_MISSING: str = "astroquery_missing"
STATUS_EMPTY_NAME: str = "empty_name"


def raw_status_to_ned_status(raw_status: str | None) -> str:
    """
    Map :func:`NedProvider.fetch_metadata` ``raw_status`` to the user-facing ``ned_status`` column.

    Tokens are verbose on purpose: CSV consumers see *how* coordinates were obtained
    (direct NED name hit vs static alias map vs pattern rule vs cone search, etc.).
    """
    if raw_status is None:
        return "NOT_FOUND_NED"
    if raw_status == STATUS_OK:
        return "OK_NED_SPARC_LABEL_DIRECT"
    if raw_status == STATUS_OK_VIA_ALIAS_MAP:
        return "OK_NED_VIA_STATIC_ALIAS_MAP"
    if raw_status == STATUS_OK_VIA_PATTERN:
        return "OK_NED_VIA_PATTERN_RULE"
    if raw_status == STATUS_OK_VIA_PREFIX:
        return "OK_NED_VIA_PREFIX_TRY_LIST"
    if raw_status == STATUS_OK_VIA_SIMBAD:
        return "OK_NED_VIA_SIMBAD_IDENTIFIER"
    if raw_status == STATUS_OK_VIA_NED_CONE:
        return "OK_NED_VIA_CONE_ARCMIN_SEARCH"
    if raw_status == STATUS_ASTROQUERY_MISSING:
        return "ASTROQUERY_MISSING"
    if raw_status == STATUS_EMPTY_NAME:
        return "ERROR_EMPTY_GALAXY_NAME"
    if raw_status.startswith("error:"):
        return "ERROR_NED_QUERY"
    return "NOT_FOUND_NED"


def _default_alias_map_path() -> Path:
    """Project-root ``sources/ned/ned_alias_map.json`` (same convention as other static tables)."""
    return Path(__file__).resolve().parent.parent.parent / "sources" / "ned" / "ned_alias_map.json"


def _load_ned_alias_config(path: Path | None) -> dict[str, Any]:
    """
    Read the JSON resolution config and normalise it into a small dict with three keys:
    ``alias_map`` (dict), ``pattern_rules`` (list of {regex, ned_query_template}) and
    ``ned_prefix_attempts`` (list of str). Missing sections degrade to empty containers.

    The historical schema (flat ``{SPARC -> NED_query}``) is still accepted to avoid
    breaking pre-existing alias files in the wild.
    """
    if path is None:
        path = _default_alias_map_path()
    if not path.is_file():
        return {"alias_map": {}, "pattern_rules": [], "ned_prefix_attempts": []}

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Could not read NED alias map %s: %s", path, exc)
        return {"alias_map": {}, "pattern_rules": [], "ned_prefix_attempts": []}

    alias_map: dict[str, str] = {}
    pattern_rules: list[dict[str, Any]] = []
    prefix_attempts: list[str] = []

    sectioned = isinstance(raw, dict) and ("alias_map" in raw or "pattern_rules" in raw or "ned_prefix_attempts" in raw)

    if sectioned:
        for k, v in (raw.get("alias_map") or {}).items():
            if str(k).startswith("_"):
                continue
            if isinstance(v, str) and v.strip():
                alias_map[str(k).strip()] = v.strip()
        for entry in raw.get("pattern_rules") or []:
            if not isinstance(entry, dict):
                continue
            regex = entry.get("regex")
            template = entry.get("ned_query_template")
            if not isinstance(regex, str) or not isinstance(template, str):
                continue
            try:
                compiled = re.compile(regex)
            except re.error as exc:
                logger.warning("Invalid regex %r in pattern_rules: %s", regex, exc)
                continue
            pattern_rules.append({"regex": compiled, "template": template})
        for entry in raw.get("ned_prefix_attempts") or []:
            if isinstance(entry, str):
                prefix_attempts.append(entry)
    else:
        # Legacy flat schema: ``{SPARC_label -> NED_query}``.
        for k, v in raw.items():
            if str(k).startswith("_"):
                continue
            if isinstance(v, str) and v.strip():
                alias_map[str(k).strip()] = v.strip()

    return {
        "alias_map": alias_map,
        "pattern_rules": pattern_rules,
        "ned_prefix_attempts": prefix_attempts,
    }


def _build_candidate_queries(
    name: str,
    *,
    alias_map: dict[str, str],
    pattern_rules: list[dict[str, Any]],
    prefix_attempts: list[str],
) -> list[tuple[str, str]]:
    """
    Build the ordered list of ``(query_string, provenance_tag)`` pairs the provider will
    attempt against NED/SIMBAD. ``provenance_tag`` is one of ``alias`` / ``pattern`` /
    ``bare`` / ``prefix`` so callers can attribute the eventual hit.

    Order
    -----
    ``alias`` first (highest signal), then ``pattern`` matches, then the bare name, then
    every configured prefix. Duplicates are removed while preserving the first occurrence.
    """
    candidates: list[tuple[str, str]] = []
    seen: set[str] = set()

    def _push(query: str, tag: str) -> None:
        q = query.strip()
        if not q or q in seen:
            return
        seen.add(q)
        candidates.append((q, tag))

    if name in alias_map:
        _push(alias_map[name], "alias")

    for rule in pattern_rules:
        if rule["regex"].fullmatch(name):
            _push(rule["template"].replace("{name}", name), "pattern")
            # Only the FIRST matching pattern produces an attempt — keeps the config
            # auditable; combine multiple variants via ``ned_prefix_attempts`` instead.
            break

    _push(name, "bare")

    # Empty prefix is implicit through the bare push above. Explicit empty entries are
    # tolerated but not duplicated thanks to ``seen``.
    for prefix in prefix_attempts:
        _push(f"{prefix}{name}", "prefix")

    return candidates


def _provenance_to_status(provenance: str) -> str:
    if provenance == "alias":
        return STATUS_OK_VIA_ALIAS_MAP
    if provenance == "pattern":
        return STATUS_OK_VIA_PATTERN
    if provenance == "prefix":
        return STATUS_OK_VIA_PREFIX
    return STATUS_OK


class NedProvider:
    """
    Thin adapter around :class:`astroquery.ipac.ned.Ned` with three layered fallbacks:
    a static alias map, regex + prefix permutations, SIMBAD, and finally a NED cone
    search seeded by SIMBAD (or by an explicit ICRS supplement).

    ``fetch_metadata`` never raises for typical query failures — it returns a dict with
    ``None`` placeholders so SPARC-derived distances can still drive kinematics. Every
    per-attempt failure is logged at ``DEBUG`` level so the batch terminal stays quiet;
    a single ``WARNING`` is emitted only when every avenue is exhausted.
    """

    def __init__(
        self,
        request_timeout_sec: float = 60.0,
        *,
        alias_map_path: Path | None = None,
        cone_radius_arcmin: float = _DEFAULT_CONE_RADIUS_ARCMIN,
    ) -> None:
        self.request_timeout_sec = request_timeout_sec
        self.cone_radius_arcmin = float(cone_radius_arcmin)

        config = _load_ned_alias_config(alias_map_path)
        self._alias_map: dict[str, str] = config["alias_map"]
        self._pattern_rules: list[dict[str, Any]] = config["pattern_rules"]
        self._prefix_attempts: list[str] = config["ned_prefix_attempts"]

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_candidates(self, name: str) -> list[tuple[str, str]]:
        return _build_candidate_queries(
            name,
            alias_map=self._alias_map,
            pattern_rules=self._pattern_rules,
            prefix_attempts=self._prefix_attempts,
        )

    @staticmethod
    def _extract_redshift(table: Any) -> float | None:
        for cand in ("Redshift", "z", "Z"):
            if cand not in table.colnames:
                continue
            zval = table[cand][0]
            if zval is None or str(zval).strip() in ("", "--"):
                return None
            try:
                return float(zval)
            except (TypeError, ValueError):
                return None
        return None

    @staticmethod
    def _extract_object_type(table: Any) -> str | None:
        for cand in ("Type", "Object Type", "OBJTYPE"):
            if cand in table.colnames:
                try:
                    return str(table[cand][0])
                except (IndexError, TypeError, ValueError):
                    return None
        return None

    def _query_ned_object(self, query: str) -> dict[str, Any] | None:
        """One ``Ned.query_object`` attempt; returns parsed dict or ``None``."""
        try:
            from astroquery.ipac.ned import Ned
        except ImportError:
            return None
        try:
            Ned.TIMEOUT = self.request_timeout_sec
            table = Ned.query_object(query)
        except Exception as exc:  # noqa: BLE001 — third-party query boundary
            logger.debug("NED query_object failed for %r: %s", query, exc)
            return None
        if table is None or len(table) == 0:
            logger.debug("NED query_object returned empty table for %r", query)
            return None
        try:
            return {
                "right_ascension_degrees": float(table["RA"][0]),
                "declination_degrees": float(table["DEC"][0]),
                "redshift": self._extract_redshift(table),
                "object_type": self._extract_object_type(table),
            }
        except (KeyError, IndexError, TypeError, ValueError) as exc:
            logger.debug("Could not parse NED table for %r: %s", query, exc)
            return None

    def _query_simbad_object(self, query: str) -> dict[str, Any] | None:
        try:
            from astroquery.simbad import Simbad
        except ImportError:
            return None
        try:
            table = Simbad.query_object(query)
        except Exception as exc:  # noqa: BLE001
            logger.debug("SIMBAD query_object failed for %r: %s", query, exc)
            return None
        if table is None or len(table) == 0:
            return None
        try:
            ra = float(table["ra"][0])
            dec = float(table["dec"][0])
        except (KeyError, IndexError, TypeError, ValueError):
            return None
        z_out: float | None = None
        for zcol in ("rvz_redshift", "redshift", "z_value"):
            if zcol not in table.colnames:
                continue
            zval = table[zcol][0]
            if zval is None or str(zval).strip() in ("", "--"):
                continue
            try:
                zf = float(zval)
                if zf == zf:  # finite check without numpy
                    z_out = zf
                    break
            except (TypeError, ValueError):
                continue
        obj_type: str | None = None
        if "main_id" in table.colnames:
            try:
                obj_type = str(table["main_id"][0])
            except (IndexError, TypeError, ValueError):
                obj_type = None
        return {
            "right_ascension_degrees": ra,
            "declination_degrees": dec,
            "redshift": z_out,
            "object_type": obj_type,
        }

    def _query_ned_cone(self, ra_deg: float, dec_deg: float) -> dict[str, Any] | None:
        """``Ned.query_region`` cone search; returns the closest match, if any."""
        try:
            from astropy import units as u
            from astropy.coordinates import SkyCoord
            from astroquery.ipac.ned import Ned
        except ImportError:
            return None
        try:
            coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
            Ned.TIMEOUT = self.request_timeout_sec
            table = Ned.query_region(coord, radius=self.cone_radius_arcmin * u.arcmin)
        except Exception as exc:  # noqa: BLE001
            logger.debug(
                "NED query_region failed at (RA=%.5f, Dec=%.5f, r=%.2f'): %s",
                ra_deg, dec_deg, self.cone_radius_arcmin, exc,
            )
            return None
        if table is None or len(table) == 0:
            return None
        try:
            return {
                "right_ascension_degrees": float(table["RA"][0]),
                "declination_degrees": float(table["DEC"][0]),
                "redshift": self._extract_redshift(table),
                "object_type": self._extract_object_type(table),
            }
        except (KeyError, IndexError, TypeError, ValueError) as exc:
            logger.debug("Could not parse NED cone table: %s", exc)
            return None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fetch_metadata(
        self,
        galaxy_name: str,
        *,
        icrs_deg: tuple[float, float] | None = None,
    ) -> dict[str, Any]:
        """
        Resolve ``galaxy_name`` to NED/SIMBAD metadata using the layered fallback chain
        documented at the module level.

        Parameters
        ----------
        galaxy_name
            SPARC label as it appears in the rotmod filename.
        icrs_deg
            Optional ``(ra_deg, dec_deg)`` seed used to trigger a NED cone search when
            every name-based lookup fails. Typically wired from
            ``NED_ICRS_CONE_SEEDS_JSON`` (:mod:`src.pipeline.constants`).

        Returns
        -------
        dict
            Flat metadata with the keys
            ``right_ascension_degrees``, ``declination_degrees``, ``redshift``,
            ``distance_mpc`` (always ``None`` here — the SPARC header is the canonical
            distance source), ``distance_source`` (``"NED"`` only when NED itself
            answered, otherwise ``None``), ``object_type`` and ``raw_status``.
        """
        empty: dict[str, Any] = {
            "right_ascension_degrees": None,
            "declination_degrees": None,
            "redshift": None,
            "distance_mpc": None,
            "distance_source": None,
            "object_type": None,
            "raw_status": "not_attempted",
        }

        name = str(galaxy_name).strip()
        if not name:
            empty["raw_status"] = STATUS_EMPTY_NAME
            return empty

        try:
            from astroquery.ipac.ned import Ned  # noqa: F401 — import probe
        except ImportError as exc:
            logger.warning("astroquery not available for NED: %s", exc)
            empty["raw_status"] = STATUS_ASTROQUERY_MISSING
            return empty

        candidates = self._build_candidates(name)

        # ----- Stage 1: NED.query_object permutation -----
        for query, provenance in candidates:
            hit = self._query_ned_object(query)
            if hit is None:
                continue
            empty.update(hit)
            empty["distance_source"] = "NED"
            empty["raw_status"] = _provenance_to_status(provenance)
            return empty

        # ----- Stage 2: SIMBAD permutation (same identifier list) -----
        simbad_hit: dict[str, Any] | None = None
        for query, _ in candidates:
            simbad_hit = self._query_simbad_object(query)
            if simbad_hit is not None:
                break

        # ----- Stage 3: NED cone search seeded by SIMBAD or by ``icrs_deg`` -----
        seed_ra: float | None = None
        seed_dec: float | None = None
        if simbad_hit is not None:
            seed_ra = simbad_hit["right_ascension_degrees"]
            seed_dec = simbad_hit["declination_degrees"]
        elif icrs_deg is not None:
            try:
                seed_ra = float(icrs_deg[0])
                seed_dec = float(icrs_deg[1])
            except (TypeError, ValueError, IndexError):
                seed_ra = seed_dec = None

        if seed_ra is not None and seed_dec is not None:
            cone_hit = self._query_ned_cone(seed_ra, seed_dec)
            if cone_hit is not None:
                empty.update(cone_hit)
                empty["distance_source"] = "NED"
                empty["raw_status"] = STATUS_OK_VIA_NED_CONE
                return empty

        if simbad_hit is not None:
            empty.update(simbad_hit)
            empty["distance_source"] = None
            empty["raw_status"] = STATUS_OK_VIA_SIMBAD
            return empty

        logger.warning(
            "Could not resolve %r in NED nor SIMBAD after %d name candidate(s); marking unresolved.",
            name,
            len(candidates),
        )
        empty["raw_status"] = STATUS_UNRESOLVED
        return empty
