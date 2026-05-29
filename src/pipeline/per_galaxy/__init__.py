"""
Per-galaxy processing chain.

Each module here owns one well-defined responsibility in the pipeline that runs *for a
single galaxy*:

* :mod:`galaxy_processor` — top-level orchestrator for **Phase 3 (validation grid)**:
  SPARC → NED → Planck → physics grid → persistence + plotting.
* :mod:`discovery_processor` — **Phase 1 (Discovery census)**: a single engine pass per galaxy
  at the physical Planck pixel scale (``f_ion = 1.0``); rows feed ``discovery_*.csv``
  and the automatic pre-filter for grid entry via :func:`discovery_audits_clean`.
* :mod:`scenario_grid_runner` — runs the (Δx, f_ion) grid using an injected
  :class:`src.physics.PhysicsEngine`.
* :mod:`scenario_selector` — picks the "best" scenario per galaxy from the grid output.
* :mod:`metrics_builder` — assembles per-scenario metric dictionaries and the structured
  :data:`audit_flags` contract (single source of truth for reliability decisions).
* :mod:`morphology` — SPARC catalogue helpers (Hubble class, photometric LSB tagging).
"""
