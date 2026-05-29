# A Falsifiable Magnetohydrodynamic Test for Galactic Rotation-Curve Anomalies: The E-G Model

> **Status: working draft.** Structural template for a future submission to a
> peer-reviewed venue (e.g. MNRAS or A&A). The numerical results that populate
> §Methodology and any future Results section are produced by the open-source
> pipeline in this repository; they are not embedded here so that the draft
> stays decoupled from any individual execution. The Portuguese version is at
> [`Draft_Paper_EG_Model_PT.md`](Draft_Paper_EG_Model_PT.md).

---

## Abstract

The flat rotation profiles of spiral and dwarf galaxies are conventionally
interpreted as evidence of a non-baryonic, collisionless dark-matter halo
that provides the missing centripetal acceleration. We test an alternative
contribution to the radial force budget: a magnetohydrodynamic
(MHD) magnetic-pressure gradient acting on the warm ionised component of
the interstellar medium (WIM). For each galaxy in the SPARC catalogue we
(i) compute the kinematic acceleration deficit
`Δa(r) = V_obs²(r)/r − V_bar²(r)/r`; (ii) invert the gas-only rotation
component into a volumetric ionised-gas density via a cylindrical-shell
decomposition; (iii) calibrate the Planck linear-polarisation Stokes amplitude
at the galaxy's sky position into a magnetic-field strength via a
CMB-energy-density dimensional anchor; and (iv) evaluate the magnetic-pressure
acceleration `a_B(r) = −∇(B²/2μ₀)/ρ_ion(r)` shell by shell. We report the
outer-disk modulus ratio `|a_B|/|Δa|` and the Pearson correlation between the
radial profiles of `a_B` and `Δa` over a scenario grid spanning three gradient
scales (50, 100, 500 pc) and three ionisation fractions (0.1, 1, 10 %). The
pipeline is exercised against synthetic NFW negative controls and reports an
explicit reliability audit per galaxy and per scenario. We frame the result
as a falsification test of the MHD-pressure ansatz, not as a refutation of
dark matter as a phenomenon.

**Keywords:** galaxies: kinematics and dynamics — galaxies: magnetic fields —
magnetohydrodynamics (MHD) — ISM: magnetic fields — methods: data analysis.

---

## 1. Introduction

Galactic rotation curves of spiral and dwarf galaxies exhibit a robust
empirical anomaly: at radii where the visible (baryonic) mass distribution
predicts a Keplerian fall-off `V ∝ r^{−1/2}`, spectroscopic measurements
consistently show flat velocity profiles (e.g. Rubin & Ford 1970; Bosma 1981).
Within the standard cosmological framework (ΛCDM), this anomaly is closed by
a spherical, collisionless dark-matter halo whose gravity supplies the
missing centripetal acceleration.

Galactic magnetic fields with characteristic strengths of order one to ten
microgauss are observationally established in spiral and dwarf galaxies
through diffuse synchrotron emission, Faraday rotation of background
polarised sources, and the linear-polarisation maps of the Planck satellite
(Planck Collaboration 2016). The dynamical contribution of these fields to
the rotation-curve problem is, however, typically treated as a second-order
correction to a gravity-only description.

The E-G (Electrodynamic–Gravity) framework presented here is **not** a
proposal of modified gravity. It is a quantitative test of a single
falsifiable claim: that the magnetic-pressure gradient
`−∇(B²/2μ₀)`, acting on the warm ionised component of the interstellar
medium, contributes a radial acceleration of order
`Δa(r) = V_obs²(r)/r − V_bar²(r)/r` at the radii where dark-matter haloes
are conventionally invoked. The test produces order-of-magnitude estimates
of `a_B(r)` from independent observational inputs (the SPARC kinematics on
one side and the Planck polarisation maps on the other) and reports the
outcome with the same rigour for confirmation and falsification.

The paper is organised as follows. Section 2 states the hypothesis as a
two-regime decomposition of the radial force budget. Section 3 derives the
magnetic-pressure-gradient ansatz from the volume Lorentz force and discusses
the coupling to the WIM. Section 4 describes the data, the calibration
scheme and the per-shell estimator pipeline. Section 5 discusses the
falsification protocol and the explicit reliability audit. Section 6
enumerates the known limitations and the planned refinements.

---

## 2. Hypothesis — Two-Regime Decomposition

The total radial acceleration at galactocentric radius `r` is decomposed as

```text
a_total(r) = a_grav(r) + a_B(r),
```

where `a_grav(r) = V_bar²(r)/r` is the Newtonian baryonic contribution from
the SPARC mass decomposition `V_bar = √(V_gas² + V_disk² + V_bulge²)` and
`a_B(r)` is the MHD magnetic-pressure contribution defined in Section 3.

Two regimes are distinguished by the relative weight of the two terms:

- **Regime A (compact/baryonic).** In the inner few kpc of a galactic disk
  and in stellar systems, `a_grav` dominates: `|a_grav| ≫ |a_B|`. Newtonian
  gravity is recovered as the local limit and the E-G framework reduces
  trivially to it.
- **Regime B (extended/magnetised).** In the outer disk, where the kinematic
  deficit `Δa(r) = V_obs²(r)/r − V_bar²(r)/r` opens up, the question is
  whether `a_B(r)` is of order `Δa(r)`. A positive answer would identify
  the MHD pressure gradient as a viable contributor to the missing
  acceleration; a negative answer (in particular `|a_B| ≪ |Δa|` by orders
  of magnitude) would falsify this particular ansatz.

The hypothesis does **not** assert that the entirety of the rotation-curve
anomaly is electromagnetic in origin; it tests the magnitude of one specific
contribution that the conventional analysis routinely neglects.

---

## 3. Proposed Mechanism — MHD Pressure Gradient on the WIM

### 3.1 From microscopic Lorentz force to magnetic pressure

A single charged particle of charge `q` and velocity `v` in an
electromagnetic field `(E, B)` experiences the Lorentz force
`F = q(E + v × B)`. For an interstellar plasma, the relevant timescales
are long compared to the gyration period and the dynamical quantity of
interest is the bulk acceleration of a fluid element. Summing over the
charge content of a fluid parcel and assuming quasi-neutrality, the
volume Lorentz force per unit volume reduces, in the non-relativistic MHD
limit, to `f = J × B = (∇ × B) × B / μ₀`. A standard vector identity
decomposes this into a pressure-gradient term `−∇(B²/2μ₀)` and a magnetic
tension term `(B · ∇)B/μ₀` (Jackson 1999).

### 3.2 The ideal-MHD pressure-gradient ansatz

This paper retains only the pressure-gradient component and treats the
plasma as ideal (high conductivity, no resistivity, no anisotropic pressure
tensor, no cosmic-ray contribution, no turbulent pressure). The acceleration
that the ionised component of the interstellar medium experiences from this
gradient is

```text
a_B(x) ≈ − ∇u_B(x) / ρ_ion(x),     with u_B(x) = B(x)² / (2 μ₀).
```

The justification for dropping the tension term is discussed as a refinement
axis in §6.

### 3.3 The WIM as the responsive component

The mass density entering the denominator is **not** the total gas density.
It is the density of the component that responds electromagnetically on the
relevant timescale, modelled as

```text
ρ_ion(r) = f_ion · ρ_gas(r),
```

with `f_ion ∈ (0, 1]` a free ionisation fraction. The observationally
relevant component is the Warm Ionised Medium (WIM): ionised hydrogen at
`T ~ 10⁴ K` permeating galactic disks and haloes with characteristic mass
densities orders of magnitude below the neutral HI gas detected by 21-cm
surveys (Reynolds 1990; Haffner et al. 2009).

### 3.4 Scaling geometry

Gravitational acceleration from a centrally concentrated baryonic mass
decays as `1/r²`. The magnetic-pressure gradient in a filamentary or
disk geometry decays substantially slower (`B ∝ 1/r` along an idealised
straight current-carrying filament, with the corresponding pressure gradient
non-trivial over inter-filament distances). The competition between the two
contributions therefore depends explicitly on the radial geometry, and the
crossover radius is in principle a derived quantity of the model rather than
a free parameter.

---

## 4. Methodology

> **Code-of-truth note.** The class names referenced below are the
> Single-Responsibility estimators implemented in `src/physics/`. The
> .cursorrules of the repository imposes that any change to those estimators
> triggers a synchronisation suggestion against this section, keeping
> methodology and code coherent.

### 4.1 Data

- **SPARC kinematics.** The Spitzer Photometry and Accurate Rotation Curves
  catalogue (Lelli, McGaugh & Schombert 2016) provides per-radius
  decompositions of `V_obs`, `V_gas`, `V_disk` and `V_bulge` for 175 spiral
  and dwarf galaxies, along with inclination, distance and morphology.
- **Planck polarisation.** Stokes `Q` and `U` maps of Galactic dust emission
  from the Planck Legacy Archive (Planck Collaboration 2016), accessed via
  the HEALPix pixelisation (Górski et al. 2005). The pipeline uses
  memory-mapped reads of the FITS payload to scale to the full sky.
- **NED.** The NASA/IPAC Extragalactic Database is queried for celestial
  coordinates `(RA, Dec)`, redshift and distance hints not present in
  SPARC, with on-disk caching to keep the run reproducible.

### 4.2 Estimators

The pipeline composes four Single-Responsibility estimators, orchestrated
by the `PhysicsEngine`:

- `BaryonicGapCalculator` evaluates `a_obs(r) = V_obs²/r`, `a_bar(r) = V_bar²/r`
  with `V_bar² = V_gas² + V_disk² + V_bulge²`, and the per-shell
  `acceleration_gap = a_obs − a_bar`. Inclination uncertainty is propagated
  linearly into a ±1σ band on `a_obs`.
- `GasDensityEstimator` inverts the gas-only rotation component into a
  volumetric ionised-gas density. Two equivalent paths
  `M_gas^A = V_gas²·r/G` and
  `M_gas^B = max(V_bar² − V_disk² − V_bulge², 0)·r/G` are evaluated and
  reconciled with `max(A,B)` to prevent over-inflation; the volumetric
  density follows from a cylindrical-shell decomposition with shell volume
  `ΔV_i = π(r_i² − r_{i−1}²)·h`. A configurable density floor and a
  stellar-dominance gate suppress shells where the inversion is meaningless.
- `MagneticFieldCalibrator` converts the Stokes amplitude
  `S = √(Q² + U²)` at the galaxy's sky position into a Tesla-valued
  magnetic-field strength `B = √(2 μ₀ u_ref) · S / S_ref`, where `u_ref` is
  the CMB energy density (a dimensional anchor, not a physical assertion
  about the origin of the disk field) and `S_ref` is the median of `S` over
  the loaded sky map (self-referential calibration). The choice of anchor
  is an explicit sensitivity axis.
- `MagneticPressureGradientEstimator` estimates `∇u_B` as a scalar finite
  difference between the central HEALPix pixel and the mean of its angular
  neighbours, divided by a physical step `Δx`. Two parallel `Δx` values are
  tracked: a scenario knob `scenario_gradient_scale_kpc ∈ {0.05, 0.1, 0.5}`
  and a physical
  `Δx_phys = D_galaxy · θ_neighbour` derived from the catalogued distance
  and the angular HEALPix neighbour spacing.
- **Disk-frame deprojected gradient (Phase 1 refinement).** On the
  physical-Δx path, the scalar `Δu_B/Δx` is replaced by a directional gradient
  obtained from a planar least-squares fit
  `u_{B,i} = u_{B,0} + g_x\,\Delta x_i + g_y\,\Delta y_i` over the HEALPix
  neighbour ring. When the disk position angle (PA) and inclination `i` are
  both available, the sky-plane offsets are rotated so the major axis aligns
  with `+x_{disk}` and the minor leg is rescaled by `1/\cos i`. The reported
  gradient magnitude is the disk-frame norm `\sqrt{g_x^2 + g_y^2}`. When PA
  is missing, the engine falls back to the axisymmetric magnitude (with the
  `\cos i` linear correction) and raises `warn_no_pa`.
- **Magnetic-tension toy term** (`TensionTermEstimator`,
  `src/physics/diagnostics/tension.py`). The radial component of
  `(\mathbf B\cdot\nabla)\mathbf B/\mu_0` is evaluated under a logarithmic-spiral
  field-line geometry with pitch angle `p`, for which the curvature radius
  reduces to `R_c(r) = r/\sin^2 p`. The reported acceleration is then
  `a_{\mathrm{tension},r}(r) = -\cos p\cdot B^2(r)/[\mu_0\,\rho_{\mathrm{ion}}(r)\,R_c(r)]`.
  Morphologies without a coherent disk spiral (`Im`, `BCD`) force the term to
  zero and raise `warn_pure_pressure`. The engine writes the side-by-side
  per-shell columns `accel_pressure_only_m_s2` (alias of the scenario-Δx
  pressure acceleration) and
  `accel_full_lorentz_m_s2 = accel_pressure_only + a_{\mathrm{tension},r}`,
  together with a parallel `theoretical_gap_coverage_ratio_full_lorentz` so
  the impact of the tension term can be audited *at the same Δx as the
  headline ratio*.
- **Monte-Carlo uncertainty propagation**
  (`src/physics/diagnostics/uncertainty.py`). The catalogue 1-σ uncertainties
  on distance `D`, inclination `i` and asymptotic flat velocity `V_{\rm flat}`
  are propagated onto `theoretical_gap_coverage_ratio` and
  `convergence_ratio_to_wim` through a closed-form multiplicative rescaling
  `R(D_s, i_s, V_s) = R_0\,(D_0/D_s)\,(\sin i_s/\sin i_0)^2\,(V_0/V_s)^2`,
  drawing `N=300` samples with seed `42`. The result is exposed as
  `theoretical_gap_coverage_ratio_p16`, `theoretical_gap_coverage_ratio_p84`,
  `convergence_ratio_to_wim_p16`, `convergence_ratio_to_wim_p84` in the
  headline CSV.
- **Anchor-sensitivity sweep**
  (`src/physics/diagnostics/anchor_sensitivity.py`). The MHD calibration
  anchor `u_{\rm ref}` is varied between the CMB energy density (default) and
  a turbulent-equipartition pressure floor
  `u_{\rm turb} = 0.5\,\rho_{\rm gas}\,\sigma_v^2`. The resulting two ratios
  per galaxy are persisted to a dedicated
  `results/step_validation/validation_anchor_sensitivity_sweep.csv` so the headline ratio
  can be audited against the anchor choice.

The engine also reports a per-shell `required_density_for_match`, i.e. the
ionised-plasma density `ρ_req` that would close the kinematic gap given the
observed magnetic-pressure gradient:

\[
\rho_{\mathrm{req}}(r) \;=\;
\frac{\bigl|\,\nabla u_B(r)\,\bigr|}{\bigl|\,\Delta a(r)\,\bigr|}.
\]

When `Δx_phys` is available (and finite), `∇u_B` in the numerator is taken at
the *physical* HEALPix neighbour separation, and the underperforming mask
that selects the shells contributing to the median uses
`|a_B(Δx_phys)|`. When `Δx_phys` is missing the engine falls back to the
scenario-driven gradient and to the scenario-driven `a_B`. This branching is
required so that **Phase 1 (Discovery census)** (which fixes `Δx = Δx_phys` and
`f_ion = 1`) anchors the WIM convergence check to the real Planck pixel,
while **Phase 2 (Validation — scenario grid)** keeps exploring the (Δx, f_ion) prior.

The population table `validation_scientific_evidence_report.csv` (written after the grid) merges, for each
galaxy that passes the direct-evidence ratio band, the Discovery headline fields
(`required_density_for_match`, `convergence_ratio_to_wim`, `delta_r_m`, `physical_pixel_scale_kpc`)
together with boolean flags `discovery_is_reliable` and `discovery_is_elite`, alongside the
best-scenario `theoretical_gap_coverage_ratio` and `pearson_correlation`.

### 4.3 Scenario grid

Per galaxy the engine evaluates the cartesian product of
`scenario_gradient_scale_kpc ∈ {0.05, 0.1, 0.5}` (kpc) and
`scenario_ionization_fraction ∈ {0.001, 0.01, 0.1}`, yielding nine
scenarios. For each scenario, the per-shell `magnetic_pressure_acceleration`
is compared to the per-shell `acceleration_gap` via
the outer-disk modulus ratio
`theoretical_gap_coverage_ratio = |⟨a_B⟩_edge| / |⟨gap⟩_edge|` and the
Pearson correlation of the two radial profiles over reliable shells.

### 4.4 Reliability audit

Each galaxy carries an explicit `audit_flags` list that documents the
gates it failed: `fail_lat` (foreground latitude), `fail_reliable_fraction`
(too few shells survived the outlier gate), `fail_density_clamped` (gas
density clamped to the floor on most shells), `fail_stellar` (stellar
gravity dominates the kinematics), `fail_skipped_physics` (the engine could
not produce `a_B`), `fail_inclination` (catalogued inclination below a
threshold), `fail_calibration` (Stokes reference fell back to an emergency
value). Non-blocking warnings include `warn_b_band` (inferred `B` outside
the literature plausibility band), `warn_no_pa` (no curated disk position
angle available — the engine fell back to the axisymmetric gradient),
`warn_pure_pressure` (irregular / BCD morphology — the magnetic-tension toy
collapses to zero), and `warn_btfr_outlier` (galaxy lies more than the
configured threshold in dex from the McGaugh+2012 baryonic Tully–Fisher
ridge). On each **scenario row**, the derived boolean `is_reliable` is true if and only if no
`fail_*` code is present; `warn_*` codes are informational and do not demote scenario `is_reliable`.
Galaxy-level inventory uses `validation_is_reliable` in `validation_all.csv` (and in
`data_reliability_assessment` inside each `*_summary.json`). Phase 1 Discovery census tables use
`discovery_is_reliable` (`results/step_discovery/discovery_reliable.csv` / `discovery_not_reliable.csv`). The consistency audit additionally surfaces two
population-level diagnostics over the reliable sample: the median
absolute residual against the Radial Acceleration Relation
(`rar_median_abs_residual_dex`) and the fraction of BTFR outliers
(`warn_btfr_outlier`).

### 4.5 Synthetic negative controls

The pipeline is also exercised on synthetic rotation curves generated from
Navarro–Frenk–White (NFW) dark-matter haloes (Navarro, Frenk & White 1996).
By construction these curves have no real magnetic gradient; the pipeline
must therefore *fail* to find a magnetic-pressure match for them. This single
test bounds one specific class of false positives.

---

## 5. Falsification Protocol

The model is considered **falsified at the configured scenario** if either:

1. `theoretical_gap_coverage_ratio` differs from unity by more than one order
   of magnitude in the population of reliable galaxies; or
2. the Pearson correlation between the radial profiles of `a_B` and `Δa`
   is consistent with zero over reliable radii.

The model is considered **consistent at the configured scenario** if
both quantities lie in the band described in §4 simultaneously over a
non-trivial fraction of the reliable subset. Consistency in one scenario
does not preclude refinement; in particular, items in §6 must be addressed
before a near-unity coverage ratio can be claimed *statistically*
distinguishable from chance.

The protocol applies equally to confirmation and to falsification, and the
pipeline reports both outcomes with the same audit trail.

---

## 6. Limitations and Planned Refinements

A complete enumeration is in Section 5 of [`EG_MODEL_THEORY_EN.md`](../theory/EG_MODEL_THEORY_EN.md)
and in the operator-facing guides under [`docs/results_interpretation/`](../results_interpretation/). The principal items are:

1. **Disk-frame deprojection.** The scalar gradient `Δu_B / Δx` is
   sky-projected; the dynamically relevant component is the disk-plane
   radial one. A multi-neighbour ring sampling with deprojected directional
   derivative is the next step.
2. **Bayesian uncertainty propagation.** Inclination, distance and
   `V_flat` uncertainties should be propagated through the coverage ratio
   via a Monte-Carlo before any near-unity value is interpreted as
   statistically significant.
3. **Gas-mass cross-check.** Integrated `gas_mass_density × shell_volume`
   should be cross-validated against the SPARC HI mass `MHI_1e9_solar`.
4. **Anchor sensitivity.** The CMB-energy-density anchor is dimensional;
   alternative anchors (turbulent equipartition, cosmic-ray equipartition,
   fixed Galactic disk B) should be scanned and reported.
5. **Magnetic tension.** The pressure-only ansatz drops the tension term
   `(B · ∇)B/μ₀`; a toy logarithmic-spiral field-line geometry would let
   the engine report `pressure-only` and `pressure+tension` columns side
   by side.
6. **Cosmic-ray pressure.** For star-forming dwarfs the cosmic-ray pressure
   can compete with thermal pressure and shift the plasma-β diagnostic.
7. **RAR/BTFR self-consistency.** Galaxies on which the magnetic-pressure
   model is claimed to fit must also satisfy the empirical Radial Acceleration
   Relation (McGaugh, Lelli & Schombert 2016) and the Baryonic Tully–Fisher
   Relation (McGaugh et al. 2000).
8. **Scenario knob retirement.** Once the disk-frame deprojection and the
   Bayesian propagation are in place, the legacy
   `scenario_gradient_scale_kpc` grid becomes a sensitivity sweep rather
   than the primary result.

These items are non-blocking for the falsifiability statement of §5 but are
necessary before any claim of a near-unity coverage ratio can be promoted
to a positive cosmological result.

---

## 7. Items to Revisit — Adjacent Predictions

Three lines of investigation sit at the boundary of the present analysis and
deserve a dedicated follow-up. None of them is required for the falsifiability
test of §5, but a positive outcome of that test would make each of them
considerably more interesting.

### 7.1 Plasma Lensing

If extended halos contain a substantial warm-ionised reservoir, photons
traversing them experience a frequency-dependent refractive deflection that is
absent from purely gravitational lensing. In the cold-plasma limit the index of
refraction reads

\[
n(\omega) \;=\; \sqrt{1 - \omega_p^{\,2}/\omega^{2}}, \qquad
\omega_p \;=\; \sqrt{\dfrac{n_e\, e^{2}}{\varepsilon_0\, m_e}},
\]

so that the deflection angle scales as `Δψ ∝ ∇⊥ n_e / ν²` and is therefore
**chromatic**, in contrast with the achromatic gravitational deflection of
General Relativity. For a halo with `n_e ∼ 10⁻³ cm⁻³` and a 10 kpc gradient
length, this predicts measurable frequency-dependent shifts of background radio
sources at GHz–MHz bands, accessible with VLBI multi-frequency interferometry
and pulsar timing arrays. A null result in the cleanest sight-lines (away from
the Galactic plane and from intervening clusters) would tighten the bounds on
the column-integrated electron density along the relevant lines of sight and
therefore on the WIM mass budget invoked in §3.3.

Relevant references: Pen & King (2012); Cordes et al. (2017); Vedantham (2019).

### 7.2 Primordial Magnetogenesis

The MHD ansatz of §3 takes the existence of a coherent magnetic field of order
`B ∼ 1`–`10 μG` as a *boundary condition* sampled by the Planck polarisation
maps. Its **origin** is left open by the present analysis. Three established
mechanisms are commonly invoked to seed cosmological magnetic fields:

- The **Biermann battery**, which generates a magnetic field whenever electron
  density and electron temperature gradients are non-parallel,
  \[
  \dfrac{\partial \mathbf{B}}{\partial t} \;=\; -\,\dfrac{c\, k_B}{e\, n_e}\,
  \nabla T_e \times \nabla n_e.
  \]
- **Inflationary** mechanisms that break the conformal invariance of
  electromagnetism during the inflationary epoch.
- **Phase-transition** mechanisms (electroweak or QCD) that produce a stochastic
  primordial field through bubble-wall collisions.

The seed strengths predicted by these channels span many decades; subsequent
amplification by the galactic dynamo (§7.3) brings them up to the μG band
sampled by the present pipeline. Establishing which seed channel is consistent
with the spectrum and coherence length of the field reconstructed by the E-G
calibrator would close one of the largest open conceptual gaps of the model.

Relevant references: Subramanian (2016); Durrer & Neronov (2013); Widrow (2002).

### 7.3 Self-Excited Dynamo

The use of Planck Stokes amplitudes as a **steady-state** probe of the disk
magnetic field is justified only if the field is sustained against ohmic
diffusion on timescales much shorter than the age of the disk. The relevant
equation is the magnetohydrodynamic induction equation,

\[
\dfrac{\partial \mathbf{B}}{\partial t} \;=\; \nabla \times (\mathbf{v} \times \mathbf{B})
\;+\; \eta\, \nabla^{2} \mathbf{B},
\]

whose first term advects and stretches the field with the velocity flow `v`
while the second term damps it with magnetic diffusivity `η`. In a differentially
rotating disk this admits a self-excited `αΩ` dynamo: differential rotation `Ω`
generates the toroidal component, while small-scale helical turbulence (the
`α`-effect) regenerates the poloidal component. Galactic-disk values
`τ_amp ∼ 10⁸`–`10⁹ yr ≪ τ_ohm` certify that the μG-scale field is in dynamical
equilibrium with the rotation curve. Within the E-G model this is the
self-consistency hook that lets the same `V_obs(r)` profile both *enter* the
acceleration gap (left-hand side) and *sustain* the field that closes it
(right-hand side); a quantitative test of that coupling — for instance by
correlating `pearson_correlation` with proxies of `Ω(r)` shear — is a natural
next step.

Relevant references: Beck (2015); Brandenburg & Subramanian (2005);
Ruzmaikin, Shukurov & Sokoloff (1988).

---

## 8. References

- Alfvén, H., 1970, *Cosmology in the Plasma Universe* (Nobel Lecture).
- Beck, R., 2015, *A&A Rev.*, 24, 4 (galactic magnetic fields).
- Bosma, A., 1981, *AJ*, 86, 1825.
- Brandenburg, A., Subramanian, K., 2005, *Phys. Rep.*, 417, 1 (astrophysical
  magnetic fields and the αΩ dynamo).
- Cordes, J. M., Wasserman, I., Hessels, J. W. T., et al., 2017, *ApJ*, 842, 35
  (plasma lensing).
- Durrer, R., Neronov, A., 2013, *A&A Rev.*, 21, 62 (cosmological magnetic fields).
- Górski, K. M., Hivon, E., Banday, A. J., et al., 2005, *ApJ*, 622, 759.
- Haffner, L. M., Dettmar, R.-J., Beckman, J. E., et al., 2009, *Rev. Mod.
  Phys.*, 81, 969.
- Jackson, J. D., 1999, *Classical Electrodynamics*, 3rd edn (Wiley).
- Lelli, F., McGaugh, S. S., Schombert, J. M., 2016, *AJ*, 152, 157
  (SPARC catalogue).
- McGaugh, S. S., Schombert, J. M., Bothun, G. D., de Blok, W. J. G., 2000,
  *ApJ*, 533, L99.
- McGaugh, S. S., Lelli, F., Schombert, J. M., 2016, *Phys. Rev. Lett.*,
  117, 201101 (Radial Acceleration Relation).
- Navarro, J. F., Frenk, C. S., White, S. D. M., 1996, *ApJ*, 462, 563.
- Pen, U.-L., King, L., 2012, *MNRAS*, 421, L132 (plasma gravitational lensing).
- Peratt, A. L., 1986, *IEEE Trans. Plasma Sci.*, 14, 763.
- Planck Collaboration, 2016, *A&A*, 594, A1; A&A 596, A105.
- Reynolds, R. J., 1990, *ApJ*, 349, L17.
- Rubin, V. C., Ford, W. K., 1970, *ApJ*, 159, 379.
- Ruzmaikin, A. A., Shukurov, A. M., Sokoloff, D. D., 1988, *Magnetic Fields
  of Galaxies* (Kluwer).
- Subramanian, K., 2016, *Rep. Prog. Phys.*, 79, 076901 (origin of cosmic
  magnetic fields).
- Vedantham, H. K., 2019, *MNRAS*, 488, 4757 (chromatic refraction by
  ionised media).
- Widrow, L. M., 2002, *Rev. Mod. Phys.*, 74, 775 (origin of galactic and
  extragalactic magnetic fields).
