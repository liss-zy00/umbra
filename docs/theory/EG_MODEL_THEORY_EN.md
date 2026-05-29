# The E-G Model — Conceptual and Mathematical Theory

> **Status: Foundational document.** This file describes the *theoretical proposition*
> behind the pipeline and the mathematical translation from physical principles to code.
> It is intentionally **timeless**: it contains no execution results, no per-galaxy
> values, and no runtime artefact paths. Concrete numbers from any specific pipeline
> run live in `results/` and the operator-facing [`README.md`](../../README.md).
> A Portuguese-language mirror of this document is kept at
> [`EG_MODEL_THEORY_PT.md`](EG_MODEL_THEORY_PT.md).

---

## 1. Introduction — The Hypothesis of Transitional States

Galactic rotation curves are the cleanest empirical anomaly in modern astrophysics.
At galactocentric radii where the visible (baryonic) mass distribution predicts a
Keplerian fall-off \(V \propto r^{-1/2}\), spectroscopic measurements consistently
show **flat** velocity profiles. The standard interpretation closes this gap with a
spherical, non-baryonic **dark-matter halo** whose gravity supplies the missing
centripetal acceleration

\[
\Delta a(r) \;=\; \underbrace{\frac{V_{\mathrm{obs}}^{2}(r)}{r}}_{a_{\mathrm{obs}}}
              \;-\; \underbrace{\frac{V_{\mathrm{bar}}^{2}(r)}{r}}_{a_{\mathrm{bar}}}
              \;>\; 0
\quad\text{(at large }r\text{)},
\qquad
V_{\mathrm{bar}} \;=\; \sqrt{V_{\mathrm{gas}}^{2}+V_{\mathrm{disk}}^{2}+V_{\mathrm{bulge}}^{2}} .
\tag{1.1}
\]

The **E-G hypothesis** does not contest Newtonian gravity in its domain of
validity. It adopts a complementary stance: galactic dynamics may be
**state-dependent on scale**, with gravity and electromagnetism behaving as
co-acting, mutually consistent fundamental forces whose relative weight changes
with the geometry of the system.

* **State A — local / compact regime.** In gravitationally bound, non-magnetised
  bodies the Newtonian potential dominates. The acceleration field of a point
  mass falls off as \(1/r^{2}\); galactic centres, stellar systems and planetary
  orbits are well described in this limit.
* **State B — extended / filamentary regime.** In magnetised plasma structures
  (current sheets, galactic-disk magnetic fields, the cosmic-web filamentary
  network observed in large-scale-structure surveys) the *spatial structure of the
  magnetic field* permits stresses that decay more slowly than \(1/r^{2}\). For
  an idealised current-carrying filament the field magnitude scales as
  \(B \propto 1/r\), and the corresponding magnetic-pressure gradient does not
  vanish over the radii where dark-matter halos are conventionally invoked.

The E-G test asks the **falsifiable** question: *can a magnetised-fluid pressure
gradient, anchored to genuine polarisation measurements, supply an acceleration
of the right order of magnitude and the right radial trend at the radii where
the gap (1.1) opens up?* A "yes" does not refute dark matter; a "no" — i.e. an
inferred magnetic field that is dimensionally implausible or whose spatial
gradient is too weak by orders of magnitude — would directly falsify this
particular MHD-pressure ansatz. The pipeline is built so that both outcomes are
reported with the same scientific rigour (Section 5).

---

## 2. The MHD Bridge — From Lorentz Force to Magnetic Pressure

The conceptual bridge of the E-G model is the passage from the **microscopic**
Lorentz force on charged particles to the **macroscopic** magnetic stress on a
fluid element. This translation is the core physics that the pipeline operator
implements; without it, the proposal would be metaphysical.

### 2.1 Microscopic origin

A single charged particle of charge \(q\) and velocity \(\mathbf{v}\) in an
electromagnetic field \((\mathbf{E},\mathbf{B})\) experiences

\[
\mathbf{F}_{\mathrm{Lorentz}} \;=\; q\,\big(\mathbf{E} + \mathbf{v}\times\mathbf{B}\big).
\tag{2.1}
\]

This is the textbook Lorentz force. For an interstellar plasma, however,
(2.1) is not the right level of description: the relevant timescales are long
compared to the particle gyration period, electron and ion species are
strongly coupled by collisions and electrostatic feedback, and the dynamical
quantity of interest is not the trajectory of one particle but the bulk
acceleration of a fluid element.

### 2.2 Continuum form: the volume Lorentz force

Summing (2.1) over the charge content of a fluid parcel and assuming
quasi-neutrality (\(\rho_{e} \to 0\), so the electric term vanishes in the
co-moving frame), the **volume Lorentz force per unit volume** reduces to

\[
\mathbf{f}_{\mathrm{Lorentz}} \;=\; \mathbf{J}\times\mathbf{B}
\;=\; \frac{1}{\mu_{0}}\big(\nabla\times\mathbf{B}\big)\times\mathbf{B},
\tag{2.2}
\]

where Ampère's law (with displacement current neglected, the **non-relativistic
MHD limit**) was used to write the current density in terms of the magnetic
field. A standard vector identity decomposes (2.2) into two physically distinct
contributions:

\[
\frac{1}{\mu_{0}}(\nabla\times\mathbf{B})\times\mathbf{B}
\;=\;
-\,\nabla\!\left(\frac{B^{2}}{2\mu_{0}}\right)
\;+\;
\frac{1}{\mu_{0}}(\mathbf{B}\cdot\nabla)\mathbf{B} .
\tag{2.3}
\]

The first term is a **pressure** — it acts as the gradient of the scalar
**magnetic energy density**

\[
u_{B}(\mathbf{x}) \;=\; \frac{B(\mathbf{x})^{2}}{2\mu_{0}} \;=\; P_{B}(\mathbf{x}),
\tag{2.4}
\]

and points from regions of higher magnetic energy toward regions of lower
magnetic energy. The second term is the **magnetic tension**, anisotropic in the
direction of \(\mathbf{B}\); it acts as a restoring stress along curved field
lines, playing the role of an "elastic" component.

### 2.3 The ideal-MHD pressure-gradient ansatz

The pipeline adopts the **ideal-MHD pressure-gradient approximation**: along the
direction in which (2.3) is sampled, the tension term is dropped. Dividing the
remaining pressure-gradient force by the inertial mass density of the magnetised
fluid yields the **magnetic-pressure acceleration**

\[
\boxed{\;
\mathbf{a}_{B}(\mathbf{x})
\;\approx\;
-\,\frac{\nabla u_{B}(\mathbf{x})}{\rho_{\mathrm{ion}}(\mathbf{x})}
\;=\;
-\,\frac{\nabla\big(B^{2}/2\mu_{0}\big)}{\rho_{\mathrm{ion}}(\mathbf{x})} .
\;}
\tag{2.5}
\]

Equation (2.5) is the **central operational equation** of the E-G test. The
denominator is the *responsive* mass density — the ionised gas component that
couples electromagnetically to the field. Its modelling is the subject of
Section 3.

The sign convention is physical: a magnetic-energy *decrease* outward
(\(\Delta u_{B} < 0\)) produces a **positive** outward radial acceleration in
the simplified one-dimensional reading of (2.5). Hence the conceptual analogy
to a "wind" of magnetic pressure pushing the magnetised plasma outward where the
field is weakening.

The justification for dropping the tension term is itself a falsifiable claim
(see Section 5.4); it is the simplest first-order MHD ansatz that lets a
data-only pipeline produce an order-of-magnitude estimate of \(\mathbf{a}_{B}\)
without requiring a global three-dimensional reconstruction of the galactic
magnetic field topology.

---

## 3. Mathematical Implementation — From Theory to Code

The pipeline is the literal numerical translation of (1.1) and (2.5) into a
sequence of estimators, each with a single physical responsibility. This section
states what each estimator computes and *why*; the source code in
`src/physics/` matches this section line-by-line.

### 3.1 The kinematic acceleration deficit

The first estimator (`BaryonicGapCalculator`) reads the SPARC rotation table
and produces SI accelerations:

\[
a_{\mathrm{obs}}(r) = \frac{V_{\mathrm{obs}}^{2}(r)}{r},\qquad
a_{\mathrm{bar}}(r) = \frac{V_{\mathrm{bar}}^{2}(r)}{r},\qquad
\mathrm{gap}(r) = a_{\mathrm{obs}}(r) - a_{\mathrm{bar}}(r).
\tag{3.1}
\]

The baryonic synthesis explicitly composes the SPARC components

\[
V_{\mathrm{bar}}(r) \;=\; \sqrt{V_{\mathrm{gas}}^{2}(r) + V_{\mathrm{disk}}^{2}(r) + V_{\mathrm{bulge}}^{2}(r)} ,
\tag{3.2}
\]

which is the rotational analogue of the baryonic velocity used in the
RAR/BTFR literature. Inclination uncertainty propagates linearly into a
symmetric ±1σ band on \(a_{\mathrm{obs}}\).

### 3.2 The gas-density estimator

The denominator of (2.5) requires a volumetric **ionised** mass density. The
estimator (`GasDensityEstimator`) uses the spherical-shell inversion of the
SPARC gas-only velocity, evaluated through **two equivalent paths** that must
agree by construction:

\[
M_{\mathrm{gas}}^{(A)}(r)
\;=\;
\frac{V_{\mathrm{gas}}^{2}(r)\,r}{G},
\qquad
M_{\mathrm{gas}}^{(B)}(r)
\;=\;
\frac{\max\!\big(V_{\mathrm{bar}}^{2}-V_{\mathrm{disk}}^{2}-V_{\mathrm{bulge}}^{2},\,0\big)\,r}{G}.
\tag{3.3}
\]

By the SPARC orthogonality (3.2) the two paths are mathematically identical;
the estimator evaluates both so that any leak between catalogue components is
auditable as a non-zero per-shell residual. Path B explicitly subtracts the
**stellar potential** from the total baryonic potential — a shell where stars
already explain the dynamics yields a vanishing residual and the conservative
\(\max(A,B)\) choice prevents the estimator from inflating the inertial reservoir.

The volumetric density follows from a cylindrical-shell decomposition. With
shell thickness \(h(R)\) configured from the catalogue extent and shell volume
\(\Delta V_{i} = \pi (r_{i}^{2}-r_{i-1}^{2})\,h\), the shell mass and density
are

\[
\Delta M_{i} \;=\; M_{\mathrm{gas}}(r_{i}) - M_{\mathrm{gas}}(r_{i-1}),\qquad
\rho_{\mathrm{gas}}(r_{i}) \;=\; \frac{\Delta M_{i}}{\Delta V_{i}} .
\tag{3.4}
\]

The effective ionised density used in (2.5) is

\[
\rho_{\mathrm{ion}}(r) \;=\; f_{\mathrm{ion}} \cdot \rho_{\mathrm{gas}}(r) ,
\tag{3.5}
\]

where \(f_{\mathrm{ion}} \in (0, 1]\) is a scenario knob exploring the unknown
ionisation fraction of the warm-ionised medium that responds to the magnetic
gradient on the relevant timescale.

A **density floor** prevents shells where the kinematic inversion yields
non-finite or non-positive densities from corrupting (2.5); a **stellar
dominance gate** removes shells where stellar gravity supplies most of the
squared-velocity budget, since the gas-only inversion is not a meaningful
density tracer there.

### 3.3 The Planck calibrator — Stokes amplitudes to Tesla

The Planck satellite measures linear polarisation through the Stokes
parameters \(Q(\hat n)\) and \(U(\hat n)\) on the sky. Their amplitude

\[
S(\hat n) \;=\; \sqrt{Q(\hat n)^{2} + U(\hat n)^{2}}
\tag{3.6}
\]

is dimensionful in map units, *not* in Tesla in the interstellar medium. The
calibrator (`MagneticFieldCalibrator`) converts (3.6) to a magnetic-field
strength via a **dimensional bridge**:

\[
B(\hat n) \;=\; \sqrt{2\,\mu_{0}\,u_{\mathrm{ref}}}\;\cdot\;\frac{S(\hat n)}{S_{\mathrm{ref}}} ,
\tag{3.7}
\]

where \(u_{\mathrm{ref}}\) is a fixed reference energy density (in SI J/m³)
and \(S_{\mathrm{ref}}\) is a calibration scale derived from the map itself.

Two important interpretive points:

1. \(u_{\mathrm{ref}}\) is the **CMB energy density**
   \(u_{\mathrm{CMB}} = a_{\mathrm{rad}} T_{\mathrm{CMB}}^{4}\) used as a
   *dimensional anchor*. This is **not** an assertion that the Galactic
   magnetic field equals the CMB field. It is a phenomenological choice:
   the CMB provides a cosmologically uniform energy-density scale with the
   correct units, the polarisation amplitude is normalised relative to itself,
   and the resulting \(B\) carries the appropriate Tesla dimensions for use
   in (2.4) and (2.5).
2. \(S_{\mathrm{ref}}\) is the **median** of (3.6) over the loaded sky map,
   which makes the calibration **self-referential** to the map content and
   independent of the absolute map scaling.

Sensitivity to the choice of \(u_{\mathrm{ref}}\) is an explicit refinement axis
(Section 5.3).

### 3.4 The magnetic-pressure gradient

The gradient \(\nabla u_{B}\) in (2.5) is estimated via a **scalar finite
difference** between the central HEALPix pixel pointing to the galaxy and the
mean of its angular neighbours:

\[
\Delta u_{B} \;=\; u_{B}^{(\mathrm{neighbours})} - u_{B}^{(\mathrm{centre})},
\qquad
\frac{\nabla u_{B}\,\cdot\,\hat e}{1} \;\approx\; \frac{\Delta u_{B}}{\Delta x} ,
\tag{3.8}
\]

where \(\Delta x\) is the physical separation between the centre and its
neighbours, and \(\hat e\) is the (unknown) unit vector linking those two
positions in the disk frame. The pipeline tracks **two parallel** \(\Delta x\)
choices: an exploratory scenario knob (50–500 pc, varied across the grid) and
a **physical** \(\Delta x_{\mathrm{phys}} = D_{\mathrm{galaxy}}\,\theta_{\mathrm{neighbour}}\)
computed from the catalogued distance and the angular HEALPix neighbour
spacing. The two values produce two parallel acceleration columns; the
honest reading of the data is the physical column, and the comparison between
them quantifies the order-of-magnitude geometric mismatch.

### 3.5 The integrated MHD acceleration

Combining (3.4)–(3.8), the pipeline evaluates (2.5) shell by shell:

\[
a_{B}(r) \;=\; -\,\frac{\Delta u_{B}}{\Delta x \cdot \rho_{\mathrm{ion}}(r)} .
\tag{3.9}
\]

Per-shell sanity diagnostics derived from (3.9) include:

* the **plasma-β** ratio of thermal to magnetic pressure,
  \(\beta = 2\mu_{0} n k_{B} T / B^{2}\), which marks the regime where the
  ideal-MHD pressure-gradient assumption is locally consistent (\(\beta < 1\))
  versus the regime where thermal pressure swamps magnetic pressure
  (\(\beta \gg 1\), where the ansatz is *out of regime*);
* the **Alfvén velocity** \(v_{A} = B / \sqrt{\mu_{0}\,\rho_{\mathrm{ion}}}\),
  the speed at which transverse magnetic perturbations propagate;
* a **plausibility band** on the inferred \(B\) so that values falling far
  outside literature-consistent disk fields surface as audit warnings.

The comparison metric of the test is the outer-disk modulus ratio

\[
\theta_{\mathrm{coverage}}
\;=\;
\frac{|\langle a_{B}\rangle_{\mathrm{edge}}|}{|\langle\mathrm{gap}\rangle_{\mathrm{edge}}|}
\tag{3.10}
\]

evaluated on the radii where dark-matter halos are conventionally invoked.
Values close to unity indicate that the magnetic-pressure proxy is *of the
right order of magnitude* to close the kinematic deficit at those radii;
values far from unity directly falsify the ansatz at the configured scenario.
A complementary Pearson correlation along reliable radii independently checks
whether the **shape** of the two profiles agrees, not just their integrated
amplitude.

---

## 4. Assumptions and Supporting Evidence

The E-G ansatz sits on five empirical pillars from established astrophysics.
None of them are private to this work; the pipeline only assembles them in a
single falsifiable test.

### 4.1 Galactic magnetic fields exist

Diffuse synchrotron emission, Faraday rotation of background polarised
sources, and direct polarisation maps from Planck and earlier CMB experiments
all confirm coherent magnetic fields permeating spiral and dwarf galaxies at
typical strengths \(\mathcal O(1\text{–}10\,\mu\mathrm{G})\). The *existence*
of a Galactic \(B(\mathbf{x})\) is therefore not an assumption but an
observational input; what (2.5) quantifies is its dynamical *contribution*.

### 4.2 Ideal magnetohydrodynamics applies on galactic scales

For typical interstellar conditions (high conductivity, low resistivity, large
length scales) the **ideal-MHD** approximation
\(\mathbf{E} + \mathbf{v}\times\mathbf{B}=0\) holds to high accuracy. Hannes
Alfvén's discovery that magnetic field lines are *frozen into* a perfectly
conducting fluid is the licence to write (2.2) as a divergence of magnetic
stresses and to retain only the pressure component on the largest scales.
Resistive corrections, anisotropic pressure tensors, and cosmic-ray pressure
are real but second-order; treating them is the subject of the refinement
program in Section 5.

### 4.3 The Universe is filamentary on the relevant scales

Large-scale-structure surveys (SDSS, eBOSS, DESI) and N-body simulations
consistently reveal the **cosmic web**: a network of filaments, sheets and
voids in which dark and baryonic matter are arranged. In an idealised
straight current-carrying filament the magnetic field magnitude scales as
\(B \propto 1/r\), so the corresponding magnetic-pressure gradient is
non-trivial over distances comparable to the inter-filament separation. The
E-G hypothesis proposes that the *outer* parts of galactic disks live in this
geometry rather than the spherical-halo geometry assumed by ΛCDM — a
geometric, not an exotic-particle, alternative.

### 4.4 The Stokes parameters carry magnetic-field information

The Planck satellite's polarisation maps \((Q,U)\) measure the linear
polarisation of dust thermal emission. For optically thin emission with random
intrinsic polarisation aligned by an ordered field, the polarisation angle
traces the projected magnetic-field orientation and the polarisation
*fraction* depends monotonically on the field strength. Equation (3.7) is the
minimal **dimensional** reading of that information: it does not assume any
specific dust grain model, only that polarisation amplitude is monotonic in
\(B\) and uniform over the calibration support.

### 4.5 Newtonian gravity is the local limit

The kinematic gap (1.1) is constructed assuming Newtonian gravity from the
SPARC mass decomposition. The E-G hypothesis does **not** modify this. In the
local regime — stellar interiors, planetary systems, the inner few kpc of a
galaxy where the baryonic component dominates — Newton is exact within
observational accuracy and the gap is small or zero. The MHD-pressure term
becomes relevant only at the radii where the deficit has opened up.

---

## 5. Refinement Parameters (Inherent Limitations)

The list below catalogues the **known weaknesses** of the present formulation.
Each item is an axis along which a refined version of the model could be
falsified more sharply; none of them invalidate the underlying conceptual
structure of Sections 1–4. Addressing them is a matter of extending the
pipeline (providers, diagnostics, and statistical reporting) while keeping
this theoretical core unchanged.

### 5.1 Sky-plane vs disk-plane geometric mismatch

The gradient in (3.8) is sampled along the **sky-projected** axis defined by
HEALPix neighbours, whereas the *physical* gradient that contributes to the
radial acceleration in (2.5) is the **disk-plane** radial component. The
projected scalar \(\Delta u_{B}/\Delta x\) is therefore an upper bound on the
radial component of the true gradient; resolving the directional difference
would require multi-neighbour reconstruction or a disk-frame deprojection of
the polarisation map.

### 5.2 Indirect gas density from the spherical-shell inversion

Equation (3.4) inverts a rotation velocity into a volumetric density under
the implicit assumption that the enclosed mass is spherically symmetric. An
HI disk explicitly violates this assumption: gas is concentrated in a thin
disk, not in a sphere. Path-A vs Path-B cross-checks (3.3) and the stellar
dominance gate prevent the estimator from inflating mass under stellar
contamination, but they do not turn the residual gas density into a true ISM
density. A cross-validation against catalogued HI mass profiles is the
natural next step.

### 5.3 Phenomenological status of the Stokes anchor

Equation (3.7) provides the correct **dimensions** for \(B\); it does not
provide a first-principles magnitude. A different choice for
\(u_{\mathrm{ref}}\) — say, equipartition with a fixed turbulent pressure —
would rescale \(B\) and therefore \(a_{B}\) by a constant factor. A
sensitivity scan over physically motivated anchor choices is required before
the *absolute* coverage ratio (3.10) can be interpreted with confidence.

### 5.4 Ideal-MHD ansatz: missing tension and non-thermal pressure

The pipeline retains only the magnetic-pressure-gradient component of (2.3),
dropping the magnetic-tension term \(({\bf B}\cdot\nabla){\bf B}/\mu_{0}\)
and treating the plasma as ideal (no resistivity, no anisotropic pressure
tensor, no cosmic-ray contribution, no turbulent pressure). Where the
magnetic-field topology is highly ordered the tension can compete with the
pressure; cosmic rays in actively star-forming dwarfs can carry a non-thermal
pressure comparable to the thermal one. The plasma-β diagnostic surfaces the
regime where the ansatz is locally out of regime, but does not correct it.

### 5.5 Inclination and distance uncertainties as symmetric bands

Inclination uncertainty is propagated into a symmetric ±1σ band on
\(a_{\mathrm{obs}}\) — but the coverage ratio (3.10) itself is computed from
the central value, not from the band edges. SPARC distances carry
non-trivial errors that linearly affect the physical \(\Delta x_{\mathrm{phys}}\)
in (3.8). Both should be propagated through (3.10) using a Bayesian
Monte-Carlo before any near-unity coverage value can be claimed
*statistically* distinguishable from the no-effect null.

### 5.6 Falsifiability scope of synthetic controls

The model is exercised against synthetic-NFW negative controls — dark-matter
rotation curves with no real magnetic gradient — to ensure the pipeline
**cannot** match those by construction. This single test bounds one specific
class of false positives. It does not exclude alternative non-MHD mechanisms
(modified inertia, non-spherical baryonic structure, beam-smearing artefacts
in real \(V_{\mathrm{obs}}\)) that could conceivably reproduce the observed
gap shape without invoking magnetic pressure. Stronger falsifiability
campaigns — anchor-sensitivity sweeps, cross-checks against the empirical
RAR/BTFR, the magnetic-tension correction — are implemented or tracked in the
codebase and validation artefacts under `results/` (see repository README and
`docs/results_interpretation/` for how to read them).

---

## 6. Reading Map

* **Programmer entering the codebase.** Section 3 explains *what
  `calculate_magnetic_pressure_acceleration` is computing* — it evaluates
  (3.9) using `MagneticFieldCalibrator` (Section 3.3),
  `GasDensityEstimator` (Section 3.2), and `MagneticPressureGradientEstimator`
  (Section 3.4), composed by the `PhysicsEngine` orchestrator. The kinematic
  gap is in Section 3.1 (`BaryonicGapCalculator`).
* **Physicist evaluating the proposition.** Sections 1–2 set up the
  hypothesis, the MHD bridge from (2.1)–(2.2) to the pressure-gradient
  ansatz (2.5), and the empirical pillars (Section 4). Section 5 lists the
  honest weaknesses; the software stack and result guides describe what is
  already exercised in code versus what remains open.
* **Reviewer asking what the test concludes for a specific dataset.** That
  question is answered by the runtime artefacts under `results/`
  and the operator-facing summary in [`README.md`](../../README.md). Those
  files are *not* part of this document by design — the theory must remain
  reproducible across runs and decoupled from any single execution.
