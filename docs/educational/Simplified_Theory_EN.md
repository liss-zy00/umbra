# The E-G Model — A Technical, Plain-Language Exposition

> **Audience.** Readers with general scientific literacy but without a
> background in plasma physics or galactic dynamics. The text is technical
> in vocabulary but linear in exposition: each section introduces only what
> the next one needs. Metaphors are deliberately avoided; analogies, when
> used, are physical.
> **Portuguese version:** [`Simplified_Theory_PT.md`](Simplified_Theory_PT.md).
> **Formal version (with full mathematics):** [`../theory/EG_MODEL_THEORY_EN.md`](../theory/EG_MODEL_THEORY_EN.md).

---

## 1. The observational problem

Spiral and dwarf galaxies rotate in a way that does not match the gravity of
their visible matter. If only the stars, the dust and the gas that telescopes
can see were present, Newtonian gravity would predict that stars in the outer
parts of a galaxy slow down with distance, just as outer planets in the Solar
System orbit more slowly than inner ones. What is actually measured is the
opposite: the rotation speed stays roughly **flat** out to the largest radii
that can be reliably observed.

The numerical statement of the discrepancy uses two accelerations at the same
radius `r`:

- the **observed** centripetal acceleration, computed from the rotation speed:
  `a_obs = V_obs² / r`;
- the **baryonic** acceleration, computed from the visible mass profile only:
  `a_bar = V_bar² / r`, with `V_bar² = V_gas² + V_disk² + V_bulge²`.

In the outer disk, `a_obs > a_bar`. The difference

```text
gap(r) = a_obs(r) − a_bar(r)
```

is the **acceleration gap**. The standard interpretation in cosmology
attributes the gap to a halo of non-baryonic *dark matter*: a substance that
has gravity but does not interact with light, has never been detected
directly, and is invoked specifically to fill this kind of missing
acceleration.

## 2. What space is actually made of

The interpretation of the gap depends on what assumptions are made about the
medium that fills a galaxy. Two facts from established astrophysics matter
here:

1. **Most of the visible Universe is plasma.** Plasma is a gas in which a
   substantial fraction of the atoms have lost one or more electrons and
   therefore carry electric charge. The interstellar and intergalactic media
   are predominantly ionised hydrogen with traces of heavier elements; they
   conduct electricity and respond to magnetic fields.
2. **Galaxies have measurable magnetic fields.** Several independent
   observational channels — Faraday rotation of background polarised sources,
   diffuse synchrotron emission, and the linear polarisation maps from the
   Planck satellite — converge on the same conclusion: spiral and dwarf
   galaxies are permeated by coherent magnetic fields with strengths of order
   one to ten microgauss.

The plasma component most relevant to this discussion is the **Warm Ionised
Medium (WIM)**: hydrogen ionised by ultraviolet radiation from massive stars,
at a temperature of roughly ten thousand Kelvin and at very low density
(of order `10⁻²⁷ kg/m³` in galactic haloes). It coexists with the neutral
gas (`HI`, of order `10⁻²² kg/m³`) that radio telescopes detect through its
21-cm line, but it is the WIM — not the neutral gas — that responds directly
to magnetic forces.

## 3. Magnetic energy as a pressure

A magnetic field of strength `B` stores energy in space. The amount of
energy per unit volume is

```text
u_B(x) = B(x)² / (2 μ₀)
```

with `μ₀` the vacuum permeability. This quantity behaves, in the equations
of magnetohydrodynamics (MHD), exactly like a **pressure** — it has the same
units as a pressure (joules per cubic metre) and, when its spatial
distribution is non-uniform, it pushes magnetised fluid from regions of
higher `u_B` toward regions of lower `u_B`. The driving quantity is
therefore not `u_B` itself but its **gradient** `∇u_B`.

The acceleration that a parcel of ionised fluid experiences from this
gradient is

```text
a_B(x) ≈ − ∇u_B(x) / ρ_ion(x)
```

where `ρ_ion` is the mass density of the responsive (ionised) component.
The minus sign means the acceleration points from high-`u_B` toward low-`u_B`.
This is **not** the textbook microscopic Lorentz force on a single charged
particle; it is the macroscopic force density obtained when the same
underlying physics is averaged over a fluid parcel. The full derivation is
in `../theory/EG_MODEL_THEORY_EN.md`, §2.

## 4. Why scale matters

Gravity from a centrally concentrated mass decays as `1/r²`. Magnetic
pressure does not decay the same way: in a structure organised in filaments
or in a disk geometry, the relevant magnetic energy gradient decays much
more slowly than `1/r²`. The exact decay depends on the geometry — an
idealised straight current-carrying filament has `B ∝ 1/r`, so its
energy-density gradient decays as `1/r³` along the filament axis but much
more slowly across the filamentary sheet.

The practical consequence is that at large galactocentric radii, where
gravity from visible matter has become weak, a magnetic-pressure gradient
that started smaller in the inner disk may end up comparable in size or
larger. The transition radius depends on the actual geometry of the field
and on the actual density of the ionised component.

## 5. The E-G hypothesis

The E-G model proposes that galactic dynamics are described by two
co-existing regimes, not by two competing theories of gravity:

- **Regime A** dominates wherever the baryonic gravitational potential is
  strong: stellar interiors, planetary systems, the innermost few kpc of a
  galactic disk. Newtonian gravity is exact within observational error and
  the gap is small or zero.
- **Regime B** is relevant where the gap has opened up: the outer parts of
  galactic disks and the diffuse haloes. There, the magnetic-pressure
  gradient defined above acts on the ionised gas component (the WIM) and
  can contribute, in principle, an acceleration of the right magnitude to
  reduce — or fully replace — the contribution attributed to dark matter.

The hypothesis is **falsifiable**. If the pipeline that computes both `gap`
and `a_B` from real data finds that `|a_B|` is many orders of magnitude
smaller than `|gap|`, the model fails. If `|a_B|` lies in the same range as
`|gap|` and the two profiles share the same radial shape over reliable
radii, the model has passed one independent observational test — though
not enough to displace the dark-matter interpretation by itself.

## 6. What the pipeline actually does

The repository implements this test on the public SPARC catalogue (175
spiral and dwarf galaxies with well-measured rotation curves and
mass decompositions) and on Planck polarisation maps (the linear-polarisation
Stokes parameters `Q` and `U` of galactic dust emission).

The steps, in order, are:

1. **Compute the gap.** For each galaxy in the catalogue, evaluate
   `a_obs(r)`, `a_bar(r)` and `gap(r)` from the SPARC rotation table.
2. **Estimate the responsive density.** Invert the gas-only velocity curve
   into a volumetric gas-density profile `ρ_gas(r)`. Multiply by a free
   ionisation fraction `f_ion ∈ (0, 1]` to obtain `ρ_ion(r) = f_ion · ρ_gas(r)`.
3. **Calibrate the magnetic field.** Convert the Planck Stokes amplitude
   `S(n̂) = √(Q² + U²)` at the galaxy's sky position into a field strength
   `B(n̂)` using a dimensional anchor (the CMB energy density, normalised
   by the map median; see `../theory/EG_MODEL_THEORY_EN.md`, §3.3).
4. **Estimate the gradient.** Compare `u_B` at the central pixel pointing
   at the galaxy with the average `u_B` of neighbouring HEALPix pixels;
   divide by the physical separation to obtain `∇u_B`.
5. **Compute `a_B`.** Combine the gradient with `ρ_ion` to obtain `a_B(r)`
   per radial shell.
6. **Compare.** Report the modulus ratio `|⟨a_B⟩_edge| / |⟨gap⟩_edge|` and
   the Pearson correlation between the radial profiles of `a_B` and `gap`
   over the radii that pass the reliability gates.

The pipeline does this for a grid of nine scenarios — three gradient scales
(50 pc, 100 pc, 500 pc) crossed with three ionisation fractions
(0.1 %, 1 %, 10 %) — and reports the outcome of every scenario, not just
the most favourable one.

## 7. Reading the result honestly

A coverage ratio close to one means that the magnetic-pressure proxy is of
the right order of magnitude to close the kinematic deficit at the radii
where dark matter is conventionally invoked. A ratio that is orders of
magnitude smaller or larger falsifies the ansatz at the configured scenario.

Two cautions are built into the reporting:

- **Audit flags.** Each galaxy is annotated with a list of conditions under
  which the result should not be trusted: low galactic latitude (foreground
  contamination), too few reliable radial shells, density floor dominating
  the inversion, stellar gravity already explaining the dynamics,
  inclination too low, calibration falling back to an emergency value.
- **Synthetic negative controls.** The pipeline is exercised against
  synthetic rotation curves generated from a Navarro–Frenk–White dark-matter
  halo (a curve that, by construction, has no real magnetic gradient). The
  pipeline must *not* find a match for these controls.

The full list of known limitations — sky-plane vs disk-plane mismatch,
indirect gas-density inversion, phenomenological calibration, dropped
magnetic-tension term, distance and inclination uncertainties — lives in
[`EG_MODEL_THEORY_EN.md`](../theory/EG_MODEL_THEORY_EN.md), Section 5, and in the
[`docs/results_interpretation/`](../results_interpretation/) guides that track CSV columns and audit flags.

## 8. What the project is not claiming

The E-G model does not claim:

- that Newtonian gravity is wrong;
- that the standard cosmological model (ΛCDM) is wrong as a whole;
- that gravitational lensing observations (such as the Bullet Cluster) are
  explained by magnetic effects;
- that the Cosmic Microwave Background fluctuations are explained without a
  cold dark matter component.

The model is restricted, by design, to the **galactic rotation-curve
problem** — the specific empirical anomaly summarised in Section 1. Its
falsification or confirmation in that single domain is the only question
the pipeline addresses.
