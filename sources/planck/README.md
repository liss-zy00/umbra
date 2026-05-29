# `sources/planck/` — Planck thermal-dust polarisation maps

**Origin.** Planck Legacy Archive, ESA. The Planck Collaboration thermal-dust
polarisation reconstructions are released as full-sky HEALPix maps in FITS
format. Mirrored at <http://pla.esac.esa.int>.

The pipeline samples the **Stokes Q/U amplitude** at each galaxy's sky
position to anchor the magnetic-field calibration. The angular separation
between the central pixel and its neighbours is projected onto the galaxy
distance to obtain the physical gradient scale `Δx_phys` used by the
magnetic-pressure-gradient estimator.

## Files

These FITS maps are **not versioned in git** (they exceed GitHub's 100 MB limit).
Download them from the [Planck Legacy Archive](http://pla.esac.esa.int) and place
them in this directory before running the pipeline.

| File | Content |
| --- | --- |
| `COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits` | GNILC component-separated thermal-dust IQU map at NSIDE=2048 (preferred when available). |
| `COM_CompMap_QU-thermaldust-commander_2048_R3.00_full.fits` | Commander pipeline thermal-dust QU map, used as a fallback. |

The provider iterates through the candidate list and keeps the first FITS that
loads successfully — see [`FITS_CANDIDATES`](../../src/pipeline/constants.py).

## Consumers

- [`src/providers/planck.py`](../../src/providers/planck.py) — `PlanckProvider` loads the FITS and samples Stokes Q/U per galaxy.

## Citation

> Planck Collaboration, *Planck 2018 results. IV. Diffuse component
> separation*, A&A 641, A4 (2020); *Planck 2018 results. XII. Galactic
> astrophysics from polarisation*, A&A 641, A12 (2020).

## Note on `healpy`

The provider raises `RuntimeError` at initialisation if a FITS file is found
but `healpy` cannot be imported in the active environment. Activate the
`dark-matter` conda env (or run via `conda run --no-capture-output -n
dark-matter python src/main.py ...`) to ensure `healpy` is available.
