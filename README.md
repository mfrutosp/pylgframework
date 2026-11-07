# pylgframework

A software framework for the Local Group.

The barycentre of the Local Group is routinely adopted as a point estimate. `pylgframework` replaces the
single-line centre-of-mass calculation with a Monte Carlo that propagates the observational uncertainty of the
member galaxies, and returns the barycentre together with its uncertainty region.

> Frutos Plaza, M. (2026). *Membership Systematics of the Local Group Barycentre: Uncertainty Region, Dependence
> on Member Galaxies, and Consistency with the CMB Dipole.* (manuscript; data and code at [doi:10.5281/zenodo.23440771](https://doi.org/10.5281/zenodo.23440771))

## Status

[![tests](https://github.com/mfrutosp/pylgframework/actions/workflows/tests.yml/badge.svg)](https://github.com/mfrutosp/pylgframework/actions/workflows/tests.yml)

## Installation

Requires Python ≥ 3.11

```bash
pip install pylgframework             # core: numpy, pandas, astropy
pip install "pylgframework[rebound]"  # adapter for REBOUND simulations
```

## Installation (for developers)

```bash
git clone https://github.com/mfrutosp/pylgframework
cd pylgframework
pip install -e .                 # core: numpy, pandas, astropy
pip install -e ".[rebound]"      # adapter for REBOUND simulations
pip install -e ".[test]"         # development and tests
```

## Quick start

```python
from pylgframework import barycentre, compare, fiducials

res = barycentre(fiducials.WU, n=1000, seed=2000)
print(res)                       # position, velocity, principal axes, CMB-frame speed, assumptions

res.median_position              # kpc, Galactocentric
res.axes                         # 1-sigma principal axes of the barycentre cloud (kpc)
res.vcmb("4")                    # speed w.r.t. the CMB (16th, 50th, 84th percentile, km/s)
res.cloud                        # DataFrame with every realization

# Paired audit of two literature fiducials (identical random draws)
compare(fiducials.WU, fiducials.VDM2012, n=1000, seed=3000)
```

Available fiducials: `WU` (Wu et al. 2026), `SAWALA` (Sawala et al. 2025), `VDM2019`, `VDM2012` and
`VDM2012_3B` (van der Marel et al. 2012, four-body and faithful three-body variants). Each one carries its
reference and notes (`fiducials.WU.reference`).

Memberships: `"2"` (MW+M31), `"3M"` (+M33), `"3L"` (+LMC), `"4"` (all four).

## Assumptions

* Observables (distance modulus, proper motions, radial velocity, mass) are drawn from Gaussians truncated at
  ±2σ (the effective standard deviation is about 0.88σ).
* Sky positions are held fixed; the MW is at rest at the origin of the Galactocentric frame
  (R0 = 8.122 kpc, v⊙ = (12.9, 245.6, 7.78) km/s).
* The output is the uncertainty **induced by the declared errors**, not a complete physical uncertainty region.
* The uncertainty model is fixed across fiducials unless `sigmas=` is supplied.

## Reproducibility

Realization `i` of an ensemble with seed `s` uses `numpy.random.default_rng([s, i])`. This scheme is part of the
public interface. The test suite reproduces Table 1 and Sects. 3 and 4 of the manuscript:

```bash
pytest            # fast tests
```

An example using REBOUND adapter is here [tests/test_rebound_simulation.py](https://github.com/mfrutosp/pylgframework/blob/main/tests/test_rebound_simulation.py).

## References

The bibliography of the manuscript, including the sources of every fiducial and of the CMB dipole, is in
[`docs/references.bib`](https://github.com/mfrutosp/pylgframework/blob/main/docs/references.bib).

## License

MIT.
