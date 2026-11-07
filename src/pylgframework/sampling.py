"""Truncated-Gaussian sampling of the observables.

The draw scheme is a public contract: realization ``i`` of an ensemble with
seed ``s`` uses ``numpy.random.default_rng([s, i])`` and consumes the random
numbers in the order given by the ``sigmas`` dictionary. Changing it breaks
reproducibility of every published ensemble (and the paired design of
``compare``).
"""
import numpy as np

from .constants import TRUNC

# 1-sigma uncertainties of the 20 observables (Wu et al. 2026, Table 1).
# Masses in Msun; mu in mag; proper motions in uas/yr; vlos in km/s.
DEFAULT_SIGMAS = {
    "MilkyWay": dict(m=2.0e11, c=2.0),
    "M31": dict(m=4.0e11, c=2.5, mu=0.032, pmra=8.1, pmdec=6.6, vlos=1.0),
    "M33": dict(m=1.0e11, c=2.0, mu=0.070, pmra=5.7, pmdec=5.7, vlos=1.8),
    "LMC": dict(m=0.5e11, c=2.0, mu=0.026, pmra=20.0, pmdec=47.0, vlos=3.4),
}


def trunc_normal(rng, mean, sigma, k=TRUNC):
    """Gaussian truncated at +-k sigma, by rejection."""
    while True:
        z = rng.standard_normal()
        if abs(z) <= k:
            return float(mean + sigma * z)


def sample_obs(rng, base, sigmas=None):
    """One realization with the same structure as ``base`` (sky positions fixed)."""
    sigmas = DEFAULT_SIGMAS if sigmas is None else sigmas
    obs = {}
    for gal, sig in sigmas.items():
        obs[gal] = dict(base[gal])
        for key, s in sig.items():
            obs[gal][key] = trunc_normal(rng, base[gal][key], s)
    for gal in base:                      # galaxies without sigmas stay fixed
        obs.setdefault(gal, dict(base[gal]))
    return obs


def sample_ensemble(base, n, seed, sigmas=None):
    """List of ``n`` realizations; realization i uses ``default_rng([seed, i])``."""
    return [sample_obs(np.random.default_rng([seed, i]), base, sigmas) for i in range(n)]
