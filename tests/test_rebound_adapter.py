import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from pylgframework import barycentre, fiducials
from pylgframework.constants import KMS, NAMES
from pylgframework.barycentre import phase_space
from pylgframework.rebound_adapter import observables_from_sim


def test_sim_roundtrip_reproduces_the_nominal_barycentre():
    base = fiducials.WU.obs
    m, r, v = phase_space([base])
    sim = rebound.Simulation()
    sim.G = 4.30091727e-6 * KMS**2          # kpc^3 / (Msun Gyr^2); the adapter only reads state
    for j in range(4):
        sim.add(m=m[0, j], x=r[0, j, 0], y=r[0, j, 1], z=r[0, j, 2],
                vx=v[0, j, 0] * KMS, vy=v[0, j, 1] * KMS, vz=v[0, j, 2] * KMS)
    obs = observables_from_sim(sim, base)
    a = barycentre(base, n=1, seed=0, sigmas={})
    b = barycentre(obs, n=1, seed=0, sigmas={})
    np.testing.assert_allclose(b.nominal_position, a.nominal_position, atol=1e-3)
    np.testing.assert_allclose(b.nominal_velocity, a.nominal_velocity, atol=1e-3)
