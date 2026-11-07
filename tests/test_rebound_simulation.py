"""REBOUND simulation of the Local Group used as a worked example of ``pylgframework``.

The simulation is built as in the dynamical code of the manuscript (IAS15, G in
kpc^3 Msun^-1 Gyr^-2, four bodies MW, M31, M33, LMC, 5 kpc softening, moved to the centre of
mass), without dynamical friction or mass transfer, which are not part of this package.

Usage pattern::

    sim = build_local_group_sim(fiducials.WU)          # your own REBOUND simulation
    obs = observables_from_sim(sim, fiducials.WU.obs)  # recover the heliocentric observables
    res = barycentre(obs, n=1000, seed=2000)           # barycentre with its uncertainty
    print(res)                                         # replaces sim.com()
"""
import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from pylgframework import barycentre, fiducials
from pylgframework.barycentre import phase_space
from pylgframework.constants import KMS, NAMES
from pylgframework.rebound_adapter import observables_from_sim

G = 4.30091727e-6 * KMS**2          # kpc^3 / (Msun Gyr^2)
SOFTENING = 5.0                      # kpc


def build_local_group_sim(fiducial, rotation=None):
    """Four-body Local Group simulation from a fiducial, in the barycentric frame."""
    m, r, v = phase_space([fiducial.obs], fiducial.zero_mass)
    sim = rebound.Simulation()
    sim.integrator = "ias15"
    sim.G = G
    sim.softening = SOFTENING
    sim.collision = "none"
    for j, name in enumerate(NAMES):
        rj, vj = r[0, j], v[0, j] * KMS                          # km/s -> kpc/Gyr
        if rotation is not None:
            rj, vj = rotation @ rj, rotation @ vj
        sim.add(m=m[0, j], x=rj[0], y=rj[1], z=rj[2], vx=vj[0], vy=vj[1], vz=vj[2], name=name)
    sim.move_to_com()
    return sim


@pytest.fixture
def sim():
    return build_local_group_sim(fiducials.WU)


def test_simulation_is_in_the_barycentric_frame(sim):
    com = sim.com()
    assert np.linalg.norm([com.x, com.y, com.z]) < 1e-9
    assert np.linalg.norm([com.vx, com.vy, com.vz]) < 1e-9


def test_adapter_recovers_the_nominal_barycentre(sim):
    """sim.com() is a point; the adapter + barycentre() gives the same point, with uncertainty."""
    obs = observables_from_sim(sim, fiducials.WU.obs)
    res = barycentre(obs, n=1, seed=0, sigmas={})
    ref = barycentre(fiducials.WU, n=1, seed=0, sigmas={})
    np.testing.assert_allclose(res.nominal_position, ref.nominal_position, atol=1e-3)
    np.testing.assert_allclose(res.nominal_velocity, ref.nominal_velocity, atol=1e-3)


def test_uncertainty_from_the_simulation_matches_the_fiducial(sim):
    obs = observables_from_sim(sim, fiducials.WU.obs)
    a = barycentre(obs, n=300, seed=2000)
    b = barycentre(fiducials.WU, n=300, seed=2000)
    np.testing.assert_allclose(a.median_position, b.median_position, atol=0.05)
    np.testing.assert_allclose(a.axes, b.axes, atol=0.05)


def test_rotated_simulation_axes_are_handled():
    th = np.deg2rad(37.0)
    R = np.array([[np.cos(th), -np.sin(th), 0.0], [np.sin(th), np.cos(th), 0.0], [0.0, 0.0, 1.0]])
    sim_r = build_local_group_sim(fiducials.WU, rotation=R)
    obs = observables_from_sim(sim_r, fiducials.WU.obs, rotation=R)
    res = barycentre(obs, n=1, seed=0, sigmas={})
    ref = barycentre(fiducials.WU, n=1, seed=0, sigmas={})
    np.testing.assert_allclose(res.nominal_position, ref.nominal_position, atol=1e-3)


def test_centre_of_mass_is_ballistic_during_integration(sim):
    """Internal forces conserve momentum: the barycentre stays at the origin (bookkeeping check)."""
    sim.integrate(2.0)
    com = sim.com()
    assert np.linalg.norm([com.x, com.y, com.z]) < 1e-6
    assert np.linalg.norm([com.vx, com.vy, com.vz]) < 1e-6
