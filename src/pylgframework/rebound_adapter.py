"""Optional adapter for REBOUND simulations (requires ``pip install pylgframework[rebound]``)."""
import numpy as np

from .constants import KMS, NAMES
from .frames import galactocentric_to_observables


def observables_from_sim(sim, template, names=NAMES):
    """Heliocentric observables of the particles of a REBOUND simulation.

    The simulation must be in kpc, Msun and Gyr, in the Galactocentric frame with the MW
    (particle 0) at rest at the origin, and particles ordered as ``names``. Masses and
    concentrations come from the particles and from ``template`` (a dict like ``Fiducial.obs``).

    The result is a dict that can be passed to ``barycentre(..., sigmas=...)``. Uncertainties
    are defined on the observables (distance modulus, proper motions, radial velocity), not on
    Cartesian coordinates, so the correlations of the frame transformation are preserved.
    """
    ps = sim.particles
    pos = np.array([[p.x, p.y, p.z] for p in ps])
    vel = np.array([[p.vx, p.vy, p.vz] for p in ps]) / KMS          # kpc/Gyr -> km/s
    obs = {}
    ob = galactocentric_to_observables(pos[1:], vel[1:])
    for j, nm in enumerate(names):
        d = dict(m=float(ps[j].m), c=float(template[nm]["c"]))
        if j > 0:
            d.update({k: float(ob[k][j - 1]) for k in ob})
        obs[nm] = d
    return obs
