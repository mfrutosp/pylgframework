"""Optional adapter for REBOUND simulations (requires ``pip install pylgframework[rebound]``)."""
import numpy as np

from .constants import KMS, NAMES
from .frames import galactocentric_to_observables


def observables_from_sim(sim, template, names=NAMES, rotation=None):
    """Heliocentric observables of the particles of a REBOUND simulation at its CURRENT state.

    Call it on the initial state (t = 0): the observables are the present-day values.

    The simulation must use kpc, Msun and Gyr, with the particles ordered as ``names``
    (particle 0 is the Milky Way). Positions and velocities are taken relative to particle 0,
    so a simulation moved to its centre of mass (``sim.move_to_com()``) is handled directly.
    If the simulation axes are rotated with respect to the Galactocentric frame
    (``x_sim = rotation @ x_galactocentric``), pass that matrix as ``rotation``.

    Masses come from the particles and concentrations from ``template`` (a dict like
    ``Fiducial.obs``). The result can be passed to ``barycentre(..., sigmas=...)``.
    Uncertainties are defined on the observables (distance modulus, proper motions, radial
    velocity), not on Cartesian coordinates, so the correlations of the frame transformation
    are preserved.
    """
    ps = sim.particles
    pos = np.array([[p.x, p.y, p.z] for p in ps])
    vel = np.array([[p.vx, p.vy, p.vz] for p in ps]) / KMS          # kpc/Gyr -> km/s
    pos, vel = pos - pos[0], vel - vel[0]                           # MW at rest at the origin
    if rotation is not None:
        R = np.asarray(rotation, float)
        pos, vel = pos @ R, vel @ R                                 # rows: (R.T @ x) for each x
    obs = {}
    ob = galactocentric_to_observables(pos[1:], vel[1:])
    for j, nm in enumerate(names):
        d = dict(m=float(ps[j].m), c=float(template[nm]["c"]))
        if j > 0:
            d.update({k: float(ob[k][j - 1]) for k in ob})
        obs[nm] = d
    return obs
