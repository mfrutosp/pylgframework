"""Heliocentric observables <-> Galactocentric phase space (vectorized)."""
import astropy.units as u
import numpy as np
from astropy.coordinates import CartesianDifferential, Galactocentric, ICRS, SkyCoord

from .constants import SUN_R0, SUN_V


def galactocentric_frame():
    return Galactocentric(galcen_distance=SUN_R0 * u.kpc,
                          galcen_v_sun=CartesianDifferential(list(SUN_V) * (u.km / u.s)))


def mu_to_kpc(mu):
    """Distance modulus -> distance in kpc."""
    return 10.0 ** ((np.asarray(mu, float) + 5.0) / 5.0) / 1.0e3


def helio_to_galactocentric(ra, dec, mu, pmra, pmdec, vlos):
    """Arrays (or scalars) of heliocentric observables -> position (K,3) kpc, velocity (K,3) km/s.

    ra, dec in deg; mu in mag; pmra (=mu_alpha*) and pmdec in uas/yr; vlos in km/s.
    """
    a = lambda x: np.atleast_1d(np.asarray(x, float))
    c = SkyCoord(ra=a(ra) * u.deg, dec=a(dec) * u.deg, distance=mu_to_kpc(a(mu)) * u.kpc,
                 pm_ra_cosdec=a(pmra) * u.uas / u.yr, pm_dec=a(pmdec) * u.uas / u.yr,
                 radial_velocity=a(vlos) * u.km / u.s, frame="icrs")
    g = c.transform_to(galactocentric_frame())
    return g.cartesian.xyz.to_value(u.kpc).T, g.velocity.d_xyz.to_value(u.km / u.s).T


def galactocentric_to_observables(pos, vel):
    """Inverse transform: position (kpc) and velocity (km/s) -> dict of heliocentric observables."""
    pos, vel = np.atleast_2d(pos), np.atleast_2d(vel)
    c = SkyCoord(x=pos[:, 0] * u.kpc, y=pos[:, 1] * u.kpc, z=pos[:, 2] * u.kpc,
                 v_x=vel[:, 0] * u.km / u.s, v_y=vel[:, 1] * u.km / u.s, v_z=vel[:, 2] * u.km / u.s,
                 frame=galactocentric_frame())
    i = c.transform_to(ICRS())
    d_pc = i.distance.to_value(u.pc)
    return dict(ra=i.ra.deg, dec=i.dec.deg, mu=5.0 * np.log10(d_pc) - 5.0,
                pmra=i.pm_ra_cosdec.to_value(u.uas / u.yr),
                pmdec=i.pm_dec.to_value(u.uas / u.yr),
                vlos=i.radial_velocity.to_value(u.km / u.s))
