"""Barycentre of the Local Group with propagated observational uncertainty."""
from dataclasses import dataclass
from typing import Any, Optional

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import ICRS, SkyCoord

from .constants import KMS, MEMBERSHIPS, NAMES, V_MW_CMB
from .fiducials import Fiducial, get as _get_fiducial
from .frames import galactocentric_frame, helio_to_galactocentric
from .sampling import DEFAULT_SIGMAS, sample_ensemble


def _resolve(source):
    """-> (name, obs dict, zero_mass tuple)."""
    if isinstance(source, str):
        source = _get_fiducial(source)
    if isinstance(source, Fiducial):
        return source.name, source.obs, source.zero_mass
    if isinstance(source, dict):
        return "custom", source, ()
    raise TypeError("source must be a Fiducial, a fiducial name or an observables dict")


def phase_space(obs_list, zero_mass=()):
    """Masses (K,4), positions (K,4,3) kpc, velocities (K,4,3) km/s (MW at rest at the origin)."""
    K = len(obs_list)
    m = np.array([[o[nm]["m"] for nm in NAMES] for o in obs_list], float)
    r, v = np.zeros((K, 4, 3)), np.zeros((K, 4, 3))
    for j, nm in enumerate(NAMES):
        if nm == "MilkyWay":
            continue
        g = lambda k: [o[nm][k] for o in obs_list]
        r[:, j], v[:, j] = helio_to_galactocentric(g("ra"), g("dec"), g("mu"),
                                                   g("pmra"), g("pmdec"), g("vlos"))
    for nm in zero_mass:
        m[:, NAMES.index(nm)] = 0.0
    return m, r, v


def _bc(m, r, v, idx):
    mm = m[:, idx]
    M = mm.sum(1)[:, None]
    return (mm[..., None] * r[:, idx]).sum(1) / M, (mm[..., None] * v[:, idx]).sum(1) / M


def build_cloud(m, r, v, membership="4"):
    """DataFrame with the barycentre of every realization and the membership diagnostics."""
    idx = MEMBERSHIPS[membership]
    rb, vb = _bc(m, r, v, idx)
    c = SkyCoord(x=rb[:, 0] * u.kpc, y=rb[:, 1] * u.kpc, z=rb[:, 2] * u.kpc,
                 frame=galactocentric_frame())
    icrs, gal = c.transform_to(ICRS()), c.transform_to("galactic")
    xyz = icrs.cartesian.xyz.to_value(u.kpc)
    df = pd.DataFrame({
        "bc0_x": rb[:, 0], "bc0_y": rb[:, 1], "bc0_z": rb[:, 2],
        "bc0_vx": vb[:, 0], "bc0_vy": vb[:, 1], "bc0_vz": vb[:, 2],
        "bc0_ra": icrs.ra.deg, "bc0_dec": icrs.dec.deg,
        "bc0_dist": icrs.distance.to_value(u.kpc),
        "bc0_l": gal.l.deg, "bc0_b": gal.b.deg,
        "bc0_icrs_x": xyz[0], "bc0_icrs_y": xyz[1], "bc0_icrs_z": xyz[2]})
    bcs = {t: _bc(m, r, v, i) for t, i in MEMBERSHIPS.items()}
    (r2, v2), (r3, v3), (r4, v4), (r5, v5) = bcs["2"], bcs["3M"], bcs["4"], bcs["3L"]
    n = lambda x: np.linalg.norm(x, axis=1)
    df["shift_M33_kpc"], df["dv_M33_kms"] = n(r3 - r2), n(v3 - v2)
    df["shift_LMC_kpc"], df["dv_LMC_kms"] = n(r5 - r2), n(v5 - v2)
    df["shift_M33_LMC_kpc"], df["dv_M33_LMC_kms"] = n(r4 - r3), n(v4 - v3)
    df["shift_2to4_kpc"], df["dv_2to4_kms"] = n(r4 - r2), n(v4 - v2)
    for T in (5.0, 10.0):
        df[f"drift_2to4_{T:g}Gyr_kpc"] = n((r4 - r2) + (v4 - v2) * KMS * T)
    for tag, (_, vv) in bcs.items():
        df[f"vcmb_{tag}_kms"] = n(vv + V_MW_CMB[None, :])
    return df


@dataclass
class Result:
    """Barycentre cloud and its summary statistics."""
    cloud: pd.DataFrame
    nominal_position: np.ndarray          # kpc, Galactocentric, unsampled central values
    nominal_velocity: np.ndarray          # km/s
    fiducial: str
    n: int
    seed: int
    membership: str
    sigmas: Any
    masses: Optional[np.ndarray] = None       # (K, 4) Msun, after ``zero_mass``
    pos: Optional[np.ndarray] = None          # (K, 4, 3) kpc, Galactocentric
    vel: Optional[np.ndarray] = None          # (K, 4, 3) km/s, MW rest frame

    @property
    def positions(self):
        return self.cloud[["bc0_x", "bc0_y", "bc0_z"]].to_numpy()

    @property
    def velocities(self):
        return self.cloud[["bc0_vx", "bc0_vy", "bc0_vz"]].to_numpy()

    @property
    def median_position(self):
        return np.median(self.positions, axis=0)

    @property
    def median_velocity(self):
        return np.median(self.velocities, axis=0)

    @property
    def covariance(self):
        return np.cov(self.positions.T)

    @property
    def principal_axes(self):
        """(semiaxes 1-sigma, descending; unit vectors as columns)."""
        w, U = np.linalg.eigh(self.covariance)
        o = np.argsort(w)[::-1]
        return np.sqrt(w[o]), U[:, o]

    @property
    def axes(self):
        return self.principal_axes[0]

    def vcmb(self, membership="4"):
        """(16th, 50th, 84th percentile) of the barycentre speed w.r.t. the CMB, km/s."""
        return tuple(np.percentile(self.cloud[f"vcmb_{membership}_kms"], [16, 50, 84]))

    def relative_to_barycentre(self, body, membership=None):
        """Position and velocity of ``body`` relative to the barycentre of ``membership``.

        One row per realization. ``body`` is one of ``NAMES`` (``"MilkyWay"``, ``"M31"``,
        ``"M33"``, ``"LMC"``); ``membership`` defaults to the one of this result and may be
        any of ``"2"``, ``"3M"``, ``"3L"``, ``"4"``. ``body`` does not need to be a member.

        Columns: ``distance_kpc``, ``speed_kms``, and the radial (``v_rad_kms``, positive
        when receding) and tangential (``v_tan_kms``) components with respect to the barycentre.
        """
        if self.pos is None:
            raise ValueError("this Result does not carry the member phase space")
        j = NAMES.index(body)
        rb, vb = _bc(self.masses, self.pos, self.vel, MEMBERSHIPS[membership or self.membership])
        dr, dv = self.pos[:, j] - rb, self.vel[:, j] - vb
        d = np.linalg.norm(dr, axis=1)
        sp = np.linalg.norm(dv, axis=1)
        vr = (dr * dv).sum(1) / d
        return pd.DataFrame({"distance_kpc": d, "speed_kms": sp, "v_rad_kms": vr,
                             "v_tan_kms": np.sqrt(np.maximum(sp**2 - vr**2, 0.0))})

    def body_budget(self, body, column="speed_kms"):
        """Median and 16th/84th percentiles of ``column`` for every membership definition."""
        rows = {}
        for tag in MEMBERSHIPS:
            q = np.percentile(self.relative_to_barycentre(body, tag)[column], [16, 50, 84])
            rows[tag] = dict(p16=q[0], median=q[1], p84=q[2])
        return pd.DataFrame(rows).T

    def summary(self):
        a = self.axes
        c = self.cloud
        q = self.vcmb(self.membership)
        return "\n".join([
            f"Local Group barycentre | fiducial: {self.fiducial} | members: {self.membership} "
            f"| N = {self.n}, seed = {self.seed}",
            f"  median position (Galactocentric, kpc): {np.round(self.median_position, 1)}",
            f"  median velocity (MW rest frame, km/s): {np.round(self.median_velocity, 1)}",
            f"  1-sigma principal axes (kpc)          : {np.round(a, 1)}",
            f"  median ICRS: RA = {c.bc0_ra.median():.1f} deg, Dec = {c.bc0_dec.median():+.1f} deg, "
            f"d = {c.bc0_dist.median():.0f} kpc",
            f"  speed w.r.t. the CMB (16/50/84 %)      : {np.round(q, 1)} km/s",
            "Assumptions:",
            "  * Gaussian sampling truncated at +-2 sigma (effective std ~0.88 sigma);",
            "  * sky positions held fixed; MW at rest at the origin of the Galactocentric frame;",
            "  * the output is the uncertainty INDUCED by the declared errors, not a complete",
            "    physical uncertainty region.",
        ])

    def __repr__(self):
        return self.summary()


def barycentre(source, n=1000, seed=2000, sigmas: Optional[dict] = None, membership="4"):
    """Barycentre of the Local Group with Monte Carlo uncertainty.

    Parameters
    ----------
    source : Fiducial, str or dict
        A fiducial (``fiducials.WU``), its name, or a dict of heliocentric observables
        with the structure of ``Fiducial.obs``.
    n : number of Monte Carlo realizations.
    seed : realization ``i`` uses ``default_rng([seed, i])`` (reproducible one by one).
    sigmas : 1-sigma uncertainties; defaults to Table 1 of Wu et al. (2026).
    membership : ``"2"`` (MW+M31), ``"3M"`` (+M33), ``"3L"`` (+LMC) or ``"4"`` (all).
    """
    name, base, zero = _resolve(source)
    sig = DEFAULT_SIGMAS if sigmas is None else sigmas
    obs_list = sample_ensemble(base, n, seed, sig)
    m, r, v = phase_space(obs_list, zero)
    df = build_cloud(m, r, v, membership)
    m0, r0, v0 = phase_space([base], zero)
    (rn, vn) = _bc(m0, r0, v0, MEMBERSHIPS[membership])
    return Result(df, rn[0], vn[0], name, n, seed, membership, sig, m, r, v)
