"""Paired comparison between two fiducials (identical random draws)."""
from dataclasses import dataclass

import numpy as np

from .barycentre import Result, barycentre
from .constants import MEMBERSHIPS


@dataclass
class Comparison:
    a: Result
    b: Result
    offset_kpc: float                 # |nominal_b - nominal_a|
    offset_principal: np.ndarray      # offset in the principal-axis frame of ``a`` (kpc)
    offset_sigma: float               # |component along the major axis of a| / a1
    dvcmb: dict                       # median v_CMB(b) - median v_CMB(a), per membership (km/s)

    def body_offset(self, body, membership=None, column="speed_kms"):
        """Paired difference (b - a) of a body's motion relative to the barycentre.

        Realization i of both fiducials uses the same random draw, so the difference is
        free of Monte Carlo noise. Returns the (16th, 50th, 84th) percentiles.
        """
        a = self.a.relative_to_barycentre(body, membership)[column].to_numpy()
        b = self.b.relative_to_barycentre(body, membership)[column].to_numpy()
        return tuple(np.percentile(b - a, [16, 50, 84]))

    def __repr__(self):
        d = ", ".join(f"{k}: {v:+.0f}" for k, v in self.dvcmb.items())
        return (f"Comparison({self.a.fiducial} -> {self.b.fiducial}): offset = {self.offset_kpc:.1f} kpc "
                f"({self.offset_sigma:.1f} sigma along the major axis of {self.a.fiducial}); "
                f"dv_CMB [km/s] = {d}")


def compare(a, b, n=1000, seed=3000, sigmas=None, membership="4"):
    """Paired experiment: realization i of both fiducials uses the same random draw,
    so differences measure the effect of the adopted central values, not Monte Carlo noise.
    """
    ra = barycentre(a, n=n, seed=seed, sigmas=sigmas, membership=membership)
    rb = barycentre(b, n=n, seed=seed, sigmas=sigmas, membership=membership)
    d = rb.nominal_position - ra.nominal_position
    axes, U = ra.principal_axes
    proj = d @ U
    dv = {t: float(rb.cloud[f"vcmb_{t}_kms"].median() - ra.cloud[f"vcmb_{t}_kms"].median())
          for t in MEMBERSHIPS}
    return Comparison(ra, rb, float(np.linalg.norm(d)), proj, float(abs(proj[0]) / axes[0]), dv)
