"""Regression against the values published in Frutos Plaza (2026), Table 1 and Sect. 3.

Tolerances reflect the rounding of the published numbers (and astropy version differences).
"""
import numpy as np
import pytest

from pylgframework import barycentre, compare, fiducials

pytestmark = pytest.mark.slow

# fiducial: (a1, a2, a3, D_BC, RA, Dec, vCMB(2), vCMB(4), Delta_2to4)   [Table 1, seed 3000]
TABLE1 = {
    "wu2026":     (53.6, 7.3, 0.2, 448, 13.4, 38.8, 621, 629, 34),
    "sawala2025": (53.6, 7.3, 0.2, 448, 13.4, 38.8, 619, 627, 34),
    "vdM2019":    (53.6, 7.3, 0.2, 448, 13.4, 38.8, 619, 626, 34),
    "vdm2012":    (48.8, 5.3, 0.2, 381, 11.6, 39.9, 604, 614, 15),
    "vdm2012_3b": (49.6, 5.5, 0.0, 401, 11.4, 40.2, 604, 602, 21),
}


@pytest.mark.parametrize("name", sorted(TABLE1))
def test_table1(name):
    a1, a2, a3, d, ra, dec, v2, v4, d24 = TABLE1[name]
    r = barycentre(name, n=1000, seed=3000)
    c = r.cloud
    np.testing.assert_allclose(r.axes, [a1, a2, a3], atol=0.15)
    assert abs(c.bc0_dist.median() - d) < 1.0
    assert abs(c.bc0_ra.median() - ra) < 0.1
    assert abs(c.bc0_dec.median() - dec) < 0.1
    assert abs(c.vcmb_2_kms.median() - v2) < 1.0
    assert abs(c.vcmb_4_kms.median() - v4) < 1.0
    assert abs(c.shift_2to4_kpc.median() - d24) < 1.0


def test_main_text_ensemble():
    """Sect. 3.1 and 3.3 (seed 2000)."""
    r = barycentre(fiducials.WU, n=1000, seed=2000)
    np.testing.assert_allclose(r.axes, [54.8, 7.7, 0.2], atol=0.15)
    assert abs(np.linalg.norm(r.median_position) - 452) < 1.5
    assert abs(r.cloud.vcmb_2_kms.median() - 622) < 1.0
    assert abs(r.cloud.vcmb_4_kms.median() - 630) < 1.0
    assert abs(r.cloud.shift_M33_kpc.median() - 53) < 1.0
    assert abs(r.cloud.shift_LMC_kpc.median() - 27) < 1.0


def test_fiducial_audit():
    """Sect. 4: 68 kpc (1.2 sigma) and 17 km/s between the modern and the HST-era fiducial."""
    c = compare(fiducials.WU, fiducials.VDM2012, n=1000, seed=3000)
    assert abs(c.offset_kpc - 67.8) < 0.5
    assert abs(c.offset_sigma - 1.2) < 0.1
    assert abs(c.dvcmb["2"] + 17) < 1.0
