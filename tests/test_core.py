import numpy as np
import pytest

from pylgframework import barycentre, compare, fiducials
from pylgframework.frames import galactocentric_to_observables, helio_to_galactocentric
from pylgframework.sampling import sample_ensemble


def test_realizations_are_reproducible_one_by_one():
    a = sample_ensemble(fiducials.WU.obs, 5, seed=7)
    b = sample_ensemble(fiducials.WU.obs, 8, seed=7)
    assert a == b[:5]                      # realization i depends only on (seed, i)


def test_truncation_is_respected():
    base = fiducials.WU.obs
    from pylgframework.sampling import DEFAULT_SIGMAS
    for o in sample_ensemble(base, 200, seed=1):
        for g, sig in DEFAULT_SIGMAS.items():
            for k, s in sig.items():
                assert abs(o[g][k] - base[g][k]) <= 2.0 * s + 1e-9


def test_sky_positions_are_fixed():
    base = fiducials.WU.obs
    for o in sample_ensemble(base, 20, seed=1):
        for g in ("M31", "M33", "LMC"):
            assert o[g]["ra"] == base[g]["ra"] and o[g]["dec"] == base[g]["dec"]


def test_transform_roundtrip():
    o = fiducials.WU.obs["M31"]
    r, v = helio_to_galactocentric(o["ra"], o["dec"], o["mu"], o["pmra"], o["pmdec"], o["vlos"])
    back = galactocentric_to_observables(r, v)
    for k in ("ra", "dec", "mu", "pmra", "pmdec", "vlos"):
        assert back[k][0] == pytest.approx(o[k], rel=1e-6, abs=1e-6)


def test_three_body_variant_degenerates_exactly():
    c = barycentre("vdm2012_3b", n=50, seed=3000).cloud
    assert (c.shift_LMC_kpc == 0).all() and (c.dv_LMC_kms == 0).all()
    np.testing.assert_allclose(c.shift_2to4_kpc, c.shift_M33_kpc)
    np.testing.assert_allclose(c.vcmb_4_kms, c.vcmb_3M_kms)


def test_paired_comparison_uses_identical_draws():
    # Same fiducial against itself: zero offset, zero velocity difference.
    c = compare(fiducials.WU, fiducials.WU, n=50, seed=3000)
    assert c.offset_kpc == pytest.approx(0.0, abs=1e-9)
    assert all(abs(v) < 1e-9 for v in c.dvcmb.values())


def test_custom_sigmas_change_the_cloud():
    from pylgframework import DEFAULT_SIGMAS
    tight = {g: {k: v / 10 for k, v in s.items()} for g, s in DEFAULT_SIGMAS.items()}
    a = barycentre(fiducials.WU, n=200, seed=1)
    b = barycentre(fiducials.WU, n=200, seed=1, sigmas=tight)
    assert b.axes[0] < 0.2 * a.axes[0]


def test_summary_states_the_assumptions():
    s = barycentre(fiducials.WU, n=20, seed=1).summary()
    assert "truncated" in s and "INDUCED" in s


def test_unknown_fiducial():
    with pytest.raises(KeyError):
        fiducials.get("nope")
