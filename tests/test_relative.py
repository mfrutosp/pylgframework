import numpy as np
import pytest

from pylgframework import barycentre, compare, fiducials


def test_milky_way_speed_relative_to_barycentre_is_the_barycentre_speed():
    """MW is at rest at the origin, so |v_MW - v_BC| = |v_BC| (exact identity)."""
    r = barycentre(fiducials.WU, n=50, seed=1)
    sp = r.relative_to_barycentre("MilkyWay", "4")["speed_kms"].to_numpy()
    np.testing.assert_allclose(sp, np.linalg.norm(r.velocities, axis=1), rtol=1e-12)


def test_components_are_consistent():
    r = barycentre(fiducials.WU, n=50, seed=1)
    d = r.relative_to_barycentre("LMC", "2")
    np.testing.assert_allclose(d.v_rad_kms**2 + d.v_tan_kms**2, d.speed_kms**2, rtol=1e-9)


def test_body_budget_has_all_memberships():
    t = barycentre(fiducials.WU, n=50, seed=1).body_budget("LMC")
    assert list(t.index) == ["2", "3M", "3L", "4"]
    assert (t.p16 <= t["median"]).all() and (t["median"] <= t.p84).all()


@pytest.mark.slow
def test_lmc_speed_relative_to_the_barycentre():
    """Worked example of the manuscript (Sect. 5): LMC speed w.r.t. the barycentre."""
    t = barycentre(fiducials.WU, n=1000, seed=2000).body_budget("LMC")
    assert abs(t.loc["2", "median"] - 271) < 2 and abs(t.loc["4", "median"] - 253) < 2


@pytest.mark.slow
def test_lmc_speed_paired_fiducial_offset():
    c = compare(fiducials.WU, fiducials.VDM2012, n=1000, seed=3000)
    lo, med, hi = c.body_offset("LMC", "4")
    assert 15 < med < 25 and lo < med < hi
