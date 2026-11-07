"""Literature fiducials with explicit provenance.

Each fiducial stores HELIOCENTRIC observables in the format used by the
sampler (``mu`` distance modulus, ``pmra``/``pmdec`` in uas/yr, ``vlos`` in
km/s, ``m`` = M200 in Msun, ``c`` = NFW concentration). RA/Dec are the usual
catalogued centres and are not sampled.

Usage::

    from pylgframework import fiducials
    fiducials.WU
    fiducials.get("vdm2012")
"""
from dataclasses import dataclass, field
from typing import Dict, Tuple

_LMC = dict(m=1.5e11, c=10.0, ra=80.894, dec=-69.756, mu=18.477,
            pmra=1910.0, pmdec=229.0, vlos=262.2)


@dataclass(frozen=True)
class Fiducial:
    name: str
    obs: Dict[str, Dict[str, float]]
    reference: str
    notes: str = ""
    zero_mass: Tuple[str, ...] = field(default_factory=tuple)   # bodies removed AFTER sampling

    def __repr__(self):
        return f"Fiducial({self.name!r}, {self.reference})"


WU = Fiducial(
    name="wu2026",
    reference="Wu, Huang, Zhang, Sun & Shao (2026), ApJ 1001, L19, Table 1, "
              "doi:10.3847/2041-8213/ae5799",
    notes="Follows the framework of Sawala et al. (2025).",
    obs={
        "MilkyWay": dict(m=1.0e12, c=10.0),
        "M31": dict(m=1.3e12, c=15.0, ra=10.6847, dec=41.2690, mu=24.407,
                    pmra=45.9, pmdec=-20.5, vlos=-301.0),
        "M33": dict(m=3.0e11, c=10.0, ra=23.4621, dec=30.6602, mu=24.670,
                    pmra=30.6, pmdec=19.0, vlos=-179.2),
        "LMC": dict(_LMC),
    })

SAWALA = Fiducial(
    name="sawala2025",
    reference="Sawala et al. (2025), Nat. Astron. 9, 1206, doi:10.1038/s41550-025-02563-1 "
              "(fiducial of the companion notebook)",
    notes="Masses converted from the 1e10 Msun units of the notebook.",
    obs={
        "MilkyWay": dict(m=1.0e12, c=10.0),
        "M31": dict(m=1.3e12, c=10.0, ra=10.6847, dec=41.2690, mu=24.407,
                    pmra=48.9, pmdec=-36.9, vlos=-301.0),
        "M33": dict(m=3.0e11, c=10.0, ra=23.4621, dec=30.6602, mu=24.67,
                    pmra=31.0, pmdec=-29.0, vlos=-179.2),
        "LMC": dict(_LMC),
    })

VDM2019 = Fiducial(
    name="vdM2019",
    reference="Sawala et al. (2025) fiducial with the M31 proper motion of "
              "van der Marel et al. (2019), ApJ 872, 24, doi:10.3847/1538-4357/ab001b",
    notes="Identical to SAWALA except for the M31 proper motion (49, -38) +- 11 uas/yr, "
          "the combined Gaia DR2 + HST value of that paper.",
    obs={
        "MilkyWay": dict(m=1.0e12, c=10.0),
        "M31": dict(m=1.3e12, c=10.0, ra=10.6847, dec=41.2690, mu=24.407,
                    pmra=49.0, pmdec=-38.0, vlos=-301.0),
        "M33": dict(m=3.0e11, c=10.0, ra=23.4621, dec=30.6602, mu=24.67,
                    pmra=31.0, pmdec=-29.0, vlos=-179.2),
        "LMC": dict(_LMC),
    })

_VDM2012_OBS = {
    "MilkyWay": dict(m=1.5e12, c=10.0),                      # Paper III canonical model
    # M31: (v_W, v_N) = (-125.2 +- 30.8, -73.8 +- 28.4) km/s at D = 770 kpc
    #      -> pmra = +34.3, pmdec = -20.2 uas/yr
    "M31": dict(m=1.5e12, c=10.0, ra=10.6847, dec=41.2690, mu=24.43,
                pmra=34.3, pmdec=-20.2, vlos=-301.0),
    # M33: Paper II Sec. 6.1 Galactocentric velocity -> pmra = +23.2, pmdec = +7.5 uas/yr
    #      at D = 794 kpc; matches the VLBA maser PM of Brunthaler et al. (2005)
    "M33": dict(m=1.5e11, c=10.0, ra=23.4621, dec=30.6602, mu=24.50,
                pmra=23.2, pmdec=7.5, vlos=-180.0),
    "LMC": dict(_LMC),                                       # NOT in vdM2012 (see notes)
}

VDM2012 = Fiducial(
    name="vdm2012",
    reference="van der Marel et al. (2012), ApJ 753, 8 (doi:10.1088/0004-637X/753/1/8) "
              "and 9 (doi:10.1088/0004-637X/753/1/9), Papers II and III",
    notes="Four-body variant. vdM2012 did NOT include the LMC; its entry is inherited "
          "from WU for comparability.",
    obs=_VDM2012_OBS)

VDM2012_3B = Fiducial(
    name="vdm2012_3b",
    reference="van der Marel et al. (2012), ApJ 753, 8 (doi:10.1088/0004-637X/753/1/8) "
              "and 9 (doi:10.1088/0004-637X/753/1/9), Papers II and III",
    notes="Faithful three-body variant (MW+M31+M33). The LMC is sampled (so draws stay "
          "paired with the other fiducials) and its mass is set to zero afterwards.",
    obs=_VDM2012_OBS, zero_mass=("LMC",))

ALL = {f.name: f for f in (WU, SAWALA, VDM2019, VDM2012, VDM2012_3B)}


def get(name):
    """Fiducial by name (``wu2026``, ``sawala2025``, ``vdM2019``, ``vdm2012``, ``vdm2012_3b``)."""
    try:
        return ALL[name]
    except KeyError:
        raise KeyError(f"unknown fiducial {name!r}; available: {sorted(ALL)}") from None
