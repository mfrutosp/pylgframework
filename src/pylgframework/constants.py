"""Single source of truth for the constants shared by the whole package."""
import numpy as np

KMS = 1.022712165                      # 1 km/s in kpc/Gyr

# Solar Galactocentric frame (Sawala et al. 2025; Wu et al. 2026)
SUN_R0 = 8.122                         # kpc
SUN_V = (12.9, 245.6, 7.78)            # km/s

# Truncation of the Gaussian sampling, in sigmas
TRUNC = 2.0

# Solar dipole w.r.t. the CMB rest frame:
# Planck 2015 results. V. LFI calibration (Planck Collaboration 2016, A&A 594, A5).
V_SUN_CMB = 369.82                     # km/s
L_CMB, B_CMB = 264.021, 48.253         # deg (Galactic l, b)

_U_SUN_CMB = V_SUN_CMB * np.array([
    np.cos(np.deg2rad(B_CMB)) * np.cos(np.deg2rad(L_CMB)),
    np.cos(np.deg2rad(B_CMB)) * np.sin(np.deg2rad(L_CMB)),
    np.sin(np.deg2rad(B_CMB))])
# Velocity of the MW (origin of the Galactocentric frame) w.r.t. the CMB
V_MW_CMB = _U_SUN_CMB - np.array(SUN_V)

# Body order is part of the public contract (indices used everywhere)
NAMES = ("MilkyWay", "M31", "M33", "LMC")

# Membership definitions: tag -> indices into NAMES
MEMBERSHIPS = {
    "2": (0, 1),            # MW + M31
    "3M": (0, 1, 2),        # MW + M31 + M33
    "3L": (0, 1, 3),        # MW + M31 + LMC
    "4": (0, 1, 2, 3),      # MW + M31 + M33 + LMC
}
