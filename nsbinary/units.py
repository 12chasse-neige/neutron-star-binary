"""Physical constants and conversions between SI and geometrized units.

The structure solver integrates with G = c = 1 and lengths in kilometres.
Tabulated equations of state and published observables use SI or cgs, so
conversions happen only when data enter or leave the solver.
"""

import lal

G = lal.G_SI                          # m^3 kg^-1 s^-2
C = lal.C_SI                          # m s^-1
MSUN = lal.MSUN_SI                    # kg
RSUN_KM = lal.MRSUN_SI / 1000         # km

PRESSURE_TO_GEOM = G / C**4 * 1e6    # Pa -> km^-2
DENSITY_TO_GEOM = G / C**2 * 1e6     # kg m^-3 -> km^-2
