"""Cold tabulated barotrope adapted to the kilometre-based structure equations.

LAL reconstructs pseudo-enthalpy and baryonic rest density from a two-column
P(E) table. The APR table is coarse enough that its native cubic interpolation
produces a small negative-compressibility interval near the inner crust, so the
frozen table is resampled monotonically in ln E versus ln P before use.
"""

from pathlib import Path

import lalsimulation as lalsim
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from .units import PRESSURE_TO_GEOM

DENSE_TABLE_POINTS = 24000


class EquationOfState:
    """Adapt LAL thermodynamics to energies and pressures in km^-2.

    ``table`` selects a tabulated barotrope by path; ``name`` selects a LAL
    built-in EOS. This adapter supplies no particle fractions or phase labels.
    """

    def __init__(self, name="APR", table=None):
        self.name = name
        self.raw = (lalsim.SimNeutronStarEOSFromFile(str(Path(table).resolve()))
                    if table else lalsim.SimNeutronStarEOSByName(name))
        self.hmax = lalsim.SimNeutronStarEOSMaxPseudoEnthalpy(self.raw)

    def at_h(self, pseudo_enthalpy):
        """Return (P [km^-2], total energy density [km^-2], rest rho [kg m^-3])."""
        # Trial ODE evaluations can step beyond the terminal surface. Keep the
        # library query positive; the event, not this floor, defines the surface.
        pseudo_enthalpy = max(float(pseudo_enthalpy), 1e-14)
        pressure = lalsim.SimNeutronStarEOSPressureOfPseudoEnthalpy(pseudo_enthalpy, self.raw)
        energy_density = lalsim.SimNeutronStarEOSEnergyDensityOfPseudoEnthalpy(pseudo_enthalpy, self.raw)
        rest_mass_density = lalsim.SimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(pseudo_enthalpy, self.raw)
        return pressure * PRESSURE_TO_GEOM, energy_density * PRESSURE_TO_GEOM, rest_mass_density

    def h_at_rho(self, rest_mass_density_g_cm3):
        """Invert the rest-density label, converting kg/m^3 to g/cm^3."""
        return brentq(lambda pseudo_enthalpy: self.at_h(pseudo_enthalpy)[2] / 1000 - rest_mass_density_g_cm3,
                      1e-10, self.hmax, xtol=1e-14)

    def dedp(self, pressure_geom):
        """Dimensionless dE/dP; its reciprocal is c_s^2/c^2."""
        return lalsim.SimNeutronStarEOSEnergyDensityDerivOfPressure(
            pressure_geom / PRESSURE_TO_GEOM, self.raw)


def monotone_table(raw, path, points=DENSE_TABLE_POINTS):
    """Resample log E(log P) with PCHIP, retaining every original node exactly.

    Both input and output columns stay in geometric m^-2 as expected by LAL.
    ``points`` counts the added logarithmic grid, not the final union of nodes.
    """
    dense_pressure = np.unique(np.r_[np.geomspace(raw[0, 0], raw[-1, 0], points), raw[:, 0]])
    pchip = PchipInterpolator(np.log(raw[:, 0]), np.log(raw[:, 1]))
    dense_energy_density = np.exp(pchip(np.log(dense_pressure)))
    dense_energy_density[np.searchsorted(dense_pressure, raw[:, 0])] = raw[:, 1]
    np.savetxt(path, np.column_stack([dense_pressure, dense_energy_density]), fmt="%.18e")
    return len(dense_pressure)


def load_apr(data_dir):
    """Load the frozen APR barotrope, rebuilding its dense table if needed.

    The dense table is a deterministic resampling of the original table; it is
    rebuilt only when absent so normal runs read saved local data.
    """
    data_dir = Path(data_dir)
    original = data_dir / "LALSimNeutronStarEOS_APR.dat"
    dense = data_dir / "APR_monotone_dense.dat"
    if not original.exists():
        raise FileNotFoundError(f"Missing APR table: {original}")
    if not dense.exists():
        raw = np.loadtxt(original)
        if not np.all(np.diff(raw, axis=0) > 0):
            raise ValueError("Original EOS table must be strictly monotone.")
        monotone_table(raw, dense)
    return EquationOfState("APR", dense)
