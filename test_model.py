"""Physics checks beyond comparing the model to its own implementation."""
import unittest
import numpy as np
from ns_model import integrate, RSUN_KM


class UniformDensityEOS:
    """Incompressible analytic TOV benchmark, not an astrophysical EOS."""
    hmax = 2.0
    energy_geom = 8e-4  # km^-2

    def at_h(self, h):
        return self.energy_geom*np.expm1(h), self.energy_geom, 1.0

    def dedp(self, pressure):
        return 0.0


class AnalyticTOVTests(unittest.TestCase):
    def test_uniform_density_schwarzschild_interior(self):
        eos = UniformDensityEOS()
        compactness = 0.15
        surface_lapse = np.sqrt(1-2*compactness)
        pc_over_e = (1-surface_lapse)/(3*surface_lapse-1)
        star = integrate(eos, np.log1p(pc_over_e), surface_h=1e-12, rtol=1e-10)
        expected_r = np.sqrt(3*compactness/(4*np.pi*eos.energy_geom))
        self.assertAlmostEqual(star.radius_km/expected_r, 1, delta=2e-8)
        self.assertAlmostEqual(star.mass_msun*RSUN_KM/(compactness*expected_r), 1, delta=2e-8)
        r = np.linspace(0.01, 0.99, 50)*expected_r
        lapse = np.sqrt(1-2*compactness*(r/expected_r)**2)
        expected_p_over_e = (lapse-surface_lapse)/(3*surface_lapse-lapse)
        calculated = np.expm1(star.solution.sol(r)[1])
        np.testing.assert_allclose(calculated, expected_p_over_e, rtol=2e-7, atol=1e-10)

    def test_hartle_newtonian_uniform_density_limit(self):
        eos = UniformDensityEOS()
        compactness = 1e-4
        lapse = np.sqrt(1-2*compactness)
        hc = np.log1p((1-lapse)/(3*lapse-1))
        star = integrate(eos, hc, surface_h=1e-14, rtol=1e-10, max_step_km=0.001)
        ratio = star.inertia_geom_km3/(star.mass_msun*RSUN_KM*star.radius_km**2)
        self.assertAlmostEqual(ratio, 0.4, delta=1e-4)

    def test_invalid_enthalpy_rejected(self):
        with self.assertRaises(ValueError):
            integrate(UniformDensityEOS(), -1)


if __name__ == "__main__":
    unittest.main()
