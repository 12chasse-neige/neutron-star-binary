"""Physical identities and independent limiting solutions for the new dynamics."""
import json
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from model import Binary, build_reservoir

ROOT = Path(__file__).resolve().parent


class PhysicalLimits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = json.loads((ROOT/"config.json").read_text())
        cls.res = build_reservoir(cls.cfg, ROOT.parent/"results")

    def test_two_independent_inertia_integrals(self):
        self.assertLess(self.res.metadata["cylindrical_vs_spherical_Iv_relative_error"], 3e-6)

    def test_isolated_lag_energy_and_angular_momentum(self):
        b = Binary(self.cfg, self.res, torque_enabled=False, gw_enabled=False)
        y0 = np.array([5., .2, 6., 6.001, 0., 0., 0., 0., 0.])
        tau = b.fp/(2*b.B*y0[3])
        s = solve_ivp(b.rhs_time, (0, 6*tau), y0, rtol=2e-11, atol=2e-12,
                      dense_output=True)
        y = s.y
        j = b.fp*y[2]+b.fv*y[3]
        energy = .5*(b.fp*y[2]**2+b.fv*y[3]**2)+y[4]
        self.assertLess(np.max(abs(j-j[0])), 1e-11)
        self.assertLess(np.max(abs(energy-energy[0])), 2e-10)
        exact_heat = .5*b.fp*b.fv*((y0[3]-y0[2])**2-(y[3,-1]-y[2,-1])**2)
        self.assertAlmostEqual(y[4,-1]/exact_heat, 1, delta=1e-5)
        # Independent linearization of the lag relaxation, for a small initial lag.
        lag_at_tau = s.sol(tau)[3]-s.sol(tau)[2]
        self.assertAlmostEqual(lag_at_tau/(.001/np.e), 1, delta=3e-4)

    def test_closed_hamiltonian_with_quadrupole(self):
        b = Binary(self.cfg, self.res, gw_enabled=False)
        y0 = np.array([8., -.5, 8.02, 8., 0., 0., 0., 0., 0.])
        s = solve_ivp(b.rhs_time, (0, 500), y0, rtol=1e-11, atol=1e-11, max_step=.5)
        n, phase, op, ov, heat = s.y[:5]
        ell, u, *_ = b.quantities(n, phase, op, ov)
        energy = -.5*ell*n+u+.5*(b.fp*op**2+b.fv*ov**2)+heat
        j = ell+b.fp*op+b.fv*ov
        self.assertLess(np.ptp(energy), 2e-8)
        self.assertLess(np.ptp(j), 2e-9)

    def test_point_particle_chirp(self):
        b = Binary(self.cfg, self.res, torque_enabled=False)
        d, s = b.evolve("paper", sample_count=300)
        n0, n1 = 2*np.pi*np.array([self.cfg["initial_orbital_frequency_hz"], self.cfg["final_orbital_frequency_hz"]])
        expected = b.ell_coeff/(8*b.flux_coeff)*(n0**(-8/3)-n1**(-8/3))
        self.assertAlmostEqual(s["duration_s"]/expected, 1., delta=1e-9)
        self.assertLess(abs(s["final_gw_phase_difference_from_point_particle_rad"]), 1e-5)
        self.assertAlmostEqual(s["heat_J"], 0.)

    def test_no_drag_and_infinite_drag_limits(self):
        vals = []
        for drag in [0., 1., 1e12]:
            cfg = {**self.cfg, "core_drag_ratio": drag, "inner_crust_drag_ratio": drag}
            vals.append(build_reservoir(cfg, ROOT.parent/"results").B)
        self.assertEqual(vals[0], 0.)
        self.assertGreater(vals[1], .1)
        self.assertLess(vals[2], 1e-10)


if __name__ == "__main__":
    unittest.main()
