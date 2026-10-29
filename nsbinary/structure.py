"""Cold spherical neutron star: TOV, Hartle frame dragging and tidal Love number.

The integration uses G = c = 1 and lengths in km; public results use SI/cgs.
An equation of state supplies the background; density thresholds only label its
radial layers.

The integration state is [m, h, w, z, y, ell], where w = Omega - omega_fd,
z = dw/dr, y = r H'/H and ell is proper radial distance. The code's h is
pseudo-enthalpy, not Planck's constant.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .units import DENSITY_TO_GEOM, PRESSURE_TO_GEOM, RSUN_KM, C, G

LAYER_NAMES = ("core", "inner crust", "outer crust")


@dataclass
class Star:
    """Dense ODE solution and surface observables, with the original EOS.

    Gravitational mass is stored as a length. Densities and boundary pressures
    are converted only at the public export boundary, never inside the ODE.
    """
    eos: object
    central_h: float
    surface_h: float
    solution: object
    radius_km: float

    @property
    def surface(self):
        return self.solution.sol(self.radius_km)

    @property
    def mass_msun(self):
        return self.surface[0] / RSUN_KM

    @property
    def inertia_geom_km3(self):
        """Match w(R), w'(R) to vacuum w=Omega-2J/r^3; return J/Omega."""
        _, _, w, z, _, _ = self.surface
        r = self.radius_km
        return r**4 * z / (6 * (w + r * z / 3))

    @property
    def love_k2(self):
        """Exterior l=2 matching for a smooth EOS without density jumps.

        ``c`` below is compactness, not the SI speed of light ``C``. A finite
        surface density or internal discontinuity needs extra jump conditions.
        """
        compactness = self.surface[0] / self.radius_km
        y = self.surface[4]
        c = compactness
        numerator = 8 / 5 * c**5 * (1 - 2*c)**2 * (2 + 2*c*(y - 1) - y)
        denominator = (2*c*(6 - 3*y + 3*c*(5*y - 8))
                       + 4*c**3*(13 - 11*y + c*(3*y - 2) + 2*c*c*(1 + y))
                       + 3*(1 - 2*c)**2*(2 - y + 2*c*(y - 1))*np.log1p(-2*c))
        return numerator / denominator

    @property
    def tidal_lambda(self):
        return 2 / 3 * self.love_k2 / (self.surface[0] / self.radius_km)**5

    def boundary(self, rho_g_cm3):
        """Locate an assumed density threshold, not a computed phase change.

        Invert rho_b(h) first and then h(r), using the dense solution so the
        result does not depend on whether an export grid hits the boundary.
        """
        h_target = self.eos.h_at_rho(rho_g_cm3)
        if not self.surface_h < h_target < self.central_h:
            raise ValueError("Layer threshold must lie inside the star.")
        r = brentq(lambda r: self.solution.sol(r)[1] - h_target,
                   self.solution.t[0], self.radius_km, xtol=1e-10)
        p, energy, rho = self.eos.at_h(h_target)
        state = self.solution.sol(r)
        return {"radius_km": float(r), "enclosed_mass_msun": float(state[0] / RSUN_KM),
                "pressure_pa": float(p / PRESSURE_TO_GEOM),
                "rest_mass_density_g_cm3": float(rho / 1000),
                "energy_equivalent_density_g_cm3": float(energy / DENSITY_TO_GEOM / 1000)}

    def profile(self, n=2401, extra_radii=()):
        """Export SI/cgs profiles with extra sampling in the thin crust."""
        # A second dense grid resolves the thin crust; include boundaries exactly.
        r = np.unique(np.r_[0, np.linspace(self.solution.t[0], self.radius_km, n),
                            np.linspace(0.85*self.radius_km, self.radius_km, n//3),
                            extra_radii])
        state = self.solution.sol(np.maximum(r, self.solution.t[0]))
        # Dense output begins at r0. Insert exact regular central values rather
        # than extrapolating the numerical solution into its singular endpoint.
        state[:, 0] = [0, self.central_h, 1, 0, 2, 0]
        therm = np.array([self.eos.at_h(h) for h in state[1]])
        p, energy, rho = therm.T
        dm_dr = 4 * np.pi * r**2 * energy / RSUN_KM
        return {"radius_km": r, "enclosed_mass_msun": state[0] / RSUN_KM,
                "pressure_pa": p / PRESSURE_TO_GEOM,
                "rest_mass_density_g_cm3": rho / 1000,
                "energy_equivalent_density_g_cm3": energy / DENSITY_TO_GEOM / 1000,
                "mass_per_radius_msun_per_km": dm_dr,
                "proper_radius_km": state[5], "pseudo_enthalpy": state[1]}


def integrate(eos, central_h, surface_h=1e-10, rtol=2e-9,
              max_step_km=0.025, start_radius_km=1e-4):
    """Integrate m, h, w, w', y and proper radial distance in areal radius r.

    w = Omega - omega is arbitrarily normalized at the centre. Matching to the
    exterior w = Omega - 2J/r^3 determines I = J/Omega. y = r H'/H is the l=2
    static tidal metric perturbation; no density jump is imposed here.
    """
    if not 0 < surface_h < central_h < eos.hmax:
        raise ValueError("Require 0 < surface_h < central_h < EOS maximum h.")
    pc, ec, _ = eos.at_h(central_h)
    r0 = start_radius_km
    # Leading regular-centre series avoids 1/r singularities in all six ODEs.
    # The Hartle amplitude is arbitrary; tidal y starts at its l=2 value.
    initial = [4*np.pi/3 * ec*r0**3,
               central_h - 2*np.pi*(ec/3 + pc)*r0**2,
               1 + 8*np.pi/5*(ec + pc)*r0**2,
               16*np.pi/5*(ec + pc)*r0, 2, r0]

    def rhs(r, state):
        # a is the positive gravitational acceleration in enthalpy units:
        # dh/dr = -a. log_j_prime couples the TOV solution to Hartle rotation.
        m, h, w, z, y, _ = state
        p, e, _ = eos.at_h(max(h, surface_h * 0.01))
        f = 1 - 2*m/r
        if f <= 0:
            raise RuntimeError("Trapped surface encountered.")
        a = (m + 4*np.pi*r**3*p) / (r*r*f)
        log_j_prime = -4*np.pi*r*(e + p)/f
        tidal_f = (1 - 4*np.pi*r*r*(e - p))/f
        tidal_q = (4*np.pi*(5*e + 9*p + (e+p)*eos.dedp(p))/f
                   - 6/(r*r*f) - 4*a*a)
        # Ordered derivatives: m, h, w, z, tidal y, proper radial distance.
        # Only E, not rho_b, sources enclosed gravitational mass.
        return [4*np.pi*r*r*e, -a, z,
                -(4/r + log_j_prime)*z - 4*log_j_prime*w/r,
                -(y*y + y*tidal_f + r*r*tidal_q)/r, 1/np.sqrt(f)]

    def surface_event(r, state):
        return state[1] - surface_h

    surface_event.terminal = True
    surface_event.direction = -1
    # Separate absolute tolerances reflect the different scales of the states;
    # dense output is needed for precise layer roots and arbitrary export grids.
    solution = solve_ivp(rhs, (r0, 60), initial, method="DOP853", rtol=rtol,
                         atol=[1e-12, 1e-14, 1e-11, 1e-12, 1e-10, 1e-11],
                         max_step=max_step_km, events=surface_event, dense_output=True)
    if not solution.success or len(solution.t_events[0]) != 1:
        raise RuntimeError(f"Surface not found: {solution.message}")
    return Star(eos, central_h, surface_h, solution, float(solution.t_events[0][0]))


def solve_target(eos, target_mass_msun=1.4, **options):
    """Shoot on h_c for the first target crossing in the supported NS bracket.

    A sampled downturn stops the search before the unstable continuation. This
    is not a full turning-point/stability analysis or a low-mass branch finder.
    """
    if target_mass_msun <= 0:
        raise ValueError("Target mass must be positive.")
    grid = np.linspace(0.06, min(0.58, 0.95*eos.hmax), 28)
    previous_h = grid[0]
    previous_mass = integrate(eos, previous_h, **options).mass_msun
    if previous_mass >= target_mass_msun:
        raise ValueError("Target is below the supported NS shooting bracket.")
    for h in grid[1:]:
        mass = integrate(eos, h, **options).mass_msun
        if mass >= target_mass_msun:
            root = brentq(lambda hc: integrate(eos, hc, **options).mass_msun - target_mass_msun,
                          previous_h, h, xtol=2e-12)
            return integrate(eos, root, **options)
        if mass < previous_mass:
            break
        previous_h, previous_mass = h, mass
    raise ValueError("Target not reached on the supported stable branch.")


def build_summary(star, cfg):
    """Computed surface observables and density-defined radial layers.

    Shell Delta m is gravitational mass, while shell labels use rest density.
    Proper thickness comes from the sixth ODE state, not coordinate subtraction.
    """
    cc = star.boundary(cfg["core_crust_rest_mass_density_g_cm3"])
    drip = star.boundary(cfg["neutron_drip_rest_mass_density_g_cm3"])
    pc, ec, rhoc = star.eos.at_h(star.central_h)
    r, m = star.radius_km, star.mass_msun
    bounds = [0, cc["radius_km"], drip["radius_km"], r]
    masses = [0, cc["enclosed_mass_msun"], drip["enclosed_mass_msun"], m]
    proper = [0] + [float(star.solution.sol(x)[5]) for x in bounds[1:]]
    layers = [{"name": name, "inner_radius_km": bounds[i], "outer_radius_km": bounds[i+1],
               "coordinate_thickness_km": bounds[i+1] - bounds[i],
               "proper_radial_thickness_km": proper[i+1] - proper[i],
               "gravitational_mass_contribution_msun": masses[i+1] - masses[i],
               "mass_fraction": (masses[i+1] - masses[i]) / m}
              for i, name in enumerate(LAYER_NAMES)]
    # I_geom is a length cubed: km^3 -> m^3, then c^2/G -> kg m^2.
    inertia = star.inertia_geom_km3 * 1e9 * C*C/G
    return {
        "eos": star.eos.name,
        "mass_msun": m, "radius_km": r, "compactness_GM_Rc2": star.surface[0]/r,
        "central_pseudo_enthalpy": star.central_h,
        "central_pressure_pa": pc/PRESSURE_TO_GEOM,
        "central_rest_mass_density_g_cm3": rhoc/1000,
        "central_energy_equivalent_density_g_cm3": ec/DENSITY_TO_GEOM/1000,
        "moment_of_inertia_kg_m2": inertia,
        "moment_of_inertia_g_cm2": inertia*1e7,
        "moment_of_inertia_geometric_msun3": star.inertia_geom_km3/RSUN_KM**3,
        "love_k2": star.love_k2, "tidal_deformability_lambda": star.tidal_lambda,
        "layer_thresholds_rest_mass_g_cm3": {
            "core_crust": cfg["core_crust_rest_mass_density_g_cm3"],
            "neutron_drip": cfg["neutron_drip_rest_mass_density_g_cm3"]},
        "core_crust_boundary": cc, "neutron_drip_boundary": drip, "layers": layers,
        "surface": {"pseudo_enthalpy_cutoff": star.surface_h,
                    "pressure_pa": star.eos.at_h(star.surface_h)[0]/PRESSURE_TO_GEOM},
    }
