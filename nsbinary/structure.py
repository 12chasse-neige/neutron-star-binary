"""Cold spherical neutron star: TOV, Hartle frame dragging and tidal Love number.

The integration uses G = c = 1 and lengths in km; public results use SI/cgs.
An equation of state supplies the background; density thresholds only label its
radial layers.

The integration state is
[mass, pseudo_enthalpy, frame_dragging, frame_dragging_derivative, tidal_y,
proper_radius], where frame_dragging = Omega - omega_fd (Hartle's bar-omega),
frame_dragging_derivative = d(bar-omega)/dr, tidal_y = r H'/H and proper_radius
is the proper radial distance.
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
    central_pseudo_enthalpy: float
    surface_pseudo_enthalpy: float
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
        """Match bar-omega(R), bar-omega'(R) to vacuum; return J/Omega."""
        _, _, frame_dragging, frame_dragging_derivative, _, _ = self.surface
        radius_km = self.radius_km
        return radius_km**4 * frame_dragging_derivative / (6 * (frame_dragging + radius_km * frame_dragging_derivative / 3))

    @property
    def love_k2(self):
        """Exterior l=2 matching for a smooth EOS without density jumps.

        ``compactness`` is not the SI speed of light ``C``. A finite surface
        density or internal discontinuity needs extra jump conditions.
        """
        compactness = self.surface[0] / self.radius_km
        tidal_y = self.surface[4]
        numerator = 8 / 5 * compactness**5 * (1 - 2*compactness)**2 * (2 + 2*compactness*(tidal_y - 1) - tidal_y)
        denominator = (2*compactness*(6 - 3*tidal_y + 3*compactness*(5*tidal_y - 8))
                       + 4*compactness**3*(13 - 11*tidal_y + compactness*(3*tidal_y - 2) + 2*compactness*compactness*(1 + tidal_y))
                       + 3*(1 - 2*compactness)**2*(2 - tidal_y + 2*compactness*(tidal_y - 1))*np.log1p(-2*compactness))
        return numerator / denominator

    @property
    def tidal_lambda(self):
        return 2 / 3 * self.love_k2 / (self.surface[0] / self.radius_km)**5

    def boundary(self, threshold_rest_mass_density_g_cm3):
        """Locate an assumed density threshold, not a computed phase change.

        Invert rho_b(h) first and then h(r), using the dense solution so the
        result does not depend on whether an export grid hits the boundary.
        """
        target_pseudo_enthalpy = self.eos.h_at_rho(threshold_rest_mass_density_g_cm3)
        if not self.surface_pseudo_enthalpy < target_pseudo_enthalpy < self.central_pseudo_enthalpy:
            raise ValueError("Layer threshold must lie inside the star.")
        radius_km = brentq(
            lambda radius_sample: (self.solution.sol(radius_sample)[1]
                                   - target_pseudo_enthalpy),
            self.solution.t[0], self.radius_km, xtol=1e-10)
        pressure, energy_density, rest_mass_density = self.eos.at_h(target_pseudo_enthalpy)
        state = self.solution.sol(radius_km)
        return {"radius_km": float(radius_km), "enclosed_mass_msun": float(state[0] / RSUN_KM),
                "pressure_pa": float(pressure / PRESSURE_TO_GEOM),
                "rest_mass_density_g_cm3": float(rest_mass_density / 1000),
                "energy_equivalent_density_g_cm3": float(energy_density / DENSITY_TO_GEOM / 1000)}

    def profile(self, profile_points=2401, extra_radii=()):
        """Export SI/cgs profiles with extra sampling in the thin crust."""
        # A second dense grid resolves the thin crust; include boundaries exactly.
        radius_km = np.unique(np.r_[
            0, np.linspace(self.solution.t[0], self.radius_km, profile_points),
            np.linspace(0.85*self.radius_km, self.radius_km, profile_points//3),
            extra_radii])
        state = self.solution.sol(np.maximum(radius_km, self.solution.t[0]))
        # Dense output begins at r0. Insert exact regular central values rather
        # than extrapolating the numerical solution into its singular endpoint.
        state[:, 0] = [0, self.central_pseudo_enthalpy, 1, 0, 2, 0]
        therm = np.array([self.eos.at_h(pseudo_enthalpy) for pseudo_enthalpy in state[1]])
        pressure, energy_density, rest_mass_density = therm.T
        dm_dr = 4 * np.pi * radius_km**2 * energy_density / RSUN_KM
        return {"radius_km": radius_km, "enclosed_mass_msun": state[0] / RSUN_KM,
                "pressure_pa": pressure / PRESSURE_TO_GEOM,
                "rest_mass_density_g_cm3": rest_mass_density / 1000,
                "energy_equivalent_density_g_cm3": energy_density / DENSITY_TO_GEOM / 1000,
                "mass_per_radius_msun_per_km": dm_dr,
                "proper_radius_km": state[5], "pseudo_enthalpy": state[1]}


def integrate(eos, central_pseudo_enthalpy, surface_pseudo_enthalpy=1e-10, rtol=2e-9,
              max_step_km=0.025, start_radius_km=1e-4):
    """Integrate mass, enthalpy, frame dragging, tidal y and proper radius.

    frame_dragging = Omega - omega is arbitrarily normalized at the centre.
    Matching to the exterior vacuum solution determines I = J/Omega.
    tidal_y = r H'/H is the l=2 static tidal metric perturbation; no density
    jump is imposed here.
    """
    if not 0 < surface_pseudo_enthalpy < central_pseudo_enthalpy < eos.hmax:
        raise ValueError(
            "Require 0 < surface_pseudo_enthalpy < central_pseudo_enthalpy "
            "< EOS maximum pseudo-enthalpy.")
    central_pressure, central_energy_density, _ = eos.at_h(central_pseudo_enthalpy)
    start_radius = start_radius_km
    # Leading regular-centre series avoids 1/r singularities in all six ODEs.
    # The Hartle amplitude is arbitrary; tidal_y starts at its l=2 value.
    initial_state = [
        4*np.pi/3 * central_energy_density*start_radius**3,
        central_pseudo_enthalpy - 2*np.pi*(
            central_energy_density/3 + central_pressure)*start_radius**2,
        1 + 8*np.pi/5*(
            central_energy_density + central_pressure)*start_radius**2,
        16*np.pi/5*(
            central_energy_density + central_pressure)*start_radius,
        2, start_radius]

    def rhs(radius_km, state):
        # gravity_acceleration is positive in enthalpy units, dh/dr = -a.
        # dlog_j_dr couples the TOV solution to Hartle rotation.
        gravitational_mass, pseudo_enthalpy, frame_dragging, frame_dragging_derivative, tidal_y, _ = state
        pressure, energy_density, _ = eos.at_h(max(pseudo_enthalpy, surface_pseudo_enthalpy * 0.01))
        lapse_squared = 1 - 2*gravitational_mass/radius_km
        if lapse_squared <= 0:
            raise RuntimeError("Trapped surface encountered.")
        gravity_acceleration = (gravitational_mass + 4*np.pi*radius_km**3*pressure) / (radius_km*radius_km*lapse_squared)
        dlog_j_dr = -4*np.pi*radius_km*(energy_density + pressure)/lapse_squared
        tidal_f = (1 - 4*np.pi*radius_km*radius_km*(energy_density - pressure))/lapse_squared
        tidal_q = (4*np.pi*(5*energy_density + 9*pressure + (energy_density+pressure)*eos.dedp(pressure))/lapse_squared
                   - 6/(radius_km*radius_km*lapse_squared) - 4*gravity_acceleration*gravity_acceleration)
        # Ordered derivatives: mass, enthalpy, frame dragging, its radial
        # derivative, tidal_y, proper radial distance.
        # Only E, not rho_b, sources enclosed gravitational mass.
        return [4*np.pi*radius_km*radius_km*energy_density, -gravity_acceleration, frame_dragging_derivative,
                -(4/radius_km + dlog_j_dr)*frame_dragging_derivative - 4*dlog_j_dr*frame_dragging/radius_km,
                -(tidal_y*tidal_y + tidal_y*tidal_f + radius_km*radius_km*tidal_q)/radius_km, 1/np.sqrt(lapse_squared)]

    def surface_event(radius_km, state):
        return state[1] - surface_pseudo_enthalpy

    surface_event.terminal = True
    surface_event.direction = -1
    # Separate absolute tolerances reflect the different scales of the states;
    # dense output is needed for precise layer roots and arbitrary export grids.
    solution = solve_ivp(rhs, (start_radius, 60), initial_state,
                         method="DOP853", rtol=rtol,
                         atol=[1e-12, 1e-14, 1e-11, 1e-12, 1e-10, 1e-11],
                         max_step=max_step_km, events=surface_event, dense_output=True)
    if not solution.success or len(solution.t_events[0]) != 1:
        raise RuntimeError(f"Surface not found: {solution.message}")
    return Star(eos, central_pseudo_enthalpy, surface_pseudo_enthalpy, solution, float(solution.t_events[0][0]))


def solve_target(eos, target_mass_msun=1.4, **options):
    """Shoot on h_c for the first target crossing in the supported NS bracket.

    A sampled downturn stops the search before the unstable continuation. This
    is not a full turning-point/stability analysis or a low-mass branch finder.
    """
    if target_mass_msun <= 0:
        raise ValueError("Target mass must be positive.")
    grid = np.linspace(0.06, min(0.58, 0.95*eos.hmax), 28)
    previous_pseudo_enthalpy = grid[0]
    previous_mass = integrate(eos, previous_pseudo_enthalpy, **options).mass_msun
    if previous_mass >= target_mass_msun:
        raise ValueError("Target is below the supported NS shooting bracket.")
    for pseudo_enthalpy in grid[1:]:
        mass = integrate(eos, pseudo_enthalpy, **options).mass_msun
        if mass >= target_mass_msun:
            root = brentq(lambda central_pseudo_enthalpy: integrate(eos, central_pseudo_enthalpy, **options).mass_msun - target_mass_msun,
                          previous_pseudo_enthalpy, pseudo_enthalpy, xtol=2e-12)
            return integrate(eos, root, **options)
        if mass < previous_mass:
            break
        previous_pseudo_enthalpy, previous_mass = pseudo_enthalpy, mass
    raise ValueError("Target not reached on the supported stable branch.")


def build_summary(star, cfg):
    """Computed surface observables and density-defined radial layers.

    Shell Delta m is gravitational mass, while shell labels use rest density.
    Proper thickness comes from the sixth ODE state, not coordinate subtraction.
    """
    core_crust_boundary = star.boundary(cfg["core_crust_rest_mass_density_g_cm3"])
    neutron_drip_boundary = star.boundary(cfg["neutron_drip_rest_mass_density_g_cm3"])
    central_pressure, central_energy_density, central_rest_mass_density = star.eos.at_h(star.central_pseudo_enthalpy)
    radius_km, mass_msun = star.radius_km, star.mass_msun
    layer_radii = [0, core_crust_boundary["radius_km"], neutron_drip_boundary["radius_km"], radius_km]
    layer_masses = [0, core_crust_boundary["enclosed_mass_msun"],
                    neutron_drip_boundary["enclosed_mass_msun"], mass_msun]
    proper_radii = [0] + [float(star.solution.sol(boundary_radius_km)[5])
                          for boundary_radius_km in layer_radii[1:]]
    layers = [{"name": name, "inner_radius_km": layer_radii[i], "outer_radius_km": layer_radii[i+1],
               "coordinate_thickness_km": layer_radii[i+1] - layer_radii[i],
               "proper_radial_thickness_km": proper_radii[i+1] - proper_radii[i],
               "gravitational_mass_contribution_msun": layer_masses[i+1] - layer_masses[i],
               "mass_fraction": (layer_masses[i+1] - layer_masses[i]) / mass_msun}
              for i, name in enumerate(LAYER_NAMES)]
    # I_geom is a length cubed: km^3 -> m^3, then c^2/G -> kg m^2.
    moment_of_inertia = star.inertia_geom_km3 * 1e9 * C*C/G
    return {
        "eos": star.eos.name,
        "mass_msun": mass_msun, "radius_km": radius_km,
        "compactness_GM_Rc2": star.surface[0]/radius_km,
        "central_pseudo_enthalpy": star.central_pseudo_enthalpy,
        "central_pressure_pa": central_pressure/PRESSURE_TO_GEOM,
        "central_rest_mass_density_g_cm3": central_rest_mass_density/1000,
        "central_energy_equivalent_density_g_cm3": central_energy_density/DENSITY_TO_GEOM/1000,
        "moment_of_inertia_kg_m2": moment_of_inertia,
        "moment_of_inertia_g_cm2": moment_of_inertia*1e7,
        "moment_of_inertia_geometric_msun3": star.inertia_geom_km3/RSUN_KM**3,
        "love_k2": star.love_k2, "tidal_deformability_lambda": star.tidal_lambda,
        "layer_thresholds_rest_mass_g_cm3": {
            "core_crust": cfg["core_crust_rest_mass_density_g_cm3"],
            "neutron_drip": cfg["neutron_drip_rest_mass_density_g_cm3"]},
        "core_crust_boundary": core_crust_boundary, "neutron_drip_boundary": neutron_drip_boundary, "layers": layers,
        "surface": {"pseudo_enthalpy_cutoff": star.surface_pseudo_enthalpy,
                    "pressure_pa": star.eos.at_h(star.surface_pseudo_enthalpy)[0]/PRESSURE_TO_GEOM},
    }
