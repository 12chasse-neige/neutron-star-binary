"""Two rigid effective rotators derived from Antonelli & Pizzochero (2017).

Two radial microphysics zones (core, inner crust) feed the straight-vortex line
integrals; the outer crust belongs to the normal rotator. This is the rigid
projection of their differential model, not a solution of its radial PDE.
All dynamics use SI, with energies and angular momenta divided by the stellar I.

Review build_reservoir first, then Binary.quantities (physical RHS), then
Binary.evolve (solver and exported traces).
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import G, c, h, neutron_mass
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

MSUN = 1.988409870698051e30
KAPPA = h / (2 * neutron_mass)


def gauss_interval(lo, hi, order):
    """Map Gauss-Legendre nodes/weights to a finite physical interval."""
    nodes, weights = leggauss(order)
    return lo + (nodes + 1) * (hi - lo) / 2, weights * (hi - lo) / 2


@dataclass
class Reservoir:
    """Canonical reservoir calibrated to the saved Hartle total inertia.

    total_I: kg m^2; I_v_fraction: I_v/I; B_eff: Y times the
    inertia-weighted free mobility. ``columns`` holds unscaled Newtonian line
    columns and normalized quadrature weights. Inertia fractions and B_eff are
    unchanged by the common rescaling.
    """
    total_I: float
    I_v_fraction: float
    B_eff: float
    metadata: dict
    columns: dict

    @property
    def I_p_fraction(self):
        return 1 - self.I_v_fraction


def build_reservoir(cfg, background):
    """Evaluate AP a, b, d, B_radial_free and their inertia-weighted cylindrical mean.

    Uses Newtonian integrals of baryon rest density on the saved TOV areal-radius
    profile, then multiplies all inertias by one factor to match the Hartle I.
    This hybrid convention is explicit; it is not a relativistic two-fluid EOS.
    """
    background = Path(background)
    star = json.loads((background / "summary.json").read_text())
    profile = np.genfromtxt(background / "profiles.csv", delimiter=",", names=True)
    # The structure export uses km and g/cm^3; every line and volume integral
    # below uses m and kg/m^3. Use rest density, not energy-equivalent density.
    radius = profile["radius_km"] * 1000
    rest_mass_density = PchipInterpolator(radius, profile["rest_mass_density_g_cm3"] * 1000)
    core_radius = star["core_crust_boundary"]["radius_km"] * 1000
    drip_radius = star["neutron_drip_boundary"]["radius_km"] * 1000
    stellar_radius = star["radius_km"] * 1000
    cylindrical_quadrature_order = cfg["cylindrical_quadrature_order"]
    line_quadrature_order = cfg["line_quadrature_order"]
    # Zone tuple: inner radius, outer radius, neutron mass fraction,
    # epsilon_n = 1 - m*_n/m_n, eta/(kappa*rho_n). No neutron reservoir beyond
    # the neutron-drip radius.
    zones = [(0., core_radius, cfg["core_neutron_fraction"], cfg["core_entrainment_epsilon_n"], cfg["core_drag_ratio"]),
             (core_radius, drip_radius, cfg["inner_crust_neutron_fraction"], cfg["inner_crust_entrainment_epsilon_n"], cfg["inner_crust_drag_ratio"])]
    for _, _, neutron_fraction, entrainment_epsilon, drag_ratio in zones:
        if not (0 <= neutron_fraction < 1 and entrainment_epsilon < 1-neutron_fraction and drag_ratio >= 0):
            raise ValueError("Require 0 <= x_n < 1, epsilon_n < x_p, drag >= 0.")
    if not 0 <= cfg["unpinned_fraction"] <= 1:
        raise ValueError("Unpinned fraction must lie in [0,1].")
    # Integrate all material for I_N; the outer crust contributes to the normal
    # component even though it contributes nothing to the superfluid line columns.
    newtonian_total_inertia = 0.
    zone_inertias = []
    for lo, hi in [(0, core_radius), (core_radius, drip_radius), (drip_radius, stellar_radius)]:
        radius_nodes, radius_weights = gauss_interval(
            lo, hi, 2*cylindrical_quadrature_order)
        zone_inertias.append(8*np.pi/3 * np.sum(
            radius_weights * radius_nodes**4 * rest_mass_density(radius_nodes)))
        newtonian_total_inertia += zone_inertias[-1]
    I_v_radial = sum(zone_inertias[i]*z[2]/(1-z[3]) for i, z in enumerate(zones))
    # Splitting at the core radius resolves the tangent-to-core boundary. Avoid
    # endpoint 0/0.
    x_zones, x_zone_weights = zip(*(gauss_interval(lo, hi, cylindrical_quadrature_order) for lo, hi in [(0, core_radius), (core_radius, drip_radius)]))
    x, x_weights = np.concatenate(x_zones), np.concatenate(x_zone_weights)
    a, b, d = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
    for j, x_line in enumerate(x):
        for lo, hi, neutron_fraction, entrainment_epsilon, drag_ratio in zones:
            if x_line >= hi:
                continue
            # Positive-z segment inside this spherical shell. Double it for
            # reflection symmetry; splitting shells prevents averaging a jump
            # in entrainment or drag across an unresolved quadrature interval.
            z_lower = np.sqrt(max(lo*lo-x_line*x_line, 0))
            z_upper = np.sqrt(hi*hi-x_line*x_line)
            z_nodes, z_weights = gauss_interval(
                z_lower, z_upper, line_quadrature_order)
            line_column = 2*neutron_fraction*np.sum(
                z_weights*rest_mass_density(np.hypot(x_line, z_nodes)))
            a[j] += line_column
            b[j] += line_column/(1-entrainment_epsilon)
            d[j] += line_column*drag_ratio
    # AP free-vortex response. All three columns have units kg/m^2, so the two
    # B columns are dimensionless. Large d suppresses radial motion rather than
    # forcing arbitrarily fast material corotation.
    B_radial_free = d*b/(a*a+d*d)
    B_azimuthal_free = a*b/(a*a+d*d)
    inertia_weights = 2*np.pi * x_weights*x**3*b
    I_v_cylindrical = np.sum(inertia_weights)
    # I_v_fraction comes from the independently evaluated spherical integral;
    # the cylindrical weights are normalized separately so B_free is not biased
    # by quadrature error.
    I_v_fraction = I_v_radial/newtonian_total_inertia
    B_free = float(np.sum(inertia_weights*B_radial_free)/I_v_cylindrical)
    B_eff = cfg["unpinned_fraction"] * B_free
    if not 0 < I_v_fraction < 1:
        raise ValueError("The effective inertia matrix must be positive.")
    metadata = {
        "stellar_mass_msun": star["mass_msun"], "radius_m": stellar_radius,
        "core_radius_m": core_radius, "neutron_drip_radius_m": drip_radius,
        "total_I_kg_m2": star["moment_of_inertia_kg_m2"],
        "newtonian_rest_density_I_kg_m2": newtonian_total_inertia,
        "common_Hartle_rescaling": star["moment_of_inertia_kg_m2"]/newtonian_total_inertia,
        "effective_v_inertia_fraction": I_v_fraction, "effective_p_inertia_fraction": 1-I_v_fraction,
        "core_material_inertia_fraction": zone_inertias[0]/newtonian_total_inertia,
        "inner_crust_material_inertia_fraction": zone_inertias[1]/newtonian_total_inertia,
        "outer_crust_material_inertia_fraction": zone_inertias[2]/newtonian_total_inertia,
        "mean_B_free": B_free, "mean_B_effective": B_eff,
        "kappa_m2_per_s": KAPPA,
    }
    columns = {"cylindrical_radius_m": x, "a_kg_m2": a, "b_kg_m2": b,
               "d_kg_m2": d, "B_radial_free": B_radial_free, "B_azimuthal_free": B_azimuthal_free,
               "inertia_weight": inertia_weights/I_v_cylindrical}
    return Reservoir(star["moment_of_inertia_kg_m2"], I_v_fraction, B_eff, metadata, columns)


class Binary:
    """Circular Hamiltonian plus quadrupole torque and leading GW radiation.

    omega_kepler denotes the Kepler frequency inferred from orbital angular
    momentum L.
    Keeping U(L, gamma) in the Hamiltonian adds its small frequency correction;
    the full system then conserves energy and angular momentum when GW and
    mutual friction are disabled. This avoids treating heat as an extra tidal
    torque, which would double-count the spin-orbit energy transfer.
    """

    def __init__(self, cfg, reservoir, rigid=False):
        """Set fixed masses/quadrupole; ``rigid`` is the single-rotator control."""
        self.cfg, self.reservoir = cfg, reservoir
        self.I = reservoir.total_I
        self.I_v_fraction = 0. if rigid else reservoir.I_v_fraction
        self.I_p_fraction = 1-self.I_v_fraction
        self.B_eff = 0. if rigid else reservoir.B_eff
        neutron_star_mass = reservoir.metadata["stellar_mass_msun"]*MSUN
        mass_ratio = cfg["mass_ratio_companion_over_ns"]
        if mass_ratio <= 0 or cfg["ellipticity"] < 0:
            raise ValueError("Require positive mass ratio and nonnegative ellipticity.")
        self.total_mass, self.reduced_mass = neutron_star_mass*(1+mass_ratio), neutron_star_mass*mass_ratio/(1+mass_ratio)
        self.companion_mass_fraction = mass_ratio/(1+mass_ratio)
        self.ellipticity = cfg["ellipticity"]
        # Divide budgets by I to avoid enormous SI state magnitudes. The
        # equations remain SI; torque/I is an angular acceleration and flux/I
        # is specific power. Restore I only when exporting physical budgets.
        self.L_over_I_coefficient = self.reduced_mass*(G*self.total_mass)**(2/3)/self.I
        self.gw_flux_over_I_coefficient = 32/5 * G**(7/3)/c**5*self.reduced_mass**2*self.total_mass**(4/3)/self.I

    def quantities(self, omega_kepler, gamma, omega_p, omega_v):
        """Evaluate physical rates and budgets at scalar or array states.

        Return order: L_over_I, U_over_I, N_over_I, omega_orb, gw_flux_over_I,
        omega_kepler_dot, omega_p_dot, omega_v_dot, P_heat_over_I. Here
        omega_kepler is the Kepler angular frequency, not Hz, and gamma =
        2(theta - f). Positive N_over_I spins up the charged component.
        """
        L_over_I = self.L_over_I_coefficient*omega_kepler**(-1/3)  # L/I
        U_over_I = -.75*self.companion_mass_fraction*self.ellipticity*omega_kepler*omega_kepler*np.cos(gamma)  # U/I
        N_over_I = -1.5*self.companion_mass_fraction*self.ellipticity*omega_kepler*omega_kepler*np.sin(gamma)  # N/I
        omega_orb = omega_kepler-6*U_over_I/L_over_I  # d(E_Kep+U)/dL
        gw_flux_over_I = self.gw_flux_over_I_coefficient*omega_kepler**(10/3)
        # Circular reaction: L proportional to n^(-1/3). Removing orbital
        # angular momentum raises n. F_J = F_GW/Omega_orb preserves dH/dt = -F_GW
        # for this Hamiltonian, including its small conservative correction.
        omega_kepler_dot = 3*omega_kepler/L_over_I*(N_over_I+gw_flux_over_I/omega_orb)
        drag_rate = 2*self.B_eff*omega_v
        omega_v_dot = -drag_rate*(omega_v-omega_p)
        # Friction exchanges spin momentum internally:
        # I_p_fraction*dOmega_p/dt + I_v_fraction*dOmega_v/dt = N/I. Heating
        # follows from the same exchange and is NOT another orbital torque.
        omega_p_dot = (N_over_I-self.I_v_fraction*omega_v_dot)/self.I_p_fraction
        if self.I_v_fraction == 0:
            omega_v_dot = omega_p_dot
        P_heat_over_I = self.I_v_fraction*drag_rate*(omega_v-omega_p)**2 if self.I_v_fraction else np.zeros_like(omega_p)
        return L_over_I, U_over_I, N_over_I, omega_orb, gw_flux_over_I, omega_kepler_dot, omega_p_dot, omega_v_dot, P_heat_over_I

    def evolve(self, initialization="locked", sample_count=5000):
        """Integrate in log(n) and return the exported trace columns.

        ``paper`` starts from the original frequencies and rounded phase;
        ``locked`` starts near a stable rigid equilibrium for conditional
        locking comparisons. ``sample_count`` adds dense log-frequency samples;
        adaptive nodes are kept as well, so it does not set ODE resolution.
        """
        cfg = self.cfg
        if initialization == "paper":
            omega_kepler_initial = 2*np.pi*cfg["initial_orbital_frequency_hz"]
            omega_p_initial = omega_v_initial = 2*np.pi*cfg["initial_spin_frequency_hz"]
            gamma_initial = -2*cfg["initial_f_minus_theta_rad"]
        elif initialization == "locked":
            # Use the rigid point-particle chirp to select a stable phase, then
            # match both spins to the corrected Hamiltonian orbital frequency.
            # No constraint subsequently forces the components to remain locked.
            omega_kepler_initial = 2*np.pi*cfg["locked_start_orbital_frequency_hz"]
            omega_kepler_dot_initial = 3*self.gw_flux_over_I_coefficient*omega_kepler_initial**(11/3)/self.L_over_I_coefficient
            sin_gamma = -2*omega_kepler_dot_initial/(3*self.companion_mass_fraction*self.ellipticity*omega_kepler_initial*omega_kepler_initial) if self.ellipticity else 0.
            if abs(sin_gamma) >= 1:
                raise ValueError("Locked initialization is above the rigid breaking frequency.")
            gamma_initial = np.arcsin(sin_gamma)
            omega_p_initial = omega_v_initial = self.quantities(omega_kepler_initial, gamma_initial, omega_kepler_initial, omega_kepler_initial)[3]
        else:
            raise ValueError("Initialization must be 'paper' or 'locked'.")
        omega_kepler_final = 2*np.pi*cfg["final_orbital_frequency_hz"]
        if not 0 < omega_kepler_initial < omega_kepler_final or omega_p_initial <= 0:
            raise ValueError("Require positive initial spin and increasing orbital frequency.")
        # Variables: gamma, Omega_p, Omega_v, t, Q/I, E_GW/I, J_GW/I,
        # W_star/I, W_orbit/I, GW phase difference from the point-particle
        # chirp.
        initial_state = np.array([gamma_initial, omega_p_initial, omega_v_initial, 0, 0, 0, 0, 0, 0, 0.])

        def rhs(log_omega_kepler, state):
            omega_kepler = np.exp(log_omega_kepler)
            gamma, omega_p, omega_v = state[:3]
            (L_over_I, U_over_I, N_over_I, omega_orb, gw_flux_over_I,
             omega_kepler_dot, omega_p_dot, omega_v_dot,
             P_heat_over_I) = self.quantities(omega_kepler, gamma, omega_p, omega_v)
            if np.real(omega_kepler_dot) <= 0 or np.real(omega_v) <= 0:
                raise RuntimeError("Frequency ceased to be monotone or vortex rotation reversed.")
            # The independent variable is log(omega_kepler), so this factor is
            # dt/dlog(omega_kepler), not an integration step. The final term
            # accumulates 2 integral(omega_orb/omega_kepler_dot - omega_kepler/
            # omega_kepler_dot_point_particle) d(omega_kepler), comparing
            # time-domain GW phases at equal Kepler frequency and initial phase.
            dt_dlog_omega_kepler = omega_kepler/omega_kepler_dot
            omega_kepler_dot_point_particle = 3*gw_flux_over_I/L_over_I
            return np.array([2*(omega_p-omega_orb)*dt_dlog_omega_kepler,
                             omega_p_dot*dt_dlog_omega_kepler,
                             omega_v_dot*dt_dlog_omega_kepler,
                             dt_dlog_omega_kepler,
                             P_heat_over_I*dt_dlog_omega_kepler,
                             gw_flux_over_I*dt_dlog_omega_kepler,
                             gw_flux_over_I/omega_orb*dt_dlog_omega_kepler,
                             omega_p*N_over_I*dt_dlog_omega_kepler,
                             omega_orb*N_over_I*dt_dlog_omega_kepler,
                             2*omega_kepler*(omega_orb/omega_kepler_dot
                                              - omega_kepler/omega_kepler_dot_point_particle)])

        atol = np.array([cfg["phase_absolute_tolerance_rad"], 1e-10, 1e-10,
                         1e-5, 1e-9, 1e-5, 1e-7, 1e-9, 1e-9, 1e-7])

        def jacobian(log_omega_kepler, state):
            # Complex-step derivatives also handle identically zero diagnostic
            # variables without finite-difference step escalation.
            perturbation = 1e-24j*np.eye(len(state))
            return np.column_stack([np.imag(rhs(log_omega_kepler, state+delta))/1e-24
                                    for delta in perturbation])

        # Before capture the circulation is fast compared with friction; an
        # explicit high-order solver is much faster. Locked drag sweeps include
        # stiff relaxation and use an implicit solver.
        method = "DOP853" if initialization == "paper" else "Radau"
        options = {"jac": jacobian} if method == "Radau" else {}
        solution = solve_ivp(rhs, (np.log(omega_kepler_initial), np.log(omega_kepler_final)), initial_state, method=method,
                             rtol=cfg["relative_tolerance"], atol=atol,
                             dense_output=True, max_step=.015, **options)
        if not solution.success:
            raise RuntimeError(solution.message)
        # Combine the adaptive mesh with dense log-frequency sampling for export.
        log_omega_kepler = np.unique(np.r_[np.linspace(np.log(omega_kepler_initial), np.log(omega_kepler_final), sample_count), solution.t])
        omega_kepler = np.exp(log_omega_kepler)
        state = solution.sol(log_omega_kepler)
        (gamma, omega_p, omega_v, time, heat_over_I, gw_energy_over_I,
         gw_angular_momentum_over_I, star_torque_work_over_I,
         orbit_torque_work_over_I,
         gw_phase_difference) = state
        (L_over_I, U_over_I, N_over_I, omega_orb, gw_flux_over_I,
         omega_kepler_dot, omega_p_dot, omega_v_dot,
         P_heat_over_I) = self.quantities(omega_kepler, gamma, omega_p, omega_v)
        spin_energy_over_I = .5*(self.I_p_fraction*omega_p*omega_p+self.I_v_fraction*omega_v*omega_v)
        orbital_energy_over_I = -.5*L_over_I*omega_kepler
        lag_relaxation_time = np.full_like(omega_v, np.inf) if self.B_eff == 0 else self.I_p_fraction/(2*self.B_eff*omega_v)
        # Inertia-weighted rms of x*Bx, then multiply by signed lag and Y.
        # This is rms over radius of the mobility-averaged drift, not an rms
        # over individual pinned/mobile vortices.
        vortex_radial_velocity_rms = np.sqrt(np.sum(reservoir_weights(self.reservoir)*
                            (self.reservoir.columns["cylindrical_radius_m"]*
                             self.reservoir.columns["B_radial_free"])**2))
        vortex_radial_velocity_rms *= cfg["unpinned_fraction"]*(omega_v-omega_p)
        out = {"time_s": time, "kepler_orbital_frequency_hz": omega_kepler/(2*np.pi),
               "orbital_frequency_hz": omega_orb/(2*np.pi), "crust_frequency_hz": omega_p/(2*np.pi),
               "canonical_superfluid_frequency_hz": omega_v/(2*np.pi), "gamma_rad": gamma,
               "canonical_lag_rad_s": omega_v-omega_p, "heat_J": heat_over_I*self.I,
               "mutual_friction_power_W": P_heat_over_I*self.I, "gw_energy_J": gw_energy_over_I*self.I,
               "star_torque_work_J": star_torque_work_over_I*self.I, "orbit_torque_work_J": orbit_torque_work_over_I*self.I,
               "spin_energy_change_J": (spin_energy_over_I-spin_energy_over_I[0])*self.I,
               "orbital_binding_energy_loss_J": (orbital_energy_over_I[0]-orbital_energy_over_I)*self.I,
               "interaction_energy_change_J": (U_over_I-U_over_I[0])*self.I,
               "instantaneous_heat_over_gw_power": P_heat_over_I/gw_flux_over_I,
               "lag_relaxation_time_s": lag_relaxation_time,
               "signed_rms_vortex_radial_velocity_m_s": vortex_radial_velocity_rms,
               "core_neutron_transport_frequency_hz": (omega_v-cfg["core_entrainment_epsilon_n"]*omega_p)/(1-cfg["core_entrainment_epsilon_n"])/(2*np.pi),
               "crust_neutron_transport_frequency_hz": (omega_v-cfg["inner_crust_entrainment_epsilon_n"]*omega_p)/(1-cfg["inner_crust_entrainment_epsilon_n"])/(2*np.pi),
               "gw_phase_difference_from_point_particle_rad": gw_phase_difference,
               "vortex_areal_density_m2": 2*omega_v/KAPPA}
        if self.I_v_fraction == 0:
            out["signed_rms_vortex_radial_velocity_m_s"] *= 0
        return out


def summarize_trace(trace):
    """Final scalar values of one evolved scenario in physical SI units."""
    heat = float(trace["heat_J"][-1])
    orbital_loss = float(trace["orbital_binding_energy_loss_J"][-1])
    return {
        "duration_s": float(trace["time_s"][-1]),
        "final_orbital_frequency_hz": float(trace["orbital_frequency_hz"][-1]),
        "final_crust_frequency_hz": float(trace["crust_frequency_hz"][-1]),
        "final_superfluid_frequency_hz": float(trace["canonical_superfluid_frequency_hz"][-1]),
        "heat_J": heat,
        "gw_energy_J": float(trace["gw_energy_J"][-1]),
        "star_torque_work_J": float(trace["star_torque_work_J"][-1]),
        "orbital_binding_energy_loss_J": orbital_loss,
        "heat_over_orbital_binding_energy_loss": heat/orbital_loss,
        "max_absolute_canonical_lag_rad_s": float(np.max(np.abs(trace["canonical_lag_rad_s"]))),
        "max_absolute_vortex_radial_velocity_m_s": float(np.max(np.abs(trace["signed_rms_vortex_radial_velocity_m_s"]))),
        "final_gw_phase_difference_from_point_particle_rad": float(trace["gw_phase_difference_from_point_particle_rad"][-1]),
    }


def reservoir_weights(reservoir):
    """Normalized canonical-inertia weights for cylindrical diagnostics."""
    return reservoir.columns["inertia_weight"]
