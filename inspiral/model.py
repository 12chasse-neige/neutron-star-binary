"""Two rigid effective rotators derived from Antonelli & Pizzochero (2017).

Two radial microphysics zones (core, inner crust) feed the straight-vortex line
integrals. The outer crust belongs to the normal rotator. This is the rigid
projection of their differential model, not a solution of its radial PDE.
All dynamics use SI, with energies and angular momenta divided by stellar I.
See MODEL.md for derivations and state tables. Review build_reservoir first,
then Binary.quantities (physical RHS), then Binary.evolve (solver/diagnostics).
"""

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.constants import G, c, h, neutron_mass

MSUN = 1.988409870698051e30
KAPPA = h / (2 * neutron_mass)


def gauss_interval(lo, hi, order):
    """Map Gauss-Legendre nodes/weights to a finite physical interval."""
    nodes, weights = leggauss(order)
    return lo + (nodes + 1) * (hi - lo) / 2, weights * (hi - lo) / 2


@dataclass
class Reservoir:
    """Canonical reservoir calibrated to the saved Hartle total inertia.

    total_I: kg m^2; fv: I_v/I; B: Y times the inertia-weighted free mobility.
    columns holds *unscaled* Newtonian line columns and normalized quadrature
    weights. Fractions and B are unchanged by the common Hartle rescaling.
    """
    total_I: float
    fv: float
    B: float
    metadata: dict
    columns: dict

    @property
    def fp(self):
        return 1 - self.fv


def build_reservoir(cfg, background, order_multiplier=1):
    """Evaluate AP a,b,d,Bx and their inertia-weighted cylindrical mean.

    Use Newtonian integrals of baryon rest density on the saved TOV areal-radius
    profile, then multiply all inertias by one factor to match the Hartle I.
    This hybrid convention is explicit; it is not a relativistic two-fluid EOS.
    """
    background = Path(background)
    summary_path = background / "summary.json"
    profile_path = background / "profiles.csv"
    star = json.loads(summary_path.read_text())
    profile = np.genfromtxt(profile_path, delimiter=",", names=True)
    # APR exports km and g/cm^3; every line/volume integral below uses m and
    # kg/m^3. Use rest density, not energy-equivalent density, in this reduction.
    radius = profile["radius_km"] * 1000
    rho = PchipInterpolator(radius, profile["rest_mass_density_g_cm3"] * 1000)
    rc = star["core_crust_boundary"]["radius_km"] * 1000
    rd = star["neutron_drip_boundary"]["radius_km"] * 1000
    outer = star["radius_km"] * 1000
    nq = cfg["cylindrical_quadrature_order"] * order_multiplier
    nz = cfg["line_quadrature_order"] * order_multiplier
    # Zone tuple: inner radius, outer radius, neutron mass fraction,
    # epsilon_n=1-m*_n/m_n, eta/(kappa*rho_n). No neutron reservoir beyond rd.
    zones = [(0., rc, cfg["core_neutron_fraction"], cfg["core_entrainment_epsilon_n"], cfg["core_drag_ratio"]),
             (rc, rd, cfg["inner_crust_neutron_fraction"], cfg["inner_crust_entrainment_epsilon_n"], cfg["inner_crust_drag_ratio"])]
    for _, _, xn, en, drag in zones:
        if not (0 <= xn < 1 and en < 1-xn and drag >= 0):
            raise ValueError("Require 0 <= x_n < 1, epsilon_n < x_p, drag >= 0.")
    if not 0 <= cfg["unpinned_fraction"] <= 1:
        raise ValueError("Unpinned fraction must lie in [0,1].")
    # Integrate all material for I_N; the outer crust contributes to the normal
    # component even though it contributes nothing to the superfluid line columns.
    total_newton = 0.
    zone_inertias = []
    for lo, hi in [(0, rc), (rc, rd), (rd, outer)]:
        rr, ww = gauss_interval(lo, hi, 2*nq)
        zone_inertias.append(8*np.pi/3 * np.sum(ww * rr**4 * rho(rr)))
        total_newton += zone_inertias[-1]
    iv_radial = sum(zone_inertias[i]*z[2]/(1-z[3]) for i, z in enumerate(zones))
    # Splitting at rc resolves the tangent-to-core boundary. Avoid endpoint 0/0.
    xs, xw = zip(*(gauss_interval(lo, hi, nq) for lo, hi in [(0, rc), (rc, rd)]))
    x, wx = np.concatenate(xs), np.concatenate(xw)
    a, b, d = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
    for j, xx in enumerate(x):
        for lo, hi, xn, en, drag in zones:
            if xx >= hi:
                continue
            # Positive-z segment inside this spherical shell. Double it for
            # reflection symmetry; splitting shells prevents averaging a jump
            # in entrainment or drag across an unresolved quadrature interval.
            zlo = np.sqrt(max(lo*lo-xx*xx, 0))
            zhi = np.sqrt(hi*hi-xx*xx)
            zz, wz = gauss_interval(zlo, zhi, nz)
            column = 2*xn*np.sum(wz*rho(np.hypot(xx, zz)))
            a[j] += column
            b[j] += column/(1-en)
            d[j] += column*drag
    # AP free-vortex response. All three columns have units kg/m^2, so bx and
    # btheta are dimensionless. Large d suppresses radial motion, rather than
    # forcing arbitrarily fast material corotation.
    bx = d*b/(a*a+d*d)
    btheta = a*b/(a*a+d*d)
    iw = 2*np.pi * wx*x**3*b
    iv_cyl = np.sum(iw)
    # fv from the independently evaluated spherical integral; the cylindrical
    # weights are normalized separately so B is not biased by quadrature error.
    fv = iv_radial/total_newton
    Bfree = float(np.sum(iw*bx)/iv_cyl)
    B = cfg["unpinned_fraction"] * Bfree
    if not 0 < fv < 1:
        raise ValueError("The effective inertia matrix must be positive.")
    metadata = {
        "stellar_mass_msun": star["mass_msun"], "radius_m": outer,
        "core_radius_m": rc, "neutron_drip_radius_m": rd,
        "total_I_kg_m2": star["moment_of_inertia_kg_m2"],
        "newtonian_rest_density_I_kg_m2": total_newton,
        "common_Hartle_rescaling": star["moment_of_inertia_kg_m2"]/total_newton,
        "effective_v_inertia_fraction": fv, "effective_p_inertia_fraction": 1-fv,
        "core_material_inertia_fraction": zone_inertias[0]/total_newton,
        "inner_crust_material_inertia_fraction": zone_inertias[1]/total_newton,
        "outer_crust_material_inertia_fraction": zone_inertias[2]/total_newton,
        "mean_B_free": Bfree, "mean_B_effective": B,
        "cylindrical_vs_spherical_Iv_relative_error": abs(iv_cyl/iv_radial-1),
        "profile_sha256": hashlib.sha256(profile_path.read_bytes()).hexdigest(),
        "background_summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        "kappa_m2_per_s": KAPPA,
    }
    columns = {"cylindrical_radius_m": x, "a_kg_m2": a, "b_kg_m2": b,
               "d_kg_m2": d, "B_radial_free": bx, "B_azimuthal_free": btheta,
               "inertia_weight": iw/iv_cyl}
    return Reservoir(star["moment_of_inertia_kg_m2"], fv, B, metadata, columns)


class Binary:
    """Circular Hamiltonian plus quadrupole torque and leading GW radiation.

    n denotes the Kepler frequency inferred from orbital angular momentum L.
    Keeping U(L,gamma) in the Hamiltonian adds its small frequency correction;
    the full system then conserves energy and angular momentum when GW and
    mutual friction are disabled. This avoids treating heat as an extra tidal
    torque, which would double-count the spin-orbit energy transfer.
    """

    def __init__(self, cfg, reservoir, rigid=False, torque_enabled=True, gw_enabled=True):
        """Set fixed masses/quadrupole and optional limiting-control switches.

        A rigid control uses all I in p and B=0; its ov is a copy of the spin.
        Disabling torque or radiation supports independent conservation tests.
        """
        self.cfg, self.reservoir = cfg, reservoir
        self.I = reservoir.total_I
        self.fv = 0. if rigid else reservoir.fv
        self.fp = 1-self.fv
        self.B = 0. if rigid else reservoir.B
        m = reservoir.metadata["stellar_mass_msun"]*MSUN
        q = cfg["mass_ratio_companion_over_ns"]
        if q <= 0 or cfg["ellipticity"] < 0:
            raise ValueError("Require positive mass ratio and nonnegative ellipticity.")
        self.mt, self.mu = m*(1+q), m*q/(1+q)
        self.ratio = q/(1+q)
        self.eps = cfg["ellipticity"] if torque_enabled else 0.
        self.gw_enabled = gw_enabled
        # Divide budgets by I to avoid enormous SI state magnitudes. The
        # equations remain SI; torque/I is an angular acceleration and flux/I
        # is specific power. Restore I only when exporting physical budgets.
        self.ell_coeff = self.mu*(G*self.mt)**(2/3)/self.I
        self.flux_coeff = 32/5 * G**(7/3)/c**5*self.mu**2*self.mt**(4/3)/self.I

    def quantities(self, n, phase, op, ov):
        """Evaluate physical rates and budgets at scalar or array states.

        Return order: L/I, U/I, N/I, Omega_orb, F_GW/I, dn/dt,
        dOmega_p/dt, dOmega_v/dt, P_heat/I. Here n is a Kepler angular
        frequency, not Hz, and phase=2(theta-f). Positive N spins up p.
        """
        ell = self.ell_coeff*n**(-1/3)  # L/I
        u = -.75*self.ratio*self.eps*n*n*np.cos(phase)  # U/I
        torque = -1.5*self.ratio*self.eps*n*n*np.sin(phase)  # N/I
        orb = n-6*u/ell  # d(E_Kep+U)/dL
        flux = self.flux_coeff*n**(10/3) if self.gw_enabled else 0.
        # Circular reaction: L proportional to n^(-1/3). Removing orbital
        # angular momentum raises n. F_J=F_GW/Omega_orb preserves dH/dt=-F_GW
        # for this Hamiltonian, including its small conservative correction.
        ndot = 3*n/ell*(torque+flux/orb)
        k = 2*self.B*ov
        vdot = -k*(ov-op)
        # Friction exchanges spin momentum internally: fp*pdot+fv*vdot=N/I.
        # Heating follows from the same exchange and is NOT another orbit torque.
        pdot = (torque-self.fv*vdot)/self.fp
        if self.fv == 0:
            vdot = pdot
        heat = self.fv*k*(ov-op)**2 if self.fv else np.zeros_like(op)
        return ell, u, torque, orb, flux, ndot, pdot, vdot, heat

    def rhs_time(self, t, y):
        """Time-domain interface for conservation/limiting-solution tests.

        y=[n,gamma,op,ov,Q/I,E_GW/I,J_GW/I,W_star/I,W_orbit/I].
        This differs from evolve's log-frequency state: n is a state here,
        elapsed time is not, and no phase-difference diagnostic is integrated.
        """
        n, phase, op, ov = y[:4]
        ell, u, torque, orb, flux, ndot, pdot, vdot, heat = self.quantities(n, phase, op, ov)
        return np.array([ndot, 2*(op-orb), pdot, vdot, heat, flux, flux/orb,
                         op*torque, orb*torque])

    def evolve(self, initialization="locked", accuracy=1., sample_count=5000):
        """Integrate in log(n), then export traces and physical audit summaries.

        ``paper`` tests capture from the original frequencies and rounded phase;
        ``locked`` starts near a stable rigid equilibrium for conditional
        locking comparisons. accuracy>1 tightens tolerances and maximum steps.
        sample_count adds dense log-frequency samples; adaptive nodes are kept
        as well, so neither this value nor CSV decimation sets ODE resolution.
        """
        cfg = self.cfg
        if initialization == "paper":
            n0 = 2*np.pi*cfg["initial_orbital_frequency_hz"]
            op0 = ov0 = 2*np.pi*cfg["initial_spin_frequency_hz"]
            phase0 = -2*cfg["initial_f_minus_theta_rad"]
        elif initialization == "locked":
            # Use the rigid point-particle chirp to select a stable phase, then
            # match both spins to the corrected Hamiltonian orbital frequency.
            # No constraint subsequently forces the components to remain locked.
            n0 = 2*np.pi*cfg["locked_start_orbital_frequency_hz"]
            ndot0 = 3*self.flux_coeff*n0**(11/3)/self.ell_coeff
            sine = -2*ndot0/(3*self.ratio*self.eps*n0*n0) if self.eps else 0.
            if abs(sine) >= 1:
                raise ValueError("Locked initialization is above the rigid breaking frequency.")
            phase0 = np.arcsin(sine)
            op0 = ov0 = self.quantities(n0, phase0, n0, n0)[3]
        else:
            raise ValueError("Initialization must be 'paper' or 'locked'.")
        n1 = 2*np.pi*cfg["final_orbital_frequency_hz"]
        if not 0 < n0 < n1 or op0 <= 0:
            raise ValueError("Require positive initial spin and increasing orbital frequency.")
        # Variables: phase, Omega_p, Omega_v, t, Q/I, E_GW/I, J_GW/I,
        # W_star/I, W_orbit/I, phase_difference_from_point_particle.
        initial = np.array([phase0, op0, ov0, 0, 0, 0, 0, 0, 0, 0.])
        def rhs(logn, y):
            n = np.exp(logn)
            phase, op, ov = y[:3]
            ell, u, torque, orb, flux, ndot, pdot, vdot, heat = self.quantities(n, phase, op, ov)
            if np.real(ndot) <= 0 or np.real(ov) <= 0:
                raise RuntimeError("Frequency ceased to be monotone or vortex rotation reversed.")
            dt = n/ndot
            ndot_pp = 3*flux/ell
            # dt here is dt/d(log n), not an integration step. The final term
            # accumulates 2*integral(Omega_orb/ndot - n/ndot_pp) dn, comparing
            # time-domain GW phases at equal Kepler frequency and initial phase.
            return np.array([2*(op-orb)*dt, pdot*dt, vdot*dt, dt,
                             heat*dt, flux*dt, flux/orb*dt, op*torque*dt,
                             orb*torque*dt, 2*n*(orb/ndot-n/ndot_pp)])
        atol = np.array([cfg["phase_absolute_tolerance_rad"], 1e-10, 1e-10,
                         1e-5, 1e-9, 1e-5, 1e-7, 1e-9, 1e-9, 1e-7])/accuracy
        def jacobian(logn, y):
            # Complex-step derivatives also handle identically zero diagnostic
            # variables without finite-difference step escalation.
            perturbation = 1e-24j*np.eye(len(y))
            return np.column_stack([np.imag(rhs(logn, y+delta))/1e-24 for delta in perturbation])
        # Before capture the circulation is fast compared with friction; an
        # explicit high-order solver is much faster. Locked drag sweeps include
        # stiff relaxation and use an implicit solver.
        method = "DOP853" if initialization == "paper" else "Radau"
        options = {"jac": jacobian} if method == "Radau" else {}
        solution = solve_ivp(rhs, (np.log(n0), np.log(n1)), initial, method=method,
                             rtol=cfg["relative_tolerance"]/accuracy, atol=atol,
                             dense_output=True, max_step=.015/np.sqrt(accuracy), **options)
        if not solution.success:
            raise RuntimeError(solution.message)
        # Audit before CSV thinning, combining the adaptive mesh with dense
        # sampling so rapid phase motion is not judged from coarse exports.
        logn = np.unique(np.r_[np.linspace(np.log(n0), np.log(n1), sample_count), solution.t])
        n = np.exp(logn)
        y = solution.sol(logn)
        phase, op, ov, time, heat, egw, jgw, wstar, worb, dephase = y
        ell, u, torque, orb, flux, ndot, pdot, vdot, power = self.quantities(n, phase, op, ov)
        spin = .5*(self.fp*op*op+self.fv*ov*ov)
        eorb = -.5*ell*n
        # Closed budgets (all divided by I): Delta H+Q+E_GW=0 and
        # Delta(L+J_spin)+J_GW=0. Keep U: it mediates the torque-work mismatch.
        energy_residual = eorb-eorb[0]+spin-spin[0]+u-u[0]+heat+egw
        j_residual = ell-ell[0]+self.fp*(op-op[0])+self.fv*(ov-ov[0])+jgw
        # Spin-only audit has no cancellation against the huge orbital binding energy.
        spin_residual = spin-spin[0]+heat-wstar
        tau = np.full_like(ov, np.inf) if self.B == 0 else self.fp/(2*self.B*ov)
        # Inertia-weighted rms of x*Bx, then multiply by signed lag and Y.
        # This is rms over radius of the mobility-averaged drift, not an rms
        # over individual pinned/mobile vortices (which would weight Y differently).
        radial_rms = np.sqrt(np.sum(reservoir_weights(self.reservoir)*
                            (self.reservoir.columns["cylindrical_radius_m"]*
                             self.reservoir.columns["B_radial_free"])**2))
        radial_rms *= cfg["unpinned_fraction"]*(ov-op)
        out = {"time_s": time, "kepler_orbital_frequency_hz": n/(2*np.pi),
               "orbital_frequency_hz": orb/(2*np.pi), "crust_frequency_hz": op/(2*np.pi),
               "canonical_superfluid_frequency_hz": ov/(2*np.pi), "gamma_rad": phase,
               "canonical_lag_rad_s": ov-op, "heat_J": heat*self.I,
               "mutual_friction_power_W": power*self.I, "gw_energy_J": egw*self.I,
               "star_torque_work_J": wstar*self.I, "orbit_torque_work_J": worb*self.I,
               "spin_energy_change_J": (spin-spin[0])*self.I,
               "orbital_binding_energy_loss_J": (eorb[0]-eorb)*self.I,
               "interaction_energy_change_J": (u-u[0])*self.I,
               "instantaneous_heat_over_gw_power": power/flux,
               "lag_relaxation_time_s": tau,
               "signed_rms_vortex_radial_velocity_m_s": radial_rms,
               "core_neutron_transport_frequency_hz": (ov-cfg["core_entrainment_epsilon_n"]*op)/(1-cfg["core_entrainment_epsilon_n"])/(2*np.pi),
               "crust_neutron_transport_frequency_hz": (ov-cfg["inner_crust_entrainment_epsilon_n"]*op)/(1-cfg["inner_crust_entrainment_epsilon_n"])/(2*np.pi),
               "gw_phase_difference_from_point_particle_rad": dephase,
               "vortex_areal_density_m2": 2*ov/KAPPA}
        if self.fv == 0:
            out["signed_rms_vortex_radial_velocity_m_s"] *= 0
        # Measure capture at a fixed frequency window, not by assuming that
        # the supplied initial phase guarantees capture after changing physics.
        capture_window = (n/(2*np.pi) >= .8) & (n/(2*np.pi) <= 1.)
        captured = bool(np.any(capture_window) and np.ptp(phase[capture_window]) < 2*np.pi
                        and np.mean(abs(op[capture_window]-orb[capture_window])/orb[capture_window]) < .02)
        # The first stored sample at/above 1 Hz supplies a fixed phase reference.
        # A full 2pi excursion is an operational slip, not a quasistatic torque
        # limit or the first moment of nonzero crust-orbit frequency mismatch.
        ref = np.searchsorted(n/(2*np.pi), 1.)
        slip = np.flatnonzero((np.arange(len(n)) >= ref) & (abs(phase-phase[ref]) >= 2*np.pi)) if captured else []
        diagnostic = {
            "initialization": initialization, "rigid": bool(self.fv == 0),
            "captured_in_0p8_to_1_hz_window": captured,
            "mean_B": self.B, "effective_v_fraction": self.fv,
            "duration_s": float(time[-1]), "heat_J": float(heat[-1]*self.I),
            "heat_over_orbital_binding_energy_loss": float(heat[-1]/(eorb[0]-eorb[-1])),
            "heat_over_gw_radiated_energy": float(heat[-1]/egw[-1]),
            "max_instantaneous_heat_over_gw_power": float(np.max(power/flux)),
            "star_torque_work_J": float(wstar[-1]*self.I),
            "orbit_torque_work_J": float(worb[-1]*self.I),
            "spin_energy_change_J": float((spin[-1]-spin[0])*self.I),
            "orbital_binding_energy_loss_J": float((eorb[0]-eorb[-1])*self.I),
            "heat_fraction_of_net_star_work": float(heat[-1]/wstar[-1]) if wstar[-1] > 0 else None,
            "first_phase_slip_gw_hz": float(orb[slip[0]]/np.pi) if len(slip) else None,
            "final_crust_frequency_hz": float(op[-1]/(2*np.pi)),
            "final_superfluid_frequency_hz": float(ov[-1]/(2*np.pi)),
            "max_absolute_canonical_lag_rad_s": float(np.max(abs(ov-op))),
            "max_absolute_vortex_radial_velocity_m_s": float(np.max(abs(radial_rms))),
            "max_energy_balance_error_over_orbital_loss": float(np.max(abs(energy_residual))/(eorb[0]-eorb[-1])),
            "max_angular_momentum_balance_error_over_orbital_loss": float(np.max(abs(j_residual))/(ell[0]-ell[-1])),
            "max_spin_balance_error_over_work": float(np.max(abs(spin_residual))/max(np.max(abs(wstar)), 1)),
            "minimum_mutual_friction_power_W": float(np.min(power)*self.I),
            "final_gw_phase_difference_from_point_particle_rad": float(dephase[-1]),
            "solver_function_evaluations": solution.nfev,
            "solver_method": method,
            "stored_samples": len(logn)
        }
        return out, diagnostic


def reservoir_weights(reservoir):
    """Normalized canonical-inertia weights for cylindrical diagnostics."""
    return reservoir.columns["inertia_weight"]
