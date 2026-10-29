"""Radial baryon fractions for a separate, composition-resolved BSk24 star.

Run with .venv/bin/python plot_composition.py. The original APR run is untouched.
The superfluid curve is an explicit fully-paired-neutron assumption, not a
calculation of pairing gaps, finite-temperature density, or entrainment.
COMPOSITION.md explains archive fields, fraction definitions and phase handling.
"""

import hashlib
import io
import json
import os
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from ns_model import RSUN_KM, PRESSURE_TO_GEOM

ARCHIVE = ROOT / "data" / "PCP_BSK24_nparam_1.zip"
OUTPUT = ROOT / "results" / "composition_bsk24"
EXPECTED_MD5 = "57a9af0939caa9c3621baa28471b2dfd"
SOURCE = "https://zenodo.org/records/14509138"
COLORS = ["#397faa", "#d6b078", "#d76a59"]


def load_table():
    """Verify the archive and reconstruct thermodynamics and baryon fractions.

    Return nb [fm^-3], P and total E [MeV/fm^3], phase codes, five fraction
    columns, the supplied mass-radius sequence and the original archive bytes.
    Phase codes: 0 core, 2 inner crust, 1 outer crust. Fraction-column order is
    free neutrons, bound baryons, free protons, all protons, all neutrons.
    """
    raw = ARCHIVE.read_bytes()
    if hashlib.md5(raw).hexdigest() != EXPECTED_MD5:
        raise ValueError("BSk24 archive differs from the published checksum")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        nb = np.loadtxt(io.BytesIO(z.read("eos.nb")), skiprows=2)
        thermo = np.loadtxt(io.BytesIO(z.read("eos.thermo")), skiprows=1)
        comp = np.loadtxt(io.BytesIO(z.read("eos.compo")))
        mn_mev = float(z.read("eos.thermo").splitlines()[0].split()[0])
        mr = np.loadtxt(io.BytesIO(z.read("eos.mr")))
    # CompOSE q1=P/nb [MeV], q7=E/(nb*m_n)-1. E includes rest energy.
    pressure = nb * thermo[:, 3]  # MeV/fm^3
    energy = nb * mn_mev * (1 + thermo[:, 9])
    assert np.array_equal(comp[:, 1], np.arange(1, len(nb) + 1))
    assert np.array_equal(thermo[:, 1], comp[:, 1])
    assert np.all(np.diff(pressure) > 0) and np.all(np.diff(energy) > 0)
    # Use charge neutrality to include protons bound inside nuclei. In this
    # archive the nuclear quadruple abundance is 1 in the crust, so it cannot
    # be used directly as a per-baryon nuclear number fraction.
    yp = comp[:, 6] + comp[:, 8]  # Ye + Ymu = TOTAL protons / baryons
    yn_free, yp_free = comp[:, 10], comp[:, 12]
    bound = 1 - yn_free - yp_free
    bound_p = yp - yp_free
    bound_n = bound - bound_p
    assert np.min(bound_n) > -2e-8 and np.min(bound_p) > -2e-8
    assert np.max(np.abs(yn_free + bound_n + yp - 1)) < 2e-8
    fractions = np.column_stack((yn_free, np.clip(bound, 0, 1), yp_free, yp, 1-yp))
    return nb, pressure, energy, comp[:, 3].astype(int), fractions, mr, raw


def make_solver(nb, pressure, energy):
    """Build an independent BSk24 TOV solver, not an APR composition overlay.

    Return a shooting function and two log-space density/pressure interpolants.
    All PCHIP interpolants forbid extrapolation; the surface is the lowest
    tabulated pressure, not the APR pseudo-enthalpy cutoff or exactly zero P.
    """
    # 1 MeV/fm^3 = 1.602176634e32 Pa. Convert to G=c=1, length in km.
    mev_to_geom = 1.602176634e32 * PRESSURE_TO_GEOM
    p, e = pressure * mev_to_geom, energy * mev_to_geom
    logp = np.log(p)
    energy_of_logp = PchipInterpolator(logp, np.log(e), extrapolate=False)
    nb_of_logp = PchipInterpolator(logp, np.log(nb), extrapolate=False)
    logp_of_logn = PchipInterpolator(np.log(nb), logp, extrapolate=False)

    def solve(nc, rtol=2e-9, step=0.025):
        """Integrate [m (km), log P (geometric km^-2)] for central nb=nc."""
        pc = np.exp(logp_of_logn(np.log(nc)))
        ec = np.exp(energy_of_logp(np.log(pc)))
        r0 = 1e-4
        # Regular-centre TOV expansion in log pressure. Central energy includes
        # rest energy; replacing it by baryon density would change the solution.
        initial = [4*np.pi*ec*r0**3/3,
                   np.log(pc) - 2*np.pi*(ec+pc)*(ec/3+pc)*r0*r0/pc]

        def rhs(r, state):
            mass, lp = state
            # Trial solver stages may cross the terminal event. Clip library
            # queries to the tabulated domain; the event defines the endpoint.
            lp = np.clip(lp, logp[0], logp[-1])
            pp, ee = np.exp(lp), np.exp(energy_of_logp(lp))
            return [4*np.pi*r*r*ee,
                    -(ee+pp)/pp*(mass+4*np.pi*r**3*pp)/(r*(r-2*mass))]

        def surface(r, state):
            return state[1] - logp[0]

        surface.terminal = True
        surface.direction = -1
        sol = solve_ivp(rhs, (r0, 30), initial, method="DOP853", rtol=rtol,
                        atol=[1e-12, 1e-10], max_step=step,
                        events=surface, dense_output=True)
        if not sol.success or len(sol.t_events[0]) != 1:
            raise RuntimeError(sol.message)
        return sol

    return solve, logp_of_logn, nb_of_logp


def radial_fractions(n, nb, phase, fractions):
    """Interpolate inside phases; retain a step at each phase boundary.

    The source brackets phase transitions rather than resolving their exact
    coexistence pressure. Use geometric density midpoints, report both brackets,
    and hold the last same-phase composition to the bracket midpoint.
    """
    drip = np.sqrt(nb[phase == 1][-1] * nb[phase == 2][0])
    cc = np.sqrt(nb[phase == 2][-1] * nb[phase == 0][0])
    out = np.empty((len(n), fractions.shape[1]))
    for region, mask in [(0, n >= cc), (2, (n >= drip) & (n < cc)), (1, n < drip)]:
        src = phase == region
        for j in range(fractions.shape[1]):
            out[mask, j] = np.interp(np.log(n[mask]), np.log(nb[src]), fractions[src, j])
    return out, drip, cc


def draw(radius, frac, rcc, rdrip, surface):
    """Plot two partitions of local baryon number, with a crust enlargement.

    The top grouping counts bound neutrons in the lattice; the lower grouping
    counts them as neutrons. Neither plot shows mass or inertia fractions.
    """
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.8), sharey=True,
                             gridspec_kw={"width_ratios": [1.25, 1]})
    fig.subplots_adjust(left=.07, right=.98, bottom=.16, top=.79, hspace=.30, wspace=.13)
    fig.suptitle("What is the star made of at each radius?", x=.07, y=.975,
                 ha="left", fontsize=22, weight="bold")
    fig.text(.07, .925, rf"BSk24  |  $M=1.4\,M_\odot$  |  $R={surface:.2f}$ km"
             "     •     Percentages count baryons locally", fontsize=12, color="#455263")
    fig.text(.07, .879, "Superfluid illustration: assume ALL unbound neutrons are paired. "
             "Lattice = baryons bound in nuclear clusters.", fontsize=11,
             bbox={"boxstyle": "round,pad=0.5", "fc": "#edf3f8", "ec": "none"})
    top_labels = ["Unbound neutrons (assumed superfluid)",
                  "Lattice: bound neutrons + protons", "Unbound protons"]
    for col in range(2):
        ax = axes[0, col]
        ax.stackplot(radius, *(100*frac[:, i] for i in range(3)),
                     colors=COLORS, labels=top_labels, linewidth=0)
        ax = axes[1, col]
        ax.stackplot(radius, 100*frac[:, 4], 100*frac[:, 3],
                     colors=[COLORS[0], COLORS[2]],
                     labels=["All neutrons (bound + unbound)", "All protons (bound + unbound)"],
                     linewidth=0)
        for ax in axes[:, col]:
            ax.set(xlim=(0 if col == 0 else rcc-.13, surface), ylim=(0, 100))
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.set_axisbelow(False)
            ax.grid(axis="y", color="white", alpha=.27, linewidth=.7)
            for boundary in (rcc, rdrip):
                ax.axvline(boundary, color="#394655", ls="--", lw=1, alpha=.8)
            ax.set_xlabel("Radius from centre (km)")
        axes[0, col].set_title("Whole star" if col == 0 else "Crust enlarged", loc="left", weight="bold", pad=25)
        for lo, hi, text in [(0, rcc, "Core"), (rcc, rdrip, "Inner crust"),
                            (rdrip, surface, "Outer crust")]:
            if col == 1 and lo == 0:
                continue
            if col == 0 and lo != 0:
                continue
            axes[0, col].text((lo+hi)/2, 1.015, text,
                              transform=axes[0, col].get_xaxis_transform(),
                              ha="center", fontsize=10, color="#455263")
    axes[0, 0].set_ylabel("Dynamical grouping\nBaryon fraction (%)")
    axes[1, 0].set_ylabel("Particle identity\nBaryon fraction (%)")
    axes[0, 0].legend(loc="lower left", fontsize=9, framealpha=.94)
    axes[1, 0].legend(loc="lower left", fontsize=9, framealpha=.94)
    axes[0, 1].text((rdrip+surface)/2, 50, "100%\nbound in nuclei", ha="center", va="center", fontsize=11)
    fig.text(.07, .065, "Bound neutrons count as neutrons AND as part of the lattice. "
             "Electrons and muons carry no baryon number and are excluded from these percentages.", fontsize=10, color="#455263")
    fig.text(.07, .037, "Pairing gaps, temperature and entrainment are not solved. "
             "BSk24 composition example; this is a separate model from the workspace's APR star.", fontsize=10, color="#455263")
    fig.text(.07, .009, "Source: Pearson et al. (2018, 2019), CompOSE PCP(BSK24), Zenodo 14509138. "
             "Surface at the lowest tabulated pressure; transitions located within table brackets.", fontsize=9, color="#66717e")
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUTPUT / f"radial_composition.{ext}", dpi=190, facecolor="white")
    plt.close(fig)


def main():
    """Shoot to 1.4 solar masses, map phase brackets and export checked fractions."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    nb, pressure, energy, phase, fractions, mr, raw = load_table()
    solve, lp_of_ln, ln_of_lp = make_solver(nb, pressure, energy)
    nc = brentq(lambda x: solve(x).y[0, -1]/RSUN_KM-1.4, .2, .9, xtol=1e-11)
    star = solve(nc)
    surface = float(star.t_events[0][0])
    _, drip, cc = radial_fractions(np.array([nc]), nb, phase, fractions)

    def radius_at_n(n):
        """Invert the same BSk24 P(nb) relation on the dense TOV solution."""
        lp = float(lp_of_ln(np.log(n)))
        return brentq(lambda r: star.sol(r)[1]-lp, star.t[0], surface)

    rcc, rdrip = radius_at_n(cc), radius_at_n(drip)
    # Include tabulated composition knots and both sides of phase boundaries.
    knots = [radius_at_n(n) for n in nb if nb[0] < n < nc]
    radius = np.unique(np.r_[0, np.linspace(star.t[0], surface, 2200),
                              np.linspace(rcc-.13, surface, 1600), knots,
                              rcc-1e-8, rcc+1e-8, rdrip-1e-8, rdrip+1e-8])
    state = star.sol(np.maximum(radius, star.t[0]))
    n = np.exp(ln_of_lp(np.clip(state[1], float(lp_of_ln(np.log(nb[0]))),
                                    float(lp_of_ln(np.log(nb[-1]))))))
    n[0] = nc
    frac, _, _ = radial_fractions(n, nb, phase, fractions)
    # Hold nc fixed for convergence; re-shooting would mix integration error
    # with a change of central state. The supplied eos.mr comparison is separate
    # and may differ because of its interpolation and surface prescription.
    fine = solve(nc, rtol=4e-10, step=.0125)
    stable_mr = mr[:np.argmax(mr[:, 1])+1]
    ref_radius = float(np.interp(1.4, stable_mr[:, 1], stable_mr[:, 0]))
    validation = {
        "target_mass_msun": float(star.y[0, -1]/RSUN_KM),
        "target_mass_absolute_error_msun": float(abs(star.y[0, -1]/RSUN_KM-1.4)),
        "refined_radius_change_m": float(abs(fine.t[-1]-surface)*1000),
        "refined_mass_change_msun": float(abs(fine.y[0, -1]-star.y[0, -1])/RSUN_KM),
        "published_eos_mr_radius_km": ref_radius,
        "radius_difference_from_eos_mr_m": float((surface-ref_radius)*1000),
        "maximum_group_fraction_sum_error": float(np.max(abs(frac[:, :3].sum(axis=1)-1))),
        "maximum_identity_fraction_sum_error": float(np.max(abs(frac[:, 3:].sum(axis=1)-1))),
        "all_fractions_in_unit_interval": bool(np.all((frac >= -2e-8) & (frac <= 1+2e-8))),
        "density_decreases": bool(np.all(np.diff(n) <= 1e-10)),
    }
    # Validate before replacing composition CSV/summary/figures. These checks
    # test fraction bookkeeping and the chosen TOV numerics, not pairing physics.
    assert validation["target_mass_absolute_error_msun"] < 1e-7
    assert validation["refined_radius_change_m"] < 1
    assert abs(validation["radius_difference_from_eos_mr_m"]) < 100
    assert validation["maximum_group_fraction_sum_error"] < 2e-8
    assert validation["all_fractions_in_unit_interval"] and validation["density_decreases"]
    validation["passed"] = True
    columns = np.column_stack((radius, n, frac))
    np.savetxt(OUTPUT / "radial_composition.csv", columns, delimiter=",", comments="",
               header="radius_km,baryon_density_fm-3,unbound_neutron_assumed_sf_fraction,"
                      "lattice_bound_baryon_fraction,unbound_proton_fraction,total_proton_fraction,total_neutron_fraction")
    summary = {"eos": "PCP(BSK24)", "source": SOURCE,
               "archive_md5": hashlib.md5(raw).hexdigest(),
               "archive_sha256": hashlib.sha256(raw).hexdigest(),
               "radius_km": surface, "central_baryon_density_fm-3": nc,
               "core_crust_radius_km": rcc, "neutron_drip_radius_km": rdrip,
               "central_neutron_percent": float(frac[0, 4]*100),
               "central_proton_percent": float(frac[0, 3]*100),
               "phase_transition_density_brackets_fm-3": {
                   "core_crust": [float(nb[phase == 2][-1]), float(nb[phase == 0][0])],
                   "neutron_drip": [float(nb[phase == 1][-1]), float(nb[phase == 2][0])]},
               "phase_transition_radius_brackets_km": {
                   "core_crust": [radius_at_n(nb[phase == 0][0]), radius_at_n(nb[phase == 2][-1])],
                   "neutron_drip": [radius_at_n(nb[phase == 2][0]), radius_at_n(nb[phase == 1][-1])]},
               "surface_pressure_mev_fm-3": float(pressure[0]),
               "assumptions": ["Local baryon number fractions; not volume or stellar mass fractions.",
                   "All unbound neutrons assumed paired; no temperature/gap/entrainment calculation.",
                   "Lattice fraction is bound baryons, not nuclear volume filling factor.",
                   "Total proton fraction follows charge neutrality Yp=Ye+Ymu.",
                   "Composition interpolated within phases, transitions at geometric density bracket midpoints.",
                   "Star ends at the lowest tabulated pressure; no atmosphere extrapolation.",
                   "Separate BSk24 star; existing APR data and configuration are unchanged."],
               "validation": validation}
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    draw(radius, frac, rcc, rdrip, surface)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
