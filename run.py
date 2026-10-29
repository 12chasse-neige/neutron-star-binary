#!/usr/bin/env python3
"""Build the APR background, classify layers and export numerical evidence.

Review sequence: prepare_eos -> solve_target (ns_model.py) -> make_summary ->
validate -> exports/plots/report. MODEL.md defines equations and units.
Binary-context inputs label the result; this script evolves no binary orbit.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import lalsimulation as lalsim
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import PchipInterpolator

from ns_model import (EquationOfState, solve_target, integrate, G, C, MSUN,
                      RSUN_KM, PRESSURE_TO_GEOM, DENSITY_TO_GEOM)

COLORS = {"core": "#335d8b", "inner crust": "#e6a34b", "outer crust": "#b4d7c3"}


def write_csv(path, columns):
    """Preserve dictionary column order and the full supplied sampling grid."""
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def monotone_table(raw, path, points):
    """Resample log E(log P) with PCHIP; preserve all original nodes exactly.

    ``points`` counts the added logarithmic grid, not the final union of nodes.
    Both input and output columns remain geometric m^-2, as expected by LAL.
    This changes interpolation between nodes, not the upstream table itself.
    """
    dense_p = np.unique(np.r_[np.geomspace(raw[0, 0], raw[-1, 0], points), raw[:, 0]])
    pchip = PchipInterpolator(np.log(raw[:, 0]), np.log(raw[:, 1]))
    dense_e = np.exp(pchip(np.log(dense_p)))
    dense_e[np.searchsorted(dense_p, raw[:, 0])] = raw[:, 1]
    np.savetxt(path, np.column_stack([dense_p, dense_e]), fmt="%.18e")
    return len(dense_p)


def prepare_eos(name):
    """Verify the frozen APR input and rebuild its dense integration table.

    If absent, copy the installed table and record its provenance. Existing
    original-table metadata supplies the integrity hash. The returned metadata
    records this run's library version and derived-table hash separately.
    """
    if name != "APR":
        raise ValueError("This reproduction project supports APR; do not silently substitute APR4.")
    target = ROOT / "data" / "LALSimNeutronStarEOS_APR.dat"
    target.parent.mkdir(exist_ok=True)
    if not target.exists():
        installed = next(Path(lalsim.__file__).parent.parent.rglob(target.name))
        shutil.copyfile(installed, target)
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    metadata = {"table": target.name, "sha256": digest,
                "distribution": "lalsuite", "distribution_version": importlib.metadata.version("lalsuite"),
                "lalsimulation_version": lalsim.__version__,
                "columns": ["pressure_geometric_m^-2", "energy_density_geometric_m^-2"],
                "upstream": "https://git.ligo.org/lscsoft/lalsuite",
                "documentation": "https://lscsoft.docs.ligo.org/lalsuite/dev/lalsimulation/group___l_a_l_sim_neutron_star_e_o_s__c.html"}
    metadata_path = target.with_suffix(".provenance.json")
    if metadata_path.exists():
        saved = json.loads(metadata_path.read_text())
        if saved["sha256"] != digest:
            raise ValueError("Frozen EOS checksum has changed.")
    else:
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    # The native cubic spline overshoots around a coarse inner-crust join.
    # PCHIP in log(P),log(E) preserves the monotone original nodes. A dense
    # sampling lets LAL reconstruct h and rho thermodynamically, while making
    # its subsequent cubic interpolation negligible. Keep every original node.
    raw = np.loadtxt(target)
    if not np.all(np.diff(raw, axis=0) > 0):
        raise ValueError("Original EOS table must be strictly monotone.")
    dense = target.with_name("APR_monotone_dense.dat")
    rows = monotone_table(raw, dense, 24000)
    metadata["original_table_rows"] = len(raw)
    metadata["integration_table"] = dense.name
    metadata["integration_table_rows"] = rows
    metadata["integration_table_sha256"] = hashlib.sha256(dense.read_bytes()).hexdigest()
    metadata["interpolation"] = "Monotone PCHIP log(E) versus log(P), densely sampled with every original node retained; LAL reconstructs h and rho."
    dense.with_suffix(".provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return EquationOfState(name, dense), metadata


def make_summary(star, cfg, provenance):
    """Convert observables and classify shells without changing the solution.

    Shell Delta m is gravitational mass, while shell labels use rest density.
    Proper thickness comes from the sixth ODE state, not coordinate subtraction.
    The prescribed epsilon*I is contextual; it does not enter the TOV solver.
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
              for i, name in enumerate(COLORS)]
    # I_geom is a length cubed: km^3 -> m^3, then c^2/G -> kg m^2.
    inertia = star.inertia_geom_km3 * 1e9 * C*C/G
    context = dict(cfg["binary_context"])
    context["companion_mass_msun"] = m * context["mass_ratio_companion_over_ns"]
    eps = context["intrinsic_equatorial_ellipticity"]
    return {"model": "cold spherical TOV + first-order Hartle rotation + static quadrupolar tides",
            "eos": provenance, "binary_context": context,
            "mass_msun": m, "radius_km": r, "compactness_GM_Rc2": star.surface[0]/r,
            "central_pseudo_enthalpy": star.central_h,
            "central_pressure_pa": pc/PRESSURE_TO_GEOM,
            "central_pressure_dyn_cm2": 10*pc/PRESSURE_TO_GEOM,
            "central_rest_mass_density_g_cm3": rhoc/1000,
            "central_energy_equivalent_density_g_cm3": ec/DENSITY_TO_GEOM/1000,
            "moment_of_inertia_kg_m2": inertia,
            "moment_of_inertia_g_cm2": inertia*1e7,
            "moment_of_inertia_geometric_msun3": star.inertia_geom_km3/RSUN_KM**3,
            "prescribed_equatorial_delta_I_kg_m2": eps*inertia,
            "love_k2": star.love_k2, "tidal_deformability_lambda": star.tidal_lambda,
            "paper_reference_I_msun3": 30.18, "paper_reference_lambda": 248,
            "I_fractional_difference_from_paper": star.inertia_geom_km3/RSUN_KM**3/30.18-1,
            "lambda_fractional_difference_from_paper": star.tidal_lambda/248-1,
            "layer_thresholds_rest_mass_g_cm3": {
                "core_crust": cfg["core_crust_rest_mass_density_g_cm3"],
                "neutron_drip": cfg["neutron_drip_rest_mass_density_g_cm3"]},
            "core_crust_boundary": cc, "neutron_drip_boundary": drip, "layers": layers,
            "surface": {"pseudo_enthalpy_cutoff": star.surface_h,
                        "pressure_pa": star.eos.at_h(star.surface_h)[0]/PRESSURE_TO_GEOM},
            "boundary_interpretation": "Configurable rest-mass-density classification, not a predicted nuclear phase transition."}


def validate(star, cfg, profile):
    """Audit the equations numerically; do not fit the paper's observables.

    Refinements hold central enthalpy fixed so they measure discretization,
    rather than hiding changes by re-shooting to the same target mass. The LAL
    comparison is an independent ODE solver using the same derived EOS; it does
    not independently establish the EOS or its low-density physical accuracy.
    """
    fine = integrate(star.eos, star.central_h, surface_h=star.surface_h,
                     rtol=cfg["relative_tolerance"]/5, max_step_km=cfg["max_step_km"]/2,
                     start_radius_km=1e-5)
    lower_surface = integrate(star.eos, star.central_h, surface_h=star.surface_h/10,
                              rtol=cfg["relative_tolerance"], max_step_km=cfg["max_step_km"])
    pc = star.eos.at_h(star.central_h)[0]/PRESSURE_TO_GEOM
    ref_r, ref_m, ref_k2 = lalsim.SimNeutronStarTOVODEIntegrateWithTolerance(pc, star.eos.raw, 1e-9)
    dedp = np.array([star.eos.dedp(x*PRESSURE_TO_GEOM) for x in profile["pressure_pa"]])
    # A log-h grid checks the EOS between the radial samples, including the
    # thin outer layers and the inner-crust interpolation interval.
    sample_h = np.geomspace(star.surface_h, star.central_h, 12000)
    sample_eos = np.array([star.eos.at_h(h) for h in sample_h])
    sample_dedp = np.array([star.eos.dedp(p) for p in sample_eos[:, 0]])
    mass_integral = simpson(profile["mass_per_radius_msun_per_km"], x=profile["radius_km"])
    lo = integrate(star.eos, star.central_h*0.999, rtol=2e-8, max_step_km=0.04).mass_msun
    hi = integrate(star.eos, star.central_h*1.001, rtol=2e-8, max_step_km=0.04).mass_msun
    checks = {
        "target_mass_relative_error": abs(star.mass_msun/cfg["target_mass_msun"] - 1),
        "lal_reference_radius_relative_error": abs(star.radius_km/(ref_r/1000) - 1),
        "lal_reference_mass_relative_error": abs(star.mass_msun/(ref_m/MSUN) - 1),
        "lal_reference_k2_relative_error": abs(star.love_k2/ref_k2 - 1),
        "refined_radius_relative_change": abs(fine.radius_km/star.radius_km - 1),
        "refined_mass_relative_change": abs(fine.mass_msun/star.mass_msun - 1),
        "refined_inertia_relative_change": abs(fine.inertia_geom_km3/star.inertia_geom_km3 - 1),
        "refined_lambda_relative_change": abs(fine.tidal_lambda/star.tidal_lambda - 1),
        "surface_cutoff_radius_change_m": abs(lower_surface.radius_km-star.radius_km)*1000,
        "mass_integral_relative_error": abs(mass_integral/star.mass_msun - 1),
        "mass_nondecreasing": bool(np.all(np.diff(profile["enclosed_mass_msun"]) >= -1e-11)),
        "pressure_decreasing": bool(np.all(np.diff(profile["pressure_pa"]) < 0)),
        "density_decreasing": bool(np.all(np.diff(profile["rest_mass_density_g_cm3"]) < 0)),
        "stable_branch_dM_dhc_positive": bool(hi > lo),
        "maximum_sound_speed_squared_over_c2_in_star": float(np.max(1/dedp)),
        "minimum_sound_speed_squared_over_c2_in_star": float(np.min(1/dedp)),
        "eos_scan_density_increasing_with_h": bool(np.all(np.diff(sample_eos[:, 2]) > 0)),
        "eos_scan_positive_compressibility": bool(np.all(sample_dedp > 0)),
        "eos_scan_causal": bool(np.all(1/sample_dedp <= 1)),
        "solution_finite": bool(all(np.all(np.isfinite(v)) for v in profile.values())),
    }
    cc = cfg["core_crust_rest_mass_density_g_cm3"]
    nd = cfg["neutron_drip_rest_mass_density_g_cm3"]
    checks["refined_core_boundary_change_m"] = abs(fine.boundary(cc)["radius_km"]-star.boundary(cc)["radius_km"])*1000
    checks["refined_neutron_drip_change_m"] = abs(fine.boundary(nd)["radius_km"]-star.boundary(nd)["radius_km"])*1000
    # Native coarse-spline differences describe a different interpolation and
    # are recorded as diagnostics, not used as convergence acceptance targets.
    native_eos = EquationOfState("APR", ROOT / "data" / "LALSimNeutronStarEOS_APR.dat")
    native_r, native_m, native_k2 = lalsim.SimNeutronStarTOVODEIntegrateWithTolerance(pc, native_eos.raw, 1e-9)
    checks["native_coarse_interpolant_radius_difference_m_at_same_pc"] = (star.radius_km*1000-native_r)
    checks["native_coarse_interpolant_mass_fractional_difference_at_same_pc"] = star.mass_msun/(native_m/MSUN)-1
    checks["native_coarse_interpolant_k2_fractional_difference_at_same_pc"] = star.love_k2/native_k2-1
    with tempfile.TemporaryDirectory(prefix="ns-eos-refinement-") as temporary:
        dense_path = Path(temporary) / "APR_refined.dat"
        original = np.loadtxt(ROOT / "data" / "LALSimNeutronStarEOS_APR.dat")
        monotone_table(original, dense_path, 48000)
        refined_eos = EquationOfState("APR", dense_path)
        refined_star = integrate(refined_eos, star.central_h, surface_h=star.surface_h,
                                 rtol=cfg["relative_tolerance"], max_step_km=cfg["max_step_km"])
    checks["doubled_eos_sampling_radius_relative_change"] = abs(refined_star.radius_km/star.radius_km-1)
    checks["doubled_eos_sampling_mass_relative_change"] = abs(refined_star.mass_msun/star.mass_msun-1)
    checks["doubled_eos_sampling_inertia_relative_change"] = abs(refined_star.inertia_geom_km3/star.inertia_geom_km3-1)
    checks["doubled_eos_sampling_lambda_relative_change"] = abs(refined_star.tidal_lambda/star.tidal_lambda-1)
    limits = {"target_mass_relative_error": 1e-8, "lal_reference_radius_relative_error": 2e-5,
              "lal_reference_mass_relative_error": 2e-5, "lal_reference_k2_relative_error": 2e-4,
              "refined_radius_relative_change": 2e-6, "refined_mass_relative_change": 2e-6,
              "refined_inertia_relative_change": 2e-6, "refined_lambda_relative_change": 2e-5,
              "surface_cutoff_radius_change_m": 0.01, "mass_integral_relative_error": 2e-5,
              "refined_core_boundary_change_m": 0.05, "refined_neutron_drip_change_m": 0.05,
              "doubled_eos_sampling_radius_relative_change": 2e-6,
              "doubled_eos_sampling_mass_relative_change": 2e-6,
              "doubled_eos_sampling_inertia_relative_change": 2e-6,
              "doubled_eos_sampling_lambda_relative_change": 2e-5}
    passed = (all(checks[k] < limit for k, limit in limits.items())
              and all(checks[k] for k in ["mass_nondecreasing", "pressure_decreasing", "density_decreasing",
                                          "stable_branch_dM_dhc_positive", "solution_finite",
                                          "eos_scan_density_increasing_with_h",
                                          "eos_scan_positive_compressibility", "eos_scan_causal"])
              and 0 < checks["minimum_sound_speed_squared_over_c2_in_star"]
              and checks["maximum_sound_speed_squared_over_c2_in_star"] <= 1)
    return {"passed": bool(passed), "checks": checks, "acceptance_limits": limits,
            "reference": "Independent LAL TOV/Love ODE solver, same frozen EOS and central pressure",
            "convergence": "Five times tighter relative tolerance, half maximum step, ten times smaller starting radius; fixed central enthalpy",
            "scope": "Numerical checks do not quantify EOS or transition-density uncertainty."}


def plot_results(profile, summary, out):
    """Draw classified shells and physical profiles; log limits affect plots only."""
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": "#29313d", "text.color": "#29313d",
                         "axes.titleweight": "bold", "figure.facecolor": "white"})
    r = profile["radius_km"]
    R = summary["radius_km"]
    cc = summary["core_crust_boundary"]["radius_km"]
    nd = summary["neutron_drip_boundary"]["radius_km"]
    layers = summary["layers"]
    fig = plt.figure(figsize=(14, 10.5), layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[0.65, 1, 1])
    full = fig.add_subplot(gs[0, 0])
    zoom = fig.add_subplot(gs[0, 1])
    for ax, low in [(full, 0), (zoom, cc)]:
        for layer in layers:
            a, b = layer["inner_radius_km"], layer["outer_radius_km"]
            ax.barh(0.4, b-a, left=a, height=0.38, color=COLORS[layer["name"]], edgecolor="white")
        ax.set(ylim=(0, 0.8), xlim=(low, R), yticks=[], xlabel="Areal radius r [km]")
        ax.spines[["left", "bottom"]].set_visible(False)
        ax.grid(axis="x", alpha=0.15)
    full.set_title("Spherical star · radial layers", loc="left", pad=14)
    full.text(cc/2, 0.4, f"Core\n0–{cc:.2f} km", ha="center", va="center", color="white", weight="bold")
    zoom.set_title("Crust enlarged · same radial coordinates", loc="left", pad=14)
    zoom.text((cc+nd)/2, 0.4, f"Inner crust\n{nd-cc:.3f} km thick", ha="center", va="center", fontsize=9)
    zoom.text((nd+R)/2, 0.4, f"Outer crust\n{R-nd:.3f} km thick", ha="center", va="center", fontsize=9)
    zoom.set_xticks([cc, nd, R], [f"{cc:.3f}\ncore edge", f"{nd:.3f}\nneutron drip", f"{R:.3f}\nsurface"])
    axes = [fig.add_subplot(gs[i, j]) for i in [1, 2] for j in [0, 1]]
    fields = [("enclosed_mass_msun", "Enclosed gravitational mass", r"m(r) [$M_\odot$]", False, 1),
              ("pressure_pa", "Hydrostatic pressure", "P(r) [Pa]", True, 1),
              ("rest_mass_density_g_cm3", "Matter density", r"Density [g cm$^{-3}$]", True, 1),
              ("mass_per_radius_msun_per_km", "Mass in a radial shell", r"dm/dr [$M_\odot$ km$^{-1}$]", False, 1)]
    for ax, (field, title, ylabel, log, scale) in zip(axes, fields):
        for layer in layers:
            ax.axvspan(layer["inner_radius_km"], layer["outer_radius_km"], color=COLORS[layer["name"]], alpha=0.14)
        ax.plot(r, profile[field]*scale, color="#273e59", lw=2)
        for x in [cc, nd]:
            ax.axvline(x, color="#78828e", ls=":", lw=1)
        ax.set(xlim=(0, R), xlabel="Areal radius r [km]", ylabel=ylabel, title=title)
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.grid(alpha=0.18)
        if log:
            ax.set_yscale("log")
        else:
            ax.set_ylim(bottom=0)
    # The dilute numerical surface is excluded from the displayed log range only.
    axes[1].set_ylim(1e23, summary["central_pressure_pa"]*2)
    axes[2].set_ylim(1e7, summary["central_energy_equivalent_density_g_cm3"]*2)
    axes[1].text(0.03, 0.05, r"Dilute surface below $10^{23}$ Pa omitted from display", transform=axes[1].transAxes, fontsize=8, color="#5a626c")
    axes[2].plot(r, profile["energy_equivalent_density_g_cm3"], color="#b35f37", lw=1.4, ls="--", label=r"Total energy / $c^2$")
    axes[2].lines[0].set_label("Rest mass density")
    axes[2].legend(frameon=False, loc="lower left", fontsize=9)
    axes[0].text(0.05, 0.94, f"Core contains {100*layers[0]['mass_fraction']:.2f}% of M\nCrust contains {100*(1-layers[0]['mass_fraction']):.2f}% of M",
                 transform=axes[0].transAxes, va="top", fontsize=10)
    thresholds = summary["layer_thresholds_rest_mass_g_cm3"]
    fig.suptitle(f"APR neutron star · {summary['mass_msun']:.2f} solar masses · R = {R:.3f} km\n"
                 f"Density-defined layers: core edge at {thresholds['core_crust']:.1e}, neutron drip at {thresholds['neutron_drip']:.1e} g/cm³",
                 fontsize=17)
    fig.savefig(out / "structure_profiles.png", dpi=200)
    fig.savefig(out / "structure_profiles.svg")
    plt.close(fig)


def report(summary, validation, sensitivity):
    """Render the saved observables and assumptions as a Markdown review report."""
    s = summary
    rows = "\n".join(f"| {x['name'].title()} | {x['inner_radius_km']:.5f}–{x['outer_radius_km']:.5f} | {x['coordinate_thickness_km']:.5f} | {x['gravitational_mass_contribution_msun']:.8f} | {100*x['mass_fraction']:.5f}% |" for x in s["layers"])
    sens_rows = "\n".join(f"| {x['threshold_g_cm3']:.2e} | {x['core_radius_km']:.5f} | {x['total_crust_thickness_km']:.5f} |" for x in sensitivity)
    return rf"""# APR neutron-star structure: simulation results

The governing equations and implementation map are in [MODEL.md](../MODEL.md).
In geometric units, with total energy density $\mathcal E$,

$$
\frac{{dm}}{{dr}}=4\pi r^2\mathcal E,\qquad
\frac{{dh}}{{dr}}=-\frac{{m+4\pi r^3P}}{{r(r-2m)}},\qquad
I=\frac J\Omega,\qquad \Lambda=\frac{{2k_2}}{{3C_\star^5}}.
$$

Generated by `run.py`; exact inputs are in `config_used.json`. Numerical checks: **{'PASS' if validation['passed'] else 'FAIL'}**.

## Connection to the probing paper

Miao & Yang, *Probing Intrinsic Ellipticity in Compact Star Binaries*, local arXiv:2511.08895v2 (20 Aug 2026), Fig. 2 and Sec. IV.
The Fig. 2 configuration is $M=1.4\,M_\odot$, $q=M'/M=3$, $\epsilon=5\times10^{{-5}}$,
initial orbital frequency 0.4 Hz, initial spin frequency 0.5 Hz, and $f-\theta=0.3986\,\mathrm{{rad}}$.
Here $M'=4.2\,M_\odot$. Sec. IV adopts APR and quotes $I=30.18$ $M_\odot^3$
(geometric units) and $\Lambda=248$ for a 1.4-solar-mass star. Its injection $q=5$
is a separate example; it is not substituted for Fig. 2's $q=3$.

This calculation builds the **spherical equilibrium background**, rather than evolving
the orbital/spin tracks of Fig. 2. The ellipticity, binary mass ratio and phase are retained
as context; they do not enter the spherical TOV equations. The prescribed equatorial
inertia difference $\epsilon I$ is {s['prescribed_equatorial_delta_I_kg_m2']:.5e} $\mathrm{{kg\,m^2}}$.
No magnetic/elastic mechanism is solved to produce that deformation.

## Global results

- Gravitational mass: {s['mass_msun']:.10f} solar masses.
- Areal/circumferential radius: **{s['radius_km']:.5f} km**.
- Central pressure: **{s['central_pressure_pa']:.6e} Pa** = {s['central_pressure_dyn_cm2']:.6e} $\mathrm{{dyn\,cm^{{-2}}}}$.
- Central rest mass density: {s['central_rest_mass_density_g_cm3']:.6e} $\mathrm{{g\,cm^{{-3}}}}$.
- Central energy-equivalent density $\mathcal E/c^2$: {s['central_energy_equivalent_density_g_cm3']:.6e} $\mathrm{{g\,cm^{{-3}}}}$.
- Relativistic moment of inertia: {s['moment_of_inertia_geometric_msun3']:.6f} $M_\odot^3$ = {s['moment_of_inertia_g_cm2']:.6e} $\mathrm{{g\,cm^2}}$.
- Difference from the paper's quoted $I$: {100*s['I_fractional_difference_from_paper']:+.4f}%.
- Tidal deformability $\Lambda$: {s['tidal_deformability_lambda']:.5f}; $k_2$={s['love_k2']:.7f}.
- Difference from the paper's quoted $\Lambda$: {100*s['lambda_fractional_difference_from_paper']:+.4f}%.

The EOS family and benchmark observables agree closely with the paper. The exact EOS
table, interpolation and crust matching used by the authors are not specified there;
the small residual difference is not fitted away or claimed to have a known cause.
The APR table is used, not an APR4 piecewise-polytropic approximation.
A monotone PCHIP interpolant in $\ln P,\ln\mathcal E$ passes through every original
table point. It is densely sampled before LAL reconstructs enthalpy and rest mass
density. This removes a spurious negative-compressibility interval in LAL's native
coarse-table cubic spline near the inner crust. The original data are preserved;
both the original and derived tables and their checksums are saved in `data/`.
The comparison to LAL below uses the same derived EOS. Differences from its
native coarse interpolant at fixed central pressure are also recorded separately.

## Radial layers

These are **density-defined regions of the numerical solution**. The core-crust boundary
is set at rest mass density {s['layer_thresholds_rest_mass_g_cm3']['core_crust']:.3e} $\mathrm{{g\,cm^{{-3}}}}$ and neutron drip at {s['layer_thresholds_rest_mass_g_cm3']['neutron_drip']:.3e} $\mathrm{{g\,cm^{{-3}}}}$.
The neutron-drip scale follows Chamel & Haensel (2008); the exact core transition is
an additional configurable modeling choice. Neither boundary is inferred from $\epsilon$,
computed by a nuclear stability calculation, or specified by the probing paper.

| Region | Radial interval (km) | Coordinate thickness (km) | Gravitational mass contribution (solar masses) | Fraction of total $M$ |
|---|---:|---:|---:|---:|
{rows}

At the core edge $P$ = {s['core_crust_boundary']['pressure_pa']:.6e} Pa and
$m$ = {s['core_crust_boundary']['enclosed_mass_msun']:.8f} solar masses.
At neutron drip $P$ = {s['neutron_drip_boundary']['pressure_pa']:.6e} Pa and
$m$ = {s['neutron_drip_boundary']['enclosed_mass_msun']:.8f} solar masses.

Radii are Schwarzschild areal coordinates. A shell's proper radial thickness is
$\int dr/\sqrt{{1-2Gm/(rc^2)}}$; this different quantity is also saved in `layers.csv`.
Layer masses are differences in enclosed **gravitational mass**, not baryonic masses
or proper-volume integrals of rest mass density. The radial mass distribution uses
total energy density including internal energy. `profiles.csv` distinguishes both densities.

## Boundary sensitivity

Only the classification threshold is changed below; the EOS and equilibrium star stay fixed.
This illustrates uncertainty in assigning the core edge; it is **not** a crust-EOS uncertainty study.

| Core threshold ($\mathrm{{g\,cm^{{-3}}}}$) | Core radius (km) | Total coordinate crust thickness (km) |
|---|---:|---:|
{sens_rows}

## Numerical validation

The implementation is checked against LAL's independent TOV/Love integrator at the
same central pressure. Resolution is checked with a smaller step, tighter tolerance
and smaller central starting radius. Lowering the numerical surface enthalpy cutoff
by ten tests the unintegrated dilute atmosphere. The surface pressure is small but
finite ({s['surface']['pressure_pa']:.4e} Pa), not exactly zero.

Full metrics and explicit acceptance limits are in `validation.json`. Additional checks
cover target mass, monotonic profiles, mass quadrature, positive stable-branch slope,
and causal sound speed in the density range sampled by this star. A 12,000-point
enthalpy scan checks EOS stability between radial samples, and doubling the dense
EOS sampling checks that resampling errors are negligible. Passing these does
not validate APR at all densities or determine systematic physical uncertainties.

## Scope and sources

- Cold, isotropic, nonrotating spherical hydrostatic equilibrium; Hartle first-order
  rotation computes $I$ without altering the TOV density or pressure profile.
- A barotrope alone does not identify nuclear species, pasta, superfluid regions or
  the chemical crust-core instability. The LAL two-column APR table includes a low-density
  extension, but carries no composition/phase flags; this project makes no unverified
  claim about its precise crust splice.
- The figure's slow rotation and $\epsilon=5\times10^{{-5}}$ motivate a spherical background; magnetic
  stresses, finite temperature, accretion, tides and deformation corrections are omitted.
- The pressure plot shows $P\geq10^{{23}}\,\mathrm{{Pa}}$ and the density plot $\rho\geq10^7$ $\mathrm{{g\,cm^{{-3}}}}$ for
  readability; the CSV includes the full profile to the numerical surface.
- [Probing paper](https://arxiv.org/abs/2511.08895v2), Fig. 2 and Sec. IV.
- [Akmal, Pandharipande & Ravenhall (1998)](https://arxiv.org/abs/nucl-th/9804027).
- [Chamel & Haensel (2008)](https://arxiv.org/abs/0812.3955), crust structure and neutron drip.
- [LAL EOS documentation](https://lscsoft.docs.ligo.org/lalsuite/dev/lalsimulation/group___l_a_l_sim_neutron_star_e_o_s__c.html).
- [Hartle (1967)](https://doi.org/10.1086/149400), relativistic slow rotation.
- [Hinderer (2008), including corrected expressions](https://arxiv.org/abs/0711.2420), tidal Love numbers.

![Calculated structure](structure_profiles.png)
"""


def main():
    """Run the pipeline; diagnostic exports precede the final validation gate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config.json")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    if not 0 < cfg["neutron_drip_rest_mass_density_g_cm3"] < cfg["core_crust_rest_mass_density_g_cm3"]:
        raise ValueError("Require positive neutron-drip density below the core-crust density.")
    eos, provenance = prepare_eos(cfg["eos"])
    star = solve_target(eos, cfg["target_mass_msun"], surface_h=cfg["surface_pseudo_enthalpy"],
                        rtol=cfg["relative_tolerance"], max_step_km=cfg["max_step_km"])
    summary = make_summary(star, cfg, provenance)
    profile = star.profile(cfg["profile_points"], [summary[x]["radius_km"] for x in ["core_crust_boundary", "neutron_drip_boundary"]])
    print(f"Solved: M={star.mass_msun:.9f} Msun, R={star.radius_km:.6f} km", flush=True)
    validation = validate(star, cfg, profile)
    sensitivity = []
    for threshold in cfg["core_crust_sensitivity_g_cm3"]:
        b = star.boundary(threshold)
        sensitivity.append({"threshold_g_cm3": threshold, "core_radius_km": b["radius_km"],
                            "total_crust_thickness_km": star.radius_km-b["radius_km"]})
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    summary["runtime"] = {"python": sys.version.split()[0], **{x: importlib.metadata.version(x) for x in ["numpy", "scipy", "matplotlib", "lalsuite"]}}
    for filename, obj in [("summary.json", summary), ("validation.json", validation), ("config_used.json", cfg)]:
        (out / filename).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    write_csv(out / "profiles.csv", profile)
    write_csv(out / "layers.csv", {k: [x[k] for x in summary["layers"]] for k in summary["layers"][0]})
    write_csv(out / "boundary_sensitivity.csv", {k: [x[k] for x in sensitivity] for k in sensitivity[0]})
    # On failure retain machine-readable diagnostics. Do not generate a fresh
    # report/figure, but note that an older report at this path may still exist.
    if not validation["passed"]:
        raise RuntimeError(f"Numerical validation failed; inspect {out / 'validation.json'}")
    plot_results(profile, summary, out)
    (out / "REPORT.md").write_text(report(summary, validation, sensitivity))
    print(json.dumps({"layers": summary["layers"], "central_pressure_pa": summary["central_pressure_pa"],
                      "validation_passed": validation["passed"], "results": str(out)}, indent=2))


if __name__ == "__main__":
    main()
