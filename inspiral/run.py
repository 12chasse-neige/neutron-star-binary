#!/usr/bin/env python3
"""Run the two-component inspiral scenarios, drag sweep and numerical audits.

MODEL.md maps the equations to model.py. This runner always consumes the root
results/{summary.json,profiles.csv}; --output only selects its own destination.
The saved background is read, not regenerated or mixed with BSk24 composition.
"""
import argparse
import csv
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT.parent/".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

from model import Binary, build_reservoir, G, c


def save_json(path, obj):
    """Write reviewable metadata; convert NumPy scalars and reject NaN/Inf."""
    def convert(value):
        return value.item() if isinstance(value, np.generic) else str(value)
    path.write_text(json.dumps(obj, indent=2, default=convert, allow_nan=False)+"\n")


def save_csv(path, columns, max_rows=None):
    """Export aligned columns, optionally thinning by index after all audits.

    The retained rows span the full trace, but their spacing is not generally
    uniform in time or frequency because adaptive solver nodes are included.
    """
    count = len(next(iter(columns.values())))
    rows = np.unique(np.linspace(0, count-1, min(count, max_rows), dtype=int)) if max_rows else range(count)
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*(np.asarray(v)[rows] for v in columns.values())))


def plots(output, traces, sweep, reservoir):
    """Draw dynamics and budgets; plotting floors do not alter saved data."""
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 180})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.8), constrained_layout=True)
    data = traces["fiducial_paper"]
    x = 2*data["kepler_orbital_frequency_hz"]
    ax = axes[0,0]
    ax.plot(x, x/2, "k--", lw=1.3, label="Orbit")
    ax.plot(x, data["crust_frequency_hz"], color="#bd562c", label="Crust / charged matter")
    ax.plot(x, data["canonical_superfluid_frequency_hz"], color="#276a9e", label="Superfluid canonical rotation")
    ref = traces["rigid_locked"]
    ax.plot(2*ref["kepler_orbital_frequency_hz"], ref["crust_frequency_hz"], color="#777777", lw=1.1, label="Rigid star, initially locked")
    ax.set(xlabel="GW frequency (Hz)", ylabel="Rotation frequency (Hz)",
           title="A  Crust locking can outlast whole-star locking", xlim=(.8,50))
    ax.legend(fontsize=8, loc="upper left")
    ax = axes[0,1]
    ax.plot(x, data["canonical_lag_rad_s"], color="#276a9e", label="Canonical lag")
    ax.axhline(0, color=".6", lw=.6)
    ax.set(xlabel="GW frequency (Hz)", ylabel=r"$\Omega_v-\Omega_p$ (rad s$^{-1}$)",
           title="B  Spin-up drives vortices inward", xlim=(.8,50))
    twin=ax.twinx()
    twin.plot(x, data["signed_rms_vortex_radial_velocity_m_s"], color="#818a36", alpha=.7)
    twin.set_ylabel("Signed rms radial vortex speed (m/s)", color="#6c752c")
    ax = axes[1,0]
    for key, label, color in [("orbital_binding_energy_loss_J", "Orbital binding energy released", "#333333"),
                              ("star_torque_work_J", "Net work delivered to star", "#bd562c"),
                              ("heat_J", "Vortex heat", "#276a9e")]:
        ax.semilogy(x, np.maximum(data[key], 1e25), label=label, color=color)
    ax.set(xlabel="GW frequency (Hz)", ylabel="Cumulative energy (J)",
           title="C  Heat is small compared with orbital loss", ylim=(1e34,1e46), xlim=(.8,50))
    ax.legend(fontsize=8)
    ax = axes[1,1]
    drag=np.array([s["core_drag_ratio"] for s in sweep])
    frac=np.array([s["heat_over_orbital_binding_energy_loss"] for s in sweep])
    ax.loglog(drag, frac, "o-", color="#276a9e", label="Two-layer integrations")
    ax.axvline(1e-4, ls=":", color="#bd562c", label="Fiducial drag")
    ax.set(xlabel=r"Core drag ratio $\eta/(\kappa\rho_n)$", ylabel="Heat / orbital binding energy released",
           title="D  Drag sweep, each started in resonance")
    ax.legend(fontsize=8)
    fig.suptitle("Two-layer inspiral with entrainment and vortex dissipation\n"
                 r"$M=1.4\,M_\odot$, $M'=4.2\,M_\odot$, $\epsilon=5\times10^{-5}$; illustrative microphysics", fontsize=14)
    for extension in ["png", "svg"]:
        fig.savefig(output/f"inspiral_dissipation.{extension}")
    plt.close(fig)
    fig, axes = plt.subplots(1,3,figsize=(13,3.8),constrained_layout=True)
    t=data["time_s"]/86400
    axes[0].plot(t,data["orbital_frequency_hz"],"k--", label="Orbit")
    axes[0].plot(t,data["crust_frequency_hz"],color="#bd562c",label="Crust")
    axes[0].plot(t,data["canonical_superfluid_frequency_hz"],color="#276a9e",label="Canonical superfluid")
    axes[0].set(xlabel="Time from paper initial conditions (days)",ylabel="Frequency (Hz)", title="Inspiral and capture")
    axes[0].legend(fontsize=8)
    axes[1].semilogy(x, np.maximum(data["instantaneous_heat_over_gw_power"], 1e-12),color="#276a9e")
    axes[1].set(xlabel="GW frequency (Hz)",ylabel="Heating power / GW luminosity",title="Instantaneous dissipation",ylim=(1e-10,1e-3))
    xx=reservoir.columns["cylindrical_radius_m"]/1000
    axes[2].semilogy(xx,reservoir.columns["B_radial_free"],color="#276a9e")
    axes[2].axvline(reservoir.metadata["core_radius_m"]/1000,color=".5",ls=":")
    axes[2].set(xlabel="Cylindrical radius (km)",ylabel=r"$\mathcal{B}^x(x)$",title="Straight-vortex line average")
    fig.savefig(output/"capture_and_vortices.png")
    plt.close(fig)


def make_report(output, cfg, reservoir, summaries, sweep, checks, bounds):
    """Explain saved scenario results and model assumptions in Markdown."""
    s = summaries["fiducial_paper"]
    rigid = summaries["rigid_locked"]
    largest=max(sweep,key=lambda row:row["heat_J"])
    rows=[]
    for row in sweep:
        slip = f'{row["first_phase_slip_gw_hz"]:.2f}' if row["first_phase_slip_gw_hz"] else ">50 or no capture"
        rows.append(f'| {row["core_drag_ratio"]:.0e} | {row["mean_B"]:.3e} | {row["heat_J"]:.3e} | '
                    f'{100*row["heat_over_orbital_binding_energy_loss"]:.6g}% | {slip} |')
    text=rf"""# Two-component inspiral, spin-orbit locking and vortex heating

The [model guide](../../inspiral/MODEL.md) defines the two canonical rotators,
line integrals and circular Hamiltonian. The heat and spin budgets obey

$$
\omega=\Omega_v-\Omega_p,\qquad
\dot\Omega_v=-2\mathcal B_{{\rm eff}}\Omega_v\omega,\qquad
P_{{\rm heat}}=2I_v\mathcal B_{{\rm eff}}\Omega_v\omega^2,
$$

$$
\Delta E_{{\rm spin}}+Q=W_\star,\qquad
\Delta(E_K+U+E_{{\rm spin}})+Q+E_{{\rm GW}}=0.
$$

## Result

For the specified illustrative microphysics, **vortex dissipation does not remove
a large fraction of the orbital energy**. Starting with the probing paper's
Fig. 2 inputs and evolving from GW frequency 0.8 to 50 Hz gives:

| Quantity | Two-layer result |
|---|---:|
| Duration | {s['duration_s']/86400:.6f} days |
| Orbital binding energy released | {s['orbital_binding_energy_loss_J']:.6e} J |
| Net quadrupole work delivered to star | {s['star_torque_work_J']:.6e} J |
| Increase in stellar rotational energy | {s['spin_energy_change_J']:.6e} J |
| Vortex heat | {s['heat_J']:.6e} J = {s['heat_J']*1e7:.6e} erg |
| Heat / orbital binding energy released | {s['heat_over_orbital_binding_energy_loss']:.6e} ({100*s['heat_over_orbital_binding_energy_loss']:.6f}%) |
| Heat / net work delivered to star | {100*s['heat_fraction_of_net_star_work']:.3f}% |
| Maximum instantaneous heat / GW luminosity | {s['max_instantaneous_heat_over_gw_power']:.6e} |
| Captured in the 0.8-1.0 Hz orbital window | {s['captured_in_0p8_to_1_hz_window']} |
| First full phase slip after capture | {s['first_phase_slip_gw_hz']:.3f} Hz GW |
| Final crust / canonical superfluid frequency | {s['final_crust_frequency_hz']:.5f} / {s['final_superfluid_frequency_hz']:.5f} Hz |

Positive net stellar work equals rotational-energy increase plus heat, within the
recorded numerical error. Heat is not an additional orbital torque. It can be a
substantial part of the much smaller energy budget delivered to the star while
remaining tiny relative to the energy lost by the orbit. No temperature increase
is inferred: heat capacity, cooling, pairing and thermal feedback are not evolved.

![Dynamics and energy](inspiral_dissipation.png)

## Meaning of two layers

There are **two effective rotational components**: p is the crust plus tightly
coupled charged matter; v is the canonical superfluid reservoir. The v component
is not the inertia of microscopic vortex cores. Core and inner-crust zones have
separate prescribed microphysics, and the outer crust is entirely in p. Straight
vortices thread both superfluid zones. This implements the rigid/body-averaged
reduction in Antonelli & Pizzochero Sec. 3.1, retaining the nonlinear factor
$2\Omega_v$; it does not solve their differential-rotation PDE.

The saved APR TOV profile is used for density, radius and layer boundaries. Its
Newtonian baryon-density inertia integrals are multiplied by a common factor
{reservoir.metadata['common_Hartle_rescaling']:.6f} to match the saved Hartle total
$I$ = {reservoir.total_I:.6e} $\mathrm{{kg\,m^2}}$. The effective fractions are
$I_v/I$ = {reservoir.fv:.6f} and $I_p/I$ = {reservoir.fp:.6f}. This is a Newtonian
two-fluid dynamics model calibrated to the relativistic total inertia, rather than
a relativistic entrained stellar perturbation calculation. The outer rotator is
not a geometrically isolated shell: it includes charged matter throughout the core.

## Physical inputs and controlled comparisons

The stellar mass is $M=1.4\,M_\odot$, the companion-to-star mass ratio is $q=3$,
and the fixed total-star ellipticity is $\epsilon=5\times10^{{-5}}$.
The quadrupole difference is $\Delta I=\epsilon I$, carried by the p component;
it is not changed to $\epsilon I_p$.
The original initial frequencies are orbital 0.4 Hz and spin 0.5 Hz, with
$f-\theta=0.3986\,\mathrm{{rad}}$. Both components initially corotate. The fiducial paper-start
run explicitly tests capture; no locking constraint is imposed during evolution.

The rigid control started at the identical original phase does not capture in
the 0.8-1.0 Hz orbital window (capture flag:
{summaries['rigid_paper']['captured_in_0p8_to_1_hz_window']}). The small circular
backreaction and phase accumulation can alter capture relative to the paper's
neglected-backreaction example. The reported 0.3986 phase is rounded and does not
guarantee capture under changed equations or constants. No initial phase is tuned
to force a rigid reproduction of Fig. 2. The separate rigid, initially locked
trajectory is the appropriate control for conditional locking evolution.

| Assumption | Core | Inner crust |
|---|---:|---:|
| Superfluid neutron mass fraction | {cfg['core_neutron_fraction']} | {cfg['inner_crust_neutron_fraction']} |
| Entrainment $\varepsilon_n=1-m_n^*/m_n$ | {cfg['core_entrainment_epsilon_n']} | {cfg['inner_crust_entrainment_epsilon_n']} |
| Drag ratio $\eta/(\kappa\rho_n)$ | {cfg['core_drag_ratio']} | {cfg['inner_crust_drag_ratio']} |
| Freely moving vortex fraction $Y$ | {cfg['unpinned_fraction']} | {cfg['unpinned_fraction']} |

These constants are **assumed sensitivity inputs**, not microphysics extracted
from the APR barotrope or a fit to either paper. No finite critical pinning force
or thermal creep law is specified. The $Y=0$ case is a permanently pinned limiting
control, with zero vortex motion and zero heating; finite-threshold unpinning and
glitches are outside this model. A suppressed-mobility $Y=10^{{-3}}$ case is also included.

The drag sweep starts each trajectory already in the stable resonant island at
orbital frequency {cfg['locked_start_orbital_frequency_hz']} Hz with corotating
components. These are **conditional-on-capture** comparisons, not capture
probabilities. The crust-to-core drag ratio is fixed at
{cfg['drag_sweep_crust_over_core']}. The sweep reaches both sides of the drag law:
more microscopic drag does not always mean more radial vortex motion.

| Core drag ratio | Mean $\mathcal B_{{\rm eff}}$ | Heat (J) | Heat / orbital loss | First full phase slip (GW Hz) |
|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

The largest integrated heat in this sampled sweep is {largest['heat_J']:.4e} J,
or {100*largest['heat_over_orbital_binding_energy_loss']:.5f}% of orbital loss.
This is a maximum over the sampled inputs, not a universal bound.

## Locking, vortex motion and orbital phase

In the spin-up portion of resonance $\Omega_p$ exceeds $\Omega_v$. The negative lag
drives vortices inward, increasing their areal density and spinning up the
superfluid. Entrainment changes the effective inertias and line-averaged mobility;
it produces no heat by itself. Drag produces nonnegative heat for either lag sign.
With the assumed inputs, the crust can remain locked while the superfluid lags,
reducing the instantaneous inertia load and extending crust locking. This effect
is distinct from increasing the star's ellipticity.

The ideal rigid, initially locked control first completes a phase slip at
{rigid['first_phase_slip_gw_hz']:.3f} Hz GW; the analytic quasistatic rigid torque
limit is {bounds['analytic_rigid_breaking_gw_hz']:.3f} Hz GW. The paper quotes
about 14.2 Hz for its finite-amplitude illustrative trajectory. A full phase slip
is an operational endpoint and need not equal the torque-limit frequency or the
first departure from frequency tracking. Initial phases and capture amplitudes
matter. The phase-slip diagnostic here uses a fixed reference at orbital 1 Hz.

The fiducial time-domain GW phase difference from a point particle at matched
frequency is {s['final_gw_phase_difference_from_point_particle_rad']:.3f} rad over
the entire 0.8-50 Hz GW interval. It includes conservative spin-up, changed locking
and dissipation together; it is not an isolated dissipative phase shift, the
frequency-domain phase in the paper's Eq. (6), or a detectability estimate.
Small fractional energy transfer can accumulate an appreciable phase over days.

![Capture and vortex mobility](capture_and_vortices.png)

## Energy-transfer bound at fixed ellipticity

The maximum quadrupole torque power at the initial orbital frequency, divided by
the GW luminosity, is {bounds['max_available_torque_power_over_gw_flux']:.6f}.
It decreases approximately as $n^{{-1/3}}$. Thus this channel cannot remove
order-unity orbital energy for the specified masses, ellipticity and frequency
range, regardless of how efficiently vortices dissipate the delivered work.
The absolute-torque integrated upper bound, allowing either sign of the actual
torque, its feedback on the chirp, and the small interaction-frequency correction, is
{bounds['absolute_torque_work_bound_J']:.4e} J,
{100*bounds['absolute_torque_work_bound_fraction']:.4f}% of orbital binding energy
released. This is a bound on orbital torque work; initial spin and interaction
energy must also be included when bounding total heat. See [MODEL.md](../../inspiral/MODEL.md) for the derivation.
This does not bound other orbital dissipation channels such as dynamical tides,
resonant stellar modes, magnetic interactions or vortex-flux-tube cutting, which
are absent from the model.

## Validation and reproducibility

Validation status: **{checks['status']}**. The run checks energy and angular
momentum balance, nonnegative heating, convergence under tighter integration
tolerances, and doubled line/cylindrical quadrature. Independent unit tests cover
isolated lag relaxation and its exact heat budget, the closed quadrupole
Hamiltonian, point-particle inspiral, and zero/infinite drag limits.
Detailed values and thresholds are in validation.json. Input checksums, versions,
parameters and all scenario summaries are in summary.json and config_used.json.
Trace CSV files are decimated to at most 15,000 rows for convenient inspection;
the energy diagnostics use the full adaptive solution plus a dense output grid.

Run from the `numerical-simulation` folder with
`.venv/bin/python inspiral/run.py`; use `--config` and `--output` to select other
paths. Rerun tests from that folder with
`.venv/bin/python -m unittest discover -s inspiral -p 'test_inspiral.py' -v`.

## Scope and sources

Large lags and a strongly radius-dependent $\mathcal B^x$ may invalidate the rigid projection;
the full AP PDE can develop differential rotation and vortex pile-up that this
two-component truncation cannot resolve. The nonaxisymmetric quadrupole is added
as a small external perturbation to an axisymmetric spherical background.
Turbulence, line bending, magnetic stress evolution, superfluid phase transitions,
finite-temperature feedback and dynamical tides are not simulated. A fixed
quadrupole in a magnetar does not establish which dissipation mechanisms are
active in its core. Results are conditional model calculations, not a measured or
unique neutron-star microphysical prediction. See [MODEL.md](../../inspiral/MODEL.md)
for equations and conventions.

- Antonelli & Pizzochero, *Axially symmetric equations for differential pulsar
  rotation with superfluid entrainment*, MNRAS 464, 721 (2017),
  [arXiv:1603.02838v4](https://arxiv.org/abs/1603.02838v4),
  [DOI](https://doi.org/10.1093/mnras/stw2376): Eqs. (2)-(7), (11), (13),
  (15), (18)-(22), (26)-(27), and the rigid reduction in Sec. 3.1.
- Miao & Yang, *Probing Intrinsic Ellipticity in Compact Star Binaries*,
  [arXiv:2511.08895v2](https://arxiv.org/abs/2511.08895v2),
  Secs. II-III and Fig. 2. The circular orbital backreaction is implemented with
  a Hamiltonian and radiation flux, as explained in [MODEL.md](../../inspiral/MODEL.md).
- Existing local APR background: ../summary.json and ../profiles.csv. Those
  files and their original equilibrium outputs are unchanged by this simulation.
"""
    (output/"REPORT.md").write_text(text)


def main():
    """Build the reservoir, run controls/sweep, validate and then publish reports."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"config.json")
    parser.add_argument("--output", type=Path, default=ROOT.parent/"results"/"inspiral")
    args=parser.parse_args()
    cfg=json.loads(args.config.read_text())
    output=args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    reservoir=build_reservoir(cfg, ROOT.parent/"results")
    save_json(output/"config_used.json",cfg)
    save_csv(output/"vortex_line_integrals.csv",reservoir.columns)
    summaries,traces={},{}
    def run(name, this_cfg=cfg, rigid=False, init="locked", save=True):
        """Run one scenario; rebuild line integrals when microphysics changes.

        Keep every summary for balance checks. Sweep traces need not be retained
        for plotting, but their full-resolution diagnostics are still evaluated.
        """
        start=time.monotonic()
        res=reservoir if this_cfg is cfg else build_reservoir(this_cfg, ROOT.parent/"results")
        data,summary=Binary(this_cfg,res,rigid=rigid).evolve(init)
        summary["wall_time_s"]=time.monotonic()-start
        summaries[name]=summary
        if save:
            save_csv(output/f"{name}.csv",data,max_rows=15000)
            traces[name]=data
        print(f'{name}: heat={summary["heat_J"]:.5e} J; heat/orbit={summary["heat_over_orbital_binding_energy_loss"]:.5e}; {summary["wall_time_s"]:.1f}s',flush=True)
        return data,summary
    # Same original phase for both paper-start models: capture is an outcome.
    # Locked controls isolate later evolution conditional on being in resonance.
    run("fiducial_paper",init="paper")
    run("rigid_paper",rigid=True,init="paper")
    run("fiducial_locked")
    run("rigid_locked",rigid=True)
    run("permanently_pinned_locked",{**cfg,"unpinned_fraction":0.})
    run("suppressed_mobility_locked",{**cfg,"unpinned_fraction":1e-3})
    # Each sweep member has its own locked initialization, not the endpoint of
    # the previous trajectory. Fix the crust/core drag ratio, not the mean B.
    sweep=[]
    for drag in cfg["drag_sweep_core_ratios"]:
        local={**cfg,"core_drag_ratio":drag,
               "inner_crust_drag_ratio":drag*cfg["drag_sweep_crust_over_core"]}
        _,s=run(f"drag_{drag:g}",local,save=False)
        sweep.append({"core_drag_ratio":drag, **s})
    save_csv(output/"drag_sweep.csv",{key:[row[key] for row in sweep] for key in sweep[0]})
    # Convergence of the heating observable, not a repeated comparison of the
    # same integrator with itself at unchanged settings.
    refined_res=build_reservoir(cfg,ROOT.parent/"results",order_multiplier=2)
    _,refined=Binary(cfg,reservoir).evolve("locked",accuracy=5,sample_count=1200)
    _,fine_paper=Binary(cfg,reservoir).evolve("paper",accuracy=5,sample_count=1200)
    quantities={
        "cylindrical_vs_spherical_Iv_error":reservoir.metadata["cylindrical_vs_spherical_Iv_relative_error"],
        "doubled_quadrature_B_relative_change":abs(refined_res.B/reservoir.B-1),
        "doubled_quadrature_Iv_relative_change":abs(refined_res.fv/reservoir.fv-1),
        "locked_heat_relative_change_tighter_tolerance":abs(refined["heat_J"]/summaries["fiducial_locked"]["heat_J"]-1),
        "paper_heat_relative_change_tighter_tolerance":abs(fine_paper["heat_J"]/summaries["fiducial_paper"]["heat_J"]-1),
        "max_energy_balance_error":max(s["max_energy_balance_error_over_orbital_loss"] for s in summaries.values()),
        "max_angular_momentum_balance_error":max(s["max_angular_momentum_balance_error_over_orbital_loss"] for s in summaries.values()),
        "max_spin_energy_balance_error":max(s["max_spin_balance_error_over_work"] for s in summaries.values()),
        "minimum_heating_power_W":min(s["minimum_mutual_friction_power_W"] for s in summaries.values()),
    }
    limits={"cylindrical_vs_spherical_Iv_error":3e-6,
            "doubled_quadrature_B_relative_change":1e-5,
            "doubled_quadrature_Iv_relative_change":1e-5,
            "locked_heat_relative_change_tighter_tolerance":3e-5,
            "paper_heat_relative_change_tighter_tolerance":3e-5,
            "max_energy_balance_error":1e-8,
            "max_angular_momentum_balance_error":1e-8,
            "max_spin_energy_balance_error":1e-6}
    passed={key:bool(quantities[key]<limit) for key,limit in limits.items()}
    passed["heating_nonnegative"]=quantities["minimum_heating_power_W"]>=-1e-10
    passed["paper_capture_preserved_under_refinement"]=fine_paper["captured_in_0p8_to_1_hz_window"]==summaries["fiducial_paper"]["captured_in_0p8_to_1_hz_window"]
    validation={"status":"PASS" if all(passed.values()) else "FAIL", "values":quantities,"limits":limits,"passed":passed,
                "refined_locked_summary":refined,"refined_paper_summary":fine_paper}
    save_json(output/"validation.json",validation)
    # Bound the available orbital torque work without assuming its sign. The
    # lambdas use I-normalized quantities; multiply the final integral by I.
    # The denominator must remain positive for the bound to apply. The default
    # small quadrupole satisfies this; new physical inputs require rechecking it.
    binary=Binary(cfg,reservoir)
    n0,n1=2*np.pi*np.array([cfg["initial_orbital_frequency_hz"],cfg["final_orbital_frequency_hz"]])
    tau_amplitude=lambda n:1.5*binary.ratio*binary.eps*n*n
    flux=lambda n:binary.flux_coeff*n**(10/3)
    alpha=tau_amplitude(n0)*n0/flux(n0)
    ell=lambda n:binary.ell_coeff*n**(-1/3)
    delta=lambda n:4.5*binary.ratio*binary.eps*n/ell(n)
    bound=quad(lambda n:tau_amplitude(n)*ell(n)*n*(1+delta(n))**2/(3*(flux(n)-n*tau_amplitude(n)*(1+delta(n)))),n0,n1,epsabs=1e-8,epsrel=1e-10)[0]*binary.I
    orbital_loss=.5*(ell(n1)*n1-ell(n0)*n0)*binary.I
    bounds={"max_available_torque_power_over_gw_flux":alpha,
            "maximum_quadrupole_frequency_fractional_correction":delta(n1),
            "absolute_torque_work_bound_J":bound,
            "absolute_torque_work_bound_fraction":bound/orbital_loss,
            "analytic_rigid_breaking_gw_hz":c**3/(np.pi*G*binary.mt)*(5/64*(1+cfg["mass_ratio_companion_over_ns"]))**(.6)*cfg["ellipticity"]**(.6)}
    save_json(output/"summary.json",{"reservoir":reservoir.metadata,"scenarios":summaries,"bounds":bounds,
               "runtime":{"python":platform.python_version(),**{p:importlib.metadata.version(p) for p in ["numpy","scipy","matplotlib"]}}})
    # Keep diagnostics for failed runs, but do not produce a fresh success
    # report/plot. A report from a previous run can still exist at this path.
    if not all(passed.values()):
        raise RuntimeError(f"Validation failed: {passed}")
    plots(output,traces,sweep,reservoir)
    make_report(output,cfg,reservoir,summaries,sweep,validation,bounds)
    print(f"Validated results: {output}",flush=True)


if __name__ == "__main__":
    main()
