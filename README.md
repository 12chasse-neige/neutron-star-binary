# Neutron-star structure and two-component inspiral

This project connects an APR stellar background to an illustrative spin–orbit
resonance model. A separate BSk24 calculation supplies composition plots.
The three calculations have distinct inputs and assumptions.

| Calculation | Entry point | Equations | Saved results |
|---|---|---|---|
| APR equilibrium, inertia and tides | [`run.py`](run.py) | [Structure model](MODEL.md) | [Structure report](results/REPORT.md) |
| Two-component inspiral and vortex heating | [`inspiral/run.py`](inspiral/run.py) | [Inspiral model](inspiral/MODEL.md) | [Inspiral report](results/inspiral/REPORT.md) |
| BSk24 composition and its own TOV solution | [`plot_composition.py`](plot_composition.py) | [Composition model](COMPOSITION.md) | [Composition guide](results/composition_bsk24/README.md) |

The resonance calculation extends the aligned circular mechanism in
[Miao & Yang](https://arxiv.org/abs/2511.08895v2) with a rigid projection of the
entrained-superfluid model of
[Antonelli & Pizzochero](https://arxiv.org/abs/1603.02838v4).
APR does not supply the assumed superfluid microphysics or predict the
prescribed ellipticity.

## Run

Use Python 3.11 or newer from this directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run.py
.venv/bin/python plot_composition.py
.venv/bin/python inspiral/run.py
```

The inspiral reads `results/summary.json` and `results/profiles.csv`, so generate
the APR background first. The BSk24 calculation is independent. Once dependencies
are installed, all three use saved local data without network access.

Both configurable runners accept `--config` and `--output`:

```sh
.venv/bin/python run.py --config config.json --output results
.venv/bin/python inspiral/run.py --config inspiral/config.json --output results/inspiral
```

Defaults resolve relative to each script; explicit relative paths resolve from
the working directory. The inspiral `--output` changes its destination, **not**
its background source, which remains the root `results/` directory. Composition
paths and target mass are fixed in `plot_composition.py`.

Runs overwrite matching output files. APR also rebuilds its dense EoS table and
provenance in `data/`, even with a different output directory. Runners save
numerical diagnostics before checking final validation, and write a new report
and plots only after validation passes. Read any retained report alongside the
current `validation.json` and `config_used.json`.

Run the existing physics tests separately:

```sh
.venv/bin/python -m unittest -v
.venv/bin/python -m unittest discover -s inspiral -p 'test_inspiral.py' -v
```

The runners add convergence and conservation checks beyond these unit tests.

## Physical inputs

The APR configuration targets $M=1.4\,M_\odot$. Its `binary_context` records

$$
q=\frac{M'}M=3,\qquad \epsilon=5\times10^{-5},\qquad
f_{\mathrm{orb},0}=0.4\,\mathrm{Hz},\qquad
f_{\mathrm{spin},0}=0.5\,\mathrm{Hz},\qquad f_0-\theta_0=0.3986.
$$

Only the EoS and target mass determine the spherical background. Binary context
does not deform its density profile. `inspiral/config.json` separately sets the
parameters actually evolved; edits to root `binary_context` do not propagate.

APR layers are a rest-mass-density classification:

| Region | Default condition on $\rho_b$ |
|---|---|
| Core | $\rho_b\geq1.4\times10^{14}\,\mathrm{g\,cm^{-3}}$ |
| Inner crust | $4\times10^{11}\leq\rho_b<1.4\times10^{14}\,\mathrm{g\,cm^{-3}}$ |
| Outer crust | $\rho_b<4\times10^{11}\,\mathrm{g\,cm^{-3}}$ to the numerical surface |

These thresholds are assumptions, not calculated nuclear transitions.
`core_crust_sensitivity_g_cm3` reclassifies the same star without changing its
EoS. BSk24 instead uses tabulated phase labels and reports transition brackets.

## Equations at a glance

The structure solver uses $G=c=1$, lengths in kilometres, total energy density
$\mathcal E$ and pressure $P$ in $\mathrm{km^{-2}}$, and dimensionless
pseudo-enthalpy $dh=dP/(\mathcal E+P)$:

$$
\frac{dm}{dr}=4\pi r^2\mathcal E,\qquad
\frac{dh}{dr}=-\frac{m+4\pi r^3P}{r(r-2m)}.
$$

Hartle slow rotation gives $I=J/\Omega$. The static quadrupolar response gives
$\Lambda=2k_2/(3C_\star^5)$, where $C_\star=GM/(Rc^2)$.
[MODEL.md](MODEL.md) defines the full state vector, central expansion, surface
matching and units.

The inspiral uses SI units, with charged-matter rotation $\Omega_p$ and canonical
superfluid rotation $\Omega_v$. For $\omega=\Omega_v-\Omega_p$,

$$
\dot\Omega_v=-2\mathcal B_{\mathrm{eff}}\Omega_v\omega,\qquad
I_p\dot\Omega_p=N-I_v\dot\Omega_v,\qquad
P_{\mathrm{heat}}=2I_v\mathcal B_{\mathrm{eff}}\Omega_v\omega^2.
$$

The imposed quadrupole is $\Delta I=\epsilon I$ with **total** $I$, although its
orientation follows the charged component. Orbital reaction and heating share
one energy budget; heat is not added as another orbital torque.

## Suggested review order

1. Read [data/README.md](data/README.md) for data contents and provenance.
2. Read [MODEL.md](MODEL.md), then `EquationOfState`, `integrate` and
   `solve_target` in [`ns_model.py`](ns_model.py).
3. Review `prepare_eos`, `make_summary` and `validate` in [`run.py`](run.py):
   interpolation, layer classification and numerical evidence are separate steps.
4. Read [inspiral/MODEL.md](inspiral/MODEL.md), then `build_reservoir`,
   `Binary.quantities` and `Binary.evolve` in [`inspiral/model.py`](inspiral/model.py).
5. Review scenario choices and acceptance limits in
   [`inspiral/run.py`](inspiral/run.py). Paper-start and initially locked cases
   address capture and subsequent locking separately.
6. Review [COMPOSITION.md](COMPOSITION.md) and
   [`plot_composition.py`](plot_composition.py) independently of APR.

The scripts have detailed comments at state definitions, physical conversions,
model reductions and numerical checks.

## Outputs and scope

| Location | Contents |
|---|---|
| `data/` | Frozen APR table, derived dense table, provenance, BSk24 archive |
| `results/profiles.csv` | Areal radius, gravitational mass, pressure, rest/energy-equivalent densities, proper radial distance |
| `results/layers.csv` | Boundaries, coordinate/proper thicknesses, gravitational mass contributions |
| `results/boundary_sensitivity.csv` | Alternative core-threshold classifications |
| `results/summary.json`, `validation.json` | APR observables, provenance, versions, measured errors and limits |
| `results/structure_profiles.png`, `.svg` | APR radial structure |
| `results/composition_bsk24/` | Independent composition figures, radial CSV and summary |
| `results/inspiral/` | Scenario traces, line integrals, drag sweep, summaries, balances and plots |

Each configurable runner saves `config_used.json`. Inspiral also hashes its
background files. Numerical convergence tests the chosen equations; it does not
establish pairing gaps, entrainment, pinning, magnetic deformation or the physical
accuracy of the rigid two-component reduction. Sources and further limits appear
in the model documents and generated reports.
