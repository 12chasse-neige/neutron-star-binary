# Neutron-star structure and two-component inspiral

Numerical core for two related calculations:

1. **Neutron-star background** — cold spherical TOV equilibrium for APR, Hartle
   slow rotation (moment of inertia) and linear static quadrupolar tides
   ($k_2$, $\Lambda$). Density thresholds classify the radial layers.
2. **Two-component inspiral** — a rigid crust/charged rotator coupled to a
   canonical superfluid reservoir through straight-vortex line integrals,
   entrainment and vortex drag, driven by a circular quadrupole resonance and
   leading gravitational-wave reaction.

The repository keeps only models and computation code. It contains no tests,
plotting, report generation or validation runs.

## Layout

```
numerical-simulation/
├── config/
│   ├── structure.json        APR background inputs and solver settings
│   └── inspiral.json         inspiral scenario and microphysics inputs
├── data/
│   ├── LALSimNeutronStarEOS_APR.dat   frozen two-column APR table
│   ├── APR_monotone_dense.dat         monotone resampling used by the solver
│   └── README.md                      units and provenance
├── docs/
│   ├── structure.md          structure equations, state and exports
│   └── inspiral.md           inspiral equations, reservoir and state
├── nsbinary/
│   ├── units.py              constants and SI/geometrized conversions
│   ├── eos.py                LAL barotrope adapter, dense-table resampling
│   ├── structure.py          Star, integrate, solve_target, build_summary
│   ├── inspiral.py           Reservoir, build_reservoir, Binary, summarize_trace
│   └── io.py                 CSV/JSON writers
├── scripts/
│   ├── compute_structure.py  solve the background, write numeric results
│   └── compute_inspiral.py   evolve scenarios, write traces and summaries
└── results/                  generated numeric outputs only
```

## Install and run

Python 3.11 or newer:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt

.venv/bin/python scripts/compute_structure.py
.venv/bin/python scripts/compute_inspiral.py
```

The inspiral reads the saved background (`results/summary.json` and
`results/profiles.csv`), so run the structure script first. Both scripts accept
`--config`; the structure script also accepts `--data` and `--output`, and the
inspiral script accepts `--background` and `--output`. Defaults resolve from the
repository root. Once dependencies are installed, all runs use saved local data
without network access.

## Outputs

| Path | Contents |
|---|---|
| `results/summary.json` | APR mass, radius, central state, $I$, $k_2$, $\Lambda$, layer boundaries |
| `results/profiles.csv` | Radial profiles: mass, pressure, rest and energy density, $dm/dr$, proper radius |
| `results/layers.csv` | Shell boundaries, thicknesses and gravitational-mass contributions |
| `results/inspiral/vortex_line_integrals.csv` | Cylindrical line columns, mobilities and inertia weights |
| `results/inspiral/<scenario>.csv` | Trace of one scenario: frequencies, lag, heat, work, GW phase |
| `results/inspiral/drag_sweep.csv` | Final heat and budget ratios over the core-drag sweep |
| `results/inspiral/summary.json` | Reservoir metadata, scenario summaries and sweep table |

## Equations at a glance

The structure solver uses $G=c=1$, lengths in kilometres, total energy density
$\mathcal E$ and pressure $P$ in $\mathrm{km^{-2}}$, and dimensionless
pseudo-enthalpy $dh=dP/(\mathcal E+P)$:

$$
\frac{dm}{dr}=4\pi r^2\mathcal E,\qquad
\frac{dh}{dr}=-\frac{m+4\pi r^3P}{r(r-2m)}.
$$

Hartle slow rotation gives $I=J/\Omega$; the static quadrupolar response gives
$\Lambda=2k_2/(3C_\star^5)$ with $C_\star=GM/(Rc^2)$.
[`docs/structure.md`](docs/structure.md) defines the full state, central
expansion, surface matching and units.

The inspiral uses SI units with charged-matter rotation $\Omega_p$ and canonical
superfluid rotation $\Omega_v$. For $\omega=\Omega_v-\Omega_p$,

$$
\dot\Omega_v=-2\mathcal B_{\mathrm{eff}}\Omega_v\omega,\qquad
I_p\dot\Omega_p=N-I_v\dot\Omega_v,\qquad
P_{\mathrm{heat}}=2I_v\mathcal B_{\mathrm{eff}}\Omega_v\omega^2.
$$

The quadrupole is $\Delta I=\epsilon I$ with the **total** $I$, although its
orientation follows the charged component. Orbital reaction and heating share one
energy budget; heat is not added as another orbital torque.
[`docs/inspiral.md`](docs/inspiral.md) derives the line integrals, rigid
projection, circular Hamiltonian and phase diagnostic.

## Scope

- The background is cold, isotropic and spherical; rotation and tides are
  first-order perturbations that do not change its density profile.
- Core and neutron-drip densities are configurable **labels**, not computed
  nuclear phase transitions.
- The inspiral is a rigid two-component reduction. It does not resolve the full
  differential-rotation PDE, vortex pile-up, finite-strength pinning, thermal
  creep, pairing gaps or temperature evolution.
- APR supplies neither the superfluid fractions/entrainment/drag nor the
  prescribed ellipticity; those are assumed inputs in `config/inspiral.json`.

Data provenance, units and the monotone resampling are described in
[`data/README.md`](data/README.md).
