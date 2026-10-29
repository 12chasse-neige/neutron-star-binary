# BSk24 composition example

[`plot_composition.py`](plot_composition.py) solves a separate
$1.4\,M_\odot$ BSk24 star and maps its own composition onto radius. It does not
attach BSk24 abundances to the APR background. Run it from the project root:

```sh
.venv/bin/python plot_composition.py
```

The archive path, target mass and output path are fixed in the script. Outputs
are PNG, SVG, PDF, a radial CSV and `summary.json` in `results/composition_bsk24/`.
See the [saved composition guide](results/composition_bsk24/README.md) for the
existing numerical results.

## Archive and thermodynamics

`load_table` reads `data/PCP_BSK24_nparam_1.zip` and verifies its published MD5.
It loads `eos.nb`, `eos.thermo`, `eos.compo` and the supplied `eos.mr` sequence.
Row indices are checked before thermodynamic and composition arrays are combined.
With baryon number density $n_b$ in $\mathrm{fm^{-3}}$, the archive conventions give

$$
P=n_bq_1,\qquad
\mathcal E=n_bm_nc^2(1+q_7).
$$

The numerical neutron rest energy comes from the `eos.thermo` header. In the
NumPy table, $q_1$ and $q_7$ are columns 3 and 9 (zero-based). Both $P$ and
$\mathcal E$ are in $\mathrm{MeV\,fm^{-3}}$, and

$$
1\,\mathrm{MeV\,fm^{-3}}=1.602176634\times10^{32}\,\mathrm{Pa}.
$$

PCHIP interpolates $\ln\mathcal E$ and $\ln n_b$ against $\ln P$, and $\ln P$
against $\ln n_b$. Extrapolation is disabled.

## Independent TOV solution

In $G=c=1$ and kilometre units, `make_solver` integrates $[m,\ln P]$:

$$
\frac{dm}{dr}=4\pi r^2\mathcal E,\qquad
\frac{d\ln P}{dr}=-\frac{\mathcal E+P}{P}
\frac{m+4\pi r^3P}{r(r-2m)}.
$$

The regular-centre expansion starts at $r_0=10^{-4}\,\mathrm{km}$. Brent's
method shoots on central $n_b$ for the target mass. The surface event is
$P(R)=P_{\min}$, the lowest tabulated pressure; there is no extrapolated
atmosphere. Unlike the APR solver, this calculation does not integrate Hartle
rotation or tidal perturbations.

## Baryon fractions

All plotted percentages are **local baryon number fractions**, not mass,
inertia or volume fractions. Charge neutrality supplies total protons:

$$
\begin{aligned}
Y_p&=Y_e+Y_\mu,&Y_n&=1-Y_p,\\
Y_{\rm bound}&=1-Y_{n,\rm free}-Y_{p,\rm free},\\
Y_{p,\rm bound}&=Y_p-Y_{p,\rm free},&
Y_{n,\rm bound}&=Y_{\rm bound}-Y_{p,\rm bound}.
\end{aligned}
$$

The upper plots partition baryons into free neutrons, bound nucleons and free
protons. The lower plots partition the same baryons into all neutrons and all
protons. Bound neutrons count in both the lattice grouping and neutron identity;
these are different partitions. Electrons and muons have zero baryon number and
do not enter the denominator. The archive's nuclear abundance field is not used
as a per-baryon cluster fraction.

The label “superfluid” assumes every unbound neutron is paired:
$Y_{n,\rm paired}=Y_{n,\rm free}$. Actual superfluid density requires pairing gaps,
temperature and entrainment. This archive alone does not determine it or proton
superconductivity. The homogeneous core has no nuclear lattice.

## Phase handling and checks

Archive phase codes are 0 (core), 2 (inner crust) and 1 (outer crust).
`radial_fractions` interpolates linearly in $\ln n_b$ **within each phase**.
A transition between adjacent density samples is placed at

$$
n_{\rm tr}=\sqrt{n_-n_+},\qquad P(r_{\rm tr})=P(n_{\rm tr}).
$$

The last same-phase abundance is held to this midpoint rather than smoothed
across a phase change. Both density and radius brackets are saved; the midpoint
is a plotting prescription, not a resolved coexistence calculation. The radial
grid includes points on both sides of each midpoint so the steps remain visible.

Checks cover fraction closure, nonnegative bound neutron/proton fractions,
target mass and radius convergence under tighter TOV settings. The supplied
mass–radius sequence is also compared with the independently integrated star.
Agreement under refinement does not remove interpolation or surface-prescription
differences from that sequence.

Sources: [CompOSE PCP(BSK24)](https://compose.obspm.fr/eos/253) and the frozen
[Zenodo version](https://zenodo.org/records/14509138). The archive contains its
source README and full model references, including Pearson et al. (2018, 2019).
