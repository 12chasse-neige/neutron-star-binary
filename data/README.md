# EoS inputs and provenance

The APR and BSk24 datasets serve separate calculations. Neither supplies the
entrainment, drag, pairing or pinning assumptions used in the inspiral model.

| File | Contents | Consumer |
|---|---|---|
| `LALSimNeutronStarEOS_APR.dat` | Frozen two-column APR pressure–energy table | `run.prepare_eos` |
| `LALSimNeutronStarEOS_APR.provenance.json` | Upstream version, units and original SHA-256 | Original-table integrity check |
| `APR_monotone_dense.dat` | Deterministic monotone resampling of APR | `EquationOfState` / LAL |
| `APR_monotone_dense.provenance.json` | Resampling method, row counts and both hashes | Reproducibility record |
| `PCP_BSK24_nparam_1.zip` | Composition-resolved BSk24 archive with source README | `plot_composition.py` |

## APR units and reconstruction

The original APR table was copied without modification from the `lalapps/data`
directory of `lalsuite==7.26.16`. Its columns are $P_{\rm geom}$ and
$\mathcal E_{\rm geom}$ in $\mathrm{m^{-2}}$, with $G=c=1$.
The second column is total energy density including rest energy. In SI,

$$
P_{\rm SI}=\frac{c^4}{G}P_{\rm geom},\qquad
\mathcal E_{\rm SI}=\frac{c^4}{G}\mathcal E_{\rm geom}.
$$

These metre-based table units differ from the solver's $\mathrm{km^{-2}}$.
LAL reconstructs dimensionless pseudo-enthalpy and baryonic rest density from the
barotrope; $\rho_b$ must not be identified with $\mathcal E/c^2$.
The runner compares the original file with its saved SHA-256 before solving.

The native coarse cubic interpolation has a small negative-compressibility
artifact near the inner crust. `run.monotone_table` instead constructs

$$
\ln\mathcal E=\operatorname{PCHIP}(\ln P),
$$

combining 24,000 logarithmically spaced pressures with every original pressure.
Original energy values are restored exactly at those nodes. LAL reads the dense
table and reconstructs its enthalpy and rest density. The original and derived
files retain separate checksums; each APR run rebuilds the derived file.
Validation also repeats the calculation with 48,000 added sampling points.

The table contains a low-density extension but no nuclear species or phase
labels. Its precise microscopic crust splice is not established here. Configured
core and neutron-drip density thresholds classify the numerical star; they are
not information extracted from these two columns. This is APR, not the APR4 fit.

## BSk24 archive

`plot_composition.py` verifies the published archive MD5
`57a9af0939caa9c3621baa28471b2dfd`. It uses `eos.nb`, `eos.thermo` and
`eos.compo` for thermodynamics and abundances, and `eos.mr` for a comparison
mass–radius sequence. Its reconstruction is

$$
P=n_bq_1,\qquad \mathcal E=n_bm_nc^2(1+q_7).
$$

The resulting star uses BSk24 for both structure and composition. See
[COMPOSITION.md](../COMPOSITION.md) for units, phase interpolation and the
explicit fully-paired-neutron illustration.

## Sources and attribution

- [LALSuite source and licensing](https://git.ligo.org/lscsoft/lalsuite).
- [LAL tabulated-EoS implementation](https://lscsoft.docs.ligo.org/lalsuite/dev/lalsimulation/_l_a_l_sim_neutron_star_e_o_s_tabular_8c_source.html).
- [Akmal, Pandharipande & Ravenhall (1998)](https://arxiv.org/abs/nucl-th/9804027).
- [CompOSE PCP(BSK24)](https://compose.obspm.fr/eos/253) and
  [frozen Zenodo archive](https://zenodo.org/records/14509138).

Upstream data retain their authorship and licensing. This copy does not assign
a new license to them. The BSk24 archive preserves its own README and references.
