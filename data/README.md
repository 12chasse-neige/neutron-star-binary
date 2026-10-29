# EoS inputs and provenance

The APR barotrope supplies the neutron-star background; it does not supply the
entrainment, drag, pairing or pinning assumptions used in the inspiral model.

| File | Contents | Consumer |
|---|---|---|
| `LALSimNeutronStarEOS_APR.dat` | Frozen two-column APR pressure–energy table | `nsbinary.eos.load_apr` |
| `APR_monotone_dense.dat` | Deterministic monotone resampling of APR | `EquationOfState` / LAL |

Both files are local and unchanged at run time. If the dense table is absent it
is rebuilt deterministically from the original table.

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

The native coarse cubic interpolation has a small negative-compressibility
artifact near the inner crust. `nsbinary.eos.monotone_table` instead constructs

$$
\ln\mathcal E=\operatorname{PCHIP}(\ln P),
$$

combining 24,000 logarithmically spaced pressures with every original pressure.
Original energy values are restored exactly at those nodes. LAL reads the dense
table and reconstructs its enthalpy and rest density. The original file is kept
as the frozen upstream data; the dense file is the resampling consumed by the
solver.

The table contains a low-density extension but no nuclear species or phase
labels. Its precise microscopic crust splice is not established here. Configured
core and neutron-drip density thresholds classify the numerical star; they are
not information extracted from these two columns. This is APR, not the APR4 fit.

## Sources and attribution

- [LALSuite source and licensing](https://git.ligo.org/lscsoft/lalsuite).
- [LAL tabulated-EoS implementation](https://lscsoft.docs.ligo.org/lalsuite/dev/lalsimulation/_l_a_l_sim_neutron_star_e_o_s_tabular_8c_source.html).
- [Akmal, Pandharipande & Ravenhall (1998)](https://arxiv.org/abs/nucl-th/9804027).
- [Chamel & Haensel (2008)](https://arxiv.org/abs/0812.3955), crust structure and neutron drip.

Upstream data retain their authorship and licensing. This copy does not assign
a new license to them.
