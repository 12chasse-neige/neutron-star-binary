# Two-component inspiral: equations and implementation

[`model.py`](model.py) implements a rigid projection of the straight-vortex
entrainment model of [Antonelli & Pizzochero (AP)](https://arxiv.org/abs/1603.02838v4),
combined with the aligned circular quadrupole resonance of
[Miao & Yang](https://arxiv.org/abs/2511.08895v2).
[`run.py`](run.py) chooses scenarios, checks convergence and writes reports.
This is an extension under explicit approximations, not a direct reproduction
of a published two-fluid binary calculation.

## 1. Components, geometry and notation

There are two effective rotational components and two superfluid microphysics
zones. These are different divisions of the star:

| Symbol | Meaning | Code |
|---|---|---|
| $\Omega_p$ | Rigid crust plus charged matter throughout the star | `op` |
| $\Omega_v$ | Spatially uniform canonical superfluid rotation | `ov` |
| $\Omega_n(r)$ | Neutron transport rotation, reconstructed per zone | Transport-frequency CSV columns |
| $\omega=\Omega_v-\Omega_p$ | Canonical lag | `ov-op` |
| $I_v$, $I_p=I-I_v$ | Effective canonical inertias | `fv*I`, `fp*I` |
| $x$, $z$ | Cylindrical radius and axial coordinate | Line quadrature |
| $r_{\rm cc}$, $R_d$ | Core edge and neutron-drip radius | Saved APR boundaries |
| $\epsilon$ | Prescribed equatorial ellipticity | `eps` |
| $\varepsilon_n$ | Neutron entrainment parameter | Zone configuration |
| $\gamma=2(\theta-f)$ | Figure–orbit phase; $f$ is orbital angle | `phase` |
| $n$, $\Omega_{\rm orb}$ | Kepler and Hamiltonian orbital angular frequencies | `n`, `orb` |

All dynamics use SI. Angular frequencies are in $\mathrm{rad\,s^{-1}}$;
frequencies in Hz divide these by $2\pi$. The canonical rotation is neither
neutron transport rotation nor the angular velocity of vortex lines:

$$
\begin{aligned}
\Omega_v&=(1-\varepsilon_n)\Omega_n+\varepsilon_n\Omega_p,\\
\Omega_n&=\frac{\Omega_v-\varepsilon_n\Omega_p}{1-\varepsilon_n},\qquad
\varepsilon_n=1-\frac{m_n^*}{m_n},\\
\rho_n&=x_n\rho_b,\qquad
\rho_n\varepsilon_n=\rho_p\varepsilon_p.
\end{aligned}
$$

A uniform $\Omega_v$ permits different $\Omega_n$ in the core and inner crust.
Their $x_n$, $\varepsilon_n$ and drag ratios are prescribed piecewise constants.
The outer crust has $x_n=0$. The charged component interpenetrates the core;
it is not just a geometrical outer shell. $I_v$ describes an angular-momentum
reservoir, not the inertia of microscopic vortex cores.

The input checks require $0\leq x_n<1$, $\varepsilon_n<x_p=1-x_n$ and
nonnegative drag, giving a locally positive inertia matrix. All neutrons counted
by $x_n$ are assumed superfluid. APR supplies neither these fractions nor pairing,
entrainment or drag. The separate BSk24 composition is not imported.

## 2. Straight-vortex line integrals

`build_reservoir` reads the saved APR rest density, boundaries and Hartle inertia.
For a line at $x<R_d$, let $z_d(x)=\sqrt{R_d^2-x^2}$,
$r=\sqrt{x^2+z^2}$ and $\kappa=h_{\rm P}/(2m_n)$, where $h_{\rm P}$ is Planck's
constant. Reflection symmetry supplies a factor of two:

$$
\begin{aligned}
a(x)&=2\int_0^{z_d}\rho_n(r)\,dz,\\
b(x)&=2\int_0^{z_d}\frac{\rho_n(r)}{1-\varepsilon_n(r)}\,dz,\\
d(x)&=\frac2\kappa\int_0^{z_d}\eta(r)\,dz
=2\int_0^{z_d}\mathcal R(r)\rho_n(r)\,dz,\qquad
\mathcal R=\frac{\eta}{\kappa\rho_n}.
\end{aligned}
$$

$a,b,d$ have units $\mathrm{kg\,m^{-2}}$. The free-vortex coefficients and
velocities are

$$
\begin{aligned}
\mathcal B^x&=\frac{db}{d^2+a^2},&
\mathcal B^\theta&=\frac{ab}{d^2+a^2},\\
v_L^x\big|_{\rm free}&=x\mathcal B^x\omega,&
v_L^\theta\big|_{\rm free}&=x(\Omega_p+\mathcal B^\theta\omega).
\end{aligned}
$$

A constant $Y\in[0,1]$ represents the mobile fraction of vortices, so the
population-averaged radial drift is $Yx\mathcal B^x\omega$. $Y$ is an assumed
mobility, not a pinning-force or thermal-creep calculation.

The Newtonian inertia integrals are

$$
\begin{aligned}
I_{v,N}&=\frac{8\pi}{3}\int_0^{R_d}\frac{r^4\rho_n}{1-\varepsilon_n}\,dr
=2\pi\int_0^{R_d}x^3b(x)\,dx,\\
I_N&=\frac{8\pi}{3}\int_0^R r^4\rho_b\,dr,\\
\langle\mathcal B^x\rangle&=
\frac{2\pi}{I_{v,N}}\int_0^{R_d}x^3b(x)\mathcal B^x(x)\,dx,\qquad
\mathcal B_{\rm eff}=Y\langle\mathcal B^x\rangle.
\end{aligned}
$$

Both axial and cylindrical Gaussian quadratures split at the core interface to
resolve changes in microphysics and lines tangent to the core. The radial and
cylindrical evaluations of $I_{v,N}$ are independent consistency checks.
The code derives the inertia fraction from the radial integral but normalizes
the drag average by its cylindrical integral to avoid propagating their small
quadrature mismatch into $\mathcal B_{\rm eff}$.

The background uses areal radius and $\rho_b$ with a **Newtonian volume element**.
A common factor $s_I=I_{\rm Hartle}/I_N$ sets $I=I_{\rm Hartle}$ and
$I_v=s_I I_{v,N}$. It cancels from the inertia fractions and drag average.
This convention does not supply relativistic entrainment, drag or two-fluid energy.

For uniform microphysics,

$$
\mathcal B^x=\frac{\mathcal R}{(1+\mathcal R^2)(1-\varepsilon_n)}.
$$

Radial mobility vanishes for both zero and infinite microscopic drag. At very
large drag the vortices follow the charged matter and cannot migrate radially.
Entrainment changes velocities and inertias without itself dissipating energy.

## 3. Rigid projection and heat

The AP vorticity equation is

$$
\partial_t\Omega_v=-Y\mathcal B^x
(2\Omega_v+x\partial_x\Omega_v)(\Omega_v-\Omega_p).
$$

Projecting onto uniform $\Omega_v$ while retaining its nonlinear factor gives

$$
\begin{aligned}
\dot\Omega_v&=-2\mathcal B_{\rm eff}\Omega_v\omega,\\
I_p\dot\Omega_p&=N-I_v\dot\Omega_v,\\
J_{\rm spin}&=I_p\Omega_p+I_v\Omega_v,\qquad
\tau_{\rm lag}=\frac{I_p}{2I\mathcal B_{\rm eff}\Omega_v}.
\end{aligned}
$$

$N$ is signed: positive torque spins up the charged component.
$\tau_{\rm lag}$ is the local small-lag relaxation time, not a constant assumed
throughout the inspiral. There is no intrinsic magnetic braking.
The canonical-variable energy and power satisfy

$$
\begin{aligned}
E_{\rm spin}&=\tfrac12(I_p\Omega_p^2+I_v\Omega_v^2),\\
\dot E_{\rm spin}&=\Omega_pN-P_{\rm heat},\\
P_{\rm heat}&=-I_v\omega\dot\Omega_v
=2I_v\mathcal B_{\rm eff}\Omega_v\omega^2\geq0
\quad(\Omega_v>0).
\end{aligned}
$$

The same power follows from vortex drag. For a free vortex,

$$
|\mathbf v_L-\mathbf v_p|^2=\frac{x^2\omega^2b^2}{a^2+d^2},\qquad
\int_{-z_d}^{z_d}\eta|\mathbf v_L-\mathbf v_p|^2dz
=\kappa d|\mathbf v_L-\mathbf v_p|^2,\qquad
n_v=\frac{2\Omega_v}{\kappa}.
$$

Integrate over vortex positions and multiply by $Y$ for the mobile population.
Pinned vortices comove with the crust and dissipate no drag heat in this
idealization. Permanent $Y=0$ pinning is a limiting control; it does not model
finite critical lag, unpinning events or glitches.

At fixed isolated-star angular momentum, the exact heat released is

$$
Q=\frac{I_pI_v}{2I}(\omega_i^2-\omega_f^2).
$$

This identity is tested independently. During driven spin-up, $\omega<0$
corresponds to inward vortex drift and increasing $n_v$.
The exported signed speed is defined operationally as

$$
v_{\rm diagnostic}=Y\omega
\left[\sum_j w_j(x_j\mathcal B^x_j)^2\right]^{1/2},\qquad
\sum_j w_j=1,
$$

where $w_j$ are canonical-inertia quadrature weights. It is a signed rms of the
mobility-averaged drift over cylindrical radius, not a number-weighted rms speed
of individual pinned and mobile vortices. A rigid control sets this diagnostic
to zero and copies its single spin into the canonical-frequency column.

## 4. Quadrupole and circular orbital reaction

The fixed equatorial quadrupole is $\Delta I=\epsilon I$ with total $I$.
Its orientation follows $\theta$, the charged component's figure angle:

$$
U(r,\gamma)=-\frac34\frac{GM'\Delta I}{r^3}\cos\gamma,\qquad
N=-\frac{\partial U}{\partial\theta}
=-\frac32\frac{GM'\Delta I}{r^3}\sin\gamma.
$$

Normalizing the same $\epsilon$ to $I_p$ would change the imposed quadrupole.
In the rigid limit with negligible orbital reaction,

$$
\ddot\gamma+3\epsilon\frac{GM'}{r^3}\sin\gamma
+2\dot\Omega_{\rm orb}=0.
$$

To include reaction, the code uses orbital angular momentum $L$ as the circular
Hamiltonian action. It allows the circular radius to change with $L$:

$$
\begin{aligned}
M_t&=M+M',&\mu&=\frac{MM'}{M_t},\\
r(L)&=\frac{L^2}{\mu^2GM_t},&n(L)&=\sqrt{\frac{GM_t}{r(L)^3}},\\
L(n)&=\mu(GM_t)^{2/3}n^{-1/3},&E_K(L)&=-\frac{GM_t\mu}{2r(L)},\\
H&=E_K+U+E_{\rm spin},&
\Omega_{\rm orb}&=\frac{\partial H}{\partial L}=n-\frac{6U}{L}.
\end{aligned}
$$

Removing angular momentum increases the Keplerian frequency; a fixed
$\mu r^2$ orbital inertia would give the wrong reaction on a circular sequence.
The quadrupole potential is evaluated on the unperturbed circular relation.
The small $-6U/L$ correction keeps the selected Hamiltonian energy accounting
consistent; this remains a leading-order circular reduction.

With leading quadrupole GW luminosity,

$$
\begin{aligned}
\mathcal F_{\rm GW}&=\frac{32}{5}\frac{G^{7/3}}{c^5}
\mu^2M_t^{4/3}n^{10/3},\\
\dot L&=-N-\frac{\mathcal F_{\rm GW}}{\Omega_{\rm orb}},\\
\dot n&=\frac{3n}{L}\left(N+\frac{\mathcal F_{\rm GW}}{\Omega_{\rm orb}}\right),\\
\dot\gamma&=2(\Omega_p-\Omega_{\rm orb}).
\end{aligned}
$$

These equations imply

$$
\dot H+\mathcal F_{\rm GW}+P_{\rm heat}=0,\qquad
\frac{d}{dt}(L+J_{\rm spin})+
\frac{\mathcal F_{\rm GW}}{\Omega_{\rm orb}}=0.
$$

Together with the spin-only balance, these are integrated numerical audits.
The code records both
$W_\star=\int\Omega_pN\,dt$ and $W_{\rm orbit}=\int\Omega_{\rm orb}N\,dt$.
Their difference reflects the changing interaction energy and circular radius.
Heat enters only through internal friction; adding it again to the orbital
torque would double-count energy transfer. Integrated heat is not a temperature.

The point-particle reference sets $U=N=0$. Post-Newtonian terms, eccentricity,
dynamical tides and radiation from stellar rotation are omitted.
Plot abscissae generally use $n/\pi$ as GW frequency, while phase-slip summaries
use $\Omega_{\rm orb}/\pi$. Both orbital frequencies are exported; their small
difference is the conservative correction above.

## 5. Initial conditions and capture diagnostics

`paper` starts with $n_0/(2\pi)=0.4\,\mathrm{Hz}$,
$\Omega_{p,0}/(2\pi)=\Omega_{v,0}/(2\pi)=0.5\,\mathrm{Hz}$ and
$\gamma_0=-2(0.3986)$. Capture is determined by evolution.
`locked` starts at $n_0/(2\pi)=0.6\,\mathrm{Hz}$ with

$$
\sin\gamma_*=-\frac{2\dot n_{\rm GW}}{3\epsilon(M'/M_t)n^2},\qquad
\cos\gamma_*>0,\qquad \Omega_{p,0}=\Omega_{v,0}=\Omega_{{\rm orb},0}.
$$

This is the stable rigid equilibrium phase with the point-particle chirp;
finite backreaction and the developing two-component lag can cause a transient.
Every drag-sweep case starts independently this way, so it tests evolution
conditional on capture rather than capture probability.

Capture requires a phase range below $2\pi$ and mean fractional charged-spin
mismatch below $0.02$ in the Kepler orbital-frequency window $0.8$–$1.0\,\mathrm{Hz}$.
The phase reference is the first stored sample at or above $1\,\mathrm{Hz}$.
The first subsequent absolute excursion of $2\pi$ is a full phase slip.
A missing slip means no detected slip by the endpoint, or no established capture;
read the capture flag too. This criterion is distinct from the first loss of
frequency tracking or the disappearance of a quasistatic fixed point.

The ideal rigid torque-limit GW frequency is

$$
f_{\rm br}=\frac{c^3}{\pi G(M+M')}
\left[\frac{5(1+q)}{64}\right]^{3/5}\epsilon^{3/5}.
$$

The instantaneous two-component load depends on lag; changing a friction time
in this rigid formula does not calculate its breaking point.
The phase diagnostic is the time-domain accumulated GW phase at matched Kepler
frequency relative to the point-particle trajectory, with equal initial phase:

$$
\Delta\Phi(n)=2\int_{n_0}^n
\left(\frac{\Omega_{\rm orb}(\tilde n)}{\dot n(\tilde n)}
-\frac{\tilde n}{\dot n_{\rm pp}(\tilde n)}\right)d\tilde n.
$$

It combines conservative spin-up, altered locking and heating. It is not an
isolated dissipative waveform term or the paper's stationary-phase correction.

## 6. Bound from available orbital torque

For fixed masses and imposed quadrupole,

$$
A(n)=\frac32\frac{M'}{M_t}\Delta I\,n^2=\max|N|,\qquad
\alpha(n)=\frac{A(n)n}{\mathcal F_{\rm GW}(n)}\propto n^{-1/3}.
$$

Ignoring the small conservative frequency correction initially,

$$
\dot n\geq\frac{3[\mathcal F_{\rm GW}-nA]}{L}>0,\qquad
|W_{\rm orbit}|\leq\int_{n_0}^{n_1}
\frac{A(n)nL(n)}{3[\mathcal F_{\rm GW}(n)-nA(n)]}\,dn.
$$

This permits either torque sign and uses the smallest allowed chirp. The
Hamiltonian correction obeys

$$
\left|\frac{\Omega_{\rm orb}}n-1\right|
\leq\delta(n)=\frac92\frac{M'}{M_t}\frac{\Delta I\,n}{L(n)}.
$$

The code evaluates the corrected bound

$$
|W_{\rm orbit}|\leq\int_{n_0}^{n_1}
\frac{A(n)n[1+\delta(n)]^2L(n)}
{3\{\mathcal F_{\rm GW}(n)-A(n)n[1+\delta(n)]\}}\,dn,
$$

which requires a positive denominator over the interval. This condition holds
for the supplied inputs; arbitrary new configurations need it checked before
interpreting the bound. Initial spin and interaction energy also belong in any
bound on total heat: net positive orbital work alone does not generally bound it.
The interaction contribution is limited by $|U_i|+|U_f|$.
This bound applies only to the imposed quadrupole channel, not all stellar
heating mechanisms.

## 7. Numerical state and review limits

`Binary.quantities` divides torques, energies and angular momenta by the saved
stellar $I$. Its `ell`, `u`, `torque` and `flux` mean $L/I$, $U/I$, $N/I$ and
$\mathcal F_{\rm GW}/I$. Physical SI outputs restore $I$.
The integration uses $s=\ln n$, with $dt/ds=n/\dot n$:

| Index | `evolve` state | Meaning |
|---|---|---|
| 0–2 | `phase`, `op`, `ov` | $\gamma,\Omega_p,\Omega_v$ |
| 3 | `time` | Elapsed seconds |
| 4–6 | `heat`, `egw`, `jgw` | $Q/I,E_{\rm GW}/I,J_{\rm GW}/I$ |
| 7–8 | `wstar`, `worb` | $W_\star/I,W_{\rm orbit}/I$ |
| 9 | `dephase` | $\Delta\Phi$ |

`rhs_time` is a separate test interface whose first four entries are
$[n,\gamma,\Omega_p,\Omega_v]$, followed by the five accumulated budgets;
its ordering must not be substituted for the table above.

Monotone $\dot n>0$ and positive $\Omega_v$ are enforced. Paper-start circulation
uses explicit DOP853; locked trajectories use implicit Radau for stiff drag
relaxation. A complex-step Jacobian differentiates the smooth RHS without
finite-difference cancellation in zero diagnostic components.
Diagnostics combine adaptive solver nodes with a dense uniform log-frequency
grid. CSV export retains at most 15,000 rows; audits precede that decimation.

Validation checks all three budgets, nonnegative heat, tighter evolution
tolerances, doubled line/cylindrical quadrature and preservation of the capture
flag. Tests independently cover isolated relaxation, the closed Hamiltonian,
point-particle chirp and zero/infinite drag. See `validation.json` for actual
errors and limits.

The rigid projection cannot resolve differential canonical rotation, vortex
pile-up, bending, turbulence, finite-strength unpinning or thermal creep.
Large lag or strongly variable mobility may invalidate the approximation.
Pairing, temperature, cooling, magnetic stress evolution and changes in the
quadrupole are not evolved. APR boundaries remain assumed density labels.
Convergence verifies the selected equations, not these physical assumptions.
