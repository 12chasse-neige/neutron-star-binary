# APR equilibrium model

This document maps [`nsbinary/structure.py`](../nsbinary/structure.py) and [`scripts/compute_structure.py`](../scripts/compute_structure.py) to their equations. The background is cold, spherical and isotropic. First-order slow rotation and linear static tides do not alter its hydrostatic density profile. Input provenance and interpolation are in [data/README.md](../data/README.md).

## State and units

The metric convention in $G=c=1$ is

$$
ds^2=-e^{\nu(r)}dt^2+(1-2m/r)^{-1}dr^2+r^2d\Omega^2.
$$

Here $r$ is areal radius and $m$ is enclosed gravitational mass in length units. $\mathcal E$ includes rest energy, while $\rho_b$ is baryonic rest-mass density; these obey $\rho_b\neq\mathcal E/c^2$ in general.

| State index | Code variable | Meaning |
|---|---|---|
| 0 | `gravitational_mass` | Gravitational mass, $\mathrm{km}$ |
| 1 | `pseudo_enthalpy` | Dimensionless pseudo-enthalpy |
| 2 | `frame_dragging` | Normalized $\bar\omega=\Omega-\omega_{\rm fd}$ |
| 3 | `frame_dragging_derivative` | $d\bar\omega/dr$ |
| 4 | `tidal_y` | $rH'/H$ for static quadrupolar tides |
| 5 | `proper_radius` | Proper radial distance $\ell$, $\mathrm{km}$ |

`EquationOfState.at_h` returns $P,\mathcal E$ in $\mathrm{km^{-2}}$, but $\rho_b$ in $\mathrm{kg\,m^{-3}}$. SI pressure and energy density convert through

$$
P_{\rm geom}=10^6\frac{G}{c^4}P_{\rm SI},\qquad
\mathcal E_{\rm geom}=10^6\frac{G}{c^4}\mathcal E_{\rm SI}.
$$

`DENSITY_TO_GEOM` uses $10^6G/c^2$ for mass density instead. Divide $\mathrm{kg\,m^{-3}}$ by $1000$ to export $\mathrm{g\,cm^{-3}}$. Also $1\,\mathrm{Pa}=10\,\mathrm{dyn\,cm^{-2}}$.

## Hydrostatic equations and the centre

LAL reconstructs pseudo-enthalpy and rest density from the tabulated barotrope. The solver evolves

$$
dh=\frac{dP}{\mathcal E+P},\qquad
\frac{dm}{dr}=4\pi r^2\mathcal E,\qquad
\frac{dh}{dr}=-\frac{m+4\pi r^3P}{r^2(1-2m/r)},\qquad
\frac{d\ell}{dr}=(1-2m/r)^{-1/2}.
$$

The regular centre has $m(0)=0$, $h(0)=h_c$, $\bar\omega(0)=1$ and $y(0)=2$. `integrate` begins at $r_0>0$ with

$$
\begin{aligned}
m(r_0)&=\frac{4\pi}{3}\mathcal E_c r_0^3,\\
h(r_0)&=h_c-2\pi\left(\frac{\mathcal E_c}{3}+P_c\right)r_0^2,\\
\bar\omega(r_0)&=1+\frac{8\pi}{5}(\mathcal E_c+P_c)r_0^2,\\
\bar\omega'(r_0)&=\frac{16\pi}{5}(\mathcal E_c+P_c)r_0,\qquad
y(r_0)=2,\qquad \ell(r_0)=r_0.
\end{aligned}
$$

The last two values retain their leading regular terms. Defaults are $r_0=10^{-4}\,\mathrm{km}$ and a surface cutoff $h_s=10^{-10}$. DOP853 terminates at the decreasing event $h(R)=h_s$, approximating $P=0$. It uses dense output and component-specific absolute tolerances.

`solve_target` scans $0.06\leq h_c\leq\min(0.58,0.95h_{\max})$ for the first crossing of the requested mass before a sampled mass decrease. Brent's method solves $M(h_c)=M_{\rm target}$. This is a supported ordinary neutron-star branch, not a global maximum-mass search or an arbitrary low-mass solver.

## Hartle slow rotation

For $j=e^{-(\nu+\lambda)/2}$ and $e^\lambda=(1-2m/r)^{-1}$,

$$
\frac{j'}j=-\frac{4\pi r(\mathcal E+P)}{1-2m/r},\qquad
\bar\omega''+\left(\frac4r+\frac{j'}j\right)\bar\omega'
+\frac4r\frac{j'}j\bar\omega=0.
$$

Matching to the vacuum solution $\bar\omega=\Omega-2J/r^3$ gives

$$
I_{\rm geom}=\frac J\Omega
=\frac{R^4\bar\omega'(R)}{6[\bar\omega(R)+R\bar\omega'(R)/3]}.
$$

The arbitrary central normalization cancels. If $I_{\rm geom}$ is in $\mathrm{km^3}$, then $I_{\rm SI}=(c^2/G)10^9I_{\rm geom}$ in $\mathrm{kg\,m^2}$. The solar-mass geometric value divides by $(GM_\odot/c^2)^3$, with all lengths in the same unit. Second-order rotational distortion is omitted.

## Static quadrupolar tides

For $y=rH'/H$, define

$$
\begin{aligned}
f_r&=1-2m/r,& a_r&=\frac{m+4\pi r^3P}{r^2f_r},\\
F&=\frac{1-4\pi r^2(\mathcal E-P)}{f_r},&
Q&=\frac{4\pi[5\mathcal E+9P+(\mathcal E+P)d\mathcal E/dP]}{f_r}
-\frac6{r^2f_r}-4a_r^2.
\end{aligned}
$$

The equation is $ry'+y^2+Fy+r^2Q=0$, with $y(0)=2$. With $C_\star=m(R)/R$ and $y_R=y(R)$, `Star.love_k2` evaluates

$$
\begin{aligned}
k_2&=\frac{\frac85C_\star^5(1-2C_\star)^2
[2+2C_\star(y_R-1)-y_R]}{D},\\
D&=2C_\star[6-3y_R+3C_\star(5y_R-8)]\\
&\quad+4C_\star^3[13-11y_R+C_\star(3y_R-2)+2C_\star^2(1+y_R)]\\
&\quad+3(1-2C_\star)^2[2-y_R+2C_\star(y_R-1)]\ln(1-2C_\star),\\
\Lambda&=\frac{2k_2}{3C_\star^5}.
\end{aligned}
$$

No finite-density jump condition is imposed. A self-bound surface or explicit sharp phase transition needs corresponding tidal matching conditions.

## Layers and exports

`Star.boundary` inverts $\rho_b(h)$ for each configured threshold, then solves $h(r)=h_{\rm threshold}$ using dense output. It rejects thresholds outside the star. The regions are $[0,r_{\rm cc}]$, $[r_{\rm cc},r_{\rm drip}]$ and $[r_{\rm drip},R]$. These are density labels, not predicted phase transitions.

For shell boundaries $r_1,r_2$,

$$
\Delta r=r_2-r_1,\qquad
\Delta\ell=\int_{r_1}^{r_2}\frac{dr}{\sqrt{1-2m(r)/r}},\qquad
\Delta M=m(r_2)-m(r_1).
$$

$\Delta M$ is gravitational mass, not proper-volume baryonic mass. `Star.profile` combines a whole-star grid, a denser crust grid and exact boundary radii. The exported $dm/dr$ is converted to $M_\odot\,\mathrm{km^{-1}}$.

## Physical scope

The background assumes cold, isotropic, nonrotating spherical hydrostatic equilibrium. First-order Hartle rotation yields $I$ without changing the TOV density or pressure profile; linear static tides yield $k_2$ and $\Lambda$ without density jumps. Second-order rotation, finite temperature, magnetic fields, accretion and dynamical tides are omitted.

Density thresholds only label radial layers; they are not computed nuclear phase transitions. The model does not predict a magnetic mountain, crustal strain, species fractions, pairing gaps or the permanent quadrupole $\Delta I=\epsilon I$ imposed by the inspiral calculation.

Sources: [Akmal, Pandharipande & Ravenhall (1998)](https://arxiv.org/abs/nucl-th/9804027); [Chamel & Haensel (2008)](https://arxiv.org/abs/0812.3955); [Hartle (1967)](https://doi.org/10.1086/149400); [Hinderer (2008)](https://arxiv.org/abs/0711.2420); [LAL EOS documentation](https://lscsoft.docs.ligo.org/lalsuite/dev/lalsimulation/group___l_a_l_sim_neutron_star_e_o_s__c.html).
