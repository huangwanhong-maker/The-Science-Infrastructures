# fluids: Water and Interfaces as Computers

*Status: draft (2026-11-06): the nine simulators below are written, each with its check program (59 checks), on physis's decision 0022. Changes from the plan: the tank's forcing conserves volume, since a plunger adds no water, so the echo state property holds for every departure but the mean level; the rotating tank is a channel by images, since a periodic box admits modes the walls forbid, and its transition is stated per azimuthal wavenumber, Θ_m = 16/(π² + (mL/R)²).* Modes: analogue, reservoir, programmed. **Depth:** the deepest family for the atmosphere and the ocean. Water obeys the Navier–Stokes equations, and a tank stands for its target by similarity (theory note §1.6). The greenest family: water, heat, soap and a turntable, all reusable. Its readings, mostly by camera, carry the cost of conversion.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Water in a tank (Fernando and Sojakka 2003) | interfering waves as a reservoir | water, motors, a camera |
| Convection cells | rolls, hexagons, spiral defect chaos; Lorenz's chaos as their truncation | water or oil, heat |
| Rotating, stratified tanks (Fultz et al. 1959; Hide and Mason 1975) | baroclinic waves, vacillation, geostrophic turbulence: the general circulation in a tank | water, heat, a turntable |
| Faraday waves; walking droplets (Couder and Fort 2006; Bush 2015) | parametric instability, quasipatterns; pilot-wave dynamics | water or silicone oil on a shaker |
| Capillary and foam networks | imbibition; coarsening and topological rearrangements | soap, water, air; a foam coarsens for good, but is cheaply remade |
| Microfluidic bubble logic (Prakash and Gershenfeld 2007) | bubbles as bits; logic, memory, timing | channels in polymer, air, water |
| Phillips's hydraulic machine (Phillips 1950) | the circular flow of income, in water | tanks and valves |

## 2. The Simulators

1. **The tank.** A damped linear wave equation of the surface on a grid (physis `continuum`), driven by motors as point sources, read by a camera as a thresholded image. Checks:
   - mode frequencies $\omega_{mn} = c\pi\sqrt{(m/L_x)^2 + (n/L_y)^2}$;
   - the damping rate;
   - the balance of energy.
2. **The convection cell.** The linear stability of Boussinesq convection (physis `continuum`), and Lorenz's three-mode truncation (Saltzman 1962; Lorenz 1963) mapped from the cell's Rayleigh number, Prandtl number and wavenumber $a$: $r = Ra/Ra_c^{\text{free}}$, $\sigma = Pr$, $b = 4/(1 + a^2)$. Checks:
   - $Ra_c = 1707.76$ at $a = 3.117$ between rigid plates, and $27\pi^4/4$ at $a = \pi/\sqrt2$ between free ones (Chandrasekhar 1961);
   - the map returns Lorenz's $b = 8/3$ at $a^2 = 1/2$.
3. **The rotating tank.** Two-layer quasi-geostrophic flow in a periodic channel (physis `continuum`), relaxed towards an imposed shear and damped by Ekman friction. The tank's rotation, temperature contrast, depth, gap, viscosity and expansion are mapped to its thermal Rossby and Taylor numbers, and so to the model's deformation wavenumber and forcing. Checks:
   - Phillips's growth rates $\sigma = kU\sqrt{(k_d^2 - K^2)/(k_d^2 + K^2)}$ for $K < k_d$;
   - energy and enstrophy conserved when forcing and friction are switched off;
   - the sequence of regimes against Hide and Mason's diagram, qualitatively.
4. **Faraday waves.** Each surface mode is a damped Mathieu oscillator, $\omega_0^2 = (gk + \sigma k^3/\rho)\tanh kh$, forced at $\Omega$. Floquet multipliers (physis `dynamics`) give the tongues of instability. Checks:
   - the subharmonic response at $\Omega = 2\omega_0$;
   - the threshold against Benjamin and Ursell's tongues at small damping.
5. **Walking droplets.** The stroboscopic model of Oza, Rosales and Bush (2013), $\kappa\ddot x + \dot x = \beta\int_{-\infty}^t J_1(|x(t) - x(s)|)\,\frac{x(t) - x(s)}{|x(t) - x(s)|}\,e^{-(t - s)}\,ds$. Checks:
   - walking begins at $\beta = 2$, at the steady speed $u = \frac{1}{\sqrt2}\sqrt{2\beta - 1 - \sqrt{1 + 4\beta}}$;
   - the model's interference is not taken as shown, after Andersen et al. (2015).
6. **Capillary networks.** Washburn's imbibition, $L^2 = \gamma r\cos\theta\, t/(2\mu)$, in each tube, and Hagen–Poiseuille conductances $\pi r^4/(8\mu L)$ in a Kirchhoff network (physis `discretisation`). Checks: Washburn's law; conductances in series and in parallel.
7. **Foams.** Two-dimensional coarsening by a Potts model (physis `statistical_mechanics`; Glazier, Anderson and Grest 1990). Checks:
   - von Neumann's law, $\dot A \propto (n - 6)$, in the mean over cells of $n$ sides;
   - the mean area growing linearly in time.
8. **Bubble logic.** A discrete-event model of bubbles at junctions, each bubble raising its channel's resistance. Provisional. Checks: the truth tables of AND, OR and NOT; bubbles conserved.
9. **Phillips's machine.** Tanks and linear valves conjugate to an income–expenditure model. Checks:
   - the steady income equals autonomous spending over $1 - c$, the multiplier;
   - the time constant;
   - the conjugacy of levels and flows to income and spending.

## 3. Tasks and their Lowering

- **Phenomena and flows** (simulators 2 to 5): the similarity map from target to tank, by its dimensionless groups.
- **Filters** (simulator 1).
- **Functions** (simulator 8).
- **Analogues** of models (simulator 9).

## 4. What the Foundations Must Gain

- physis `continuum`: linear waves in a basin, the linear stability of convection, two-layer quasi-geostrophic flow;
- physis `dynamics`: Floquet multipliers;
- physis `statistical_mechanics`: the Potts model;
- physis `discretisation`: Kirchhoff networks.

## 5. References

- Fernando and Sojakka (2003); Chandrasekhar (1961); Fultz et al. (1959); Hide and Mason (1975); Couder and Fort (2006); Bush (2015); Andersen et al. (2015); von Neumann (1952); Phillips (1950): in the theory note.
- Benjamin, T. B. and Ursell, F. (1954), The stability of the plane free surface of a liquid in vertical periodic motion, *Proceedings of the Royal Society A* 225, 505–515.
- Glazier, J. A., Anderson, M. P. and Grest, G. S. (1990), Coarsening in the two-dimensional soap froth and the large-Q Potts model: a detailed comparison, *Philosophical Magazine B* 62, 615–645.
- Lorenz, E. N. (1963), Deterministic nonperiodic flow, *Journal of the Atmospheric Sciences* 20, 130–141.
- Oza, A. U., Rosales, R. R. and Bush, J. W. M. (2013), A trajectory equation for walking droplets: hydrodynamic pilot-wave theory, *Journal of Fluid Mechanics* 737, 552–570.
- Phillips, N. A. (1954), Energy transformations and meridional circulations associated with simple baroclinic waves in a two-level, quasi-geostrophic model, *Tellus* 6, 273–286.
- Prakash, M. and Gershenfeld, N. (2007), Microfluidic bubble logic, *Science* 315, 832–835.
- Saltzman, B. (1962), Finite amplitude free convection as an initial value problem, *Journal of the Atmospheric Sciences* 19, 329–341.
- Washburn, E. W. (1921), The dynamics of capillary flow, *Physical Review* 17, 273–283.
