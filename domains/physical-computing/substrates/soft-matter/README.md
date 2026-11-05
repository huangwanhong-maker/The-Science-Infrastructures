# soft-matter: Ionic Membranes, Hydrogels, Liquid Crystals and Active Matter

*Status: planned (ADR 0076).* Modes: reservoir, analogue, learning. **Depth:** genuine many-body emergence (active matter) and continuum instabilities (liquid crystals). Provisional until condensed-matter is founded (ADR 0076 §2.5).

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Ionic membranes and iontronic memristors (Kamsma et al. 2023) | volatile memory, excitability, spiking in water | water, salts, polymer channels; reusable |
| Hydrogels and self-oscillating gels (Yoshida et al. 1996; Yashin and Balazs 2006) | chemistry turned into motion; swelling waves | alginate and agarose are biodegradable, polyacrylamide's monomer is toxic; the gel's BZ reagents are consumed |
| Liquid crystals; active nematics (Sanchez et al. 2012) | the Freedericksz transition, electroconvection, defect turbulence; active turbulence | low voltage, reusable; liquid-crystal monomers are persistent and toxic (Su et al. 2019); active nematics consume ATP |
| Active colloids (Bricard et al. 2013; Palacci et al. 2013) | flocking, living crystals, symmetry breaking | hydrogen peroxide as fuel, or hexadecane as solvent |

## 2. The Simulators

1. **The iontronic memristor.** A volatile conductance, $\dot g = (g_\infty(V) - g)/\tau$, with $g_\infty$ from the channel's accumulation of ions (provisional, after Kamsma et al.). Checks:
   - the loop is pinched;
   - the lobes' area vanishes at both low and high frequency and is greatest near $\omega\tau \sim 1$.
2. **The self-oscillating gel.** The Oregonator ([`chemical`](../chemical/)) drives the swelling, $\dot\phi = (\phi_{eq}(v) - \phi)/\tau_g$: the gel is a low-pass mechanical readout of its chemistry. Check: the amplitude's transfer $1/\sqrt{1 + \omega^2\tau_g^2}$. Yashin and Balazs's lattice model of the gel comes later.
3. **The liquid-crystal cell.** The director's angle obeys $\gamma\theta_t = K\theta_{zz} + \varepsilon_0\Delta\varepsilon E^2 \sin\theta\cos\theta$, in one elastic constant. Checks:
   - the Freedericksz threshold $V_c = \pi\sqrt{K/(\varepsilon_0\Delta\varepsilon)}$, independent of the thickness;
   - the pitchfork amplitude above the threshold, growing as $\sqrt{V - V_c}$.
4. **Active nematics.** Two-dimensional active nematohydrodynamics: a $Q$-tensor coupled to Stokes flow. Later, on physis's fluids. Check: the length of the active instability, $\ell_a = \sqrt{K/\alpha}$ (after Ramaswamy and Simha 2002).
5. **Active colloids.** Vicsek's model (physis `statistical_mechanics`) as their normal form, with the parameters of Quincke rollers, provisional. Checks:
   - the order parameter tends to 1 at zero noise, and to $O(1/\sqrt N)$ at full noise;
   - the transition against Vicsek et al. (1995), qualitatively.

## 3. Tasks and their Lowering

- **Filters** (simulators 1 and 2).
- **Phenomena** (simulators 3 to 5).
- **Learning:** networks of iontronic memristors.

## 4. What the Foundations Must Gain

- physis `statistical_mechanics`: Vicsek's model;
- physis `continuum`: the active nematic, after its fluids.

## 5. References

- Yoshida et al. (1996); Sanchez et al. (2012); Bricard et al. (2013); Palacci et al. (2013); Su et al. (2019): in the theory note.
- de Gennes, P.-G. and Prost, J. (1993), *The Physics of Liquid Crystals*, 2nd edition, Oxford University Press.
- Kamsma, T. M., Boon, W. Q., ter Rele, T., Spitoni, C. and van Roij, R. (2023), Iontronic neuromorphic signaling with conical microfluidic memristors, *Physical Review Letters* 130, 268401.
- Ramaswamy, S. and Simha, R. A. (2002), Hydrodynamic fluctuations and instabilities in ordered suspensions of self-propelled particles, *Physical Review Letters* 89, 058101.
- Vicsek, T., Czirók, A., Ben-Jacob, E., Cohen, I. and Shochet, O. (1995), Novel type of phase transition in a system of self-driven particles, *Physical Review Letters* 75, 1226–1229.
- Yashin, V. V. and Balazs, A. C. (2006), Pattern formation and shape changes in self-oscillating polymer gels, *Science* 314, 798–801.
