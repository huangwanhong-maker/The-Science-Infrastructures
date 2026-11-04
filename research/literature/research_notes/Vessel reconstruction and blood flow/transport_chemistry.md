# Transport and Chemical Reactions Carried by Blood in Vessels: Oxygen, Mass Transport, Coagulation and Thrombosis, Drugs and Contrast, Perfusion

Notes for a team implementing advection–diffusion–reaction (ADR) models coupled to simulated blood flow (0D/1D networks to 3D), with reaction networks for coagulation and oxygen binding, verified against published reference results. Compiled 2026-11-04.

**Verification legend.**
- **[V]** means the bibliographic record (authors, year, title, venue, volume, issue, pages, DOI) was checked against the Crossref API on 2026-11-04. The source link is the DOI.
- **[W]** means it was checked only by web search or a repository or publisher listing.
- **[U]** means it was not verified and must be checked before citing.
- **†** marks a content summary (an equation, model scope or result) taken from my knowledge of the paper and not re-read in full this session. Check it against the paper before using it as a test fixture.
- Without †, the content claim comes from the page linked.

**Classical vs latest.** Within each question, "Classical" means work up to 2022. "Latest" means 2023–2026. A few 2021–2022 papers that the latest work builds on directly are marked "bridge".

---

## 1. Oxygen Transport: Haemoglobin Binding, the Krogh Cylinder, Arteries and Their Wall, the Microcirculation, and High-Schmidt Boundary Layers

### Takeaway
A complete oxygen module needs three pieces:
- **A saturation law.** Hill or Adair for theory. Severinghaus 1979 or Dash–Bassingthwaighte for physiological accuracy, with pH, CO₂, temperature and 2,3-DPG effects.
- **A large-vessel ADR equation.** Oxyhaemoglobin is carried with the plasma, so the slope dS/dP multiplies the convective term (Moore & Ethier 1997; Tarbell 2003). It needs boundary-layer resolution, because Sc ≈ 10³ for O₂.
- **A microvascular scheme.** Either the Krogh cylinder for verification, or network methods: Secomb's Green's functions, 3D–1D coupling, or porous continua.

Closed-form fixtures exist: the Hill and Severinghaus curves, the Krogh–Erlang profile, and Graetz/Lévêque mass-transfer solutions. Network reference results are published (e.g., Celaya-Alcala 2021), but the code licences are unclear. Recent work (2023–2026) resolves individual red cells (Takeishi 2026; Zhang 2025) and coupled coronary or cortical networks (Wang 2025; Sharma & Secomb 2025).

### Cited Findings

#### 1.1 Haemoglobin binding (classical)
- **Hill (1910).** "The possible effects of the aggregation of the molecules of haemoglobin on its dissociation curves." *J Physiol* 40 (Suppl., Proceedings of the Physiological Society, 22 Jan 1910): iv–vii.
  - Crossref resolves doi:10.1113/jphysiol.1910.sp001386 to the Proceedings as a whole [V for the Proceedings]. The communication's exact title and pages are [U].
  - This is the origin of the Hill equation S = Pⁿ/(Pⁿ + P₅₀ⁿ), with n ≈ 2.7 for adult human blood †. — [Source](https://doi.org/10.1113/jphysiol.1910.sp001386)
- **Adair GS (with Bock AV, Field H) (1925).** "The hemoglobin system. VI. The oxygen dissociation curve of hemoglobin." *J Biol Chem* 63(2):529–545. doi:10.1016/S0021-9258(18)85018-9 [V]
  - Crossref lists all three names and gives the short title "The hemoglobin system".
  - It introduced four-step successive binding: S = (K₁P + 2K₁K₂P² + 3K₁K₂K₃P³ + 4K₁K₂K₃K₄P⁴) / [4(1 + K₁P + K₁K₂P² + K₁K₂K₃P³ + K₁K₂K₃K₄P⁴)] †. — [Source](https://doi.org/10.1016/S0021-9258(18)85018-9)
- **Severinghaus JW (1979).** "Simple, accurate equations for human blood O₂ dissociation computations." *J Appl Physiol* 46(3):599–602. doi:10.1152/jappl.1979.46.3.599 [V]
  - Standard-condition curve: S = [23400/(P³ + 150P) + 1]⁻¹, with P in mmHg.
  - It is corrected to a "virtual" PO₂ by temperature, pH and base excess or PCO₂ factors †. The exact correction coefficients must be taken from the paper. — [Source](https://doi.org/10.1152/jappl.1979.46.3.599)
- **Kelman GR (1966).** "Digital computer subroutine for the conversion of oxygen tension into saturation." *J Appl Physiol* 21(4):1375–1376. doi:10.1152/jappl.1966.21.4.1375 [V]
  - A rational-polynomial fit with temperature, pH and base-excess corrections †. It is a long-standing comparison baseline. — [Source](https://doi.org/10.1152/jappl.1966.21.4.1375)
- **Dash RK, Bassingthwaighte JB (2004).** "Blood HbO₂ and HbCO₂ dissociation curves at varied O₂, CO₂, pH, 2,3-DPG and temperature levels." *Ann Biomed Eng* 32(12):1676–1693. doi:10.1007/s10439-004-7821-6 [V]
  - Their erratum, *Ann Biomed Eng* 38(4):1683–1701 (2010), doi:10.1007/s10439-010-9948-y [V], is a full 19-page corrected version. Implement from the 2010 text. — [Source](https://doi.org/10.1007/s10439-004-7821-6); [Erratum](https://doi.org/10.1007/s10439-010-9948-y)
- **Dash RK, Korman B, Bassingthwaighte JB (2016).** "Simple accurate mathematical models of blood HbO₂ and HbCO₂ dissociation curves at varied physiological conditions: evaluation and comparison with other models." *Eur J Appl Physiol* 116(1):97–113. doi:10.1007/s00421-015-3228-3 [V]
  - It compares Hill, Adair, Kelman, Severinghaus and their own models †. It is the natural cross-model fixture paper. — [Source](https://doi.org/10.1007/s00421-015-3228-3)

#### 1.2 The Krogh cylinder (classical)
- **Krogh A (1919).** "The number and distribution of capillaries in muscles with calculations of the oxygen pressure head necessary for supplying the tissue." *J Physiol* 52(6):409–415. doi:10.1113/jphysiol.1919.sp001839 [V]
  - A capillary of radius r_c sits in a coaxial tissue cylinder of radius R, with uniform consumption M and Krogh coefficient K = Dα. Steady radial diffusion gives the Krogh–Erlang profile, a closed-form fixture (my derivation, consistent with the classical result †):
    P(r) = P_c + (M/4K)(r² − r_c²) − (M R²/2K) ln(r/r_c).
  - The no-flux condition is applied at r = R. — [Source](https://doi.org/10.1113/jphysiol.1919.sp001839)

#### 1.3 Large arteries and the wall (classical)
- **Moore JA, Ethier CR (1997).** "Oxygen mass transfer calculations in large arteries." *J Biomech Eng* 119(4):469–475. doi:10.1115/1.2798295 [V]
  - Arterial O₂ transport with oxyhaemoglobin carriage, entering through the slope of the dissociation curve, and wall-side consumption †. — [Source](https://doi.org/10.1115/1.2798295)
- **Tarbell JM (2003).** "Mass transport in arteries and the localization of atherosclerosis." *Annu Rev Biomed Eng* 5:79–118. doi:10.1146/annurev.bioeng.5.040202.121529 [V]
  - Reviews O₂, albumin, LDL and ATP transport. Covers fluid-phase versus wall-limited transport (Damköhler-type arguments) and the links to wall shear stress (WSS) and atherosclerosis †. — [Source](https://doi.org/10.1146/annurev.bioeng.5.040202.121529)
- **Ethier CR (2002).** "Computational modeling of mass transfer and links to atherosclerosis." *Ann Biomed Eng* 30(4):461–471. doi:10.1114/1.1468890 [V] — [Source](https://doi.org/10.1114/1.1468890)
- **Kaazempur-Mofrad MR, Ethier CR (2001).** "Mass transport in an anatomically realistic human right coronary artery." *Ann Biomed Eng* 29(2):121–127. doi:10.1114/1.1349704 [V] — [Source](https://doi.org/10.1114/1.1349704)
- **Coppola G, Caro C (2008).** "Oxygen mass transfer in a model three-dimensional artery." *J R Soc Interface* 5(26):1067–1075. doi:10.1098/rsif.2007.1338 [V] — [Source](https://doi.org/10.1098/rsif.2007.1338)
- **Tada S (2010).** "Numerical study of oxygen transport in a carotid bifurcation." *Phys Med Biol* 55(14):3993–4010. doi:10.1088/0031-9155/55/14/004 [V] — [Source](https://doi.org/10.1088/0031-9155/55/14/004)

#### 1.4 The microcirculation (classical)
- **Hsu R, Secomb TW (1989).** "A Green's function method for analysis of oxygen delivery to tissue by microvascular networks." *Math Biosci* 96(1):61–78. doi:10.1016/0025-5564(89)90083-7 [V] — [Source](https://doi.org/10.1016/0025-5564(89)90083-7)
- **Secomb TW, Hsu R, Park EYH, Dewhirst MW (2004).** "Green's function methods for analysis of oxygen delivery to tissue by microvascular networks." *Ann Biomed Eng* 32(11):1519–1529. doi:10.1114/B:ABME.0000049036.08817.44 [V]
  - Vessels and tissue are represented by point sources and sinks. The method was applied to realistic 3D networks from skeletal muscle, brain and tumour. — [Source](https://doi.org/10.1114/B:ABME.0000049036.08817.44); [Secomb lab page](https://experts.arizona.edu/en/publications/greens-function-methods-for-analysis-of-oxygen-delivery-to-tissue/)
- **Secomb TW (2016).** "A Green's function method for simulation of time-dependent solute transport and reaction in realistic microvascular geometries." *Math Med Biol* 33(4):475–494. doi:10.1093/imammb/dqv031 [V]
  - The lab's "Green's function method" page distributes the method. The code licence was not confirmed. — [Source](https://doi.org/10.1093/imammb/dqv031); [Lab page](https://sites.arizona.edu/secomb/greens-function-method)
- **Pries AR, Secomb TW, Gaehtgens P, Gross JF (1990).** "Blood flow in microvascular networks. Experiments and simulation." *Circ Res* 67(4):826–834. doi:10.1161/01.RES.67.4.826 [V]
  - Empirical phase-separation (plasma-skimming) laws for haematocrit at bifurcations †. These are needed because O₂ delivery depends on red-cell flux. — [Source](https://doi.org/10.1161/01.RES.67.4.826)
- **Pries AR, Secomb TW (2005).** "Microvascular blood viscosity in vivo and the endothelial surface layer." *Am J Physiol Heart Circ Physiol* 289(6):H2657–H2664. doi:10.1152/ajpheart.00297.2005 [V]
  - The in vivo viscosity law used by most network solvers †. — [Source](https://doi.org/10.1152/ajpheart.00297.2005)
- **Goldman D (2008).** "Theoretical models of microvascular oxygen transport to tissue." *Microcirculation* 15(8):795–811. doi:10.1080/10739680801938289 [V] — review. — [Source](https://doi.org/10.1080/10739680801938289)
- **Secomb TW (2017).** "Blood flow in the microcirculation." *Annu Rev Fluid Mech* 49:443–461. doi:10.1146/annurev-fluid-010816-060302 [V] — [Source](https://doi.org/10.1146/annurev-fluid-010816-060302)
- **Fang Q, Sakadžić S, Ruvinskaya L, Devor A, Dale AM, Boas DA (2008).** "Oxygen advection and diffusion in a three-dimensional vascular anatomical network." *Opt Express* 16(22):17530. doi:10.1364/OE.16.017530 [V] — [Source](https://doi.org/10.1364/OE.16.017530)
- **Linninger AA, Gould IG, Marinnan T, Hsu CY, Chojecki M, Alaraj A (2013).** "Cerebral microcirculation and oxygen tension in the human secondary cortex." *Ann Biomed Eng* 41(11):2264–2284. doi:10.1007/s10439-013-0828-0 [V] — [Source](https://doi.org/10.1007/s10439-013-0828-0)
- **Lücker A, Weber B, Jenny P (2015).** "A dynamic model of oxygen transport from capillaries to tissue with moving red blood cells." *Am J Physiol Heart Circ Physiol* 308(3):H206–H216. doi:10.1152/ajpheart.00447.2014 [V]
  - Uses discrete parachute-shaped red blood cells (RBCs). It reports that homogeneous tissue oxygenation needs capillary networks denser on the venular side than on the arteriolar side. — [Source](https://journals.physiology.org/doi/full/10.1152/ajpheart.00447.2014)
- **Köppl T, Vidotto E, Wohlmuth B (2020).** "A 3D–1D coupled blood flow and oxygen transport model to generate microvascular networks." *Int J Numer Method Biomed Eng* 36(10):e3386. doi:10.1002/cnm.3386 [V] — [Source](https://doi.org/10.1002/cnm.3386)
- **Hartung G, Badr S, Moeini M, Lesage F, Kleinfeld D, Alaraj A, Linninger A (2021).** "Voxelized simulation of cerebral oxygen perfusion elucidates hypoxia in aged mouse cortex." *PLoS Comput Biol* 17(1):e1008584. doi:10.1371/journal.pcbi.1008584 [V] — [Source](https://doi.org/10.1371/journal.pcbi.1008584)
- **Celaya-Alcala JT, Lee GV, Smith AF, Li B, Sakadžić S, Boas DA, Secomb TW (2021).** "Simulation of oxygen transport and estimation of tissue perfusion in extensive microvascular networks: application to cerebral cortex." *J Cereb Blood Flow Metab* 41(3):656–669. doi:10.1177/0271678X20927100 [V]
  - Estimates flow when the boundary conditions are incomplete. Predicted tissue PO₂ is broad, with 99% of tissue between 20 and 80 mmHg under reference conditions. This is a usable distributional fixture. — [Source](https://experts.arizona.edu/en/publications/simulation-of-oxygen-transport-and-estimation-of-tissue-perfusion/)
- **Possenti L, Cicchetti A, Rosati R, Cerroni D, Costantino ML, Rancati T, Zunino P (2021) (bridge).** "A mesoscale computational model for microvascular oxygen transfer." *Ann Biomed Eng* 49(12):3356–3373. doi:10.1007/s10439-021-02807-x [V] — [Source](https://doi.org/10.1007/s10439-021-02807-x)
- **Józsa TI, Padmos RM, Samuels N, El-Bouri WK, Hoekstra AG, Payne SJ (2021) (bridge).** "A porous circulation model of the human brain for in silico clinical trials in ischaemic stroke." *Interface Focus* 11(1):20190127. doi:10.1098/rsfs.2019.0127 [V] — [Source](https://doi.org/10.1098/rsfs.2019.0127)

#### 1.5 Latest (2023–2026)
- **Amiri FA, Zhang J (2023).** "Oxygen transport across tank-treading red blood cell: individual and joint roles of flow convection and oxygen–hemoglobin reaction." *Microvasc Res* 145:104447. doi:10.1016/j.mvr.2022.104447 [V]
  - Lattice-Boltzmann model of a single RBC. Cytoplasmic convection matters at low PO₂ and high shear rate. — [Source](https://www.sciencedirect.com/science/article/abs/pii/S002628622200139X)
- **Zhang J (2025).** "Convection effect of plasma flow on oxygen transport in capillaries: an in-depth numerical investigation." *Microcirculation* 32(3):e70011. doi:10.1111/micc.70011 [V]
  - Couples plasma flow, RBC deformation and O₂–Hb kinetics. Plasma convection's effect on overall O₂ transport can be neglected in typical microcirculatory conditions, which justifies diffusion-only plasma gaps in reduced models. — [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC12009512)
- **Sharma A, Secomb TW (2025).** "A theoretical model for oxygen transport to the cerebral cortex: effects of flow redistribution by penetrating arterioles." *Microvasc Res* 161:104836. doi:10.1016/j.mvr.2025.104836 [V] — [Source](https://doi.org/10.1016/j.mvr.2025.104836)
- **Wang H, Choy JS, Kassab GS, Lee LC (2025).** "Computer model coupling hemodynamics and oxygen transport in the coronary capillary network: pulsatile vs. non-pulsatile analysis." *Comput Methods Programs Biomed* 258:108486. doi:10.1016/j.cmpb.2024.108486 [V] — [Source](https://doi.org/10.1016/j.cmpb.2024.108486)
- **Takeishi N, Kobayashi J, Wada S, Ii S (2026).** "Diffuse interface approach to oxygen transport and metabolism under blood flow dynamics in microcirculations." *Int J Heat Mass Transf* 266:128822. doi:10.1016/j.ijheatmasstransfer.2026.128822 [V]
  - Couples O₂ transport, reaction and consumption with the motion and deformation of individual RBCs across a capillary network.
  - RBCs release more O₂ where local PO₂ is low, which keeps tissue oxygenation uniform even when RBCs are distributed unevenly.
  - The university press release says experimental validation is still planned. — [Source](https://www.eurekalert.org/news-releases/1125999); [Kyushu University](https://www.kyushu-u.ac.jp/en/researches/view/390)
- **Mathur T, Tronolone JJ, Jain A (2023), preprint.** "AngioMT: an in silico platform for digital sensing of oxygen transport through heterogenous microvascular networks." *bioRxiv*. doi:10.1101/2023.01.09.523275 [V as a preprint record]
  - A finite-element O₂ transport program, said in its abstract to be validated. — [Source](https://doi.org/10.1101/2023.01.09.523275)

#### 1.6 Fixtures available
- **Closed forms:**
  - Hill and Adair saturation curves (Hill 1910; Adair 1925).
  - Severinghaus 1979 at standard conditions.
  - The Krogh–Erlang profile (Krogh 1919).
- **Cross-model comparison tables:** Dash et al. 2016 compares the saturation models †. — [Source](https://doi.org/10.1007/s00421-015-3228-3)
- **Distributional network result:** Celaya-Alcala 2021, with 99% of tissue between 20 and 80 mmHg. — [Source](https://experts.arizona.edu/en/publications/simulation-of-oxygen-transport-and-estimation-of-tissue-perfusion/)

### Inferences
- **Large-vessel O₂ equation.** Write total O₂ as αP + [Hb]·C_Hb·S(P). Assuming Hb and plasma O₂ move together, the transport equation for P becomes
  (1 + ([Hb]C_Hb/α)·dS/dP)(∂P/∂t + u·∇P) = D∇²P.
  This is the form used in the Moore & Ethier line †. The factor (1 + …) can exceed 50 near P₅₀. So the effective Péclet number is much larger than for dissolved O₂ alone, and the boundary-layer mesh must allow for it.
- **Schmidt number and boundary-layer thickness.** These values are my assumptions, to be checked against Tarbell 2003: ν_blood ≈ 3.3×10⁻⁶ m²/s and D_O₂ ≈ 2×10⁻⁹ m²/s, giving Sc ≈ 1.7×10³. Lévêque scaling gives δ_c/δ_u ≈ Sc^(−1/3) ≈ 0.08, so roughly ten-micrometre concentration boundary layers in large arteries. Plan prism layers of about 1 µm or a wall-function/Lévêque correction (see §5).
- **Microvascular verification path.**
  1. Krogh closed form.
  2. Single capillary with discrete RBCs (Lücker 2015; Zhang 2025).
  3. A network with the Pries viscosity and phase-separation laws.
  4. Compare the tissue-PO₂ histogram with Celaya-Alcala 2021.
- **For an infrastructure with an exact-arithmetic style.** Hill and Adair are rational functions of P, and Severinghaus is a rational function after the virtual-PO₂ transform. All can be held to closed forms in tests before being coupled to flow.

### Gaps
- Hill 1910's exact title and pages inside the Proceedings were not confirmed (Crossref has only the Proceedings DOI).
- Popel AS (1989), "Theory of oxygen transport to tissue," *Crit Rev Biomed Eng* 17(3):257–321, is [U]: it has no DOI and was not checked.
- The licence of the Secomb lab's Green's-function code, and of the Linninger and Köppl codes, was not found.
- No 2023–2026 paper on O₂ transport in large arteries or the arterial wall surfaced. Crossref discovery searches returned only heat-transfer papers. The area seems quiet since about 2010–2016, or is indexed under other terms.
- The Severinghaus virtual-PO₂ correction coefficients are given here only as †. Take them from the paper.

---

## 2. Mass Transport of LDL, Albumin, Nitric Oxide and ATP; Wall Permeability (Kedem–Katchalsky); Relation to Atherosclerosis and WSS

### Takeaway
The standard model has three layers:
- **Lumen.** Navier–Stokes with ADR, at very high Schmidt number: Sc ~ 10⁴–10⁶ for albumin and LDL.
- **Membrane.** The endothelium as a Kedem–Katchalsky membrane, or the WSS-dependent three-pore model.
- **Wall.** A porous wall (Darcy plus ADR with reaction), one or several layers.

Concentration polarisation (LDL accumulating at the luminal surface because water filters through the wall) is the central coupled effect. Models span from wall-free (Wada & Karino 2002; Rappitsch & Perktold 1996), through coupled lumen–wall (Stangeby & Ethier 2002; Prosi 2005; Olgac 2008; Yang & Vafai 2006), to multiscale with glycocalyx (Vincent & Weinberg 2014). The latest work couples LDL transport to cellular plaque growth on patient-specific coronaries (Corti 2026).

For nitric oxide, the reaction is the hard part: very fast scavenging by haemoglobin. The literature is Lancaster 1994 → Vaughn 1998 → Buerk 2001 → Chen 2006, and flow-chamber models (Fadel 2009; Plata 2010).

### Cited Findings

#### 2.1 Membrane transport (classical)
- **Kedem O, Katchalsky A (1958).** "Thermodynamic analysis of the permeability of biological membranes to non-electrolytes." *Biochim Biophys Acta* 27:229–246. doi:10.1016/0006-3002(58)90330-5 [V]
  - The Kedem–Katchalsky (K–K) relations †:
    J_v = L_p(Δp − σ_d·Δπ)
    J_s = P·Δc + (1 − σ_f)·J_v·c̄
  - These are the membrane conditions used at the endothelium and the internal elastic lamina (IEL). — [Source](https://doi.org/10.1016/0006-3002(58)90330-5)

#### 2.2 Wall-free and coupled lumen–wall LDL and albumin models (classical)
- **Rappitsch G, Perktold K (1996).** "Pulsatile albumin transport in large arteries: a numerical simulation study." *J Biomech Eng* 118(4):511–519. doi:10.1115/1.2796038 [V] — [Source](https://doi.org/10.1115/1.2796038)
- **Stangeby DK, Ethier CR (2002).** "Computational analysis of coupled blood-wall arterial LDL transport." *J Biomech Eng* 124(1):1–8. doi:10.1115/1.1427041 [V] — [Source](https://doi.org/10.1115/1.1427041)
- **Wada S, Karino T (2002).** "Theoretical prediction of low-density lipoproteins concentration at the luminal surface of an artery with a multiple bend." *Ann Biomed Eng* 30(6):778–791. doi:10.1114/1.1495868 [V] — [Source](https://doi.org/10.1114/1.1495868)
- **Prosi M, Zunino P, Perktold K, Quarteroni A (2005).** "Mathematical and numerical models for transfer of low-density lipoproteins through the arterial walls: a new methodology for the model set up with applications to the study of disturbed lumenal flow." *J Biomech* 38(4):903–917. doi:10.1016/j.jbiomech.2004.04.024 [V]
  - A hierarchy of wall models (wall-free, fluid-wall, multilayer), with domain-decomposition coupling †. — [Source](https://doi.org/10.1016/j.jbiomech.2004.04.024)
- **Kaazempur-Mofrad MR, Wada S, Myers JG, Ethier CR (2005).** "Mass transport and fluid flow in stenotic arteries: axisymmetric and asymmetric models." *Int J Heat Mass Transf* 48(21–22):4510–4517. doi:10.1016/j.ijheatmasstransfer.2005.05.004 [V] — [Source](https://doi.org/10.1016/j.ijheatmasstransfer.2005.05.004)
- **Yang N, Vafai K (2006).** "Modeling of low-density lipoprotein (LDL) transport in the artery—effects of hypertension." *Int J Heat Mass Transf* 49(5–6):850–867. doi:10.1016/j.ijheatmasstransfer.2005.09.019 [V]
  - The four-layer wall: endothelium, intima, IEL and media, each with porous-medium parameters †. Its parameter table is widely reused. — [Source](https://doi.org/10.1016/j.ijheatmasstransfer.2005.09.019)
- **Ai L, Vafai K (2006).** "A coupling model for macromolecule transport in a stenosed arterial wall." *Int J Heat Mass Transf* 49(9–10):1568–1591. doi:10.1016/j.ijheatmasstransfer.2005.10.041 [V] — [Source](https://doi.org/10.1016/j.ijheatmasstransfer.2005.10.041)
- **Olgac U, Kurtcuoglu V, Poulikakos D (2008).** "Computational modeling of coupled blood-wall mass transport of LDL: effects of local wall shear stress." *Am J Physiol Heart Circ Physiol* 294(2):H909–H919. doi:10.1152/ajpheart.01082.2007 [V]
  - The three-pore endothelium: vesicular, normal-junction and leaky-junction pathways. The fraction of leaky cells depends on WSS †. — [Source](https://doi.org/10.1152/ajpheart.01082.2007)
  - Its 2009 patient-specific follow-up on a CT-derived left coronary artery found the diseased-state high-LDL region displaced distal to the plaque, consistent with downstream plaque growth. That follow-up is [W]: PubMed 19329764, full bibliographic record not checked. — [Source](https://pubmed.ncbi.nlm.nih.gov/19329764/)
- **Vincent PE, Sherwin SJ, Weinberg PD (2009).** "The effect of a spatially heterogeneous transmural water flux on concentration polarization of low density lipoprotein in arteries." *Biophys J* 96(8):3102–3115. doi:10.1016/j.bpj.2009.01.022 [V] — [Source](https://doi.org/10.1016/j.bpj.2009.01.022)
- **Vincent P, Sherwin S, Weinberg P (2010).** "The effect of the endothelial glycocalyx layer on concentration polarisation of low density lipoprotein in arteries." *J Theor Biol* 265(1):1–17. doi:10.1016/j.jtbi.2010.04.015 [V] — [Source](https://doi.org/10.1016/j.jtbi.2010.04.015)
- **Vincent PE, Weinberg PD (2014).** "Flow-dependent concentration polarization and the endothelial glycocalyx layer: multi-scale aspects of arterial mass transport and their implications for atherosclerosis." *Biomech Model Mechanobiol* 13(2):313–326. doi:10.1007/s10237-013-0512-1 [V]
  - The DOI ends in **-1**. The often-quoted "-3" does not resolve. — [Source](https://doi.org/10.1007/s10237-013-0512-1)

#### 2.3 High-Schmidt near-wall transport (classical)
- **Hansen KB, Shadden SC (2016).** "A reduced-dimensional model for near-wall transport in cardiovascular flows." *Biomech Model Mechanobiol* 15(3):713–722. doi:10.1007/s10237-015-0719-4 [V] — [Source](https://doi.org/10.1007/s10237-015-0719-4)
- **Arzani A, Gambaruto AM, Chen G, Shadden SC (2016).** "Lagrangian wall shear stress structures and near-wall transport in high-Schmidt-number aneurysmal flows." *J Fluid Mech* 790:158–172. doi:10.1017/jfm.2016.6 [V] — [Source](https://doi.org/10.1017/jfm.2016.6)
- **Hansen KB, Arzani A, Shadden SC (2019).** "Finite element modeling of near-wall mass transport in cardiovascular flows." *Int J Numer Method Biomed Eng* 35(1):e3148. doi:10.1002/cnm.3148 [V] — [Source](https://doi.org/10.1002/cnm.3148)

#### 2.4 Nitric oxide and ATP (classical)
- **Lancaster JR (1994).** "Simulation of the diffusion and reaction of endogenously produced nitric oxide." *PNAS* 91(17):8137–8141. doi:10.1073/pnas.91.17.8137 [V] — [Source](https://doi.org/10.1073/pnas.91.17.8137)
- **Vaughn MW, Kuo L, Liao JC (1998).** "Estimation of nitric oxide production and reaction rates in tissue by use of a mathematical model." *Am J Physiol* 274(6):H2163–H2176. doi:10.1152/ajpheart.1998.274.6.H2163 [V]
  - Arteriole–wall–tissue geometry with an RBC-free layer and Hb scavenging †. — [Source](https://doi.org/10.1152/ajpheart.1998.274.6.H2163)
- **Buerk DG (2001).** "Can we model nitric oxide biotransport? A survey of mathematical models for a simple diatomic molecule with surprisingly complex biological activities." *Annu Rev Biomed Eng* 3:109–143. doi:10.1146/annurev.bioeng.3.1.109 [V] — [Source](https://doi.org/10.1146/annurev.bioeng.3.1.109)
- **Chen X, Jaron D, Barbee KA, Buerk DG (2006).** "The influence of radial RBC distribution, blood velocity profiles, and glycocalyx on coupled NO/O₂ transport." *J Appl Physiol* 100(2):482–492. doi:10.1152/japplphysiol.00633.2005 [V] — [Source](https://doi.org/10.1152/japplphysiol.00633.2005)
- **Fadel AA, Barbee KA, Jaron D (2009).** "A computational model of nitric oxide production and transport in a parallel plate flow chamber." *Ann Biomed Eng* 37(5):943–954. doi:10.1007/s10439-009-9658-5 [V] — [Source](https://doi.org/10.1007/s10439-009-9658-5)
- **Plata AM, Sherwin SJ, Krams R (2010).** "Endothelial nitric oxide production and transport in flow chambers: the importance of convection." *Ann Biomed Eng* 38(9):2805–2816. doi:10.1007/s10439-010-0039-x [V] — [Source](https://doi.org/10.1007/s10439-010-0039-x)
- **John K, Barakat AI (2001).** "Modulation of ATP/ADP concentration at the endothelial surface by shear stress: effect of flow-induced ATP release." *Ann Biomed Eng* 29(9):740–751. doi:10.1114/1.1397792 [V]
  - Shear-dependent release and ectonucleotidase degradation, applied as a reactive wall boundary condition †. — [Source](https://doi.org/10.1114/1.1397792)

#### 2.5 Latest (2023–2026)
- **Corti A, De Nisco G, Wentzel JJ, Migliavacca F, Morbiducci U, Chiastra C (2026).** "Personalized multiscale modeling of coronary plaque progression: the interaction between low-density-lipoprotein transport and cellular dynamics." *Comput Methods Appl Mech Eng* 448:118427. doi:10.1016/j.cma.2025.118427 [V]
  - Combines lumen haemodynamics, a three-pore LDL endothelium and agent-based cellular dynamics on a patient-specific coronary.
  - Calibrated against 1-year follow-up, it predicted the 1-year change in wall area in two cross-sections with less than 10% error.
  - Blood LDL level changed plaque morphology, from localised to diffuse. — [Source](https://pure.eur.nl/en/publications/personalized-multiscale-modeling-of-coronary-plaque-progression-t/)
- **Fan Z, Yao J, Xu J, Liu X, Liu M, Ye X, Deng X (2024).** "Nitric oxide transport in carotid bifurcation after different stent interventions: a numerical study." *Math Model Nat Phenom* 19:7. doi:10.1051/mmnp/2023039 [V] — [Source](https://doi.org/10.1051/mmnp/2023039)

#### 2.6 Fixtures available
- **K–K and three-pore algebraic relations.** Unit-testable from Kedem & Katchalsky 1958 and Olgac 2008.
- **Four-layer parameter table.** Yang & Vafai 2006 † gives a 1D steady transmural profile with a closed form per layer (an inference: piecewise exponential solutions of the 1D ADR with first-order reaction).
- **Concentration polarisation in a straight tube with wall filtration.** Wada & Karino 2002 † and Stangeby & Ethier 2002 † are standard benchmark geometries.

### Inferences
- **Schmidt numbers (my assumptions; check against Tarbell 2003).**

  | Species | D (m²/s) | Sc |
  |---|---|---|
  | Albumin | ≈ 7×10⁻¹¹ | ≈ 5×10⁴ |
  | LDL | ≈ 5×10⁻¹² | ≈ 7×10⁵ |

  At Sc ≈ 7×10⁵, δ_c/δ_u ≈ 0.011. LDL boundary layers are therefore sub-micrometre to a few micrometres. This is why Hansen & Shadden developed reduced near-wall models, and why wall-free models replace the wall with K–K Robin conditions.
- **Order of implementation.**
  1. K–K Robin boundary condition on a fixed lumen, as in Wada & Karino.
  2. Coupled single-layer porous wall, as in Stangeby & Ethier.
  3. Multilayer wall, as in Yang & Vafai or Prosi.
  4. WSS-dependent three-pore model, as in Olgac.
  5. Coupling to plaque growth, as in Corti 2026.
- **NO is reaction-limited, not convection-limited.** Hb scavenging is very fast. The reaction term should be integrated implicitly or with exact local exponential updates; see §5.

### Gaps
- Few 2023–2026 papers on coupled LDL or albumin transport surfaced. Crossref discovery searches were noisy, and the Corti 2026 introduction itself says patient-specific coupled growth models are scarce.
- No public, licensed reference dataset for LDL wall concentration was found. Validation is usually qualitative, against plaque sites or experimental concentration-polarisation measurements.
- The 2009 Olgac patient-specific paper's full citation was not verified.

---

## 3. Coagulation and Thrombosis: Reaction Networks, Flow-Coupled Clotting, Platelets, Venous Thrombosis, Validation Data, Model Reduction and Software

### Takeaway
**The well-mixed lineage** runs Lawson/Jones & Mann 1994 → Hockin 2002 (34 ODEs, 42 rate constants, extrinsic pathway with TFPI and AT) → Butenas 2004 and Chatterjee 2010 (contact pathway) → Mitrophanov 2014 (fibrin). It is curated in BioModels (CC0) with reproduced figures, so these models are the best ODE fixtures.

**Flow-coupled lineage.**
- Kuharsky & Fogelson 2001: a well-mixed zone with mass-transfer coefficients.
- Fogelson & Tania 2005.
- Leiderman & Fogelson 2011: spatial, with a Brinkman resistance from platelet aggregates.
- Link 2018/2020 and Miyazawa 2023 I–II.
- Open-source clotFoam (2023).

**Macroscopic and device lineage.** Sorensen 1999; Wu 2017; Taylor 2016; Méndez Rojano 2022; Zhussupbekov 2022.

**Particle and multiscale lineage.** Xu 2008; Flamm 2012; Lu 2016; Shankar 2022/2023; Yazdani 2017.

**Validation.** Calibrated automated thrombography (CAT; Hemker 2003) and synthetic-plasma thrombin generation (Butenas 1999). The 2023 critical evaluation (Ranc 2023) shows the Chatterjee and Butenas schemes need modifications to fit CAT factor-deficiency data, and publishes its data and code.

**Reductions.** Hansen & Shadden 2019 automated reduction (9- and 7-species models in BioModels); Papadopoulos 2014.

**Venous thrombosis.** Thrombosis under low shear and stasis is modelled by Govindarajan 2016 and Bouchnita 2020 (the latter in BioModels), and reviewed in Panteleev 2024. Its physiology is anchored in Virchow's triad (Bagot & Arya 2008) and in valve-pocket stasis and hypoxia (Karino & Motomiya 1984; Bovill & van der Vliet 2011). Patient-specific DVT CFD has begun: Jimoh-Taiwo 2022 (femoral vein, thrombin–fibrinogen) and George & Ngoepe 2026 (valve aperture).

### Cited Findings

#### 3.1 Well-mixed reaction networks (classical)
- **Lawson JH, Kalafatis M, Stram S, Mann KG (1994).** "A model for the tissue factor pathway to thrombin. I. An empirical study." *J Biol Chem* 269(37):23357–23366. doi:10.1016/S0021-9258(17)31661-7 [V]
  - The experimental data behind the first ODE model †. — [Source](https://doi.org/10.1016/S0021-9258(17)31661-7)
- **Jones KC, Mann KG (1994).** "A model for the tissue factor pathway to thrombin. II. A mathematical simulation." *J Biol Chem* 269(37):23367–23373. doi:10.1016/S0021-9258(17)31662-9 [V]
  - BioModels: **BIOMD0000000336** (Jones1994_BloodCoagulation), curated. — [Source](https://doi.org/10.1016/S0021-9258(17)31662-9); [BioModels search](https://biomodels.org/search?query=coagulation&numResults=100&format=json)
- **Hockin MF, Jones KC, Everse SJ, Mann KG (2002).** "A model for the stoichiometric regulation of blood coagulation." *J Biol Chem* 277(21):18322–18333. doi:10.1074/jbc.M201173200 [V]
  - 34 ODEs and 42 rate constants. The extrinsic system with the TFPI and antithrombin-III anticoagulants.
  - BioModels **BIOMD0000000335** (Hockin2002_BloodCoagulation), curated by reproducing the paper's figures.
  - Also run on BioSimulations in COPASI, and simulated for Fig. 3 on JWS Online. These are ready reference trajectories. — [Source](https://www.omicsdi.org/dataset/biomodels/BIOMD0000000335); [BioSimulations](https://biosimulations.org/projects/BIOMD0000000335_copasi_A_model_for_the_stoichiometric); [JWS Online Fig. 3](https://jjj.mib.ac.uk/models/experiments/hockin2002_fig3-user)
- **Butenas S, van't Veer C, Mann KG (1999).** "'Normal' thrombin generation." *Blood* 94(7):2169–2178. doi:10.1182/blood.V94.7.2169.419k22_2169_2178 [V]
  - Synthetic-plasma thrombin generation at mean plasma factor levels †. It is the usual empirical comparison for Hockin-type models. — [Source](https://doi.org/10.1182/blood.V94.7.2169.419k22_2169_2178)
- **Butenas S, Orfeo T, Gissel MT, Brummel KE, Mann KG (2004).** "The significance of circulating factor IXa in blood." *J Biol Chem* 279(22):22875–22882. doi:10.1074/jbc.M400531200 [V]
  - Ranc 2023 describes the Butenas 2004 JBC scheme as Hockin plus FX activation by IXa.
  - BioModels lists **BIOMD0000000362** (Butenas2004_BloodCoagulation). That this entry is this paper is my inference. — [Source](https://doi.org/10.1074/jbc.M400531200); [Ranc 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10461854/)
- **Wajima T, Isbister GK, Duffull SB (2009).** "A comprehensive model for the humoral coagulation network in humans." *Clin Pharmacol Ther* 86(3):290–298. doi:10.1038/clpt.2009.87 [V]
  - BioModels **BIOMD0000000338** (aPTT test), **BIOMD0000000339** (PT test) and **BIOMD0000000340** (warfarin/heparin). — [Source](https://doi.org/10.1038/clpt.2009.87); [BioModels](https://biomodels.org/search?query=coagulation&numResults=100&format=json)
- **Chatterjee MS, Denney WS, Jing H, Diamond SL (2010).** "Systems biology of coagulation initiation: kinetics of thrombin generation in resting and activated human blood." *PLoS Comput Biol* 6(9):e1000950. doi:10.1371/journal.pcbi.1000950 [V]
  - Extends Hockin with contact-pathway initiation.
  - BioModels **MODEL1108260014** (non-curated). Its OmicsDI record still carries "copyright (c) 2005–2011 The BioModels.net Team". — [Source](https://doi.org/10.1371/journal.pcbi.1000950); [OmicsDI](https://www.omicsdi.org/dataset/biomodels/MODEL1108260014)
- **Danforth CM, Orfeo T, Mann KG, Brummel-Ziedins KE, Everse SJ (2009).** "The impact of uncertainty in a blood coagulation model." *Math Med Biol* 26(4):323–336. doi:10.1093/imammb/dqp011 [V]
  - Propagates literature ranges of the Hockin rate constants by sampling †. — [Source](https://doi.org/10.1093/imammb/dqp011)
- **Mitrophanov AY, Wolberg AS, Reifman J (2014).** "Kinetic model facilitates analysis of fibrin generation and its modulation by clotting factors: implications for hemostasis-enhancing therapies." *Mol BioSyst* 10(9):2347–2357. doi:10.1039/C4MB00263F [V]
  - BioModels **MODEL1806260001** ("Extended Hockin … with additional fibrin-associated species"). Related entries: **BIOMD0000000951** (Mitrophanov2015, pH) and MODEL1806040001 (2011, rFVIIa). — [Source](https://doi.org/10.1039/C4MB00263F); [BioModels](https://biomodels.org/search?query=coagulation&numResults=100&format=json)
- **Anand M, Rajagopal K, Rajagopal KR (2008).** "A model for the formation, growth, and lysis of clots in quiescent plasma. A comparison between the effects of antithrombin III deficiency and protein C deficiency." *J Theor Biol* 253(4):725–738. doi:10.1016/j.jtbi.2008.04.015 [V] — [Source](https://doi.org/10.1016/j.jtbi.2008.04.015)
- **Mann KG (2012) and Hemker HC, Kerdelo S, Kremers RMW (2012): the debate on model value.**
  - Mann: "Is there value in kinetic modeling of thrombin generation? Yes." *J Thromb Haemost* 10(8):1463–1469. doi:10.1111/j.1538-7836.2012.04799.x [V]
  - Hemker et al.: "Is there value in kinetic modeling of thrombin generation? No (unless…)." *J Thromb Haemost* 10(8):1470–1477. doi:10.1111/j.1538-7836.2012.04802.x [V]
  - Read both before claiming a model is "validated". — [Mann](https://doi.org/10.1111/j.1538-7836.2012.04799.x); [Hemker](https://doi.org/10.1111/j.1538-7836.2012.04802.x)
- **Other BioModels coagulation entries, useful as fixtures.** The search returned 66 models in all. All IDs are from the BioModels search. — [Source](https://biomodels.org/search?query=coagulation&numResults=100&format=json)
  - Bungay2003 (BIOMD0000000334)
  - Hockin1999 Va inactivation (BIOMD0000000365)
  - Panteleev2010 full model (BIOMD0000000740)
  - Dunster2016 nondimensional (BIOMD0000000925)
  - Bravo2012 (BIOMD0000000739)
  - Nayak2015 (BIOMD0000000611)
  - Nagashima2002 (BIOMD0000000747)
  - Pearce2021 fibrin polymerisation (BIOMD0000001054)
  - Kuharsky2001 (MODEL1109150000)
  - Anand2003 (MODEL1806130003)
  - Bouchnita2020 venous flow (MODEL2011030003)
  - Rojano2019 (MODEL1907120001)

#### 3.2 Licences of the BioModels entries
- BioModels' Terms of Use and FAQ state that all models, in every format, and their annotations are provided under the **Creative Commons CC0** public-domain dedication. Curated models pass a curation pipeline. Non-curated models are also CC0. — [BioModels Terms of Use](https://www.ebi.ac.uk/biomodels/termsofuse); [FAQ](https://www.ebi.ac.uk/biomodels/faq); [Wikipedia summary](https://en.wikipedia.org/wiki/BioModels)
- A 2022 reproducibility mailing-list thread asked whether models calibrated on more restrictively licensed data can truly be CC0. It was unresolved. — [Thread](https://lists.simtk.org/pipermail/vp-reproduce-subgroup/2022-January/000308.html)
- **Conflict.** Older OmicsDI mirror records (e.g., Chatterjee2010) show a "copyright (c) 2005–2011 BioModels.net Team" notice, while the live terms say CC0. Treat the live terms as authoritative and record the retrieval date. — [OmicsDI](https://www.omicsdi.org/dataset/biomodels/MODEL1108260014)

#### 3.3 Flow-coupled coagulation and platelet models (classical)
- **Kuharsky AL, Fogelson AL (2001).** "Surface-mediated control of blood coagulation: the role of binding site densities and platelet deposition." *Biophys J* 80(3):1050–1074. doi:10.1016/S0006-3495(01)76085-7 [V]
  - Coagulation and platelet deposition in a well-mixed reaction zone over a small injury. Flow enters through mass-transfer coefficients, and reactions occur on platelet surfaces with explicit binding-site densities †. Parameter tables are given in the paper †. — [Source](https://doi.org/10.1016/S0006-3495(01)76085-7)
- **Fogelson AL, Tania N (2005).** "Coagulation under flow: the influence of flow-mediated transport on the initiation and inhibition of coagulation." *Pathophysiol Haemost Thromb* 34(2–3):91–108. doi:10.1159/000089930 [V] — [Source](https://doi.org/10.1159/000089930)
- **Anand M, Rajagopal K, Rajagopal KR (2003).** "A model incorporating some of the mechanical and biochemical factors underlying clot formation and dissolution in flowing blood." *J Theor Med* 5(3–4):183–218. doi:10.1080/10273660412331317415 [V]
  - The journal now indexes as *Comput Math Methods Med*.
  - An intrinsic-pathway ADR network in a shear-thinning fluid †. BioModels MODEL1806130003. — [Source](https://doi.org/10.1080/10273660412331317415)
- **Panteleev MA, Ovanesov MV, Kireev DA, et al. (2006).** "Spatial propagation and localization of blood coagulation are regulated by intrinsic and protein C pathways, respectively." *Biophys J* 90(5):1489–1500. doi:10.1529/biophysj.105.069062 [V]
  - Nine authors. A reaction–diffusion model †. — [Source](https://doi.org/10.1529/biophysj.105.069062)
- **Bodnár T, Sequeira A (2008).** "Numerical simulation of the coagulation dynamics of blood." *Comput Math Methods Med* 9(2):83–104. doi:10.1080/17486700701852784 [V] — [Source](https://doi.org/10.1080/17486700701852784)
- **Leiderman K, Fogelson AL (2011).** "Grow with the flow: a spatial-temporal model of platelet deposition and blood coagulation under flow." *Math Med Biol* 28(1):47–84. doi:10.1093/imammb/dqq005 [V]
  - Two-dimensional Navier–Stokes with a Brinkman (porosity) resistance from bound platelets. Coagulation proteins are advected, diffused and react in the fluid and on platelet and subendothelial surfaces †. — [Source](https://doi.org/10.1093/imammb/dqq005)
- **Fogelson AL, Neeves KB (2015).** "Fluid mechanics of blood clot formation." *Annu Rev Fluid Mech* 47:377–403. doi:10.1146/annurev-fluid-010814-014513 [V] — the field's standard review. — [Source](https://doi.org/10.1146/annurev-fluid-010814-014513)
- **Link KG, Stobb MT, Di Paola J, Neeves KB, Fogelson AL, Sindi SS, Leiderman K (2018).** "A local and global sensitivity analysis of a mathematical model of coagulation and platelet deposition under flow." *PLoS ONE* 13(7):e0200917. doi:10.1371/journal.pone.0200917 [V] — [Source](https://doi.org/10.1371/journal.pone.0200917)
- **Link KG, Stobb MT, Sorrells MG, Bortot M, Ruegg K, Manco-Johnson MJ, et al. (2020).** "A mathematical model of coagulation under flow identifies factor V as a modifier of thrombin generation in hemophilia A." *J Thromb Haemost* 18(2):306–317. doi:10.1111/jth.14653 [V] — 11 authors. — [Source](https://doi.org/10.1111/jth.14653)
- **Du J, Fogelson AL (2018).** "A two-phase mixture model of platelet aggregation." *Math Med Biol* 35(2):225–256. doi:10.1093/imammb/dqx001 [V] — [Source](https://doi.org/10.1093/imammb/dqx001)
- **Sorensen EN, Burgreen GW, Wagner WR, Antaki JF (1999).** "Computational simulation of platelet deposition and activation: I. Model development and properties." *Ann Biomed Eng* 27(4):436–448. doi:10.1114/1.200 [V]
  - Convection–diffusion–reaction of resting and activated platelets and agonists (ADP, TxA₂, thrombin) †. — [Source](https://doi.org/10.1114/1.200)
- **Wu WT, Jamiolkowski MA, Wagner WR, Aubry N, Massoudi M, Antaki JF (2017).** "Multi-constituent simulation of thrombus deposition." *Sci Rep* 7:42720. doi:10.1038/srep42720 [V] — [Source](https://doi.org/10.1038/srep42720)
- **Taylor JO, Meyer RS, Deutsch S, Manning KB (2016).** "Development of a computational model for macroscopic predictions of device-induced thrombosis." *Biomech Model Mechanobiol* 15(6):1713–1731. doi:10.1007/s10237-016-0793-2 [V] — [Source](https://doi.org/10.1007/s10237-016-0793-2)
- **Menichini C, Xu XY (2016).** "Mathematical modeling of thrombus formation in idealized models of aortic dissection: initial findings and potential applications." *J Math Biol* 73(5):1205–1226. doi:10.1007/s00285-016-0986-4 [V] — [Source](https://doi.org/10.1007/s00285-016-0986-4)
- **Yazdani A, Li H, Humphrey JD, Karniadakis GE (2017).** "A general shear-dependent model for thrombus formation." *PLoS Comput Biol* 13(1):e1005291. doi:10.1371/journal.pcbi.1005291 [V] — [Source](https://doi.org/10.1371/journal.pcbi.1005291)
- **Xu Z, Chen N, Kamocka MM, Rosen ED, Alber M (2008).** "A multiscale model of thrombus development." *J R Soc Interface* 5(24):705–722. doi:10.1098/rsif.2007.1202 [V] — [Source](https://doi.org/10.1098/rsif.2007.1202)
- **Flamm MH, Colace TV, Chatterjee MS, et al. (2012).** "Multiscale prediction of patient-specific platelet function under flow." *Blood* 120(1):190–198. doi:10.1182/blood-2011-10-388140 [V] — nine authors. — [Source](https://doi.org/10.1182/blood-2011-10-388140)
- **Lu Y, Lee MY, Zhu S, Sinno T, Diamond SL (2016).** "Multiscale simulation of thrombus growth and vessel occlusion triggered by collagen/tissue factor using a data-driven model of combinatorial platelet signalling." *Math Med Biol*, advance access dqw015. doi:10.1093/imammb/dqw015 [V]
  - The final issue and pages are [U]. — [Source](https://doi.org/10.1093/imammb/dqw015)
- **Bannish BE, Keener JP, Fogelson AL (2014).** "Modelling fibrinolysis: a 3D stochastic multiscale model." *Math Med Biol* 31(1):17–44. doi:10.1093/imammb/dqs029 [V] — [Source](https://doi.org/10.1093/imammb/dqs029)
- **Belyaev AV, Panteleev MA, Ataullakhanov FI (2015).** "Threshold of microvascular occlusion: injury size defines the thrombosis scenario." *Biophys J* 109(2):450–456. doi:10.1016/j.bpj.2015.06.019 [V] — [Source](https://doi.org/10.1016/j.bpj.2015.06.019)
- **Zhussupbekov M, Méndez Rojano R, Wu WT, Massoudi M, Antaki JF (2021).** "A continuum model for the unfolding of von Willebrand factor." *Ann Biomed Eng* 49(9):2646–2658. doi:10.1007/s10439-021-02845-5 [V] — [Source](https://doi.org/10.1007/s10439-021-02845-5)
- **Zhussupbekov M, Méndez Rojano R, Wu WT, Antaki JF (2022) (bridge).** "von Willebrand factor unfolding mediates platelet deposition in a model of high-shear thrombosis." *Biophys J* 121(21):4033–4047. doi:10.1016/j.bpj.2022.09.040 [V] — [Source](https://doi.org/10.1016/j.bpj.2022.09.040)
- **Méndez Rojano R, Lai A, Zhussupbekov M, Burgreen GW, Cook K, Antaki JF (2022) (bridge).** "A fibrin enhanced thrombosis model for medical devices operating at low shear regimes or large surface areas." *PLoS Comput Biol* 18(10):e1010277. doi:10.1371/journal.pcbi.1010277 [V] — [Source](https://doi.org/10.1371/journal.pcbi.1010277)
- **Shankar KN, Zhang Y, Sinno T, Diamond SL (2022) (bridge).** "A three-dimensional multiscale model for the prediction of thrombus growth under flow with single-platelet resolution." *PLoS Comput Biol* 18(1):e1009850. doi:10.1371/journal.pcbi.1009850 [V]
  - Embeds a reduced coagulation cascade in a 3D flow framework. — [Source](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8827456/)

#### 3.4 Venous thrombosis (low shear, stasis)
- **Govindarajan V, Rakesh V, Reifman J, Mitrophanov AY (2016).** "Computational study of thrombus formation and clotting factor effects under venous flow conditions." *Biophys J* 110(8):1869–1885. doi:10.1016/j.bpj.2016.03.010 [V]
  - Couples a Hockin-type cascade to venous-shear flow †. — [Source](https://www.sciencedirect.com/science/article/pii/S0006349516300595)
- **Bouchnita A, Terekhov K, Nony P, Vassilevski Y, Volpert V (2020).** "A mathematical model to quantify the effects of platelet count, shear rate, and injury size on the initiation of blood coagulation under venous flow conditions." *PLoS ONE* 15(7):e0235392. doi:10.1371/journal.pone.0235392 [V]
  - Available in BioModels as **MODEL2011030003** ("Bouchnita2020 – Model for blood coagulation under venous flow condition"). — [Source](https://doi.org/10.1371/journal.pone.0235392); [BioModels](https://biomodels.org/search?query=coagulation&numResults=100&format=json)
- **Panteleev MA, Bershadsky ES, Shibeko AM, Nechipurenko DY (2024).** "Current issues in computational modeling of thrombosis, fibrinolysis, and thrombolysis." *Computer Research and Modeling* 16(4):975–995. doi:10.20537/2076-7633-2024-16-4-975-995 [V]
  - Reviews arterial, venous and microcirculatory thrombosis, fibrinolysis and thrombolysis, and the basic model types and their computational complexity. — [Source](https://crm-en.ics.org.ru/journal/article/3504/)
- **Bouchnita A, Yadav K, Llored JP, Gurovich A, Volpert V (2023).** "Thrombin generation thresholds for coagulation initiation under flow." *Axioms* 12(9):873. doi:10.3390/axioms12090873 [V] — [Source](https://doi.org/10.3390/axioms12090873)
- **Virchow's triad: stasis, hypercoagulability and vessel-wall injury.**
  - **Bagot CN, Arya R (2008).** "Virchow and his triad: a question of attribution." *Br J Haematol* 143(2):180–190. doi:10.1111/j.1365-2141.2008.07323.x [V]
  - It traces how the "triad" came to be attributed to Virchow †. Cite this rather than Virchow's 1856 text, which has no DOI. — [Source](https://doi.org/10.1111/j.1365-2141.2008.07323.x)
- **Mechanisms of deep vein thrombosis (DVT).**
  - **Mackman N (2012).** "New insights into the mechanisms of venous thrombosis." *J Clin Invest* 122(7):2331–2336. doi:10.1172/JCI60229 [V]
  - A one-page erratum appeared at 122(9):3368, doi:10.1172/JCI66121 [V]. — [Source](https://doi.org/10.1172/JCI60229)
- **The valve pocket as the nidus.**
  - **Karino T, Motomiya M (1984).** "Flow through a venous valve and its implication for thrombus formation." *Thromb Res* 36(3):245–257. doi:10.1016/0049-3848(84)90224-X [V]
    - Flow visualisation of the vortex trapped in the valve sinus †. — [Source](https://doi.org/10.1016/0049-3848(84)90224-X)
  - **Bovill EG, van der Vliet A (2011).** "Venous valvular stasis–associated hypoxia and thrombosis: what is the link?" *Annu Rev Physiol* 73:527–545. doi:10.1146/annurev-physiol-012110-142305 [V]
    - Stasis in the valve pocket leads to local haemoglobin desaturation and hypoxia of the endothelium, and hence to a procoagulant state †.
    - This couples the oxygen module (§1) directly to thrombosis. — [Source](https://doi.org/10.1146/annurev-physiol-012110-142305)
  - **Lurie F, Kistner RL, Eklof B, Kessler D (2003).** "Mechanism of venous valve closure and role of the valve in circulation: a new concept." *J Vasc Surg* 38(5):955–961. doi:10.1016/S0741-5214(03)00711-0 [V]
    - Its 2002 precursor: *J Vasc Surg* 35(4):713–717, doi:10.1067/mva.2002.121123 [V]. These describe the valve-cycle phases that a moving-valve simulation must reproduce †. — [Source](https://doi.org/10.1016/S0741-5214(03)00711-0)
- **Valve and DVT models.**
  - **Ariane M, Wen W, Vigolo D, Brill A, Nash FGB, Barigou M, Alexiadis A (2017).** "Modelling and simulation of flow and agglomeration in deep veins valves using discrete multi physics." *Comput Biol Med* 89:96–103. doi:10.1016/j.compbiomed.2017.07.020 [V]
    - A mesh-free particle model of elastic valve leaflets and agglomerating particles. The follow-up is Baksamawi et al. (2021), *ChemEngineering* 5(3):40, doi:10.3390/chemengineering5030040 [V]. — [Source](https://doi.org/10.1016/j.compbiomed.2017.07.020)
  - **Jimoh-Taiwo Q, Haffejee R, Ngoepe M (2022) (bridge).** "A mechano-chemical computational model of deep vein thrombosis." *Front Phys* 10:886193. doi:10.3389/fphy.2022.886193 [V]
    - CFD of a 3D patient-specific common femoral vein: thrombin–fibrinogen reactions, pulsatile flow and clot growth. Per its abstract, it couples the mechanical and biochemical factors that earlier DVT models treated separately. — [Source](https://doi.org/10.3389/fphy.2022.886193)
  - **George S, Ngoepe M (2026).** "A computational fluid dynamics model of venous thrombogenesis: valve aperture effects in deep vein thrombosis." *J Precis Med Health Dis* 7:100055. doi:10.1016/j.premed.2026.100055 [V] — [Source](https://doi.org/10.1016/j.premed.2026.100055)
  - **Ballard MS, Stong C, Biesinger J (2026), preprint.** "A computational investigation of the effect of venous valve morphology on flow conditions conducive to thrombus formation." Research Square, doi:10.21203/rs.3.rs-10681787/v1 [V as a preprint record] — [Source](https://doi.org/10.21203/rs.3.rs-10681787/v1)
  - **Assefa AM, Palaiodimou L, Bourantas G, Sakellarios A (2025).** "Personalized computational fluid dynamics analysis of cerebral venous hemodynamics in a case of deep cerebral vein thrombosis." *J Pers Med* 15(12):570. doi:10.3390/jpm15120570 [V]
    - Patient-specific steady CFD of the venous sinuses from MR venography: outflow redistribution, pressure and WSS. — [Source](https://doi.org/10.3390/jpm15120570)

#### 3.5 Validation data and model evaluation (classical)
- **Hemker HC, Giesen P, Al Dieri R, Regnault V, de Smedt E, Wagenvoord R, Lecompte T, Béguin S (2003).** "Calibrated automated thrombin generation measurement in clotting plasma." *Pathophysiol Haemost Thromb* 33(1):4–15. doi:10.1159/000071636 [V]
  - The CAT method. Its readouts (lag time, peak, time-to-peak, endogenous thrombin potential) are the usual comparison quantities †. — [Source](https://doi.org/10.1159/000071636)
- **Chelle P, Morin C, Montmartin A, Piot M, Cournil M, Tardy-Poncet B (2018).** "Evaluation and calibration of in silico models of thrombin generation using experimental data from healthy and haemophilic subjects." *Bull Math Biol* 80(8):1989–2025. doi:10.1007/s11538-018-0440-4 [V]
  - The correct DOI ends **0440-4**. — [Source](https://doi.org/10.1007/s11538-018-0440-4)
- **Ratto N, Bouchnita A, Chelle P, Marion M, Panteleev M, Nechipurenko D, Tardy-Poncet B, Volpert V (2021).** "Patient-specific modelling of blood coagulation." *Bull Math Biol* 83(5):50. doi:10.1007/s11538-021-00890-8 [V] — [Source](https://doi.org/10.1007/s11538-021-00890-8)

#### 3.6 Model reduction
- **Hansen KB, Shadden SC (2019).** "Automated reduction of blood coagulation models." *Int J Numer Method Biomed Eng* 35(10):e3220. doi:10.1002/cnm.3220 [V]
  - BioModels **BIOMD0000000755** (nine-species reduced model, curated) and **MODEL1907180003** (seven-species). — [Source](https://doi.org/10.1002/cnm.3220); [BioModels](https://biomodels.org/search?query=coagulation&numResults=100&format=json)
- **Papadopoulos KP, Gavaises M, Atkin C (2014).** "A simplified mathematical model for thrombin generation." *Med Eng Phys* 36(2):196–204. doi:10.1016/j.medengphy.2013.10.012 [V] — [Source](https://doi.org/10.1016/j.medengphy.2013.10.012)
- **Dunster2016 nondimensional coagulation model** in BioModels (BIOMD0000000925). The paper's citation was not verified. — [BioModels](https://biomodels.org/search?query=coagulation&numResults=100&format=json)

#### 3.7 Latest (2023–2026)
- **Ranc A, Bru S, Mendez S, Giansily-Blaizot M, Nicoud F, Méndez Rojano R (2023).** "Critical evaluation of kinetic schemes for coagulation." *PLoS ONE* 18(8):e0290531. doi:10.1371/journal.pone.0290531 [V]. Findings, from the page linked:
  - **Models and data.** It evaluates the Chatterjee 2010 (intrinsic) and Butenas 2004 (extrinsic) schemes against CAT in platelet-poor plasma. The deficiency series covered FXII, FXI and FVIII at 100, 50, 15, 5, 1 and <1%; the extrinsic runs used 1 pM tissue factor.
  - **Intrinsic model defect.** The original intrinsic model generates thrombin at 0% FXII, because FVIIa activates FIX and FX without TF. Those reactions were removed.
  - **Extrinsic model defect.** The extrinsic model is delayed and too weak. The authors set initial FVIIa–TF to 1 pM and raised the IXa→FX k_cat 60-fold, from 8.0×10⁻⁴ to 4.8×10⁻² s⁻¹.
  - **Remaining failure.** Both the original and modified intrinsic models fail on endogenous thrombin potential (ETP), peak thrombin and time-to-peak for FVIII < 1%.
  - **Parameter caution.** Table 1 gives initial concentrations. There is an internal discrepancy in the FXII activation rate k₁ (5×10³ in the text vs 5×10⁻⁴ to 5×10⁻³ in Table 2).
  - **Code and data.** Implemented in Python as mass-action ODEs. Data and code are on GitHub, `rodrigomrxvi/CriticalEvaluationCoagulation`. Fibrin reactions 50–57 were dropped to save time. — [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC10461854/)
- **Miyazawa K, Fogelson AL, Leiderman K (2023).** "Inhibition of platelet-surface-bound proteins during coagulation under flow I: TFPI." *Biophys J* 122(1):99–113. doi:10.1016/j.bpj.2022.11.023 [V] — [Source](https://doi.org/10.1016/j.bpj.2022.11.023)
- **Miyazawa K, Fogelson AL, Leiderman K (2023).** "… II: Antithrombin and heparin." *Biophys J* 122(1):230–240. doi:10.1016/j.bpj.2022.10.038 [V]
  - Adds antithrombin inhibition at activated platelet surfaces, with unfractionated and low-molecular-weight heparin.
  - The model is sensitive to AT levels, especially at low TF density. Inactivating no single enzyme substantially suppresses thrombin generation. — [Source](https://par.nsf.gov/servlets/purl/10430377); [record](https://ouci.dntb.gov.ua/works/4Kv2ryq7)
- **Montgomery D, Municchi F, Leiderman K (2023).** "clotFoam: an open-source framework to simulate blood clot formation under arterial flow." *SoftwareX* 23:101483. doi:10.1016/j.softx.2023.101483 [V]
  - Built on OpenFOAM. A continuum platelet model (advection, diffusion, aggregation in a dynamic fluid) and a simplified coagulation model with reactive boundary conditions for wall-bound species.
  - Code: github.com/d-montgomery/clotFoam. The code's licence is not stated on the abstract page; the arXiv paper is CC BY 4.0. — [Source](https://arxiv.org/abs/2304.09180)
- **Shankar KN, Diamond SL, Sinno T (2023).** "Development of a parallel multiscale 3D model for thrombus growth under flow." *Front Phys* 11:1256462. doi:10.3389/fphy.2023.1256462 [V] — [Source](https://doi.org/10.3389/fphy.2023.1256462)
- **Chung D, Bakshi S, van der Graaf PH (2023).** "A review of quantitative systems pharmacology models of the coagulation cascade: opportunities for improved usability." *Pharmaceutics* 15(3):918. doi:10.3390/pharmaceutics15030918 [V]
  - Notes that BioModels provides the code for many coagulation models, and that the field suffers from unclear connections between models and irreproducible code. — [Source](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10054658/)
- **Grande Gutiérrez N, Mukherjee D, Bark D (2024).** "Decoding thrombosis through code: a review of computational models." *J Thromb Haemost* 22(1):35–47. doi:10.1016/j.jtha.2023.08.021 [V]
  - Focuses on platelet-based thrombosis models, identifying thrombus-prone regions in patient-specific vasculature, and thrombus deformation and embolisation.
  - Concludes there is no one-size-fits-all model. — [Source](https://profiles.wustl.edu/en/publications/decoding-thrombosis-through-code-a-review-of-computational-models/)
- **Tuna R, Yi W, Crespo Cruz E, Romero J, Ren Y, Guan J, et al. (2024).** "Platelet biorheology and mechanobiology in thrombosis and hemostasis: perspectives from multiscale computation." *Int J Mol Sci* 25(9):4800. doi:10.3390/ijms25094800 [V]
  - Eleven authors. Also posted as a preprint, doi:10.20944/preprints202402.0799.v2. — [Source](https://doi.org/10.3390/ijms25094800)
- **Bouchnita A, Volpert V (2024).** "Intraplatelet calcium signaling regulates thrombus growth under flow: insights from a multiscale model." *Computation* 12(5):99. doi:10.3390/computation12050099 [V] — [Source](https://doi.org/10.3390/computation12050099)
- **Panteleev et al. (2024)**, the review in §3.4. [V]
- **Xu Z, Shi Y, Wu X, Wang S, He F, Hao P, et al. (2026).** "Multiscale insights into thrombus growth and detachment under non-physiological blood flow." *Comput Methods Programs Biomed* 281:109345. doi:10.1016/j.cmpb.2026.109345 [V]
  - Dissipative particle dynamics with coarse-grained cells: platelet activation, adhesion and fibrin formation. Validated against microfluidic experiments.
  - Flow intensity shifts the cellular-to-fibrin ratio.
  - Two detachment mechanisms: shear-driven boundary fragmentation, and pressure-gradient-driven internal layer separation. — [Source](https://pubmed.ncbi.nlm.nih.gov/41962341/)
- **Preprints and unverified venues, 2025–2026.** Both are [W].
  - "Platelet plug microstructure and flow modulate fibrin gelation dynamics: insights from computational simulations" (arXiv 2604.07844; also PMC13082137): a preformed discrete platelet plug, a reduced thrombin-generation model adapted from clotFoam, and fibrin polymerisation. Fibrin larger than monomers is treated as non-transported. — [arXiv](https://arxiv.org/html/2604.07844); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13082137/)
  - "Mesoscale simulations of thrombin activation and fibrin formation in microvascular and in vitro settings" (arXiv 2512.12859). Its methods were not confirmed. — [arXiv](https://arxiv.org/pdf/2512.12859)

#### 3.8 Fixtures available
- **ODE trajectories (best).** BIOMD0000000335 (Hockin 2002, curated by reproducing the paper's figures, with a COPASI run on BioSimulations and JWS Online's Fig. 3 experiment), BIOMD0000000336 (Jones & Mann 1994), BIOMD0000000362 (Butenas 2004), BIOMD0000000338–340 (Wajima 2009), BIOMD0000000755 (Hansen & Shadden nine-species), BIOMD0000000951 (Mitrophanov 2015). All are CC0.
- **Experimental CAT data with code.** Ranc 2023 GitHub repository, and Chelle 2018 (healthy and haemophilic).
- **Flow-coupled.** Bouchnita 2020 is in BioModels (non-curated). Link 2018 gives sensitivity indices for the Kuharsky–Fogelson / Leiderman–Fogelson lineage †. clotFoam provides an open 3D OpenFOAM implementation.

### Inferences
- **Stiffness and time integration.** Coagulation networks combine nanomolar-to-micromolar species with rate constants spanning many orders of magnitude, and a long lag before a rapid burst. Treat them as stiff:
  - BDF (CVODE) or Rosenbrock (ROS3P, MATLAB ode23s) for the well-mixed ODEs.
  - In flow coupling, Strang splitting: explicit or stabilised transport, then an implicit local reaction solve per cell or node. See §5.
  - These are my recommendations, not taken from a single source.
- **A fixture path.**
  1. Import the curated SBML (e.g., BIOMD0000000335).
  2. Reproduce the thrombin time course against BioSimulations/COPASI output to a stated tolerance.
  3. Compare readouts (lag, peak, ETP) with Butenas 1999 and Hemker-style CAT data.
  4. Repeat Ranc 2023's deficiency series as a falsification test.
  5. Only then couple to flow, using Kuharsky–Fogelson as the 0D-with-mass-transfer intermediate before the full spatial Leiderman–Fogelson model.
- **Licence compatibility.** BioModels SBML (CC0) is compatible with an MIT repository. clotFoam derives from OpenFOAM, which is GPL by my knowledge (not verified this session). Copying its code into an MIT repository would raise copyleft issues. Treat it as an external comparison only.
- **Venous models.** Venous models differ from arterial ones mainly in shear (low), residence time (long) and the role of fibrin over platelets. Govindarajan 2016 and Bouchnita 2020 are the citable starting points. Fibrin-enhanced low-shear models (Méndez Rojano 2022) bridge to stasis regions.

### Gaps
- **Venous thrombosis.** The physiology and valve-pocket sources are now in §3.4 [V]. Still missing:
  - Welsh et al. (2019, *J Clin Invest*), "Hemodynamic regulation of perivalvular endothelial gene expression prevents deep venous thrombosis." Crossref returned only Faculty Opinions records, so it is [U].
  - Any venous-valve model coupling a full coagulation network with a moving valve.
  - Published venous clot-growth reference trajectories.
  - The O₂-desaturation-in-the-valve-pocket → endothelial-activation chain (Bovill & van der Vliet 2011) has not, in what I found, been modelled as a coupled ADR problem. It is a natural research target that joins §1 and §3.
- No reference flow-coupled clotting benchmark with published numerical trajectories, comparable to the curated ODE fixtures, was found. Link 2018 and the clotFoam paper are the closest. Their supplementary data was not inspected.
- The parameter tables of Kuharsky & Fogelson 2001 and Leiderman & Fogelson 2011 exist in the papers † but were not re-read here.
- The licences of the Diamond-lab and Antaki-lab codes were not checked.
- 2024–2026 work by Fogelson and Leiderman after Miyazawa 2023 was not found in my searches.

---

## 4. Drug and Contrast-Agent Transport: Drug-Eluting Stents, Indicator Dilution, Tracer Kinetics (CT/MR Perfusion), and Porous-Medium Perfusion Models

### Takeaway
**Drug-eluting stents (DES).** The models couple:
- release from a polymer coating, by diffusion or dissolution;
- luminal advection–diffusion, in which flow patterns dictate deposition (Kolachalama 2009; Balakrishnan 2005);
- transport in a multilayer porous wall with reversible binding.

Saturable, nonlinear binding to specific receptors and to the extracellular matrix (ECM) is essential: McGinty & Pontrelli 2016; Tzafriri 2012; Bozsak 2014, who found paclitaxel convection-dominated and sirolimus binding-dominated.

**Contrast agents.** Indicator-dilution theory (Stewart 1897 / Hamilton 1932 → Meier–Zierler 1954 central-volume theorem) underlies both DSC deconvolution (Østergaard 1996, SVD; Wu 2003, oSVD) and DCE tracer kinetics (Tofts 1991/1999, extended Tofts, AATH). CT contrast at whole-body scale is a 0D compartment (PBPK) model (Bae 1998).

**Perfusion.** Multi-compartment Darcy continua (Michler 2013; Hyde 2013/2014) or 1D–3D coupling (D'Angelo & Quarteroni 2008), now combined with coronary trees and the whole heart (Papamanolis 2021; Di Gregorio 2021; Zingaro 2023; Montino Pelagi 2024) and with deforming hierarchical vasculature (Hohl 2026).

**Verification resources.** Strong, openly benchmarked material exists for tracer kinetics:
- the OSIPI DCE challenge (Shalom 2024) and lexicon (Dickie 2024);
- QIBA-style digital reference objects (Bosca & Jackson 2016);
- the Kudo 2013 digital phantom for CT and MR perfusion;
- the Parker 2006 population arterial input function (AIF).

### Cited Findings

#### 4.1 Drug-eluting stents and arterial uptake (classical)
- **Hwang CW, Wu D, Edelman ER (2001).** "Physiological transport forces govern drug distribution for stent-based delivery." *Circulation* 104(5):600–605. doi:10.1161/hc3101.092214 [V] — [Source](https://doi.org/10.1161/hc3101.092214)
- **Balakrishnan B, Tzafriri AR, Seifert P, Groothuis A, Rogers C, Edelman ER (2005).** "Strut position, blood flow, and drug deposition: implications for single and overlapping drug-eluting stents." *Circulation* 111(22):2958–2965. doi:10.1161/CIRCULATIONAHA.104.512475 [V] — [Source](https://doi.org/10.1161/CIRCULATIONAHA.104.512475)
- **Pontrelli G, de Monte F (2007).** "Mass diffusion through two-layer porous media: an application to the drug-eluting stent." *Int J Heat Mass Transf* 50(17–18):3658–3669. doi:10.1016/j.ijheatmasstransfer.2006.11.003 [V]
  - Analytical two-layer solutions †. A closed-form fixture. — [Source](https://doi.org/10.1016/j.ijheatmasstransfer.2006.11.003)
- **Kolachalama VB, Tzafriri AR, Arifin DY, Edelman ER (2009).** "Luminal flow patterns dictate arterial drug deposition in stent-based delivery." *J Control Release* 133(1):24–30. doi:10.1016/j.jconrel.2008.09.075 [V]
  - A corrigendum is published: *J Control Release* 146(1):160 (2010), doi:10.1016/j.jconrel.2010.05.007 [V]. — [Source](https://doi.org/10.1016/j.jconrel.2008.09.075)
- **Zunino P, D'Angelo C, Petrini L, Vergara C, Capelli C, Migliavacca F (2009).** "Numerical simulation of drug eluting coronary stents: mechanics, fluid dynamics and drug release." *Comput Methods Appl Mech Eng* 198(45–46):3633–3644. doi:10.1016/j.cma.2008.07.019 [V] — [Source](https://doi.org/10.1016/j.cma.2008.07.019)
- **Vairo G, Cioffi M, Cottone R, Dubini G, Migliavacca F (2010).** "Drug release from coronary eluting stents: a multidomain approach." *J Biomech* 43(8):1580–1589. doi:10.1016/j.jbiomech.2010.01.033 [V] — [Source](https://doi.org/10.1016/j.jbiomech.2010.01.033)
- **Tzafriri AR, Groothuis A, Price GS, Edelman ER (2012).** "Stent elution rate determines drug deposition and receptor-mediated effects." *J Control Release* 161(3):918–926. doi:10.1016/j.jconrel.2012.05.039 [V] — [Source](https://doi.org/10.1016/j.jconrel.2012.05.039)
- **Bozsak F, Chomaz JM, Barakat AI (2014).** "Modeling the transport of drugs eluted from stents: physical phenomena driving drug distribution in the arterial wall." *Biomech Model Mechanobiol* 13(2):327–347. doi:10.1007/s10237-013-0546-4 [V]
  - A multilayer wall with reversible binding. Paclitaxel transport is convection-dominated; sirolimus transport is binding-dominated. — [Source](https://www.researchgate.net/publication/259607779_Modeling_the_transport_of_drugs_eluted_from_stents_Physical_phenomena_driving_drug_distribution_in_the_arterial_wall)
- **Bozsak F, Gonzalez-Rodriguez D, Sternberger Z, Belitz P, Bewley T, Chomaz JM, Barakat AI (2015).** "Optimization of drug delivery by drug-eluting stents." *PLoS ONE* 10(6):e0130182. doi:10.1371/journal.pone.0130182 [V]
  - Paclitaxel stents perform optimally with very fast (hours) or very slow (months) release, at lower concentrations. Sirolimus stents perform optimally only with slow release. — [Source](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4470631/)
- **McGinty S (2014).** "A decade of modelling drug release from arterial stents." *Math Biosci* 257:80–90. doi:10.1016/j.mbs.2014.06.016 [V] — review. — [Source](https://doi.org/10.1016/j.mbs.2014.06.016)
- **McGinty S, Pontrelli G (2016).** "On the role of specific drug binding in modelling arterial eluting stents." *J Math Chem* 54(4):967–976. doi:10.1007/s10910-016-0618-7 [V]
  - Simulating ECM-bound and specific-receptor-bound phases separately is important if receptor saturation is to be shown. — [Source](https://link.springer.com/article/10.1007/s10910-016-0618-7)
- **Lucena RM, Mangiavacchi N, Pontes J, Anjos G, McGinty S (2018).** "On the transport through polymer layer and porous arterial wall in drug-eluting stents." *J Braz Soc Mech Sci Eng* 40(12):572. doi:10.1007/s40430-018-1487-8 [V]
  - Only 52–54% of the sirolimus mass diffuses into the wall. — [Source](https://link.springer.com/article/10.1007/s40430-018-1487-8)
- **McKittrick CM, McKee S, Kennedy S, Oldroyd K, Wheel M, Pontrelli G, et al. (2019).** "Combining mathematical modelling with in vitro experiments to predict in vivo drug-eluting stent performance." *J Control Release* 303:151–161. doi:10.1016/j.jconrel.2019.03.012 [V] — [Source](https://doi.org/10.1016/j.jconrel.2019.03.012)
- **Escuer J, Aznar I, McCormick C, Peña E, McGinty S, Martínez MA (2021).** "Influence of vessel curvature and plaque composition on drug transport in the arterial wall following drug-eluting stent implantation." *Biomech Model Mechanobiol* 20(2):767–786. [V for authors, venue, volume and pages; the DOI was truncated in my check (10.1007/s10237-020-0…) and must be looked up.] — [Crossref search](https://api.crossref.org/works?query.bibliographic=Escuer+vessel+curvature+plaque+composition+drug+transport)
- **Escuer J, Fensterseifer Schmidt A, Peña E, Martínez MA, McGinty S (2022) (bridge).** "Mathematical modelling of endovascular drug delivery: balloons versus stents." *Int J Pharm* 620:121742. doi:10.1016/j.ijpharm.2022.121742 [V] — [Source](https://doi.org/10.1016/j.ijpharm.2022.121742)

#### 4.2 Drug-eluting stents: latest (2023–2026)
- **Salvi S, Jain A, Pontrelli G, McGinty S (2023).** "Modeling dual drug delivery from eluting stents: the influence of non-linear binding competition and non-uniform drug loading." *Pharm Res* 40(1):215–230. doi:10.1007/s11095-022-03419-3 [V]
  - Competition for a limited density of non-specific binding sites, and large differences in forward and backward rate coefficients, determine uptake. — [Source](https://link.springer.com/article/10.1007/s11095-022-03419-3)
- No peer-reviewed 2024–2026 DES transport model surfaced in two web searches and the Crossref discovery queries. A 2023 Preprints.org 2D finite-volume PLGA-strut model with reversible binding exists but is not peer-reviewed. — [Preprints.org](https://www.preprints.org/manuscript/202305.0090)

#### 4.3 Indicator dilution and tracer kinetics (classical)
- **Stewart GN (1897).** "Researches on the circulation time and on the influences which affect it. IV. The output of the heart." *J Physiol* 22(3):159–183. doi:10.1113/jphysiol.1897.sp000684 [V]
  - Crossref gives the title without "IV. The output of the heart" †. — [Source](https://doi.org/10.1113/jphysiol.1897.sp000684)
- **Hamilton WF, Moore JW, Kinsman JM, Spurling RG (1932).** "Studies on the circulation. IV. Further analysis of the injection method, and of changes in hemodynamics under physiological and pathological conditions." *Am J Physiol* 99(3):534–551. doi:10.1152/ajplegacy.1932.99.3.534 [V]
  - Crossref has only the main title; the subtitle is †.
  - Stewart–Hamilton: F = m/∫₀^∞ c(t) dt, with the recirculation tail extrapolated †. — [Source](https://doi.org/10.1152/ajplegacy.1932.99.3.534)
- **Meier P, Zierler KL (1954).** "On the theory of the indicator-dilution method for measurement of blood flow and volume." *J Appl Physiol* 6(12):731–744. doi:10.1152/jappl.1954.6.12.731 [V]
  - The central-volume theorem V = F·t̄, where t̄ is the mean of the transit-time density h(t) †. — [Source](https://doi.org/10.1152/jappl.1954.6.12.731)
- **Zierler KL (1962).** "Theoretical basis of indicator-dilution methods for measuring flow and volume." *Circ Res* 10(3):393–407. doi:10.1161/01.RES.10.3.393 [V] — [Source](https://doi.org/10.1161/01.RES.10.3.393)
- **Zierler K (2000).** "Indicator dilution methods for measuring blood flow, volume, and other properties of biological systems: a brief history and memoir." *Ann Biomed Eng* 28(8):836–848. doi:10.1114/1.1308496 [V] — [Source](https://doi.org/10.1114/1.1308496)
- **Tofts PS, Kermode AG (1991).** "Measurement of the blood-brain barrier permeability and leakage space using dynamic MR imaging. 1. Fundamental concepts." *Magn Reson Med* 17(2):357–367. doi:10.1002/mrm.1910170208 [V] — [Source](https://doi.org/10.1002/mrm.1910170208)
- **Tofts PS, Brix G, Buckley DL, Evelhoch JL, Henderson E, Knopp MV, et al. (1999).** "Estimating kinetic parameters from dynamic contrast-enhanced T₁-weighted MRI of a diffusable tracer: standardized quantities and symbols." *J Magn Reson Imaging* 10(3):223–232. doi:10.1002/(SICI)1522-2586(199909)10:3<223::AID-JMRI2>3.0.CO;2-S [V]
  - Thirteen authors.
  - The (extended) Tofts model †: C_t(t) = v_p·C_p(t) + K^trans ∫₀^t C_p(τ) e^(−k_ep(t−τ)) dτ, with k_ep = K^trans/v_e. — [Source](https://doi.org/10.1002/(SICI)1522-2586(199909)10:3%3C223::AID-JMRI2%3E3.0.CO;2-S)
- **Østergaard L, Weisskoff RM, Chesler DA, Gyldensted C, Rosen BR (1996).** "High resolution measurement of cerebral blood flow using intravascular tracer bolus passages. Part I: Mathematical approach and statistical analysis." *Magn Reson Med* 36(5):715–725. doi:10.1002/mrm.1910360510 [V]
  - Tissue curve = CBF·(C_a ⊗ R)(t). R is recovered by truncated-SVD deconvolution. CBV is the ratio of areas, and MTT = CBV/CBF (central volume) †. — [Source](https://doi.org/10.1002/mrm.1910360510)
- **Østergaard L, Sorensen AG, Kwong KK, Weisskoff RM, Gyldensted C, Rosen BR (1996).** "… Part II: Experimental comparison and preliminary results." *Magn Reson Med* 36(5):726–736. doi:10.1002/mrm.1910360511 [V] — [Source](https://doi.org/10.1002/mrm.1910360511)
- **St Lawrence KS, Lee TY (1998).** "An adiabatic approximation to the tissue homogeneity model for water exchange in the brain: I. Theoretical derivation." *J Cereb Blood Flow Metab* 18(12):1365–1377. doi:10.1097/00004647-199812000-00011 [V] — the AATH model. — [Source](https://doi.org/10.1097/00004647-199812000-00011)
- **Wu O, Østergaard L, Weisskoff RM, Benner T, Rosen BR, Sorensen AG (2003).** "Tracer arrival timing-insensitive technique for estimating flow in MR perfusion-weighted imaging using singular value decomposition with a block-circulant deconvolution matrix." *Magn Reson Med* 50(1):164–174. doi:10.1002/mrm.10522 [V] — oSVD. — [Source](https://doi.org/10.1002/mrm.10522)
- **Parker GJM, Roberts C, Macdonald A, Buonaccorsi GA, Cheung S, Buckley DL, et al. (2006).** "Experimentally-derived functional form for a population-averaged high-temporal-resolution arterial input function for dynamic contrast-enhanced MRI." *Magn Reson Med* 56(5):993–1000. doi:10.1002/mrm.21066 [V]
  - Ten authors. A closed-form AIF: two Gaussians plus a sigmoid-modulated exponential †. A standard synthetic input. — [Source](https://doi.org/10.1002/mrm.21066)
- **Sourbron SP, Buckley DL (2011).** "On the scope and interpretation of the Tofts models for DCE-MRI." *Magn Reson Med* 66(3):735–745. doi:10.1002/mrm.22861 [V] — [Source](https://doi.org/10.1002/mrm.22861)
- **Sourbron SP, Buckley DL (2013).** "Classic models for dynamic contrast-enhanced MRI." *NMR Biomed* 26(8):1004–1027. doi:10.1002/nbm.2940 [V] — the standard review of compartment models. — [Source](https://doi.org/10.1002/nbm.2940)
- **Bae KT, Heiken JP, Brink JA (1998).** "Aortic and hepatic contrast medium enhancement at CT. Part I. Prediction with a computer model." *Radiology* 207(3):647–655. doi:10.1148/radiology.207.3.9609886 [V]
  - A whole-body compartment (PBPK) model of iodinated contrast †. A natural 0D network test case. — [Source](https://doi.org/10.1148/radiology.207.3.9609886)
- **Bae KT (2010).** "Intravenous contrast medium administration and scan timing at CT: considerations and approaches." *Radiology* 256(1):32–61. doi:10.1148/radiol.10090908 [V] — [Source](https://doi.org/10.1148/radiol.10090908)
- **Eslami P, Seo JH, Rahsepar AA, George R, Lardo AC, Mittal R (2015).** "Computational study of computed tomography contrast gradients in models of stenosed coronary arteries." *J Biomech Eng* 137(9):091002. doi:10.1115/1.4030891 [V]
  - 3D CFD with contrast advection–diffusion. — [Source](https://doi.org/10.1115/1.4030891)
- **Kudo K, Christensen S, Sasaki M, Østergaard L, Shirato H, Ogasawara K, et al. (2013).** "Accuracy and reliability assessment of CT and MR perfusion analysis software using a digital phantom." *Radiology* 267(1):201–211. doi:10.1148/radiol.12112618 [V] — [Source](https://doi.org/10.1148/radiol.12112618)
- **Bosca RJ, Jackson EF (2016).** "Creating an anthropomorphic digital MR phantom—an extensible tool for comparing and evaluating quantitative imaging algorithms." *Phys Med Biol* 61(2):974–982. doi:10.1088/0031-9155/61/2/974 [V] — [Source](https://doi.org/10.1088/0031-9155/61/2/974)

#### 4.4 Tracer kinetics: latest (2023–2026)
- **Dickie BR, Ahmed Z, Arvidsson J, Bell LC, Buckley DL, Debus C, et al. (2024).** "A community-endorsed open-source lexicon for contrast agent-based perfusion MRI: a consensus guidelines report from the ISMRM Open Science Initiative for Perfusion Imaging (OSIPI)." *Magn Reson Med* 91(5):1761–1773. doi:10.1002/mrm.29840 [V] — 16 authors. — [Source](https://doi.org/10.1002/mrm.29840)
- **Shalom ER, et al. (2024).** "The ISMRM Open Science Initiative for Perfusion Imaging (OSIPI): results from the OSIPI–Dynamic Contrast-Enhanced challenge." *Magn Reson Med* 91(5):1803–1821. doi:10.1002/mrm.29909 [W; the full author list was not taken]
  - Across 10 submissions the OSIPI-gold score ranged from 28% to 78% (median 59%).
  - Accuracy ranged 0.54–0.92, repeatability 0.64–0.86 and reproducibility 0.65–1.00.
  - There was high inter-software variability in K^trans. Manual AIF selection worsened reproducibility, and a detailed standard operating procedure was critical. — [White Rose eprint](https://eprints.whiterose.ac.uk/id/eprint/207709/)
- **The OSIPI DCE-DSC-MRI code collection.**
  - Python implementations of T₁ mapping, signal-to-concentration conversion and pharmacokinetic models, tested with Pytest and GitHub Actions.
  - It is now closed to contributions, and development has moved to the `osipy` library. An ISMRM 2022 abstract reported 74 contributions collected and 14 tested, with DSC behind DCE.
  - The repository's licence was not checked. — [GitHub](https://github.com/OSIPI/DCE-DSC-MRI_CodeCollection); [ISMRM 2022 abstract, van Houdt et al., doi:10.58530/2022/2691 [V]](https://archive.ismrm.org/2022/2691.html)
- **OSIPI deep-learning DCE challenge (2024).** Kargar S, Kershaw L, Kazerooni AF, Bell L, Van der Heijden R, Mutsaerts HJ, et al. "The ISMRM OSIPI: a challenge for reproducible DCE-MRI AI-based analysis." ISMRM 2024 abstract 0245, doi:10.58530/2024/0245 [V]. Uses uterine DCE data. — [Source](https://cds.ismrm.org/protected/24MProceedings/PDFfiles/0245_En1mDaR8b.html)
- **Chung KJ, De Sarno D, Lee TY (2023).** "CT perfusion stroke lesion threshold calibration between deconvolution algorithms." *Sci Rep* 13:21458. doi:10.1038/s41598-023-48700-6 [V] — [Source](https://doi.org/10.1038/s41598-023-48700-6)

#### 4.5 Porous-medium and multiscale perfusion (classical)
- **D'Angelo C, Quarteroni A (2008).** "On the coupling of 1D and 3D diffusion-reaction equations: application to tissue perfusion problems." *Math Models Methods Appl Sci* 18(8):1481–1504. doi:10.1142/S0218202508003108 [V] — [Source](https://doi.org/10.1142/S0218202508003108)
- **Chapelle D, Gerbeau JF, Sainte-Marie J, Vignon-Clementel IE (2010).** "A poroelastic model valid in large strains with applications to perfusion in cardiac modeling." *Comput Mech* 46(1):91–101. doi:10.1007/s00466-009-0452-x [V] — [Source](https://doi.org/10.1007/s00466-009-0452-x)
- **Cookson AN, Lee J, Michler C, Chabiniok R, Hyde E, Nordsletten DA, et al. (2012).** "A novel porous mechanical framework for modelling the interaction between coronary perfusion and myocardial mechanics." *J Biomech* 45(5):850–855. doi:10.1016/j.jbiomech.2011.11.026 [V] — [Source](https://doi.org/10.1016/j.jbiomech.2011.11.026)
- **Michler C, Cookson AN, Chabiniok R, Hyde E, Lee J, Sinclair M, et al. (2013).** "A computationally efficient framework for the simulation of cardiac perfusion using a multi-compartment Darcy porous-media flow model." *Int J Numer Method Biomed Eng* 29(2):217–232. doi:10.1002/cnm.2520 [V]
  - Multi-compartment Darcy flow with inter-compartment exchange, homogenised from cryomicrotome vessel data in porcine myocardium. — [Source](https://doi.org/10.1002/cnm.2520); [context](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10453666/)
- **Hyde ER, Michler C, Lee J, Cookson AN, Chabiniok R, Nordsletten DA, Smith NP (2013).** "Parameterisation of multi-scale continuum perfusion models from discrete vascular networks." *Med Biol Eng Comput* 51(5):557–570. doi:10.1007/s11517-012-1025-2 [V] — [Source](https://doi.org/10.1007/s11517-012-1025-2)
- **Hyde ER, Cookson AN, Lee J, Michler C, Goyal A, Sochi T, et al. (2014).** "Multi-scale parameterisation of a myocardial perfusion model using whole-organ arterial networks." *Ann Biomed Eng* 42(4):797–811. doi:10.1007/s10439-013-0951-y [V]
  - Thirteen authors. This is myocardial, not liver: the brief's "liver/heart" pairing conflates it with the liver work below. — [Source](https://doi.org/10.1007/s10439-013-0951-y)
- **Liver.**
  - Debbaut C, Vierendeels J, Casteleyn C, et al. (2012). "Perfusion characteristics of the human hepatic microcirculation based on three-dimensional reconstructions and computational fluid dynamic analysis." *J Biomech Eng* 134(1):011003. doi:10.1115/1.4005545 [V] — [Source](https://doi.org/10.1115/1.4005545)
  - Debbaut C, Vierendeels J, Siggers JH, Repetto R, Monbaliu D, Segers P (2014). "A 3D porous media liver lobule model: the importance of vascular septa and anisotropic permeability for homogeneous perfusion." *Comput Methods Biomech Biomed Eng* 17(12):1295–1310. doi:10.1080/10255842.2012.744399 [V] — [Source](https://doi.org/10.1080/10255842.2012.744399)
  - Ricken T, Werner D, Holzhütter HG, König M, Dahmen U, Dirsch O (2015). "Modeling function–perfusion behavior in liver lobules including tissue, blood, glucose, lactate and glycogen by use of a coupled two-scale PDE–ODE approach." *Biomech Model Mechanobiol* 14(3):515–536. doi:10.1007/s10237-014-0619-z [V] — [Source](https://doi.org/10.1007/s10237-014-0619-z)
- **Cattaneo L, Zunino P (2014).** "A computational model of drug delivery through microcirculation to compare different tumor treatments." *Int J Numer Method Biomed Eng* 30(11):1347–1371. doi:10.1002/cnm.2661 [V] — [Source](https://doi.org/10.1002/cnm.2661)
- **Papamanolis L, Kim HJ, Jaquet C, Sinclair M, Schaap M, Danad I, et al. (2021).** "Myocardial perfusion simulation for coronary artery disease: a coupled patient-specific multiscale model." *Ann Biomed Eng* 49(5):1432–1447. doi:10.1007/s10439-020-02681-z [V]
  - 1D coronary flow coupled to a single-compartment porous myocardium. It cites the cost and parameterisation difficulty of multi-compartment models. — [Source](https://link.springer.com/article/10.1007/s10439-020-02681-z)
- **Di Gregorio S, Fedele M, Pontone G, Corno AF, Zunino P, Vergara C, Quarteroni A (2021).** "A computational model applied to myocardial perfusion in the human heart: from large coronaries to microvasculature." *J Comput Phys* 424:109836. doi:10.1016/j.jcp.2020.109836 [V] — [Source](https://doi.org/10.1016/j.jcp.2020.109836)

#### 4.6 Perfusion: latest (2023–2026)
- **Zingaro A, Vergara C, Dede' L, Regazzoni F, Quarteroni A (2023).** "A comprehensive mathematical model for cardiac perfusion." *Sci Rep* 13:14220. doi:10.1038/s41598-023-41312-0 [V]
  - Multi-compartment Darcy perfusion inside an electrophysiology–mechanics–haemodynamics heart model. — [Source](https://www.nature.com/articles/s41598-023-41312-0)
- **Alves JR, Berg LA, Gaio ED, Rocha BM, de Queiroz RAB, dos Santos RW (2023).** "A hybrid model for cardiac perfusion: coupling a discrete coronary arterial tree model with a continuous porous-media flow model of the myocardium." *Entropy* 25(8):1229. doi:10.3390/e25081229 [V]
  - Includes contrast-agent transport with fibrosis–contrast interaction and recirculation, to simulate infarction and gadolinium late-enhancement MRI. A direct template for "contrast in porous perfusion". — [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC10453666/)
- **Montino Pelagi G, Regazzoni F, Huyghe JM, Baggiano A, Alì M, Bertoluzza S, et al. (2024).** "Modeling cardiac microcirculation for the simulation of coronary flow and 3D myocardial perfusion." *Biomech Model Mechanobiol* 23(6):1863–1888. doi:10.1007/s10237-024-01873-z [V] — [Source](https://doi.org/10.1007/s10237-024-01873-z)
- **Wang Y, Yin X (2024).** "Modelling coronary flow and myocardial perfusion by integrating a structured-tree coronary flow model and a hyperelastic left ventricle model." *Comput Methods Programs Biomed* 243:107928. doi:10.1016/j.cmpb.2023.107928 [V] — [Source](https://doi.org/10.1016/j.cmpb.2023.107928)
- **Hohl J, Ebrahem A, Jessen E, ten Eikelder MFP, Schillinger D (2026).** "A multi-compartment homogenized perfusion model for deforming hierarchical vasculature." *Biomech Model Mechanobiol* 25(1):10. doi:10.1007/s10237-025-02026-6 [V] — [Source](https://doi.org/10.1007/s10237-025-02026-6)
- **p-Brain (2026).** "A modular open-source framework for automated quantitative DCE-MRI of cerebral perfusion, microvasculature, and blood–brain barrier permeability." *Magn Reson Med* 96(6). doi:10.1002/mrm.70606 [W]
  - Crossref returned no authors for the journal record. The preprint is by Tireli ED, Cramer SP, Lindberg U, Tireli D, Vestergaard MB, Larsson HBW (doi:10.64898/2026.02.15.705995). — [Source](https://doi.org/10.1002/mrm.70606)

#### 4.7 Fixtures available
- **Closed forms.**
  - Two-layer DES diffusion (Pontrelli & de Monte 2007 †).
  - Stewart–Hamilton and central-volume identities on synthetic gamma-variate curves.
  - The Tofts / extended Tofts convolution with the Parker AIF (closed-form AIF).
- **Benchmarks.**
  - The OSIPI DCE challenge data and scores (Shalom 2024), and the OSIPI code collection's Pytest reference values.
  - The digital MR phantom (Bosca & Jackson 2016).
  - The CT and MR perfusion digital phantom (Kudo 2013).
- **Whole-body 0D contrast.** Bae 1998 aortic and hepatic enhancement curves †.

### Inferences
- **Contrast transport in a 0D/1D/3D stack.** The natural verification ladder:
  1. 0D compartments (Bae 1998). Fluxes are flows between compartments; a closed linear ODE system.
  2. 1D advection along a network.
  3. 3D lumen ADR (Eslami 2015).
  4. Porous tissue uptake (Alves 2023).

  At each level, check that the simulated curves satisfy the Meier–Zierler central-volume theorem, V = F·MTT, on the simulated residue function. That is an exact, model-independent invariant.
- **Deconvolution is ill-posed.** It needs regularisation: truncated SVD, oSVD, or Bayesian methods. Use the Kudo phantom and the OSIPI data to set tolerances, not exact agreement.
- **DES models.** The minimal credible model is a dissolution- or diffusion-release coating, luminal ADR, and a two-phase (ECM plus receptor) saturable-binding wall. Binding is stiff (fast on, slow off), so an implicit local reaction solve is needed.

### Gaps
- No 2024–2026 peer-reviewed DES transport model was found; the field may have moved to drug-coated balloons.
- Kety SS (1951), "The theory and applications of the exchange of inert gas at the lungs and tissues," *Pharmacol Rev* 3(1):1–41, is [U]: it has no DOI and was not checked.
- 1D network models of contrast or drug transport in arterial trees (e.g., on ADAN or Müller–Toro networks) were not located. Discovery searches returned only hydrology papers.
- The QIBA DCE digital reference object's citation (Barboriak and colleagues), the OSIPI code collection's licence, and the full author list of Shalom 2024 were not verified.

---

## 5. Numerical Methods for Advection-Dominated Transport in Blood Flow, and Stiff Reaction Networks

### Takeaway
**Transport discretisation.** At Péclet numbers of 10⁵–10⁹, Galerkin finite elements oscillate. The standard remedies are SUPG (Brooks & Hughes 1982), GLS (Hughes, Franca & Hulbert 1989) and their ADR comparisons (Codina 1998), together with boundary-layer meshing or reduced near-wall models (Hansen & Shadden 2016; Hansen, Arzani & Shadden 2019). Finite-volume alternatives exist (clotFoam on OpenFOAM), as do lattice-Boltzmann advection–diffusion schemes (Ginzburg 2005; Chopard 2009; Krüger 2017).

**Stiff chemistry.** Coagulation chemistry is stiff. Use BDF (CVODE in SUNDIALS) or Rosenbrock methods (ROS3P; Shampine & Reichelt's ode23s), coupled to transport by Strang or IMEX splitting.

**Residence time.** It is a useful scalar for low-shear thrombosis and can be computed by an Eulerian ADR equation (Esmaily-Moghadam 2013).

### Cited Findings

#### 5.1 Stabilised finite elements for advection–diffusion (classical)
- **Brooks AN, Hughes TJR (1982).** "Streamline upwind/Petrov-Galerkin formulations for convection dominated flows with particular emphasis on the incompressible Navier-Stokes equations." *Comput Methods Appl Mech Eng* 32(1–3):199–259. doi:10.1016/0045-7825(82)90071-8 [V] — [Source](https://doi.org/10.1016/0045-7825(82)90071-8)
- **Hughes TJR, Franca LP, Hulbert GM (1989).** "A new finite element formulation for computational fluid dynamics: VIII. The Galerkin/least-squares method for advective-diffusive equations." *Comput Methods Appl Mech Eng* 73(2):173–189. doi:10.1016/0045-7825(89)90111-4 [V] — [Source](https://doi.org/10.1016/0045-7825(89)90111-4)
- **Codina R (1998).** "Comparison of some finite element methods for solving the diffusion-convection-reaction equation." *Comput Methods Appl Mech Eng* 156(1–4):185–210. doi:10.1016/S0045-7825(97)00206-5 [V]
  - Covers the reaction-dominated regime, which matters for NO and binding †. — [Source](https://doi.org/10.1016/S0045-7825(97)00206-5)
- **Hansen KB, Arzani A, Shadden SC (2019).** "Finite element modeling of near-wall mass transport in cardiovascular flows." *Int J Numer Method Biomed Eng* 35(1):e3148. doi:10.1002/cnm.3148 [V] — [Source](https://doi.org/10.1002/cnm.3148)
- **Esmaily-Moghadam M, Hsia TY, Marsden AL (2013).** "A non-discrete method for computation of residence time in fluid mechanics simulations." *Phys Fluids* 25(11):110802. doi:10.1063/1.4819142 [V] — [Source](https://doi.org/10.1063/1.4819142)
- **Arthurs CJ, Khlebnikov R, Melville A, Marčan M, Gomez A, Dillon-Murphy D, et al. (2021).** "CRIMSON: an open-source software framework for cardiovascular integrated modelling and simulation." *PLoS Comput Biol* 17(5):e1008881. doi:10.1371/journal.pcbi.1008881 [V]
  - 34 authors. Whether it ships a scalar-transport solver was not checked. — [Source](https://doi.org/10.1371/journal.pcbi.1008881)

#### 5.2 Lattice-Boltzmann advection–diffusion (classical)
- **Ginzburg I (2005).** "Equilibrium-type and link-type lattice Boltzmann models for generic advection and anisotropic-dispersion equation." *Adv Water Resour* 28(11):1171–1195. doi:10.1016/j.advwatres.2005.03.004 [V] — [Source](https://doi.org/10.1016/j.advwatres.2005.03.004)
- **Chopard B, Falcone JL, Latt J (2009).** "The lattice Boltzmann advection-diffusion model revisited." *Eur Phys J Spec Top* 171(1):245–249. doi:10.1140/epjst/e2009-01035-5 [V] — [Source](https://doi.org/10.1140/epjst/e2009-01035-5)
- **Krüger T, Kusumaatmaja H, Kuzmin A, Shardt O, Silva G, Viggen EM (2017).** *The Lattice Boltzmann Method: Principles and Practice.* Springer (Graduate Texts in Physics). doi:10.1007/978-3-319-44649-3 [V] — [Source](https://doi.org/10.1007/978-3-319-44649-3)

#### 5.3 Stiff integrators and splitting (classical)
- **Strang G (1968).** "On the construction and comparison of difference schemes." *SIAM J Numer Anal* 5(3):506–517. doi:10.1137/0705041 [V] — [Source](https://doi.org/10.1137/0705041)
- **Hairer E, Wanner G (1996).** *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*, 2nd ed. Springer. doi:10.1007/978-3-642-05221-7 [V] — [Source](https://doi.org/10.1007/978-3-642-05221-7)
- **Shampine LF, Reichelt MW (1997).** "The MATLAB ODE Suite." *SIAM J Sci Comput* 18(1):1–22. doi:10.1137/S1064827594276424 [V]
  - ode15s (NDF/BDF) and ode23s (modified Rosenbrock) †. — [Source](https://doi.org/10.1137/S1064827594276424)
- **Ascher UM, Ruuth SJ, Spiteri RJ (1997).** "Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations." *Appl Numer Math* 25(2–3):151–167. doi:10.1016/S0168-9274(97)00056-1 [V] — [Source](https://doi.org/10.1016/S0168-9274(97)00056-1)
- **Lang J, Verwer J (2001).** "ROS3P—an accurate third-order Rosenbrock solver designed for parabolic problems." *BIT Numer Math* 41(4):731–738. doi:10.1023/A:1021900219772 [V] — [Source](https://doi.org/10.1023/A:1021900219772)
- **Hindmarsh AC, Brown PN, Grant KE, Lee SL, Serban R, Shumaker DE, Woodward CS (2005).** "SUNDIALS: suite of nonlinear and differential/algebraic equation solvers." *ACM Trans Math Softw* 31(3):363–396. doi:10.1145/1089014.1089020 [V]
  - SUNDIALS is BSD-licensed by my knowledge; not verified this session. — [Source](https://doi.org/10.1145/1089014.1089020)

#### 5.4 Latest (2023–2026)
- **Montgomery, Municchi & Leiderman (2023), clotFoam** (see §3.7). An open finite-volume implementation on OpenFOAM coupling flow, a continuum platelet model and reactive wall boundary conditions. — [Source](https://arxiv.org/abs/2304.09180)
- **Eggington L, Entwistle J, Xu X, Halliday I, Spendlove J (2026), preprint.** "Benchmarking lattice Boltzmann advection–diffusion models for residence time prediction in blood flow." SSRN, doi:10.2139/ssrn.7584003 [V as a preprint record; not peer-reviewed] — [Source](https://doi.org/10.2139/ssrn.7584003)
- **Key K, Abdelmalik MRA, Elgeti S, Hughes TJR, Baidoo F (2023), preprint.** "Finite element and isogeometric stabilized methods for the advection-diffusion-reaction equation." SSRN, doi:10.2139/ssrn.4473578 [V as a preprint record] — [Source](https://doi.org/10.2139/ssrn.4473578)
- **Bridge: Shankar et al. (2022, 2023).** A reduced coagulation cascade inside a 3D multiscale flow solver with single-platelet resolution, and its parallel version. These show that reduced networks are the practical choice inside 3D flow. — [Source](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8827456/); [2023](https://doi.org/10.3389/fphy.2023.1256462)

### Inferences
- **A scheme fitted to this infrastructure's calx foundation.** calx already holds the θ-method for large linear systems, RK and Störmer–Verlet, and conjugate gradients.
  1. Transport by SUPG/GLS finite elements or upwind finite volumes, at the θ-method's implicit level.
  2. Reactions as a separate stiff ODE per node, by BDF or Rosenbrock (ROS3P is 3rd order and L-stable for parabolic problems).
  3. Strang splitting between them, with a convergence test that halving Δt reduces the splitting error fourfold.

  This keeps "the law apart from the field and from the integrator", as physis's `continuum` already does for heat conduction. Rosenbrock methods need only the Jacobian and linear solves, with no Newton iteration, and are easy to write on calx's dense and sparse systems.
- **Verification tests for transport.**
  - 1D steady advection–diffusion with exact exponential boundary layers (Pe up to 10⁶).
  - The Graetz/Lévêque wall-flux problem in a tube. The Lévêque solution is standard, but Lévêque 1928 is [U].
  - Rotating-cone or Gaussian-hill pure advection for numerical diffusion.
  - Conservation to rounding error in finite-volume or locally conservative finite-element schemes.
- **Verification tests for reactions.**
  - Mass-action conservation laws: total prothrombin plus thrombin plus complexes is constant, absent production.
  - Reproduction of the BioModels SBML trajectories.
  - Order-of-accuracy tests on a linear test equation y′ = λy with λ ≪ 0, for stability.

### Gaps
- No peer-reviewed 2023–2026 benchmark comparing SUPG, finite volume and LBM specifically for haemodynamic mass transport was found. The 2026 SSRN preprint on LBM residence time is the nearest.
- Lévêque MA (1928), "Les lois de la transmission de chaleur par convection," *Ann Mines* 13:201–299, is [U]: it has no DOI and was not verified.
- The licences of OpenFOAM and SUNDIALS are stated from my knowledge only.
