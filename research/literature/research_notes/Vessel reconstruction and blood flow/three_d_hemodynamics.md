# Three-Dimensional Simulation of Blood Flow in Image-Derived Vessels

Scope: finite-element/finite-volume CFD and the lattice Boltzmann method (LBM) on geometries reconstructed from CT/MRI; blood rheology; fluid–structure interaction (FSI); haemodynamic indices; verification and validation; clinical uses (arterial and venous); latest work (2023–2026).

Verification note (2026-11-04): every DOI below was resolved through doi.org content negotiation (Crossref/DataCite metadata), and authors, year, volume, issue and pages were checked against that record unless marked **[unverified]**. Where the online-first year differs from the print year, the print year and volume are given. Licences were read from the repositories' licence files or the GitLab/GitHub licence APIs on 2026-11-04. Each section separates **Classical** work from **Latest (2023–2026)** work.

---

## 1. Image-Based CFD Practice: Geometry, Mesh, Inflow, Outflow, Backflow and Transition

### Takeaway
The standard pipeline is: segment the lumen, build a surface, mesh it, prescribe an inlet flow waveform, couple the outlets to reduced (0D Windkessel or 1D) models, and solve the incompressible Navier–Stokes equations with backflow stabilisation. The largest real-world errors come from segmentation, inflow rates and outflow strategy rather than from the solver. In stenoses and some aneurysms, resolving transitional flow needs low-dissipation numerics and fine time steps (roughly 800k tetrahedra and about 5000 steps per cycle for aneurysms).

### Cited Findings

**Classical reviews and frameworks**
- Steinman DA (2002). "Image-based computational fluid dynamics modeling in realistic arterial geometries." *Annals of Biomedical Engineering* 30(4):483–497. doi:10.1114/1.1467679. This is the founding review of image-based CFD of arterial flow and its links to atherosclerosis. — [doi](https://doi.org/10.1114/1.1467679)
- Taylor CA, Hughes TJR, Zarins CK (1998). "Finite element modeling of blood flow in arteries." *Computer Methods in Applied Mechanics and Engineering* 158(1–2):155–196. doi:10.1016/S0045-7825(98)80008-X. It describes an integrated system (geometry, meshing, stabilised FE solver) for vascular flow. — [doi](https://doi.org/10.1016/S0045-7825(98)80008-X)
- Taylor CA, Figueroa CA (2009). "Patient-specific modeling of cardiovascular mechanics." *Annual Review of Biomedical Engineering* 11:109–134. doi:10.1146/annurev.bioeng.10.061807.160521. It covers model construction, material properties, boundary conditions and solution methods, with free full text at PMC4581431. — [doi](https://doi.org/10.1146/annurev.bioeng.10.061807.160521); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4581431)
- Taylor CA, Steinman DA (2010). "Image-based modeling of blood flow and vessel wall dynamics: applications, methods and future directions." *Annals of Biomedical Engineering* 38(3):1188–1203. doi:10.1007/s10439-010-9901-0. It is the position paper of the 6th International Bio-Fluid Mechanics Symposium (2008) and reviews progress from mid-2005 to early 2009, including a SWOT analysis. — [doi](https://doi.org/10.1007/s10439-010-9901-0); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3137382)
- Quarteroni A, Veneziani A, Vergara C (2016). "Geometric multiscale modeling of the cardiovascular system, between theory and practice." *CMAME* 302:193–252. doi:10.1016/j.cma.2016.01.007. It reviews 3D FSI, 1D and 0D (lumped) models and their coupling. — [doi](https://doi.org/10.1016/j.cma.2016.01.007)

**Outflow boundary conditions (resistance, RCR Windkessel, coupled 0D/1D)**
- Vignon-Clementel IE, Figueroa CA, Jansen KE, Taylor CA (2006). "Outflow boundary conditions for three-dimensional finite element modeling of blood flow and pressure in arteries." *CMAME* 195(29–32):3776–3796. doi:10.1016/j.cma.2005.04.014.
  - The paper derives outlet conditions for any downstream model in which pressure is an explicit function of flow or velocity at the interface (resistance, impedance, RCR, 1D).
  - It uses a Dirichlet-to-Neumann / variational multiscale approach within a stabilised semi-discrete FE formulation.
  - It notes that constant-pressure or prescribed-velocity outlets strongly distort solutions.
  - [doi](https://doi.org/10.1016/j.cma.2005.04.014); [SCOREC record](https://scorec.rpi.edu/reports/view_report.php?id=433)
- Kim HJ, Vignon-Clementel IE, Coogan JS, Figueroa CA, Jansen KE, Taylor CA (2010). "Patient-specific modeling of blood flow and pressure in human coronary arteries." *Annals of Biomedical Engineering* 38(10):3195–3209. doi:10.1007/s10439-010-0083-6. Coronary outlets are coupled to lumped models that include the ventricular pressure (closed-loop heart model, aorta, three-element Windkessels). This is the basis of FFR-CT boundary conditions. — [doi](https://doi.org/10.1007/s10439-010-0083-6)
- Pirola S, Cheng Z, Jarral OA, O'Regan DP, Pepper JR, Athanasiou T, Xu XY (2017). "On the choice of outlet boundary conditions for patient-specific analysis of aortic flow using computational fluid dynamics." *Journal of Biomechanics* 60:15–21. doi:10.1016/j.jbiomech.2017.06.005. Of five outlet strategies compared in a normal aorta, a well-tuned three-element Windkessel at all outlets gave the best agreement with in vivo flow and physiological pressure waveforms. Zero-pressure "open" outlets at multiple branches failed to reproduce physiological features. — [doi](https://doi.org/10.1016/j.jbiomech.2017.06.005)
- Chnafa C, Brina O, Pereira VM, Steinman DA (2018). "Better than nothing: a rational approach for minimizing the impact of outflow strategy on cerebrovascular simulations." *AJNR* 39(2):337–343. doi:10.3174/ajnr.A5484. With 70 MCA-aneurysm cases in 0D and 10 in 3D, outlet flow rates differed by up to 70% between outflow strategies. Murray-law and "splitting" methods gave near-physiological flows, and only splitting was insensitive to arbitrary truncation of the model. — [doi](https://doi.org/10.3174/ajnr.A5484)
- Murray CD (1926). "The physiological principle of minimum work: I. The vascular system and the cost of blood volume." *PNAS* 12(3):207–214. doi:10.1073/pnas.12.3.207. This is the cubic-diameter flow-split law used to distribute flow or resistance among outlets. — [doi](https://doi.org/10.1073/pnas.12.3.207)

**Backflow stabilisation**
- Esmaily Moghadam M, Bazilevs Y, Hsia T-Y, Vignon-Clementel IE, Marsden AL (2011). "A comparison of outlet boundary treatments for prevention of backflow divergence with relevance to blood flow simulations." *Computational Mechanics* 48(3):277–291. doi:10.1007/s00466-011-0599-0.
  - It compares three methods: (1) a boundary stabilisation term, (2) constraining outlet velocity to be normal, and (3) Lagrange-multiplier constraints on the outlet profile.
  - Stabilisation had the least effect on local and global haemodynamics, improved stability, and added little to no cost.
  - [doi](https://doi.org/10.1007/s00466-011-0599-0)

**Inflow profiles and rates**
- Morbiducci U, Ponzini R, Gallo D, Bignardi C, Rizzo G (2013). "Inflow boundary conditions for image-based computational hemodynamics: impact of idealized versus measured velocity profiles in the human aorta." *Journal of Biomechanics* 46(1):102–109. doi:10.1016/j.jbiomech.2012.10.012.
  - Measured axial (1D) PC-MRI profiles capture disturbed shear adequately.
  - Idealised (flat or parabolic) profiles distort bulk flow and especially helical structures in the aorta.
  - [doi](https://doi.org/10.1016/j.jbiomech.2012.10.012)
- Hoi Y, Wasserman BA, Xie YJ, Najjar SS, Ferrucci L, Lakatta EG, et al. (2010). "Characterization of volumetric flow rate waveforms at the carotid bifurcations of older adults." *Physiological Measurement* 31(3):291–302. doi:10.1088/0967-3334/31/3/002. It is a source of representative carotid inflow waveforms. Metadata was verified; the content was not read. — [doi](https://doi.org/10.1088/0967-3334/31/3/002)
- Womersley JR (1955). "Method for the calculation of velocity, rate of flow and viscous drag in arteries when the pressure gradient is known." *Journal of Physiology* 127(3):553–563. doi:10.1113/jphysiol.1955.sp005276. Its oscillatory pipe-flow solution is used both as an inlet profile and as an analytic verification case. — [doi](https://doi.org/10.1113/jphysiol.1955.sp005276); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1365740)

**Mesh, time step and numerics**
- Valen-Sendstad K, Steinman DA (2014). "Mind the gap: impact of computational fluid dynamics solution strategy on prediction of intracranial aneurysm hemodynamics and rupture status indicators." *AJNR* 35(3):536–543. doi:10.3174/ajnr.A3793.
  - It studied 12 MCA aneurysms at high versus normal resolution.
  - Commercial solvers' default low-order stabilisation adds artificial dissipation.
  - Typical strategies may suffice for highly reduced rupture indices but cannot be relied on for WSS magnitudes.
  - [doi](https://doi.org/10.3174/ajnr.A3793)
- Khan MO, Valen-Sendstad K, Steinman DA (2015). "Narrowing the expertise gap for predicting intracranial aneurysm hemodynamics: impact of solver numerics versus mesh and time-step resolution." *AJNR* 36(7):1310–1316. doi:10.3174/ajnr.A4263.
  - Solver numerics (discretisation scheme and settings) matter at least as much as mesh and time-step resolution.
  - With an optimised second-order, minimally dissipative solver, differences were negligible down to about 800,000 tetrahedra and about 5000 time steps per cardiac cycle.
  - [doi](https://doi.org/10.3174/ajnr.A4263); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7965263)

**Turbulence and transition in stenoses and aneurysms**
- Varghese SS, Frankel SH, Fischer PF (2007). "Direct numerical simulation of stenotic flows. Part 1. Steady flow." *Journal of Fluid Mechanics* 582:253–280. doi:10.1017/S0022112007005848. "Part 2. Pulsatile flow," *JFM* 582:281–318. doi:10.1017/S0022112007005836. They used spectral-element DNS of axisymmetric and eccentric 75% (area) stenoses, and a public DNS database exists at Purdue. — [Part 1](https://doi.org/10.1017/S0022112007005848); [Part 2](https://doi.org/10.1017/S0022112007005836); [database](https://engineering.purdue.edu/CFDLAB/DNS_stenoticflow_database/main.html)
- Lee SE, Lee S-W, Fischer PF, Bassiouny HS, Loth F (2008). "Direct numerical simulation of transitional flow in a stenosed carotid bifurcation." *Journal of Biomechanics* 41(11):2551–2561. doi:10.1016/j.jbiomech.2008.03.038.
  - It is a patient-specific spectral-element DNS.
  - Velocity fluctuations were at 100–300 Hz, with shedding strongest in systolic deceleration.
  - WSS in the stenosis reached 25–45 Pa.
  - One catalogue lists an extra author (C. D. Bertram); the DOI record lists five.
  - [doi](https://doi.org/10.1016/j.jbiomech.2008.03.038); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3279123)
- Valen-Sendstad K, Mardal K-A, Mortensen M, Reif BAP, Langtangen HP (2011). "Direct numerical simulation of transitional flow in a patient-specific intracranial aneurysm." *Journal of Biomechanics* 44(16):2826–2832. doi:10.1016/j.jbiomech.2011.08.015. Turbulent structures in an MCA aneurysm changed both the frequency of WSS direction changes and its magnitude, with a maximum of 41.5 Pa. — [doi](https://doi.org/10.1016/j.jbiomech.2011.08.015)

**Latest (2023–2026)**
- Taylor CA, Petersen K, Xiao N, Sinclair M, Bai Y, Lynch SR, Updegrove A, Schaap M (2023). "Patient-specific modeling of blood flow in the coronary arteries." *CMAME* 417(Part B):116414. doi:10.1016/j.cma.2023.116414. It is a recent review from the HeartFlow group. The DOI metadata misspells Updegrove as "UpdePac", and the content was not read. — [doi](https://doi.org/10.1016/j.cma.2023.116414)

### Inferences
- For a C++ solver, the outlet interface should be a pluggable "pressure as a function of flow history" operator, P = f(Q, state). Resistance, RCR and coronary models then become 0D ODE subsystems coupled implicitly or semi-implicitly. This is what Vignon-Clementel 2006, Kim 2010 and Pirola 2017 imply.
- Backflow stabilisation (Esmaily Moghadam 2011) is cheap and should be the default at every Neumann or Windkessel outlet.
- Default commercial-style first-order stabilised settings will suppress the high-frequency content that DNS studies find. A minimally dissipative second-order scheme is the safer reference.

### Gaps
- No single verified paper gave a quantitative mesh rule (for example, cells per diameter) for voxel/LBM grids in aneurysms; see §2 for the Randles group's 20 µm statement.
- Valen-Sendstad et al. (2015) on estimating inlet flow rates for aneurysm CFD was not verified, so it is omitted.

---

## 2. Lattice Boltzmann Methods for Blood in Vascular Geometries

### Takeaway
LBM suits voxelised segmentations because the fluid domain is a set of lattice sites taken directly from the image mask, with no unstructured meshing. It is also local and memory-bandwidth-bound, so it maps well to GPUs. Accuracy hinges on the wall treatment: simple bounce-back is first-order and staircased, while interpolated Bouzidi–Firdaouss–Lallemand (BFL) is second-order and cheap. Head-to-head comparisons with finite-volume solvers agree well in momentum-driven (MCA) aneurysms and less well in shear-driven (ICA) ones. Most major open LBM codes are GPL or AGPL; HemeLB (LGPL-3.0) is the main exception.

### Cited Findings

**Classical theory**
- Chen S, Doolen GD (1998). "Lattice Boltzmann method for fluid flows." *Annual Review of Fluid Mechanics* 30:329–364. doi:10.1146/annurev.fluid.30.1.329. It reviews the method for single- and multiphase flow and notes the ease of complex boundaries. — [doi](https://doi.org/10.1146/annurev.fluid.30.1.329)
- Krüger T, Kusumaatmaja H, Kuzmin A, Shardt O, Silva G, Viggen EM (2017). *The Lattice Boltzmann Method: Principles and Practice*. Springer, Graduate Texts in Physics. doi:10.1007/978-3-319-44649-3. The eBook was released 7 November 2016 with copyright 2017. The book has example code on GitHub. — [doi](https://doi.org/10.1007/978-3-319-44649-3); [Springer](https://www.springer.com/978-3-319-44647-9)
- d'Humières D, Ginzburg I, Krafczyk M, Lallemand P, Luo L-S (2002). "Multiple-relaxation-time lattice Boltzmann models in three dimensions." *Philosophical Transactions of the Royal Society A* 360(1792):437–451. doi:10.1098/rsta.2001.0955.
  - It gives D3Q15 and D3Q19 MRT models with better stability than BGK, tested on a 3D diagonally driven cavity at Re 500 and 2000.
  - The Crossref record lists only d'Humières; MEDLINE lists all five authors.
  - [doi](https://doi.org/10.1098/rsta.2001.0955)
- Ginzburg I, Verhaeghe F, d'Humières D (2008). "Two-relaxation-time lattice Boltzmann scheme: about parametrization, velocity, pressure and mixed boundary conditions." *Communications in Computational Physics* 3(2):427–478. doi:10.4208/cicp.2008.v3.p427.
  - It introduces the TRT "magic" parametrisation, under which permeabilities and boundary positions become viscosity-independent.
  - It gives second- and third-order velocity and pressure boundary schemes.
  - [doi](https://doi.org/10.4208/cicp.2008.v3.p427)
- Bouzidi M, Firdaouss M, Lallemand P (2001). "Momentum transfer of a Boltzmann-lattice fluid with boundaries." *Physics of Fluids* 13(11):3452–3459. doi:10.1063/1.1399290. It combines bounce-back with first- or second-order interpolation for curved walls, and computes forces by momentum exchange. OpenLB cites it as the source of its Bouzidi boundary. — [doi](https://doi.org/10.1063/1.1399290)
- Nash RW, Carver HB, Bernabeu MO, Hetherington J, Groen D, Krüger T, Coveney PV (2014). "Choice of boundary condition for lattice-Boltzmann simulation of moderate-Reynolds-number flow in complex domains." *Physical Review E* 89(2):023303. doi:10.1103/PhysRevE.89.023303.
  - In HemeLB on Poiseuille, Womersley and curved-pipe flow, BFL and Guo–Zheng–Shi converge at second order while simple bounce-back degrades to first order.
  - BFL beats GZS in unsteady flow and is much cheaper.
  - Junk–Yang is unstable at larger Re and is not recommended.
  - The arXiv v1 (2012) reported different conclusions, so cite the published version.
  - [doi](https://doi.org/10.1103/PhysRevE.89.023303); [arXiv](https://arxiv.org/abs/1211.0205)

**HemeLB and large-scale vascular LBM**
- Mazzeo MD, Coveney PV (2008). "HemeLB: A high performance parallel lattice-Boltzmann code for large scale fluid flow in complex geometries." *Computer Physics Communications* 178(12):894–914. doi:10.1016/j.cpc.2008.02.013. It is a C++ code built for sparse vascular geometries. — [doi](https://doi.org/10.1016/j.cpc.2008.02.013)
- Groen D, Hetherington J, Carver HB, Nash RW, Bernabeu MO, Coveney PV (2013). "Analysing and modelling the performance of the HemeLB lattice-Boltzmann simulation environment." *Journal of Computational Science* 4(5):412–422. doi:10.1016/j.jocs.2013.03.002. A later HemeLB performance study reports more than a 25× peak-performance gain between 2007 and 2014. — [doi](https://doi.org/10.1016/j.jocs.2013.03.002); [arXiv 1410.4713](https://arxiv.org/pdf/1410.4713)
- HemeLB licence: **LGPL-3.0**. The README says "open source under the LGPL license". — [GitHub hemelb-codes/hemelb](https://github.com/hemelb-codes/hemelb)
- HemeLB outlets: later HemeLB work imposes outlet pressure via density (ghost density) and drives it with two- or three-element Windkessels (P computed from Q). — [arXiv 2411.11797](https://arxiv.org/pdf/2411.11797); [arXiv 2202.11770](https://arxiv.org/pdf/2202.11770)
- HARVEY (Randles A, Draeger EW, Oppelstrup T, Krauss L, Gunnels JA). "Massively parallel models of the human circulatory system." SC15, Gordon Bell finalist.
  - It is a C/C++ LBM code that ran the systemic arterial tree on 1,572,864 Blue Gene/Q cores.
  - A companion paper states that a resolution of at least 20 µm is needed to study key macroscopic risk factors in patient-specific data.
  - The lab reports porting to GPU exascale systems (Frontier, Aurora).
  - [OSTI](https://www.osti.gov/servlets/purl/1241975); [Randles lab](https://randleslab.pratt.duke.edu/research/high-performance-computing)

**Rheology in LBM**
- Boyd J, Buick JM, Green S (2007). "Analysis of the Casson and Carreau-Yasuda non-Newtonian blood models in steady and oscillatory flows using the lattice Boltzmann method." *Physics of Fluids* 19(9):093103. doi:10.1063/1.2772250.
  - It is a 2D LBM study with shear-rate-dependent local relaxation time.
  - Non-Newtonian effects were significant in steady flow.
  - In oscillatory flow the differences were largest at low Re and low Womersley number, and larger for Casson.
  - [doi](https://doi.org/10.1063/1.2772250)

**Accuracy versus FV/FEM solvers**
- Jain K, Jiang J, Strother C, Mardal K-A (2016). "Transitional hemodynamics in intracranial aneurysms — comparative velocity investigations with high resolution lattice Boltzmann simulations, normal resolution ANSYS simulations, and MR imaging." *Medical Physics* 43(11):6186–6198. doi:10.1118/1.4964793.
  - In two basilar aneurysms, LBM DNS, ANSYS Fluent (laminar) and MR agreed well for the laminar case.
  - They differed clearly in the transitional case, but the differences became minimal after ensemble averaging and downsampling.
  - [doi](https://doi.org/10.1118/1.4964793)
- Jain K (2020). "Efficacy of the FDA nozzle benchmark and the lattice Boltzmann method for the analysis of biomedical flows in transitional regime." *Medical & Biological Engineering & Computing* 58(8):1817–1830. doi:10.1007/s11517-020-02188-8. LBM-based DNS at throat Re 2000 and 3500 gave acceptable agreement with PIV for jet breakdown location, mean velocity, shear stress and pressure. — [doi](https://doi.org/10.1007/s11517-020-02188-8); [arXiv](https://arxiv.org/abs/2005.07119)

**Open LBM codes and their licences** (relevant to an MIT-licensed infrastructure; checked 2026-11-04)

| Code | Licence | GPU | Source |
|---|---|---|---|
| HemeLB | **LGPL-3.0** | via HemePure-GPU | [GitHub](https://github.com/hemelb-codes/hemelb) |
| HemePure-GPU | **LGPLv3** (per Zacharoudiou 2023; the GitHub API was rate-limited, so this is not re-checked) | CUDA C++ | [Pure QUB](https://pure.qub.ac.uk/en/publications/development-and-performance-of-a-hemelb-gpu-code-for-human-scale-/) |
| Palabos | **AGPL-3.0** ("freely available under the terms of an open-source AGPLv3 license"; GitLab API agpl-3.0) | GPU-capable pre-release announced | [palabos.unige.ch](https://palabos.unige.ch/) |
| OpenLB | **GPL-2.0-or-later** (GitLab API; LICENSE is GPL v2) | GPU since 1.5 (2022-04, CUDA); AMD HIP since 1.9 (2025-12) | [openlb.net/download](https://www.openlb.net/download/) |
| waLBerla | **GPL-3.0-only** (GitLab API; website "GPLv3") | CUDA and HIP; code generation via lbmpy | [walberla.net](https://www.walberla.net/) |
| lbmpy | **GPL-3.0-or-later** (GitLab API) | code generation for CPU and GPU | i10git.cs.fau.de/pycodegen/lbmpy |
| TCLB | **GPL-3.0** (GitHub API) | CUDA | github.com/CFD-GO/TCLB |
| FluidX3D | **Custom, non-OSI**: free for "public research, education or personal use"; no commercial use, no military use, no training AI models on the code | OpenCL; FP32 arithmetic with FP16S/FP16C storage | [LICENSE.md](https://github.com/ProjectPhysX/FluidX3D) |

- FluidX3D's README lists SRT/BGK (default) and TRT collision, mid-grid bounce-back (stationary and moving) and equilibrium inflow/outflow boundaries. It says FP16 storage nearly halves memory and nearly doubles speed "without noticeably affecting accuracy for most setups". It does not mention interpolated bounce-back. — [GitHub](https://github.com/ProjectPhysX/FluidX3D)

**Latest (2023–2026)**
- Zacharoudiou I, McCullough JWS, Coveney PV (2023). "Development and performance of a HemeLB GPU code for human-scale blood flow simulation." *Computer Physics Communications* 282:108548. doi:10.1016/j.cpc.2022.108548.
  - It is a CUDA C++ port that keeps the CPU code's strong scaling, with scaling analysis "up to tens of thousands of GPUs".
  - The code is released as HemePure-GPU under LGPLv3.
  - [doi](https://doi.org/10.1016/j.cpc.2022.108548); [arXiv](https://arxiv.org/abs/2202.11770)
- Benemerito I, McCullough J, Narracott A, Coveney PV, Marzo A (2025). "Comparison of Navier-Stokes and lattice Boltzmann solvers for subject-specific modelling of intracranial aneurysms." *Computers in Biology and Medicine* 197(Part B):111050. doi:10.1016/j.compbiomed.2025.111050.
  - It compares Ansys CFX (FV) with HemeLB (LBM) on ICA and MCA aneurysms across 462 boundary-condition profiles.
  - The MCA aneurysm (momentum-driven) showed good quantitative agreement for velocity and WSS metrics.
  - The ICA aneurysm (shear-driven) varied more.
  - Both identified the same locations of extreme rupture-risk metrics.
  - [doi](https://doi.org/10.1016/j.compbiomed.2025.111050); [Pure QUB](https://pure.qub.ac.uk/en/publications/comparison-of-navier-stokes-and-lattice-boltzmann-solvers-for-sub/)
- Huang F, Noël R, Berg P, Hosseini SA (2022 preprint). "Simulation of the FDA nozzle benchmark: a lattice Boltzmann study." arXiv:2204.10566. The journal version is **[unverified]**. — [arXiv](https://arxiv.org/abs/2204.10566)

### Inferences
- **Why LBM suits voxelised segmentations.** A binary or level-set mask from CT/MRI maps one-to-one to fluid, wall and solid lattice sites. Sub-voxel wall distance q, taken from a signed distance field or the level set, feeds BFL interpolation directly. No surface-to-volume meshing step (TetGen/Gmsh) is needed, which also removes the AGPL TetGen dependency (§8). The method is explicit and local, so it maps naturally onto GPU kernels and a task-graph scheduler.
- **Accuracy recipe for a C++ LBM.** Use D3Q19 with TRT (magic parameter 1/4) or MRT, BFL walls, and Windkessel pressure outlets as in HemeLB. Verify against Poiseuille and Womersley (Nash 2014 test set) and the FDA nozzle (Jain 2020).
- **Known LBM caveats** (standard results in Krüger et al. 2017, not quoted from a fetched page):
  - LBM is weakly compressible, with compressibility error of order Ma².
  - Viscosity is tied to the relaxation time, so very low viscosity or high Re needs fine lattices or regularised/MRT collision.
  - Pulsatile flow needs small time steps because Δt ∝ Δx² at fixed τ.
  - Plan diffusive scaling for convergence studies.
- **Licensing.** Only HemeLB (LGPL) can be linked from MIT code with modest obligations. Palabos, OpenLB, waLBerla, lbmpy and TCLB (GPL/AGPL) can be used as external cross-validation tools at a process boundary, but cannot be incorporated. FluidX3D cannot be used for anything commercial.

### Gaps
- The search found no peer-reviewed 2024–2025 report of interactive, minutes-scale GPU-LBM simulation of a patient-specific aneurysm.
- The quantitative error metrics in Benemerito 2025 were not read (abstract only).
- HARVEY's public availability and licence were not found.
- The Sailfish GPU LBM repository was not found at its former GitHub location, so its licence is unknown.

---

## 3. Blood Rheology: Models, Parameters, Haematocrit, Fåhræus–Lindqvist, and When It Matters

### Takeaway
Blood is shear-thinning, viscoelastic and thixotropic, but in large arteries (shear rates above about 100 s⁻¹) a Newtonian viscosity of about 3.5 mPa·s is usually adequate. Several high-resolution aneurysm and coronary studies find rheology to be a second-order effect next to segmentation, inflow and numerics. Non-Newtonian models matter more in low-shear regions (aneurysm domes, recirculation, veins, the venous IVC), for OSI and for transport. The answer is parameter-set dependent: different published Carreau–Yasuda fits give materially different OSI.

### Cited Findings

**Classical**
- Cho YI, Kensey KR (1991). "Effects of the non-Newtonian viscosity of blood on flows in a diseased arterial vessel. Part 1: Steady flows." *Biorheology* 28(3–4):241–262. doi:10.3233/BIR-1991-283-415.
  - It is an FE study in a coronary casting comparing several constitutive models.
  - Non-Newtonian treatment notably changed pressure drop at Re of about 100 or below.
  - It gave higher WSS in separation and recirculation zones.
  - [doi](https://doi.org/10.3233/BIR-1991-283-415)
- Cho–Kensey Carreau parameters, widely reproduced: μ₀ = 0.056 Pa·s, μ∞ = 0.00345 Pa·s, λ = 3.313 s, n = 0.3568 — **[unverified against the original paper]**. The values appear in many secondary sources, but no fetched page confirmed the mapping. A 2011 conference paper lists n ≈ 0.357, μ₀ ≈ 0.056, μ∞ ≈ 0.00345, λ ≈ 3.313 s without attributing them. — [search context](https://bura.brunel.ac.uk/bitstream/2438/6835/1/MNF2011.pdf)
- A widely used Carreau–Yasuda set (verified as used in Assi et al. 2024): μ∞ = 0.0035 Pa·s, μ₀ = 0.16 Pa·s, λ = 8.2 s, n = 0.2128, a = 0.64. It is commonly attributed to Gijsen et al.; that attribution was not verified here. — [Assi 2024](https://doi.org/10.3389/fbioe.2024.1302063)
- Johnston BM, Johnston PR, Corney S, Kilpatrick D (2004). "Non-Newtonian blood flow in human right coronary arteries: steady state simulations." *Journal of Biomechanics* 37(5):709–720. doi:10.1016/j.jbiomech.2003.09.016. It compares five non-Newtonian models with Newtonian in four reconstructed RCAs. Later citing work summarises it as showing little difference in WSS patterns. — [doi](https://doi.org/10.1016/j.jbiomech.2003.09.016)
- Quemada D (1978). "Rheology of concentrated disperse systems III. General features of the proposed non-Newtonian model. Comparison with experimental data." *Rheologica Acta* 17(6):643–653. doi:10.1007/BF01522037. Part II (the model): *Rheol Acta* 17:632–642 (DOI not verified).
  - The model is η = η_p(1 − ½ k(γ̇) H)⁻², with H the haematocrit.
  - k is an intrinsic viscosity interpolating between k₀ and k∞ via a critical shear rate γ̇_c (exponent p = 1/2 for blood).
  - Check the exact exponent placement against the PDF.
  - [doi](https://doi.org/10.1007/BF01522037)
- Fåhræus R, Lindqvist T (1931). "The viscosity of the blood in narrow capillary tubes." *American Journal of Physiology* 96:562–568 (no DOI checked). Apparent viscosity falls with tube diameter below roughly 0.3 mm, and tube haematocrit is lower than feed haematocrit. — [Wikipedia summary](https://en.wikipedia.org/wiki/F%C3%A5hr%C3%A6us%E2%80%93Lindqvist_effect)
- Pries AR, Neuhaus D, Gaehtgens P (1992). "Blood viscosity in tube flow: dependence on diameter and hematocrit." *American Journal of Physiology* 263(6):H1770–H1778 (no DOI checked).
  - It gives empirical fits of relative apparent viscosity as a function of tube diameter and haematocrit.
  - The viscosity–haematocrit relation is linear for diameters of 6 µm or less and over-proportional for 9 µm or more.
  - [index](https://data.zhenzhi.com/details/1/999364089)

**Standard constitutive forms** (no parameter values implied)
- Carreau–Yasuda: μ(γ̇) = μ∞ + (μ₀ − μ∞)[1 + (λγ̇)^a]^((n−1)/a). Carreau is a = 2.
- Cross: μ = μ∞ + (μ₀ − μ∞)/(1 + (λγ̇)^m).
- Casson: √τ = √τ_y + √(η_c γ̇), which has a yield stress.
- In each case γ̇ = √(2 D:D) from the strain-rate tensor D.

**When non-Newtonian matters**
- Khan MO, Steinman DA, Valen-Sendstad K (2017). "Non-Newtonian versus numerical rheology: practical impact of shear-thinning on the prediction of stable and unstable flows in intracranial aneurysms." *International Journal for Numerical Methods in Biomedical Engineering* 33(7):e2836. doi:10.1002/cnm.2836.
  - With a modified-Cross model, time-averaged and instantaneous WSS, the spectral content of instabilities and dome-averaged metrics were minimally affected.
  - Solution strategy mattered more.
  - Point-wise normalisation of non-Newtonian by Newtonian WSS artificially amplifies differences in low-WSS regions.
  - [doi](https://doi.org/10.1002/cnm.2836)
- Kandangwa P, Torii R, Gatehouse PD, Sherwin SJ, Weinberg PD (2022). "Influence of right coronary artery motion, flow pulsatility and non-Newtonian rheology on wall shear stress metrics." *Frontiers in Bioengineering and Biotechnology* 10:962687. doi:10.3389/fbioe.2022.962687.
  - With Carreau–Yasuda using four published parameter sets, the Cho & Kensey set gave Spearman correlations of 0.99 with Newtonian for all metrics.
  - The Gijsen set gave a minimum correlation of 0.95 (OSI) and point-wise OSI differences from about +25% to −250%.
  - Non-Newtonian effects were "important" for two of the four parameter sets.
  - [doi](https://doi.org/10.3389/fbioe.2022.962687)

**Latest (2023–2026)**
- Lynch S, Nama N, Figueroa CA (2022, late-2022 boundary). "Effects of non-Newtonian viscosity on arterial and venous flow and transport." *Scientific Reports* 12:20568. doi:10.1038/s41598-022-19867-1.
  - Power-law and Carreau–Yasuda models were implemented in CRIMSON for a thoracic aortic aneurysm and an IVC–iliac venous model.
  - Newtonian gave larger in-plane velocity and vorticity, while non-Newtonian gave larger WSS in both arterial and venous models.
  - Newtonian overestimated near-wall and luminal mass transport.
  - [doi](https://doi.org/10.1038/s41598-022-19867-1); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9709089)
- Sasmal C (2026, arXiv:2609.20855, submitted 31 Aug 2026). "How does human blood rheology influence arterial and cardiovascular device hemodynamics? A comprehensive review, discussion, and future directions." It reviews generalised Newtonian, viscoelastic and thixo-elasto-viscoplastic models. Only the abstract was read. — [arXiv](https://arxiv.org/abs/2609.20855)

### Inferences
- Implement rheology as a pointwise viscosity law μ(γ̇, H), a strategy object, so that Newtonian, Carreau–Yasuda, Cross, Casson (regularised) and Quemada with haematocrit are interchangeable.
  - In LBM this sets a local relaxation time; clamp τ to keep it stable.
  - In FEM it enters the viscous term with Picard or Newton linearisation.
- Report results for both the Newtonian and the Cho–Kensey set, and avoid point-wise ratios in low-WSS zones (Khan 2017).
- Fåhræus–Lindqvist and Pries fits matter only for vessels below about 0.3 mm, which is outside CT/MRI-resolvable 3D geometry. They belong in 0D/1D downstream models, not in the 3D domain.

### Gaps
- The Cho–Kensey parameter values and the Gijsen attribution were not verified against the original texts.
- No verified, quantitative guideline (for example, a threshold Womersley or shear rate) was found for when non-Newtonian models change clinical indices.

---

## 4. Fluid–Structure Interaction: Methods, Wall Models, and When FSI Changes WSS

### Takeaway
Two families dominate:
- **Coupled momentum method (CMM).** A fixed fluid mesh plus a membrane wall model gives cheap, robust, small-deformation FSI.
- **Arbitrary Lagrangian–Eulerian (ALE) FSI.** Full nonlinear solids, often monolithic, are needed for large motions such as a dissection flap.

Partitioned schemes suffer the added-mass instability for blood-like density ratios. Rigid-wall models usually reproduce high-WSS regions and global pressures and flows. Compliance changes low and oscillatory WSS, OSI and recirculation, sometimes enormously (AnXplore: sac OSI −36% to +674%). Mobile structures (dissection flaps) need FSI.

### Cited Findings

**Classical**
- Figueroa CA, Vignon-Clementel IE, Jansen KE, Hughes TJR, Taylor CA (2006). "A coupled momentum method for modeling blood flow in three-dimensional deformable arteries." *CMAME* 195(41–43):5685–5706. doi:10.1016/j.cma.2005.11.011.
  - The fluid and wall are strongly coupled, with the wall a linear membrane enhanced with transverse shear.
  - It is motivated by the difficulty of ALE on large realistic models and is described as "remarkably robust".
  - [doi](https://doi.org/10.1016/j.cma.2005.11.011)
- Bazilevs Y, Calo VM, Zhang Y, Hughes TJR (2006). "Isogeometric fluid–structure interaction analysis with applications to arterial blood flow." *Computational Mechanics* 38(4–5):310–322. doi:10.1007/s00466-006-0084-3. It couples a NURBS-based incompressible fluid with a nonlinear elastic solid at large displacements, and builds analysis-suitable NURBS from patient images. — [doi](https://doi.org/10.1007/s00466-006-0084-3)
- Bazilevs Y, Calo VM, Hughes TJR, Zhang Y (2008). "Isogeometric fluid-structure interaction: theory, algorithms, and computations." *Computational Mechanics* 43(1):3–37. doi:10.1007/s00466-008-0315-x. It presents a fully coupled monolithic ALE formulation with a hyperelastic solid, including the tangent matrix with shape derivatives. — [doi](https://doi.org/10.1007/s00466-008-0315-x)
- Causin P, Gerbeau J-F, Nobile F (2005). "Added-mass effect in the design of partitioned algorithms for fluid–structure problems." *CMAME* 194(42–44):4506–4527. doi:10.1016/j.cma.2004.12.005. It explains why loosely coupled schemes go unstable for blood-flow density ratios and thin walls, and why strongly coupled schemes then need many subiterations. — [doi](https://doi.org/10.1016/j.cma.2004.12.005)
- Holzapfel GA, Gasser TC, Ogden RW (2000). "A new constitutive framework for arterial wall mechanics and a comparative study of material models." *Journal of Elasticity* 61(1–3):1–48. doi:10.1023/A:1010835316564. It is the HGO two-layer, fibre-reinforced hyperelastic model with residual stress via an opened sector, using three parameters per layer in its specific form. — [doi](https://doi.org/10.1023/A:1010835316564)
- Brown AG, Shi Y, Marzo A, Staicu C, Valverde I, Beerbaum P, et al. (2012). "Accuracy vs. computational time: translating aortic simulations to the clinic." *Journal of Biomechanics* 45(3):516–523. doi:10.1016/j.jbiomech.2011.11.041.
  - FSI took longest to reach periodicity.
  - Boundary pressures and flow rates were comparable across rigid and FSI approaches.
  - The WSS comparison was not read.
  - [doi](https://doi.org/10.1016/j.jbiomech.2011.11.041)
- Alimohammadi M, Sherwood JM, Karimpour M, Agu O, Balabani S, Díaz-Zuccarini V (2015). "Aortic dissection simulation models for clinical support: fluid-structure interaction vs. rigid wall models." *BioMedical Engineering OnLine* 14:34. doi:10.1186/s12938-015-0032-6.
  - In one type-B dissection (two-way FSI, Windkessel outlets, Carreau–Yasuda, SST turbulence), high-WSS regions changed little.
  - Some low and oscillatory WSS regions appeared only with FSI.
  - The authors judge FSI justified for interventional planning.
  - [doi](https://doi.org/10.1186/s12938-015-0032-6)
- Bäumler K, Vedula V, Sailer AM, Seo J, Chiu P, Mistelbauer G, et al. (2020). "Fluid–structure interaction simulations of patient-specific aortic dissection." *Biomechanics and Modeling in Mechanobiology* 19(5):1607–1628. doi:10.1007/s10237-020-01294-8.
  - It uses two-way FSI with prestress, external tissue support and intercostal/lumbar tethering, with boundary conditions from 4D flow MRI, and agrees with 4D flow.
  - Lowering flap stiffness increased flap displacement from 1.4 to 13.4 mm, reduced the high-TAWSS area, and cut abdominal true-lumen flow by up to 20%.
  - [doi](https://doi.org/10.1007/s10237-020-01294-8)

**Latest (2023–2026)**
- Goetz A, Jeken-Rico P, Pelissier U, Chau Y, Sédat J, Hachem E (2024). "AnXplore: a comprehensive fluid-structure interaction study of 101 intracranial aneurysms." *Frontiers in Bioengineering and Biotechnology* 12:1433811. doi:10.3389/fbioe.2024.1433811.
  - 101 real aneurysm sacs were mounted on one idealised toroidal artery and simulated with rigid and compliant walls.
  - Sac-averaged OSI changed by −36% to +674% relative to the rigid wall, and near-dome recirculation was notably affected.
  - In a flow-diverter case, systolic sac velocity was 73% higher with the compliant wall.
  - The dataset and parameters are public.
  - [doi](https://doi.org/10.3389/fbioe.2024.1433811); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11243300)
- svMultiPhysics (SimVascular): a C++ re-implementation of the Fortran svFSI solver. It is an FE multiphysics code for fluid, solid, FSI (ALE and CMM in the svFSI lineage), diffusion and electrophysiology, with PETSc/Trilinos backends. Licence **BSD-3-Clause**. — [GitHub](https://github.com/SimVascular/svMultiPhysics)
- Zhu C, Vedula V, Parker D, Wilson N, Shadden S, Marsden A (2022). "svFSI: a multiphysics package for integrated cardiac modeling." *Journal of Open Source Software* 7(78):4118. doi:10.21105/joss.04118. — [doi](https://doi.org/10.21105/joss.04118)

### Inferences
- **Staged FSI path for a C++ code.** Rigid-wall first; then CMM, a small-strain membrane on a fixed mesh that is cheap and fits both FEM and, in principle, LBM via moving-wall bounce-back; then full ALE or immersed methods only for flaps and valves.
- **Partitioned coupling** needs Robin–Neumann or other added-mass-aware interface conditions, or a monolithic solve (Causin 2005).
- **Wall-sensitive indices.** OSI, RRT and low-WSS area are the most sensitive to wall motion (AnXplore, Alimohammadi). TAWSS and pressure drop are least sensitive. Report which indices were computed with rigid walls.

### Gaps
- No systematic, multi-patient quantitative meta-analysis of rigid-versus-FSI WSS differences was found. The evidence is case series (AnXplore is the largest).
- Prestress algorithms (for example, Hsu & Bazilevs 2011) were not verified here.

---

## 5. Haemodynamic Indices: Definitions and Formulas

### Takeaway
Most studies use TAWSS, OSI, RRT, transWSS, ECAP and helicity descriptors, all functionals of the wall shear stress vector τ_w(x, t) over a cycle T. TAWSS, OSI and RRT are strongly correlated, and RRT has been proposed as a single robust low-and-oscillatory metric. transWSS captures multidirectional flow that OSI misses. Indices are sensitive to WSS computation at the wall (gradient recovery), so normalise by a parent-artery value when comparing groups.

### Cited Findings

**Definitions**
- WSS vector: τ_w = μ (∇u + ∇uᵀ)·n projected tangentially, τ_w = σ_v·n − ((σ_v·n)·n) n.
- TAWSS = (1/T) ∫₀ᵀ |τ_w| dt.
- OSI = ½ [1 − |∫₀ᵀ τ_w dt| / ∫₀ᵀ |τ_w| dt], which ranges from 0 (unidirectional) to 0.5 (purely oscillatory).
- RRT = 1/[(1 − 2·OSI)·TAWSS] = 1/|(1/T)∫₀ᵀ τ_w dt|.
- These standard forms are as given in recent secondary sources (the AortaCFD tutorial attributes RRT to Himburg 2004, and a PMC table states the identity). — [AortaCFD hemodynamics](https://jiewangnk.github.io/AortaCFD-web/hemodynamics/); [PMC table](https://pmc.ncbi.nlm.nih.gov/articles/PMC8365882/table/T2)
- Ku DN, Giddens DP, Zarins CK, Glagov S (1985). "Pulsatile flow and atherosclerosis in the human carotid bifurcation. Positive correlation between plaque location and low oscillating shear stress." *Arteriosclerosis* 5(3):293–302. doi:10.1161/01.ATV.5.3.293. It is the origin of OSI. In a laser-Doppler model, the outer sinus wall had a mean WSS of −0.5 dyn/cm² oscillating between −7 and +4 dyn/cm², co-located with plaque. — [doi](https://doi.org/10.1161/01.ATV.5.3.293)
- He X, Ku DN (1996). "Pulsatile flow in the human left coronary artery bifurcation: average conditions." *Journal of Biomechanical Engineering* 118(1):74–82. doi:10.1115/1.2795948. It is the commonly cited form of the OSI definition, though the formula was not seen in the fetched text. — [doi](https://doi.org/10.1115/1.2795948)
- Himburg HA, Grzybowski DM, Hazel AL, LaMack JA, Li X-M, Friedman MH (2004). "Spatial comparison between wall shear stress measures and porcine arterial endothelial permeability." *American Journal of Physiology – Heart and Circulatory Physiology* 286(5):H1916–H1922. doi:10.1152/ajpheart.00897.2003.
  - In vivo permeability fell with increasing TAWSS and rose slightly with OSI.
  - RRT is usually attributed to this paper, although RRT was not seen in the abstract.
  - [doi](https://doi.org/10.1152/ajpheart.00897.2003)
- Lee S-W, Antiga L, Steinman DA (2009). "Correlations among indicators of disturbed flow at the normal carotid bifurcation." *Journal of Biomechanical Engineering* 131(6):061013. doi:10.1115/1.3127252.
  - In 50 MRI-based normal carotids, TAWSS, OSI and RRT were strongly inter-correlated.
  - RRT is recommended as a robust single metric of low and oscillating shear (per an abstract summary).
  - [doi](https://doi.org/10.1115/1.3127252)
- Peiffer V, Sherwin SJ, Weinberg PD (2013). "Computation in the rabbit aorta of a new metric – the transverse wall shear stress – to quantify the multidirectional character of disturbed blood flow." *Journal of Biomechanics* 46(15):2651–2658. doi:10.1016/j.jbiomech.2013.08.003.
  - transWSS = (1/T)∫₀ᵀ |τ_w · (n × τ̄/|τ̄|)| dt, the cycle average of WSS components perpendicular to the mean WSS vector τ̄.
  - Its patterns differ from those of TAWSS and OSI.
  - [doi](https://doi.org/10.1016/j.jbiomech.2013.08.003)
- Di Achille P, Tellides G, Figueroa CA, Humphrey JD (2014). "A haemodynamic predictor of intraluminal thrombus formation in abdominal aortic aneurysms." *Proceedings of the Royal Society A* 470(2172):20140163. doi:10.1098/rspa.2014.0163. It is the source usually cited for ECAP = OSI/TAWSS (endothelial cell activation potential). The bibliographic record was verified; the formula was not seen in a fetched text. — [doi](https://doi.org/10.1098/rspa.2014.0163)
- Helicity:
  - The helicity density is H_k = u·ω, with vorticity ω = ∇×u.
  - The local normalised helicity is LNH = u·ω/(|u||ω|), in [−1, 1].
  - Morbiducci et al. (2009) quantified aortic helical flow from 4D PC-MRI with a Lagrangian Helical Flow Index (not LNH): Morbiducci U, Ponzini R, Rizzo G, Cadioli M, Esposito A, De Cobelli F, et al. "In vivo quantification of helical blood flow in human aorta by time-resolved three-dimensional cine phase contrast magnetic resonance imaging." *Annals of Biomedical Engineering* 37(3):516–531. doi:10.1007/s10439-008-9609-6.
  - LNH is cited in later work as a standard helicity descriptor.
  - [doi](https://doi.org/10.1007/s10439-008-9609-6); [PMC5316512](https://pmc.ncbi.nlm.nih.gov/articles/PMC5316512)
- Gallo D, Steinman DA, Bijari PB, Morbiducci U (2012). "Helical flow in carotid bifurcation as surrogate marker of exposure to disturbed shear." *Journal of Biomechanics* 45(14):2398–2404. doi:10.1016/j.jbiomech.2012.07.007. Helicity-based bulk-flow descriptors, selected by AIC in regression, predict exposure to disturbed WSS. — [doi](https://doi.org/10.1016/j.jbiomech.2012.07.007)
- Low-TAWSS masking: some tools mask OSI where TAWSS is below 0.5 Pa because OSI is numerically unreliable there. — [AortaCFD](https://jiewangnk.github.io/AortaCFD-web/hemodynamics/)

### Inferences
- Compute τ_w from a recovered gradient: superconvergent patch recovery in FEM, or the non-equilibrium stress tensor in LBM, −(1 − 1/(2τ)) Σᵢ cᵢcᵢ fᵢ^neq. Then accumulate running time integrals of τ_w and |τ_w| per wall node. All of TAWSS, OSI, RRT, transWSS and ECAP follow from three accumulators (∫τ_w, ∫|τ_w| and a second pass for transWSS), which allows streaming computation on a GPU.
- Normalise sac or region WSS by a parent-vessel reference before comparing across teams or patients (Valen-Sendstad 2018).

### Gaps
- Exact definitions of Gallo's h1–h4 helicity descriptors were not retrieved.
- The original OSI formula text in He & Ku (1996) and the ECAP formula in Di Achille (2014) were not seen in fetched full text; the standard forms are given.

---

## 6. Validation, Verification and Reproducibility

### Takeaway
Verification uses analytic solutions (Poiseuille, Womersley 1955) and DNS databases. Validation uses the FDA benchmark nozzle (multi-lab PIV, throat Re 500–6500, public data) and the FDA blood pump. The FDA nozzle's transitional cases (Re 2000–3500) separate low-dissipation solvers from the rest. Multi-team aneurysm challenges show that pressure is robust across solvers, but flow patterns, WSS magnitudes, segmentation (volume ±20%, ostium ±30%) and inflow choices dominate real-world variability. Credibility in a regulatory context follows ASME V&V 40-2018 and the FDA's November 2023 guidance.

### Cited Findings

**Analytic and DNS verification**
- Womersley (1955), §1, gives the oscillatory pipe flow. The Nash et al. (2014) LBM test set is Poiseuille, Womersley and steady curved-pipe flow, and the Purdue DNS stenosis database is in §1. — [Womersley](https://doi.org/10.1113/jphysiol.1955.sp005276); [Nash](https://doi.org/10.1103/PhysRevE.89.023303)

**FDA benchmark nozzle and blood pump (classical)**
- Hariharan P, Giarra M, Reddy V, Day SW, Manning KB, Deutsch S, et al. (2011). "Multilaboratory particle image velocimetry analysis of the FDA benchmark nozzle model to support validation of computational fluid dynamics simulations." *Journal of Biomechanical Engineering* 133(4):041002. doi:10.1115/1.4003440.
  - Three laboratories measured PIV and pressure at throat Re 500–6500 (laminar, transitional, turbulent).
  - Interlaboratory uncertainty was about 10% at most locations for the laminar and turbulent cases.
  - It reached about 60% for jet size and velocity in the transitional Re 2000 case.
  - [doi](https://doi.org/10.1115/1.4003440); [PDF on NCIHub](https://ncihub.cancer.gov/app/site/collections/2999/Hariharan_PIV_Analysis_of_Nozzle_2011.pdf)
- Stewart SFC, Paterson EG, Burgreen GW, Hariharan P, Giarra M, Reddy V, et al. (2012). "Assessment of CFD performance in simulations of an idealized medical device: results of FDA's first computational interlaboratory study." *Cardiovascular Engineering and Technology* 3(2):139–160. doi:10.1007/s13239-012-0087-5.
  - The submissions' results varied widely.
  - CFD used for device safety therefore requires careful experimental validation.
  - The submissions used laminar or RANS models; later work points to LES and transition modelling.
  - [doi](https://doi.org/10.1007/s13239-012-0087-5)
- Malinauskas RA, Hariharan P, Day SW, Herbertson LH, Buesen M, Steinseifer U, et al. (2017). "FDA benchmark medical device flow models for CFD validation." *ASAIO Journal* 63(2):150–160. doi:10.1097/MAT.0000000000000499.
  - It covers the nozzle and a centrifugal blood pump, with more than 20 CFD groups.
  - Discrepancies concentrated in separated flow.
  - 57% of pump pressure-head predictions fell within 1 SD of measurement.
  - Only 37% of 52 submissions reported haemolysis.
  - [doi](https://doi.org/10.1097/MAT.0000000000000499)
- Raben JS, Hariharan P, Robinson R, Malinauskas R, Vlachos PP (2014 preprint). "Time-resolved particle image velocimetry measurements with wall shear stress and uncertainty quantification for the FDA benchmark nozzle model." arXiv:1405.3125. It reports time-resolved PIV with WSS and uncertainty at Re 500, 2000, 5000 and 8000. The journal version is **[unverified]**. — [arXiv](https://arxiv.org/abs/1405.3125)
- Fehn N, Wall WA, Kronbichler M (2018 preprint). "Modern discontinuous Galerkin methods for the simulation of transitional and turbulent flows in biomedical engineering: a comprehensive LES study of the FDA benchmark nozzle model." arXiv:1811.06752. It finds jet-breakdown location to be the main challenge and to be sensitive to numerics. The journal version is **[unverified]**. — [arXiv](https://arxiv.org/abs/1811.06752)
- Data access:
  - The FDA CFD wiki (ncihub.org/wiki/FDA_CFD; pump at nciphub.org/wiki/FDA_CFD/ComputationalRoundRobin2Pump) hosts geometry, CAD and PIV data. The ncihub.cancer.gov host was unreachable from here on 2026-11-04.
  - The FDA Regulatory Science Tool catalogue lists the nozzle and pump validation data (pressure, velocity), CAD and descriptions on GitHub, with Zenodo DOI 10.5281/zenodo.6629049.
  - A separate RST entry gives haemolysis (plasma-free haemoglobin) data, dated 05/05/2026.
  - [FDA RST search result](https://cdrh-rst.fda.gov/node/1036); [Zenodo](https://zenodo.org/records/6629049/files/CFD%20Tool.pdf)

**Cerebral aneurysm CFD challenges**
- Steinman DA, Hoi Y, Fahy P, Morris L, Walsh MT, Aristokleous N, et al. (56 authors) (2013). "Variability of computational fluid dynamics solutions for pressure and flow in a giant aneurysm: the ASME 2012 Summer Bioengineering Conference CFD Challenge." *Journal of Biomechanical Engineering* 135(2):021016. doi:10.1115/1.4023382.
  - The article number is 021016 per the DOI record; one index lists 021015.
  - In Phase I, 25 groups mostly agreed on peak-systolic pressure drop within 8%.
  - In Phase II, median predictions came within a few mmHg of measurement once sub-millimetre reconstruction differences were accounted for.
  - Sac flow patterns varied, and only a few groups with high temporal resolution captured instabilities.
  - [doi](https://doi.org/10.1115/1.4023382)
- Janiga G, Berg P, Sugiyama S, Kono K, Steinman DA (2015). "The Computational Fluid Dynamics Rupture Challenge 2013—Phase I: prediction of rupture status in intracranial aneurysms." *AJNR* 36(3):530–536. doi:10.3174/ajnr.A4157.
  - 21 of 26 groups (81%) identified which of two similar MCA aneurysms had ruptured, against 39 of 43 neurosurgeons (91%).
  - Nobody located the rupture site.
  - [doi](https://doi.org/10.3174/ajnr.A4157)
- Berg P, Roloff C, Beuing O, Voss S, Sugiyama S, Aristokleous N, et al. (50 authors) (2015). "The Computational Fluid Dynamics Rupture Challenge 2013—Phase II: variability of hemodynamic simulations in two intracranial aneurysms." *Journal of Biomechanical Engineering* 137(12):121008. doi:10.1115/1.4031794.
  - 26 groups from 15 countries took part, and about 80% produced similar results.
  - Only two cases deviated by more than 35% from the mean.
  - Steady CFD agreed with PIV in a silicone model.
  - [doi](https://doi.org/10.1115/1.4031794)
- Berg P, Voß S, Saalfeld S, Janiga G, Bergersen AW, Valen-Sendstad K, et al. (51 authors) (2018). "Multiple Aneurysms AnaTomy CHallenge 2018 (MATCH): Phase I: Segmentation." *Cardiovascular Engineering and Technology* 9(4):565–581. doi:10.1007/s13239-018-00376-0.
  - 26 groups segmented 3D-DSA of five aneurysms, with no consensus on software or algorithms.
  - Aneurysm volume differed by up to 20% and ostium area by up to 30%.
  - Surfaces were systematically over-segmented by about 10% against a 2D reference, and morphology parameters varied by up to 25%.
  - [doi](https://doi.org/10.1007/s13239-018-00376-0)
- Valen-Sendstad K, Bergersen AW, Shimogonya Y, Goubergrits L, Bruening J, Pallares J, et al. (55 authors) (2018). "Real-world variability in the prediction of intracranial aneurysm wall shear stress: the 2015 International Aneurysm CFD Challenge." *Cardiovascular Engineering and Technology* 9(4):544–564. doi:10.1007/s13239-018-00374-2.
  - 26 teams submitted 28 datasets covering the whole pipeline.
  - Segmentation, model extent and inflow rate drove a large spread in sac-average WSS, which normalising by parent-artery WSS reduced substantially.
  - Maximum WSS and low-shear area were more variable, and team experience did not predict variability.
  - The authors recommend standardising extents, inflow rates and blood properties.
  - [doi](https://doi.org/10.1007/s13239-018-00374-2); [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6290689/)

**Credibility frameworks**
- ASME V&V 40-2018, "Assessing Credibility of Computational Modeling through Verification and Validation: Application to Medical Devices". It is a risk-informed framework covering the question of interest, the context of use, model risk and credibility goals. — [UMN summary](https://dmd.umn.edu/2023/vv-40)
- FDA final guidance, "Assessing the Credibility of Computational Modeling and Simulation in Medical Device Submissions". It was issued on 16 November 2023, with the Federal Register notice on 17 November 2023 (docket FDA-2021-D-0980, draft December 2021). It builds on V&V 40 and applies to physics-based and mechanistic models. — [Federal Register](https://www.govinfo.gov/content/pkg/FR-2023-11-17/html/2023-25470.htm); [FDA PDF](https://fda.gov/media/175618/download)

**Latest (2023–2026)**
- Ponnaluri SV, Hariharan P, Herbertson LH, Manning KB, Malinauskas RA, Craven BA (2023; online Nov 2022). "Results of the interlaboratory computational fluid dynamics study of the FDA benchmark blood pump." *Annals of Biomedical Engineering* 51(1):253–269. doi:10.1007/s10439-022-03105-w.
  - It reviews 24 submissions.
  - CFD of rotary pumps needs validation across all operating points and quantities of interest, and the diffuser is the hardest region.
  - Haemolysis can be predicted accurately despite inaccurate flow fields.
  - [doi](https://doi.org/10.1007/s10439-022-03105-w); [NSF PAR](https://par.nsf.gov/servlets/purl/10462950)
- Benemerito et al. (2025) FV versus LBM; see §2.

### Inferences
- **Recommended verification ladder for a new C++/GPU solver:**
  1. Poiseuille: exact, checks second-order wall convergence.
  2. Womersley at α = 4–16.
  3. Steady curved pipe (Dean flow).
  4. Purdue stenosis DNS (Varghese 2007).
  5. FDA nozzle at Re 500 (laminar), 2000 and 3500 (transitional; jet-breakdown location), and 6500, against multi-lab PIV with its reported uncertainty.
  6. A cerebral aneurysm challenge case with normalised WSS.
- Then frame each clinical use with V&V 40: the context of use and model risk set how much validation is enough.

- FDA nozzle specification, as given in FDA/NCIHub-hosted documents:
  - A 12 mm inner-diameter inlet tube leads into a 4 mm throat, 40 mm long. Measured values were 12.00 ± 0.13 mm and 4.00 ± 0.08 mm.
  - It is run in two orientations: through a gradual conical contraction, or through a sudden contraction into the throat.
  - Throat Re was 500, 2000, 3500, 5000 and 6500, over flow rates of 0.3–4.0 L/min.
  - The CFD specification is density 1056 kg/m³ and viscosity 3.5 cP.
  - The physical blood analogue (sodium iodide, glycerin and water) measured 6.9–8 cP and 1.65–1.73 g/cm³ at 21 °C, so match Re rather than properties.
  - [Herbertson et al. 2015, Artificial Organs (NCIHub PDF)](https://ncihub.cancer.gov/app/site/collections/2999/Herbertson_et_al-2015-Artificial_Organs.pdf); [FDA CPI CFD ASAIO poster 2009](https://ncihub.cancer.gov/app/site/collections/2911/CPI_CFD_ASAIO_Poster_May_2009.pdf)

### Gaps
- The diffuser cone angle of the FDA nozzle was not confirmed. The bibliographic details of Herbertson et al. (2015) were not DOI-verified.
- The current GitHub URL of the FDA "CFD and Blood-Damage Benchmarks" repository was not retrieved.

---

## 7. Clinical Applications: Coronary, Cerebral, Aortic, Venous, Fontan

### Takeaway
FFR-CT is the one image-based CFD application with randomised-trial-level and regulatory maturity: DISCOVER-FLOW, NXT, PLATFORM and PRECISE, all HeartFlow-funded. Aneurysm rupture-risk CFD is still a research tool: it discriminates rupture status in retrospective cohorts but cannot locate rupture sites, and low-versus-high WSS remains unresolved. Aortic dissection needs FSI for flap motion. Venous applications are emerging:
- iliac vein compression, where shear rate is elevated;
- portal hypertension, where virtual HVPG from CTA and CFD gives AUC 0.89;
- Fontan surgical planning, using CFD of exercise and respiration.

### Cited Findings

**Coronary FFR-CT**
- Taylor CA, Fonte TA, Min JK (2013). "Computational fluid dynamics applied to cardiac computed tomography for noninvasive quantification of fractional flow reserve: scientific basis." *JACC* 61(22):2233–2241. doi:10.1016/j.jacc.2012.11.083. Rest and hyperaemic flow and pressure are computed from CTA using allometric scaling and microvascular resistance models. — [doi](https://doi.org/10.1016/j.jacc.2012.11.083)
- Koo B-K, Erglis A, Doh J-H, Daniels DV, Jegere S, Kim H-S, et al. (2011). "Diagnosis of ischemia-causing coronary stenoses by noninvasive fractional flow reserve computed from coronary computed tomographic angiograms: results from the prospective multicenter DISCOVER-FLOW study." *JACC* 58(19):1989–1997. doi:10.1016/j.jacc.2011.06.066.
  - 103 patients and 159 vessels were studied.
  - Per vessel, FFR-CT had sensitivity 88% and specificity 82%, against 91% and 40% for CTA. These figures come from the TCTMD summary.
  - FFR-CT correlated with invasive FFR at r = 0.72.
  - HeartFlow-sponsored.
  - [doi](https://doi.org/10.1016/j.jacc.2011.06.066); [TCTMD](https://www.tctmd.com/news/discover-flow-published-noninvasive-ffr-shows-high-diagnostic-accuracy)
- Nørgaard BL, Leipsic J, Gaur S, Seneviratne S, Ko BS, Ito H, et al. (2014). "Diagnostic performance of noninvasive fractional flow reserve derived from coronary computed tomography angiography in suspected coronary artery disease: the NXT trial." *JACC* 63(12):1145–1155. doi:10.1016/j.jacc.2013.11.043.
  - 254 patients were studied.
  - The per-patient AUC was 0.90 against 0.81 for CTA (p = 0.0008), with markedly higher specificity.
  - Meeting reports gave different numbers, so do not mix them.
  - HeartFlow-funded.
  - [doi](https://doi.org/10.1016/j.jacc.2013.11.043); [PubMed](https://pubmed.ncbi.nlm.nih.gov/24486266/)
- Douglas PS, Pontone G, Hlatky MA, Patel MR, Nørgaard BL, Byrne RA, et al. (2015). "Clinical outcomes of fractional flow reserve by computed tomographic angiography-guided diagnostic strategies vs. usual care in patients with suspected coronary artery disease: the prospective longitudinal trial of FFRCT: outcome and resource impacts study" (PLATFORM). *European Heart Journal* 36(47):3359–3367. doi:10.1093/eurheartj/ehv444. In this non-randomised study of 584 patients, an FFR-CT-guided strategy reduced invasive angiograms showing no obstructive disease. — [doi](https://doi.org/10.1093/eurheartj/ehv444)
- **Latest:** Douglas PS, Nanna MG, Kelsey MD, Yow E, Mark DB, Patel MR, et al. (2023). "Comparison of an initial risk-based testing strategy vs usual testing in stable symptomatic patients with suspected coronary artery disease: the PRECISE randomized clinical trial." *JAMA Cardiology* 8(10):904–914. doi:10.1001/jamacardio.2023.2595.
  - 2103 patients were randomised to a strategy of deferred testing for minimal-risk patients and CTA with selective FFR-CT for the rest.
  - Catheterisation without obstructive CAD fell to 2.6% from 10.2%.
  - Death and MI were similar.
  - HeartFlow-funded.
  - [doi](https://doi.org/10.1001/jamacardio.2023.2595)
- Kim et al. (2010) supplies the coronary boundary conditions (§1).

**Intracranial aneurysm rupture risk**
- Xiang J, Natarajan SK, Tremmel M, Ma D, Mocco J, Hopkins LN, et al. (2011). "Hemodynamic-morphologic discriminants for intracranial aneurysm rupture." *Stroke* 42(1):144–152. doi:10.1161/STROKEAHA.110.592923. In 119 aneurysms (38 ruptured), hemodynamics was as important as morphology in discriminating rupture status. — [doi](https://doi.org/10.1161/STROKEAHA.110.592923)
- Cebral JR, Mut F, Weir J, Putman C (2011). "Quantitative characterization of the hemodynamic environment in ruptured and unruptured brain aneurysms." *AJNR* 32(1):145–151. doi:10.3174/ajnr.A2419.
  - In 210 aneurysms, ruptured ones had more concentrated inflow, higher maximum WSS, more concentrated shear and lower viscous dissipation ratio.
  - Low-WSS area did not differ.
  - The group differences held for steady and pulsatile inflow.
  - [doi](https://doi.org/10.3174/ajnr.A2419)
- Meng H, Tutino VM, Xiang J, Siddiqui A (2014). "High WSS or low WSS? Complex interactions of hemodynamics with intracranial aneurysm initiation, growth, and rupture: toward a unifying hypothesis." *AJNR* 35(7):1254–1262. doi:10.3174/ajnr.A3558. It proposes that low WSS and high OSI drive an inflammatory-cell pathway, while high WSS with a positive gradient drives a mural-cell pathway. — [doi](https://doi.org/10.3174/ajnr.A3558)
- Janiga et al. (2015) found rupture status was predictable but the rupture site was not (§6).

**Aortic dissection**: Alimohammadi (2015) and Bäumler (2020), §4.

**Venous flow**
- Assi IZ, Lynch SR, Ricker BD, Ranjane SV, Williams DM, Wakefield TW, Obi AT, Figueroa CA (2024). "A comparative study of altered hemodynamics in iliac vein compression syndrome." *Frontiers in Bioengineering and Biotechnology* 12:1302063. doi:10.3389/fbioe.2024.1302063.
  - Methods: CRIMSON, rigid walls, and ultrasound-Doppler inflow at the external and internal iliac veins over respiratory cycles. An RCR Windkessel at the IVC was tuned to about 10 mmHg. Carreau–Yasuda rheology was used on about 4 million tetrahedra with a time step of 10⁻⁴ s, and velocities were validated against ultrasound within 10%.
  - Shear rate in the left common iliac vein was 550 ± 103 s⁻¹ in patients against 75 ± 37 s⁻¹ in controls.
  - The left-to-right shear-rate ratio was 6.56 against 1.43 (p = 0.00008).
  - Peak shear rate was about 1653 against 256 s⁻¹, compared with a platelet-activation threshold of about 1000 s⁻¹.
  - With 4 patients and 4 controls, the authors describe the work as hypothesis-generating.
  - [doi](https://doi.org/10.3389/fbioe.2024.1302063)
- A Taiwan group built 31 patient-specific CTA-based CFD models of iliac vein compression (Bioengineering 2023, MDPI 10(6):688). It proposes the tilt angle of the left common iliac vein as an auxiliary severity indicator. Authors and DOI are **[unverified]**. — [PMC10295234](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10295234/)
- Lynch, Nama and Figueroa (2022): in an IVC–iliac model, non-Newtonian viscosity raised WSS and lowered transport (§3).
- Qi X, An W, Liu F, Qi R, Wang L, Liu Y, et al. (2019). "Virtual hepatic venous pressure gradient with CT angiography (CHESS 1601): a prospective multicenter study for the noninvasive diagnosis of portal hypertension." *Radiology* 290(2):370–377. doi:10.1148/radiol.2018180425.
  - 102 cirrhotic patients at 3 centres were studied (NCT02842697).
  - Virtual HVPG from 3D reconstruction and CFD correlated with invasive HVPG at R = 0.61.
  - The rule-out threshold was 7.3 mmHg and the rule-in threshold 13.0 mmHg for clinically significant portal hypertension (HVPG of 10 mmHg or more).
  - The validation AUC was 0.89 (95% CI 0.81–0.96, n = 73).
  - [doi](https://doi.org/10.1148/radiol.2018180425)

**Fontan circulation**
- Marsden AL, Vignon-Clementel IE, Chan FP, Feinstein JA, Taylor CA (2007). "Effects of exercise and respiration on hemodynamic efficiency in CFD simulations of the total cavopulmonary connection." *Annals of Biomedical Engineering* 35(2):250–263. doi:10.1007/s10439-006-9224-3.
  - Two patient-specific models were run with a respiration model.
  - Efficiency was high at rest, fell with exercise, and energy loss became significant even at moderate exercise.
  - [doi](https://doi.org/10.1007/s10439-006-9224-3)

### Inferences
- The clinical evidence ladder separates applications whose output is a pressure or flow ratio (FFR-CT, virtual HVPG) from WSS-based risk indices (aneurysm, atherosclerosis).
  - Pressure and flow quantities are robust across solvers (Steinman 2013) and have reached prospective and randomised trials.
  - WSS-based risk indices remain retrospective and variable (Valen-Sendstad 2018).
- Venous models need respiration-driven inflow, low-shear rheology, and often collapsible-wall considerations. These are good tests of an infrastructure's 0D coupling and rheology plug-ins.

### Gaps
- No verified 2023–2026 randomised evidence was found for aneurysm CFD or venous CFD.
- No recent Fontan virtual-surgery outcome studies were verified.
- NICE or other HTA assessments of FFR-CT were not checked.

---

## 8. Latest Work (2023–2026): GPUs, ML Surrogates, Real-Time Models, Open Solvers and Datasets

### Takeaway
Since 2023:
- GPU ports of production vascular solvers have matured (HemeLB-GPU; OpenLB GPU, with HIP from 2025; waLBerla GPU code generation).
- ML surrogates moved from point-cloud networks to SE(3)-equivariant mesh networks, GNN reduced-order models and neural-field operators, trained on CFD of tens to tens of thousands of geometries.
- Large synthetic aneurysm CFD datasets appeared: AneuG-Flow at NeurIPS 2025 (14,000 steady cases) and AneumoBench (about 85k simulations, preprint).
- Automated 0D/1D reduced-order models from SimVascular give 1–10% errors at a fraction of 3D cost.

Licensing splits the solver landscape. BSD/MIT/Apache options (svMultiPhysics, SimVascular, svZeroDSolver, MFEM, Nek5000/nekRS, deal.II, vmtk, MONAI, neuraloperator) can sit behind MIT interfaces. LGPL (HemeLB, DOLFINx, Oasis) requires care. GPL/AGPL (OpenFOAM, CRIMSON, OpenLB, waLBerla, Palabos, lifex, TetGen, Gmsh, VaMPy) cannot be incorporated.

### Cited Findings

**GPU solvers**: HemeLB-GPU (2023), OpenLB GPU and HIP, and waLBerla CUDA/HIP; see §2.

**ML surrogates and physics-informed learning (precursors, before 2023)**
- Raissi M, Yazdani A, Karniadakis GE (2020). "Hidden fluid mechanics: learning velocity and pressure fields from flow visualizations." *Science* 367(6481):1026–1030. doi:10.1126/science.aaw4741. It infers velocity and pressure from concentration fields with PINNs and includes an aneurysm example. — [doi](https://doi.org/10.1126/science.aaw4741)
- Kissas G, Yang Y, Hwuang E, Witschey WR, Detre JA, Perdikaris P (2020). "Machine learning in cardiovascular flows modeling: predicting arterial blood pressure from non-invasive 4D flow MRI data using physics-informed neural networks." *CMAME* 358:112623. doi:10.1016/j.cma.2019.112623. It uses 1D-model-constrained PINNs, and the outputs also calibrate Windkessel parameters. — [doi](https://doi.org/10.1016/j.cma.2019.112623)
- Pfaff T, Fortunato M, Sanchez-Gonzalez A, Battaglia PW (2021). "Learning mesh-based simulation with graph networks" (MeshGraphNets). ICLR 2021, arXiv:2010.03409. It is an encode–process–decode GNN on mesh graphs, one to two orders of magnitude faster than the source simulators. — [arXiv](https://arxiv.org/abs/2010.03409)
- Li G, Wang H, Zhang M, Tupin S, Qiao A, Liu Y, Ohta M, Anzai H (2021). "Prediction of 3D cardiovascular hemodynamics before and after coronary artery bypass surgery via deep learning." *Communications Biology* 4:99. doi:10.1038/s42003-020-01638-1. It reports about 600× speed-up over CFD, under 1 s per new geometry, and about 90% accuracy. — [doi](https://doi.org/10.1038/s42003-020-01638-1)

**ML surrogates (2023–2026)**
- Pegolotti L, Pfaller MR, Rubio NL, Ding K, Brugarolas Brufau R, Darve E, Marsden AL (2024). "Learning reduced-order models for cardiovascular simulations with graph neural networks." *Computers in Biology and Medicine* 168:107676. doi:10.1016/j.compbiomed.2023.107676. A GNN trained on 3D CFD predicts centreline pressure and flow, with errors below 2% (pressure) and 3% (flow) given enough training data, per the preprint. It outperforms physics-based 1D models. — [doi](https://doi.org/10.1016/j.compbiomed.2023.107676)
- Suk J, de Haan P, Lippe P, Brune C, Wolterink JM (2024). "Mesh neural networks for SE(3)-equivariant hemodynamics estimation on the artery wall." *Computers in Biology and Medicine* 173:108328. doi:10.1016/j.compbiomed.2024.108328.
  - It is an end-to-end SE(3)-equivariant graph convolution on triangular wall meshes, applied to WSS in synthetic coronaries.
  - Reported error figures differ between listings (a 7.6% approximation error versus an NMAE of 0.4%); check the PDF.
  - Code: github.com/sukjulian/coronary-mesh-convolution.
  - [doi](https://doi.org/10.1016/j.compbiomed.2024.108328)
- Suk J, Nannini G, Rygiel P, Brune C, Pontone G, Redaelli A, Wolterink JM (2025). "Deep vectorised operators for pulsatile hemodynamics estimation in coronary arteries from a steady-state prior." *Computer Methods and Programs in Biomedicine* 271:108958. doi:10.1016/j.cmpb.2025.108958. It is a discretisation-independent neural-field operator (message passing plus self-attention) mapping steady to pulsatile haemodynamics in 74 CCTA stenotic coronaries with pulsatile CFD ground truth. — [doi](https://doi.org/10.1016/j.cmpb.2025.108958)
- Preprints, **[not peer-reviewed as verified]**:
  - Garnier P, et al., including Goetz A and Hachem E (Dec 2025). "Graph deep learning for intracranial aneurysm blood flow simulation and risk assessment." arXiv:2512.09013. — [arXiv](https://arxiv.org/abs/2512.09013)
  - Rygiel P, Suk J, Brune C, Yeung KK, Wolterink JM (Jul 2025). "Wall shear stress estimation in abdominal aortic aneurysms: towards generalisable neural surrogate models." arXiv:2507.22817. — [arXiv](https://arxiv.org/abs/2507.22817)
  - Flouris K, et al. (Jul 2025). "Localized FNO for spatiotemporal hemodynamic upsampling in aneurysm MRI." arXiv:2507.13789. It is a Fourier-neural-operator variant (Laplacian eigenvectors) predicting WSS from 4D flow MRI. — [arXiv](https://arxiv.org/abs/2507.13789)

**Real-time and reduced-order modelling**
- Pfaller MR, Pham J, Verma A, Pegolotti L, Wilson NM, Parker DW, Yang W, Marsden AL (2022). "Automated generation of 0D and 1D reduced-order models of patient-specific blood flow." *International Journal for Numerical Methods in Biomedical Engineering* 38(10):e3639. doi:10.1002/cnm.3639.
  - The SimVascular pipeline builds 0D and 1D models automatically from 3D geometry.
  - Across 72 public models, flow and pressure errors were typically 1–10%.
  - 0D needed about a third of the 1D runtime.
  - svZeroDSolver is **BSD-3-Clause**.
  - [doi](https://doi.org/10.1002/cnm.3639)

**Datasets**
- Wilson NM, Ortiz AK, Johnson AB (2013). "The Vascular Model Repository: a public resource of medical imaging data and blood flow simulation results." *Journal of Medical Devices* 7(4):040923. doi:10.1115/1.4025983.
  - It covers the aorta and cerebral, coronary, aortofemoral and pulmonary arteries plus single-ventricle anatomy, distributed as DICOM and VTK.
  - Later counts vary: 144 projects and 111 results in one presentation; "over 300" projects in a Stanford posting; 206 subjects in a 2023 paper.
  - An NIH R01 (LM013120) is expanding it toward an "AI-ready" VMR.
  - [doi](https://doi.org/10.1115/1.4025983); [NIH poster](https://datascience.nih.gov/sites/default/files/3-Pegolotti%2C%20Luca%2C%20An%20AI-Ready%20Vascular%20Model%20Repository%20for%20Modeling.pdf)
- AneuX morphology database:
  - Zenodo doi:10.5281/zenodo.6678442 (June 2022) holds 750 intracranial aneurysm geometries from 3D rotational angiography, merging AneuX, @neurIST and Aneurisk, under a **CC BY-NC 4.0** licence.
  - The associated paper is Juchler N, Schilling S, Bijlenga P, Kurtcuoglu V, Hirsch S (2022). "Shape trumps size: image-based morphological analysis reveals that the 3D shape discriminates intracranial aneurysm disease status better than aneurysm size." *Frontiers in Neurology* 13:809391. doi:10.3389/fneur.2022.809391.
  - It has geometry only, not CFD.
  - [doi](https://doi.org/10.3389/fneur.2022.809391); [ZHAW record](https://digitalcollection.zhaw.ch/handle/11475/25379)
- AnXplore (2024): 101 rigid and FSI aneurysm simulations, public (§4).
- AneuG-Flow (NeurIPS 2025 Datasets & Benchmarks): "A large-scale synthetic dataset of diverse intracranial aneurysm geometries and hemodynamics".
  - It has 14,000 steady laminar CFD cases (more than 3M cells each) and pulsatile cases, given as 200 in the abstract and 730 on Hugging Face and in the paper.
  - Shapes come from a generative model trained on about 109–116 real MCA-bifurcation aneurysms.
  - The data include velocity, pressure and WSS fields; the PointNet and graph U-Net baseline reaches 4.67% relative L2 on normalised WSS.
  - Hosted at huggingface.co/datasets/whding123/AneuG-Flow. The licence is **[not stated in the indexed page]**.
  - [NeurIPS](https://proceedings.neurips.cc/paper_files/paper/2025/hash/e2b8ff0035bc9f572a7deefbcbea85bc-Abstract-Datasets_and_Benchmarks_Track.html); [HF](https://huggingface.co/datasets/whding123/AneuG-Flow)
- Aneumo, retitled in v2 as "AneumoBench: a source-linked benchmark for synthetic-geometry transfer in aneurysm CFD": Li X, et al. (2025 preprint, arXiv:2505.14717).
  - v1 abstract: 427 real geometries, 10,660 synthetic shapes and 8 steady flow rates, giving 85,280 CFD results, plus segmentation masks.
  - Counts differ between versions.
  - Code and data: github.com/Xigui-Li/Aneumo.
  - [arXiv](https://arxiv.org/abs/2505.14717)

**Open-source solvers and their licences** (checked 2026-11-04)

| Software | Kind | Licence | MIT-compatibility note |
|---|---|---|---|
| svMultiPhysics | C++ FE: fluid, solid, FSI, EP | **BSD-3-Clause** | permissive |
| svFSI (Fortran, predecessor) | FE multiphysics | MIT-style permission text with BSD disclaimer (GitHub "Other") | permissive |
| SimVascular | pipeline: segmentation to simulation | "provided under the BSD license"; the text is MIT-style, with "additional license information" in the source | permissive; check bundled third-party code |
| svZeroDSolver | 0D/1D | **BSD-3-Clause** | permissive |
| CRIMSON GUI | pipeline | **AGPL-3.0-or-later + Commons Clause** (no right to "Sell") | not OSI open source; incompatible |
| CRIMSON Flowsolver | FE (SimVascular-derived) | **GPL-3.0-or-later + Commons Clause** | not OSI open source; incompatible |
| FEniCS / DOLFINx | FE framework | **LGPL-3.0-or-later** | dynamic linking only; changes to it stay LGPL |
| deal.II | FE library | **Apache-2.0 WITH LLVM-exception, or LGPL-2.1-or-later** (dual) | permissive under the Apache option |
| OpenFOAM | FV | **GPL-3.0** | incompatible for incorporation |
| MFEM | FE library (GPU) | **BSD-3-Clause** | permissive |
| Nek5000 / nekRS | spectral element (nekRS GPU) | BSD-3-style (UChicago Argonne) | permissive |
| lifex | FE cardiac (deal.II-based) | **GPL-3.0-or-later** | incompatible |
| Oasis | FEniCS NS solver | **LGPL-3.0** | as LGPL |
| OasisMove | moving-domain Oasis | **GPL-3.0** | incompatible |
| VaMPy | vascular pipeline (Oasis) | **GPL-3.0** (some files LGPL/MIT) | incompatible |
| Kratos | multiphysics | core BSD-style (CIMNE); per-application licences | permissive core |
| vmtk | vessel modelling | **BSD** | permissive; optional TetGen |
| TetGen 1.6 | tetrahedral mesher | **AGPL-3.0-or-later** or commercial (WIAS) | incompatible unless licensed commercially |
| Gmsh | mesher | **GPL-2.0-or-later** (GitLab API) | incompatible for incorporation |
| MONAI | medical-image DL | **Apache-2.0** | permissive |
| neuraloperator | FNO library | **MIT** | permissive |

Sources:
- [svMultiPhysics](https://github.com/SimVascular/svMultiPhysics); [SimVascular LICENSE.md](https://github.com/SimVascular/SimVascular)
- [CRIMSON GUI LICENSE](https://github.com/carthurs/CRIMSONGUI); [CRIMSON Flowsolver LICENSE](https://github.com/carthurs/CRIMSONFlowsolver)
- [DOLFINx](https://github.com/FEniCS/dolfinx); [deal.II LICENSE.md](https://github.com/dealii/dealii/blob/master/LICENSE.md)
- [OpenFOAM licencing](https://www.openfoam.com/documentation/licencing)
- [TetGen (WIAS)](https://www.wias-berlin.de/software/index.jsp?id=TetGen)
- [vmtk LICENSE](https://github.com/vmtk/vmtk)
- [VaMPy LICENSE](https://github.com/KVSlab/VaMPy)
- The other entries come from the GitHub or GitLab licence API.

**Software papers**
- Updegrove A, Wilson NM, Merkow J, Lan H, Marsden AL, Shadden SC (2017; online Dec 2016). "SimVascular: an open source pipeline for cardiovascular simulation." *Annals of Biomedical Engineering* 45(3):525–541. doi:10.1007/s10439-016-1762-8. — [doi](https://doi.org/10.1007/s10439-016-1762-8)
- Arthurs CJ, Khlebnikov R, Melville A, Marčan M, Gomez A, Dillon-Murphy D, et al. (34 authors) (2021). "CRIMSON: an open-source software framework for cardiovascular integrated modelling and simulation." *PLoS Computational Biology* 17(5):e1008881. doi:10.1371/journal.pcbi.1008881. The paper says "open-source", but the repository licence adds the Commons Clause. — [doi](https://doi.org/10.1371/journal.pcbi.1008881)

### Inferences
- **Implementation choice for an MIT-licensed C++/GPU infrastructure on CT/MRI geometries.** Two tracks, cross-validated against each other and against the FDA nozzle:
  - **Track A (voxel-native, GPU): an own LBM.** D3Q19 TRT/MRT, BFL walls from a signed-distance field, Windkessel outlets via density, and a pluggable μ(γ̇) law, following HemeLB's documented design (LGPL code can be consulted; reimplement rather than copy GPL codes).
  - **Track B (mesh-based, reference): an FE solver behind the infrastructure's interfaces.** Use svMultiPhysics (BSD-3) or MFEM/deal.II, with RCR/coronary 0D coupling, backflow stabilisation, and CMM/ALE FSI.
  - Meshing should avoid TetGen and Gmsh unless they run as separate processes, or use vmtk (BSD) with a permissive mesher.
- **Surrogates.** Train on the VMR (open), AnXplore (public) and AneuG-Flow (check the licence). AneuX is non-commercial (CC BY-NC), which constrains redistribution of derived models. neuraloperator (MIT) and MONAI (Apache-2.0) are licence-compatible ML stacks.

### Gaps
- No verified peer-reviewed report of interactive or real-time 3D patient-specific CFD (as opposed to 0D/1D or ML surrogates) in 2023–2026 was found.
- The licence of the AneuG-Flow data was not confirmed.
- The VMR's current licence and terms were not checked.
- The Gmsh licence was read from the GitLab API; Gmsh's GPL linking exceptions were not checked.
