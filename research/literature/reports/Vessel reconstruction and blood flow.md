# Geometry, Not Solvers, Limits Vascular Simulation

The literature agrees on one finding above all: **the reconstructed geometry limits the accuracy of image-based haemodynamics more than the flow solver does**. When 26 teams segmented the same five aneurysms, their volumes differed by up to 20% and their ostium areas by up to 30%. When 26 teams ran whole pipelines, the interquartile range of sac-averaged wall shear stress reached 56%. Yet 25 groups computing a pressure drop agreed within 8%. That asymmetry explains the clinical record. Pressure-type outputs have reached a randomised trial (FFR-CT) and a prospective venous study (a virtual hepatic venous pressure gradient, AUC 0.89); shear-type risk indices remain retrospective. For veins the evidence is lopsided. Imaging methods and labelled datasets are overwhelmingly arterial: the main open-licensed source of large-vein labels (venae cavae, iliac, portal–splenic, brachiocephalic, pulmonary and jugular veins) is TotalSegmentator. Venous mechanics, by contrast, is mature in one dimension: a collapsible tube law, valves, gravity and the calf muscle pump, exact Riemann solutions through collapse, and a closed-loop network of 2,185 arteries and 189 veins. But no public venous benchmark with reference data exists. The infrastructure should therefore build in this order:

1. an explicit coordinate chain from DICOM's LPS frame to SI metres;
2. classical, testable enhancement, segmentation, distance transforms, topology-preserving skeletons and surfaces, scored by metrics that see topology and calibre;
3. a vascular graph that carries cross-sectional area as well as inscribed radius;
4. a reduced-order solver whose single tube-law family covers arteries, veins and dural sinuses, verified against exact solutions and the Boileau and PWDB benchmarks;
5. a voxel-native lattice Boltzmann solver, verified against Poiseuille, Womersley and the FDA nozzle;
6. transport with stiff chemistry, held to closed forms and public-domain reaction networks;
7. only then, a forward measurement operator for 4D flow MRI, with assimilation built on it.

Licences constrain the plan as much as the science does. TetGen (AGPL), Gmsh and CGAL's meshing (GPL), nearly every open lattice Boltzmann code (GPL or AGPL), CRIMSON (AGPL with the Commons Clause) and vesselFM (GPL code, non-commercial weights) cannot enter an MIT core. Those methods must be reimplemented. ITK, VTK, VMTK built without TetGen, svMultiPhysics, svZeroDSolver, openBF and nnU-Net can serve as references.

## 1. Reconstruction error outweighs solver error in every multi-team comparison

### 1.1 Multi-team comparisons isolate geometry as the dominant error

The community challenges are the cleanest evidence, because they hold the patient fixed and let the methods vary.

- **Segmentation.** In MATCH Phase I, 26 groups segmented the same 3D rotational angiograms of five aneurysms with no common software. **Aneurysm volumes differed by up to 20% and ostium areas by up to 30%**, surfaces were systematically oversegmented by about 10% against a high-resolution 2D reference, and morphological parameters varied by up to 25% ([Berg et al. 2018](https://doi.org/10.1007/s13239-018-00376-0)).
- **Whole pipelines.** In the 2015 International Aneurysm CFD Challenge, 28 datasets came from 26 teams, each free to choose segmentation, extent, inflow and solver. Variability in segmentation, model extent and inflow rate drove **interquartile ranges of sac-averaged WSS of up to 56%**. Normalising by the parent artery's WSS reduced the spread substantially, and team experience did not predict it ([Valen-Sendstad et al. 2018](https://doi.org/10.1007/s13239-018-00374-2)).
- **Repeat scans.** Repeated MRI of the same carotid bifurcations reproduced the lumen to about 0.15 mm (5%). Yet time-averaged WSS was reproducible only to about 37%, a figure "determined largely by the precision of the lumen boundary extraction" ([Thomas et al. 2003](https://doi.org/10.1114/1.1540102)).
- **Phantoms.** An earlier phantom study already found about 15% mean error between WSS fields computed on MR-derived and true geometries ([Moore et al. 1999](https://doi.org/10.1114/1.163)).
- **Outflow strategy.** Outflow strategy alone moved the computed outlet flows of MCA aneurysms by **up to 70%**. Only a "splitting" rule was insensitive to arbitrary truncation of the model ([Chnafa et al. 2018](https://doi.org/10.3174/ajnr.A5484)).

Pressure, by contrast, is robust. In the ASME 2012 challenge, 25 groups mostly agreed on the **peak-systolic pressure drop within 8%** of each other, while flow patterns in the sac varied ([Steinman et al. 2013](https://doi.org/10.1115/1.4023382)).

Solver numerics are a real but second-order effect, and they bite where shear matters. Across 12 MCA aneurysms, normal-resolution solution strategies correlated with high-resolution ones at R² > 0.9 for most reduced indices but only **0.23 for OSI** ([Valen-Sendstad & Steinman 2014](https://doi.org/10.3174/ajnr.A3793)). With a minimally dissipative second-order solver, differences became negligible at about **800,000 tetrahedra and 5,000 time steps per cycle** ([Khan et al. 2015](https://doi.org/10.3174/ajnr.A4263)).

### 1.2 Pressure-type outputs have earned clinical evidence; shear-type indices have not

The clinical evidence follows the same divide.

- **FFR-CT (coronary).** The only image-based CFD application with randomised-trial evidence; all the studies below were funded by HeartFlow.
  - In DISCOVER-FLOW (103 patients, 159 vessels), per-vessel sensitivity and specificity were 88% and 82%, against 91% and 40% for CTA alone ([Koo et al. 2011](https://doi.org/10.1016/j.jacc.2011.06.066); figures from the [TCTMD summary](https://www.tctmd.com/news/discover-flow-published-noninvasive-ffr-shows-high-diagnostic-accuracy)).
  - In NXT, the per-patient AUC was **0.90 against 0.81** for CTA alone ([Nørgaard et al. 2014](https://doi.org/10.1016/j.jacc.2013.11.043)).
  - In the randomised PRECISE trial (2,103 patients), a strategy with selective FFR-CT cut catheterisations without obstructive disease from **10.2% to 2.6%** ([Douglas et al. 2023](https://doi.org/10.1001/jamacardio.2023.2595)).
- **Venous pressure.** In a prospective multicentre study of 102 cirrhotic patients, a virtual hepatic venous pressure gradient computed by CFD on CTA correlated with the invasive gradient at R = 0.61. Its validation AUC was **0.89** for clinically significant portal hypertension ([Qi et al. 2019](https://doi.org/10.1148/radiol.2018180425)).
- **Aneurysm rupture risk.** Remains a research tool. In the 2013 Rupture Challenge, 81% of CFD groups identified which aneurysm had ruptured, against 91% of neurosurgeons, and none located the rupture site ([Janiga et al. 2015](https://doi.org/10.3174/ajnr.A4157)). Whether low or high WSS drives rupture remains unresolved ([Meng et al. 2014](https://doi.org/10.3174/ajnr.A3558)).

Diagnosis is outside the infrastructure's scope, but the lesson for its design is direct. **Geometry and its uncertainty must be first-class objects of the infrastructure, not preprocessing.** The uncertainty in geometry is the largest and the least often propagated ([Sankaran et al. 2016](https://doi.org/10.1016/j.jbiomech.2016.01.002); [Maher et al. 2021](https://doi.org/10.1016/j.cma.2021.114038)).

## 2. Hessian filters, topology-aware losses and one vessel foundation model define segmentation in 2026

### 2.1 Classical enhancement is cheap, exactly specified and still benchmarked

The Hessian-eigenvalue filters remain the backbone of classical vessel enhancement. Each costs one 3×3 symmetric eigenproblem per voxel per scale.

**Frangi.** Sort the eigenvalues of the Gaussian-derivative Hessian by magnitude, |λ1| ≤ |λ2| ≤ |λ3|. The filter combines three ratios:

- R_A = |λ2|/|λ3| separates plates from lines;
- R_B = |λ1|/√|λ2λ3| separates blobs from lines;
- S = √(Σλi²) measures "structureness".

They combine as V = (1 − e^(−R_A²/2α²))·e^(−R_B²/2β²)·(1 − e^(−S²/2c²)), with α = β = 0.5 and c by default half the maximum of S ([Frangi et al. 1998](https://doi.org/10.1007/BFb0056195); formula as implemented in [scikit-image v0.24](https://raw.githubusercontent.com/scikit-image/scikit-image/v0.24.0/skimage/filters/ridges.py)).

**Implementations disagree** on scale normalisation and sign. scikit-image applies σ² scaling in `sato` but not in `frangi`. A from-scratch implementation should therefore be checked against ITK (Apache-2.0) or the Lamy benchmark's C++ code, not against scikit-image alone.

**Jerman.** Jerman's filter targets Frangi's non-uniform response across the lumen and at bifurcations. It regularises the third eigenvalue to λρ, then responds with V = λ2²(λρ − λ2)·27/(λ2 + λρ)³, saturating to 1 where λ2 ≥ λρ/2 > 0 ([Jerman et al. 2016](https://doi.org/10.1109/TMI.2016.2550102); [author's code](https://raw.githubusercontent.com/timjerman/JermanEnhancementFilter/master/vesselness3D.m)). The author's MATLAB code carries no licence, so it should be reimplemented from the paper.

**Optimally oriented flux (OOF).** OOF replaces the Hessian by the flux of the image gradient through a sphere. By the divergence theorem, each of the six distinct entries of the flux matrix is a convolution with a closed-form Fourier kernel, Ψ(u) = 4πr·u_i u_j·e^(−2(π‖u‖σ)²)·‖u‖⁻²·(cos 2πr‖u‖ − sin(2πr‖u‖)/(2πr‖u‖)). Radii are sampled every half voxel and the eigenvalues normalised by 4πr². Because the response depends only on where the sphere's surface touches a boundary, OOF is **robust to closely adjacent structures**, where Hessian filters blur ([Law & Chung 2008](https://home.cse.ust.hk/~achung/eccv08_law_chung.pdf)). RORPO offers a morphological, non-Hessian alternative ([Merveille et al. 2018](https://doi.org/10.1109/TPAMI.2017.2672972)).

**The Lamy benchmark.** It compares seven filters in one C++ framework on IRCAD liver CT, the Bullitt brain MRA set and synthetic VascuSynth data, region by region and by diameter class. Its two published conclusions are that **bifurcations "are often missed by vesselness filters"** and that default parameters are "often not enough" ([Lamy et al. 2022](https://doi.org/10.1109/TMI.2022.3192679)). The per-filter ranking could not be read.

All of this maps onto foundations the infrastructure already has. Frangi, Sato and Jerman need separable Gaussian derivatives and a 3×3 eigensolver. OOF needs, for each radius, six FFT-domain multiplications by Ψ, then the same eigensolver; both are in calx.

### 2.2 Classical segmentation supplies baselines and centrelines without training data

The classical toolbox remains the training-free baseline that deep networks are measured against. Its canonical review is [Lesage et al. 2009](https://doi.org/10.1016/j.media.2009.07.011).

- **Mixture models (TOF-MRA).** For time-of-flight angiograms, models of the intensity histogram fitted by EM give a threshold or a posterior map cheaply ([Wilson & Noble 1999](https://doi.org/10.1109/42.811277); [Hassouna et al. 2006](https://doi.org/10.1016/j.media.2004.11.009)).
- **Velocity coherence (phase-contrast MRA).** For phase-contrast angiograms and therefore 4D flow, segmentation can also use the coherence of flow directions among neighbouring voxels, which no TOF method can ([Chung et al. 2004](https://doi.org/10.1109/TMI.2004.836877)).
- **Minimal paths.** These solve the eikonal equation |∇U| = P by fast marching from a seed, with the potential P taken, for example, as 1/(ε + vesselness). Tracing back along the gradient of U gives a globally optimal centreline between two points ([Cohen & Kimmel 1997](https://doi.org/10.1023/A:1007922224810)). Lifting the path to (x, y, z, r) optimises the radius jointly with the path ([Li & Yezzi 2007](https://doi.org/10.1109/TMI.2007.903696)).
- **Ridge traversal.** Steps along the tube direction and re-centres in the normal plane ([Aylward & Bullitt 2002](https://doi.org/10.1109/42.993126)).
- **Graph cuts.** Globally optimal for their energy ([Boykov & Jolly 2001](https://doi.org/10.1109/ICCV.2001.937505)), but their boundary term shrinks thin vessels unless flux terms are added. The widely used Boykov–Kolmogorov max-flow code reportedly has had research-only licensing, which the notes did not verify. An MIT core should therefore write its own push–relabel solver.

### 2.3 Deep networks win on accuracy, but licences decide what can be borrowed

**Baseline: nnU-Net.** nnU-Net (Apache-2.0) remains the reference baseline. Its 2024 "revisited" study argues that well-configured convolutional U-Nets are still not beaten by transformer or Mamba variants ([Isensee et al. 2021](https://doi.org/10.1038/s41592-020-01008-z); [Isensee et al. 2024](https://arxiv.org/abs/2404.09556)).

**Topology-aware losses.** Vessel-specific gains come mainly from losses that see topology. All four have permissive code that can sit behind an MIT interface:

| Loss | Code licence | Source |
|---|---|---|
| clDice | MIT | [Shit et al. 2021](https://doi.org/10.1109/CVPR46437.2021.01629) |
| Betti matching | MIT | [Stucki et al. 2023](https://arxiv.org/abs/2211.15272); [3D version, 2024](https://arxiv.org/abs/2407.04683) |
| Skeleton Recall | Apache-2.0 | [Kirchhoff et al. 2024](https://doi.org/10.1007/978-3-031-72980-5_13) |
| cbDice | Apache-2.0 | [Shi et al. 2024](https://arxiv.org/abs/2407.01517) |

**vesselFM, the one vessel foundation model.** It is the only foundation model aimed at vessels across modalities. It is trained on roughly 70% domain-randomised synthetic data, 20% real data from 17 annotated sources and 10% flow-matching samples. Its zero-shot Dice is compared below with the general promptable models ([Wittmann et al. 2025](https://arxiv.org/html/2411.17386v2)):

| Zero-shot Dice | 7T TOF-MRA (SMILE-UHURA) | Hepatic-vessel CT (MSD Task08) |
|---|---|---|
| vesselFM | **74.7%** | **29.7%** |
| VISTA3D | 5.1% | 23.8% |
| MedSAM-2 | 3.9% | 14.5% |
| SAM-Med3D | 2.1% | 7.9% |

Fine-tuning on three patches lifts vesselFM's MSD score only to 45.0%. Part of the shortfall comes from MSD's reference labels, which leave the aorta unannotated. Three further limits apply:

- **Licences.** The code is **GPL-3.0** and the weights are **Open RAIL++-M, research and non-commercial only**.
- **A non-comparable checkpoint.** An August 2026 checkpoint was trained on a set that includes the evaluation data, so results obtained with it are no longer comparable with the paper's ([vesselFM repository](https://github.com/bwittmann/vesselFM)).
- **What it does not do.** It segments binary masks only and leaves connectivity repair to future work.

vesselFM is therefore an external comparison for the infrastructure, never a borrowed component.

**Licences of the other general models.** These vary by task and must be read individually.

- **TotalSegmentator.** The code is Apache-2.0, and its weights for the `total`, `lung_vessels`, `liver_vessels` and `headneck_bones_vessels` tasks are Apache-2.0. Its `coronary_arteries`, `heartchambers_highres` and `brain_structures` tasks are free only for non-commercial use, and `brain_aneurysm` is CC BY-NC 4.0 ([TotalSegmentator README](https://github.com/wasserth/TotalSegmentator)).
- **nnInteractive.** The weights are CC BY-NC-SA 4.0 ([nnInteractive](https://github.com/MIC-DKFZ/nnInteractive)).
- **VISTA3D.** Its README is internally inconsistent about the licences ([VISTA3D](https://github.com/Project-MONAI/VISTA/tree/main/vista3d)).

### 2.4 Veins are the least served target of vessel segmentation

Almost every method above was developed and benchmarked on arteries: coronary CTA, the Circle of Willis, TOF-MRA and the aorta. Venous labels exist mostly as classes inside general anatomical models.

| Source | Modality | Venous structures labelled | Licence | Role for the infrastructure |
|---|---|---|---|---|
| TotalSegmentator v2, `total` task | CT, 1,228 cases | superior and inferior vena cava, portal and splenic vein, iliac veins, brachiocephalic veins, pulmonary vein | data CC BY 4.0; weights Apache-2.0 | first reference and baseline ([data](https://doi.org/10.5281/zenodo.10047292); [README](https://github.com/wasserth/TotalSegmentator)) |
| TotalSegmentator `lung_vessels`, `headneck_bones_vessels` | CT | lung veins; internal jugular veins | weights Apache-2.0 | baseline |
| TotalSegmentator `brain_structures` | CT | dural venous sinuses | non-commercial | comparison only |
| MSD Task08 HepaticVessel | portal-venous CT | hepatic vessels and tumour (aorta not labelled) | CC BY-SA 4.0 per AWS, which links a CC BY deed | benchmark ([AWS registry](https://registry.opendata.aws/msd/)) |
| HVA-CT (2026) | CT | hepatic venous annotations | CC BY-SA 4.0 (content unchecked) | candidate ([Zenodo](https://doi.org/10.5281/zenodo.19850108)) |
| 3D-IRCADb-01 | CT, 20 cases | liver vessels (which venous systems: unverified) | CC BY-NC-ND 4.0 | fetch-only ([IRCAD](https://www.ircad.fr/research/data-sets/liver-segmentation-3d-ircadb-01/)) |
| KiPA22 | CTA, 130 cases | renal arteries and veins | application form, no licence found | fetch-only ([KiPA22](https://kipa22.grand-challenge.org/)) |
| HiPaS | CT | pulmonary arteries and veins | dataset licence unchecked | candidate ([Chu et al. 2025](https://arxiv.org/abs/2404.07671)) |

**No dataset found offers vessel labels for phase-contrast MRA or 4D flow under an open licence.** The one strong venous geometric reference located is indirect: synthetic trees optimised under different principles gave branching exponents between 2.0 and 3.0, validated against a human liver portal-venous corrosion cast ([Jessen et al. 2024](https://doi.org/10.1109/tbme.2023.3334758)).

The inference for a veins-first programme is that it should start from contrast CT in the portal-venous phase. TotalSegmentator's vein classes serve there as reference and baseline. Segmentation from MR venography and venous 4D flow should be treated as research.

### 2.5 Metrics must see topology and calibre

Volumetric Dice is dominated by large vessels and blind to connectivity, so vessel work pairs it with three other kinds of metric:

- **Skeleton overlap.** clDice is the harmonic mean of topological precision |S_P ∩ V_L|/|S_P| and topological sensitivity |S_L ∩ V_P|/|S_L|, computed from the skeletons S of the prediction P and the label L.
- **Topology.** Betti-number errors, or the Betti matching error.
- **Surface distances.** HD95 and ASSD.

The challenges add task-specific protocols:

- **VESSEL12** scored ROC curves on points labelled by three annotators, reported separately for small, medium and large vessels ([VESSEL12](https://vessel12.grand-challenge.org/Details/)).
- **TopCoW** released 250 MRA and CTA scans with 13 Circle-of-Willis classes and graph-edge labels. Its best solutions reached nearly 90% class-average Dice ([TopCoW](https://arxiv.org/abs/2312.17670); [TopBrain proposal](https://doi.org/10.5281/zenodo.15084013)).

Human agreement depends strongly on calibre. Inter-operator Dice for the whole aorta on contrast CT is about **96%** ([arXiv:2002.03463](https://arxiv.org/abs/2002.03463)), while agreement on small vessels is far lower. **Dice must therefore be reported per diameter class** before a method is compared with human agreement.

Every one of these metrics is computable in C++ from masks, a skeletoniser and a homology routine. morphe already holds homology over ℤ/2 and persistent homology.

## 3. Surfaces, centrelines and graphs inherit every segmentation error and add their own

### 3.1 A watertight surface comes from an interior isosurface and non-shrinking smoothing

The method depends on what the segmentation provides.

- **Probability map.** Take the 0.5 isosurface of the probability field. Marching Cubes 33 guarantees topology for the trilinear interpolant ([Lewiner et al. 2003](https://doi.org/10.1080/10867651.2003.10487582)). MC33's own case table had errors that break that guarantee, later corrected ([Custodio et al. 2013](https://doi.org/10.1016/j.cag.2013.04.004)).
- **Hard label map.** Constrained surface nets place one vertex per surface cube and relax it while keeping it inside its cube, which avoids staircasing ([Gibson 1998](https://doi.org/10.1007/BFb0056277)). The 2024 parallel SurfaceNets is "typically one to two orders of magnitude faster" than sequential discrete isosurfacing and handles multiple labels ([Schroeder et al. 2024](https://arxiv.org/abs/2401.14906)).
- **Smoothing.** Taubin's λ|μ filter alternates a shrinking Laplacian step with an inflating one. It acts as a low-pass filter on the mesh Laplacian's eigenbasis and does not shrink the surface ([Taubin 1995](https://doi.org/10.1145/218380.218473)).

Non-shrinking smoothing matters for small vessels. Plain Laplacian smoothing removes radius, and a one-voxel loss is a large relative error in a vessel two to four voxels across. MC33's guarantees cover the interpolant only, not the segmentation. A spurious bridge between an artery and its companion vein, or a broken thin vein, must be prevented upstream by topology-aware losses or removed by explicit editing.

### 3.2 Centrelines carry radius, but an inscribed sphere misreads a vein

Three families of centreline coexist.

- **Voxel thinning.** Lee–Kashyap–Chu; Palágyi–Kuba; Bertrand–Couprie isthmus thinning. It preserves topology by construction and scales to huge volumes, but it leaves spurs that need pruning ([Lee et al. 1994](https://doi.org/10.1006/cgip.1994.1042); [Bertrand & Couprie 2015](https://doi.org/10.1016/j.gmod.2015.05.001)).
- **Voronoi centrelines (Antiga; VMTK).** Minimal paths of an eikonal equation on the Voronoi diagram of a surface, which carry the maximal-inscribed-sphere radius exactly ([Antiga et al. 2003](https://doi.org/10.1109/TMI.2003.812261); [VMTK tutorial](http://vmtk.github.io/tutorials/Centerlines.html)).
  - They need watertight, reasonably smooth surfaces and labelled endpoints.
  - VMTK's tutorial calls the inverse radius the "wave speed", while Antiga's formulation minimises ∫ds/R. An implementer must fix the convention from the source code.
  - A 2026 alternative solves the eikonal equation by finite elements on the mesh itself. It "considerably outperforms VMTK" on the rough surfaces that machine-learning segmentations produce, and needs no outlet labelling ([Sanchez Van Moer & Shadden 2026](https://doi.org/10.1002/cnm.70183)).
- **Image-space minimal paths and trackers.** These need no prior segmentation (§2.2). A CNN orientation tracker reached 93.7% overlap and 0.21 mm mean distance on the CAT08 coronary standard ([Wolterink et al. 2019](https://doi.org/10.1016/j.media.2018.10.005); [Schaap et al. 2009](https://doi.org/10.1016/j.media.2009.06.003)).

**Radius is reliable only when the lumen spans several voxels.**

- **Diameter.** Accurate diameters in TOF and PC angiography need the lumen to span **at least three pixels** ([Hoogeveen et al. 1998](https://doi.org/10.1002/jmri.1880080608)).
- **Flow volume.** Accurate flow volumes need **five to six voxels** across ([Dyverfeldt et al. 2015](https://doi.org/10.1186/s12968-015-0174-5)).
- **Quantisation (an inference).** A surface accurate to ±0.5 voxel puts the relative radius error near 0.5/R_vox: 25% at two voxels, 10% at five.

**Veins need more than an inscribed radius.** The maximal inscribed sphere measures the *minimum* half-width of a lumen, so for elliptical or partly collapsed cross-sections it differs systematically from the radius of a circle of equal area. Veins are often elliptical or partly collapsed at low transmural pressure. A one-dimensional venous model, moreover, is written in area A and a reference area A0 (§4.3). The vascular graph should therefore store, at each centreline station:

- the cross-sectional area and its contour, or a fitted ellipse;
- the inscribed radius beside them.

### 3.3 Vascular graphs have loops, and Murray's cube law fails in most territories

A vascular graph collapses skeleton junctions into nodes and branches into edges carrying length, radius, tortuosity and orientation. Two facts shape its design.

- **Scale.** Whole mouse-brain graphs exceed **10⁸ segments** ([Kirst et al. 2020](https://doi.org/10.1016/j.cell.2020.01.028)), and VesselGraph (MIT) is the standard learning benchmark on them ([Paetzold et al. 2021](https://arxiv.org/abs/2108.13233)).
- **Loops.** The cortical microvasculature forms interconnected loops, not a tree ([Blinder et al. 2013](https://doi.org/10.1038/nn.3426)). The graph's first Betti number, β1 = E − V + C, counts independent loops. Comparing the graph's β0 and β1 with those of the mask, computed by cubical homology, detects topological errors that skeletonisation introduced.

Murray's minimum-work law, r0³ = r1³ + r2³, is widely used to split flow among outlets and to repair radii. The measurements say it should be estimated per territory, not assumed:

| Territory | Exponent n in r0ⁿ = r1ⁿ + r2ⁿ | Source |
|---|---|---|
| Human cerebral arteries (internal carotid system) | 2.9 ± 0.7 (N = 157) | [Rossitti & Löfgren 1993](https://doi.org/10.1161/01.str.24.3.371) |
| Human epicardial coronaries | 2.15 for inlet flow; 2.38 for microvascular resistance | [Taylor et al. 2022](https://doi.org/10.3389/fphys.2022.871912) |
| Coronary meta-analysis, 1,070 trees | **2.39 (95% CI 2.24–2.54)** | [Taylor et al. 2024](https://doi.org/10.1152/ajpheart.00142.2024) |
| Human pulmonary arteries | 2.31 ± 0.60 (n = 7) | [Altieri Correa et al. 2025](https://doi.org/10.1115/1.4068886) |
| Synthetic trees against a portal-venous cast | 2.0–3.0, depending on the constraint | [Jessen et al. 2024](https://doi.org/10.1109/tbme.2023.3334758) |

Two consequences follow. A fixed n = 3 builds in bias. Exponent fits should use only branches several voxels across, because radius error grows as calibre falls.

### 3.4 Meshes need per-patient convergence studies, and most meshers are copyleft

**Pipelines.** The established open pipelines are:

- VMTK: level sets, Voronoi centrelines, meshing with boundary layers ([Antiga et al. 2008](https://doi.org/10.1007/s11517-008-0420-1));
- SimVascular ([Updegrove et al. 2017](https://doi.org/10.1007/s10439-016-1762-8));
- CRIMSON ([Arthurs et al. 2021](https://doi.org/10.1371/journal.pcbi.1008881)).

**Convergence.** Grid convergence must be checked per patient. Among five aneurysm models, four converged monotonically with 12–16% uncertainty on the finest grids. One complex dome oscillated, so "each patient-specific model required its own grid convergence study" ([Hodis et al. 2012](https://doi.org/10.1016/j.jbiomech.2012.07.030); procedure by [Roache 1994](https://doi.org/10.1115/1.2910291)).

**Hexahedral meshing.** Meshing directly from centrelines is a newer route. It meshed 60 cerebral networks with 92% of vessels and 83% of bifurcations needing no manual repair ([Decroocq et al. 2023](https://doi.org/10.1016/j.media.2023.102912)).

**Licences.** These shape what the infrastructure can link:

- **Copyleft, to run only as external processes:**
  - TetGen 1.5 and later is AGPLv3, and the TetGen 1.4.3 that VMTK bundles forbids commercial distribution ([WIAS](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1); [VMTK's TetGen licence](https://raw.githubusercontent.com/vmtk/vmtk/master/vtkVmtk/Utilities/tetgen1.4.3/LICENSE));
  - Gmsh is GPL-2.0-or-later ([Gmsh](https://gmsh.info/LICENSE.txt));
  - every CGAL package that a vessel pipeline would use (meshing, polygon-mesh processing, skeletonisation, isosurfacing) is GPL, only the kernel being LGPL ([CGAL packages](https://doc.cgal.org/latest/Manual/packages.html)).
- **Compatible with an MIT core:**
  - VTK and Geogram (BSD-3);
  - ITK (Apache-2.0);
  - fTetWild (MPL-2.0, file-level copyleft only) ([fTetWild](https://raw.githubusercontent.com/wildmeshing/fTetWild/master/LICENSE.MPL2)).

### 3.5 Imaging physics fixes which vessels can be reconstructed at all

**Modalities.** The resolving limit of each modality bounds the vessels it can reconstruct:

- **CT angiography** needs more than 250 HU of intravascular opacification and reconstructs slices thinner than 1 mm, 0.5–0.6 mm at best ([SCCT 2016](https://doi.org/10.1016/j.jcct.2016.10.002)). Photon-counting CT reaches 0.2 mm slices ([Mergen et al. 2022](https://doi.org/10.1097/rli.0000000000000897)).
- **4D flow MRI** typically runs at 1.5–3 mm isotropic and 30–40 ms ([Dyverfeldt et al. 2015](https://doi.org/10.1186/s12968-015-0174-5)).
- **7T TOF**, with prospective motion correction, reaches 0.15 mm isotropic ([Mattern et al. 2018](https://doi.org/10.1002/mrm.27033)).

Combined with the three-voxel and five-to-six-voxel criteria (§3.2), simple arithmetic gives the reliable ranges:

| Modality | Resolution | Reliable diameter from | Accurate flow from |
|---|---|---|---|
| Conventional CTA | 0.5 mm | about 1.5 mm | — |
| Photon-counting CT | 0.2 mm | about 0.6 mm | — |
| 4D flow MRI | 1.5–3 mm | — | about 8–18 mm: the aorta and great vessels (by size also the central veins; typical calibres were not verified in the notes) |

Radii read from 4D flow images of intracranial vessels are therefore unreliable.

**Coordinates.**

- **DICOM.** Image Position (Patient) gives the centre of the first voxel, and Image Orientation (Patient) the row and column direction cosines, in a right-handed LPS frame. That frame's x increases to the patient's left, y posteriorly and z toward the head ([DICOM PS3.3 §C.7.6.2](https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.7.6.2.html)).
- **NIfTI.** RAS, through either a quaternion qform with its sign factor qfac or a general affine sform. "To transform between DICOM and NIFTI, you just have to negate the x- and y-coordinates" ([nifti1.h](https://raw.githubusercontent.com/NIFTI-Imaging/nifti_clib/master/niftilib/nifti1.h)).

**The pixel-spacing trap** follows from these definitions. Stepping along a row, in direction X, uses the *second* Pixel Spacing value, and stepping down a column uses the first. The voxel-to-LPS affine is therefore A = [X·PS₂, Y·PS₁, N·Δk | IPP]. The slice direction N and step Δk should come from the Image Positions of successive slices, never from the nominal Slice Thickness or the implementation-specific Slice Location.

**What the notes did not cover.**

- Patient-position variants (head first or feet first).
- DICOM's enhanced multi-frame objects.

Both matter for veins. Venous calibre depends on posture, and images are acquired supine while venous physiology is often studied upright.

## 4. Veins need collapse, valves and gravity, which reduced-order models already carry

### 4.1 Lumped models are a graph of blocks, and steady flow is a graph Laplacian

**Lumped (0D) models.** These are ODE or DAE systems on an electrical-analogue graph:

- **Outlets.** Two-, three- and four-element Windkessels close them ([Westerhof et al. 2009](https://doi.org/10.1007/s11517-008-0359-2)).
- **Heart.** A time-varying elastance, E(t) = P(t)/(V(t) − V0), drives them ([Suga & Sagawa 1974](https://doi.org/10.1161/01.RES.35.1.117)).
- **Closed loops.** These join heart, arteries, veins and the pulmonary circuit ([Shi et al. 2011](https://doi.org/10.1186/1475-925X-10-33)).

The current open reference implementation is svZeroDSolver, a BSD-3 C++ solver that defines models as JSON "blocks" and ships a calibrator ([Menon et al. 2025](https://doi.org/10.21105/joss.07595)).

**Steady Poiseuille flow.** On a vascular graph with conductances g = πr⁴/(8μL), steady flow is exactly a weighted graph-Laplacian problem, Kirchhoff's current law at the nodes. physis's discretisation already solves flows with conductances, so the first verification costs almost nothing.

**Microvessels.** For microvessels, a constant viscosity is wrong by factors that depend on diameter and haematocrit. The in vivo law of Pries and colleagues is the standard correction, and it belongs in the reduced models rather than in 3D, because vessels below about 0.3 mm are not resolved by CT or MRI ([Pries et al. 1994](https://doi.org/10.1161/01.RES.75.5.904)).

### 4.2 One-dimensional models are a 2×2 balance law, closed by a tube law

**Governing equations.** On each edge, area A and flow Q obey ∂A/∂t + ∂Q/∂z = 0 and ∂Q/∂t + ∂(αQ²/A)/∂z + (A/ρ)∂P/∂z = −2μ(γ+2)Q/(ρA). They are closed by a tube law P(A), elastic or with a Kelvin–Voigt viscous term ([openBF, Benemerito et al. 2024](https://doi.org/10.1088/1361-6579/ad9663)). Junctions impose conservation of mass and continuity of total pressure through the characteristics, and outlets are closed by Windkessels or structured trees.

**Numerical schemes.** Six schemes were compared on the same cases: discontinuous Galerkin, locally conservative Galerkin, Galerkin least-squares, finite volumes, MacCormack and a simplified trapezium rule. **All agreed well** ([Boileau et al. 2015](https://doi.org/10.1002/cnm.2732)). Their errors against 3D were at most **1.2% in pressure, 2.6% in flow and 4.3% in area**. Against the in vitro 37-artery silicone network they were 4.0% in pressure and 25.6% in flow ([Jin & Alastruey 2021](https://doi.org/10.1098/rsif.2020.0881)).

**Taper.** openBF's largest discrepancies against the population database come from heavily tapered vessels that it does not taper. Image-derived edges therefore need A0(z) and β(z), and with them a well-balanced scheme ([Müller et al. 2013](https://doi.org/10.1016/j.jcp.2013.01.050)).

### 4.3 The venous tube law keeps the system hyperbolic through collapse

Veins differ from arteries in three ways a solver must carry:

- a tube law valid through collapse;
- valves that make flow one-directional;
- hydrostatic gravity with external (muscle) pressure.

The classical mechanics of collapsible tubes comes from [Shapiro 1977](https://doi.org/10.1115/1.3426281), the review by [Grotberg & Jensen 2004](https://doi.org/10.1146/annurev.fluid.36.050802.121918), and the roll waves of gravity-driven collapsible flow ([Brook, Falle & Pedley 1999](https://doi.org/10.1017/S0022112099006084)).

**The venous tube law.** The law now in use, from the arterial–venous network ADAVN after Müller and Toro, is p_tm = K[(A/A0)^m − (A/A0)^n] + p0, plus a viscoelastic term ([Müller et al. 2023](https://doi.org/10.3389/fphys.2023.1162391)). Its exponents encode the wall's behaviour:

- m = 10 describes passive stiffening;
- n = −3/2 describes collapse;
- the dural sinuses, which do not collapse, use m = 1/2, n = 0;
- the reference pressure p0 is about 5 mmHg.

**A short derivation from that law.** Write α = A/A0. The wave speed is c = √((K/ρ)(mα^m − nα^n)). Because n < 0, both terms are positive for every α > 0, so **the system stays strictly hyperbolic through collapse**. The speed is smallest at α* = (n²/m²)^(1/(m−n)) ≈ 0.72, which is where flow limitation (speed index u/c → 1) should appear. The choice (m, n) = (1/2, 0) recovers the elastic √A law of arteries. **One tube-law class with per-edge (K, m, n) therefore covers arteries, veins and sinuses in a single code.**

**Valves and collapse at boundaries.** Both are lumped elements.

- **Valves.** The valve model of [Mynard et al. 2012](https://doi.org/10.1002/cnm.1466) serves heart valves and venous valves alike, as a junction whose loss depends on a state variable driven by the pressure difference. ADAVN uses it for both.
- **Collapse at the sinuses.** Where dural sinuses meet cerebral veins, Starling-resistor elements obey dζ/dt = K_so(1 − ζ)Δp when opening and K_sc·ζ·Δp when closing. The "venous waterfall" is then driven by intracranial pressure ([Müller et al. 2023, eqs. 41–42](https://doi.org/10.3389/fphys.2023.1162391)).
- **Pressure shielding.** A valve substantially reduces the dynamic pressure on the vein below it when gravity acts ([Zervides et al. 2008](https://doi.org/10.1186/1475-925X-7-8)).

**The calf muscle pump.** Two 1D models treat it.

- [Keijsers et al. 2015](https://doi.org/10.1002/cnm.2714) model contraction as raised extravascular pressure on deep veins. In their model, **proximal valves increase effective venous return by 53%**, and the superficial system cuts refilling time by 37%.
- [Ferlito et al. 2026](https://doi.org/10.1109/TBME.2026.3701594) extend this to deep and superficial veins joined by perforators, collapsible mechanics, dynamic valves, gravity and a 5 s contraction cycle. They compare competent with incompetent valves and external compression. Only the abstract was available, so their equations and code status are unknown.

**Whole circulations with veins.** Three models exist.

- Müller and Toro's global closed loop characterised the head and neck veins from MRI and validated them against phase-contrast flows ([Müller & Toro 2014](https://doi.org/10.1002/cnm.2622)).
- Mynard and Smolich modelled the entire adult circulation, including portal veins ([Mynard & Smolich 2015](https://doi.org/10.1007/s10439-015-1313-8)).
- ADAVN couples **2,185 arteries to 189 veins** and publishes its network as a CSV file ([Müller et al. 2023](https://doi.org/10.3389/fphys.2023.1162391)).

**Posture.** Gravity enters the momentum equation as a source term along each edge's axis. The image-derived graph must therefore carry each edge's direction in the body frame and a posture rotation; this is an inference from the cited models, which treat posture explicitly. Further, the area seen in a supine image is the area at the imaged transmural pressure, not A0. A0 must therefore be inferred through the tube law from an assumed or measured venous pressure; that too is an inference, not a published recipe.

### 4.4 A verification ladder exists for arteries; for veins only exact solutions anchor it

**The arterial ladder** climbs from exact to population-scale:

- **Analytic solutions:**
  - Poiseuille and Kirchhoff on a graph;
  - Womersley's oscillatory flow ([Womersley 1955](https://doi.org/10.1113/jphysiol.1955.sp005276)).
- **Exact Riemann solutions:**
  - collapsible tubes with discontinuous properties ([Toro & Siviglia 2013](https://doi.org/10.4208/cicp.210611.240212a));
  - veins including collapse ([Spiller et al. 2017](https://doi.org/10.1016/j.amc.2017.01.024));
  - a general momentum-correction coefficient ([Spilimbergo et al. 2024](https://doi.org/10.4208/cicp.OA-2023-0250)).
- **Boileau's six cases.** Their official data location was not found. openBF's repository holds reference solutions for the carotid, the iliac bifurcation and ADAN56 under Apache-2.0 ([openBF tests](https://github.com/INSIGNEO/openBF/tree/master/test)).
- **The in vitro 37-artery network.**
- **The Pulse Wave Database.** 4,374 virtual adults on a 116-segment network ([Charlton et al. 2019](https://doi.org/10.1152/ajpheart.00218.2019)). It is dedicated to the public domain under **ODC-PDDL** ([Zenodo](https://zenodo.org/records/3275625)), the only reference dataset found that an MIT repository can redistribute without restriction.

**For veins there is no Boileau-style dataset.** The anchors are:

- the exact Riemann solutions with collapse;
- Brook, Falle and Pedley's steady states and roll waves;
- the published waveforms of Müller and Toro (against PC-MRI), Keijsers and Ferlito, as semi-quantitative targets.

This gap is also an opportunity. A public venous verification suite would be new.

### 4.5 The same outlet object should close a 1D edge and a 3D boundary

**Coupling.** A 3D region is closed by lumped or 1D models through an outlet that returns pressure P(Q) and dP/dQ, coupled implicitly and modularly ([Esmaily Moghadam et al. 2013](https://doi.org/10.1016/j.jcp.2012.07.035)). The theory is reviewed by [Quarteroni et al. 2016](https://doi.org/10.1016/j.cma.2016.01.007).

**Automatic reduced models.** Pfaller and colleagues built 0D and 1D models automatically from 3D geometry, detecting junctions and stenoses. On 72 public models their terminal-outlet errors were typically **1–10%** against rigid-wall 3D, and 0D cost about a third of 1D ([Pfaller et al. 2022](https://doi.org/10.1002/cnm.3639)). That is the closest published analogue of the target pipeline: geometry, then centreline graph, then reduced model, then comparison with 3D.

**Learned components.** Since 2023 they have entered as closures and emulators, not replacements:

- graph networks predict centreline pressure and flow within about 3% ([Pegolotti et al. 2024](https://doi.org/10.1016/j.compbiomed.2023.107676));
- hybrid models predict bifurcation losses better than classical loss models ([Rubio et al. 2025](https://doi.org/10.1016/j.compbiomed.2024.109420));
- physics-informed neural operators support Bayesian inference from sparse 4D flow in Fontan patients ([Ryan et al. 2026](https://doi.org/10.1007/s10439-026-04219-1)).

Each needs a physics solver to generate its training data, so the solver comes first. **No learned surrogate specific to venous networks was found.**

## 5. Lattice Boltzmann suits voxel masks, and validation sets the bar for every 3D solver

### 5.1 Image-based CFD practice is settled; its weak points are inflow and outflow

The standard pipeline is to segment, build a surface and mesh it, prescribe an inlet waveform, couple the outlets to 0D or 1D models, and solve the incompressible Navier–Stokes equations with backflow stabilisation. Each boundary choice has a measured effect:

- **Outlets.** Constant-pressure and prescribed-velocity outlets distort solutions strongly ([Vignon-Clementel et al. 2006](https://doi.org/10.1016/j.cma.2005.04.014)). In a normal aorta, well-tuned three-element Windkessels at every outlet reproduced in vivo flows and physiological pressures, while zero-pressure outlets did not ([Pirola et al. 2017](https://doi.org/10.1016/j.jbiomech.2017.06.005)).
- **Backflow.** Of three outlet treatments compared, a boundary stabilisation term had the least effect on haemodynamics and added almost no cost ([Esmaily Moghadam et al. 2011](https://doi.org/10.1007/s00466-011-0599-0)). It should be the default.
- **Inflow profile.** Idealised flat or parabolic profiles distort helical flow in the aorta, where measured profiles do not ([Morbiducci et al. 2013](https://doi.org/10.1016/j.jbiomech.2012.10.012)).
- **Transitional flow.** Spectral-element DNS found velocity fluctuations at 100–300 Hz in a stenosed carotid ([Lee et al. 2008](https://doi.org/10.1016/j.jbiomech.2008.03.038)). The public DNS database of 75% stenoses ([Varghese et al. 2007](https://doi.org/10.1017/S0022112007005848); [Purdue database](https://engineering.purdue.edu/CFDLAB/DNS_stenoticflow_database/main.html)) is a verification target.

### 5.2 Lattice Boltzmann maps the mask to the lattice, if the wall is interpolated

**Why it suits voxels.** The lattice Boltzmann method takes the fluid domain directly from the image mask, with no unstructured meshing. It thereby sidesteps the copyleft meshers of §3.4 entirely. It is local and bound by memory bandwidth, so it maps well to GPUs.

**Wall treatment.** Accuracy hinges on how the wall is treated. In HemeLB's tests on Poiseuille, Womersley and curved-pipe flow, the interpolated Bouzidi–Firdaouss–Lallemand (BFL) wall **converges at second order while simple bounce-back degrades to first**, and BFL is cheaper than the Guo–Zheng–Shi alternative ([Nash et al. 2014](https://doi.org/10.1103/PhysRevE.89.023303)). The sub-voxel wall distance BFL needs comes directly from a signed distance field. The two-relaxation-time collision with its "magic" parameter makes wall positions independent of viscosity ([Ginzburg et al. 2008](https://doi.org/10.4208/cicp.2008.v3.p427)).

**Accuracy against finite volumes.** Across 462 boundary-condition profiles on two aneurysms, HemeLB and a commercial finite-volume solver agreed quantitatively in the momentum-driven MCA aneurysm but less in the shear-driven ICA aneurysm. Both located the same extreme rupture-risk metrics ([Benemerito et al. 2025](https://doi.org/10.1016/j.compbiomed.2025.111050)).

**Licences.** These rule out most open codes as components:

| Code | Licence | Usable from MIT code? |
|---|---|---|
| HemeLB | LGPL-3.0 | yes, with modest obligations ([HemeLB](https://github.com/hemelb-codes/hemelb)) |
| Palabos | AGPL-3.0 | no ([Palabos](https://palabos.unige.ch/)) |
| OpenLB | GPL-2.0-or-later | no ([OpenLB](https://www.openlb.net/download/)) |
| waLBerla | GPL-3.0 | no ([waLBerla](https://www.walberla.net/)) |
| FluidX3D | custom; forbids commercial use | no ([FluidX3D](https://github.com/ProjectPhysX/FluidX3D)) |

An own D3Q19 lattice Boltzmann solver, with TRT or MRT collision, BFL walls and Windkessel outlets imposed through density as HemeLB does, is therefore both scientifically apt and legally necessary. Its known costs are weak compressibility of order Ma², a viscosity tied to the relaxation time, and Δt ∝ Δx² under diffusive scaling.

### 5.3 Rheology matters in slow flow, which is where veins live

**Large arteries.** In large arteries, where shear rates exceed about 100 s⁻¹, a Newtonian viscosity near 3.5 mPa·s is usually adequate. In high-resolution aneurysm simulations, a shear-thinning model barely changed WSS or the spectral content of instabilities, while solution strategy mattered more. Point-wise ratios of non-Newtonian to Newtonian WSS also exaggerate differences in low-shear zones ([Khan et al. 2017](https://doi.org/10.1002/cnm.2836)).

**The answer depends on the parameters.** With four published Carreau–Yasuda sets in a right coronary artery, the Cho–Kensey set correlated with Newtonian at 0.99. The Gijsen set produced point-wise OSI differences from **+25% to −250%** ([Kandangwa et al. 2022](https://doi.org/10.3389/fbioe.2022.962687)).

**Veins.** In an IVC–iliac venous model, non-Newtonian viscosity raised WSS, and the Newtonian model overestimated near-wall and luminal transport ([Lynch et al. 2022](https://doi.org/10.1038/s41598-022-19867-1)). This is the evidence that matters most for veins.

**Design.** Rheology should be a strategy object, μ(γ̇, H), so that the following are interchangeable:

- Newtonian;
- Carreau–Yasuda, μ = μ∞ + (μ0 − μ∞)[1 + (λγ̇)^a]^((n−1)/a);
- Cross;
- regularised Casson;
- Quemada with haematocrit.

Results should be reported for at least two parameter sets.

### 5.4 Rigid walls get pressure and high shear right, and get oscillatory shear wrong

**Two FSI families.** Fluid–structure interaction comes in two main forms:

- **Coupled momentum method.** A membrane wall on a fixed fluid mesh; cheap, robust, small-deformation ([Figueroa et al. 2006](https://doi.org/10.1016/j.cma.2005.11.011)).
- **ALE.** Full arbitrary Lagrangian–Eulerian FSI, needed for large motions.

**Added mass.** Partitioned schemes become unstable at blood-like density ratios because of added mass, so they need strongly coupled or Robin-type interfaces ([Causin et al. 2005](https://doi.org/10.1016/j.cma.2004.12.005)).

**What compliance changes.** It changes OSI, low WSS and recirculation, sometimes enormously:

- **Aneurysms.** Across 101 aneurysms mounted on one parent artery, the compliant wall changed sac-averaged OSI by **−36% to +674%** relative to the rigid one ([AnXplore, Goetz et al. 2024](https://doi.org/10.3389/fbioe.2024.1433811)).
- **Dissection.** In aortic dissection, lowering flap stiffness raised flap displacement from 1.4 to 13.4 mm and cut abdominal true-lumen flow by up to 20% ([Bäumler et al. 2020](https://doi.org/10.1007/s10237-020-01294-8)).

**Order of work.** Rigid walls first, then the coupled momentum method, then ALE, kept for flaps and valves. Every report should say which indices were computed with rigid walls.

### 5.5 Haemodynamic indices follow from three running integrals per wall node

**Definitions.** The standard indices are functionals of the wall shear stress vector τ_w over a cycle T:

| Index | Definition | Source |
|---|---|---|
| TAWSS | (1/T)∫\|τ_w\| dt | standard |
| OSI | ½[1 − \|∫τ_w dt\|/∫\|τ_w\| dt] | [Ku et al. 1985](https://doi.org/10.1161/01.ATV.5.3.293) |
| RRT | 1/[(1 − 2·OSI)·TAWSS] | attributed to [Himburg et al. 2004](https://doi.org/10.1152/ajpheart.00897.2003) |
| transWSS | the cycle mean of shear transverse to the mean WSS direction | [Peiffer et al. 2013](https://doi.org/10.1016/j.jbiomech.2013.08.003) |
| ECAP | OSI/TAWSS | — |

**Correlations.** TAWSS, OSI and RRT are strongly inter-correlated, and RRT has been proposed as a single robust metric ([Lee, Antiga & Steinman 2009](https://doi.org/10.1115/1.3127252)).

**Computation.** All of them follow from running integrals of τ_w and |τ_w| per wall node, plus a second pass for transWSS, so they can be accumulated while the solver streams. The WSS itself comes from recovered gradients in finite elements, or from the non-equilibrium stress in lattice Boltzmann.

### 5.6 Validation is a ladder of analytic cases, the FDA nozzle and community challenges

**The FDA benchmark nozzle.** It is the standard validation case for transitional flow. A 12 mm tube leads into a 4 mm throat, measured by PIV in three laboratories at throat Reynolds numbers of 500–6,500 ([Hariharan et al. 2011](https://doi.org/10.1115/1.4003440)).

- **Measurement uncertainty.** Interlaboratory uncertainty was about **10%** in the laminar and turbulent cases, but about **60%** for jet size and velocity at Re 2,000.
- **CFD performance.** Discrepancies among more than 20 CFD groups concentrated in separated flow ([Malinauskas et al. 2017](https://doi.org/10.1097/MAT.0000000000000499)).
- **Lattice Boltzmann.** LBM-based DNS matched the PIV acceptably at Re 2,000 and 3,500 ([Jain 2020](https://doi.org/10.1007/s11517-020-02188-8)).
- **Data.** The geometry, CAD and data are distributed through the FDA's Regulatory Science Tools ([Zenodo 6629049](https://zenodo.org/records/6629049/files/CFD%20Tool.pdf)).

**Credibility.** For anything near a clinical context, ASME V&V 40-2018 and the FDA's final guidance of November 2023 set a risk-informed credibility framework for physics-based models ([Federal Register](https://www.govinfo.gov/content/pkg/FR-2023-11-17/html/2023-25470.htm)).

**Venous 3D work.** Venous applications are emerging, but all are small:

- **Iliac vein compression.** In four patients and four controls, shear rate in the left common iliac vein was **550 ± 103 s⁻¹ against 75 ± 37 s⁻¹**, with peaks near 1,650 s⁻¹, above a platelet-activation threshold of about 1,000 s⁻¹. The authors call this hypothesis-generating ([Assi et al. 2024](https://doi.org/10.3389/fbioe.2024.1302063)).
- **Fontan circulation.** Its efficiency falls with exercise, with respiration modelled ([Marsden et al. 2007](https://doi.org/10.1007/s10439-006-9224-3)).

Venous 3D models need respiration-driven inflow, low-shear rheology and often collapsible walls. They are therefore the hardest test of the 0D coupling and rheology plug-ins.

## 6. Transport and chemistry can be held to closed forms and public-domain reaction networks

### 6.1 Oxygen needs a saturation law, a slope-weighted transport equation and thin boundary layers

**Saturation laws.** The oxygen–haemoglobin saturation laws range from rational functions to full physiological models:

- **Hill**, S = Pⁿ/(Pⁿ + P50ⁿ), and Adair's four-step binding: closed forms for tests.
- **Severinghaus.** At standard conditions, **S = [23400/(P³ + 150P) + 1]⁻¹**, with P in mmHg. It is corrected to a "virtual" PO2 for temperature, pH and CO2 ([Severinghaus 1979](https://doi.org/10.1152/jappl.1979.46.3.599)).
- **Dash–Bassingthwaighte.** This model adds CO2, pH, 2,3-DPG and temperature, and should be implemented from its corrected 2010 text ([Dash & Bassingthwaighte 2010](https://doi.org/10.1007/s10439-010-9948-y)). Their 2016 comparison of all these models is the natural cross-model fixture ([Dash et al. 2016](https://doi.org/10.1007/s00421-015-3228-3)).

**Large vessels.** Oxyhaemoglobin moves with the plasma, so in large vessels the transport equation for PO2 takes the form (1 + ([Hb]C_Hb/α)·dS/dP)(∂P/∂t + u·∇P) = D∇²P. This is the form used in the Moore–Ethier line, reconstructed in the notes ([Moore & Ethier 1997](https://doi.org/10.1115/1.2798295); [Tarbell 2003](https://doi.org/10.1146/annurev.bioeng.5.040202.121529)). The factor in brackets can exceed 50 near P50. With a Schmidt number for oxygen near 1.7 × 10³ (the notes' estimate), concentration boundary layers are a few percent of the velocity layer's thickness, and the near-wall mesh must allow for them.

**The microcirculation.** Several methods treat it:

- **Krogh cylinder.** It gives an exact fixture, P(r) = P_c + (M/4K)(r² − r_c²) − (MR²/2K)·ln(r/r_c) ([Krogh 1919](https://doi.org/10.1113/jphysiol.1919.sp001839)).
- **Green's functions.** These handle realistic networks ([Secomb et al. 2004](https://doi.org/10.1114/B:ABME.0000049036.08817.44)), with the phase-separation laws for haematocrit at bifurcations ([Pries et al. 1990](https://doi.org/10.1161/01.RES.67.4.826)).
- **A distributional fixture.** In extensive cortical networks, **99% of tissue lies between 20 and 80 mmHg** ([Celaya-Alcala et al. 2021](https://experts.arizona.edu/en/publications/simulation-of-oxygen-transport-and-estimation-of-tissue-perfusion/)).
- **Plasma convection.** Recent red-cell-resolved work finds plasma convection negligible in typical capillaries, which justifies diffusion-only plasma gaps in reduced models ([Zhang 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12009512)).

### 6.2 Wall transport is a three-layer problem with concentration polarisation at its centre

The standard model has three layers:

- **Lumen.** The lumen is solved at very high Schmidt numbers: for LDL the notes estimate about 7 × 10⁵, giving sub-micrometre to micrometre boundary layers.
- **Endothelium.** The endothelium is a Kedem–Katchalsky membrane, J_v = L_p(Δp − σΔπ) and J_s = PΔc + (1 − σ)J_v·c̄ ([Kedem & Katchalsky 1958](https://doi.org/10.1016/0006-3002(58)90330-5)). Alternatively it is a three-pore model whose leaky junctions depend on WSS ([Olgac et al. 2008](https://doi.org/10.1152/ajpheart.01082.2007)).
- **Wall.** The wall is porous, with advection, diffusion and reaction.

**Concentration polarisation.** LDL accumulates at the luminal surface because water filters through the wall. This is the central coupled effect.

**Patient-specific growth.** The latest work couples LDL transport to agent-based plaque growth on a patient-specific coronary. Calibrated against one-year follow-up, it predicted the change in wall area within 10% in two cross-sections ([Corti et al. 2026](https://pure.eur.nl/en/publications/personalized-multiscale-modeling-of-coronary-plaque-progression-t/)).

**Numerics.** Reduced near-wall models exist because boundary layers this thin are expensive to resolve ([Hansen & Shadden 2016](https://doi.org/10.1007/s10237-015-0719-4)).

### 6.3 Coagulation has curated public-domain fixtures; venous thrombosis has an unmodelled link

**The well-mixed coagulation lineage** gives the best fixtures in the whole field.

- **Hockin 2002.** The model of [Hockin et al. 2002](https://doi.org/10.1074/jbc.M201173200) has **34 ODEs and 42 rate constants**. It is curated in BioModels as BIOMD0000000335 by reproducing the paper's figures, with a COPASI run on BioSimulations ([BioSimulations](https://biosimulations.org/projects/BIOMD0000000335_copasi_A_model_for_the_stoichiometric)).
- **Other curated models** in the same repository include:
  - Jones–Mann 1994 (BIOMD0000000336);
  - Butenas 2004 (BIOMD0000000362);
  - Wajima 2009 (BIOMD0000000338–340);
  - a nine-species automatic reduction (BIOMD0000000755) ([Hansen & Shadden 2019](https://doi.org/10.1002/cnm.3220)).
- **Licence.** BioModels releases every model under **CC0** ([BioModels terms](https://www.ebi.ac.uk/biomodels/termsofuse)). The models are therefore compatible with an MIT repository and can sit in its test fixtures.

**Falsification data.** These come from calibrated automated thrombography ([Hemker et al. 2003](https://doi.org/10.1159/000071636)). A 2023 critical evaluation found that both standard schemes needed repair against factor-deficiency data, and it publishes its code and data ([Ranc et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10461854/)):

- **The intrinsic model** generated thrombin with no factor XII, because FVIIa activated FIX and FX without tissue factor.
- **The extrinsic model** was too slow and needed the IXa→FX k_cat raised **60-fold**.
- **Both** still failed at FVIII below 1%.

The debate on what such validation means should be read before claiming it ([Hemker et al. 2012](https://doi.org/10.1111/j.1538-7836.2012.04802.x)).

**Flow-coupled lineage.** It runs:

1. from a well-mixed zone with mass-transfer coefficients ([Kuharsky & Fogelson 2001](https://doi.org/10.1016/S0006-3495(01)76085-7));
2. to spatial models with a Brinkman resistance from platelet aggregates ([Leiderman & Fogelson 2011](https://doi.org/10.1093/imammb/dqq005));
3. to the open OpenFOAM-based clotFoam, which inherits OpenFOAM's GPL and so serves only as an external comparison ([Montgomery et al. 2023](https://arxiv.org/abs/2304.09180)).

**Venous thrombosis.** It differs in low shear, long residence time and the dominance of fibrin over platelets. Its models and anatomy so far:

- **Models.**
  - A Hockin-type cascade under venous flow ([Govindarajan et al. 2016](https://www.sciencedirect.com/science/article/pii/S0006349516300595)).
  - A model of platelet count, shear and injury size under venous flow, available in BioModels ([Bouchnita et al. 2020](https://doi.org/10.1371/journal.pone.0235392)).
  - Fibrin-enhanced low-shear models ([Méndez Rojano et al. 2022](https://doi.org/10.1371/journal.pcbi.1010277)).
- **Patient-specific deep-vein work.**
  - The common femoral vein, with thrombin–fibrinogen reactions ([Jimoh-Taiwo et al. 2022](https://doi.org/10.3389/fphy.2022.886193)).
  - The effect of valve aperture ([George & Ngoepe 2026](https://doi.org/10.1016/j.premed.2026.100055)).
- **The valve pocket as the site of origin.** A vortex is trapped in the valve sinus ([Karino & Motomiya 1984](https://doi.org/10.1016/0049-3848(84)90224-X)). Stasis there desaturates haemoglobin and makes the endothelium hypoxic and procoagulant ([Bovill & van der Vliet 2011](https://doi.org/10.1146/annurev-physiol-012110-142305)).

**No model was found that couples flow in the valve pocket, oxygen desaturation and the coagulation cascade as one advection–diffusion–reaction problem.** The infrastructure's oxygen and coagulation modules, on one venous geometry with a valve, would address this directly; it is a natural research target.

### 6.4 Contrast, drugs and perfusion share an exact invariant and open benchmarks

**Indicator dilution.** It rests on the Meier–Zierler central-volume theorem, V = F·t̄, where t̄ is the mean transit time ([Meier & Zierler 1954](https://doi.org/10.1152/jappl.1954.6.12.731)). Every simulated contrast curve, whether 0D, 1D, 3D or porous, must satisfy it on its residue function. **That gives an exact, model-independent invariant** for the transport code.

**Tracer kinetics.** Several standard models build on that theory:

- **Tofts.** The extended Tofts model, C_t = v_p·C_p + K^trans ∫C_p(τ)e^(−k_ep(t−τ))dτ ([Tofts et al. 1999](https://doi.org/10.1002/(SICI)1522-2586(199909)10:3%3C223::AID-JMRI2%3E3.0.CO;2-S)), driven by Parker's closed-form population arterial input function ([Parker et al. 2006](https://doi.org/10.1002/mrm.21066)), gives a synthetic fixture.
- **DSC deconvolution.** This is ill-posed and needs SVD regularisation ([Østergaard et al. 1996](https://doi.org/10.1002/mrm.1910360510)).
- **Whole-body contrast.** It is a 0D compartment network ([Bae et al. 1998](https://doi.org/10.1148/radiology.207.3.9609886)).

**Software variability.** In the OSIPI DCE challenge, ten submissions scored between **28% and 78%** on its gold-standard score, with high variability in K^trans between packages ([Shalom et al. 2024](https://eprints.whiterose.ac.uk/id/eprint/207709/)). A digital phantom for CT and MR perfusion offers a second reference ([Kudo et al. 2013](https://doi.org/10.1148/radiol.12112618)).

**Perfusion.** It is modelled as multi-compartment Darcy flow ([Michler et al. 2013](https://doi.org/10.1002/cnm.2520)). A 1D coronary model coupled to a porous myocardium "compare[d] overall well" with ¹⁵O-water PET in five patients ([Papamanolis et al. 2021](https://doi.org/10.1007/s10439-020-02681-z)).

**Drug-eluting stents.** Models of drug release show that binding cannot be ignored. Paclitaxel is convection-dominated and sirolimus binding-dominated, and saturable two-phase binding must be modelled ([Bozsak et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4470631/)). No peer-reviewed model from 2024–2026 was found.

### 6.5 Stabilised transport plus split stiff reactions is the standard numerical recipe

**Transport.** At Péclet numbers of 10⁵ to 10⁹, Galerkin finite elements oscillate. SUPG ([Brooks & Hughes 1982](https://doi.org/10.1016/0045-7825(82)90071-8)) and finite volumes are the remedies.

**Reactions.** Coagulation has a long lag followed by a burst, with rate constants spanning many orders of magnitude, so it is stiff. It needs BDF or Rosenbrock methods; ROS3P is third-order and L-stable ([Lang & Verwer 2001](https://doi.org/10.1023/A:1021900219772)).

**Coupling.** Strang splitting couples the two ([Strang 1968](https://doi.org/10.1137/0705041)), with a convergence test that halving Δt cuts the splitting error fourfold.

**Residence time.** For low-shear thrombosis, residence time is computed as an Eulerian advection–diffusion field ([Esmaily-Moghadam et al. 2013](https://doi.org/10.1063/1.4819142)).

**Fit with the infrastructure.** The recipe keeps the law apart from the field and the integrator, as physis's `continuum` already does for heat conduction.
- **What calx already holds:** the θ-method, conjugate gradients and small dense solves.
- **What it still needs:** a Rosenbrock or BDF integrator.

## 7. Imaging and simulation meet through an explicit measurement operator

### 7.1 4D flow MRI must be compared with simulation in the space of the measurement

**The 2023 consensus update** sets the acquisition targets ([Bissell et al. 2023](https://doi.org/10.1186/s12968-023-00942-z)):

| Parameter | Recommendation |
|---|---|
| Acquired voxel, adult vessels | aim for 2 mm, never coarser than 2.5 mm |
| Temporal resolution | 30 ms, never coarser than 50 ms |
| VENC | about 10% above the peak velocity, never more than 25% above |
| VENC for venous and Fontan flow | 50–80 cm/s |

It also fixes a processing order: background-phase correction, anti-aliasing, segmentation, visualisation, quantification. Reports must state both the acquired and the reconstructed resolution.

**Derived quantities carry their own errors.**

- **WSS from the image.** It is biased low and depends on resolution. A learned estimator trained on CFD of 37 aortas reached 4.3% relative error against CFD but 7.1% on noisy synthetic 4D flow ([WSSNet, Ferdian et al. 2022](https://doi.org/10.3389/fcvm.2021.769927)).
- **Relative pressure.** It needs a Poisson, work–energy or virtual work–energy estimator ([Nolte et al. 2021](https://doi.org/10.1016/j.media.2021.102195)).
- **Image quality.** An image can be screened by its compatibility with the Navier–Stokes equations ([Garay et al. 2022](https://doi.org/10.1002/cnm.3603)).

**The decisive validation result.** In vitro, raw 4D flow and CFD agreed at only **r² = 0.63**. Agreement rose to **0.97** once MR artefacts were corrected and the CFD was degraded to the MRI's spatial and temporal resolution ([Puiseux et al. 2019](https://doi.org/10.1002/nbm.4063)).

**The measurement operator H.** The interface between simulation and image should therefore be an explicit forward operator H. It maps a simulated velocity field to the image the scanner would produce, in four steps:

1. average over the acquired voxel and cardiac bin;
2. wrap at the VENC;
3. add noise that scales with VENC and SNR;
4. add background-phase offsets.

Bloch simulators such as CMRsim are its high-fidelity form ([Weine et al. 2024](https://doi.org/10.1002/mrm.30010); [Puiseux et al. 2021](https://doi.org/10.1371/journal.pone.0248816)).

**Open data are scarce.** The one licence confirmed for an open aortic 4D flow atlas is **CC BY-NC-SA 3.0 US** ([Northwestern atlas](https://prism.northwestern.edu/records/r29kc-x5336)). It belongs in the data cache, never in the repository. No open venous 4D flow dataset was found.

### 7.2 Assimilation splits into ensemble filters, adjoints and Bayesian calibration on surrogates

**Sequential filtering.** The reduced-order unscented Kalman filter works in parameter space. It "does not require any adjoint or tangent problems" and runs its 2p + 1 sigma-point simulations in parallel ([Moireau & Chapelle 2011](https://doi.org/10.1051/cocv/2010006); [Bertoglio et al. 2012](https://doi.org/10.1002/cnm.1476)). Run on 1D models of an in vitro network with real measurements, it estimated terminal resistances and wall parameters from at most one pressure or flow signal per vessel. That study also documented the filter's limits in clinically realistic setups ([Caiazzo et al. 2017](https://doi.org/10.1002/cnm.2843)).

The infrastructure already holds the filter's analysis step: tyche conditions a Gaussian one observation at a time. ROUKF reduces to that step plus a scheduler of forward runs. CRIMSON's implementation cannot be reused (AGPL with the Commons Clause), so ROUKF should be reimplemented from Moireau and Chapelle.

**Variational assimilation.** It needs a differentiable solver and suits whole-state or inflow-field estimation:

- transient Navier–Stokes recovered from noisy, low-resolution PC-MRI in cerebral aneurysms ([Funke et al. 2019](https://doi.org/10.1002/cnm.3152));
- a discontinuous Galerkin solver fitted to 4D flow and checked against PIV and 7T in vivo data ([Töger et al. 2020](https://doi.org/10.1002/mrm.28269)).

**Bayesian calibration on surrogates.** This is the current frontier. Richter and colleagues run 3D once, build an optimised 0D surrogate, and obtain Windkessel posteriors by sequential Monte Carlo. Optimising the surrogate cut its median error against 3D by **nearly an order of magnitude** across 72 public models ([Richter et al. 2025](https://doi.org/10.1098/rsta.2024.0223)).

**Geometry as a parameter.** Bayesian inverse Navier–Stokes problems now make the vessel boundary itself a parameter inferred jointly with the flow. One reconstructed and segmented MR velocimetry at SNR ≈ 6 with negligible error, which the authors equate to a 27-fold shorter scan ([Kontogiannis et al. 2022](https://doi.org/10.1017/jfm.2022.503); [Kontogiannis et al. 2025](https://doi.org/10.1088/1361-6420/ad9cb7)). This is the formal answer to §1's finding that geometry is the dominant error.

### 7.3 Physics-informed networks are inverse estimators, not faster forward solvers

**Where PINNs have earned their place.** In vascular flow they work as inverse solvers and filters:

- recovering velocity and pressure from a transported scalar ([Raissi et al. 2020](https://doi.org/10.1126/science.aaw4741));
- recovering near-wall WSS from sparse data with unknown boundary conditions ([Arzani et al. 2021](https://doi.org/10.1063/5.0055600)), though its code carries no licence;
- predicting pressure from 4D flow with 1D physics ([Kissas et al. 2020](https://doi.org/10.1016/j.cma.2019.112623)).

**Learned surrogates.** They report few-percent errors in distribution. Shape-model-driven networks reached 6.0% in pressure and 4.0% in velocity, about 4,000 times faster than the solver ([Pajaziti et al. 2023](https://doi.org/10.1371/journal.pcbi.1011055)).

**Independent assessments temper this.**

- Weak baselines and reporting biases produce overoptimism in machine learning for fluid PDEs ([McGreivy & Hakim 2024](https://doi.org/10.1038/s42256-024-00897-5)).
- PINNs rarely beat finite elements on forward problems ([Grossmann et al. 2024](https://doi.org/10.1093/imamat/hxae011)).
- Generalisation to unseen anatomies remains open, as a benchmark of geometric deep learning on coronary meshes shows ([Nannini et al. 2025](https://doi.org/10.1016/j.compbiomed.2025.110477)).

**Placement in the infrastructure.**

- Learned components belong in the Python works.
- Their errors must be carried in the likelihood of any calibration.
- MIT-licensed code can be ported: 4DFlowNet, WSSNet, neuraloperator.

### 7.4 Digital twins and credibility are defined; atlases supply templates, not priors

**Definitions.** NASEM defines a digital twin as a computational model of a physical counterpart, "dynamically updated through bidirectional data flows" ([NASEM 2024](https://www.nationalacademies.org/publications/26894)). For imaging-constrained haemodynamics that bidirectional flow is the assimilation layer of §7.2. The cardiovascular framing comes from [Corral-Acero et al. 2020](https://doi.org/10.1093/eurheartj/ehaa159) and [Niederer et al. 2021](https://doi.org/10.1038/s43588-021-00072-5).

**Credibility for patient-specific models.** Building on V&V 40, the credibility case adds two analyses:

- an analysis of whether the validation evidence applies to each new patient ([Pathmanathan et al. 2017](https://doi.org/10.1115/1.4037671));
- an account of how far the model is personalised ([Galappaththige et al. 2022](https://doi.org/10.1371/journal.pcbi.1010541)).

Because V&V 40 is risk-informed, an open infrastructure should record a **context of use** with each assimilated model, not one global accuracy figure.

**Atlases.** Normal-anatomy atlases supply labelled geometry and terminology, and no work was found that uses them directly as priors for haemodynamics. The working prior is a statistical shape model built from a cohort ([Bruse et al. 2016](https://doi.org/10.1186/s12880-016-0142-z)).

- **BodyParts3D.** The licence recorded in the integration notes needs correcting. Those notes, quoting the 2011 archive README, give CC BY-SA 2.1 Japan ([DBCLS archive README](https://dbarchive.biosciencedbc.jp/data/bodyparts3d/20110915/README_e.html)). In fact **DBCLS relicensed BodyParts3D under CC BY 4.0 on 2025-02-25** (README_e.html, §3), so derivative meshes no longer need to stay share-alike. The OBJ files' headers still carry the earlier CC BY-SA 2.1 Japan notice. Provenance should record both the header text and the superseding licence.
- **Human Reference Atlas.** Its 3D reference objects are CC BY 4.0, the brain excepted ([CCF library](https://hubmapconsortium.github.io/ccf/pages/ccf-3d-reference-library.html)).

Either atlas can serve as a template that is registered and deformed to a patient. A shape model's Gaussian then serves as the prior on geometry in the Bayesian inverse problem.

## 8. Build coordinates, geometry and reduced-order veins before voxel flow and assimilation

### 8.1 Three principles fix the order

The evidence fixes the order of work by three principles.

1. **Dependence.** Every later stage consumes the coordinate frame and the reconstructed geometry, and geometry is where the error lives (§1).
2. **Verifiability.** Each stage should enter with exact solutions before it touches data. Every stage below has at least one closed form or exact invariant, so the infrastructure's habit of holding code to closed forms carries over intact.
3. **Leverage for veins.** The defining venous phenomena (collapse, valves, gravity, the muscle pump) are formulated and verified in one dimension. Venous 3D models, in turn, need the 0D outlets and the rheology plug-ins anyway (§5.6). Reduced-order veins therefore come before 3D.

### 8.2 The build order, with verification at each stage

| Stage | Placement | First content | Exact verification | Benchmarks and data (licence) |
|---|---|---|---|---|
| 0 | medical-imaging; eikon | DICOM IPP/IOP/Pixel Spacing → LPS affine with slice direction from successive IPPs; NIfTI qform/sform → RAS; explicit LPS↔RAS and mm↔m transformations; grids that carry frame, unit and provenance | DICOM PS3.3 eq. C.7.6.2.1-1; nifti1.h; fixtures generated by reference readers | IXI TOF headers (CC BY-SA 3.0); TotalSegmentator v2 (CC BY 4.0) |
| 1 | eikon | Gaussian derivatives; 3×3 eigensolver; Frangi, Sato, Jerman; OOF by FFT; exact Euclidean distance transform ([Maurer et al. 2003](https://doi.org/10.1109/TPAMI.2003.1177156)); fast marching; hysteresis region growing; EM mixtures; push–relabel graph cut; topology-preserving thinning; MC33 and surface nets; Taubin smoothing | synthetic tubes, tori and bifurcations of known radius; analytic Hessian of a Gaussian tube; distance and arrival time of a point source | Lamy C++ benchmark cases; ITK outputs (Apache-2.0); IXI raw TOF |
| 2 | eikon, with morphe | Dice per diameter class, clDice, β0/β1 errors and Betti matching, HD95/ASSD, Schaap overlap | hand-countable complexes; known-genus surfaces | TopCoW (attribution; commercial use needs permission); AVT (CC BY 4.0) |
| 3 | anatomy, `vessels` | vascular graph: Voronoi/eikonal and on-mesh centrelines; per-station area, contour, ellipse and inscribed radius; tortuosity; β1; fitted Murray exponent; artery/vein labels; valve sites; geometric uncertainty | analytic tubes, ellipses and bifurcations | TotalSegmentator vein classes; MSD Task08; HVA-CT; TopCoW graph edges |
| 4 | physis fluids; physiology `circulation` | Poiseuille–Kirchhoff (already in physis); 0D block DAEs (RC, RCR, elastance heart, Mynard valves, Starling resistor); 1D balance law with a (K, m, n) tube-law family, well-balanced MUSCL finite volumes, junctions, gravity and external pressure; outlet objects returning P(Q) and dP/dQ | RC/RCR closed forms; Womersley; exact Riemann problems with collapse (Toro–Siviglia, Spiller, Spilimbergo); the graph-Laplacian steady state | Boileau cases via openBF references (Apache-2.0); AN37 in vitro; ADAN56; PWDB (ODC-PDDL); ADAVN network; Müller–Toro, Keijsers and Ferlito waveforms; Pfaller's 72 models |
| 5 | physis fluids; physiology rheology | D3Q19 TRT/MRT lattice Boltzmann with BFL walls from a signed distance field; Windkessel outlets through density; μ(γ̇, H) strategies; streaming WSS-index accumulators | Poiseuille at second order; Womersley at α 4–16; Dean flow | Purdue stenosis DNS; FDA nozzle at Re 500–6,500; SBC2012 and 2015 challenge cases; AnXplore |
| 6 | physis; physiology | ADR by SUPG or finite volumes; Strang splitting with Rosenbrock/BDF reactions; residence time; O2 saturation laws; SBML import of coagulation; tracer kinetics | Krogh–Erlang; Graetz–Lévêque; exponential layers at Pe up to 10⁶; V = F·t̄; Tofts with Parker's AIF; mass-action conservation | BioModels BIOMD0000000335 and others (CC0); Ranc 2023 CAT data; Bouchnita 2020 venous model; Celaya-Alcala distribution; OSIPI DCE; Kudo phantom |
| 7 | medical-imaging; tyche | 4D flow operator H (voxel and bin averaging, VENC wrap, noise, offsets); ROUKF on Gaussian conditioning; SMC Windkessel calibration on 0D surrogates; geometry as an inferred parameter; context-of-use records | twin experiments: H applied to the infrastructure's own simulations | Northwestern aorta atlas (CC BY-NC-SA, cache only); Puiseux in vitro protocol; Vascular Model Repository (terms unchecked) |

Fabrica follows from stage 3 onward, drawing each stage's fields through the visualization layer's scene model. It shows the vascular graph, the 1D waves on it, the voxel flow, and then the transported species.

### 8.3 Placement between the foundations and the domains

The placement follows the repository's rule that general methods go to foundations and named objects to domains.

**physis.** It takes the general mechanics:

- the hyperbolic balance law with an abstract tube law and its Riemann solvers;
- junction and outlet interfaces;
- lattice Boltzmann;
- advection–diffusion–reaction with stabilisation;
- the generalised-Newtonian interface;
- WSS functionals.

**physiology `circulation`.** It takes what names blood and its vessels:

- the venous and arterial tube-law parameters;
- the Carreau–Yasuda, Quemada and Pries parameter sets;
- the elastance heart, venous valves, Starling resistors, the muscle pump and posture;
- the oxygen saturation laws;
- the coagulation networks.

**anatomy `vessels`.** It owns the vascular graph as an anatomical object, with named vessels and valve sites.

**medical-imaging.** It owns modality physics and the operators H.

**eikon.** It owns every image algorithm that names no organ, with morphe supplying homology.

**Order of the domains.** Physiology, placed at height 2 by ADR 0061, will need anatomy's vessels and medical-imaging's operators. Its declared dependencies are therefore a placement decision to record before work begins.

### 8.4 Licences decide what is borrowed, linked or reimplemented

| Rule | Code | Data |
|---|---|---|
| Borrow or port, keeping notices | ITK, OpenVDB, nnU-Net, MONAI, TotalSegmentator code, openBF (Apache-2.0); VTK, Geogram, svMultiPhysics, svZeroDSolver, MFEM (BSD-3); SimVascular, svOneDSolver (MIT-style); VesselGraph, clDice, 4DFlowNet, WSSNet, neuraloperator (MIT); VMTK (BSD) without its bundled TetGen | small values in fixtures, with attribution: PWDB (ODC-PDDL); BioModels (CC0); TotalSegmentator CT and MRI, AVT, HRA objects, BodyParts3D since 2025-02-25 (CC BY 4.0) |
| Link dynamically, never copy | HemeLB (LGPL-3.0); Netgen (LGPL-2.1); MMG (LGPL-3.0+); DOLFINx, pyadjoint (LGPL-3.0); Verdandi, DeepXDE (LGPL-2.1); fTetWild (MPL-2.0) | share-alike, cached and never committed: MSD (CC BY-SA 4.0); IXI (CC BY-SA 3.0); HVA-CT (CC BY-SA 4.0) |
| External process only, for cross-verification | TetGen ≥ 1.5 (AGPL-3.0); Gmsh (GPL-2.0+); CGAL meshing, mesh processing and skeletons (GPL-3.0+); OpenLB, waLBerla, Palabos, OpenFOAM, clotFoam, Oasis, VaMPy, ClearMap, VesselVio (GPL/AGPL) | non-commercial or gated, fetched by the user: 3D-IRCADb (CC BY-NC-ND 4.0); Northwestern 4D flow atlas (CC BY-NC-SA 3.0 US); AneuX (CC BY-NC 4.0); TopCoW (commercial use by permission); ImageCAS, ASOCA, Parse2022, KiPA22, VESSEL12 (agreements); COSTA (academic only) |
| Reimplement from the paper; compare only | CRIMSON (AGPL-3.0 + Commons Clause); vesselFM (GPL-3.0, non-commercial weights); FluidX3D (non-commercial); nnInteractive weights (CC BY-NC-SA 4.0); Jerman's MATLAB, PINN-wss, STARFiSh (no licence) | — |

### 8.5 What waits

Four things should wait.

- **Full ALE fluid–structure interaction.** Compliance mainly moves oscillatory indices (§5.4), which the clinically proven pressure outputs do not need.
- **Training foundation models.** vesselFM shows the cost, and its licences show the constraint. Topology-aware fine-tuning of nnU-Net with the permissive losses is the practical route.
- **Learned surrogates and PINNs.** They need the verified solvers to make their training data and to bound their errors (§4.5, §7.3).
- **Patient-scale venous 3D assimilation.** It has no open venous 4D flow data to assimilate.

## Conclusion

**The venous problem inverts the arterial one.** Arteries are data-rich and mechanics-light: thousands of labelled scans, multi-team benchmarks and clinical trials, governed by elastic walls and Murray-like splits. Veins are mechanics-rich and data-poor: collapse, valves, gravity and muscle pumping are well formulated, but labelled venous images are scattered, venous 4D flow data are absent, and no venous reference benchmark exists. An infrastructure that holds its code to closed forms therefore has a distinctive contribution to make for veins. One tube-law family is exactly verified through collapse by Riemann solutions. A public-domain venous verification suite, published under the same terms as PWDB, would fill a gap the literature names but has not closed. **Licensing and science point the same way for 3D:** the copyleft of the meshers favours a voxel-native lattice Boltzmann solver, and the evidence favours it too.

**One abstraction serves both validation and assimilation.** Because geometry dominates the error, the reconstructed vessel should be treated as an uncertain parameter, not a fixed input, and the comparison "in the space of the measurement" (Puiseux) and the assimilation of that measurement (ROUKF, SMC, Bayesian inverse Navier–Stokes) share one forward operator H, which removes much duplicated machinery. **The venous valve pocket is a research target that only an integrated programme can pose:** anatomy places the valve, flow creates stasis and residence time, oxygen transport desaturates the pocket, and the coagulation network responds to the hypoxic endothelium. The literature describes each link but has not coupled them; a library holding image geometry, reduced and 3D flow, oxygen and coagulation under one set of conventions is what that coupling needs.
