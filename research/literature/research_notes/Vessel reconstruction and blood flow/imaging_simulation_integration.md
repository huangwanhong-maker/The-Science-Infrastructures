# Imaging–Simulation Integration for Vascular Flow: 4D Flow MRI, Data Assimilation, Physics-Informed Learning and Cardiovascular Digital Twins

*How these notes were checked (2026-11-04).* Every bibliographic entry below that carries a DOI was checked against Crossref's metadata on 2026-11-04: authors, year, title, venue, volume, issue, pages or article number, and DOI. A link to `doi.org` is given as the source. Where a finding describes what a paper *reports*, it comes from the abstract Crossref returned, or from a page that was fetched; claims that rest on a title alone are worded as such. Licences were read from each repository's LICENSE file on GitHub (raw files, 2026-11-04) unless stated otherwise. Items not verified are marked **[unverified]**. "Classical" means up to 2022; "Latest" means 2023–2026.

## 1. 4D Flow MRI and Phase-Contrast Velocimetry: Acquisition, Errors, Derived Quantities, Validation and Open Data

### Takeaway
There are two consensus statements, both verified: Dyverfeldt et al. 2015 (JCMR 17:72) and the 2023 update by Bissell et al. (JCMR 25:40). Each derived quantity has its own error behaviour. Image-based wall shear stress (WSS) is biased low and depends on resolution. Relative pressure needs a PPE, work–energy (WERP) or virtual work–energy (vWERP) estimator and is sensitive to noise. A comparison with CFD agrees only after the MRI artefacts are corrected and the CFD is degraded to the MRI's resolution. The simulation-to-image interface should therefore be an explicit forward measurement operator, which an MR simulator such as CMRsim can supply. Open 4D flow data are scarce and often under non-commercial terms: the one licence confirmed is CC BY-NC-SA 3.0 US.

### Cited Findings

#### Consensus statements and reviews
- **The 2015 consensus statement.** Dyverfeldt P, Bissell M, Barker AJ, Bolger AF, Carlhäll CJ, Ebbers T, Francois CJ, Frydrychowicz A, Geiger J, Giese D, Hope MD, Kilner PJ, et al. (17 authors). "4D flow cardiovascular magnetic resonance consensus statement." *J Cardiovasc Magn Reson* 2015;17(1):72. doi:10.1186/s12968-015-0174-5 — [DOI](https://doi.org/10.1186/s12968-015-0174-5)
- **The 2023 update, confirmed as asked.** Bissell MM, Raimondi F, Ait Ali L, Allen BD, Barker AJ, Bolger A, Burris N, Carhäll CJ, Collins JD, Ebbers T, Francois CJ, Frydrychowicz A, et al. (41 authors). "4D Flow cardiovascular magnetic resonance consensus statement: 2023 update." *J Cardiovasc Magn Reson* 2023;25(1):40. doi:10.1186/s12968-023-00942-z — [DOI](https://doi.org/10.1186/s12968-023-00942-z). It is open access as PMC10357639. Its content below was read from the full text. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10357639/fullTextXML)
  - **The update's focus.** It builds on the 2015 statement, with recommended acquisition parameters for clinical use drawn from centres using 4D flow clinically, a processing workflow, quality assurance, validation, and reporting standards.
  - **Table 1, acquired spatial resolution** ("aim for" / "do not exceed"):

    | Target | Aim for | Do not exceed |
    |---|---|---|
    | Adults, whole heart | 2.5 mm³ | 3 mm³ |
    | Adults, vessels | 2 mm³ | 2.5 mm³ |
    | Paediatric, whole heart | 2 mm³ | 2.5 mm³ |
    | Paediatric, vessels | 1.5 mm³ | 2 mm³ |
    | Neonates | 0.75–1 mm³ | 1.5 mm³ |

  - **Table 1, other parameters.**
    - Acquired temporal resolution: 30 ms, not exceeding 50 ms.
    - VENC: the maximum expected velocity; "10% higher than maximum expected velocity, do not exceed 25%."
    - ECG gating: retrospective.
    - Clinical acquisitions: kept to 5–10 min.
  - **VENC guidance.** VENC is to be as low as possible for velocity-to-noise ratio (VNR) while avoiding aliasing. Without a prior 2D PC scan and with stenosis suspected, start at 250 cm/s. Otherwise:
    - aorta and pulmonary artery: 150 cm/s;
    - false lumen of a dissection: 50–150 cm/s;
    - venous and Fontan flow: 50–80 cm/s;
    - intracardiac: 100–150 cm/s.
  - **The processing workflow.**
    1. Background phase offset correction (eddy currents, concomitant gradients).
    2. Anti-aliasing (phase unwrapping) if required. Regions with "incorrigible aliasing" are to be excluded.
    3. Segmentation.
    4. Visualisation.
    5. Quantification, with an internal consistency check.
  - **Multi-VENC.** Multi-VENC acquisition and reconstruction "can generate 4D Flow CMR data with the favorable VNR of a low-VENC acquisition but without velocity aliasing." Physics-informed neural networks are mentioned as early work for denoising, resolution and automatic unwrapping.
  - **Table 2: three ways to validate.**
    - In vivo: lacks a gold standard.
    - Phantoms: realistic phantoms are hard to make.
    - Simulations: "very controllable, different sources of error can be separated, underlying numerical velocity data serves as ground truth," but "fidelity uncertain, computational cost."

    The statement adds that "the generation of accurate numerical velocity data using computational fluid dynamics is a separate field of research." It also notes that reference data are often unavailable in vivo for WSS, turbulence, kinetic energy and vorticity.
  - **Table 3: reporting.** Acquired *and* reconstructed spatial and temporal resolution must both be reported. The example given: acquired 2.5 × 2.5 × 2.5 mm³, reconstructed 1.0 × 1.0 × 1.0 mm³, effective temporal resolution 35 ms ≈ TR × flow encodings × segments. Every correction applied must also be reported: eddy currents, gradient non-uniformity, intravoxel dephasing, concomitant gradients, aliasing, filtering.
- **A standard review.** Markl M, Frydrychowicz A, Kozerke S, Hope M, Wieben O. "4D flow MRI." *J Magn Reson Imaging* 2012;36(5):1015–1036. doi:10.1002/jmri.23632 — [DOI](https://doi.org/10.1002/jmri.23632). *Caveat:* Crossref also lists a second DOI for the same title without pages (10.1002/jmri.23556); cite 10.1002/jmri.23632.
- **Basic facts of phase-contrast flow measurement.** Lotz J, Meier C, Leppert A, Galanski M. "Cardiovascular flow measurement with phase-contrast MR imaging: basic facts and implementation." *RadioGraphics* 2002;22(3):651–671. doi:10.1148/radiographics.22.3.g02ma11651 — [DOI](https://doi.org/10.1148/radiographics.22.3.g02ma11651)

#### Sources of acquisition error (classical)
- **Noise and background phase.** Walker PG, Cranney GB, Scheidegger MB, et al. "Semiautomated method for noise reduction and background phase error correction in MR phase velocity data." *J Magn Reson Imaging* 1993;3(3):521–530. doi:10.1002/jmri.1880030315 — [DOI](https://doi.org/10.1002/jmri.1880030315)
- **Partial volume.** Tang C, Blatter DD, Parker DL. "Accuracy of phase-contrast flow measurements in the presence of partial-volume effects." *J Magn Reson Imaging* 1993;3(2):377–385. doi:10.1002/jmri.1880030213 — [DOI](https://doi.org/10.1002/jmri.1880030213)
- **Concomitant (Maxwell) gradient terms.** Bernstein MA, Zhou XJ, Polzin JA, King KF, Ganin A, Pelc NJ, Glover GH. "Concomitant gradient terms in phase contrast MR: analysis and correction." *Magn Reson Med* 1998;39(2):300–308. doi:10.1002/mrm.1910390218 — [DOI](https://doi.org/10.1002/mrm.1910390218)
- **Background phase offsets across sites and vendors.** Gatehouse PD, Rolf MP, Graves MJ, Hofman MB, Totman J, et al. "Flow measurement by cardiovascular magnetic resonance: a multi-centre multi-vendor study of background phase offset errors that can compromise the accuracy of derived regurgitant or shunt flow measurements." *J Cardiovasc Magn Reson* 2010;12(1):5. doi:10.1186/1532-429X-12-5 — [DOI](https://doi.org/10.1186/1532-429x-12-5)
- **Aliasing when velocity exceeds VENC: phase unwrapping.** Loecher M, Schrauben E, Johnson KM, Wieben O. "Phase unwrapping in 4D MR flow with a 4D single-step Laplacian algorithm." *J Magn Reson Imaging* 2016;43(4):833–842. doi:10.1002/jmri.25045 — [DOI](https://doi.org/10.1002/jmri.25045)
- **Aliasing and velocity-to-noise: dual-VENC.** Schnell S, Ansari SA, Wu C, Garcia J, Murphy IG, Rahman OA, Rahsepar AA, Aristova M, Collins JD, Carr JC, Markl M. "Accelerated dual-venc 4D flow MRI for neurovascular applications." *J Magn Reson Imaging* 2017;46(1):102–114. doi:10.1002/jmri.25595 — [DOI](https://doi.org/10.1002/jmri.25595)
- **Moving valve planes: retrospective valve tracking.** Westenberg JJM, Roes SD, Ajmone Marsan N, et al. "Mitral valve and tricuspid valve blood flow: accurate quantification with 3D velocity-encoded MR imaging with retrospective valve tracking." *Radiology* 2008;249(3):792–800. doi:10.1148/radiol.2492080146 — [DOI](https://doi.org/10.1148/radiol.2492080146)
- **Acquisition parameters in an aortic phantom.** Montalba C, Urbina J, Sotelo J, Andia ME, Tejos C, Irarrazaval P, Hurtado DE, Valverde I, Uribe S. "Variability of 4D flow parameters when subjected to changes in MRI acquisition parameters using a realistic thoracic aortic phantom." *Magn Reson Med* 2018;79(4):1882–1892. doi:10.1002/mrm.26834 — [DOI](https://doi.org/10.1002/mrm.26834)
- **The noise in a turbulent acquisition, simulated.** Petersson S, Dyverfeldt P, Gårdhagen R, Karlsson M, Ebbers T. "Simulation of phase contrast MRI of turbulent flow." *Magn Reson Med* 2010;64(4):1039–1046. doi:10.1002/mrm.22494 — [DOI](https://doi.org/10.1002/mrm.22494)

#### Wall shear stress from 4D flow
- Stalder AF, Russe MF, Frydrychowicz A, Bock J, Hennig J, Markl M. "Quantitative 2D and 3D phase contrast MRI: optimized analysis of blood flow and vessel wall parameters." *Magn Reson Med* 2008;60(5):1218–1231. doi:10.1002/mrm.21778 — [DOI](https://doi.org/10.1002/mrm.21778)
- Petersson S, Dyverfeldt P, Ebbers T. "Assessment of the accuracy of MRI wall shear stress estimation using numerical simulations." *J Magn Reson Imaging* 2012;36(1):128–138. doi:10.1002/jmri.23610 — [DOI](https://doi.org/10.1002/jmri.23610)
- van Ooij P, Potters WV, Guédon A, et al. "Wall shear stress estimated with phase contrast MRI in an in vitro and in vivo intracranial aneurysm." *J Magn Reson Imaging* 2013;38(4):876–884. doi:10.1002/jmri.24051 — [DOI](https://doi.org/10.1002/jmri.24051)
- Potters WV, van Ooij P, Marquering H, vanBavel E, Nederveen AJ. "Volumetric arterial wall shear stress calculation based on cine phase contrast MRI." *J Magn Reson Imaging* 2015;41(2):505–516. doi:10.1002/jmri.24560 — [DOI](https://doi.org/10.1002/jmri.24560)
- Cibis M, Potters WV, Gijsen FJ, Marquering H, van Ooij P, vanBavel E, Wentzel JJ, Nederveen AJ. "The effect of spatial and temporal resolution of cine phase contrast MRI on wall shear stress and oscillatory shear index assessment." *PLoS One* 2016;11(9):e0163316. doi:10.1371/journal.pone.0163316 — [DOI](https://doi.org/10.1371/journal.pone.0163316)
- **The bias, stated by WSSNet's authors.** Image-based WSS from 4D flow MRI "underestimates true WSS values, and the accuracy is dependent on spatial resolution." Their learned estimator, trained on CFD of 37 patient-specific aortas, reached the following on validation and test cases:
  - against CFD: mean absolute error 0.55 ± 0.60 Pa, relative error 4.34 ± 4.14%, Pearson r 0.92;
  - on noisy synthetic 4D flow at 2.4 mm: 0.99 ± 0.91 Pa, 7.13 ± 6.27%, r 0.79.

  Ferdian E, Dubowitz DJ, Mauger CA, Wang A, Young AA. "WSSNet: aortic wall shear stress estimation using deep learning on 4D flow MRI." *Front Cardiovasc Med* 2022;8:769927. doi:10.3389/fcvm.2021.769927 — [DOI](https://doi.org/10.3389/fcvm.2021.769927)

#### Pressure from velocity
- **Early relative-pressure mapping.** Tyszka JM, Laidlaw DH, Asa JW, Silverman JM. "Three-dimensional, time-resolved (4D) relative pressure mapping using magnetic resonance imaging." *J Magn Reson Imaging* 2000;12(2):321–329. doi:10.1002/1522-2586(200008)12:2<321::AID-JMRI15>3.0.CO;2-2 — [DOI](https://doi.org/10.1002/1522-2586(200008)12:2%3C321::AID-JMRI15%3E3.0.CO;2-2)
- **Early relative-pressure estimation.** Ebbers T, Wigström L, Bolger AF, Engvall J, Karlsson M. "Estimation of relative cardiovascular pressures using time-resolved three-dimensional phase contrast MRI." *Magn Reson Med* 2001;45(5):872–879. doi:10.1002/mrm.1116 — [DOI](https://doi.org/10.1002/mrm.1116)
- **Nasiraei-Moghaddam et al.: the effects of noise, resolution and flow rate on pressure in stenoses.** Their abstract reports that pressures were computed in an axisymmetric geometry. Phantoms with 50%, 75% and 90% area stenoses were both imaged and simulated by CFD, and Gaussian noise was added to the CFD velocities. Nasiraei-Moghaddam A, Behrens G, Fatouraee N, Agarwal R, Choi ET, Amini AA. "Factors affecting the accuracy of pressure measurements in vascular stenoses from phase-contrast MRI." *Magn Reson Med* 2004;52(2):300–309. doi:10.1002/mrm.20152 — [DOI](https://doi.org/10.1002/mrm.20152)
- **Krittian et al.: pressure computed directly by finite elements (the user's citation, verified).** Krittian SBS, Lamata P, Michler C, Nordsletten DA, Bock J, Bradley CP, Pitcher A, Kilner PJ, Markl M, Smith NP. "A finite-element approach to the direct computation of relative cardiovascular pressure from time-resolved MR velocity data." *Med Image Anal* 2012;16(5):1029–1037. doi:10.1016/j.media.2012.04.003 — [DOI](https://doi.org/10.1016/j.media.2012.04.003)
- **The work–energy relative pressure (WERP).** Donati F, Figueroa CA, Smith NP, Lamata P, Nordsletten DA. "Non-invasive pressure difference estimation from PC-MRI using the work-energy equation." *Med Image Anal* 2015;26(1):159–172. doi:10.1016/j.media.2015.08.012 — [DOI](https://doi.org/10.1016/j.media.2015.08.012)
- **A review of the estimators, with new approaches.** Bertoglio C, Nuñez R, Galarce F, Nordsletten D, Osses A. "Relative pressure estimation from velocity measurements in blood flows: state-of-the-art and new approaches." *Int J Numer Method Biomed Eng* 2018;34(2):e2925. doi:10.1002/cnm.2925 — [DOI](https://doi.org/10.1002/cnm.2925)
- **The virtual work–energy relative pressure (vWERP).**
  - Marlevi D, et al. "Estimation of cardiovascular relative pressure using virtual work-energy." *Sci Rep* 2019;9(1):1375. doi:10.1038/s41598-018-37714-0 — [DOI](https://doi.org/10.1038/s41598-018-37714-0)
  - Marlevi D, Balmus M, Hessenthaler A, et al. "Non-invasive estimation of relative pressure for intracardiac flows using virtual work-energy." *Med Image Anal* 2021;68:101948. doi:10.1016/j.media.2020.101948 — [DOI](https://doi.org/10.1016/j.media.2020.101948)
- **Validation of 4D-flow-based relative pressure maps in the aorta.** Nolte D, Urbina J, Sotelo J, Sok L, Montalba C, Valverde I, Osses A, Uribe S, Bertoglio C. "Validation of 4D Flow based relative pressure maps in aortic flows." *Med Image Anal* 2021;74:102195. doi:10.1016/j.media.2021.102195 — [DOI](https://doi.org/10.1016/j.media.2021.102195)
- **A quality check by Navier–Stokes compatibility.** Garay J, Mella H, Sotelo J, Cárcamo C, Uribe S, Bertoglio C, Mura J. "Assessment of 4D flow MRI's quality by verifying its Navier–Stokes compatibility." *Int J Numer Method Biomed Eng* 2022;38(6):e3603. doi:10.1002/cnm.3603 — [DOI](https://doi.org/10.1002/cnm.3603)

#### Validation against CFD, and simulating the MRI
- **Cibis et al., as the user cited.** Two verified papers match "Cibis et al.":
  - the 2016 *PLoS One* paper on the effect of resolution on WSS and OSI (above);
  - Cibis M, Jarvis K, Markl M, et al. "The effect of resolution on viscous dissipation measured with 4D flow MRI in patients with Fontan circulation: evaluation using computational fluid dynamics." *J Biomech* 2015;48(12):2984–2989. doi:10.1016/j.jbiomech.2015.07.039 — [DOI](https://doi.org/10.1016/j.jbiomech.2015.07.039)
- **CFD against 4D PC-MRI in the Circle of Willis.** Berg P, Stucht D, Janiga G, Beuing O, Speck O, Thévenin D. "Cerebral blood flow in a healthy Circle of Willis and two intracranial aneurysms: computational fluid dynamics versus four-dimensional phase-contrast magnetic resonance imaging." *J Biomech Eng* 2014;136(4):041003. doi:10.1115/1.4026108 — [DOI](https://doi.org/10.1115/1.4026108)
- **Puiseux et al. 2019, in vitro: the comparison works only in the space of the measurement.**
  - Between raw 4D flow and CFD over the whole domain, r² was only 0.63.
  - It rose to r² = 0.97 when the MR artefacts were removed by post-processing and the CFD was down-sampled to the MRI's spatial and temporal resolution.

  Puiseux T, Sewonu A, Meyrignac O, Rousseau H, Nicoud F, Mendez S, Moreno R. "Reconciling PC-MRI and CFD: an in-vitro study." *NMR Biomed* 2019;32(5):e4063. doi:10.1002/nbm.4063 — [DOI](https://doi.org/10.1002/nbm.4063)
- **Simulating the 4D flow acquisition itself from a CFD field.** Puiseux T, Sewonu A, Moreno R, Mendez S, Nicoud F. "Numerical simulation of time-resolved 3D phase-contrast magnetic resonance imaging." *PLoS One* 2021;16(3):e0248816. doi:10.1371/journal.pone.0248816 — [DOI](https://doi.org/10.1371/journal.pone.0248816)
- **Inflow boundary conditions taken from 4D flow MRI.**
  - Youssefi P, Gomez A, Arthurs C, Sharma R, Jahangiri M, Figueroa CA. "Impact of patient-specific inflow velocity profile on hemodynamics of the thoracic aorta." *J Biomech Eng* 2018;140(1):011002. doi:10.1115/1.4037857 — [DOI](https://doi.org/10.1115/1.4037857)
  - Pirola S, Jarral OA, O'Regan DP, et al. "Computational study of aortic hemodynamics for patients with an abnormal aortic valve: the importance of secondary flow at the ascending aorta inlet." *APL Bioeng* 2018;2(2):026101. doi:10.1063/1.5011960 — [DOI](https://doi.org/10.1063/1.5011960)

#### Latest (2023–2026)
- **The 2023 consensus update** (above).
- **CMRsim, an open-source Python package for cardiovascular MR simulation.** It takes dynamic digital phantoms with motion and flow and solves the Bloch equations numerically or analytically on GPUs, through TensorFlow. Among its demonstrations is phase-contrast imaging of turbulent flow downstream of a stenosis. Weine J, McGrath C, Dirix P, Buoso S, Kozerke S. "CMRsim–A python package for cardiovascular MR simulations incorporating complex motion and flow." *Magn Reson Med* 2024;91(6):2621–2637. doi:10.1002/mrm.30010 — [DOI](https://doi.org/10.1002/mrm.30010)
- **Dirix et al. 2022: synthetic patient-specific 4D flow MRI paired with ground-truth CFD, for training inference methods.**
  - The scan time is held fixed, so changing the resolution or the time averaging rescales the SNR.
  - Velocity and turbulent kinetic energy are reconstructed by a Bayesian method.

  This is the basis of ETH's synthetic cohorts. Dirix P, Buoso S, Peper ES, Kozerke S. "Synthesis of patient-specific multipoint 4D flow MRI data of turbulent aortic flow downstream of stenotic valves." *Sci Rep* 2022;12(1):16004. doi:10.1038/s41598-022-20121-x — [DOI](https://doi.org/10.1038/s41598-022-20121-x)
- **Inlet velocity profiles generated from data, for the boundary conditions of an aneurysmal ascending aorta.** Saitta S, Maga L, Armour C, Votta E, et al. "Data-driven generation of 4D velocity profiles in the aneurysmal ascending aorta." *Comput Methods Programs Biomed* 2023;233:107468. doi:10.1016/j.cmpb.2023.107468 — [DOI](https://doi.org/10.1016/j.cmpb.2023.107468)

#### CT for geometry, and perfusion imaging as data on transport
- **CT angiography to FFR by CFD, the basis of HeartFlow.** Taylor CA, Fonte TA, Min JK. "Computational fluid dynamics applied to cardiac computed tomography for noninvasive quantification of fractional flow reserve: scientific basis." *J Am Coll Cardiol* 2013;61(22):2233–2241. doi:10.1016/j.jacc.2012.11.083 — [DOI](https://doi.org/10.1016/j.jacc.2012.11.083)
- **The uncertainty from geometry.**
  - Sankaran S, Kim HJ, Choi G, Taylor CA. "Uncertainty quantification in coronary blood flow simulations: impact of geometry, boundary conditions and blood viscosity." *J Biomech* 2016;49(12):2540–2547. doi:10.1016/j.jbiomech.2016.01.002 — [DOI](https://doi.org/10.1016/j.jbiomech.2016.01.002)
  - Maher GD, Fleeter CM, Schiavazzi DE, Marsden AL. "Geometric uncertainty in patient-specific cardiovascular modeling with convolutional dropout networks." *Comput Methods Appl Mech Eng* 2021;386:114038. doi:10.1016/j.cma.2021.114038 — [DOI](https://doi.org/10.1016/j.cma.2021.114038)
- **Tracer kinetics, the observation models for DCE-MRI.**
  - Tofts PS, Brix G, Buckley DL, et al. "Estimating kinetic parameters from dynamic contrast-enhanced T1-weighted MRI of a diffusable tracer: standardized quantities and symbols." *J Magn Reson Imaging* 1999;10(3):223–232. doi:10.1002/(SICI)1522-2586(199909)10:3<223::AID-JMRI2>3.0.CO;2-S — [DOI](https://doi.org/10.1002/(SICI)1522-2586(199909)10:3%3C223::AID-JMRI2%3E3.0.CO;2-S)
  - Sourbron SP, Buckley DL. "Classic models for dynamic contrast-enhanced MRI." *NMR Biomed* 2013;26(8):1004–1027. doi:10.1002/nbm.2940 — [DOI](https://doi.org/10.1002/nbm.2940)
- **Perfusion as porous flow.**
  - Michler C, Cookson AN, Chabiniok R, et al. "A computationally efficient framework for the simulation of cardiac perfusion using a multi-compartment Darcy porous-media flow model." *Int J Numer Method Biomed Eng* 2013;29(2):217–232. doi:10.1002/cnm.2520 — [DOI](https://doi.org/10.1002/cnm.2520)
  - Hyde ER, Cookson AN, Lee J, et al. "Multi-scale parameterisation of a myocardial perfusion model using whole-organ arterial networks." *Ann Biomed Eng* 2014;42(4):797–811. doi:10.1007/s10439-013-0951-y — [DOI](https://doi.org/10.1007/s10439-013-0951-y)
- **Papamanolis et al. 2021: from CT through PET to the myocardium.** Their abstract reports:
  - vessels were segmented from coronary CTA and extended to the arterioles by space-filling synthetic trees;
  - a 1D coronary model was coupled to a single-compartment Darcy myocardium;
  - in five patients, at rest and in hyperaemia, the simulations "compare overall well" with [¹⁵O]H₂O PET.

  Papamanolis L, Kim HJ, Jaquet C, Sinclair M, Schaap M, Danad I, van Diemen P, Knaapen P, Najman L, Talbot H, Taylor CA, Vignon-Clementel I. "Myocardial perfusion simulation for coronary artery disease: a coupled patient-specific multiscale model." *Ann Biomed Eng* 2021;49(5):1432–1447. doi:10.1007/s10439-020-02681-z — [DOI](https://doi.org/10.1007/s10439-020-02681-z)
- **Józsa et al. 2021: a brain model mapped by ASL.** A three-compartment porous model of the brain's microcirculation, whose mapping to the large arteries comes from vessel-encoded ASL MRI and clustering, for in silico trials in stroke. Józsa TI, Padmos RM, Samuels N, El-Bouri WK, Hoekstra AG, Payne SJ. "A porous circulation model of the human brain for in silico clinical trials in ischaemic stroke." *Interface Focus* 2021;11(1):20190127. doi:10.1098/rsfs.2019.0127 — [DOI](https://doi.org/10.1098/rsfs.2019.0127)
- **CFD combined with ASL.** Schollenberger J, Osborne NH, Hernandez-Garcia L, et al. "A combined computational fluid dynamics and arterial spin labeling MRI modeling strategy to quantify patient-specific cerebral hemodynamics in cerebrovascular occlusive disease." *Front Bioeng Biotechnol* 2021;9:722445. doi:10.3389/fbioe.2021.722445 — [DOI](https://doi.org/10.3389/fbioe.2021.722445)
- **Physics-informed networks for perfusion.** van Herten RLM, Chiribiri A, Breeuwer M, Veta M, Scannell CM. "Physics-informed neural networks for myocardial perfusion MRI quantification." *Med Image Anal* 2022;78:102399. doi:10.1016/j.media.2022.102399 — [DOI](https://doi.org/10.1016/j.media.2022.102399)

#### Open 4D flow datasets and their licences
- **Northwestern's "4D Flow Healthy Aorta Atlases"** (Scott MB, Barker AJ, Markl M), v2.0.0, published 2020-02-26. It is a 37.9 MB archive under **CC BY-NC-SA 3.0 US (non-commercial, share-alike)**, DOI 10.18131/g3-knrk-5191, with the companion code at github.com/alexbarks/atlas-tool. — [Northwestern PRISM record](https://prism.northwestern.edu/records/r29kc-x5336)
- **Records found whose licence could not be confirmed** (search snippets only). — [search results: Zenodo 12515230](https://zenodo.org/records/12515230), [Zenodo 4882572](https://zenodo.org/records/4882572), [Stanford AS4DF](https://purl.stanford.edu/dz488kx6180), [Auckland figshare](https://auckland.figshare.com/articles/dataset/4DFlowNet_-_high_resolution_aortic_CFD_dataset/24424888)
  - Ohio State: free-running 3D cine and 4D flow, about 6.6 GB, Zenodo 12515230 (v2, June 2024);
  - in vitro 4D flow of silicone tubes and an aneurysm phantom, Zenodo 4882572 (2021);
  - Stanford's AS4DF: a compliant 3D-printed thoracic aorta phantom, with 2D PC, cine and 4D flow at three temporal resolutions;
  - the aortic CFD dataset of 4DFlowNet, on the University of Auckland's figshare (the fetch returned HTTP 403);
  - ETH Zurich's synthetic RACLETTE collection (189 virtual subjects for aortic 4D flow and 25 for 5D flow), announced only through a LinkedIn link.
- **A restricted record.** Zenodo 10658225 (a 4D flow pilot in TEVAR patients) is listed as restricted. — [Zenodo](https://zenodo.org/records/10658225)
- **The Vascular Model Repository** holds images, models and simulation results; its licence was not checked. Wilson NM, Ortiz AK, Johnson AB. "The Vascular Model Repository: a public resource of medical imaging data and blood flow simulation results." *J Med Devices* 2013;7(4):040923. doi:10.1115/1.4025983 — [DOI](https://doi.org/10.1115/1.4025983)

### Inferences
- **The forward measurement operator H.** The layer needs an explicit H that maps a simulated velocity field to a 4D flow image. H should apply:
  - averaging over the voxel and the cardiac bin (partial volume, temporal blurring);
  - wrapping at VENC;
  - noise whose variance depends on VENC and SNR;
  - background phase offsets and concomitant-gradient terms.

  Validation should then be done in the space of the image, as H(simulation) against the image, never by upsampling the image (Puiseux 2019). A full Bloch-equation H, as in CMRsim or Puiseux 2021, is the high-fidelity option; voxel-averaging plus noise is the cheap one.
- **H's default parameters can come from Table 1 of the 2023 consensus.** For adult vessels: acquired voxels of 2–2.5 mm, temporal resolution 30–50 ms, VENC 10–25% above the peak velocity. H must use the *acquired* resolution, not the reconstructed one, which the consensus requires every report to state separately; reconstructed data are usually interpolated, e.g. to 1 mm.
- **The statement itself defines simulation-based validation's place.** It lists it as "very controllable" but of "fidelity uncertain." The infrastructure's simulations should therefore carry their own VVUQ record (§4) when used as ground truth for imaging methods.
- **WSS and relative pressure are best treated as derived quantities of interest of the assimilated model**, not as direct observations, since both are biased or noisy when taken straight from the image (WSSNet; Petersson 2012; Nasiraei-Moghaddam 2004). Navier–Stokes compatibility (Garay 2022) can serve as a check of an image's quality before it is assimilated.
- **The licence of the only confirmed aortic 4D flow atlas (CC BY-NC-SA 3.0 US) is not compatible with redistributing it in an MIT repository.** Fetch it into the data cache at use time, as the infrastructure already does for Pangu-Weather's weights, and never commit it.

### Gaps
- **Pelc et al. 1991 could not be verified.** "Phase contrast cine magnetic resonance imaging," *Magn Reson Q* 1991, the usual source of the noise relation σ_v ∝ VENC/SNR, did not appear in Crossref, which returned a different Pelc 1995 paper. **[unverified]**
- **The 2023 update was read for acquisition, processing and validation (Tables 1–3), but not for its sections on individual clinical applications.** It gives no numerical accuracy figures for WSS or relative pressure in the passages extracted; those come from the primary papers above.
- **Licences of the other 4D flow datasets are unknown** (Ohio State, AS4DF, 4DFlowNet's CFD set, RACLETTE). No open in vivo 4D flow dataset under CC BY or CC0 was found.
- **The Vascular Model Repository's terms of use were not checked.**

## 2. Data Assimilation in Cardiovascular Modelling: Kalman Filters, Variational and Adjoint Methods, Windkessel Calibration, Bayesian UQ and Multi-Fidelity Methods

### Takeaway
Two lineages dominate.

- **Sequential: Kalman filtering.** Moireau and Chapelle's reduced-order unscented Kalman filter (ROUKF) works in parameter space, needs no adjoint, and is applied to fluid–structure interaction and Windkessel parameters by Bertoglio, Moireau and Gerbeau (2012). It is implemented in Verdandi (LGPL-2.1) and CRIMSON (AGPL-3.0 with the Commons Clause).
- **Variational: adjoint-based.** D'Elia, Perego and Veneziani (2012) and Funke et al. (2019), the latter on dolfin-adjoint (LGPL-3.0).

**A correction to the brief:** the verified D'Elia–Perego–Veneziani 2012 paper is a *variational* method, not a Kalman filter. Since 2020 the frontier is Bayesian calibration with posteriors, made affordable by surrogates: Richter et al. 2025 calibrate Windkessel parameters by sequential Monte Carlo on optimised 0D surrogates (svZeroDSolver, BSD-3). The same period added multi-fidelity estimators and Bayesian inverse Navier–Stokes problems that infer the vessel boundary jointly with the flow.

### Cited Findings

#### Sequential (Kalman) estimation (classical)
- Moireau P, Chapelle D, Le Tallec P. "Joint state and parameter estimation for distributed mechanical systems." *Comput Methods Appl Mech Eng* 2008;197(6–8):659–677. doi:10.1016/j.cma.2007.08.021 — [DOI](https://doi.org/10.1016/j.cma.2007.08.021)
- Moireau P, Chapelle D. "Reduced-order unscented Kalman filtering with application to parameter identification in large-dimensional systems." *ESAIM Control Optim Calc Var* 2011;17(2):380–405. doi:10.1051/cocv/2010006 — [DOI](https://doi.org/10.1051/cocv/2010006)
- **Bertoglio, Moireau and Gerbeau 2012: ROUKF applied to haemodynamics.** The filter is restricted to the parameter space, "does not require any adjoint or tangent problems," and "can easily be run in parallel."
  - It estimates the wall's stiffness from the wall's displacement in an idealised 3D abdominal aortic aneurysm.
  - It gives preliminary estimates of a proximal Windkessel resistance.

  Bertoglio C, Moireau P, Gerbeau JF. "Sequential parameter estimation for fluid–structure problems: application to hemodynamics." *Int J Numer Method Biomed Eng* 2012;28(4):434–455. doi:10.1002/cnm.1476 — [DOI](https://doi.org/10.1002/cnm.1476)
- **The library behind ROUKF.** Chapelle D, Fragu M, Mallet V, Moireau P. "Fundamental principles of data assimilation underlying the Verdandi library: applications to biophysical model personalization within euHeart." *Med Biol Eng Comput* 2013;51(11):1221–1233. doi:10.1007/s11517-012-0969-6 — [DOI](https://doi.org/10.1007/s11517-012-0969-6)
- **From clinical measurements to parameter estimates in multi-scale CFD.** Pant S, Fabrèges B, Gerbeau JF, Vignon-Clementel IE. "A methodological paradigm for patient-specific multi-scale CFD simulations: from clinical measurements to parameter estimates for individual analysis." *Int J Numer Method Biomed Eng* 2014;30(12):1614–1648. doi:10.1002/cnm.2692 — [DOI](https://doi.org/10.1002/cnm.2692)
- **Caiazzo et al. 2017: ROUKF on real data.** The filter was run on 1D blood-flow models of an in vitro arterial network, using *real* experimental data. It estimated terminal resistances and wall parameters from at most one pressure or flow measurement per vessel. The study names the filter's limitations in clinically realistic setups. Caiazzo A, Caforio F, Montecinos G, Müller LO, Blanco PJ, Toro EF. "Assessment of reduced-order unscented Kalman filter for parameter identification in 1-dimensional blood flow models using experimental data." *Int J Numer Method Biomed Eng* 2017;33(8):e2843. doi:10.1002/cnm.2843 — [DOI](https://doi.org/10.1002/cnm.2843)
- Lal R, Mohammadi B, Nicoud F. "Data assimilation for identification of cardiovascular network characteristics." *Int J Numer Method Biomed Eng* 2017;33(5):e2824. doi:10.1002/cnm.2824 — [DOI](https://doi.org/10.1002/cnm.2824)
- **Arthurs et al. 2020: the CRIMSON framework.** Its parts:
  - ROUKF for the wall's material and simple lumped-parameter-network (LPN) boundary parameters;
  - a constrained least-squares augmentation (ROUKF-CLS) for complex LPNs;
  - a "Netlist" implementation of the networks;
  - a demonstration on a healthy volunteer's non-invasive anatomy, flow and pressure;
  - implementation in CRIMSON.

  Arthurs CJ, Xiao N, Moireau P, Schaeffter T, Figueroa CA. "A flexible framework for sequential estimation of model parameters in computational hemodynamics." *Adv Model Simul Eng Sci* 2020;7(1):48. doi:10.1186/s40323-020-00186-x — [DOI](https://doi.org/10.1186/s40323-020-00186-x)
- **Assimilation in a subject-specific Circle of Willis.** Gaidzik F, Pathiraja S, Saalfeld S, Stucht D, Speck O, Thévenin D, Janiga G. "Hemodynamic data assimilation in a subject-specific circle of Willis geometry." *Clin Neuroradiol* 2021;31(3):643–651. doi:10.1007/s00062-020-00959-2 — [DOI](https://doi.org/10.1007/s00062-020-00959-2)

#### Variational and adjoint methods, and boundary control (classical)
- **D'Elia et al. 2012: variational, not Kalman.** The brief attributes Kalman filtering to D'Elia, Perego and Veneziani 2012; the verified 2012 papers of this group are variational.
  - D'Elia M, Perego M, Veneziani A. "A variational data assimilation procedure for the incompressible Navier-Stokes equations in hemodynamics." *J Sci Comput* 2012;52(2):340–359. doi:10.1007/s10915-011-9547-6 — [DOI](https://doi.org/10.1007/s10915-011-9547-6)
  - D'Elia M, Mirabella L, Passerini T, Perego M, Piccinelli M, Vergara C, Veneziani A. "Applications of variational data assimilation in computational hemodynamics." In: *Modeling of Physiological Flows* (MS&A series). Springer Milan; 2012:363–394. doi:10.1007/978-88-470-1935-5_12 — [DOI](https://doi.org/10.1007/978-88-470-1935-5_12)
- **The wall's compliance by inverse FSI.** Perego M, Veneziani A, Vergara C. "A variational approach for estimating the compliance of the cardiovascular tissue: an inverse fluid-structure interaction problem." *SIAM J Sci Comput* 2011;33(3):1181–1211. doi:10.1137/100808277 — [DOI](https://doi.org/10.1137/100808277)
- **Tuning to physiological measurements.** Spilker RL, Taylor CA. "Tuning multidomain hemodynamic simulations to match physiological measurements." *Ann Biomed Eng* 2010;38(8):2635–2648. doi:10.1007/s10439-010-0011-9 — [DOI](https://doi.org/10.1007/s10439-010-0011-9)
- **Windkessel parameters by adjoint.** Ismail M, Wall WA, Gee MW. "Adjoint-based inverse analysis of windkessel parameters for patient-specific vascular models." *J Comput Phys* 2013;244:113–130. doi:10.1016/j.jcp.2012.10.028 — [DOI](https://doi.org/10.1016/j.jcp.2012.10.028)
- **CFD regularised by 3D PC-MRI.** Rispoli VC, Nielsen JF, Nayak KS, Carvalho JLA. "Computational fluid dynamics simulations of blood flow regularized by 3D phase contrast MRI." *Biomed Eng Online* 2015;14(1):110. doi:10.1186/s12938-015-0104-7 — [DOI](https://doi.org/10.1186/s12938-015-0104-7)
- **Velocity tracking.** Tiago J, Guerra T, Sequeira A. "A velocity tracking approach for the data assimilation problem in blood flow simulations." *Int J Numer Method Biomed Eng* 2017;33(10):e2856. doi:10.1002/cnm.2856 — [DOI](https://doi.org/10.1002/cnm.2856)
- **Boundary control.** Koltukluoğlu TS, Blanco PJ. "Boundary control in computational haemodynamics." *J Fluid Mech* 2018;847:329–364. doi:10.1017/jfm.2018.329 — [DOI](https://doi.org/10.1017/jfm.2018.329)
- **Funke et al. 2019 (the user's citation, verified).** Variational assimilation on a *transient* Navier–Stokes model recovers high-resolution flow from noisy, low-resolution PC-MRI. It was tested in 2D+t and 3D+t with pulsatile flow in cerebral aneurysms. Funke SW, Nordaas M, Evju Ø, Alnæs MS, Mardal KA. "Variational data assimilation for transient blood flow simulations: cerebral aneurysms as an illustrative example." *Int J Numer Method Biomed Eng* 2019;35(1):e3152. doi:10.1002/cnm.3152 — [DOI](https://doi.org/10.1002/cnm.3152)
- **Töger et al. 2020: "4D-flow HIKING".** An optimising CFD solver with adjoint-based optimisation, on a discontinuous Galerkin discretisation, fits CFD to 4D flow. It was tested:
  - in synthetic 2D tests;
  - on a pulsatile phantom, against laser particle-image velocimetry, at 1.5 T;
  - in vivo, in a volunteer's intracranial arteries at 7 T.

  Töger J, Zahr MJ, Aristokleous N, Markenroth Bloch K, Carlsson M, Persson PO. "Blood flow imaging by optimal matching of computational fluid dynamics to 4D-flow data." *Magn Reson Med* 2020;84(4):2231–2245. doi:10.1002/mrm.28269 — [DOI](https://doi.org/10.1002/mrm.28269)
- **Resistance-type boundary conditions by optimal control.** Fevola E, Ballarin F, Jiménez-Juan L, Fremes S, Grivet-Talocia S, Rozza G, Triverio P. "An optimal control approach to determine resistance-type boundary conditions from in-vivo data for cardiovascular simulations." *Int J Numer Method Biomed Eng* 2021;37(10):e3516. doi:10.1002/cnm.3516 — [DOI](https://doi.org/10.1002/cnm.3516)
- **Kontogiannis et al. 2022: a Bayesian inverse Navier–Stokes problem that reconstructs and segments at once.**
  - Its prior is a set of Gaussian random fields, and the posterior covariance is approximated by quasi-Newton.
  - It reconstructed and segmented MR velocimetry at SNR ≈ 6, with "negligible errors" against high-SNR images. The authors call this a reduction of the total scanning time by a factor of 27.

  Kontogiannis A, Elgersma SV, Sederman AJ, Juniper MP. "Joint reconstruction and segmentation of noisy velocity images as an inverse Navier–Stokes problem." *J Fluid Mech* 2022;944:A40. doi:10.1017/jfm.2022.503 — [DOI](https://doi.org/10.1017/jfm.2022.503)

#### Reduced-order assimilation
- **Merging CFD and 4D flow by proper orthogonal decomposition (POD).** Bakhshinejad A, Baghaie A, Vali A, Saloner D, Rayz VL, D'Souza RM. "Merging computational fluid dynamics and 4D Flow MRI using proper orthogonal decomposition and ridge regression." *J Biomech* 2017;58:162–173. doi:10.1016/j.jbiomech.2017.05.004 — [DOI](https://doi.org/10.1016/j.jbiomech.2017.05.004)
- **Denoising and super-resolution by POD.** Fathi MF, Bakhshinejad A, Baghaie A, Saloner D, Sacho RH, Rayz VL, D'Souza RM. "Denoising and spatial resolution enhancement of 4D flow MRI using proper orthogonal decomposition and lasso regularization." *Comput Med Imaging Graph* 2018;70:165–172. doi:10.1016/j.compmedimag.2018.07.003 — [DOI](https://doi.org/10.1016/j.compmedimag.2018.07.003)
- **Multi-fidelity data in reduced-order assimilation.** Habibi M, D'Souza RM, Dawson STM, Arzani A. "Integrating multi-fidelity blood flow data with reduced-order data assimilation." *Comput Biol Med* 2021;135:104566. doi:10.1016/j.compbiomed.2021.104566 — [DOI](https://doi.org/10.1016/j.compbiomed.2021.104566)
- **Galarce et al.: reconstruction from Doppler ultrasound with reduced models** (the parametrized-background data-weak lineage).
  - Galarce F, Lombardi D, Mula O. "Reconstructing haemodynamics quantities of interest from Doppler ultrasound imaging." *Int J Numer Method Biomed Eng* 2021;37(2):e3416. doi:10.1002/cnm.3416 — [DOI](https://doi.org/10.1002/cnm.3416)
  - Galarce F, Gerbeau JF, Lombardi D, Mula O. "Fast reconstruction of 3D blood flows from Doppler ultrasound images and reduced models." *Comput Methods Appl Mech Eng* 2021;375:113559. doi:10.1016/j.cma.2020.113559 — [DOI](https://doi.org/10.1016/j.cma.2020.113559)

#### Windkessel and lumped boundary conditions, Bayesian UQ, multi-fidelity (classical)
- **The Windkessel and LPN coupling to be calibrated.**
  - Vignon-Clementel IE, Figueroa CA, Jansen KE, Taylor CA. "Outflow boundary conditions for three-dimensional finite element modeling of blood flow and pressure in arteries." *Comput Methods Appl Mech Eng* 2006;195(29–32):3776–3796. doi:10.1016/j.cma.2005.04.014 — [DOI](https://doi.org/10.1016/j.cma.2005.04.014)
  - Kim HJ, Vignon-Clementel IE, Figueroa CA, et al. "On coupling a lumped parameter heart model and a three-dimensional finite element aorta model." *Ann Biomed Eng* 2009;37(11):2153–2169. doi:10.1007/s10439-009-9760-8 — [DOI](https://doi.org/10.1007/s10439-009-9760-8)
- **Sankaran and Marsden 2011 (the user's citation, verified).** Sankaran S, Marsden AL. "A stochastic collocation method for uncertainty quantification and propagation in cardiovascular simulations." *J Biomech Eng* 2011;133(3):031001. doi:10.1115/1.4003259 — [DOI](https://doi.org/10.1115/1.4003259)
- **A guide to UQ and sensitivity analysis.** Eck VG, Donders WP, Sturdy J, Feinberg J, Delhaas T, Hellevik LR, Huberts W. "A guide to uncertainty quantification and sensitivity analysis for cardiovascular applications." *Int J Numer Method Biomed Eng* 2016;32(8):e02755. doi:10.1002/cnm.2755 — [DOI](https://doi.org/10.1002/cnm.2755)
- **UQ for virtual surgery.** Schiavazzi DE, Arbia G, Baker C, Hlavacek AM, Hsia TY, Marsden AL, Vignon-Clementel IE; MOCHA Investigators. "Uncertainty quantification in virtual surgery hemodynamics predictions for single ventricle palliation." *Int J Numer Method Biomed Eng* 2016;32(3):e02737. doi:10.1002/cnm.2737 — [DOI](https://doi.org/10.1002/cnm.2737)
- **Schiavazzi et al. 2017 (the user's citation, verified).** Schiavazzi DE, Baretta A, Pennati G, Hsia TY, Marsden AL. "Patient-specific parameter estimation in single-ventricle lumped circulation models under uncertainty." *Int J Numer Method Biomed Eng* 2017;33(3):e02799. doi:10.1002/cnm.2799 — [DOI](https://doi.org/10.1002/cnm.2799)
- **Tran et al. (the user's "Tran et al.", verified).**
  - Tran JS, Schiavazzi DE, Ramachandra AB, Kahn AM, Marsden AL. "Automated tuning for parameter identification and uncertainty quantification in multi-scale coronary simulations." *Comput Fluids* 2017;142:128–138. doi:10.1016/j.compfluid.2016.05.015 — [DOI](https://doi.org/10.1016/j.compfluid.2016.05.015)
  - Tran JS, Schiavazzi DE, Kahn AM, Marsden AL. "Uncertainty quantification of simulated biomechanical stimuli in coronary artery bypass grafts." *Comput Methods Appl Mech Eng* 2019;345:402–428. doi:10.1016/j.cma.2018.10.024 — [DOI](https://doi.org/10.1016/j.cma.2018.10.024)
- **Gaussian-process emulators for a 1D model.** Melis A, Clayton RH, Marzo A. "Bayesian sensitivity analysis of a 1D vascular model with Gaussian process emulators." *Int J Numer Method Biomed Eng* 2017;33(12):e2882. doi:10.1002/cnm.2882 — [DOI](https://doi.org/10.1002/cnm.2882)
- **The uncertainty of a pressure drop across a coarctation.** Brüning J, Hellmeier F, Yevtushenko P, Kühne T, Goubergrits L. "Uncertainty quantification for non-invasive assessment of pressure drop across a coarctation of the aorta using CFD." *Cardiovasc Eng Technol* 2018;9(4):582–596. doi:10.1007/s13239-018-00381-3 — [DOI](https://doi.org/10.1007/s13239-018-00381-3)
- **Model discrepancy and model selection.** Paun LM, Colebank MJ, Olufsen MS, Hill NA, Husmeier D. "Assessing model mismatch and model selection in a Bayesian uncertainty quantification analysis of a fluid-dynamics model of pulmonary blood circulation." *J R Soc Interface* 2020;17(173):20200886. doi:10.1098/rsif.2020.0886 — [DOI](https://doi.org/10.1098/rsif.2020.0886)
- **Multi-fidelity estimators.**
  - Seo J, Fleeter C, Kahn AM, Marsden AL, Schiavazzi DE. "Multifidelity estimators for coronary circulation models under clinically informed data uncertainty." *Int J Uncertain Quantif* 2020;10(5):449–466. doi:10.1615/Int.J.UncertaintyQuantification.2020033068 — [DOI](https://doi.org/10.1615/int.j.uncertaintyquantification.2020033068)
  - Fleeter CM, Geraci G, Schiavazzi DE, Kahn AM, Marsden AL. "Multilevel and multifidelity uncertainty quantification for cardiovascular hemodynamics." *Comput Methods Appl Mech Eng* 2020;365:113030. doi:10.1016/j.cma.2020.113030 — [DOI](https://doi.org/10.1016/j.cma.2020.113030)
- **A review of UQ in haemodynamics.** Ninos G, Bartzis V, Merlemis N, Sarris IE. "Uncertainty quantification implementations in human hemodynamic flows." *Comput Methods Programs Biomed* 2021;203:106021. doi:10.1016/j.cmpb.2021.106021 — [DOI](https://doi.org/10.1016/j.cmpb.2021.106021)
- **Nolte and Bertoglio 2022: a review of inverse problems in blood flow.** It covers the last decade of inverse problems for fully dimensional fluid and fluid–solid models from clinical images. Nolte D, Bertoglio C. "Inverse problems in blood flow modeling: a review." *Int J Numer Method Biomed Eng* 2022;38(8):e3613. doi:10.1002/cnm.3613 — [DOI](https://doi.org/10.1002/cnm.3613)
- **Automatic 0D and 1D models from 3D geometry.** Pfaller MR, Pham J, Verma A, Pegolotti L, Wilson NM, Parker DW, Yang W, Marsden AL. "Automated generation of 0D and 1D reduced-order models of patient-specific blood flow." *Int J Numer Method Biomed Eng* 2022;38(10):e3639. doi:10.1002/cnm.3639 — [DOI](https://doi.org/10.1002/cnm.3639)

#### Latest (2023–2026)
- **The uncertainty of blood's rheology.** Ranftl S, Müller TS, Windberger U, Brenn G, von der Linden W. "A Bayesian approach to blood rheological uncertainties in aortic hemodynamics." *Int J Numer Method Biomed Eng* 2023;39(4):e3576. doi:10.1002/cnm.3576 — [DOI](https://doi.org/10.1002/cnm.3576)
- **Lumped boundary conditions by vector fitting.** Fevola E, Bradde T, Triverio P, Grivet-Talocia S. "A vector fitting approach for the automated estimation of lumped boundary conditions of 1D circulation models." *Cardiovasc Eng Technol* 2023;14(4):505–525. doi:10.1007/s13239-023-00669-z — [DOI](https://doi.org/10.1007/s13239-023-00669-z)
- **Kontogiannis et al. 2025: a Bayesian inverse Navier–Stokes problem with the boundary as a parameter.** Velocimetry is assimilated to reconstruct a 3D flow jointly with unknown Navier–Stokes parameters, *the boundary's position among them*. The method uses:
  - Gaussian priors on the parameters;
  - a stabilised Nitsche weak form;
  - a viscous signed distance field as an auxiliary variable for the geometry;
  - an adjoint-consistent cut-cell finite element method;

  and it is applied to flow-MRI. Kontogiannis A, Elgersma SV, Sederman AJ, Juniper MP. "Bayesian inverse Navier–Stokes problems: joint flow field reconstruction and parameter learning." *Inverse Problems* 2025;41(1):015008. doi:10.1088/1361-6420/ad9cb7 — [DOI](https://doi.org/10.1088/1361-6420/ad9cb7)
- **Richter et al. 2025: Bayesian Windkessel calibration on an optimised 0D surrogate.**
  - The 3D model is run once, for an initial choice of boundary conditions, and an optimised 0D surrogate is built from it.
  - Sequential Monte Carlo then gives the posterior of all the Windkessel parameters.
  - Optimising the 0D models lowered their median error against 3D "by nearly one order of magnitude in 72 publicly available vascular models."
  - The posteriors were evaluated at several signal-to-noise ratios.

  Richter J, Nitzler J, Pegolotti L, Menon K, Biehler J, Wall WA, Schiavazzi DE, Marsden AL, Pfaller MR. "Bayesian Windkessel calibration using optimized zero-dimensional surrogate models." *Philos Trans R Soc A* 2025;383(2292):20240223. doi:10.1098/rsta.2024.0223 — [DOI](https://doi.org/10.1098/rsta.2024.0223)
- **svZeroDSolver** (BSD-3-Clause), which installs both `svzerodsolver` and `svzerodcalibrator` command-line tools. — [GitHub](https://github.com/SimVascular/svZeroDSolver). Menon K, Richter J, Pfaller MR, Pham J, Mathew EM, Harold KE, Dorn NC, Verma A, Marsden AL. "svZeroDSolver: a modular package for lumped-parameter cardiovascular simulations." *J Open Source Softw* 2025;10(109):7595. doi:10.21105/joss.07595 — [DOI](https://doi.org/10.21105/joss.07595)

### Inferences
- **The two lineages divide the work by what they need.**
  - ROUKF needs only forward runs of an ensemble of 2p+1 sigma points, and is parallel. It suits a C++ solver without adjoints and a few tens of parameters (Windkessels, stiffness).
  - Adjoint 4D-Var needs a differentiable solver. It suits estimating the whole state or an inflow field (Funke 2019; Töger 2020).

  A layer should support both behind one interface: observation operator H, error covariance R, parameter prior, and a choice of estimator.
- **The infrastructure already holds the Kalman analysis step.** tyche's existing Gaussian conditioning, one observation at a time (CLAUDE.md, ADR 0054), is that step. ROUKF can be built on it, with the sigma-point propagation a scheduler of forward runs.
- **Calibration should run on surrogates, with their error kept.** Richter 2025's pattern is cheap: one 3D run, an optimised 0D model, then SMC. Multi-fidelity estimators (Seo 2020; Fleeter 2020) and GP emulators (Melis 2017) are the alternatives. Each needs the surrogate's own error kept in the likelihood, which links to model discrepancy (Paun 2020).
- **Geometry is the largest and least often propagated uncertainty** (Sankaran 2016; Maher 2021). Kontogiannis 2025 shows how the boundary can become a parameter inferred jointly with the flow.

### Gaps
- **No Kalman paper by D'Elia, Perego and Veneziani dated 2012 was found.** Their verified 2012 works are variational. A tutorial chapter by Bertagna, D'Elia, Perego and Veneziani (c. 2014) on data assimilation in cardiovascular FSI may treat Kalman filters, but was not searched or verified. **[unverified]**
- **The details of Verdandi's current distribution are unverified.** The INRIA GitLab page was behind an access challenge; the licence was read from the GitHub mirror `feelpp/verdandi` (LGPL-2.1). Maintenance status is unknown.
- **Quantitative results of Funke 2019, Tran 2017 and Fleeter 2020 were not retrieved** (abstracts truncated or missing): estimation errors and speed-ups of multi-fidelity estimators.

## 3. Physics-Informed Neural Networks and Neural Operators for Vascular Flow: Reported Accuracy and Limits

### Takeaway
PINNs act in vascular flow as *inverse solvers* and as filters for super-resolution and denoising, not as faster forward solvers.

- **Inverse solvers:**
  - Hidden fluid mechanics recovers velocity and pressure from a scalar visualised in the flow (Raissi, Yazdani & Karniadakis 2020).
  - Kissas et al. 2020 predict pressure from 4D flow with 1D PINNs.
  - Arzani, Wang & D'Souza 2021 recover near-wall WSS from sparse data with unknown boundary conditions.
- **Super-resolution:** 4DFlowNet (2020), Fathi et al. 2020, and in 2023–2024 their physics-informed, implicit-representation and ensemble successors.

Learned surrogates of CFD report errors of a few per cent on in-distribution synthetic cohorts:

| Surrogate | Reported error |
|---|---|
| SSM-driven aorta networks | 6.0% pressure, 4.0% velocity |
| WSSNet | 4.3% relative WSS error against CFD; 7.1% on synthetic MRI |
| Li 2021 (coronary bypass) | "around 90%" accuracy |

Independent assessments warn of weak baselines and reporting bias in machine learning for fluid PDEs (McGreivy & Hakim 2024) and that PINNs rarely beat finite elements on forward problems (Grossmann et al. 2024). Generalisation across geometries (GINO, DIMON, equivariant mesh networks) is the open problem of 2023–2026.

### Cited Findings

#### Foundations (classical)
- **PINNs.** Raissi M, Perdikaris P, Karniadakis GE. "Physics-informed neural networks: a deep learning framework for solving forward and inverse problems involving nonlinear partial differential equations." *J Comput Phys* 2019;378:686–707. doi:10.1016/j.jcp.2018.10.045 — [DOI](https://doi.org/10.1016/j.jcp.2018.10.045). Reference code `maziarraissi/PINNs`: **MIT**. — [LICENSE](https://raw.githubusercontent.com/maziarraissi/PINNs/master/LICENSE)
- **Hidden fluid mechanics.** The method uses the Navier–Stokes equations to extract velocity and pressure from flow visualisations (smoke, contrast agents); blood flow in an aneurysm is among its examples. Raissi M, Yazdani A, Karniadakis GE. "Hidden fluid mechanics: learning velocity and pressure fields from flow visualizations." *Science* 2020;367(6481):1026–1030. doi:10.1126/science.aaw4741 — [DOI](https://doi.org/10.1126/science.aaw4741). Code `maziarraissi/HFM`: **MIT**. — [LICENSE](https://raw.githubusercontent.com/maziarraissi/HFM/master/LICENSE)
- **Kissas et al. 2020 (the user's citation, verified).** The networks are constrained by a 1D pulsatile flow model, conserving mass and momentum. They are trained on "noisy and scattered clinical data of flow and wall displacement," including "in-vivo measurements near the aorta/carotid bifurcation of a healthy human subject." They return "physically consistent predictions for velocity, pressure and wall displacement pulse wave propagation," and post-processing yields Windkessel parameters. The abstract gives no numerical errors. — [arXiv:1905.04817](https://arxiv.org/abs/1905.04817). Kissas G, Yang Y, Hwuang E, Witschey WR, Detre JA, Perdikaris P. "Machine learning in cardiovascular flows modeling: predicting arterial blood pressure from non-invasive 4D flow MRI data using physics-informed neural networks." *Comput Methods Appl Mech Eng* 2020;358:112623. doi:10.1016/j.cma.2019.112623 — [DOI](https://doi.org/10.1016/j.cma.2019.112623)
- **A surrogate trained without simulation data.** Sun L, Gao H, Pan S, Wang JX. "Surrogate modeling for fluid flows based on physics-constrained deep learning without simulation data." *Comput Methods Appl Mech Eng* 2020;361:112732. doi:10.1016/j.cma.2019.112732 — [DOI](https://doi.org/10.1016/j.cma.2019.112732)
- **Arzani, Wang and D'Souza 2021 (the user's citation, verified).** With inlet and outlet conditions unknown, "PINN could use sparse velocity measurements away from the wall to quantify WSS with very high accuracy." Arzani A, Wang JX, D'Souza RM. "Uncovering near-wall blood flow from sparse data with physics-informed neural networks." *Phys Fluids* 2021;33(7):071905. doi:10.1063/5.0055600 — [DOI](https://doi.org/10.1063/5.0055600). Its code `amir-cardiolab/PINN-wss` (PyTorch) covers 1D advection–diffusion, a 2D stenosis, 2D and 3D aneurysms and the identification of viscosity, and has **no LICENSE file at its root** (default: all rights reserved). — [README](https://raw.githubusercontent.com/amir-cardiolab/PINN-wss/main/README.md)
- **Physics-informed CNNs for super-resolution without high-resolution labels.** Gao H, Sun L, Wang JX. "Super-resolution and denoising of fluid flow using physics-informed convolutional neural networks without high-resolution labels." *Phys Fluids* 2021;33(7):073603. doi:10.1063/5.0054312 — [DOI](https://doi.org/10.1063/5.0054312)
- **Brain haemodynamics.** Sarabian M, Babaee H, Laksari K. "Physics-informed neural networks for brain hemodynamic predictions using medical imaging." *IEEE Trans Med Imaging* 2022;41(9):2285–2303. doi:10.1109/TMI.2022.3161653 — [DOI](https://doi.org/10.1109/tmi.2022.3161653)
- **Operator learning: DeepONet and its physics-informed form.**
  - Lu L, Jin P, Pang G, Zhang Z, Karniadakis GE. "Learning nonlinear operators via DeepONet based on the universal approximation theorem of operators." *Nat Mach Intell* 2021;3(3):218–229. doi:10.1038/s42256-021-00302-5 — [DOI](https://doi.org/10.1038/s42256-021-00302-5)
  - Wang S, Wang H, Perdikaris P. "Learning the solution operator of parametric partial differential equations with physics-informed DeepONets." *Sci Adv* 2021;7(40):eabi8605. doi:10.1126/sciadv.abi8605 — [DOI](https://doi.org/10.1126/sciadv.abi8605)
- **The training pathology of PINNs.** Wang S, Teng Y, Perdikaris P. "Understanding and mitigating gradient flow pathologies in physics-informed neural networks." *SIAM J Sci Comput* 2021;43(5):A3055–A3081. doi:10.1137/20M1318043 — [DOI](https://doi.org/10.1137/20m1318043)
- **DeepXDE.** Lu L, Meng X, Mao Z, Karniadakis GE. "DeepXDE: a deep learning library for solving differential equations." *SIAM Rev* 2021;63(1):208–228. doi:10.1137/19M1274067 — [DOI](https://doi.org/10.1137/19m1274067). Licence **LGPL-2.1**. — [LICENSE](https://raw.githubusercontent.com/lululxvi/deepxde/master/LICENSE)
- **Reviews.**
  - Arzani A, Dawson STM. "Data-driven cardiovascular flow modelling: examples and opportunities." *J R Soc Interface* 2021;18(175):20200802. doi:10.1098/rsif.2020.0802 — [DOI](https://doi.org/10.1098/rsif.2020.0802)
  - Arzani A, Wang JX, Sacks MS, Shadden SC. "Machine learning for cardiovascular biomechanics modeling: challenges and beyond." *Ann Biomed Eng* 2022;50(6):615–627. doi:10.1007/s10439-022-02967-4 — [DOI](https://doi.org/10.1007/s10439-022-02967-4)

#### Super-resolution of 4D flow (classical)
- **4DFlowNet** (the user's citation, verified). Ferdian E, Suinesiaputra A, Dubowitz DJ, Zhao D, Wang A, Cowan B, Young AA. "4DFlowNet: super-resolution 4D flow MRI using deep learning and computational fluid dynamics." *Front Phys* 2020;8:138. doi:10.3389/fphy.2020.00138 — [DOI](https://doi.org/10.3389/fphy.2020.00138). Code `EdwardFerdian/4DFlowNet`: **MIT**. — [LICENSE](https://raw.githubusercontent.com/EdwardFerdian/4DFlowNet/master/LICENSE)
- **Physics-informed super-resolution and denoising.** Fathi MF, Perez-Raya I, Baghaie A, Berg P, Janiga G, Arzani A, D'Souza RM. "Super-resolution and denoising of 4D-Flow MRI using physics-informed deep neural nets." *Comput Methods Programs Biomed* 2020;197:105729. doi:10.1016/j.cmpb.2020.105729 — [DOI](https://doi.org/10.1016/j.cmpb.2020.105729)
- **Cerebrovascular 4D flow enhanced by learning on CFD.** Rutkowski DR, Roldán-Alzate A, Johnson KM. "Enhancement of cerebrovascular 4D flow MRI velocity fields using machine learning and computational fluid dynamics simulation data." *Sci Rep* 2021;11(1):10240. doi:10.1038/s41598-021-89636-z — [DOI](https://doi.org/10.1038/s41598-021-89636-z)
- **WSSNet**: accuracy given in §1. Code `EdwardFerdian/WSSNet`: **MIT**. — [LICENSE](https://raw.githubusercontent.com/EdwardFerdian/WSSNet/master/LICENSE)

#### Learned surrogates of CFD (classical)
- **The thoracic aorta.** Liang L, Mao W, Sun W. "A feasibility study of deep learning for predicting hemodynamics of human thoracic aorta." *J Biomech* 2020;99:109544. doi:10.1016/j.jbiomech.2019.109544 — [DOI](https://doi.org/10.1016/j.jbiomech.2019.109544)
- **Li et al. 2021: coronary bypass.** The network agreed with CFD at "prediction accuracy of around 90%" over more than 2 million nodes, 600-fold faster than CFD and within 1 s. Li G, Wang H, Zhang M, Tupin S, Qiao A, Liu Y, Ohta M, Anzai H. "Prediction of 3D cardiovascular hemodynamics before and after coronary artery bypass surgery via deep learning." *Commun Biol* 2021;4(1):99. doi:10.1038/s42003-020-01638-1 — [DOI](https://doi.org/10.1038/s42003-020-01638-1)
- **A synthetic database for training.** Thamsen B, Yevtushenko P, Gundelwein L, Setio AAA, Lamecker H, Kelm M, Schafstedde M, Heimann T, Kuehne T, Goubergrits L. "Synthetic database of aortic morphometry and hemodynamics: overcoming medical imaging data availability." *IEEE Trans Med Imaging* 2021;40(5):1438–1449. doi:10.1109/TMI.2021.3057496 — [DOI](https://doi.org/10.1109/tmi.2021.3057496)

#### Latest (2023–2026)
- **Physics-informed super-resolution of cardiac 4D flow.** Shone F, Ravikumar N, Lassila T, MacRaild M, Wang Y, Taylor ZA, Jimack P, Dall'Armellina E, Frangi AF. "Deep physics-informed super-resolution of cardiac 4D-flow MRI." In: *Information Processing in Medical Imaging (IPMI 2023)*, LNCS. Springer; 2023:511–522. doi:10.1007/978-3-031-34048-2_39 — [DOI](https://doi.org/10.1007/978-3-031-34048-2_39)
- **Cerebrovascular super-resolution followed by physics-informed processing.** Ferdian E, Marlevi D, Schollenberger J, Aristova M, Edelman ER, Schnell S, Figueroa CA, Nordsletten DA, Young AA. "Cerebrovascular super-resolution 4D Flow MRI – sequential combination of resolution enhancement by deep learning and physics-informed image processing to non-invasively quantify intracranial velocity, flow, and relative pressure." *Med Image Anal* 2023;88:102831. doi:10.1016/j.media.2023.102831 — [DOI](https://doi.org/10.1016/j.media.2023.102831)
- **Pajaziti et al. 2023: networks driven by a statistical shape model.**
  - 3,000 CFD runs were made on synthetic aortas drawn from a statistical shape model (SSM) of 67 real aortas.
  - On 200 test shapes the errors were 6.01% ± 3.12 for pressure and 3.99% ± 0.93 for velocity.
  - Inference took about 0.075 s, "4,000x faster than the solver."

  Pajaziti E, Montalt-Tordera J, Capelli C, Sivera R, Sauvage E, Quail M, Schievano S, Muthurangu V. "Shape-driven deep neural networks for fast acquisition of aortic 3D pressure and velocity flow fields." *PLoS Comput Biol* 2023;19(4):e1011055. doi:10.1371/journal.pcbi.1011055 — [DOI](https://doi.org/10.1371/journal.pcbi.1011055)
- **Geometry-informed neural operators (GINO).** Li Z, Kovachki N, Choy C, et al. "Geometry-informed neural operator for large-scale 3D PDEs." *Advances in Neural Information Processing Systems 36* (NeurIPS 2023):35836–35854. doi:10.52202/075280-1556 — [DOI](https://doi.org/10.52202/075280-1556)
- **Graph neural networks as reduced-order models.** Pegolotti L, Pfaller MR, Rubio NL, Ding K, Brugarolas Brufau R, Darve E, Marsden AL. "Learning reduced-order models for cardiovascular simulations with graph neural networks." *Comput Biol Med* 2024;168:107676. doi:10.1016/j.compbiomed.2023.107676 — [DOI](https://doi.org/10.1016/j.compbiomed.2023.107676)
- **Equivariant mesh networks on the artery wall.** Suk J, de Haan P, Lippe P, Brune C, Wolterink JM. "Mesh neural networks for SE(3)-equivariant hemodynamics estimation on the artery wall." *Comput Biol Med* 2024;173:108328. doi:10.1016/j.compbiomed.2024.108328 — [DOI](https://doi.org/10.1016/j.compbiomed.2024.108328)
- **Implicit neural representations, unsupervised.** Saitta S, Carioni M, Mukherjee S, Schönlieb CB, Redaelli A. "Implicit neural representations for unsupervised super-resolution and denoising of 4D flow MRI." *Comput Methods Programs Biomed* 2024;246:108057. doi:10.1016/j.cmpb.2024.108057 — [DOI](https://doi.org/10.1016/j.cmpb.2024.108057)
- **Super-resolution generalised across the circulation by ensembles.** Ericsson L, Hjalmarsson A, Akbar MU, Ferdian E, Bonini M, Hardy B, Schollenberger J, Aristova M, Winter P, Burris N, Fyrdahl A, Sigfridsson A, et al. "Generalized super-resolution 4D flow MRI – using ensemble learning to extend across the cardiovascular system." *IEEE J Biomed Health Inform* 2024;28(12):7239–7250. doi:10.1109/JBHI.2024.3429291 — [DOI](https://doi.org/10.1109/jbhi.2024.3429291)
- **DIMON: operators over geometries, with cardiac electrophysiology among its examples.** Yin M, Charon N, Brody R, Lu L, Trayanova N, Maggioni M. "A scalable framework for learning the geometry-dependent solution operators of partial differential equations." *Nat Comput Sci* 2024;4(12):928–940. doi:10.1038/s43588-024-00732-2 — [DOI](https://doi.org/10.1038/s43588-024-00732-2)
- **Limits of PINNs and of ML for fluids.**
  - Grossmann et al. compare PINNs and FEM systematically on Poisson, Allen–Cahn and semilinear Schrödinger equations, for cost and accuracy. Grossmann TG, Komorowska UJ, Latz J, Schönlieb CB. "Can physics-informed neural networks beat the finite element method?" *IMA J Appl Math* 2024;89(1):143–174. doi:10.1093/imamat/hxae011 — [DOI](https://doi.org/10.1093/imamat/hxae011)
  - McGreivy N, Hakim A. "Weak baselines and reporting biases lead to overoptimism in machine learning for fluid-related partial differential equations." *Nat Mach Intell* 2024;6(10):1256–1269. doi:10.1038/s42256-024-00897-5 — [DOI](https://doi.org/10.1038/s42256-024-00897-5)
- **Neural differential equations for 1D blood flow.** Csala H, Mohan A, Livescu D, Arzani A. "Physics-constrained coupled neural differential equations for one dimensional blood flow modeling." *Comput Biol Med* 2025;186:109644. doi:10.1016/j.compbiomed.2024.109644 — [DOI](https://doi.org/10.1016/j.compbiomed.2024.109644)
- **A benchmark of geometric deep learning on coronary meshes.** Nannini G, Suk J, Rygiel P, Saitta S, Mariani L, Maranga R, Baggiano A, Pontone G, Wolterink JM, Redaelli A. "Learning hemodynamic scalar fields on coronary artery meshes: a benchmark of geometric deep learning models." *Comput Biol Med* 2025;195:110477. doi:10.1016/j.compbiomed.2025.110477 — [DOI](https://doi.org/10.1016/j.compbiomed.2025.110477)
- **2026 works whose metadata only was seen in Crossref** (contents not read):
  - Osman D, Sel K, Spatz E, Jafari R. "Cardiovascular digital twins using a Windkessel physics informed neural network." *npj Digit Med* 2026;9(1):443. doi:10.1038/s41746-026-02610-9 — [DOI](https://doi.org/10.1038/s41746-026-02610-9)
  - Elhadidy M, Viknesh S, D'Souza RM, Arzani A. "Wall shear stress reconstruction from concentration: differentiable physics and physics-informed neural networks." *Phys Fluids* 2026;38(8):081908. doi:10.1063/5.0346659 — [DOI](https://doi.org/10.1063/5.0346659)
- **Preprints, not peer-reviewed** (search snippets only). — [search results](https://arxiv.org/abs/2512.09013v1), [arXiv 2507.13789](https://web3.arxiv.org/pdf/2507.13789), [SSRN](https://doi.org/10.2139/ssrn.5010501)
  - a graph-network surrogate of intracranial aneurysms that reproduces full-field haemodynamics in under a minute per cardiac cycle, validated on out-of-distribution aneurysms of the MATCH challenge (arXiv 2512.09013);
  - a domain-agnostic Fourier neural operator for patient-specific aneurysms (arXiv 2507.13789);
  - an E(3)-equivariant estimator of WSS in abdominal aortic aneurysms (Rygiel et al., arXiv 2507.22817);
  - Villié A, et al. "Physics-informed neural networks for enhancing 4D flow MRI: correcting artifacts and assimilating mean pressure and Reynolds stresses in turbulent blood flow," SSRN preprint 2024, doi:10.2139/ssrn.5010501.
- **Libraries.**
  - NVIDIA PhysicsNeMo (formerly Modulus): **Apache-2.0** — [LICENSE](https://raw.githubusercontent.com/NVIDIA/physicsnemo/main/LICENSE.txt)
  - `neuraloperator`: **MIT** — [LICENSE](https://raw.githubusercontent.com/neuraloperator/neuraloperator/main/LICENSE)
  - Transolver: **MIT** — [LICENSE](https://raw.githubusercontent.com/thuml/Transolver/main/LICENSE)

### Inferences
- **Inside a C++ assimilation layer, PINNs are better framed as one more inverse estimator than as a solver.** Arzani 2021 and Kontogiannis 2022/2025 address the same problem, inferring flow and WSS with unknown boundary conditions; the adjoint finite element method gives a posterior covariance as well. Where PINNs are used, they belong in the Python works, on PyTorch (BSD), DeepXDE (LGPL) or PhysicsNeMo (Apache).
- **Learned surrogates belong in calibration and UQ, with their errors carried in the likelihood.** Their reported few-per-cent errors are in-distribution, on synthetic cohorts drawn from statistical shape models. They are emulators for posteriors (Richter 2025; Seo 2020), not replacements for verified CFD in a credibility case.
- **Out-of-distribution geometry is the main limit.** Every 2023–2026 paper on generalisation (GINO, DIMON, equivariant meshes, the benchmark of Nannini 2025) treats it as open. Any surrogate should report its errors on held-out anatomies and its applicability domain, in the sense of Pathmanathan 2017 (§4).
- **All the super-resolution codes checked (4DFlowNet, WSSNet) are MIT**, so they can be vendored or ported. PINN-wss has no licence and cannot be reused without permission.

### Gaps
- **Quantitative accuracy was not retrieved** for Kissas 2020, 4DFlowNet, Fathi 2020, Pegolotti 2024, Suk 2024, Ericsson 2024 or Saitta 2024 (no abstracts in Crossref); read the papers for their error figures.
- **FNO, MeshGraphNets and Krishnapriyan were not verified in Crossref:**
  - Li Z, Kovachki N, Azizzadenesheli K, et al. "Fourier neural operator for parametric partial differential equations" (ICLR 2021; arXiv:2010.08895) **[unverified]**;
  - Pfaff T et al. "Learning mesh-based simulation with graph networks" (ICLR 2021) **[unverified]**;
  - Krishnapriyan A et al. "Characterizing possible failure modes in physics-informed neural networks" (NeurIPS 2021) **[unverified]**.
- **Licences of Kissas's code** (PredictiveIntelligenceLab) and of Pegolotti's GNN code were not found.

## 4. Cardiovascular Digital Twins: Position Papers, Roadmaps, VVUQ Frameworks and the Regulatory Context

### Takeaway
The cardiovascular digital twin is framed by Corral-Acero et al. 2020 (*Eur Heart J*) and Niederer, Sacks, Girolami & Willcox 2021 (*Nat Comput Sci*), and generalised in NASEM's *Foundational Research Gaps and Future Directions for Digital Twins* (2024, doi:10.17226/26894). NASEM defines a digital twin as a model coupled to its physical counterpart and dynamically updated by bidirectional data flows. Credibility is governed by ASME V&V 40-2018 and FDA's final guidance of November 2023 on the credibility of computational modelling in device submissions. Patient-specific models add the questions of applicability and of personalisation (Pathmanathan 2017; Galappaththige 2022; Aycock 2024). Community CFD challenges (the FDA nozzle; SBC2012; the 2015 aneurysm challenge) are the classical benchmarks of variability.

### Cited Findings

#### Position papers and frameworks (classical)
- **Corral-Acero et al. 2020: the digital twin of precision cardiology.** A position paper on the digital twin as the second pillar of precision medicine, emphasising "the synergies between mechanistic and statistical models." Corral-Acero J, Margara F, Marciniak M, Rodero C, Loncaric F, Feng Y, Gilbert A, Fernandes JF, Bukhari HA, Wajdan A, Martinez MV, Santos MS, et al. (39 authors). "The 'Digital Twin' to enable the vision of precision cardiology." *Eur Heart J* 2020;41(48):4556–4564. doi:10.1093/eurheartj/ehaa159 — [DOI](https://doi.org/10.1093/eurheartj/ehaa159)
- **Computational models in cardiology.** Niederer SA, Lumens J, Trayanova NA. "Computational models in cardiology." *Nat Rev Cardiol* 2019;16(2):100–111. doi:10.1038/s41569-018-0104-y — [DOI](https://doi.org/10.1038/s41569-018-0104-y)
- **Niederer et al. 2021 (the user's citation, verified).** Niederer SA, Sacks MS, Girolami M, Willcox K. "Scaling digital twins from the artisanal to the industrial." *Nat Comput Sci* 2021;1(5):313–320. doi:10.1038/s43588-021-00072-5 — [DOI](https://doi.org/10.1038/s43588-021-00072-5)
- **Digital twins at scale as probabilistic graphical models.** Kapteyn MG, Pretorius JVR, Willcox KE. "A probabilistic graphical model foundation for enabling predictive digital twins at scale." *Nat Comput Sci* 2021;1(5):337–347. doi:10.1038/s43588-021-00069-0 — [DOI](https://doi.org/10.1038/s43588-021-00069-0)
- **A twin of the systemic circulation by inverse analysis.** Chakshu NK, Sazonov I, Nithiarasu P. "Towards enabling a cardiovascular digital twin for human systemic circulation using inverse analysis." *Biomech Model Mechanobiol* 2021;20(2):449–465. doi:10.1007/s10237-020-01393-6 — [DOI](https://doi.org/10.1007/s10237-020-01393-6)
- **Electrophysiological twins from the ECG.** Gillette K, Gsell MAF, Prassl AJ, et al. "A framework for the generation of digital twins of cardiac electrophysiology from clinical 12-leads ECGs." *Med Image Anal* 2021;71:102080. doi:10.1016/j.media.2021.102080 — [DOI](https://doi.org/10.1016/j.media.2021.102080)
- **The health digital twin, reviewed.** Coorey G, Figtree GA, Fletcher DF, et al. "The health digital twin to tackle cardiovascular disease—a review of an emerging interdisciplinary field." *npj Digit Med* 2022;5(1):126. doi:10.1038/s41746-022-00640-7 — [DOI](https://doi.org/10.1038/s41746-022-00640-7)
- **Mathematical and numerical foundations.**
  - Quarteroni A, Manzoni A, Vergara C. "The cardiovascular system: mathematical modelling, numerical algorithms and clinical applications." *Acta Numerica* 2017;26:365–590. doi:10.1017/S0962492917000046 — [DOI](https://doi.org/10.1017/s0962492917000046)
  - Morris PD, Narracott A, von Tengg-Kobligk H, et al. "Computational fluid dynamics modelling in cardiovascular medicine." *Heart* 2016;102(1):18–28. doi:10.1136/heartjnl-2015-308044 — [DOI](https://doi.org/10.1136/heartjnl-2015-308044)
  - Alastruey J, Xiao N, Fok H, Schaeffter T, Figueroa CA. "On the impact of modelling assumptions in multi-scale, subject-specific models of aortic haemodynamics." *J R Soc Interface* 2016;13(119):20160073. doi:10.1098/rsif.2016.0073 — [DOI](https://doi.org/10.1098/rsif.2016.0073)
  - Xiao N, Alastruey J, Figueroa CA. "A systematic comparison between 1-D and 3-D hemodynamics in compliant arterial models." *Int J Numer Method Biomed Eng* 2014;30(2):204–231. doi:10.1002/cnm.2598 — [DOI](https://doi.org/10.1002/cnm.2598)

#### VVUQ and credibility (classical)
- **The applicability of validation evidence.** Pathmanathan P, Gray RA, Romero VJ, Morrison TM. "Applicability analysis of validation evidence for biomedical computational models." *J Verif Valid Uncertain Quantif* 2017;2(2):021005. doi:10.1115/1.4037671 — [DOI](https://doi.org/10.1115/1.4037671)
- **A special issue on VVUQ of cardiovascular models.** Steinman DA, Migliavacca F. "Editorial: special issue on verification, validation, and uncertainty quantification of cardiovascular models: towards effective VVUQ for translating cardiovascular modelling to clinical utility." *Cardiovasc Eng Technol* 2018;9(4):539–543. doi:10.1007/s13239-018-00393-z — [DOI](https://doi.org/10.1007/s13239-018-00393-z)
- **The diversity of personalisation.** Gray RA, Pathmanathan P. "Patient-specific cardiovascular computational modeling: diversity of personalization and challenges." *J Cardiovasc Transl Res* 2018;11(2):80–88. doi:10.1007/s12265-018-9792-2 — [DOI](https://doi.org/10.1007/s12265-018-9792-2)
- **What clinical decision support needs.** Huberts W, Heinen SGH, Zonnebeld N, et al. "What is needed to make cardiovascular models suitable for clinical decision support? A viewpoint paper." *J Comput Sci* 2018;24:68–84. doi:10.1016/j.jocs.2017.07.006 — [DOI](https://doi.org/10.1016/j.jocs.2017.07.006)
- **VVUQ for regulatory in silico trials.** Viceconti M, Pappalardo F, Rodriguez B, Horner M, Bischoff J, Musuamba Tshinanu F. "In silico trials: verification, validation and uncertainty quantification of predictive models used in the regulatory evaluation of biomedical products." *Methods* 2021;185:120–127. doi:10.1016/j.ymeth.2020.01.011 — [DOI](https://doi.org/10.1016/j.ymeth.2020.01.011)
- **Galappaththige et al. 2022: credibility for patient-specific models.** They work out how VVUQ applies to patient-specific models, the basis of digital twins, using cardiac electrophysiology as the exemplar; earlier credibility frameworks targeted devices or generic patients. Galappaththige S, Gray RA, Costa CM, Niederer S, Pathmanathan P. "Credibility assessment of patient-specific computational modeling using patient-specific cardiac modeling as an exemplar." *PLoS Comput Biol* 2022;18(10):e1010541. doi:10.1371/journal.pcbi.1010541 — [DOI](https://doi.org/10.1371/journal.pcbi.1010541)
- **Benchmarks and challenges.**
  - **The FDA nozzle.** Hariharan P, Giarra M, Reddy V, Day SW, Manning KB, Deutsch S, Stewart SFC, Myers MR, Berman MR, Burgreen GW, Paterson EG, Malinauskas RA. "Multilaboratory particle image velocimetry analysis of the FDA benchmark nozzle model to support validation of computational fluid dynamics simulations." *J Biomech Eng* 2011;133(4):041002. doi:10.1115/1.4003440 — [DOI](https://doi.org/10.1115/1.4003440)
  - **The SBC2012 challenge.** In phase I, 25 groups agreed remarkably on pressure, the majority predicting peak systolic pressure drops "within 8% of each other," while the aneurysm's flow patterns varied more. Steinman DA, Hoi Y, Fahy P, Morris L, Walsh MT, Aristokleous N, Anayiotos AS, Papaharilaou Y, Arzani A, Shadden SC, Berg P, Janiga G, et al. "Variability of computational fluid dynamics solutions for pressure and flow in a giant aneurysm: the ASME 2012 Summer Bioengineering Conference CFD challenge." *J Biomech Eng* 2013;135(2):021016. doi:10.1115/1.4023382 — [DOI](https://doi.org/10.1115/1.4023382)
  - **The 2015 aneurysm challenge.** Valen-Sendstad K, Bergersen AW, Shimogonya Y, Goubergrits L, Bruening J, Pallares J, Cito S, Piskin S, Pekkan K, Geers AJ, Larrabide I, Rapaka S, et al. "Real-world variability in the prediction of intracranial aneurysm wall shear stress: the 2015 International Aneurysm CFD Challenge." *Cardiovasc Eng Technol* 2018;9(4):544–564. doi:10.1007/s13239-018-00374-2 — [DOI](https://doi.org/10.1007/s13239-018-00374-2)
  - **Backflow at open boundaries.** Bertoglio C, Caiazzo A, Bazilevs Y, et al. "Benchmark problems for numerical treatment of backflow at open boundaries." *Int J Numer Method Biomed Eng* 2018;34(2):e2918. doi:10.1002/cnm.2918 — [DOI](https://doi.org/10.1002/cnm.2918)
  - **Numerical schemes for 1D arterial flow.** Boileau E, Nithiarasu P, Blanco PJ, et al. "A benchmark study of numerical schemes for one-dimensional arterial blood flow modelling." *Int J Numer Method Biomed Eng* 2015;31(10):e02732. doi:10.1002/cnm.2732 — [DOI](https://doi.org/10.1002/cnm.2732)

#### Latest (2023–2026)
- **NASEM 2024.** *Foundational Research Gaps and Future Directions for Digital Twins.* Washington, DC: The National Academies Press; 2024. doi:10.17226/26894. — [NAP](https://www.nationalacademies.org/publications/26894); [study page](https://www.nationalacademies.org/our-work/foundational-research-gaps-and-future-directions-for-digital-twins); [VPH Institute note](https://vph-institute.org/news/nasem-report-on-foundational-research-gaps-and-future-directions-for-digital-twins.html)
  - It is a consensus study sponsored by DoD, DOE, NIH and NSF.
  - It defines a digital twin as coupling computational models with a physical counterpart, "dynamically updated through bidirectional data flows."
  - It recommends that federal agencies build crosscutting programmes on the mathematical, statistical and computational foundations.
  - Its ISBN appears as 9780309700450 in one listing and 9780309700429 in another.
- **GPU-accelerated heart twins.** Viola et al. replicate the heart's full multiphysics "within a few hours per heartbeat" on GPUs, for synthetic cohorts; their proof of concept is left bundle branch block and cardiac resynchronisation. Viola F, Del Corso G, De Paulis R, Verzicco R. "GPU accelerated digital twins of the human heart open new routes for cardiovascular research." *Sci Rep* 2023;13(1):8230. doi:10.1038/s41598-023-34098-8 — [DOI](https://doi.org/10.1038/s41598-023-34098-8)
- **Digital twins in medicine.** Laubenbacher R, Mehrad B, Shmulevich I, Trayanova N. "Digital twins in medicine." *Nat Comput Sci* 2024;4(3):184–191. doi:10.1038/s43588-024-00607-6 — [DOI](https://doi.org/10.1038/s43588-024-00607-6)
- **Sel et al. 2024: building cardiovascular twins.** A review of building twins with "rapid bidirectional interactions" between the virtual and the physical, from multimodal and multiscale data. Sel K, Osman D, Zare F, Masoumi Shahrbabak S, Brattain L, Hahn J, Inan OT, Mukkamala R, Palmer J, Paydarfar D, Pettigrew RI, Quyyumi AA, et al. "Building digital twins for cardiovascular health: from principles to clinical impact." *J Am Heart Assoc* 2024;13(19):e031981. doi:10.1161/JAHA.123.031981 — [DOI](https://doi.org/10.1161/jaha.123.031981)
- **Twins in the era of generative AI.** Thangaraj PM, Benson SH, Oikonomou EK, Asselbergs FW, Khera R. "Cardiovascular care with digital twin technology in the era of generative artificial intelligence." *Eur Heart J* 2024;45(45):4808–4821. doi:10.1093/eurheartj/ehae619 — [DOI](https://doi.org/10.1093/eurheartj/ehae619)
- **Whole-heart emulation by latent neural ODEs.** Salvador M, Strocchi M, Regazzoni F, et al. "Whole-heart electromechanical simulations using latent neural ordinary differential equations." *npj Digit Med* 2024;7(1):90. doi:10.1038/s41746-024-01084-x — [DOI](https://doi.org/10.1038/s41746-024-01084-x)
- **Aycock et al. 2024: credibility for in silico clinical trials.** A hierarchical framework for in silico clinical trials whose submodels are of mixed kinds (physics-based, data-driven, rule-based). Aycock KI, Battisti T, Peterson A, Yao J, Kreuzer S, Capelli C, Pant S, Pathmanathan P, Hoganson DM, Levine SM, Craven BA. "Toward trustworthy medical device in silico clinical trials: a hierarchical framework for establishing credibility and strategies for overcoming key challenges." *Front Med* 2024;11:1433372. doi:10.3389/fmed.2024.1433372 — [DOI](https://doi.org/10.3389/fmed.2024.1433372)
- **FDA's final guidance and ASME V&V 40.** The guidance is titled "Assessing the Credibility of Computational Modeling and Simulation in Medical Device Submissions." — [HPN](https://www.hpnonline.com/regulatory/news/53078879/fda-issues-final-guidance-assessing-the-credibility-of-computational-modeling-and-simulation-in-medical-device-submissions); [MedTech Intelligence](https://medtechintelligence.com/news_article/fda-publishes-final-guidance-on-assessing-credibility-of-cms-in-device-submissions/amp); [NIBIB IMAG summary](https://www.imagwiki.nibib.nih.gov/node/538691); [FDA document (not opened)](https://www.fda.gov/media/175363/download)
  - It was issued final in November 2023, covered in the trade press on 2023-11-21, with an FDA webinar on 2024-01-11.
  - It follows a December 2021 draft.
  - It sets a risk-informed credibility framework drawing on the FDA-recognised ASME V&V 40-2018, applicable to any physics-based or mechanistic model.
  - FDA's webinar slides note that V&V 40 implicitly assumes validation against prospective, well-controlled bench tests.
- **A Windkessel PINN twin** (2026, metadata only): Osman et al., *npj Digit Med* 2026 (§3).

### Inferences
- **For imaging-constrained haemodynamics, the "bidirectional data flow" of NASEM's definition is exactly the assimilation layer (§2).** Its credibility case under V&V 40 and the FDA guidance needs four parts:
  - **Code verification:** benchmarks, the method of manufactured solutions, and the community cases above.
  - **Validation in the space of the measurement** (§1).
  - **Applicability analysis** for each new patient.
  - **UQ of the predicted quantities of interest**, with geometry, boundary conditions and rheology propagated.
- **V&V 40 is risk-informed**: how much evidence is needed depends on the model's influence on a decision and on its consequence. An open infrastructure should therefore record a context of use with each assimilated twin, not one global accuracy figure.
- **Variability between codes is a lower bound on credible accuracy.** The challenges show pressure drops agreeing within about 8% while flow patterns and WSS vary much more.

### Gaps
- **The text of ASME V&V 40-2018 and the FDA guidance itself were not opened** (paywalled or not fetched). The guidance's exact docket identifiers and its sections on patient-specific models are unverified.
- **No 2025–2026 consensus roadmap specific to cardiovascular digital twins was located.** An ESC or EHRA statement, or the European Virtual Human Twin (EDITH) roadmap, may exist. **[unverified]**
- **The regulatory history of FFR-CT** (e.g. FDA's de novo classification of HeartFlow) was not verified in this session.

## 5. Open-Source Tools for Integration, Their Licences, and Atlases of Normal Anatomy as Priors

### Takeaway
The permissive core that fits an MIT infrastructure:

- **Solvers and models:** SimVascular and svMultiPhysics (BSD), svZeroDSolver with its calibrator (BSD-3), vmtk (BSD; its optional TetGen dependency carries a modified MIT licence), openBF (Apache-2.0).
- **Learning:** the PINN and super-resolution codes (MIT: HFM, PINNs, 4DFlowNet, WSSNet, neuraloperator, Transolver), PhysicsNeMo (Apache-2.0).
- **UQ:** chaospy (MIT).

Copyleft must be flagged:

- **LGPL**, usable as a dynamically linked dependency: Verdandi (LGPL-2.1), DeepXDE (LGPL-2.1), dolfin-adjoint/pyadjoint (LGPL-3.0), DOLFINx (LGPL-3.0-or-later).
- **GPL-3.0:** Oasis, turtleFSI, VaMPy, morphMan.
- **AGPL-3.0 with the Commons Clause:** CRIMSON's GUI. The clause forbids selling, so CRIMSON is not open source in the OSI sense.

Atlases of normal anatomy (BodyParts3D, the Human Reference Atlas) give labelled anatomy. The priors actually used in haemodynamics are statistical shape models and synthetic virtual cohorts built from them.

### Cited Findings

#### Licences, read from each repository's LICENSE file on 2026-11-04
| Tool | Role in the integration layer | Licence | Source |
|---|---|---|---|
| SimVascular | image → model → CFD pipeline | BSD ("SimVascular is provided under the BSD License") | [LICENSE.md](https://raw.githubusercontent.com/SimVascular/SimVascular/master/LICENSE.md) |
| svMultiPhysics (svFSI successor) | 3D CFD/FSI solver | BSD-style 3-clause text (© 2016–2025 Stanford, UC Regents) | [LICENSE.txt](https://raw.githubusercontent.com/SimVascular/svMultiPhysics/main/LICENSE.txt) |
| svZeroDSolver (+ `svzerodcalibrator`) | 0D LPN solver and calibration | BSD-3-Clause | [GitHub](https://github.com/SimVascular/svZeroDSolver) |
| CRIMSON GUI | pipeline incl. ROUKF parameter estimation (Arthurs 2020) | **AGPL-3.0 + Commons Clause** (no right to "Sell"); MITK parts BSD-style | [LICENSE.txt](https://raw.githubusercontent.com/carthurs/CRIMSONGUI/master/LICENSE.txt) |
| Verdandi | data assimilation (ROUKF, etc.) | **LGPL-2.1** (GitHub mirror `feelpp/verdandi`) | [license](https://raw.githubusercontent.com/feelpp/verdandi/master/license) |
| vmtk | centrelines, meshing for vascular CFD | BSD (TetGen optional, modified MIT) | [LICENSE](https://raw.githubusercontent.com/vmtk/vmtk/master/LICENSE) |
| morphMan | morphological manipulation of vessels | **GPL-3.0** | [LICENSE](https://raw.githubusercontent.com/KVSlab/morphMan/master/LICENSE) |
| VaMPy | vascular CFD pipeline on Oasis | **GPL-3.0** (some files differ) | [LICENSE](https://raw.githubusercontent.com/KVSlab/VaMPy/master/LICENSE) |
| Oasis | FEniCS Navier–Stokes solver | **GPL-3.0** | [COPYING](https://raw.githubusercontent.com/mikaem/Oasis/master/COPYING) |
| turtleFSI | FSI solver | **GPL-3.0** | [LICENSE](https://raw.githubusercontent.com/KVSlab/turtleFSI/master/LICENSE) |
| DOLFINx | finite elements | **LGPL-3.0-or-later** (README) | [README](https://raw.githubusercontent.com/FEniCS/dolfinx/main/README.md) |
| pyadjoint (dolfin-adjoint) | automatic adjoints for 4D-Var | **LGPL-3.0** | [LICENSE](https://raw.githubusercontent.com/dolfin-adjoint/pyadjoint/main/LICENSE) |
| openBF | 1D arterial network solver | Apache-2.0 (© 2023 INSIGNEO) | [LICENSE](https://raw.githubusercontent.com/INSIGNEO/openBF/master/LICENSE) |
| chaospy | polynomial chaos and UQ | MIT | [LICENSE.txt](https://raw.githubusercontent.com/jonathf/chaospy/master/LICENSE.txt) |
| DeepXDE | PINN library | **LGPL-2.1** | [LICENSE](https://raw.githubusercontent.com/lululxvi/deepxde/master/LICENSE) |
| NVIDIA PhysicsNeMo | PINNs and operators | Apache-2.0 | [LICENSE.txt](https://raw.githubusercontent.com/NVIDIA/physicsnemo/main/LICENSE.txt) |
| neuraloperator | FNO and related | MIT | [LICENSE](https://raw.githubusercontent.com/neuraloperator/neuraloperator/main/LICENSE) |
| Transolver | transformer PDE surrogate | MIT | [LICENSE](https://raw.githubusercontent.com/thuml/Transolver/main/LICENSE) |
| HFM; PINNs (Raissi) | reference PINN code | MIT | [HFM](https://raw.githubusercontent.com/maziarraissi/HFM/master/LICENSE); [PINNs](https://raw.githubusercontent.com/maziarraissi/PINNs/master/LICENSE) |
| 4DFlowNet; WSSNet | 4D flow super-resolution; WSS | MIT | [4DFlowNet](https://raw.githubusercontent.com/EdwardFerdian/4DFlowNet/master/LICENSE); [WSSNet](https://raw.githubusercontent.com/EdwardFerdian/WSSNet/master/LICENSE) |
| PINN-wss (Arzani) | near-wall PINN | **no licence file** (all rights reserved by default) | [README](https://raw.githubusercontent.com/amir-cardiolab/PINN-wss/main/README.md) |
| Healthy Aorta 4D flow atlases (Northwestern) | normal 4D flow atlas data | **CC BY-NC-SA 3.0 US** | [PRISM](https://prism.northwestern.edu/records/r29kc-x5336) |

- **Papers for the tools.**
  - **SimVascular.** Updegrove A, Wilson NM, Merkow J, Lan H, Marsden AL, Shadden SC. "SimVascular: an open source pipeline for cardiovascular simulation." *Ann Biomed Eng* 2017;45(3):525–541. doi:10.1007/s10439-016-1762-8 — [DOI](https://doi.org/10.1007/s10439-016-1762-8)
  - **CRIMSON.** Its pipeline runs from segmentation through meshing and boundary-condition design to Navier–Stokes and FSI simulation and post-processing. Arthurs CJ, Khlebnikov R, Melville A, Marčan M, Gomez A, Dillon-Murphy D, Cuomo F, Silva Vieira M, Schollenberger J, Lynch SR, Tossas-Betancourt C, Iyer K, et al. "CRIMSON: an open-source software framework for cardiovascular integrated modelling and simulation." *PLoS Comput Biol* 2021;17(5):e1008881. doi:10.1371/journal.pcbi.1008881 — [DOI](https://doi.org/10.1371/journal.pcbi.1008881)
  - **dolfin-adjoint.** Mitusch S, Funke S, Dokken J. "dolfin-adjoint 2018.1: automated adjoints for FEniCS and Firedrake." *J Open Source Softw* 2019;4(38):1292. doi:10.21105/joss.01292 — [DOI](https://doi.org/10.21105/joss.01292)
  - **chaospy.** Feinberg J, Langtangen HP. "Chaospy: an open source tool for designing methods of uncertainty quantification." *J Comput Sci* 2015;11:46–57. doi:10.1016/j.jocs.2015.08.008 — [DOI](https://doi.org/10.1016/j.jocs.2015.08.008)
  - **svZeroDSolver** (JOSS 2025) and **Verdandi** (Chapelle 2013): §2.

#### Atlases of normal anatomy and shape priors
- **BodyParts3D.** Mitsuhashi N, Fujieda K, Tamura T, Kawamoto S, Takagi T, Okubo K. "BodyParts3D: 3D structure database for anatomical concepts." *Nucleic Acids Res* 2009;37(Database issue):D782–D785. doi:10.1093/nar/gkn613 — [DOI](https://doi.org/10.1093/nar/gkn613)
  - Its licence is **CC BY-SA 2.1 Japan (share-alike)**, per DBCLS's archive README.
  - The required credit reads: "BodyParts3D, Copyright© 2008 The Database Center for Life Science licensed by CC Attribution-Share Alike 2.1 Japan."

  — [DBCLS README](https://dbarchive.biosciencedbc.jp/data/bodyparts3d/20110915/README_e.html)
- **The Human Reference Atlas.** Börner K, Teichmann SA, Quardokus EM, Gee JC, et al. "Anatomical structures, cell types and biomarkers of the Human Reference Atlas." *Nat Cell Biol* 2021;23(11):1117–1128. doi:10.1038/s41556-021-00788-6 — [DOI](https://doi.org/10.1038/s41556-021-00788-6).
  - The 3D reference objects of HuBMAP's Common Coordinate Framework are released under **CC BY 4.0**.
  - They are initially built from the Visible Human male and female datasets of the National Library of Medicine, and distributed as GLB.
  - The brain comes from the Allen Human Reference Atlas–3D (2020); check its own terms.

  — [CCF 3D Reference Object Library](https://hubmapconsortium.github.io/ccf/pages/ccf-3d-reference-library.html) Its blood-vasculature tables appear in a preprint: Boppana A, Lee S, Malhotra R, et al. "Anatomical structures, cell types, and biomarkers of the healthy human blood vasculature." bioRxiv 2022. doi:10.1101/2022.02.28.482302 — [DOI](https://doi.org/10.1101/2022.02.28.482302)
- **Shape priors and virtual cohorts used in haemodynamics.**
  - Bruse JL, McLeod K, Biglino G, et al.; MOCHA Collaborative Group. "A statistical shape modelling framework to extract 3D shape biomarkers from medical imaging data: assessing arch morphology of repaired coarctation of the aorta." *BMC Med Imaging* 2016;16(1):40. doi:10.1186/s12880-016-0142-z — [DOI](https://doi.org/10.1186/s12880-016-0142-z)
  - Strocchi M, Augustin CM, Gsell MAF, et al. "A publicly available virtual cohort of four-chamber heart meshes for cardiac electro-mechanics simulations." *PLoS One* 2020;15(6):e0235145. doi:10.1371/journal.pone.0235145 — [DOI](https://doi.org/10.1371/journal.pone.0235145)
  - Rodero C, Strocchi M, Marciniak M, et al. "Linking statistical shape models and simulated function in the healthy adult human heart." *PLoS Comput Biol* 2021;17(4):e1008851. doi:10.1371/journal.pcbi.1008851 — [DOI](https://doi.org/10.1371/journal.pcbi.1008851)
  - Romero P, Lozano M, Martínez-Gil F, et al. "Clinically-driven virtual patient cohorts generation: an application to aorta." *Front Physiol* 2021;12:713118. doi:10.3389/fphys.2021.713118 — [DOI](https://doi.org/10.3389/fphys.2021.713118)
  - Thamsen et al. 2021 (§3), a synthetic aortic database from an SSM.
  - Pajaziti et al. 2023 (§3), an SSM of 67 aortas giving 3,000 CFD cases.
  - Hoeijmakers MJMM, Huberts W, Rutten MCM, van de Vosse FN. "The impact of shape uncertainty on aortic-valve pressure-drop computations." *Int J Numer Method Biomed Eng* 2021;37(10):e3518. doi:10.1002/cnm.3518 — [DOI](https://doi.org/10.1002/cnm.3518)
- **Geometry from images.**
  - Radl L, Jin Y, Pepe A, et al. "AVT: multicenter aortic vessel tree CTA dataset collection with ground truth segmentation masks." *Data Brief* 2022;40:107801. doi:10.1016/j.dib.2022.107801 — [DOI](https://doi.org/10.1016/j.dib.2022.107801)
  - Kong F, Wilson N, Shadden S. "A deep-learning approach for direct whole-heart mesh reconstruction." *Med Image Anal* 2021;74:102222. doi:10.1016/j.media.2021.102222 — [DOI](https://doi.org/10.1016/j.media.2021.102222)
  - Kong F, Shadden SC. "Learning whole heart mesh generation from patient images for computational simulations." *IEEE Trans Med Imaging* 2023;42(2):533–545. doi:10.1109/TMI.2022.3219284 — [DOI](https://doi.org/10.1109/tmi.2022.3219284)

### Inferences
- **For an MIT C++/Python infrastructure, the dependencies that can be embedded or ported freely** are:
  - **solvers and models:** SimVascular, svMultiPhysics, svZeroDSolver, vmtk, openBF;
  - **learning:** PhysicsNeMo, neuraloperator, Transolver;
  - **UQ:** chaospy;
  - **4D flow super-resolution:** 4DFlowNet and WSSNet.

  LGPL libraries (Verdandi, DeepXDE, pyadjoint, DOLFINx) can be called as separately installed dependencies, with no code copied in. GPL-3.0 code (Oasis, turtleFSI, VaMPy, morphMan) should be used only as external tools that read and write files, for instance for cross-verification. CRIMSON's AGPL with the Commons Clause rules out reuse of its code: re-implement ROUKF from Moireau & Chapelle 2011 instead.
- **Normal-anatomy atlases are not, in the literature found, used directly as priors for vascular CFD.** They supply labelled reference geometry and terminology. The working prior is a statistical shape model (a mean and modes, with a Gaussian on the coefficients) built from a cohort (Bruse 2016; Romero 2021; Thamsen 2021; Pajaziti 2023). A proposed design: BodyParts3D or HRA geometry as the template that is registered and deformed, and the SSM's Gaussian as the prior on geometry in the Bayesian inverse problem (cf. Kontogiannis 2025's boundary as a parameter, and Maher 2021's geometric uncertainty).

### Gaps
- **The atlases' licences were confirmed from search results quoting their pages**, not by downloading the licence files: BodyParts3D CC BY-SA 2.1 JP; HRA's 3D objects CC BY 4.0, the brain excepted.
  - Under share-alike, derivative *meshes* of BodyParts3D must stay CC BY-SA. Keep them in the data cache, not in the MIT code.
  - HRA's CC BY 4.0 permits redistribution with attribution.
- **No published work was found that uses BodyParts3D or the HRA directly as a prior for haemodynamic simulation.** This absence is from a limited search, not a systematic review.
- **Not checked:** the licences of CMRsim, Kissas's 1D PINN code, Pegolotti's GNN code, lifex-cfd, Firedrake, the AVT dataset and the Vascular Model Repository.
