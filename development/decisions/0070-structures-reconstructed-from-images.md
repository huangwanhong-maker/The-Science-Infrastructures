# 0070. Structures Reconstructed from Images: eikon, the Medical-Imaging Domain, Fabrica, and Vessels First

- **Status:** Accepted, provisionally
- **Date:** 2026-11-04
- **Amends:** [ADR 0066](0066-the-anatomy-domain-and-the-body-atlas.md) (the Body Atlas renamed Fabrica, and its scope widened), [ADR 0061](0061-the-domains-of-life.md) (physiology's dependencies), [ADR 0035](0035-the-mathematical-foundations.md) (a sixth foundation)
- **Decided by:** the user, who asked for "features that can utilize integrative imaging data to reconstruct physical organs, and phenomena on top of it … a reusable and separate module", vessels first, "as we interested in fluid dynamics"; agreed that analysis is "not a single domain. Because in fact they utilize some reusable, common algorithms", that "medical-imaging only is ok" as a domain, that physiology be founded with circulation, and that fluids belong in a foundation ("It is foundation for many things"); renamed the Body Atlas Fabrica ("visualize body components … from various medical imaging tech data integratively reconstruct them, and algorithms analysis them, visualize phenomena/prediction/simulation"; diagnosis out of scope); and asked for the stage's task: "module(s) and algorithm(s) for integratively reconstruct body structures from (heterogeneous kinds of) imaging data, and vein as the first", after a survey. Claude surveyed the literature, chose the layout, the methods and their checks, and wrote this record.

## 1. Context

The Body Atlas (ADR 0066) shows one standard body, BodyParts3D's man. The user's aim is the body of a person as images show it: structures reconstructed from CT, MR and other images together, analysed, and carrying phenomena, first the flow of blood in veins. A survey of classical and recent work (`research/literature/reports/Vessel reconstruction and blood flow.md`, six notes beside it) found:

- that the reconstructed geometry, not the flow solver, limits what can be computed on a vessel: 26 teams segmenting the same aneurysms differed by up to 20% in volume, while 25 groups' pressure drops on one geometry agreed within 8%;
- that every vesselness filter's defaults were insufficient and bifurcations their common failure (Lamy et al. 2022);
- that topology needs its own measures (clDice, Betti numbers) beside overlap;
- that labelled venous data are scarce, and that licences sort the tools: TetGen, CGAL's meshing, most lattice-Boltzmann codes and vesselFM cannot enter the MIT core.

The infrastructure had no library for images. What reconstruction needs is mostly not medical: scale space, tubular features, distances, thinning, surfaces, registration and fusion serve a river in a satellite image or a filament in a telescope's as well as a vessel in a CT.

## 2. Decision

### 2.1 A Foundation for Images: eikon

The methods of image analysis are a sixth foundation, [eikon](https://github.com/huangwanhong-maker/eikon) (Greek *εἰκών*, image), on calx, morphe and physis, in its own repository as the submodule `foundations/eikon` (its decisions 0001–0003). Its first ten modules, each checked against exact answers on phantoms (116 checks):

| Module | Content |
|---|---|
| `grid` | oriented grids, sampling, resampling, NIfTI-1 |
| `phantoms` | tubes, trees and tori with exact distances, partial volumes, noise |
| `filters` | Gaussian scale space, the γ-normalised Hessian |
| `features` | Frangi, Sato and Jerman vesselness, the answering scale and axis |
| `segmentation` | Otsu, components, growing, hysteresis, holes |
| `distance` | exact distance transform with nearest points, morphology, second-order fast marching, minimal paths |
| `skeleton` | homotopic thinning, the skeleton as morphe's multigraph with radii, pruning |
| `surfaces` | marching tetrahedra to physis's triangle surfaces, Taubin smoothing, pieces |
| `registration` | rigid registration by mutual information; ICP with Horn's quaternion |
| `fusion` | STAPLE; Dice, surface distances, clDice, Betti numbers |

The closed-form symmetric 3 × 3 eigenproblem, which names no image, went to calx (its decision 0015).

### 2.2 The Medical-Imaging Domain

`domains/medical-imaging` is a field at height 1, on anatomy and the foundations calx, morphe, physis and eikon. It holds only what is medical:

- `dicom`: DICOM files read and series assembled in the patient frame LPS, with Pixel Spacing read row spacing first, slice spacing from positions, gaps and gantry tilt reported; NIfTI's RAS carried to LPS through anatomy's frames.
- `modalities`: CT's Hounsfield units, windows and forward model; TOF-MRA's inflow signal by the spoiled gradient echo, with venous saturation.
- `reconstruction`: vessels from one image. They are detected by vesselness and delineated at the local half-maximum. They are thinned to the anatomy domain's vascular tree, with radii, and bounded by a closed surface. Presets exist for contrast CT and TOF-MRA.
- `integration`: vessels from several images of one body. Each image is registered to a reference by mutual information, the evidence is fused by STAPLE, and each branch is classified as artery or vein by what each modality shows.

On a contrast CT of a known tree, the reconstruction scores Dice 0.94 and clDice 0.99. It finds the tree's 7 junctions and 9 ends, with radii within 0.10 mm on average, and fits Murray's exponent within 0.25 of the truth's 3. On a CT and a TOF-MRA of one body, the second image moved by 6°, the images register within 0.25 mm, and arteries and veins are told apart at Dice 0.87 and 0.94. Real CT is checked against the Medical Segmentation Decathlon's hepatic vessels (§2.5).

### 2.3 Vascular Trees in Anatomy

The anatomy domain gains `vessels`: vascular trees of segments with centrelines and radii, rooted, with generations and Strahler orders, bifurcation angles, tortuosity, Murray's exponent fitted rather than assumed, and naming by an atlas. It names the objects; the medical-imaging domain fills them from images.

### 2.4 Physiology Founded; Fluids a Foundation Module

Physiology (ADR 0061) is founded with `circulation` as its first module, for the next stage: flow in the reconstructed vessels. Its declared dependencies gain anatomy; its height stays 2. The general mechanics of fluids becomes a module of physis, `fluids` (its decision 0021, next stage): reduced-order (0D, 1D with the venous tube law) and three-dimensional (lattice Boltzmann on the voxels), with transport. These serve the blood, the sea, the air and the circulation of value alike. The survey's §8 orders that work.

### 2.5 Data

The data library reads members of a remote tar archive by range requests (`data/acquisition/remotetar.py`), indexing the archive's headers once. The registry's `msd-hepatic-vessels` takes four cases of the Medical Segmentation Decathlon's Task 08 (CC BY-SA 4.0) with 1.5 mm slices, about 130 MB of the 9.35 GB archive. The cases were chosen after reading every case's spacing from its header, since the smallest files have 7.5 mm slices, where vessels are not resolved. `data/products/msd_vessels.py` unpacks them for the C++. On four real portal-venous CTs, against the experts' labels within the annotated territory, the classical chain reaches a mean Dice of 0.32 and clDice of 0.32. Two cases come near 0.5; the weakest, near 0.04, has an enhancement of only 82 HU in the parenchyma. This is in the range the survey reports for classical filters, below learned networks, and it is the baseline later work must not fall below. Two lessons came from these cases. Frangi's c must be taken within the region: over a whole crop, the edges of bone and air suppressed every vessel. And the hysteresis must be lower than on phantoms. TOF-MRA of the head (IXI, CC BY-SA 3.0) and the large veins labelled by TotalSegmentator (CC BY 4.0) are planned.

### 2.6 Fabrica

The Body Atlas is renamed **Fabrica** (`applications/fabrica`, port 8990, silver white), after Vesalius's *De humani corporis fabrica*. Its scope:

- to show the body's components;
- to reconstruct them from images of several kinds, integratively;
- to analyse them;
- to show phenomena, predictions and simulations on the reconstructed structures.

Diagnosis is out of scope: it would bring clinical responsibility the infrastructure does not take. Fabrica gains a mode, Reconstruct.

## 3. Alternatives Considered

- **Image analysis in the medical-imaging domain.** Rejected by the user: the algorithms are common to every science with images.
- **A domain of analysis.** Rejected: analysis is methods, not objects.
- **Fluids as a domain.** Rejected by the user: fluids name no object, and serve blood, sea, air and economics.
- **Learned segmentation first.** Deferred. The classical chain is the reference against which learned models are measured. The licences of open models' weights are mixed.
- **Keeping the name Body Atlas.** Rejected by the user: it narrowed the ambition to showing.

## 4. Consequences

- The foundations are six. FOUNDATION_DEPENDS in the dependency checks gains eikon on calx, morphe and physis. The root CMake builds eikon after nomos and medical-imaging after anatomy.
- eikon's GitHub repository is to be created by the user. Until then the submodule has no remote, and the pointer cannot be pushed.
- A program reconstructing a person's vessels from their images can now be composed from tested parts. Its limits are measured, not assumed: on real CT the Decathlon's labels set the baseline that later changes must not fall below.
- Next: flow and transport in the reconstructed veins (physis `fluids`, physiology `circulation`), as the survey's §8 orders them, and their visualisation in Fabrica.
