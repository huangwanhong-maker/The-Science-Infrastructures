# medical-imaging

*Status: draft (2026-11-04). Field, height 1, on anatomy.* Images of the body and the structures reconstructed from them. Founded at the user's request ("medical-imaging only is ok") for Fabrica, with vessels first ([ADR 0070](../../development/decisions/0070-structures-reconstructed-from-images.md)).

## 1. Purpose

- **The domain holds only what is medical.** Every method on images — scale space, vesselness, segmentation, distances, skeletons, surfaces, registration, fusion — is the foundation [eikon](https://github.com/huangwanhong-maker/eikon)'s, and serves every science with images. This domain knows what an image of the body is: the patient frame and how scanners place images in it; what each modality measures and how a vessel appears in it; which of eikon's methods suit a modality, and with what settings; and how modalities that show different things are combined.
- **Reconstruction is checked against a truth.** Each modality has a forward model, built on eikon's phantoms. It makes the image a known body would give, so every pipeline is held to the body it was given. On real images the reference is experts' labels; precision against them is a lower bound, since labels miss small branches.
- **Integration is the point.** A contrast CT shows arteries and veins alike. A TOF-MRA with a venous saturation band shows arteries only. Together, registered and fused, they tell one from the other, which neither does alone.
- **Diagnosis is not its purpose.** The domain reconstructs and measures structures. It does not judge them (ADR 0070).

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `dicom` | DICOM Part 10 files (uncompressed, little endian; compressed transfer syntaxes refused with their UID). A series is sorted along its normal, and its spacing taken from the positions, not Slice Thickness. Pixel Spacing is read row spacing first. Gaps and gantry tilt are reported. The grid lies in the patient frame LPS, and NIfTI's RAS is carried to it through anatomy's frames | draft |
| `modalities` | CT: Hounsfield units, windows, and a forward model with its point-spread function and noise. TOF-MRA: the spoiled gradient echo's approach to steady state, inflow and venous saturation, and a forward model with Rician noise | draft |
| `reconstruction` | Vessels from one image: confined to a region, detected by vesselness with intensity, and delineated at each voxel's local half-maximum. They are thinned to anatomy's vascular tree with radii and bounded by a closed surface. Presets for contrast CT and TOF-MRA | draft |
| `integration` | Vessels from several images: each registered to a reference by mutual information, reconstructed on its own grid, and carried across. The vessels are fused by STAPLE, and each branch is classified as artery or vein by what each modality shows | draft |
| `flow` | Phase-contrast and 4D flow MRI: velocity encoding, aliasing, the measurement operator for comparison and assimilation with simulated flow (the survey's §5) | planned |
| `atlas` | Registration of a person's structures to BodyParts3D's body, naming vessels by the atlas's centrelines | planned |
| `learned` | Segmentation by learned networks (nnU-Net, TotalSegmentator), run outside the core and measured against the classical chain | planned |

## 3. Dependencies

| Source | What medical imaging uses from it |
|---|---|
| anatomy | The patient frames LPS and RAS (`frames`); vascular trees, their measures and naming (`vessels`) |
| eikon | Every method on images: grids and NIfTI, phantoms, filters, vesselness, segmentation, distances, skeletons, surfaces, registration, fusion |
| calx, morphe, physis | Through eikon and anatomy; calx's random streams for the forward models' noise |

## 4. Checks

`node domains/medical-imaging/tests/run-tests.mjs` builds the domain's C++ and runs:

- **`check_dicom`** (11). Files written in explicit and implicit VR, read back to every voxel's position. Pixel Spacing in order, spacing from positions, gaps, tilt, compressed syntaxes refused, RAS to LPS.
- **`check_modalities`** (10). Hounsfield units. The spoiled gradient echo's limits and the Ernst angle. Inflow and saturation. The CT model's edge response and noise.
- **`check_reconstruction`** (8). A contrast CT of a known tree: Dice 0.94, clDice 0.99, its junctions and ends, radii within 0.10 mm, Murray's exponent, a closed surface of the truth's volume.
- **`check_integration`** (4). A CT and a TOF-MRA of one body moved between scans: registration within 0.25 mm; the vessels at Dice 0.91; arteries and veins at 0.87 and 0.94.
- **`check_msd`** (2). Four real portal-venous CTs from the Medical Segmentation Decathlon, held to the experts' labels at the baseline measured on 2026-11-04: mean Dice 0.32 and clDice 0.32. Two cases come near 0.5 and one near 0.04, whose parenchyma at 82 HU is weakly enhanced. The thresholds of the portal-venous preset were tuned on these cases, so a held-out evaluation is open. It skips without the data library's product (`python data/acquisition/fetch.py msd-hepatic-vessels`, then `python data/products/msd_vessels.py`).

## 5. References

- DICOM PS3.3 (Information Object Definitions) C.7.6.2 and 10.7.1.3; PS3.5 (Data Structures and Encoding); PS3.10 (Media Storage and File Format). NEMA.
- Cox, R. W. et al. (2004), A (sort of) new image data format standard: NIfTI-1, the Data Format Working Group.
- Hounsfield, G. N. (1973), Computerized transverse axial scanning (tomography), *Br. J. Radiol.* 46, 1016.
- Haacke, E. M. et al. (1999), *Magnetic Resonance Imaging: Physical Principles and Sequence Design*, Wiley.
- Antonelli, M. et al. (2022), The Medical Segmentation Decathlon, *Nat. Commun.* 13, 4128.
- The survey behind ADR 0070: `research/literature/reports/Vessel reconstruction and blood flow.md`.

## 6. Open Questions

- **An organ's mask.** The real-CT check searches the territory the experts annotated, since the Decathlon gives no liver mask. An organ segmentation (TotalSegmentator's, or a classical one) would make the search independent of the labels.
- **Thick slices.** Vessels of 1–6 mm are not resolved between slices 7.5 mm apart. Resampling, or vesselness with anisotropic scales, may recover some.
- **Veins in particular.** Open venous labels are scarce, and no open venous 4D flow data were found (the survey). A venous verification suite with reference data is a research target.
