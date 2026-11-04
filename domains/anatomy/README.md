# anatomy

*Status: draft. Field, height 0.* The structures of the body, as models and, later, as images. The domain was founded on 2026-11-04 at the user's request, for the Body Atlas ([ADR 0066](../../development/decisions/0066-the-anatomy-domain-and-the-body-atlas.md)). README §4.4 had placed it in the order beforehand, with medical imaging and histology among its modules.

## 1. Purpose

Anatomy names the parts of the body, says where each lies and what it is made of, and measures them. The domain holds this knowledge, and every claim in it can be checked against a body.

- **Names and kinds.** A part is a concept of the Foundational Model of Anatomy (FMA). Each concept has its kinds (is-a) and the wholes it belongs to (part-of).
- **Places.** A part lies in a frame of the body's own directions, so its left is the subject's left. It is cut in the planes in which imaging sees it.
- **Measures.** The foundations measure what the mesh of a part bounds: its volume, area, principal axes, extents and topology.

The first body is BodyParts3D's adult man, 2,234 meshes named by FMA ([data registry](../../data/registry/bodyparts3d.json)). Medical imaging comes next, for reconstructing vessels from scans, and physiology's flow will join it later. The domain's frames and sections are already the planes of DICOM and of CT.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `frames` | The anatomical position and directions; the patient frames LPS (DICOM, BodyParts3D) and RAS (NIfTI, MNI) and their explicit conversion; the cardinal planes; the radiological and neurological displays of a section; where a point lies relative to another | draft |
| `parts` | FMA's concepts; is-a and part-of as forests with their closures, meets and down-sets. Elements below a concept by the tree or by the source's lists (FMA's part-of is a graph, of which BodyParts3D's tree is a spanning tree). The seventeen systems of the body by an ordered table of rules; sides read from names | draft |
| `body` | BodyParts3D's body read from the data library's product: the atlas and the meshes in LPS millimetres | draft |
| `measures` | A part's topology (morphe), its area, its solid (physis: volume, centre, inertia, principal axes), its extents. A compound's measures composed from its elements'. The parts that touch | draft |
| `sections` | The body cut by a cardinal plane in the display's coordinates, each part's pieces with their areas. The vertebral level of a transverse plane, read from the body's own vertebrae | draft |
| `knowledge` | What is known of the parts: Latin names (TA2), function, blood supply, venous drainage, nerves and clinical notes, with references; a seed of fifty entries (`knowledge/structures.json`), inherited along is-a | seed |
| `vessels` | Vascular trees: segments between junctions and ends with their centrelines and radii, rooted (generations, Strahler orders), measured (length, tortuosity, radius away from the junctions, bifurcation angles), Murray's exponent fitted, named by an atlas ([ADR 0070](../../development/decisions/0070-structures-reconstructed-from-images.md)); filled from images by the medical-imaging domain | draft |
| `histology` | Tissues and their cells (README §4.5) | planned |

## 3. Dependencies

| Source | What anatomy uses from it |
|---|---|
| morphe | `topology`: surfaces given as triangles (decision 0010), the classification of surfaces; `combinatorics`: spanning forests |
| physis | `computational_geometry`: triangle surfaces, welding, sections by a plane, close pairs (decision 0019); `mechanics`: the solid a closed surface bounds, `combine`, `principal_axes` |
| calx | `verification` (the checks); `linalg` through physis's principal axes |

The domain declares no other domain (`module.json`). Physiology will stand on it for circulation, and so will clinical medicine.

## 4. Used By

- **The [Body Atlas](../../applications/body-atlas/README.md).** It shows the body by system and draws a card for each part. It cuts the body in the planes of imaging and offers guided tours.

## 5. Conventions

- **Frames.** Coordinates are LPS millimetres unless named otherwise: +x the subject's left, +y posterior, +z superior. RAS is LPS turned by π about the superior axis. `frames::convert` is the only crossing between them.
- **Displays.** A transverse section is shown radiologically, seen from the feet with the subject's left on the image's right, unless the neurological display is asked for. A coronal section is seen from the front, a sagittal one from the subject's left.
- **Sides.** A name's "left" or "right" is the body's side only far from the midline. Near it, the word may name a chamber or a branch: the right ventricle and the right coronary artery's posterior branches lie left of the midline.
- **Units.** Volumes in mm³ and areas in mm², as the meshes' millimetres give them. The database publishes its volumes in cm³.

## 6. Checks

`node domains/anatomy/tests/run-tests.mjs` builds the C++ with CMake (in WSL on Windows) and runs five programs: frames 14, parts 16, body 9, measures 21 and sections 10. It also checks the form of the knowledge base, and holds its keys to the release's concepts.

The checks that read the body skip when the product has not been built. With the product, they hold:
- the volumes to the database's published ones (fixture `tests/fixtures/bodyparts3d-volumes.tsv`, made by `tests/make_bp3d_fixtures.py`);
- the sections to Cavalieri's principle on the left kidney;
- the order of the vertebrae, and the aorta to the left of the vena cava at T10;
- the trees as spanning trees of their relations.

## 7. Open Questions

1. **The female body.** HuBMAP's Human Reference Atlas (CC BY 4.0) has male and female reference organs and a named blood vasculature, and would be a second body for the same domain.
2. **The vessels as a network.** The meshes' contacts suggest how vessels join. Whether centrelines from the meshes or from images come first is for the vein work to decide.
3. **The knowledge base**, a seed to be checked line by line against Gray's Anatomy and extended beyond its fifty entries.
4. **Duplicates and mislabels.** BodyParts3D 4.0 has duplicate copies, which the product marks, and one mislabelled pair (the flexor pollicis brevis, its sides swapped), which the checks name. Its published volumes scatter for thin parts, apparently because they were computed on meshes finer than those distributed.

## 8. References

- Mitsuhashi, N., Fujieda, K., Tamura, T., Kawamoto, S., Takagi, T. and Okubo, K. (2009), BodyParts3D: 3D structure database for anatomical concepts, *Nucleic Acids Research* 37, D782–D785.
- Rosse, C. and Mejino, J. L. V. (2003), A reference ontology for biomedical informatics: the Foundational Model of Anatomy, *Journal of Biomedical Informatics* 36, 478–500.
- FIPAT (2019), *Terminologia Anatomica*, 2nd ed.
- Standring, S. (ed.) (2020), *Gray's Anatomy*, 42nd ed., Elsevier.
- DICOM PS3.3 (2024), C.7.6.2.1.1, the patient coordinate system.
