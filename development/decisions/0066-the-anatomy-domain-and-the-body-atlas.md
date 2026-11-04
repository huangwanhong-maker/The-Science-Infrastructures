# 0066. The Anatomy Domain and the Body Atlas: the Whole Body by Its Parts, Measured, and Cut in the Planes of Imaging

- **Status:** Accepted, provisionally. The knowledge base is a seed. Imaging and vessels are planned.
- **Date:** 2026-11-04
- **Decided by:**
  - **The request.** The user asked to begin the infrastructure's first application of the life sciences beyond the brain and the genome: "Let us begin a body atlas? or better names? We visualize various component in anatomy. I think UI major color is silver white is OK." They named the next tasks: vein reconstruction from imaging data, and the simulation and visualization of flow in veins.
  - **The user's choices.** Asked, they chose the name *Body Atlas*, and agreed that the anatomy domain be founded now (CLAUDE.md §6.5).
  - **What Claude did.** It wrote the domain, the foundations' additions (physis decision 0019, morphe decision 0010), the data library's manifest and product, the application, and this record.

## 1. Context

README §4.4 placed anatomy as a field at height 0, "the structures of the body, as images and as models", with medical imaging and histology among its modules. [ADR 0061](0061-the-domains-of-life.md) left it unplaced until a work reached it. The Brain Atlas draws one organ ([ADR 0028](0028-neuroscience-and-the-brain-atlas.md)); nothing drew the body.

Four open models of the whole body were weighed:
- **BodyParts3D** (DBCLS): an adult man of 2,234 meshes, each the representation of a concept of the Foundational Model of Anatomy, with FMA's is-a and part-of trees. Relicensed CC BY 4.0 on 2025-02-25.
- **Z-Anatomy**: built on BodyParts3D, and distributed as a Blender file.
- **HuBMAP's Human Reference Atlas**: CC BY 4.0, male and female reference organs and a named blood vasculature, but not the skeleton and muscles whole.
- **The Visible Human Project's images**: images, not models.

## 2. Decision

### 2.1 The Anatomy Domain (`domains/anatomy`)

A field at height 0, C++ from the start, on calx, morphe and physis. Five modules (its [README](../../domains/anatomy/README.md), its [theory](../../domains/anatomy/theory/anatomy.md)):
- **`frames`.** The anatomical directions; LPS (DICOM, BodyParts3D) and RAS (NIfTI, MNI) with their explicit conversion; the cardinal planes; the radiological and neurological displays.
- **`parts`.** FMA's concepts, the two relations as forests with their closures and meets, and the systems of the body by an ordered table of rules: tissue kinds by is-a before organ systems by part-of, regional exceptions first.
- **`body`.** The reader of the data library's product.
- **`measures`.** Topology by morphe; area, solid, principal axes and extents by physis; compounds by `combine`; touching parts by close pairs.
- **`sections`.** The body cut in the display's coordinates; vertebral levels read from the body's own vertebrae.

A knowledge base of fifty entries (`knowledge/structures.json`) gives Latin names, function, supply, drainage, nerves and clinical notes, with references. It is a seed, inherited along is-a.

### 2.2 What Went to the Foundations

Following the rule that general methods go to the foundations:
- **physis** (decision 0019): triangle surfaces (welding, area, sections by a plane as oriented loops, close pairs of labelled points), and the solid a closed surface bounds (`solid_polyhedron`, composing with the existing `MassProperties`).
- **morphe** (decision 0010): surfaces given as lists of triangles, their 2-complexes read directly in O(F log F), split into components and classified, their orientation checked.

### 2.3 The Data

- **The manifest.** `data/registry/bodyparts3d.json` names the release of 2013-06-19 at its pinned URLs, with checksums.
- **The product.** `data/products/body_anatomy.py` quantises the meshes to 16 bits over each element's box (finer than the files' digits for parts under 650 mm) into four packs, and writes the trees and compounds as tables the domain reads.
- **What the product resolves.**
  - Seventeen headers name no concept; each such element takes the most specific compound containing it.
  - Copies are marked when they share a concept, bounds within 3 mm and volume within 3 % (13 found).

### 2.4 The Body Atlas (`applications/body-atlas`, port 8990)

- **The engine.** C++ on the domain (`engine/main.cpp`), built by `engine.py` and kept running by the server. It answers `catalogue`, `section`, `part` and `profile`. The catalogue (about 20 s) is kept on disk, keyed by the engine's sources and the product.
- **The page.** three.js 0.170.0 and three-mesh-bvh 0.8.3:
  - **Body.** The body in the round by system (shown, glass or hidden), the parts in either tree, a search, and a card for each part: what is known, where it lies (its vertebral level), its measures against the database's, what it touches, its place in both trees, and for long parts their calibre along their length. A cut by any cardinal plane, whose face is the engine's section.
  - **Sections.** A plate in the conventions of imaging, its pieces named on leader lines, and the cut in the round beside it.
  - **Tours.** Four guided tours: the body in layers; the heart and the great vessels; the way back from the foot (the venous return); a walk down the trunk in sections.
- **The style.** Silver white on pearl paper by day and silver on graphite by night, at the user's choice, with a rose-gold minor for what is chosen. Serif type, hairlines. The tissues take softened anatomical colours, and arteries and veins are checked apart in every vision by the visualization layer's simulation.

## 3. Alternatives Considered

- **The Human Reference Atlas first.** Deferred, not rejected. It lacks the whole skeleton and the muscles. Its vasculature and female body make it the natural second body for the vein work (the domain's open question 1).
- **Anatomy as a module of physiology.** Rejected. The map holds it as a field of its own, and physiology (planned, height 2) will stand on it.
- **The knowledge in the page.** Rejected. What is known of the parts names objects in nature, so it belongs to the domain, in a file its checks hold to the release's concepts.
- **Volumes held to the database's published values part by part.** Rejected as a check. Those values scatter for thin parts, apparently because they were computed on finer meshes. The checks hold the median and the parts over 100 cm³, and Cavalieri's principle, which is exact.

## 4. Consequences

- **The order.** It gains anatomy at height 0. The diagrams in CLAUDE.md §3.2, README §3 and §4.6, `domains/README.md` and `docs/architecture.md` show it, and README §4.4 marks it built.
- **Checks.**
  - The domain: 70 C++ checks, and its knowledge base's checks.
  - physis: 29. morphe: 10.
  - The application: 22, among them the engine end to end.
- **Findings about BodyParts3D 4.0, recorded in the domain's checks.**
  - Its part-of tree is a spanning tree of FMA's part-of graph: its compound lists equal the tree's closure in 1,063 of 1,368 cases, and contain it in all.
  - One pair is mislabelled: the flexor pollicis brevis, its sides swapped.
  - A name's side is the body's only far from the midline.
- **Next.** The vein work: imaging (CT volumes in LPS, segmentation, surfaces), vessels as networks with centrelines and calibres, and flow in physiology's circulation, each a module or domain of its own.
