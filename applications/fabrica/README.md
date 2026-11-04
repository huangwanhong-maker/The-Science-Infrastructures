# Fabrica

*Status: draft (2026-11-04). Serves: learning, use and research. Draws on: anatomy, medical imaging.* Decision records: [ADR 0066](../../development/decisions/0066-the-anatomy-domain-and-the-body-atlas.md) (as the Body Atlas), [ADR 0070](../../development/decisions/0070-structures-reconstructed-from-images.md) (renamed, and reconstruction).

Fabrica, after Vesalius's *De humani corporis fabrica* (1543), is a workbench of the body: it shows the body's components, reconstructs them from images of several kinds together, analyses them, and shows phenomena, predictions and simulations on the structures reconstructed. Diagnosis is outside its scope. It was the Body Atlas until the user renamed it (ADR 0070).

Fabrica shows the whole human body in the round, every part of it clickable for what it is: its name and kind, where it lies, what it measures, what it touches, and what is known of it. It cuts the body in the planes of imaging and draws each section as an anatomical plate. Guided tours run through it. The body is BodyParts3D's adult man, 2,234 meshes each named by a concept of the Foundational Model of Anatomy. Every measure on its cards is the anatomy domain's.

## 1. Concepts

- **One body, its parts and its systems.**
  - **Parts.** Every mesh is a part, a concept of FMA. The parts are browsed by FMA's part-of tree (the human body, its systems and regions) or by its is-a tree (kinds: bone organ, artery, muscle organ).
  - **Systems.** The domain sorts every part into one of seventeen systems by an ordered table of rules, tissue kinds before organ systems. So the vessels of the lung are vessels, and the larynx's cartilages respiratory.
  - **Display.** Each system is shown, drawn as glass or hidden. At first the skin is glass, the muscles and fascia hidden, so the interior is seen.
- **The card of a part.**
  - **What is known.** Its Latin name, what it does, its blood supply, venous drainage and nerves, and a clinical note, with references. These come from the domain's knowledge base, a seed of fifty entries. A part without an entry of its own takes its kind's, and says so.
  - **Where it lies.** Its vertebral level, its height above the soles, its side.
  - **What it measures.** Its volume, surface and extents along its principal axes, beside the volume the database published. Its form: one closed surface of genus 0, or what it is.
  - **What it touches.** The parts within 1.5 mm of it.
  - **In the trees.** Its chains in both trees.
  - **Calibre.** A long part, such as a vessel, has its calibre along its length: the diameter of the circle of each cross-section's area.
- **Choosing a part.** It is lit in rose gold and drawn faintly through whatever lies in front of it, so a part inside the body is found at once.
- **Cutting.** A cut by a transverse, coronal or sagittal plane keeps the half whose face turns to the view. Its face is the engine's section, drawn exactly where the domain cut the meshes.
- **Sections.** The Sections mode draws the cut as a plate, beside the cut body in the round.
  - **Conventions.** As in imaging, a transverse section is seen from the feet, the subject's left on the image's right; or, at a key, from above.
  - **Names.** The larger pieces are named in the margins on leader lines.
  - **Orientation.** The subject's directions are lettered at the edges, with a scale bar of 50 mm.
  - **Levels.** A transverse level is named by its vertebra, read from the body's own vertebrae.
- **Tours.**
  - The body in layers.
  - The heart and the great vessels.
  - The way back from the foot: the venous return from the great saphenous vein to the right atrium.
  - A walk down the trunk in sections, at T4, T7, T10, L1 and L4.

- **Reconstruct** ([ADR 0070](../../development/decisions/0070-structures-reconstructed-from-images.md)). Vessels reconstructed from images by the medical-imaging domain, on the foundation eikon.
  - **Cases.** A simulated body, imaged twice by the domain's forward models: a contrast CT, and a TOF-MRA of the patient moved between scans. Its two images are integrated: registered by mutual information, fused, and arteries told from veins by what each modality shows. Then real portal-venous CTs of the liver from the Medical Segmentation Decathlon, reconstructed within the territory its experts annotated.
  - **In the round.** The vessels found, in silver, or the arteries in red and the veins in blue. The reference, the experts' labels or the simulated truth, as rose glass. The centrelines coloured by radius.
  - **Slices.** The image's slices as in imaging, windowed. The vessels found and the reference are tinted, and the region not searched is darkened.
  - **The card.** The image's grid. The agreement with the reference: Dice, clDice with its centreline precision and sensitivity, and surface distances. For the simulated body, the registration's error and the arteries' and veins' agreement. The vascular tree: its branches, junctions, ends, length, components and loops, and Murray's exponent fitted where enough bifurcations determine it.
  - **The engine.** It reconstructs a case once and keeps it. A real CT takes about a minute.

## 2. Running

From the repository root, fetch the data, build the product and compile the engine with the atlas's recipe (about 145 MB, and two minutes), then serve the page:

```
python data/acquisition/recipe.py applications/fabrica/recipe.json
python applications/fabrica/serve.py --open
```

The page is at http://127.0.0.1:8990/. The first start measures every part once (about 20 seconds) and keeps the catalogue in `data/cache/fabrica/`. The address bar keeps the mode, the part chosen, the cut, the systems and the view, so a link reproduces a view.

- **In the round:** drag to turn, scroll to approach, right-drag to move. Click a part to read it; double-click to fly to it. Escape clears the choice and the isolation.
- **On the plate:** scroll to step through the levels (shift: five at a time); ctrl-scroll to zoom about the pointer; drag to move; double-click to fit again. Click a piece or its name to read it.
- **Anywhere:** type `/` or click *Find a part* for any of the 3,432 concepts that hold parts.

## 3. Data

| Data | Source | Dataset | Licence |
|---|---|---|---|
| BodyParts3D 4.0: 2,234 meshes of an adult man, FMA's is-a and part-of trees, the elements of each compound | The Database Center for Life Science (Mitsuhashi et al. 2009) | `bodyparts3d` | CC BY 4.0 (since 2025-02-25) |

The product (`data/cache/products/body-bp3d*`) is built from these files by `data/products/body_anatomy.py` and is never committed. It holds:
- the catalogue the page reads;
- the tables the domain reads;
- four packs of quantised meshes, raw for the engine and gzipped for the page.

## 4. Checks

```
node applications/fabrica/tests/run-tests.mjs
```

They check what is the atlas's own, and skip, saying so, when the product, the catalogue or the engine is missing:
- names are read as a reader expects;
- LPS maps to the scene as `data.js` says, by a rotation;
- every system has a colour, and arteries and veins stay apart in every vision (the visualization layer's simulation);
- the packs are as long as the product says, and decoded meshes lie in their boxes with unit normals;
- every part and vertebra a tour names exists;
- the search finds what it should.

With the engine, end to end:
- the section at T10 names its level and its loops are closed;
- seen as in imaging, the aorta lies on the image's right of the inferior vena cava, and seen from above on its left;
- the left kidney's profile integrates to its volume;
- the liver's solid parts come to over a litre;
- an unnamed plane is refused with a reason.

The science (frames, parts, systems, measures, sections, the knowledge base's form) is checked with the anatomy domain.

## 5. Structure

```
applications/fabrica/
├── index.html, styles.css    the page and its style (silver white; rose gold for what is chosen)
├── serve.py                  the local server (port 8990): the engine kept alive, the catalogue kept on disk
├── engine.py                 the engine built with g++ (natively or in WSL) and its worker
├── engine/main.cpp           the engine on the anatomy domain: catalogue, section, part, profile
├── recipe.json               the data, the product and the engine
├── src/app.js                loading, choosing, the cut, the sections, the tours, the search, the address
├── src/viewer.js             the round view: three.js, systems' materials, glass, isolation, picking by BVH, the cut's face
├── src/data.js               the product, the packs, the catalogue, the knowledge; LPS to the scene
├── src/card.js               the card of a part, and its calibre
├── src/plate.js              the section as a plate: pieces, names on leader lines, orientation, zoom
├── src/tree.js               the systems, the parts' trees, the search
├── src/tours.js              the tours
├── src/palette.js            the systems' colours and first styles
└── tests/run-tests.mjs
```

## 6. Open Questions

1. **The veins from images.** The user's next work: reconstructing veins from CT or MRI in their patient frames (the domain's planned `imaging`), and following the flow in them (physiology's circulation). The atlas's calibre profiles and LPS sections are the first steps.
2. **A female body,** from HuBMAP's Human Reference Atlas, and its named vasculature.
3. **The knowledge base,** to be checked against Gray's Anatomy and extended.
4. **Lighter meshes.** The muscles of the trunk are finely meshed (the external oblique 112,000 triangles). Simplifying them would halve the download.
