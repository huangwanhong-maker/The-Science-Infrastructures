# Brain Atlas

*Status: draft (2026-10-23). Serves: learning and use. Draws on: neuroscience.* Decision records: [ADR 0028](../../development/decisions/0028-neuroscience-and-the-brain-atlas.md), [ADR 0032](../../development/decisions/0032-the-anatomy-of-the-brain-atlas.md).

The Brain Atlas shows the human brain in the round, every part of it clickable for its anatomy: where it lies, what it does, what it connects, what supplies it with blood, and what follows when it fails. It cuts the brain and shows it in sections, colours its cortex by other divisions, maps data onto it, places the electrodes of the 10–20 system on the scalp around it, and shows the mouse brain beside it. The brain is the unbiased average of 152 adults, the MNI 152 ICBM 2009c template that the Human Brain Project's atlas uses.

<p align="center">
<img src="../../assets/gallery/brain-02-hippocampus.jpg" width="49%" alt="The deep structures through a glass cortex, with the hippocampus and its card">
<img src="../../assets/gallery/brain-03-coronal-section.jpg" width="49%" alt="A coronal cut, its face washed in the colours of the structures">
</p>

## 1. Concepts

- **One brain, many atlases, one grid.** The cortex is a surface made from the template's tissue maps, divided into the 48 gyri of the Harvard–Oxford atlas on each side. The deep structures come from FreeSurfer's segmentation of the template, refined by MASSP 2.0 (7 T MRI of the subcortex), the thalamic nuclei of Najdenovska et al. and Julich-Brain's cerebellar nuclei. The fibre bundles come from the EBRAINS atlas of the deep white matter. All lie on the template's 1 mm grid, so that the round view, the sections and the coordinates agree.
- **A tree of structures.** The neuroscience domain names 152 structures, from the brain down to the gyri, nuclei and bundles, each joined to its surfaces and labels on each side (`brain.js`), and ten functional systems that cut across the tree: the motor and sensory systems, vision, hearing, language, memory, emotion, reward, the neuromodulatory nuclei and the cerebrospinal fluid.
- **What is known of each.** Each structure's card comes from the domain's knowledge base (`knowledge.js`): its Latin name after *Terminologia Neuroanatomica*, where it lies, what it does, its connections, its blood supply, its development and its clinical significance, with references (textbooks and papers checked against Crossref). The card adds what the atlas measures: the volume of each side and their asymmetry, the centre in MNI and Talairach coordinates, the structures it borders, the systems it belongs to, and three small sections through it.
- **Seeing inside.** The cortex may be solid, glass or hidden; a cut along any axis shows the template's section on its face; the parts may be drawn apart. Names are set out beside the brain as on an engraved plate, joined to their structures by fine lines, and follow the view.
- **The modes.**
  - *Anatomy:* the brain in the round, the tree of structures, the systems, five guided tours.
  - *Sections:* the template's image in three planes, every structure washed in its colour.
  - *Atlases:* the cortex coloured by Julich-Brain's cytoarchitectonic areas or by Schaefer's parcels of Yeo's seven networks, each area with where it lies in the gross anatomy; and the mouse brain in the Allen Common Coordinate Framework, with its ontology as a tree.
  - *Data:* a pasted table of structures and values colours the brain (names as the atlas gives them, as FreeSurfer's aseg.stats gives them, or the regions of the atlas shown); pasted points, in MNI or Talairach coordinates, are drawn with what lies at each. Two examples come from the atlas itself: the volume of every structure, and the asymmetry between its sides.
  - *Head:* the template's scalp, cast from its T1 image; the electrodes of the 10–20 or 10–10 system placed on it by the arcs between the nasion, the inion and the ears (Jasper 1958; Oostenveld and Praamstra 2001); a subject's head length, breadth and height scale the brain and the scalp, and the electrodes are placed again; pasted channel values are interpolated over the scalp by spherical splines (Perrin et al. 1989) and carried radially onto the cortex. Each electrode's card names the cortex beneath it.

## 2. Running

From the repository root, fetch the data and build the anatomy (about 110 MB, and a minute to build):

```
python data/acquisition/fetch.py templateflow-mni152 massp-2 harvard-oxford mial-thalamic-nuclei ebrains-deep-white-matter schaefer-2018 ebrains-julich-brain allen-mouse-ccf
python data/products/brain_anatomy.py
python applications/brain-atlas/serve.py --open
```

The page is at http://127.0.0.1:8810/. The address bar keeps the mode, the chosen structure, the cortex's style and colouring, the cut and the view, so that a link reproduces a view.

- **In the round:** drag to turn, scroll to approach, right-drag to move. Click a part to read it; double-click to fly to it.
- **In the sections:** click or drag to set the point, scroll to step through the slices (shift: five at a time); click a section's name to make it the main one.
- **Anywhere:** type in *Find a structure* for any name, Latin name or other name a structure is known by.

## 3. Data

| Data | Source | Dataset | Licence |
|---|---|---|---|
| MNI 152 ICBM 2009c nonlinear asymmetric, 1 mm: T1 image, brain mask, tissue probability maps, FreeSurfer's segmentation | McConnell Brain Imaging Centre (Fonov et al. 2011), via TemplateFlow | `templateflow-mni152` | MNI: permission to use, copy, modify and distribute |
| MASSP 2.0, 31 subcortical structures on each side | Bazin et al. (2020), via TemplateFlow | `massp-2` | CC BY 4.0 |
| Harvard–Oxford cortical atlas, lateralised, 25 % threshold | Makris et al. (2006) and others, FSL, via TemplateFlow | `harvard-oxford` | FSL licence, non-commercial |
| Thalamic nuclei from diffusion MRI | Najdenovska et al. (2018), via TemplateFlow | `mial-thalamic-nuclei` | CC BY 4.0 |
| Deep white-matter bundles | Guevara et al. (2012), EBRAINS | `ebrains-deep-white-matter` | CC BY-NC-SA 4.0 |
| Julich-Brain 3.1, the maximum-probability maps of 207 areas | Amunts et al. (2020), EBRAINS | `ebrains-julich-brain` | CC BY-NC-SA 4.0 |
| Schaefer 2018, 100 and 400 parcels in 7 networks | Schaefer et al. (2018), via TemplateFlow | `schaefer-2018` | MIT |
| Allen Mouse Brain CCFv3, template and 2017 annotation at 50 µm, ontology | Allen Institute (Wang et al. 2020) | `allen-mouse-ccf` | Allen Institute terms of use |

The anatomy product (`data/cache/products/anatomy-mni152.json` with its three packs of meshes and the joined labels) is built from these files by `data/products/brain_anatomy.py`, and is never committed.

## 4. Checks

```
node applications/brain-atlas/tests/run-tests.mjs
```

The checks need the fetched data and the product, and are skipped, and say so, when they are missing. They verify what is the atlas's own: that every file it reads is in its dataset's manifest; that the product's labels are numbered once, each a structure of the catalogue, on the template's grid; that Schaefer's parcels lie on the side their names say and Julich-Brain's areas fall in their lobes, with soft colours; that the electrodes fall where they should on the template's scalp; that pasted names are read (the atlas's, FreeSurfer's, Latin and other names, with or without a side); and that the tours name structures that exist. The science (the tree, the sides of every key, the knowledge base) is checked with the neuroscience domain.

## 5. Structure

```
applications/brain-atlas/
├── index.html, styles.css         the page and its style
├── serve.py                       the local server (port 8810)
├── src/app.js                     the application: loading, painting, selection, the card, the plate, the cut, the modes, the address
├── src/scene/viewer.js            the round view: three.js, materials, glass, cuts and their face, picking, markers, scaling
├── src/sections/sections.js       the sections: the template's image washed in the labels' colours, at any voxel size
├── src/data/anatomy.js            the anatomy product and the template, from the cache
├── src/data/parcellations.js      Julich-Brain and Schaefer on the template's grid, with their colours
├── src/data/mouse.js              the Allen mouse brain, reoriented to RAS, with a compact label index
├── src/modes/atlases.js, mouse.js, data.js, head.js   the modes
├── src/ui/card.js, plate.js, tree.js, controls.js, tours.js, data.js   the card, the plate of names, the tree, the toolbar and search, the tours, pasted data
├── src/render/colormaps.js        colour maps (as the Earth Atlas)
├── src/workers/surface.worker.js  surfaces by surface nets, off the main thread (the scalp, the mouse)
└── tests/run-tests.mjs
```

## 6. Open Questions

1. **Surfaces of the cortex** with its layers (the pial and white surfaces of fsaverage, or BigBrain's), and parcellations defined on them.
2. **Probability maps** of Julich-Brain, so that a point is assigned with its uncertainty.
3. **The rat,** with the Waxholm Space atlas.
4. **Digitised electrode positions** and a measured head, in place of the template's provisional landmarks; a forward model, so that scalp maps lead to their sources.
5. **The fibre bundles** as streamlines from tractography, in place of the territories of a maximum-probability map.
6. **A colour of its own** for the atlas, within the warm family; it shares the orchid palette for now.
