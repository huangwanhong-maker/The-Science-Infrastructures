# Brain Atlas

*Status: draft (2026-10-22). Serves: use and learning. Draws on: neuroscience.* Decision record: [ADR 0028](../../development/decisions/0028-neuroscience-and-the-brain-atlas.md).

The Brain Atlas shows a template brain in three slices and in the round, divides it by a parcellation whose regions form a tree, and maps data onto it: values by region, points by coordinate, and brain-wave channels from electrodes on the scalp. It shows the human brain in the space of the Human Brain Project's atlas (MNI 152 ICBM 2009c) with Julich-Brain's cytoarchitectonic areas and Schaefer's functional parcels, and the mouse brain in the Allen Institute's Common Coordinate Framework with its ontology. The brain can be scaled to a subject's head, and the electrodes of the 10–20 system follow.

## 1. Concepts

- **A point and its slices.** Three orthogonal slices pass through one point, the axial, coronal and sagittal; the fourth view shows the same planes in the round, with the brain's surface. Click or drag a slice to move the point, scroll to step through slices, use the arrow keys, or click the round view. Views follow the neurological convention: the subject's left on the viewer's left; the sagittal view looks from the left.
- **A parcellation as a tree.** Regions contain regions: the cerebral cortex its lobes, an area its two hemispheric parts, a network its parcels, the isocortex its areas and their layers. The tree is browsed and searched; a chosen region is lit in every view, extracted as a surface, and measured (volume, centroid, extent). The point card names the region under the point and its ancestors.
- **Coordinates in every frame.** MNI millimetres, voxels, and Talairach coordinates by two published approximations (Brett; Lancaster et al.) for the human; RAS millimetres and the framework's own micrometres for the mouse.
- **Data by region.** A pasted table of regions and values colours the regions by value; a value given for a parent colours everything beneath it. Regions are matched by name, acronym, label or FreeSurfer key, with the hemisphere written any way. Built-in examples: the volume of each region, its hemisphere, its depth in the tree, its group.
- **Data by coordinate.** Pasted points, in MNI or Talairach coordinates (or RAS and CCF for the mouse), are drawn as spheres in the round and as rings in the slices they are near, and each is assigned to the region it falls in.
- **The head and its electrodes.** The template's scalp carries the electrodes of the 10–20 system (or the 65 positions of the 10–10 system above the landmark level), placed by the arcs between the nasion, inion and preauricular points as Jasper (1958) and Oostenveld and Praamstra (2001) define them. A subject's head length, breadth and height scale the brain, and the electrodes are placed again on the scaled scalp. Pasted channel values become a map over the scalp by spherical splines (Perrin et al. 1989), shown on the scalp and on the brain beneath.
- **Everything computed here.** The slices are drawn on the graphics card from the template and the label volume; the surfaces are extracted by surface nets in a worker; every other number comes from the neuroscience domain. The page loads only files of the data library's cache, three.js and a font.

## 2. Running

From the repository root, fetch the data (about 80 MB):

```
python data/acquisition/fetch.py templateflow-mni152 schaefer-2018 ebrains-julich-brain allen-mouse-ccf
python applications/brain-atlas/serve.py --open
```

The page is at http://127.0.0.1:8810/. The address bar keeps the atlas, the parcellation, the layout, the point and the chosen region, so a link reproduces a view.

- **Layouts:** four views; the round view alone; the three slices.
- **Gestures:** click or drag a slice to set the point; scroll to step (shift: 5 mm); ctrl-scroll to zoom; shift-drag to pan; double-click to reset. Drag the round view to turn it, scroll to approach, shift-drag to pan, click to pick a point on a plane, a surface or a marker. Arrow keys move the point in the plane; Page Up and Down across it; Escape deselects.

## 3. Data

| Data | Source | Dataset | Licence |
|---|---|---|---|
| MNI 152 ICBM 2009c nonlinear asymmetric, 1 mm, brain and head, with its mask | McConnell Brain Imaging Centre (Fonov et al. 2011), via TemplateFlow | `templateflow-mni152` | MNI: permission to use, copy, modify and distribute |
| Julich-Brain v3.1, the maximum probability maps of 207 areas, one volume per hemisphere, with the EBRAINS configuration | Forschungszentrum Jülich, EBRAINS (Amunts et al. 2020) | `ebrains-julich-brain` | EBRAINS terms, non-commercial; fetched, not redistributed |
| Schaefer 2018, 100 and 400 parcels in 7 networks | CBIG (Schaefer et al. 2018; Yeo et al. 2011), via TemplateFlow | `schaefer-2018` | MIT |
| Allen Mouse Brain CCF v3, template and 2017 annotation at 50 µm, with the ontology | Allen Institute (Wang et al. 2020) | `allen-mouse-ccf` | Allen Institute terms of use |

The Desikan–Killiany parcellation is built into the domain (names, labels, colours and lobes, from FreeSurfer's table) and is listed for use with an aparc+aseg volume of one's own; no volume of it is fetched here.

## 4. Checks

```
node applications/brain-atlas/tests/run-tests.mjs
```

The checks need the fetched data and are skipped, and say so, when a dataset is missing. They verify what is the atlas's own: that the parcellations' volumes lie on the template's grid; that the label index joins every voxel's label to a region; that the sides come out right (every parcel named left has its centroid at negative x; Julich-Brain's left-hemisphere volume holds only left-hemisphere areas); that the networks sit where anatomy puts them (the visual network most posterior); that the mouse annotation's labels are all in the ontology and the isocortex lies at the dorsal surface; and that the manifests name every file the atlas reads. The science is checked with the neuroscience domain.

## 5. Structure

```
applications/brain-atlas/
├── index.html, styles.css       the page and its style
├── serve.py                     the local server (port 8810), mounting domains/, foundations/, visualization/ and the data cache
├── src/app.js                   the application: state, loading, the panel, the address bar
├── src/data/atlases.js          the files of each atlas, fetching and inflating, the label index on the template grid
├── src/render/stage.js          the WebGL stage: four views, slice shader, surfaces, markers, gestures
├── src/render/colormaps.js      colour maps (as the Earth Atlas)
├── src/workers/surface.worker.js  surfaces by surface nets, off the main thread
├── src/ui/tree.js, panel.js, data.js   the region tree, the panel's HTML, the parsing of pasted data
└── tests/run-tests.mjs
```

## 6. Open Questions

1. **Surfaces of the cortex** (the pial and white surfaces of fsaverage) and surface-based parcellations.
2. **Probability maps** of Julich-Brain, so that a point is assigned with its uncertainty, as the domain already computes.
3. **The rat**, with the Waxholm Space atlas.
4. **Digitised electrode positions** and a measured head, in place of the template's provisional landmarks; a forward model for the electrodes.
5. **A colour of its own** for the atlas, within the warm family; it shares the orchid palette for now.
