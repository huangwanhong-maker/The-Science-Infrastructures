# 0028. Neuroscience as a Domain, and the Brain Atlas

- **Status:** Accepted; amended by [0032](0032-the-anatomy-of-the-brain-atlas.md) and [0034](0034-biosignals-on-the-brain.md)
- **Date:** 2026-10-22
- **Decided by:** the user asked for a new domain, neuroscience, beginning with a brain atlas like the Human Brain Project's, of high quality, MIT-licensed, in a git repository of its own attached as a submodule; then said that what they love is a system that visualises the brain and maps data onto it, and that setting the brain's size and shape would help, for brain-wave data among others. Claude chose the design below and wrote this record.

## 1. Context

The domains so far are those of the Earth and the sky. CLAUDE.md §1 names the life sciences and neuroscience among the domains to come, and [0026](0026-domains-as-libraries.md) places a new domain in the order by what it may depend on, with a graduation rule for a domain that gains a repository of its own.

The Human Brain Project's atlas (EBRAINS; Amunts et al. 2020; Dickscheid et al. 2026) is organised by reference spaces, parcellations whose regions form a tree, labelled and statistical maps, data features and locations. Its data are open, under licences of their own: the MNI templates (permissive), Julich-Brain (Creative Commons with a non-commercial clause), Schaefer's parcellation (MIT), the Allen Institute's mouse atlas (its terms of use).

## 2. Decision

1. **Neuroscience is a field of height 0,** depending on no other domain, declaring calx and physis as its foundations for the C++ to come. It holds the objects of the brain: reference spaces and their template grids, volumes and their affines, parcellations and maps, and the head. Its contract is that of 0026.
2. **It is a repository of its own,** `https://github.com/huangwanhong-maker/neuroscience`, MIT-licensed, attached as the submodule `domains/neuroscience`, by the graduation rule of 0026 applied at birth at the user's request. It has its own CLAUDE.md, decision log and ideas; the shared conventions are the infrastructure's. The checks of the infrastructure run it like any domain, alone with its down-set.
3. **The first content is a brain atlas** in the sense of EBRAINS, for the human, the mouse and the rat: the catalogue of spaces (Talairach, MNI 305, MNI 152 linear, 6th generation and 2009c, Colin 27, BigBrain, fsaverage, the Allen CCF v3, Waxholm Space) with their verified grids; readers of NIfTI-1 and -2 and of NRRD from bytes; the MNI–Talairach transforms of Brett and of Lancaster et al.; region trees with labels keyed by fragment, and readers of FreeSurfer's table, BIDS tables, Schaefer's networks, the Allen ontology and the EBRAINS configuration; the Desikan–Killiany parcellation built in with FreeSurfer's labels and colours, checked against its table; labelled and statistical maps with their statistics, overlaps and assignments.
4. **The head is part of the domain:** landmarks and head frames (CTF, Neuromag), the scalp as a star-shaped surface (an ellipsoid, or cast through a head mask), the 10–20 and 10–10 electrode positions constructed by arcs (Jasper 1958; Oostenveld and Praamstra 2001), scaling of a template head to a subject's measurements and Umeyama's similarity from matched points, spherical splines over the scalp (Perrin et al. 1989) and projection onto the brain. Template landmarks are provisional, placed by rule.
5. **The Brain Atlas** (`applications/brain-atlas`, port 8810) is the third application: a template brain in three slices and in the round, drawn on the graphics card from 3D textures; a parcellation as a tree, with a chosen region lit, extracted as a surface and measured; coordinates in every frame that applies; data mapped by region (pasted tables, with built-in examples), by coordinate (pasted points in MNI or Talairach coordinates) and from electrodes (pasted channel values, interpolated over the scalp and shown on the brain); the head scaled to a subject's size, the electrodes placed again. The mouse brain is shown in the Allen framework with its ontology. It shares the orchid palette and the serif type of the other atlases.
6. **Four datasets join the data library:** `templateflow-mni152`, `schaefer-2018`, `ebrains-julich-brain` and `allen-mouse-ccf`, each with its licence and citation; the application reads their files from the cache through the domain's readers, inflating gzip in the browser, with no product built in Python.
7. **Surfaces by surface nets** join the shared visualization layer (`visualization/js/isosurface.js`): a surface at a level of any volume, closed and oriented, for any application.
8. **Conventions.** World coordinates are RAS millimetres; the views follow the neurological convention; a parcellation names its version and a map its volumes; every transcribed number cites its source and is checked against it.

## 3. Alternatives Considered

- **Neuroscience inside this repository,** like the other domains. The user asked for a repository of its own; the graduation rule allows it, and the checks keep it independent either way.
- **A client of the EBRAINS services** (siibra). It would make the domain impure and tie the application to a service; the files are open and small enough to fetch into the cache.
- **Products built in Python,** as the Earth Atlas has. The atlases are already in standard formats that the domain reads; a product would be a copy.
- **Marching cubes** for the surfaces. Surface nets need no case tables, close the surface and give smoother meshes.
- **Published electrode positions** on a template head. Their licences were unclear, and positions do not scale to a subject; the construction by arcs does, and reduces to the published definitions on a sphere, which the tests check.
- **A fourth atlas, the rat,** now. Its space and parcellation are catalogued; its volumes wait for the application to grow.

## 4. Consequences

- **The order of the domains** gains a second element at height 0: metrology and neuroscience. The diagrams in CLAUDE.md, the README and the architecture document show it.
- **Checks.** The domain has 90 checks against the headers of the published templates, TemplateFlow's descriptions, the published coefficients of the MNI–Talairach transforms, FreeSurfer's table, the Allen ontology, the EBRAINS configuration, constructed maps, and the 10–20 system on a sphere and an ellipsoid; the application 23 on the fetched data (grids, the label index, sides and anatomy coming out right); the visualization layer 6 on surfaces. Every domain still passes alone.
- **The repository** `huangwanhong-maker/neuroscience` must be created on GitHub and pushed before the submodule resolves for others; until then the submodule points at a commit that exists only here.
- **Data.** About 80 MB fetched for the two atlases; Julich-Brain and the Allen volumes stay in the cache under their own terms.
- *Caveat.*
  - The template head's landmarks are placed by rule and stand in for digitised ones.
  - The MNI–Talairach transforms fit the ICBM-152 linear and 6th-generation templates; the 2009c space is treated as the same frame, to about a millimetre.
  - Brain-wave data are shown as scalp potentials above the cortex, not as sources; the inverse problem is an idea in the domain's log.
  - A few numerics (an eigen-solver, trilinear interpolation, a polar decomposition) are written in the domain provisionally and belong in calx.
