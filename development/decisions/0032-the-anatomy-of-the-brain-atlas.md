# 0032. The Anatomy of the Brain Atlas: Surfaces, a Tree of Structures and What Is Known of Each

- **Status:** Accepted; amended by [0034](0034-biosignals-on-the-brain.md)
- **Date:** 2026-10-23
- **Amends:** [0028](0028-neuroscience-and-the-brain-atlas.md), §2.5 and §2.6
- **Decided by:** the user found the first Brain Atlas rough, its interface poor, and asked that the anatomical structures be comprehensive and well rendered, and that clicking a part show its anatomy; they asked for high quality, creativity, and time taken. Claude chose the design below and wrote this record.

## 1. Context

The first Brain Atlas (0028) drew slices of the template and, in the round, a brain surface and one region at a time, extracted in the page from label volumes. It showed parcellations well, but not the brain as an anatomist knows it: the gyri, the deep nuclei, the ventricles and the fibre bundles together, each named, each with something to say. Surfaces extracted in the page at the grid's resolution were coarse, and a structure had no description beyond its name.

Open atlases on the template's grid now cover the whole brain at the level of structures: FreeSurfer's segmentation of the template, the Harvard–Oxford gyri, MASSP 2.0's subcortical nuclei from 7 T MRI (Bazin et al. 2020), the thalamic nuclei of Najdenovska et al. (2018), Julich-Brain's cerebellar nuclei, and the deep white-matter bundles of EBRAINS (Guevara et al. 2012).

## 2. Decision

1. **An anatomy product built in Python** (`data/products/brain_anatomy.py`): surfaces of the cortex, the white matter and the cerebellum from the template's tissue maps, the cortex's vertices labelled by gyrus; 130 meshes of nuclei, ventricles and bundles from the atlases above, extracted by surface nets from blurred and upsampled masks and smoothed by Taubin's method; each vertex shaded by a baked obscurance and, on the large surfaces, by sulcal depth; in three packs, read on demand; with the labels of all sources joined on the template's grid and the contacts between them. This amends 0028 §2.6: the surfaces are too costly to make well in the page, and the product is a derived work of the fetched files, rebuilt from them, never committed.
2. **The brain as a tree of structures,** in the domain (`brain.js`): 152 structures from the brain to the gyri and nuclei, each with the keys of its labels and meshes by side, a colour (the gyri in hues of their lobes), a short name, and ten functional systems that cut across the tree. The sides of every key are checked by the centroids of its voxels.
3. **What is known of each structure,** in the domain (`knowledge.js`): its Latin name after *Terminologia Neuroanatomica*, a summary, where it lies, what it does, what it connects, its blood supply (after Tatu et al.), its development and its clinical significance, with references, each paper checked against Crossref. A division of a gyrus inherits what it does not say from the gyrus. Where an atlas's region differs from the classical structure, the entry says how.
4. **The application is rebuilt around them** (amending 0028 §2.5): a viewer in the round with glass, cuts whose face shows the section, explosion, and names set out as on an engraved plate; a card for every structure; sections washed in the structures' colours; guided tours. Five modes: Anatomy, Sections, Atlases (Julich-Brain and Schaefer on the cortex, each area placed in the gross anatomy; the Allen mouse brain), Data (values by structure, by FreeSurfer's names or an atlas's regions; points with what lies there) and Head (the scalp cast from the template, the 10–20 and 10–10 systems, scaling, scalp maps carried onto the cortex, the cortex beneath each electrode).
5. **Four datasets join the data library:** `massp-2`, `harvard-oxford`, `mial-thalamic-nuclei` and `ebrains-deep-white-matter`, with their licences: CC BY 4.0, the FSL licence (non-commercial), CC BY 4.0 and CC BY-NC-SA 4.0.

## 3. Alternatives Considered

- **Surfaces made in the page,** as before. They are coarse at 1 mm and slow for 130 structures; the product makes them once, smoothly.
- **Published meshes** of the brain (for instance from BodyParts3D or a commercial model). Their licences and their spaces vary; meshes made from the atlases on the template's grid agree with the sections and the coordinates exactly.
- **Descriptions fetched from an encyclopaedia.** They would make the domain impure and their accuracy would be unchecked; the knowledge base is written from the standard references and checked by the tests for coverage, references and typography.
- **The Desikan–Killiany gyri** in place of Harvard–Oxford's. FreeSurfer's parcellation of the template's own surface is not distributed with TemplateFlow; Harvard–Oxford's lateralised maximum-probability map is, on the same grid.

## 4. Consequences

- The atlas shows the brain as a whole and in its parts, and says what each part is; the cards are a learning text in their own right.
- The product must be built once (`python data/products/brain_anatomy.py`, about a minute) after the datasets are fetched.
- The Harvard–Oxford atlas and the EBRAINS data carry non-commercial licences; they are fetched, not redistributed, as before.
- The knowledge base is long prose in a domain; it is for learning, not for clinical decisions, and the About dialog says so.
