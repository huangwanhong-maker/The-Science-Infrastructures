# 0072. A Woman's Body, from the Human Reference Atlas

- **Status:** Accepted, provisionally
- **Date:** 2026-11-04
- **Amends:** [ADR 0066](0066-the-anatomy-domain-and-the-body-atlas.md) (the anatomy domain's body), [ADR 0070](0070-structures-reconstructed-from-images.md) (Fabrica)
- **Decided by:** the user: "The body add feminine one. Also in gallery", and, for the gallery, that pictures should not show sensitive parts of either body. Claude chose the source, the way of reading it and its checks, and wrote this record.

## 1. Context

Fabrica showed one body, BodyParts3D's adult man. No open dataset holds a woman's whole body at BodyParts3D's completeness: the commercial anatomies are not openly licensed, and Z-Anatomy, the open atlas, is BodyParts3D's man. The Human Reference Atlas of HuBMAP (CC BY 4.0) places its 3D reference organs, modelled from the Visible Human Female, in one frame: the 3D Reference Organ Set for United, Female, v1.8, 911 meshes in ten system layers with each structure's UBERON identifier. It holds her skin and breasts, brain, spinal cord and eyes, heart and 110 named vessels, lungs, digestive, urinary and reproductive organs, spleen, thymus and a lymph node, her vertebral column, sternum, sacrum, coccyx and pelvis, and both knees with the rectus femoris. It does not hold the skull, ribs, arms or feet, nor most muscles.

## 2. Decision

1. **The data.** The registry's `hra-female` fetches the united GLB (319 MB) and its graph. `data/products/body_female.py` writes the product `body-hra-f` in body-bp3d's layout. Each mesh is carried by its nodes' transforms from the model's metres (y up, z anterior, +x the subject's left) into LPS millimetres, the soles at z = 0. That is a proper rotation, so triangles stay outward. Every node is a concept, so the node tree is the part-of tree; there is no is-a tree. Each concept keeps its UBERON identifier in a table of its own.
2. **The domain.** The anatomy domain reads her with the same reader. Her systems come from an ordered table of her own rules (`HRA_SYSTEM_RULES`), not FMA's: an element takes the system of the layer it lies in by part-of, and where a layer holds several systems, its name decides.
   - Circulatory layer: veins and arteries by their names; the rest is the heart.
   - Integumentary layer: the breasts' parts go with the genital system, as accessory organs of reproduction; the rest is the skin.
   Her vertebral levels come from her own vertebrae, C1 to L6. She has six lumbar vertebrae, a variant the model keeps.
3. **Fabrica.** A key in the header chooses the man or the woman; the page reloads with `?body=body-hra-f`. The engine keeps one state per body, reading and measuring each when first asked: her catalogue takes 41 s. Her parts are browsed by the Atlas's layers. The tours walk the man's FMA concepts and are hidden for her. The card credits the Atlas and names her structures as the Atlas does.
4. **The gallery.** No picture shows a sensitive part. The genital system is hidden in every picture of either body. For the woman it includes the breasts, and her skin is hidden too. In the Sections mode, whose round view looks up at the trunk from below, both are hidden as well. The man's midsagittal plate, whose skin and veins traced the external genitalia, was replaced by a section at L1. Fabrica itself shows the whole of either body when asked.

## 3. Alternatives Considered

- **A commercial female anatomy.** Rejected: not openly licensed.
- **Mapping her structures onto FMA.** Deferred: UBERON to FMA crosswalks are partial. Classifying by the Atlas's own layers is exact for what it holds. Her cards will gain the knowledge base's entries when an UBERON key is added to it.
- **Waiting for a complete female skeleton and musculature.** Rejected: the organs, the vessels and the trunk's skeleton serve most of Fabrica's purposes, imaging and flow first. The parts the Atlas lacks are stated in the registry and on the card.

## 4. Consequences

- **The checks.** The body check holds her when her product is built:
  - 911 elements and 7,060,019 triangles;
  - one part-of tree under the united model;
  - a system for every element, with 55 arteries, 55 veins and 14 heart parts;
  - her height, 1,657.8 mm, as the Atlas gives it;
  - the spleen on her left, the liver on her right, and the kidneys behind the liver.
- **A plate fix.** A plate is now framed on all the pieces cut. With the skin hidden, the largest loop was one organ.
- **Next.** The reconstruction of a person's vessels (ADR 0070) may now be registered to an atlas of the matching sex.
