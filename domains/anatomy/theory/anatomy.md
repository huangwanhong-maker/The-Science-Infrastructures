# The Anatomy Domain: Frames, Parts, Measures and Sections

*Status: draft (2026-11-04). Serves: learning and use.* The mathematics of the anatomy domain, stated with its hypotheses, with what the code computes and where it is checked ([ADR 0066](../../../development/decisions/0066-the-anatomy-domain-and-the-body-atlas.md)).

## 1. Overview

A body, for this domain, is a finite set $E$ of *elements* (meshes). Each element $e$ represents a concept $c(e)$ of the Foundational Model of Anatomy (FMA). The concepts are ordered by two relations, is-a and part-of. Each element's mesh is a triangle surface $S_e \subset \mathbb{R}^3$ in a patient frame.

The domain says four things:
- where things lie: §2, the frames;
- what a thing is and what it belongs to: §3, the parts;
- how large it is and of what form: §4, the measures;
- what a plane through the body meets: §5, the sections.

## 2. Frames of the Body

**Definition 2.1 (the anatomical position and directions).** The body stands upright, facing forward, its upper limbs at its sides. The six directions are named from the body: left and right are the subject's; anterior faces the way the face looks; superior points toward the head. *Remark.* Terminologia Anatomica (FIPAT 2019) and every anatomical text use these terms. They never depend on an observer, which is why a picture of the body must say where it is seen from.

**Definition 2.2 (patient frames).** A patient frame is a right-handed orthonormal frame named by the directions of its axes:
- **LPS:** $(e_x, e_y, e_z) = (\text{left}, \text{posterior}, \text{superior})$. This is DICOM's patient coordinate system (PS3.3, C.7.6.2.1.1) and BodyParts3D's.
- **RAS:** $(\text{right}, \text{anterior}, \text{superior})$. This is NIfTI's, MNI's and the Brain Atlas's.

**Proposition 2.3.** Both frames are right-handed, and $x_{\mathrm{RAS}} = R\,x_{\mathrm{LPS}}$ with $R = \operatorname{diag}(-1, -1, 1) \in SO(3)$, the turn by $\pi$ about the superior axis; $R^2 = I$. *Proof.* $e_L \times e_P = e_S$ in the body (facing $-e_P$ with $e_S$ up, the left is $e_S \times (-e_P) = e_L$). And $R$ has determinant $+1$. ∎ (`frames.hpp`; `check_frames.cpp`.)

**Definition 2.4 (the display of a section).** A display is a basis $(u, v, w)$ of the frame, with $u$ to the image's right, $v$ up and $w = u \times v$ toward the viewer. The **radiological** convention, set by Hounsfield's scanners, shows a transverse section from the feet ($w = -e_S$, $u = e_L$, $v = -e_P$), so the subject's left lies on the image's right. It shows a coronal section from the front ($w = -e_P$). The **neurological** convention shows them from above and from behind. A sagittal section is seen from the subject's left ($w = e_L$, $u = e_P$): anterior on the image's left.

## 3. Parts and Their Systems

**Definition 3.1.** Let $C$ be the set of concepts. A *forest* on $C$ is a partial map $p : C \rightharpoonup C$, the parent, with no cycles. Its reflexive transitive closure $\le$ is a partial order whose principal up-sets are chains: $x \le a$ iff $a = p^k(x)$ for some $k \ge 0$. The *meet* of two concepts in one tree is their least common ancestor. The *down-set* of $a$ is $\{x : x \le a\}$.

**Fact 3.2 (FMA's part-of is a graph; BodyParts3D's tree is a spanning tree).**
- FMA's part-of relation is a directed acyclic graph: a muscle of the body wall is part of the body wall and part of its region (Rosse and Mejino 2003, §3).
- BodyParts3D publishes a forest $p$ for browsing, one parent for each concept.
- Separately, it lists for each compound $a$ the elements $L(a) \subseteq E$ composing it, according to the whole graph.

The checks establish, for both relations and every compound $a$ (`check_body.cpp`):
$$\{e : c(e) \le a\} \subseteq L(a).$$
Equality holds for 2,904 of 2,905 is-a compounds and 1,063 of 1,368 part-of compounds. The elements of a concept are therefore read as
$$E(a) = \{e : c(e) \le a\} \cup \{e : e \in L(b) \text{ for some } b \le a\},$$
that is, an element lies below $a$ when its own concept does in the tree, or when a compound it is listed under does (`Atlas::element_below`).

**3.3 The systems.** A system is assigned by an ordered list of rules $(s_i, R_i, a_i)$, $i = 1, \ldots, n$: "below $a_i$ in relation $R_i$". The system of an element is $s_{i^*}$, where $i^*$ is the least $i$ with $e$ below $a_i$ in $R_i$. The order encodes anatomical judgement, so it is the domain's to state:
1. **Regional exceptions first.** Teeth are digestive, and the larynx's cartilages and ligaments respiratory. Otherwise they would fall under cartilage and ligament as skeleton.
2. **Tissue kinds by is-a.** Vessels, then nerves and the central nervous system, then bone, cartilage and ligament, then muscle and tendon. The vessels of the lung are vessels, though the lung is respiratory by part-of.
3. **Organ systems by part-of.** The heart, eye, ear, respiratory, digestive, urinary, genital, lymphoid, endocrine and integument, with a few kinds the trees leave outside them (the testis, the salivary glands).
4. **Last, the nervous system by part-of.**

On BodyParts3D 4.0 every element falls under a rule (`check_body.cpp`). Seventeen named parts are checked to fall where anatomy puts them.

**3.4 Sides.** A name's side is read from the words *left* and *right*. Far from the midline ($|x| > 100$ mm, the limbs) the side named is the body's. Of 429 such parts, all lie on their named side in LPS but one pair, which the database mislabels. Nearer the midline a name's side may be a chamber's or a branch's (the domain's README, §5).

## 4. Measures

**4.1 Topology.** A mesh is welded at $10^{-4}$ mm and split into components. Each component is classified by morphe (decision 0010): whether it is a surface, closed, its genus and its boundary circles, and whether the triangles' own orientation is consistent.

**Proposition 4.2 (the solid).** Let $S$ be a closed surface, consistently oriented outward, bounding a solid $\Omega$. For any reference point $o$,
$$|\Omega| = \sum_{t=(a,b,c)} \tfrac16 \det[a-o,\ b-o,\ c-o],$$
and likewise the first and second moments, by Tonon's formulas for the tetrahedra $(o, a, b, c)$.
- *Proof.* The divergence theorem applied to $x \mapsto (x - o)/3$ and its polynomial analogues. Each triangle's flux is the signed tetrahedron's moment.
- *Remark.* If the orientation is inward, every term changes sign; physis refuses a non-positive volume, and the domain reads the mesh reversed (`inward`).
- *Computation.* physis's `solid_polyhedron`, decision 0019. The principal axes are the eigenvectors of the inertia of unit density. The first axis, of least moment, is the axis along which the part is longest. The *extents* are the lengths of the part's projections on the axes.

**4.3 The database's volumes.** BodyParts3D publishes a volume for each element. The volumes computed here agree with them at a median ratio of 1.0004. The parts over 100 cm³ agree to 1 % in 94 % of cases. Thin parts scatter, by factors up to 40 for some small vessels. The published bounds agree with the meshes' to 0.5 mm, so the difference lies in the volumes alone. That is what one expects if the database computed its volumes on its full-resolution meshes, of which the distributed files keep 1 % of the polygons. *Caveat.* This is an inference, not documented by the database. The volumes here are exact for the meshes distributed.

**4.4 Compounds.** A compound's solid is the sum of its elements' (physis's `combine`, the parallel axis theorem). Overlapping elements are counted twice. When not all elements bound solids, the compound has no volume; the sum over those that do is reported as a lower bound.

## 5. Sections

**5.1** The section of $S$ by the plane $\{x : w \cdot x = h\}$ is computed by physis's `plane_section` with normal $w$, the display's. Each piece's loops then run counter-clockwise in $(u, v)$ round the solid, and its area $A(h)$ is the sum of their shoelace areas.

**Proposition 5.2 (Cavalieri).** For a solid bounded by a closed surface, $\int A(h)\,dh = |\Omega|$. *Proof.* Fubini's theorem on the indicator of $\Omega$. ∎

The checks hold the computed sections to the computed volume, on the left kidney, to $5.6 \times 10^{-6}$ by the midpoint rule on 400 planes. This joins the two computations (the divergence theorem and the slicing), which share no code beyond the mesh.

**5.3 Vertebral levels.** The centres $z_k$ of the vertebrae's meshes, C1 to the sacrum, descend strictly (checked). A height $z$ is named by the vertebra $k$ whose interval holds it. The intervals are bounded by the midpoints $(z_{k-1} + z_k)/2$, and at the ends by as far again as half the gap to the neighbour. *Remark.* This is the vertebra's centre, not its body's plate. Atlases name planes by the disc or the body they cross, and the two differ by up to half a vertebra.

**5.4 The calibre of a part.** Sections normal to a part's long axis give $A(s)$ along it. The diameter of the circle of that area, $2\sqrt{A(s)/\pi}$, is its calibre. For a vessel it is the first step toward the hemodynamics planned next. The right great saphenous vein's calibre runs from about 3 mm at the ankle to about 5 mm near the groin.

## References

- Cavalieri, B. (1635), *Geometria indivisibilibus continuorum*; Fubini's theorem in any text of measure theory.
- DICOM PS3.3 (2024), C.7.6.2.1.1.
- FIPAT (2019), *Terminologia Anatomica*, 2nd ed.
- Minetto, R., Volpato, N., Stolfi, J., Gregori, R. M. M. H. and da Silva, M. V. G. (2017), An optimal algorithm for 3D triangle mesh slicing, *Computer-Aided Design* 92, 1–10.
- Mitsuhashi, N. et al. (2009), BodyParts3D, *Nucleic Acids Research* 37, D782–D785.
- Rosse, C. and Mejino, J. L. V. (2003), *Journal of Biomedical Informatics* 36, 478–500.
- Tonon, F. (2004), Explicit exact formulas for the 3-D tetrahedron inertia tensor in terms of its vertex coordinates, *Journal of Mathematics and Statistics* 1, 8–11.
