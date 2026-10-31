# molecular-science

*Status: planned. Field, height 0.* Molecules as objects of a science: atoms and the bonds between them, the graphs and shapes of molecules, their energies and motions, the macromolecules of life in space, and the reactions by which molecules become others. Placed in the order with the other domains of life at the user's request (2026-10-31, [ADR 0061](../../development/decisions/0061-the-domains-of-life.md)); its work is tier 1.6 of the [README's order](../../README.md#46-the-order-of-implementation).

## 1. Purpose

A molecule is described at three levels, each a structure with its own equivalence.

- **A graph.** Atoms are labelled vertices and bonds labelled edges. Two descriptions name one constitution when their graphs are isomorphic; a canonical form gives each constitution one name, and the automorphism group gives its classes of equivalent atoms. Rings are a basis of the cycle space over ℤ/2. The isomers of a formula, as substitutions on a fixed skeleton, are orbits of the skeleton's symmetry group, counted by Pólya's theorem.
- **A shape.** A conformation is a point of ℝ^{3n} modulo the rigid motions SE(3). The energy is a function on that quotient, invariant under reflections too, and the conformers are its local minima. The point group is the set of isometries that carry the molecule to itself by permuting like atoms; a molecule is chiral exactly when its point group contains no improper rotation (no S_n, among them mirror planes and the inversion).
- **A wave function.** The state of the electrons, from which the energy comes. In Hückel's model the π-orbital energies of a conjugated hydrocarbon are α + λβ for the eigenvalues λ of its carbon skeleton's adjacency matrix; Hartree–Fock in a basis of Gaussians is the next step.

Proteins and nucleic acids are molecules too: chains whose sequence is a word in their residues. Genomics holds the word; this domain holds the chain in space, with its backbone, fold and contacts. Structural biology and biochemistry are this domain's modules (README §4.5). Chemistry, of substances in bulk, will stand on it; so will cell biology, which declares it now, and drug discovery.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `atoms` | The elements and their isotopes: standard atomic weights (CIAAW) and isotopic masses, electronegativities, covalent and van der Waals radii, valences | planned |
| `graphs` | Molecules as labelled graphs: formulae, valence and charge; SMILES read and written; canonical forms and classes of equivalent atoms; rings as a cycle basis; aromaticity; stereo descriptors (CIP); isomers enumerated by Pólya's theorem | planned |
| `geometry` | Molecules in space: internal coordinates; force fields (bonds, angles, torsions, non-bonded terms); conformers by minimisation; superposition by Kabsch's rotation and the RMSD; point groups; normal modes from the Hessian | planned |
| `dynamics` | Molecular dynamics: Störmer–Verlet and RATTLE for constrained bonds, thermostats and Langevin dynamics, averages over ensembles, free energies by thermodynamic integration | planned |
| `electronic` | Electronic structure: Hückel's model as the spectrum of the graph, then Hartree–Fock in Gaussian bases; orbitals, charges, dipoles | planned |
| `biopolymers` | Proteins and nucleic acids in space: residues and chains; PDB and mmCIF; backbone dihedrals and Ramachandran's map; secondary structure from hydrogen bonds (DSSP); contact maps; structural alignment and the TM-score; lattice models of folding | planned |
| `reactions` | Reactions between named molecules: mechanisms, rate laws by Arrhenius and Eyring, enzyme kinetics (Michaelis–Menten as a reduction of mass action), binding equilibria; their networks handed to physis's `chemistry` | planned |

## 3. Dependencies

| Source | What molecular science uses from it |
|---|---|
| morphe | `combinatorics` (graphs, canonical forms, automorphisms), `topology` (cycles over ℤ/2), `algebra` (group actions, Burnside's and Pólya's counting, characters) |
| physis | point groups and their character tables, `geometry` (rotations), `mechanics` (motion from a discrete action), `quantum`, `chemistry` (reaction networks, mass action and their Markov chains), `statistical_mechanics` (ensembles, when they are written) |
| calx | linear algebra and the symmetric eigenproblem, optimisation, Störmer–Verlet and RATTLE, special functions |

A general method this domain needs, one that names no molecule, is written in the foundation where it belongs. The integrals of Gaussian bases, for instance, need the Boys function, which belongs in calx's special functions.

## 4. Used By

Cell biology stands on it. The Chemistry Workbench's plug-in *Symmetry and generation* keeps its little chemistry (elements, idealised geometries, CIP for a single centre, checked against RDKit) in the plug-in, as the user chose, until a chemistry domain is founded ([ADR 0051](../../development/decisions/0051-the-chemistry-workbench.md)). Chemistry and drug discovery will stand on this domain when they are founded.

## 5. Open Questions

- Where the line falls between this domain and chemistry: the molecule and its reactions here; substances in bulk (phases, solutions, electrochemistry) there.
- Which references the checks are held to, by scripts that run them for the fixtures: RDKit for graphs and stereochemistry, PySCF for electronic structure, OpenMM for force fields and dynamics, DSSP for secondary structure.
- How a protein's word (genomics) and its chain in space (here) meet when a structure's sequence differs from its gene's.

## 6. References

- Atkins, P., de Paula, J. and Keeler, J. (2018), *Atkins' Physical Chemistry*, 11th ed., Oxford University Press.
- Leach, A. R. (2001), *Molecular Modelling: Principles and Applications*, 2nd ed., Prentice Hall.
- Szabo, A. and Ostlund, N. S. (1996), *Modern Quantum Chemistry*, Dover.
- Frenkel, D. and Smit, B. (2002), *Understanding Molecular Simulation*, 2nd ed., Academic Press.
- Pólya, G. (1937), Kombinatorische Anzahlbestimmungen für Gruppen, Graphen und chemische Verbindungen, *Acta Mathematica* 68, 145–254.
- Weininger, D. (1988), SMILES, a chemical language and information system. 1. Introduction to methodology and encoding rules, *Journal of Chemical Information and Computer Sciences* 28, 31–36.
- Kabsch, W. (1976), A solution for the best rotation to relate two sets of vectors, *Acta Crystallographica* A32, 922–923.
- Kabsch, W. and Sander, C. (1983), Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features, *Biopolymers* 22, 2577–2637.
- Zhang, Y. and Skolnick, J. (2004), Scoring function for automated assessment of protein structure template quality, *Proteins* 57, 702–710.
- Berman, H. M. et al. (2000), The Protein Data Bank, *Nucleic Acids Research* 28, 235–242.
