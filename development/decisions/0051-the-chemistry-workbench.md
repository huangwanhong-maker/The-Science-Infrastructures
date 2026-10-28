# 0051. The Chemistry Workbench, and the Manuscript on Chemical Structure Generation Made Exact

- **Status:** Accepted, provisionally (the workbench's name and colour are provisional; the chemistry its plug-in holds moves to a chemistry domain when one is founded)
- **Date:** 2026-10-29
- **Decided by:**
  - **The request.** The user asked Claude to "implement chemistry lab": first only the algorithm of their manuscript *An Inquiry into Chemical Structure Generation through Symmetry Breaking: Two Formalisms Drawing on Group-Equivariant Grammar and Group Field Theory* (v1.4, 2026-10-28), as a plug-in. They asked to "serious formal verify the paper's formalism are corrected", to implement it in code, and to try to visualize examples. They added: "if anything you think not make sense, or error, please tell me, we need revise paper".
  - **The user's choices.** Asked by Claude (CLAUDE.md §1 retires "laboratory"), they named it the **Chemistry Workbench**. Asked whether to found a chemistry domain now, they chose to keep the plug-in's small chemistry inside it for now.
  - **What Claude did.** It wrote the foundations' new methods, the workbench, its plug-in, the verification note and this record.

## 1. Context

The manuscript proposes two symmetry-first formalisms for chemical structure generation:
1. a typed chemical graph grammar with a group acting on rules, states and contexts, and its stochastic dynamics;
2. a group-field-inspired generator of labelled 2-complexes, read chemically through a separate readout.

Around them it states three conjectures and the contrastive methodology that keeps formal, physical and ontological claims apart. Its case studies give exact fixtures: the stereocentre of CHBrClF, the six-carbon cycle, Boulatov's tetrahedral interaction over Z₂, and a two-vertex complex.

The infrastructure had no chemistry. The foundations held graphs with canonical forms, group actions, homology and the coloured and tensorial group field theories. They had no graph rewriting, no point groups, no uncoloured stranded field theory, and no reaction kinetics; physis's `chemistry` module was planned.

## 2. Decision

### 2.1 The General Methods Go to the Foundations

- **morphe, decision 0007.** Graph rewriting by double pushouts with relabelling (Habel and Plump 2002): matches and their automorphisms, derivations, derivation spaces with counted derivations, classes of graphs and rules, label permutations, closure and equivariance.
- **physis, decision 0011.**
  - Point groups of labelled points, with their character tables through morphe and the responses they allow.
  - Stranded field theories on real fields over finite groups: Boulatov's vertex, its symmetries and gauge projection, Wick pairings with their faces, 2-complexes and classes, and the Gaussian moments by brute force.
  - Reaction kinetics, the first content of `chemistry`: stoichiometry, conservation and balance, mass action, and the Markov chain solved exactly and sampled.

### 2.2 The Chemistry Workbench

`applications/chemistry-workbench/` serves on port 8950. It is a shell whose rail is built from plug-ins: `plugins/<name>/plugin.json` declares a plug-in's views, and a plug-in with C++ has `engine/main.cpp`. `engine.py` builds each engine into `data/cache/engines/chemistry/<name>/`, judging staleness by the headers it includes followed transitively, and takes the sessions' WSL lock while compiling. The server keeps each engine alive as a JSON-lines worker. A plug-in joins by being a directory; nothing in the shell names it.

Its colour is provisional, until the user chooses: wisteria violet (#8E7CC3) with an apricot minor (#D9915E), on lilac paper by day and aubergine by night. Type is serif and lines are hairlines, as in every application.

### 2.3 The Plug-in *Symmetry and Generation*

`plugins/symmetry-generation/` executes the manuscript on the foundations and verifies it.

**Its engine** runs the manuscript's examples:
- the grammars: ethanol, glycine backbones, the stereocentre, the acid–base edit;
- the idealized molecules with their point groups;
- Frank's autocatalysis as mass action and as an exact Markov chain;
- the cycle, its filled complex, and the two-vertex complex;
- the tetrahedral field with its characters and order parameters;
- the Wick expansion to order 3;
- the checks of the manuscript's claims: 55 computations behind 34 claims.

**Its six views** are Grammar, Stereocentre, Point groups, Selection, Complexes and Field.

**Its note**, `theory/formal-verification.md`, proves what holds and states what must change.

**Its chemistry is provisional:** nine elements with covalent radii, idealized geometries, and the CIP rule for a centre with four single-atom ligands. It stays here by the user's choice until a chemistry domain is founded.

### 2.4 What the Verification Found

Four errors, seven amendments and three statements that hold more strongly than written (note §7). The four errors:
1. **Definition 4.2** cannot change a kept atom's charge. Adopt relabelling.
2. **For abelian G**, the symmetry §5.3 breaks is the GFT gauge redundancy §5.2 warns against. The gauge projection annihilates every odd order parameter.
3. **§5.3's order parameter** η lies along a flat direction of the tetrahedral interaction (V[t(−1)^a] ≡ 0). The odd character (−1)^b gives S = 4μη² + 64λη⁴, the Landau form derived rather than postulated.
4. **The two-vertex complex** is generated by the Wick expansion of the Gaussian integral, not by a chosen field configuration.

Most claims hold. Among those confirmed:
- every number of App. L;
- the generator covariance and the readout covariance;
- the factorization criterion and the data-processing bound;
- the cubic normal form and Frank's equation;
- τ = +1 is (R), proved and confirmed by RDKit;
- the Feynman expansion ⟨Vⁿ⟩ = Σ |G|^F, exactly, to order 3.

## 3. Alternatives Considered

- **Write the plug-in as one program holding all its mathematics.** This would have repeated graphs, groups, homology and field theory that the foundations hold or should hold, against the rule that a work writes no general method.
- **Found a chemistry domain now.** This is the user's decision (CLAUDE.md §6.5). They chose to wait.
- **Check point groups and CIP descriptors against values typed from memory.** The point groups are checked from generic orbits of the groups themselves. The CIP descriptors are checked against RDKit, by a script that makes the fixture (`tests/make_fixtures.py`).

## 4. Consequences

- The manuscript can be revised from the note's §7, and every revision rechecked by the engine.
- A chemistry domain, when founded (molecular science is in the README's first tier), will find the following in place:
  - its general methods in morphe and physis;
  - the plug-in's elements, geometries and stereochemistry, to move into the domain.
- Further plug-ins (reactions, spectra, crystals) join the workbench as directories.
