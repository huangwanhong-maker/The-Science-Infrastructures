# 0061. The Domains of Life: from the Molecule to the Ecosystem, Placed in the Order

- **Status:** Accepted, provisionally. Each domain is `planned`. Its modules, any dependencies beyond those declared, and whether it lives here or in a repository of its own are decided when its work begins.
- **Date:** 2026-10-31
- **Decided by:**
  - **The request.** The user asked Claude: "in domains add life science related domain(s). And also update domains in root README.md".
  - **What Claude did.** It took the domains from the README's own map (§4.4), placed them in the order, wrote each one's README and `module.json`, brought the map and the diagrams of the order up to date, and wrote this record.

## 1. Context

The README's map of the sciences (§4.4) lists eight domains of life. Three have code: neuroscience ([ADR 0028](0028-neuroscience-and-the-brain-atlas.md)), and genomics and evolution ([ADR 0043](0043-the-genome-infrastructure.md)). Five had no directory: cell biology (priority 3.3), ecology (3.7), and microbiology, physiology and developmental biology (tier 4). Cell biology stands on molecular science (1.6, in the user's tier 1). The map files molecular science under matter, but it holds structural biology and biochemistry (§4.5).

A domain without a directory is not in the order. No `module.json` can name it, so nothing can be placed above it. The map said what each domain stands on. It did not say what each holds against its neighbours, or which foundations it would use.

The map was also out of date. It counted "fifty-six domains, ten of them built" when it lists fifty-seven, of which seventeen are built, and it gave genomics and evolution their tiers although both had been built.

## 2. Decision

### 2.1 Six Domains Placed

Six directories, each with a README (purpose, scope by module, dependencies, used by, open questions, references) and a `module.json`, all `planned`:

| Domain | Kind | Height | May depend on | Foundations | Priority |
|---|---|---|---|---|---|
| [molecular-science](../../domains/molecular-science/) | field | 0 | — | calx, morphe, physis | 1.6 |
| [cell-biology](../../domains/cell-biology/) | field | 1 | genomics, molecular-science | calx, morphe, physis, tyche | 3.3 |
| [microbiology](../../domains/microbiology/) | field | 2 | genomics, evolution, cell-biology | calx, morphe, physis, tyche | 4 |
| [physiology](../../domains/physiology/) | field | 2 | cell-biology, neuroscience | calx, physis, tyche | 4 |
| [developmental-biology](../../domains/developmental-biology/) | field | 2 | cell-biology | calx, morphe, physis | 4 |
| [ecology](../../domains/ecology/) | field | 5 | geography, climate, evolution | calx, physis, tyche | 3.7 |

Placing a domain does not change its priority. Placement says where a domain stands in the order; priority says when its work begins.

### 2.2 Composition and Descent

The domains of life lie along two axes.

- **Composition.** A cell is made of molecules, a tissue of cells, an organism of tissues. Where the order follows composition, a domain stands on the domains of what its objects are made of, whose models explain its own, and never the reverse: molecular science, then cell biology, then physiology and developmental biology.
- **Descent.** Genomes (genomics), their history (evolution), and the communities of organisms in a place (ecology).

Microbiology is where the two axes meet for single-celled life: the cell and its lineage are one object of study, so it stands on both.

### 2.3 What Each Holds Against Its Neighbours

- **Genomics and cell biology.** A gene as a word, and the count of its transcripts, belong to genomics (`annotation`, `expression`). A gene as a node of a regulatory network, and a cell as a point in expression space, belong to cell biology.
- **physis, molecular science and cell biology.** A reaction network with mass action and its Markov chain names nothing in nature, and is physis's `chemistry`. A reaction between named molecules, with its rate and mechanism, is molecular science's. A pathway of the cell is cell biology's.
- **Molecular science and chemistry.** The molecule belongs here. Substances in bulk (phases, solutions, electrochemistry) belong to chemistry, which stands on molecular science when it is founded. The Chemistry Workbench keeps its little chemistry in its plug-in until a chemistry domain is founded, as the user chose ([ADR 0051](0051-the-chemistry-workbench.md)).
- **Neuroscience and physiology.** The brain, its maps and signals, and networks of spiking neurons belong to neuroscience. The heart, muscle and the regulation of the body belong to physiology, which declares neuroscience for the nervous control of organs.
- **Evolution and microbiology.** The systematics of plants, animals and fungi stays with evolution (README §4.5). The systematics of bacteria, archaea and viruses goes to microbiology: it needs evolution's trees and is bound up with their pangenomes. So microbiology declares evolution, beyond what the map proposed.
- **Microbiology, physiology and epidemiology.** A pathogen within its host belongs to microbiology (`infection`), the host's response to physiology (`immunity`), and spread between hosts to epidemiology. Which bridge couples them is left open.
- **Cell biology and developmental biology.** A cell's regulatory states belong to cell biology; how they are laid out in space and time in an embryo belongs to developmental biology.
- **Evolution and ecology.** Ecological time belongs to ecology and evolutionary time to evolution. Ecology declares evolution, for phylogenetic diversity and population genetics. It also stands on climate, which makes it, at height 5, the highest domain in the order.

### 2.4 Disciplines Held as Modules

README §4.5 gains these rows:
- cell biology: systems biology, synthetic biology, cell signalling and metabolic networks;
- microbiology: bacteriology, virology, the microbiome and the systematics of microbes;
- developmental biology: embryology, morphogenesis and regeneration;
- physiology: plant physiology, beside immunology, endocrinology and nutrition;
- ecology: behavioural ecology, beside biodiversity, biogeography and fisheries.

### 2.5 General Methods Go to the Foundations

Each README names the general methods its domain will need, to be written in a foundation when a module first needs them (CLAUDE.md §2.1):
- **calx:** linear programming, for flux balance; the Boys function, for Gaussian bases.
- **morphe:** parallel rewriting (L-systems) in `formal`, and local moves on 2-complexes (T1 transitions of epithelia) in `topology`. Canonical forms of graphs and the counting of orbits are already written there.
- **physis:** finite dynamical systems (Boolean networks) and delay equations in `dynamics`, and ensembles in `statistical_mechanics`, all planned there. Reaction–diffusion, Turing's condition and master equations are already written.
- **tyche:** point processes, clustering and principal components (planned there), and compositional data by Aitchison's log-ratios.

### 2.6 What Is Not Decided

- Each domain's first module, its first work, and the references its checks are held to.
- Whether a domain lives here or in a repository of its own (the graduation rule of [ADR 0026](0026-domains-as-libraries.md)).
- Further dependencies, each declared when a module needs it: ecology on physiology (primary production, the metabolic theory of ecology), and developmental biology on evolution (lineage trees, evolutionary developmental biology).
- Drug discovery (1.7), epidemiology, anatomy and clinical medicine are not placed. Each waits for its work, and epidemiology also waits for demography, on which it stands.

## 3. Alternatives Considered

- **Only the five domains of the map's branch of life,** with cell biology standing on genomics alone until molecular science was founded. Rejected. To found cell biology it would have changed the map's placement of it, and molecular science is the user's next domain after genomics (tier 1.6), with structural biology and biochemistry among its modules.
- **All of life and health,** with drug discovery, epidemiology, anatomy and clinical medicine too. Rejected for now. Epidemiology stands on demography, which is not founded, and clinical medicine stands on epidemiology. Drug discovery, a bridge over genomics and molecular science, is placed when its work begins.
- **One domain, biology, with these as its modules.** Rejected. The map holds each as a domain because each has objects and models of its own, and the order then says which may use which. A single domain would hide that physiology needs neuroscience and ecology needs climate.

## 4. Consequences

- The order gains six domains, and its first height 5 (ecology). The diagrams in README §3 and §4.6, `domains/README.md` §2, CLAUDE.md §3.2 and `docs/architecture.md` §2 show them.
- README §4.4 links every domain that has a directory, counts fifty-seven domains (seventeen built, nine more placed and planned), and marks genomics and evolution built. §4.5 gains the rows of §2.4. §6 and CLAUDE.md §5 list the new directories.
- No code changes. The checks see six more domains with no suites of their own, and the dependency checker verifies that the declared order is acyclic.
