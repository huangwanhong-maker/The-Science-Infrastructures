# 0025. An Umbrella of Laboratories for All the Sciences

- **Status:** Accepted. Written as the umbrella repository's 0001, and renumbered 0025 when the first laboratory's history was joined to it. The proposal of §2.4, to raise the shared libraries to this repository, awaits the user's agreement.
- **Date:** 2026-10-22
- **Decided by:** the user moved the Integrative Astronomy, Geography, Oceanography Laboratory into a folder named *The Science*. They asked that the folder become a repository for "an integrative laboratory, study, research and application infrastructure for all science domains", with the laboratory as a submodule. Claude wrote this record.

## 1. Context

The first laboratory began as a programme for astronomy, geography, oceanography and meteorology. Its own design anticipates growth:
- shared libraries for computation (calx), mathematical structures (morphe) and physics (physis), each in its own repository (its ADR 0008);
- fields organised so that a new one is added, not edited in (its ADR 0010);
- future fields, such as geophysics, already proposed.

The user now intends the programme to cover all the sciences. A single repository for all of them would be too large to clone, review or change safely. A set of unrelated repositories would lose the conventions and the map that hold them together.

## 2. Decision

1. **An umbrella repository.** *The Science Infrastructures* (`https://github.com/huangwanhong-maker/The-Science-Infrastructures`) holds what concerns all the sciences, and nothing particular to one:
   - the map of the sciences and the laboratories;
   - the conventions they share;
   - the decisions that cross laboratories;
   - the roadmap.
2. **Laboratories as submodules.** Each laboratory covers one domain or a group of closely connected domains. It lives in its own repository and is attached here as a git submodule. The first is the Integrative Astronomy, Geography, Oceanography Laboratory, at the path the user gave it, `astronomy and geometry/`.
3. **Shared conventions.**
   - The first laboratory's conventions apply to every repository of the infrastructure: SI units and explicit frames, the writing style, the status vocabulary, decision records, idea capture, the MIT licence for software, and datasets that keep their own licences.
   - A change to them is made in all repositories together.
   - A laboratory may add conventions of its own domain.
4. **The shared libraries, proposed to move here.** calx, morphe and physis are general, and every laboratory should be able to use them. Today they are submodules of the first laboratory. Moving them to this repository's `lib/` would:
   - make them visibly shared;
   - let a second laboratory use them without depending on the first;
   - keep one version of each across the infrastructure.

   This restructures the first laboratory (its ADR 0008), so it is proposed here and awaits the user's agreement.
5. **Dependencies.** Libraries ← laboratories ← bridges ← this repository. A laboratory uses another only through shared libraries and published products (data products, tools), never its internals. No laboratory depends on this repository.
6. **Numbering.** This repository's decisions are numbered from 0001. They are independent of each laboratory's and library's.

## 3. Alternatives Considered

- **One repository for all the sciences.** Simplest to search. But it would be too large to clone, to review and to change safely, and each laboratory's history, issues and licence questions would mix.
- **Growing the first laboratory into all the sciences.** Its name, structure and conventions are those of the geosciences and astronomy. Life or social sciences would sit there as guests.
- **Independent repositories without an umbrella.** They are flexible, but nothing would hold the map, the shared conventions or the cross-laboratory decisions, and the programme would drift apart.

## 4. Consequences

- **A new domain is a new laboratory,** and adding one touches only this repository's README, CLAUDE.md and roadmap.
- **Cloning with `--recurse-submodules`** brings every laboratory and library.
- **Commits go in dependency order:** libraries, then laboratories, then this repository's submodule pointers. Pushing remains the user's decision.
- **Links from this repository into a laboratory** use absolute GitHub URLs.
- **The first laboratory's documents** now say that it is part of this infrastructure.
- *Caveat.*
  - The first laboratory's path contains spaces and names only part of its scope; renaming it is an open decision.
  - Until the libraries move, a second laboratory would have to take them as submodules of its own, which risks diverging versions.
