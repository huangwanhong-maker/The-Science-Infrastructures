# foundations: The Shared Libraries

*Status: draft.* The mathematics, computation and physics on which every domain stands, written once.

## 1. Purpose

A foundation holds laws and methods that name no object in the world: a group, a quadrature rule, the heat equation, the spectral transform on a sphere. A domain holds a science's objects and their models: the planets, the ocean, the plates. The test of placement is therefore simple: content that names a thing in nature belongs to a domain (see [`domains/`](../domains/)); content that does not, and that more than one domain needs, belongs here.

Each foundation is a repository of its own, attached here as a git submodule, because it is general enough to be used outside this programme. Each has its own CLAUDE.md, decision log and ideas, and must make sense on its own.

## 2. Libraries

| Library | Scope | Depends on | Status |
|---|---|---|---|
| [calx](https://github.com/huangwanhong-maker/calx) | Computation: a retargetable stack (representation, schedule, compiler, runtime, backends) and numerical methods on floating-point numbers and arrays | — | planned |
| [morphe](https://github.com/huangwanhong-maker/morphe) | Mathematical structures and their equivalences: arithmetic, algebra, categories, topology, algebraic geometry, formal structures | — | planned |
| [physis](https://github.com/huangwanhong-maker/physis) | Geometry and theoretical physics: units, geometry, discretisation, mechanics, gravitation, thermodynamics, continuum, chemistry, radiation | calx, morphe | planned |
| inference | Estimation, filtering, data assimilation, uncertainty and learning: inference under uncertainty, common to every science. Its name is to be chosen when its first code is written. | calx, morphe | planned |

## 3. Dependencies

- calx and morphe depend on nothing in the programme; physis depends on both.
- No foundation depends on a domain, an application or this repository.
- Domains, the shared edges and applications use a foundation only through its public interface.

## 4. Working Here

- Work inside `foundations/<name>`, and read that library's CLAUDE.md first.
- Commit in the library, then commit the updated submodule pointer in this repository.
- Commit in dependency order: calx and morphe, then physis, then this repository.
- Links into a library use absolute GitHub URLs, since a relative link into a submodule does not resolve on GitHub.
