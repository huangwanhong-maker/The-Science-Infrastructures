# 0032. The C++ Plug-in Classes Completed

- **Status:** Accepted, provisionally: built at the user's request; open to revision on review. Amends [ADR 0029](0029-plugins-in-python-and-cpp.md) §2.3 and §4.
- **Date:** 2026-10-23
- **Decided by:** the user asked for everything Claude had listed as next, including the C++ classes beyond `Diagnostic`. Claude chose the design below.

## 1. Context

ADR 0029 gave C++ plug-ins only `Diagnostic`, and named `ForecastModel`, `DataAssimilation` and `Trajectories` in C++ as the next step. [ADR 0030](0030-the-api-edge.md) moved the plug-in contract out of the application into the API.

## 2. Decision

1. **The header is the API's**, `api/cpp/plugin.hpp`, as the Python classes are (`api/python/science/plugins.py`). A plug-in is compiled with `-I api/cpp`.
2. **The classes**: `Plugin`, `Diagnostic`, `ForecastModel<State>`, `DataAssimilation`, `Trajectories`, with the hooks of the JavaScript and Python classes. Three differences belong to the language: needs, outputs and what to show are members; a forecast's state is a type parameter; `Diagnostic::compute` keeps its first signature.
3. **Inputs and results.** Needs of kind `observations` reach the program as tab-separated files of the set's rows; parameters may be choices of text. The file contract gains lines for the grid and times, vectors, the start, verify, observations, sketches (polyline, marker, label, moving point), tables and notes, all documented at the head of the header.
4. **What runs where.** The program runs what calls its hooks: the loops, and the parcel step. The host replays the rest through `plugins.py`'s own classes: a forecast's verification beside persistence, and an analysis's increment and observation markers. One implementation of these serves every language, and any other language that writes the file contract.
5. **Compilers.** On Windows without a native compiler, the host compiles with `g++` in WSL and runs the Linux program with `wsl -e`, its paths translated.
6. **Examples, each checked against a closed form:** semi-Lagrangian advection (Williamson et al. 1992, case 1: a cosine bell over the poles, normalised l2 = 0.065 after one revolution on a 2° grid), parcels in a solid-body rotation, and a Cressman (1959) analysis of one observation.

## 3. Alternatives Considered

- **The verification in C++.** Faster, but a second implementation to keep equal to `plugins.py`'s; the replay keeps one, at the cost of `markers()` not being a C++ hook.

## 4. Consequences

- Plug-ins exist in all four kinds in each of JavaScript, Python and C++.
- On Windows, a C++ plug-in's first load after an idle minute costs 5–10 s while WSL's virtual machine starts.
