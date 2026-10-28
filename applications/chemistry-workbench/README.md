# Chemistry Workbench

*Status: draft (2026-10-29). Serves: research, learning.* A workbench for chemistry, built of plug-ins ([ADR 0051](../../development/decisions/0051-the-chemistry-workbench.md)).

## 1. Purpose

Chemistry joins the infrastructure through its mathematics first: symmetry, graphs and their rewriting, reaction networks, and the field theories that generate complexes. The workbench is a place where each such study is a plug-in. A plug-in states what it computes, on the foundations, and shows it in its own views.

Its first plug-in, **Symmetry and generation**, executes the user's manuscript *An Inquiry into Chemical Structure Generation through Symmetry Breaking* (v1.4) and verifies it claim by claim, with a note that proves what holds and states what must change.

## 2. Use

```
python applications/chemistry-workbench/serve.py            # http://127.0.0.1:8950/
```

The engines are built from the page on first use, or beforehand:

```
python data/acquisition/recipe.py applications/chemistry-workbench/recipe.json
# or: python applications/chemistry-workbench/engine.py build
```

g++ 11 or later is needed, natively or in WSL. Compiling takes a minute or two. The build takes the WSL lock the sessions share (`%TEMP%/science-wsl.lock`).

## 3. Plug-ins

| Plug-in | Views | Computes on |
|---|---|---|
| [Symmetry and generation](plugins/symmetry-generation/) | Grammar, Stereocentre, Point groups, Selection, Complexes, Field | morphe: graph rewriting, homology, characters. physis: point groups, stranded field theories, reaction kinetics. calx: integrators, random streams |

**The contract.** A plug-in is a directory `plugins/<name>/` with:
- `plugin.json`: its title, summary, status, activity, `views` (`id`, `title`, `module`, `about`) and, if it computes, `engine`;
- its views, ES modules exporting `render(root, ctx)`. `ctx` gives:
  - `ask(op, body)`, one request to the plug-in's engine;
  - the drawing kit (`chart`, `graph`, `complex`, and `/src/draw/molecule.js` by import);
  - `dom`, `base` (the plug-in's URL), `onDispose`, `go`;
- for C++, `engine/main.cpp`. It answers JSON lines (`{"op": …}` → `{"ok": true, "result": …}`) with `--serve`, as the server keeps it.

A plug-in joins by being there: the rail is built from the manifests, and nothing in the shell names a plug-in.

## 4. Layout

| Path | Content |
|---|---|
| `index.html`, `styles.css` | the page; wisteria violet with an apricot minor, on lilac paper by day and aubergine by night (provisional, until the user chooses) |
| `src/app.js` | the shell: rail, routing, theme, building engines from the page |
| `src/api.js`, `src/dom.js` | the server's API; elements, SVG and the hover note |
| `src/draw/` | `chart.js` (lines, bars), `graph.js` (labelled graphs: atoms by element, nonterminals as dashed boxes), `complex.js` (2-complexes and integer matrices), `molecule.js` (three.js 0.170.0: molecules, symmetry elements, a dipole) |
| `engine.py` | builds each plug-in's engine into `data/cache/engines/chemistry/<name>/`; `Worker` keeps it alive |
| `serve.py` | the server: `GET /api/plugins`, `POST /api/<plugin>/<op>`, the recipe at `/setup` |
| `recipe.json` | the build of the engines |
| `tests/` | `run-tests.mjs`; `make_fixtures.py` and `fixtures/` (RDKit's CIP descriptors) |

## 5. Checks

`node applications/chemistry-workbench/tests/run-tests.mjs` (22 checks) covers:
- the manifests, the recipe, and every page module parsing;
- when the engine is built:
  - the 55 checks of the manuscript's claims;
  - CHBrClF's CIP descriptors against RDKit;
  - the catalogue's 15 point groups, and the character tables;
  - the derivation counts of four glycines;
  - the complexes' homology;
  - Wick sums against brute-force moments;
  - the field's characters;
  - Frank's network.

The mathematics is checked in the foundations: morphe's `check_graph_rewriting`, physis's `check_point_groups`, `check_stranded` and `check_kinetics`.

## 6. Dependencies

The foundations morphe, physis and calx (C++, header-only), and the API's C++ JSON. In the page: three.js 0.170.0 from jsDelivr, and Source Serif 4. There is no domain yet: the plug-in's chemistry stays in it until a chemistry domain is founded.
