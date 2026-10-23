# The History Complex

*Status: draft.* An application of the history domain ([ADR 0041](../../development/decisions/0041-the-history-domain.md)). Histories are recorded as occasions, beings and ties, and built into a 2-complex whose topology is computed exactly. Its signed ties become the spin foam of a ℤ₂ gauge theory fitted to the record, from which the unattested signs are reconstructed. Serves research and learning.

## 1. What It Shows

A **canon** of classical facts, recorded by the domain's [codebook](../../domains/history/codebook.md), becomes a complex by the domain's construction ([theory](../../domains/history/theory/the-history-complex.md)): occasions are its vertices, the lives of the beings between occasions its edges, and each tie between two beings a chain of lenses, one 2-cell between each two of their meetings. The application shows it six ways:

- **The foam.** The history in spacetime, in three dimensions:
  - time rises, and each being is a thread through the occasions it attends, drifting back toward its own place between them;
  - each occasion is a knot where the threads of those present meet, and each tie's lenses are films between its two threads;
  - the lenses are coloured by discipline, by their spin in the foam, by their flips, or by the certainty of an inferred sign;
  - opened, each occasion shows its presences, the rungs of its ties and their triangles, frustrated ones in vermilion: the refined complex;
  - a cut at any date is the network of that moment, and the cut can be played upward through the history.
- **The loom.** Time runs across; each being is a thread; an occasion is a knot binding the threads of those present; a tie is a ribbon of lenses between two threads, coloured by its discipline, hatched where it is negative, marked where it flips. A being, a tie or an occasion is chosen to follow it; the disciplines are shown or hidden.
- **The moment.** The network of beings and live ties at a chosen date, the spin network of the moment: links signed, their certainty fading with the meetings since they were last attested, frustrated triangles shaded; beneath it, the frustration of the network year by year.
- **The complex.** It shows:
  - its counts and Euler characteristic, its homology and fundamental group;
  - the loops of H₁, told as the beings and occasions they pass: circles of encounters no recorded tie fills, where the canon is silent;
  - the closed shells of H₂, with the Bianchi identity on each;
  - the record's other complexes: the refined complex, the complex with its shells filled, Dowker's complex of co-presence with Atkin's Q-analysis, and the history in time with its persistent homology drawn as bars.

  The complex can be saved as JSON for the [Complex Atlas](../complex-atlas/).
- **The record.** Every assertion beside its sources, confidence and provenance, the problems the construction found, and the check of the beings against Wikidata.
- **The readings.** They show:
  - the spin foam: the couplings fitted to the record, what each label means, the cells most often in the foam, each tie's memory across its life against an isolated tie's, the power of each being over the reconstruction, and SU(2) kinematics;
  - the persistence of ties by discipline;
  - the frustrated triangles at occasions;
  - the disciplines' lenses, and the beings and occasions where they meet.

The canon chooser lists every canon of the domain (six: Athens, Florence and Rome, the Scientific Revolution, Europe 1785–1828, the European alliances and Vienna), each with a guided tour. It also offers **all canons, joined**: one record, glued along the beings, occasions and ties the canons share, where they agree who a being is (the theory note, Definition 1.3).

## 2. Running

From the repository root, build the engine once (g++ 11 or later, natively or in WSL) and serve the page:

```
python data/acquisition/recipe.py applications/history-complex/recipe.json
python applications/history-complex/serve.py --open
```

The page is at http://127.0.0.1:8870/. With `?fixtures` in the address the page answers from `tests/fixtures/`, the engine's own answers for the first canon, to see every view without the engine.

## 3. The Protocol

The engine is `engine/history.cpp`, header-only C++20 on the history domain (`domains/history/`), morphe (homology, π₁) and the API's JSON. The server runs it once per request: one JSON request on its standard input, one JSON answer on its standard output, an error as `{"error": "..."}`. The canon is passed by its path; the engine reads nothing else.

### 3.1 Requests

| Request | Answer |
|---|---|
| `{"op": "build", "canon": "<path>.jsonl"}` | the record, the complex, its invariants and the readings (§3.2) |
| `{"op": "build", "canons": ["<path>.jsonl", …]}` | the same, for the canons joined into one record |
| `{"op": "moment", "canon": "<path>.jsonl", "date": "<EDTF>"}` | the network of the moment at the last day of the date, with its reconstruction (§3.3) |
| `{"op": "power", "canon": "<path>.jsonl"}` | the power of each being, and of balance, over the reconstruction (§3.4) |

### 3.2 The Answer to `build`

Beings, occasions and ties are referred to by their index in the answer's arrays. Days are Julian day numbers (`jdn`, integers, noon-based); `year` is the decimal Gregorian year of a day, for drawing.

```
{
  "canon":   {"id", "title", "about", "codebook", "period", "status", "drafted_by", "drafted"},
  "sources": [{"id", "cite"}],
  "beings":  [{"id", "name", "type", "lived", "wikidata", "colour", "occasions": [o, ...]}],
  "occasions": [{"id", "date", "calendar", "place", "type", "label",
                 "first", "last", "year",
                 "present": [{"being": b, "role", "derived": bool}],
                 "sources": [id], "confidence", "drafted_by", "checked_by", "note"}],
  "ties": [{"id", "a": b, "b": b, "type", "colour", "confidence", "sources": [id], "note",
            "attested": [{"occasion": o, "sign": ±1}], "ends": o | null,
            "life": [o_first, o_last],
            "rungs": [{"occasion": o, "sign": ±1, "attested": bool, "expected": m}]}],
  "problems": [{"severity": "error" | "warning", "where": id, "message"}],
  "complex": {
    "vertices": [o, ...],
    "edges": [{"being": b, "from": o, "to": o}],
    "faces": [{"tie": t, "from": o, "to": o, "signs": [s, s], "attested": [bool, bool],
               "flip": bool, "boundary": [±(edge + 1), ...]}],
    "atlas": { the Complex Atlas's JSON of the complex, with positions }
  },
  "invariants": {
    "counts": [V, E, F], "euler", "components",
    "homology": [{"dimension", "rank", "torsion": [..], "text"}],
    "pi1": {"generators", "relators", "simplified_generators", "simplified_relators", "text"},
    "loops":  [{"steps": [{"edge": e, "being": b, "from": o, "to": o, "forward": bool}],
                "beings": [b, ...], "occasions": [o, ...]}],
    "shells": [{"beings": [b, b, b], "ties": [t, t, t], "faces": [f, f, f],
                "from": o, "to": o, "product_from": ±1, "product_to": ±1, "flips"}],
    "shell_rank"
  },
  "readings": {
    "persistence": [{"colour", "lenses", "flips", "p", "memory", "beta"}],
    "flips": [f, ...],
    "triangles": [{"occasion": o, "beings": [b, b, b], "ties": [t, t, t],
                   "product": ±1, "expected", "attested": bool}],
    "moments": [{"year", "jdn", "ties", "triangles", "frustrated", "expected"}],
    "colours": {"lenses": {colour: n},
                "carriers": [{"being": b, "colours": {colour: n}}],
                "meetings": [{"occasion": o, "colours": [colour, ...]}]}
  }
}
```

- `rungs` are a tie's values at every joint occasion within its life: attested there, or carried from the last attestation with `expected` = sign · m^k, where m is the memory of the tie's colour and k the lenses since.
- `persistence` counts only lenses with both ends attested: `p` = flips / lenses, `memory` = 1 − 2p = tanh β, `beta` = artanh(memory) (null when memory ≤ 0); its last entry, colour `all`, pools every colour.
- A loop's `steps` follow a closed walk on the 1-skeleton, `forward` when the walk goes with time; `beings` are those it passes, in order, and `occasions` those where it passes from one being to another.
- A shell is minimal: three beings pairwise tied, all present at two consecutive meetings of the three ties; `product_to` = `product_from` · (−1)^`flips` (the Bianchi identity).

The answer also carries the spin foam and the other complexes (the domain's theory note, §6–7):

```
  "foam": {
    "couplings": [{"of", "n", "k", "p", "t", "beta"}],          // by colour, "all", and "balance"
    "foams",                                                    // b₂ of K̂ over ℤ/2
    "cells": [{"kind": "lens" | "triangle", "index", "beta", "t", "prior", "spin", "support", "posterior"}],
    "triangles": [{"occasion": o, "beings": [b, b, b], "ties": [t, t, t]}],
    "ties": [{"tie": t, "posterior": [⟨u⟩ at each rung]}],
    "lives": [{"tie": t, "lenses", "wilson", "error", "chain"}],
    "moments": [{"year", "triangles", "mean", "lo", "hi"}],      // the reconstruction's frustration
    "sweeps",
    "su2": {"segments": [{"live", "channels"}], "carrying", "admissible"}
  },
  "models": {
    "refined": {"counts", "euler", "homology", "triangles", "foams"},
    "filled": {"counts", "euler", "homology"},
    "dowker": {"simplices", "homology", "q"},
    "time": {"intervals": [{"dimension", "birth": o, "death": o | null, "creator", "killer"}], "betti": [[b₀, b₁, b₂] by stage]}
  }
```

### 3.3 The Answer to `moment`

```
{"date", "jdn",
 "nodes": [b, ...],
 "links": [{"tie": t, "a": b, "b": b, "sign": ±1, "expected", "since": o, "lenses_since"}],
 "triangles": [{"beings": [b, b, b], "ties": [t, t, t], "product": ±1, "expected"}],
 "frustrated", "expected",
 "posterior": {"links": [⟨u⟩], "triangles": [⟨product⟩], "frustrated": {"mean", "lo", "hi"}}}
```

### 3.4 The Answer to `power`

```
{"beings": [{"being": b, "ties", "divergence", "rungs"}], "balance": {"divergence", "rungs"}, "sweeps"}
```

### 3.5 The Server

| Path | |
|---|---|
| `GET /api/engine` | whether the engine is built and current |
| `GET /api/canons` | the canons in `domains/history/canons/`, each its first line, and `integrated`, all of them joined, when there are two or more |
| `POST /api/build` `{"canon": id}` | the engine's `build` |
| `POST /api/moment` `{"canon": id, "date"}` | the engine's `moment` |
| `POST /api/power` `{"canon": id}` | the engine's `power` |
| `GET /api/check?canon=id` | the last check of the canon's beings against Wikidata, or `{"status": "none"}` |
| `POST /api/check/run` `{"canon": id}` | runs the check (`data/products/history_check.py`) in the background |
| `GET /setup`, `POST /setup/run` | the recipe: building the engine |

## 4. Files

| Path | |
|---|---|
| `index.html`, `styles.css`, `src/` | The page: ES modules, no build step |
| `serve.py`, `engine.py`, `recipe.json` | The server, the engine's builder, the recipe |
| `engine/history.cpp` | The engine |
| `src/foam3d.js`, `src/foamlayout.js` | The foam in three dimensions (three.js 0.170.0), and its layout |
| `src/models.js`, `src/foamreadings.js` | The other complexes and persistence; the spin foam's readings |
| `tests/` | The application's checks, and the engine's answers as fixtures |
