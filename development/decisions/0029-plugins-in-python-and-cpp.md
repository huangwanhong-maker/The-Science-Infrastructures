# 0029. Plug-ins in Python and C++: the Same Contract, Run by the Server

- **Status:** Accepted, provisionally: built at the user's request; open to revision on review. Amends [ADR 0027](0027-earth-atlas-plugins.md). Amended by [ADR 0030](0030-the-api-edge.md): the Python classes, the field products' reading and writing, and the hosting moved to `api/python/science`, the C++ classes to `api/cpp/plugin.hpp`; the Earth Atlas binds them.
- **Date:** 2026-10-22
- **Decided by:** the user asked that plug-ins be possible in Python and C++ as well as JavaScript. Claude chose the design below.

## 1. Context

ADR 0027 made a plug-in an ES module that the page loads and runs. Algorithms in Python, the language of the data library and of numpy, and in C++, the language of the computational core to come ([ADR 0003](0003-cpp-computational-core.md)), cannot run in the page. The application's server already runs Python for the page: requests on demand ([ADR 0023](0023-ocean-column-and-data-on-demand.md)) and the data recipes. ADR 0027 §3 foresaw the classes mirrored in Python, with the server as the context.

## 2. Decision

1. **One contract, three languages.** The classes, their hooks, the needs, the parameters and the results have the same names and meanings in JavaScript (`src/plugins/plugin.js`), in Python (`py/plugins.py`) and in C++ (`plugins/cpp/plugin.hpp`). The page loads, runs and shows a plug-in the same way whatever its language, and names the language beside its title.
2. **Python plug-ins run in the server's process.** `py/host.py` imports the script and asks it for its plug-in, resolves its needs from the products on disk, which `py/fields.py` reads in the page's own format, runs it on numpy arrays, and writes each field result as a product under `data/cache/plugins/<id>/`, which the page reads as it reads any product; sketches, tables and notes pass as JSON. The server exposes three requests: describe, run, and a job's state.
3. **C++ plug-ins are programs.** The host compiles one with the compiler at hand (`c++ -std=c++20 -O2`) when its binary is older than its source or the header, asks it to describe itself (`--describe` prints JSON), and runs it on files (`--run <dir>`: `input.txt` and one float32 block per input variable in; `output.txt` and one float32 block per output out). The header gives `Diagnostic`, and only the standard library is used, so any C++20 compiler builds a plug-in. The same file contract serves any language that reads and writes files.
4. **Only files under the application's `plugins/` are run by the server.** A script elsewhere is not executed on the serving computer; the page's own plug-ins keep loading from any URL, as ADR 0027 allows.
5. **Checked without a page or a server.** `py/check.py` runs the Python plug-in, writes a product and reads it back within half a quantum, and compiles and runs the C++ plug-in when a compiler exists; the Earth Atlas's suite runs it, and skips what the computer cannot do.

## 3. Alternatives Considered

- **Python and C++ in the page,** through Pyodide and WebAssembly. Possible, and it would keep every plug-in in one place; but it brings a toolchain and a 10 MB runtime, loses numpy's speed and the server's direct reach to the data, and the server already exists.
- **One process contract for all three languages.** Python in the process keeps numpy arrays without copies and tracebacks in hand; the JavaScript runtime keeps its directness. The process contract is for languages that need it.
- **A shared library for C++,** loaded into the server. Faster per call, but platform-specific loading, and a fault in a plug-in would bring the server down; a process is isolated and simple.

## 4. Consequences

- **New files:** `applications/earth-atlas/py/` (`fields.py`, `plugins.py`, `host.py`, `check.py`), `plugins/cpp/plugin.hpp`, and the same diagnostic in three languages, `vorticity.js`, `vorticity.py` and `vorticity.cpp`, checked against the same closed form.
- **The C++ classes stop at `Diagnostic`.** `ForecastModel`, `DataAssimilation` and `Trajectories` in C++ are the next step; their generic parts (verification, increments, the parcel step) can stay in the host, as the Python base classes show.
- **What a computer cannot do, the panel says:** no compiler, no C++ plug-in; no numpy, no Python plug-in.
- **Security, as ADR 0027:** a plug-in the server runs has the server's authority. Only files in the repository's `plugins/` run.
