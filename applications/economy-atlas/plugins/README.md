# plugins: Analyses of the Economy Atlas's Networks

*Status: draft (2026-10-24).* Scripts in Python and programs in C++ that the atlas's server runs on the network shown, and whose results it draws on the map ([ADR 0046](../../../development/decisions/0046-the-economy-atlas-data-and-plugins.md)). The classes belong to no application: [api/python/science/networks.py](../../../api/python/science/networks.py) and [api/cpp/network_plugin.hpp](../../../api/cpp/network_plugin.hpp). Only files directly in this directory are run.

## 1. The Contract

A plug-in is given a **network**: places, each with an id, a name, a latitude and longitude (degrees, WGS 84) and a group; flows between them, each with a value in the network's unit and, optionally, a kind; and, where the atlas computed it, its own decomposition of the net flows (the attraction of each place, the driven, local and global parts of each net leg). It answers **results**, which the atlas draws:

| Result | Drawn as |
|---|---|
| `values` | a number for each place, by id: the colour of its country, or of its circle where no country is its place; `diverging` about 0 or `sequential` |
| `arcs` | flows between places, as arcs: `[{from, to, value}]` |
| `sketch` | polylines, markers and labels in longitude and latitude |
| `table` | rows in the panel |
| `note` | text in the panel |

Its **parameters** (`number`, `integer`, `choice`, `boolean`, with a default and optional bounds) become a form; a run receives them defaulted and typed.

Both languages run as processes behind one file contract, so that a plug-in can neither stop the server nor outlive its time limit: the server writes `input.json` (`{"network", "params"}`) into a directory of the run's own, runs the plug-in with `--run <dir>`, follows the progress lines it prints (`progress\t<fraction>\t<text>`) and reads `output.json` (`{"results": [...]}`). `--describe` prints the plug-in's description. The full contract is at the top of [network_plugin.hpp](../../../api/cpp/network_plugin.hpp).

## 2. A Plug-in in Python

```python
import numpy as np

def create(api):
    class Balance(api.NetworkAnalysis):
        id = "balance"
        title = "Net inflow"
        description = "Each place's imports less its exports."
        params = [api.param("scale", "Scale", 1.0, minimum=0)]

        def run(self, network, params, context):
            context.progress(0.5, "adding up")
            b = network.balance() * params["scale"]          # aligned with network.places
            return [api.values("Net inflow", b, unit=network.unit, scale="diverging"),
                    api.note(f"{network.n} places, {len(network.value)} flows")]
    return Balance()
```

`network` has `places`, `ids`, `index`, `n`, the flows as numpy arrays `source`, `target`, `value`, and `matrix()`, `exports()`, `imports()`, `balance()`, `distances_km()` and `place(id)`. numpy is at hand; so is anything installed in the Python that runs the server.

## 3. A Plug-in in C++

```cpp
#include "network_plugin.hpp"
namespace net = science::network;
using science::Json;

struct Balance : net::Analysis {
  Json describe() const override {
    return Json::object({{"id", "balance"}, {"title", "Net inflow"}, {"description", "Imports less exports."},
                         {"params", Json::array({net::param("scale", "Scale", 1.0)})}});
  }
  std::vector<Json> run(const net::Network& n, const net::Params& p, net::Progress& progress) override {
    std::map<std::string, double> b;
    const auto balance = n.balance();
    for (std::size_t i = 0; i < n.size(); i++) b[n.places[i].id] = balance[i] * p.number("scale");
    return {net::values("Net inflow", b, n.unit, "diverging")};
  }
};
SCIENCE_NETWORK_PLUGIN(Balance)
```

The server compiles it with g++ (natively, or in WSL), C++20, against `api/cpp` and the economics domain and the foundations, so that a plug-in may use the domain itself, as [stationary.cpp](stationary.cpp) uses its stationary stocks. It is compiled again when it or a header it includes changes; the first compile takes some tens of seconds.

## 4. The Examples

| File | Language | What it does |
|---|---|---|
| [gravity.py](gravity.py) | Python | Tinbergen's gravity equation, ln F = c + a ln X + b ln M − d ln D, by least squares: each economy coloured by how far its exports lie above the line, the flows furthest above it as arcs, the coefficients |
| [concentration.py](concentration.py) | Python | The Herfindahl–Hirschman index of each economy's exports over its partners, and the most concentrated with their main partner |
| [stationary.cpp](stationary.cpp) | C++ | Where value settles if each economy passes its stock on in its export shares: the economics domain's stationary stocks of the damped rule |

## 5. Running One from the Command Line

From the repository root, with a network written as JSON (as the atlas's server writes it, or by hand):

```
python -m science.networks try applications/economy-atlas/plugins gravity.py network.json arcs=20
```

(with `api/python` on `PYTHONPATH`). A C++ plug-in compiles against `api/cpp`, and against the directories given with `-I<dir>`: the atlas's server passes the economics domain's and the foundations' (`applications/economy-atlas/engine.py`, `include_dirs()`).
