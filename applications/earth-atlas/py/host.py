"""The host of the Earth Atlas's plug-ins that run in the server (ADR 0029):
Python scripts, imported into this process and run on numpy, and C++ programs,
compiled and run as processes behind one file contract. The page cannot tell
them from its own JavaScript plug-ins: the application's server exposes

    GET  /plugin-host/describe?file=plugins/x.py   the plug-in's identity, parameters and needs, with the sources that serve each
    POST /plugin-host/run  {file, chosen, params, context}   starts a job and returns {job}
    GET  /plugin-host/job/<id>                      the job's state, progress and results

and the host.js of the page loads, runs and shows the results as for any plug-in.

The hosting itself belongs to no application: it is the package of the API,
api/python/science (hosting.py, with the classes in plugins.py and the products
in fields.py). This module binds it to the atlas's plugins/ directory, keeping
the names serve.py and check.py use. The API runs the same plug-ins at
/api/v1/plugins (api/README.md).

    python applications/earth-atlas/py/host.py describe plugins/vorticity.py       from the command line
    python applications/earth-atlas/py/host.py run plugins/vorticity.py [need=product|level ...] [name=value ...]
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE.parent
ROOT = APP.parents[1]
sys.path.insert(0, str(ROOT / "api" / "python"))
from science import fields, hosting  # noqa: E402,F401
from science import plugins as api  # noqa: E402

PLUGINS_DIR = APP / "plugins"
PRODUCTS = ROOT / "data" / "cache" / "products"     # what the server serves at /data/
OUT = ROOT / "data" / "cache" / "plugins"            # ... at /data/plugins/
UTC = dt.timezone.utc
MODEL_ORDER = hosting.MODEL_ORDER
HOST = hosting.Host(APP, PRODUCTS, OUT)
model_rank = hosting.model_rank
resolve = hosting.resolve
plugin_path = HOST.plugin_path
load_plugin = HOST.load_plugin
describe = HOST.describe
run = HOST.run


def Library() -> hosting.Library:
    """The atlas's products and those its plug-ins made before."""
    return hosting.Library(PRODUCTS, OUT)


def ProcessPlugin(source: Path) -> hosting.ProcessPlugin:
    """A C++ plug-in, its binary under the atlas's data/cache/plugins/bin/."""
    return hosting.ProcessPlugin(source, OUT)


class Jobs(hosting.Jobs):
    """Runs of the atlas's plug-ins for its server, one thread each; the page polls their state."""

    def __init__(self):
        super().__init__(HOST)


def main(argv: list[str]) -> int:
    if len(argv) < 2 or argv[0] not in ("describe", "run"):
        print(__doc__)
        return 2
    if argv[0] == "describe":
        print(json.dumps(describe(argv[1]), indent=1, ensure_ascii=False))
        return 0
    chosen, params = hosting.defaults(describe(argv[1]))
    for item in argv[2:]:
        key, _, value = item.partition("=")
        if key in chosen:
            chosen[key] = value
        else:
            params[key] = float(value) if value.replace(".", "", 1).lstrip("-").isdigit() else value
    _, results = run(argv[1], chosen, params, api.Context(time=dt.datetime.now(UTC), progress=lambda f, t="": print(f"  {f:4.0%} {t}")))
    for r in results:
        print(r.get("url") or r.get("title") or r.get("text") or r["kind"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
