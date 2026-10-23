"""Serve the Complex Atlas on this computer.

The page is served from this directory, and the code it shares from the repository's
``domains/``, ``foundations/`` and ``visualization/`` under the same paths, with the data cache
at ``/data/``. The complexes it shows are built and analysed by the engine (``engine.py``: the
foundations' C++, morphe's topology and physis's gauge fields and state sums):

    GET  /api/engine                     whether the engine is built and current, and the limits of a request
    GET  /api/catalogue                  the families of complexes, from the engine
    POST /api/build                      {"family", "params"}: a complex of a family
    POST /api/analyse                    {"complex"}: its invariants
    POST /api/gauge                      {"complex", "group", "beta", ...}: a gauge field on it, by Monte Carlo
    POST /api/turaev-viro                {"complex", "r"}: the Turaev-Viro state sum of a 3-manifold
    GET  /api/theories                   the field theories that generate spacetimes (physis decision 0006)
    POST /api/expand                     {"theory", "params", "order"}: its Feynman graphs and series
    POST /api/graph-complex              {"coloured_graph", "view"}: the 2-complex or pseudomanifold of a graph
    POST /api/large-n                    {"theory", "params"}: its large-N phase
    POST /api/condensate                 {"spins", "G", "phi", "points"}: a GFT condensate's cosmology
    POST /api/random-geometry            {"kind", "size", "D", "seed"}: a random quadrangulation or melonic graph
    POST /api/dimensions                 {"complex", "graph", ...}: its spectral and Hausdorff dimensions
    POST /api/spin-foam                  {"complex", "r", "mode", ...}: a spin foam on a 3-manifold's dual
    POST /api/cancel                     stop the computations in progress

The engine is built by the application's recipe (``recipe.json``), which the page may run:

    GET  /setup                          its one step, and whether it is done
    POST /setup/run[?again=1]            build the engine (again)

Only 127.0.0.1 is bound unless ``--host`` says otherwise. On a public address, building needs
the server's token (``--token``), and requests are held within smaller limits.

Usage::

    python applications/complex-atlas/serve.py [--host 127.0.0.1] [--port 8850] [--open]

then open http://127.0.0.1:8850/ in a browser.
"""

from __future__ import annotations

import argparse
import http.server
import json
import sys
import threading
import urllib.parse
import webbrowser
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
MOUNTS = {
    "/domains/": ROOT / "domains",
    "/foundations/": ROOT / "foundations",
    "/visualization/": ROOT / "visualization",
    "/data/": ROOT / "data" / "cache",
}
TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json",
    ".svg": "image/svg+xml",
    ".png": "image/png",
}

sys.path.insert(0, str(APP))
import engine  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
# The limits of a request: the cells of a complex, the sweeps of a Monte Carlo run, the order
# of a state sum, the wall time of each kind of request.
LIMITS = {"cells": 250_000, "sweeps": 200_000, "r": 9, "order": 6, "size": 20_000,
          "seconds": {"build": 60.0, "analyse": 120.0, "gauge": 300.0, "turaev_viro": 300.0, "theories": 30.0, "expand": 180.0,
                      "graph_complex": 30.0, "large_n": 60.0, "condensate": 60.0, "random_geometry": 60.0, "dimensions": 120.0,
                      "spin_foam": 300.0}}
RUNNING: set = set()
SLOTS = threading.BoundedSemaphore(3)
OPS = {"/api/build": "build", "/api/analyse": "analyse", "/api/gauge": "gauge", "/api/turaev-viro": "turaev_viro",
       # The generation of spacetime by field theories (physis decision 0006): theories and their expansion,
       # the complexes of Feynman graphs, large N and condensates, random geometries, dimensions, spin foams.
       "/api/theories": "theories", "/api/expand": "expand", "/api/graph-complex": "graph_complex", "/api/large-n": "large_n",
       "/api/condensate": "condensate", "/api/random-geometry": "random_geometry", "/api/dimensions": "dimensions",
       "/api/spin-foam": "spin_foam"}


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


def write_refusal(handler) -> str | None:
    if not WRITES["allowed"]:
        return "this server does not build: start it with --token to allow building with a token"
    token = WRITES["token"]
    if token and handler.headers.get("Authorization", "") != f"Bearer {token}":
        return "building on this server needs its token"
    return None


def cells_of(complex_: dict) -> int:
    """The number of cells of a complex as the page sends it."""
    try:
        n = int(complex_.get("vertices", 0)) + len(complex_.get("edges", [])) + len(complex_.get("faces", []))
        n += sum(len(layer) for layer in complex_.get("cells", []) or [])
        return n
    except (TypeError, ValueError, AttributeError):
        return 0


def refuse_beyond_limits(op: str, body: dict) -> str | None:
    """The reason a request exceeds this server's limits, or None. The engine checks everything else."""
    c = body.get("complex")
    if c is not None:
        if not isinstance(c, dict):
            return "a complex is a JSON object"
        n = cells_of(c)
        if n > LIMITS["cells"]:
            return f"{n:,} cells: this server takes at most {LIMITS['cells']:,}"
    if op == "gauge":
        try:
            sweeps = int(body.get("sweeps", 0)) + int(body.get("thermalise", 0))
        except (TypeError, ValueError):
            return "sweeps and thermalisation are whole numbers"
        if sweeps > LIMITS["sweeps"]:
            return f"{sweeps:,} sweeps: this server runs at most {LIMITS['sweeps']:,}"
    if op in ("turaev_viro", "spin_foam"):
        try:
            r = int(body.get("r", 0))
        except (TypeError, ValueError):
            return "r is a whole number"
        if r > LIMITS["r"]:
            return f"r = {r}: this server sums at most to r = {LIMITS['r']}"
    if op == "spin_foam":
        try:
            sweeps = int(body.get("sweeps", 0))
        except (TypeError, ValueError):
            return "sweeps are a whole number"
        if sweeps > LIMITS["sweeps"]:
            return f"{sweeps:,} sweeps: this server runs at most {LIMITS['sweeps']:,}"
    if op == "expand":
        try:
            order = int(body.get("order", 0))
        except (TypeError, ValueError):
            return "the order is a whole number"
        if order > LIMITS["order"]:
            return f"order {order}: this server expands at most to order {LIMITS['order']}"
    if op == "random_geometry":
        try:
            size = int(body.get("size", 0))
        except (TypeError, ValueError):
            return "the size is a whole number"
        if size > LIMITS["size"]:
            return f"size {size:,}: this server makes geometries of at most {LIMITS['size']:,}"
    return None


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **TYPES}

    def translate_path(self, path: str) -> str:
        path = urllib.parse.unquote(path.split("?", 1)[0].split("#", 1)[0])
        prefix = next((p for p in MOUNTS if path.startswith(p)), None)
        base, rel = (MOUNTS[prefix], path[len(prefix):]) if prefix else (APP, path.lstrip("/"))
        target = (base / rel).resolve()
        if target != base and base not in target.parents:
            return str(base / "__outside__")
        return str(target)

    def do_GET(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/api/engine":
            return self.reply_json({**engine.status(), "limits": LIMITS, "running": len(RUNNING)}, 200)
        if url.path == "/api/catalogue":
            c = engine.catalogue()
            return self.reply_json(c if c is not None else {"message": "the engine is not built"}, 200 if c is not None else 503)
        if url.path == "/setup":
            return self.reply_json(*self.recipe_request(url))
        if url.path == "/api/theories":
            return self.engine_request("theories", {})
        return super().do_GET()

    def do_POST(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/setup/run":
            return self.reply_json(*self.recipe_request(url))
        if url.path == "/api/cancel":
            runs = list(RUNNING)
            for run in runs:
                run.cancel()
            return self.reply_json({"cancelled": len(runs)}, 200)
        if url.path in OPS:
            return self.engine_request(OPS[url.path])
        self.send_error(404)

    def engine_request(self, op: str, body: dict | None = None) -> None:
        if body is None:
            length = int(self.headers.get("Content-Length", "0") or 0)
            if length > 32 * 1024 * 1024:
                return self.reply_json({"message": "the request is too large"}, 413)
            try:
                body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
                if not isinstance(body, dict):
                    raise ValueError("a request is a JSON object")
            except ValueError as error:
                return self.reply_json({"message": f"not a request: {error}"}, 400)
        reason = refuse_beyond_limits(op, body)
        if reason:
            return self.reply_json({"message": reason}, 400)
        state = engine.status()
        if not state["current"]:
            return self.reply_json({"message": "the engine is not built from the present sources: build it first", "engine": state}, 503)
        # A page fires several requests at once (a theory's expansion, its large-N phase, a foam):
        # each waits its turn for one of the engine's slots rather than failing.
        if not SLOTS.acquire(timeout=120):
            return self.reply_json({"message": "the engine is busy with other requests; try again in a moment"}, 409)
        try:
            run = engine.Run(json.dumps({**body, "op": op}))
            RUNNING.add(run)
            code, out, err = run.execute(timeout=LIMITS["seconds"][op])
        finally:
            RUNNING.discard(run)
            SLOTS.release()
        try:
            answer = json.loads(out) if out.strip() else None
        except ValueError:
            answer = None
        if isinstance(answer, dict) and "error" in answer:
            return self.reply_json({"message": str(answer["error"])}, 400)
        if code != 0 or answer is None:
            message = err.strip().splitlines()[-1] if err.strip() else "the engine failed"
            return self.reply_json({"message": message}, 499 if code < 0 else 500)
        # The answer is passed through as the engine wrote it.
        self.reply_text(out, 200)

    def recipe_request(self, url) -> tuple[dict, int]:
        refusal = write_refusal(self)
        policy = {"writable": refusal is None, "refusal": refusal}
        try:
            if url.path == "/setup":
                return {**recipe_job().view(), **policy}, 200
            if refusal:
                return {"message": refusal, **policy}, 403
            query = dict(urllib.parse.parse_qsl(url.query))
            return {**recipe_job().start(None, again=query.get("again") == "1"), **policy}, 200
        except RuntimeError as error:
            return {"message": str(error)}, 409

    def reply_json(self, body: dict, code: int) -> None:
        self.reply_text(json.dumps(body, ensure_ascii=False), code)

    def reply_text(self, text: str, code: int) -> None:
        data = text.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, format: str, *args) -> None:
        if args and str(args[1]).startswith(("4", "5")):
            super().log_message(format, *args)


class Server(http.server.ThreadingHTTPServer):
    request_queue_size = 64
    daemon_threads = True


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Serve the Complex Atlas locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8850)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows building the engine")
    parser.add_argument("--max-cells", type=int, help="the most cells of a complex in a request")
    parser.add_argument("--max-sweeps", type=int, help="the most Monte Carlo sweeps in a request")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
        LIMITS.update(cells=30_000, sweeps=20_000, r=7, order=4, size=4_000,
                      seconds={"build": 20.0, "analyse": 30.0, "gauge": 60.0, "turaev_viro": 60.0, "theories": 10.0, "expand": 60.0,
                               "graph_complex": 10.0, "large_n": 20.0, "condensate": 20.0, "random_geometry": 20.0,
                               "dimensions": 40.0, "spin_foam": 60.0})
    if args.max_cells is not None:
        LIMITS["cells"] = args.max_cells
    if args.max_sweeps is not None:
        LIMITS["sweeps"] = args.max_sweeps
    if not engine.status()["current"]:
        print("The engine is not built from the present sources. Build it from the page, or run:\n"
              "  python data/acquisition/recipe.py applications/complex-atlas/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if not public else f"{url} and on {args.host}"
    print(f"Complex Atlas at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
