"""Serve the Economy Atlas on this computer.

The page is served from this directory, and the code it shares from the repository's
``domains/``, ``foundations/`` and ``visualization/`` under the same paths, with the data cache
at ``/data/``. What it shows is computed by the engine (``engine.py``: the economics domain's
C++, on physis's discrete exterior calculus and morphe's complexes) from the data library's
products. The protocol is in the README, §3:

    GET  /api/engine                     whether the engine is built and current, and the limits of a request
    GET  /api/data                       which products are built: the trade flows' years, the currencies' months
    POST /api/map                        {"year"}: the world's trade in a year, its exchange, drive and circulation
    POST /api/spacetime                  {"first", "last", "currencies"?}: currencies held and exchanged
    POST /api/production                 {"table", "value_added_coefficients"?, "final_demand"?}: a Leontief model
    POST /api/circuits                   {"network"}: a transfer network's deficits, prices and circulation
    POST /api/cancel                     stop the computations in progress

The engine and the data are made by the application's recipe (``recipe.json``), which the page
may run:

    GET  /setup                          its steps, and whether each is done
    POST /setup/run[?again=1]            run them (again)

Answers for a year or a span are kept, keyed by the engine's build and the product, so that each
is computed once. Only 127.0.0.1 is bound unless ``--host`` says otherwise; on a public address
the recipe needs the server's token (``--token``), and requests are held within smaller limits.

Usage::

    python applications/economy-atlas/serve.py [--host 127.0.0.1] [--port 8910] [--open]

then open http://127.0.0.1:8910/ in a browser.
"""

from __future__ import annotations

import argparse
import hashlib
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
    ".geojson": "application/geo+json",
    ".md": "text/markdown; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
}
PRODUCTS = ROOT / "data" / "cache" / "products"
TRADE = PRODUCTS / "trade-flows.json"
CURRENCIES = PRODUCTS / "currency-spacetime.json"
KEPT = ROOT / "data" / "cache" / "economy-atlas" / "answers"

sys.path.insert(0, str(APP))
import engine  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
LIMITS = {"bytes": 1024 * 1024, "seconds": {"map": 120.0, "spacetime": 120.0, "production": 30.0, "circuits": 30.0}}
RUNNING: set = set()
SLOTS = threading.BoundedSemaphore(3)
OPS = {"/api/map": "map", "/api/spacetime": "spacetime", "/api/production": "production", "/api/circuits": "circuits"}


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


def write_refusal(handler) -> str | None:
    if not WRITES["allowed"]:
        return "this server does not run its recipe: start it with --token to allow it with a token"
    token = WRITES["token"]
    if token and handler.headers.get("Authorization", "") != f"Bearer {token}":
        return "running the recipe on this server needs its token"
    return None


def data_state() -> dict:
    """The products built, and what the page may ask of them."""
    out: dict = {"trade": None, "currencies": None}
    try:
        index = json.loads((PRODUCTS / "index.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        index = {}
    if TRADE.exists():
        entry = next((e for e in index.get("trade", []) if e["name"] == "trade-flows"), None)
        out["trade"] = entry or {"name": "trade-flows"}
    if CURRENCIES.exists():
        entry = next((e for e in index.get("currencies", []) if e["name"] == "currency-spacetime"), None)
        out["currencies"] = entry or {"name": "currency-spacetime"}
    return out


def answer_key(op: str, request: dict, product: Path | None) -> str:
    h = hashlib.sha256(json.dumps({"engine": engine.read_build() or {}, "op": op, "request": request}, sort_keys=True).encode("utf-8"))
    if product is not None:
        stat = product.stat()
        h.update(f"{product.name}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    return h.hexdigest()[:32]


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
        if url.path == "/api/data":
            return self.reply_json(data_state(), 200)
        if url.path == "/setup":
            return self.reply_json(*self.recipe_request(url))
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

    def read_body(self) -> dict | None:
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > LIMITS["bytes"]:
            self.reply_json({"message": "the request is too large"}, 413)
            return None
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            if not isinstance(body, dict):
                raise ValueError("a request is a JSON object")
            return body
        except ValueError as error:
            self.reply_json({"message": f"not a request: {error}"}, 400)
            return None

    def engine_request(self, op: str) -> None:
        body = self.read_body()
        if body is None:
            return
        toolchain = (engine.read_build() or {}).get("toolchain") or "native"
        request: dict = {"op": op}
        product = None
        if op == "map":
            if not TRADE.exists():
                return self.reply_json({"message": "the trade flows are not built: run the recipe"}, 404)
            product = TRADE
            request["flows"] = engine.path_for(toolchain, TRADE)
            if isinstance(body.get("year"), int):
                request["year"] = body["year"]
        elif op == "spacetime":
            if not CURRENCIES.exists():
                return self.reply_json({"message": "the currencies' record is not built: run the recipe"}, 404)
            product = CURRENCIES
            request["record"] = engine.path_for(toolchain, CURRENCIES)
            for key in ("first", "last"):
                if isinstance(body.get(key), str):
                    request[key] = body[key]
            if isinstance(body.get("currencies"), list):
                request["currencies"] = [c for c in body["currencies"] if isinstance(c, str)]
        elif op == "production":
            if not isinstance(body.get("table"), dict):
                return self.reply_json({"message": "a production request needs a table"}, 400)
            request["table"] = body["table"]
            for key in ("value_added_coefficients", "final_demand"):
                if isinstance(body.get(key), list):
                    request[key] = body[key]
        elif op == "circuits":
            if not isinstance(body.get("network"), dict):
                return self.reply_json({"message": "a circuits request needs a network"}, 400)
            request["network"] = body["network"]
        state = engine.status()
        if not state["current"]:
            return self.reply_json({"message": "the engine is not built from the present sources: build it first", "engine": state}, 503)
        keyed = {k: v for k, v in request.items() if k not in ("flows", "record")}
        kept = KEPT / f"{answer_key(op, keyed, product)}.json" if product is not None else None
        if kept is not None and kept.exists():
            return self.reply_text(kept.read_text(encoding="utf-8"), 200)
        if not SLOTS.acquire(timeout=120):
            return self.reply_json({"message": "the engine is busy with other requests; try again in a moment"}, 409)
        run = None
        try:
            run = engine.Run(json.dumps(request))
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
            message = err.replace("\0", "").strip().splitlines()[-1] if err.strip() else "the engine failed"
            return self.reply_json({"message": message}, 499 if code < 0 else 500)
        if kept is not None:
            try:
                KEPT.mkdir(parents=True, exist_ok=True)
                kept.write_text(out, encoding="utf-8")
            except OSError:
                pass
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
    parser = argparse.ArgumentParser(description="Serve the Economy Atlas locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8910)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows running the recipe")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
        LIMITS.update(bytes=256 * 1024, seconds={"map": 60.0, "spacetime": 60.0, "production": 10.0, "circuits": 10.0})
    if not engine.status()["current"]:
        print("The engine is not built from the present sources. Build it from the page, or run:\n"
              "  python data/acquisition/recipe.py applications/economy-atlas/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if not public else f"{url} and on {args.host}"
    print(f"Economy Atlas at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
