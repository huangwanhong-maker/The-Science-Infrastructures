"""Serve Fabrica on this computer.

The page is served from this directory, with the repository's ``domains/``, ``foundations/`` and
``visualization/`` under the same paths and the data cache at ``/data/``. The engine (``engine.py``) is kept
running here and answers the page's questions of the anatomy domain:

    GET  /api/status                     whether the product and the engine are ready
    GET  /api/catalogue[?body=…]         every part of a body measured, the trees, what touches what, the levels;
                                         computed once (about 20 s) and kept in data/cache/fabrica/
    POST /api/section | part | profile   {…}: one request to the engine
    POST /api/cases | reconstruct | slice  {…}: vessels reconstructed from images (ADR 0070)
    GET  /setup                          the recipe (recipe.json): the data, the product, the engine
    POST /setup/run[?again=1]            run it (again)

Only 127.0.0.1 is bound unless ``--host`` says otherwise; on a public address, running the recipe needs
``--token``.

Usage::

    python applications/fabrica/serve.py [--host 127.0.0.1] [--port 8990] [--open]

then open http://127.0.0.1:8990/ in a browser.
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
    ".md": "text/markdown; charset=utf-8",
    ".svg": "image/svg+xml",
    ".gz": "application/gzip",
}
PRODUCTS = ROOT / "data" / "cache" / "products"
PRODUCT = PRODUCTS / "body-bp3d.json"
# The bodies: BodyParts3D's man (the default) and the Human Reference Atlas's woman, when her product is built.
BODIES = ("body-bp3d", "body-hra-f")
KEPT = ROOT / "data" / "cache" / "fabrica"

sys.path.insert(0, str(APP))
import engine  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
WORKER = engine.Worker()
LOCK = threading.Lock()
CATALOGUE: dict[str, str] = {}


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


def catalogue_key(body: str = "body-bp3d") -> str | None:
    """A body's catalogue depends on the engine's sources and its product's build."""
    state = engine.status()
    product = PRODUCTS / f"{body}.json"
    if body not in BODIES or not state["current"] or not product.exists():
        return None
    built = json.loads(product.read_text(encoding="utf-8")).get("built_at", "")
    return body + "-" + hashlib.sha256(f"{state['sources']}|{built}".encode()).hexdigest()[:16]


def catalogue_text(body: str = "body-bp3d") -> str:
    key = catalogue_key(body)
    if key is None:
        raise RuntimeError("the product or the engine is not ready: run the recipe")
    with LOCK:
        if key in CATALOGUE:
            return CATALOGUE[key]
        path = KEPT / f"catalogue-{key}.json"
        if path.exists():
            CATALOGUE[key] = path.read_text(encoding="utf-8")
            return CATALOGUE[key]
    answer = WORKER.ask({"op": "catalogue", "body": body})
    if not answer.get("ok"):
        raise RuntimeError(str(answer.get("error", "the engine failed")))
    text = json.dumps(answer["result"], ensure_ascii=False, separators=(",", ":"))
    with LOCK:
        KEPT.mkdir(parents=True, exist_ok=True)
        for old in KEPT.glob(f"catalogue-{body}-*.json"):
            old.unlink(missing_ok=True)
        (KEPT / f"catalogue-{key}.json").write_text(text, encoding="utf-8")
        CATALOGUE[key] = text
    return text


def warm() -> None:
    try:
        catalogue_text()
    except RuntimeError:
        pass


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
        if url.path == "/api/status":
            bodies = [b for b in BODIES if (PRODUCTS / f"{b}.json").exists()]
            return self.reply_json({"product": PRODUCT.exists(), "engine": engine.status(), "catalogue": catalogue_key() is not None,
                                    "bodies": bodies}, 200)
        if url.path == "/api/catalogue":
            body = urllib.parse.parse_qs(url.query).get("body", ["body-bp3d"])[0]
            try:
                return self.reply_text(catalogue_text(body), 200)
            except RuntimeError as error:
                return self.reply_json({"message": str(error)}, 503)
        if url.path == "/setup":
            return self.reply_json(*self.recipe_request(url))
        return super().do_GET()

    def do_POST(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/setup/run":
            with LOCK:
                CATALOGUE.clear()
            WORKER.stop()
            return self.reply_json(*self.recipe_request(url))
        parts = url.path.strip("/").split("/")
        if len(parts) == 2 and parts[0] == "api" and parts[1] in ("section", "part", "profile", "cases", "reconstruct", "slice"):
            return self.engine_request(parts[1])
        self.send_error(404)

    def engine_request(self, op: str) -> None:
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > 1024 * 1024:
            return self.reply_json({"message": "the request is too large"}, 413)
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            if not isinstance(body, dict):
                raise ValueError("a request is a JSON object")
        except ValueError as error:
            return self.reply_json({"message": f"not a request: {error}"}, 400)
        if not engine.status()["current"]:
            return self.reply_json({"message": "the engine is not built from the present sources: run the recipe"}, 503)
        try:
            answer = WORKER.ask({**body, "op": op})
        except RuntimeError as error:
            return self.reply_json({"message": str(error)}, 500)
        if not answer.get("ok"):
            return self.reply_json({"message": str(answer.get("error", "the engine failed"))}, 400)
        self.reply_text(json.dumps(answer["result"], ensure_ascii=False, separators=(",", ":")), 200)

    def recipe_request(self, url) -> tuple[dict, int]:
        refusal = None
        if not WRITES["allowed"]:
            refusal = "this server does not build: start it with --token to allow building with a token"
        elif WRITES["token"] and self.headers.get("Authorization", "") != f"Bearer {WRITES['token']}":
            refusal = "building on this server needs its token"
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
    parser = argparse.ArgumentParser(description="Serve Fabrica locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8990)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows running the recipe")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
    if catalogue_key() is None:
        print("The product or the engine is not ready. Run the recipe from the page, or:\n"
              "  python data/acquisition/recipe.py applications/fabrica/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    threading.Thread(target=warm, daemon=True).start()
    url = f"http://127.0.0.1:{args.port}/"
    print(f"Fabrica at {url if not public else url + ' and on ' + args.host}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        WORKER.stop()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
