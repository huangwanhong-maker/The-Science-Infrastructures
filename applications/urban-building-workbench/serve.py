"""Serve the Urban Building Workbench on this computer.

The page is served from this directory, with the repository's ``domains/``, ``foundations/`` and
``visualization/`` under the same paths and the data cache at ``/data/``. Its plug-ins live in
``plugins/<name>/``: a manifest (``plugin.json``) naming its views, the views' modules, and, for
those with C++, an engine built by ``engine.py`` and kept running here, one process each.

    GET  /api/plugins                    the plug-ins: their manifests and whether their engines are current
    POST /api/<plugin>/<op>              {…}: one request to a plug-in's engine
    GET  /setup                          the recipe (recipe.json): building the engines, and whether it is done
    POST /setup/run[?again=1]            build them (again)

Answers to requests without parameters are kept until the engine is rebuilt. Only 127.0.0.1 is
bound unless ``--host`` says otherwise; on a public address, building needs ``--token``.

Usage::

    python applications/urban-building-workbench/serve.py [--host 127.0.0.1] [--port 8970] [--open]

then open http://127.0.0.1:8970/ in a browser.
"""

from __future__ import annotations

import argparse
import http.server
import json
import sys
import threading
import time
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
}

sys.path.insert(0, str(APP))
import engine  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
WORKERS: dict[str, engine.Worker] = {}
CACHE: dict[tuple[str, str], str] = {}
LOCK = threading.Lock()


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


def worker(name: str) -> engine.Worker:
    with LOCK:
        if name not in WORKERS:
            WORKERS[name] = engine.Worker(name)
        return WORKERS[name]


STATUS: dict[str, tuple[float, dict]] = {}


def engine_status(name: str, fresh: bool = False) -> dict:
    """An engine's status, kept for 30 s: judging it hashes every source and header it includes, which takes
    seconds, and every request to the engine asks it."""
    now = time.monotonic()
    kept = STATUS.get(name)
    if fresh or not kept or now - kept[0] > 30:
        kept = (now, engine.status(name))
        STATUS[name] = kept
    return kept[1]


def warm() -> None:
    """Ask each current engine the requests its manifest names under "warm" (without parameters), once, so that
    the first page to ask finds their answers kept (the tides of every street take seconds to fit)."""
    for p in engine.plugins():
        if not p["engine"] or not engine_status(p["name"])["current"]:
            continue
        for op in p.get("warm", []):
            try:
                answer = worker(p["name"]).ask({"op": op})
            except RuntimeError:
                continue
            if answer.get("ok"):
                CACHE[(p["name"], op)] = json.dumps(answer["result"], ensure_ascii=False)


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
        if url.path == "/api/plugins":
            out = []
            for p in engine.plugins():
                p["engineStatus"] = engine_status(p["name"], fresh=True) if p["engine"] else None
                out.append(p)
            return self.reply_json({"plugins": out}, 200)
        if url.path == "/setup":
            return self.reply_json(*self.recipe_request(url))
        return super().do_GET()

    def do_POST(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/setup/run":
            with LOCK:
                CACHE.clear()
                STATUS.clear()
                for w in WORKERS.values():
                    w.stop()
            return self.reply_json(*self.recipe_request(url))
        parts = url.path.strip("/").split("/")
        if len(parts) == 3 and parts[0] == "api":
            return self.engine_request(parts[1], parts[2])
        self.send_error(404)

    def engine_request(self, name: str, op: str) -> None:
        if name not in {p["name"] for p in engine.plugins() if p["engine"]}:
            return self.reply_json({"message": f"no plug-in {name} with an engine"}, 404)
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > 4 * 1024 * 1024:
            return self.reply_json({"message": "the request is too large"}, 413)
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            if not isinstance(body, dict):
                raise ValueError("a request is a JSON object")
        except ValueError as error:
            return self.reply_json({"message": f"not a request: {error}"}, 400)
        state = engine_status(name)
        if not state["current"]:
            return self.reply_json({"message": "the engine is not built from the present sources: build it first", "engine": state}, 503)
        key = (name, op) if not body else None
        if key and key in CACHE:
            return self.reply_text(CACHE[key], 200)
        try:
            answer = worker(name).ask({**body, "op": op})
        except RuntimeError as error:
            return self.reply_json({"message": str(error)}, 500)
        if not answer.get("ok"):
            return self.reply_json({"message": str(answer.get("error", "the engine failed"))}, 400)
        text = json.dumps(answer["result"], ensure_ascii=False)
        if key:
            CACHE[key] = text
        self.reply_text(text, 200)

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
    parser = argparse.ArgumentParser(description="Serve the Urban Building Workbench locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8970)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows building the engines")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
    stale = [p["name"] for p in engine.plugins() if p["engine"] and not engine.status(p["name"])["current"]]
    if stale:
        print(f"Engines not built from the present sources ({', '.join(stale)}). Build them from the page, or run:\n"
              "  python data/acquisition/recipe.py applications/urban-building-workbench/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    threading.Thread(target=warm, daemon=True).start()
    url = f"http://127.0.0.1:{args.port}/"
    print(f"Urban Building Workbench at {url if not public else url + ' and on ' + args.host}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        for w in WORKERS.values():
            w.stop()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
