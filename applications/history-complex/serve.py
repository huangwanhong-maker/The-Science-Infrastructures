"""Serve the History Complex on this computer.

The page is served from this directory, and the code it shares from the repository's
``domains/``, ``foundations/`` and ``visualization/`` under the same paths, with the data cache
at ``/data/``. The histories it shows are the canons of the history domain
(``domains/history/canons/*.jsonl``), built into their complexes and read by the engine
(``engine.py``: the domain's C++, with morphe's topology). The protocol is in the README, §3:

    GET  /api/engine                     whether the engine is built and current, and the limits of a request
    GET  /api/canons                     the canons, each its first line
    POST /api/build                      {"canon"}: the record, the complex, its invariants and readings
    POST /api/moment                     {"canon", "date"}: the network of the moment at a date
    POST /api/power                      {"canon"}: the power of each being over the reconstruction
    POST /api/field                      {"canon"}: the history as a Feynman graph and as a sample path
    GET  /api/check?canon=<id>           the last check of the canon's beings against Wikidata
    POST /api/check/run                  {"canon"}: check it again (data/products/history_check.py)
    POST /api/cancel                     stop the computations in progress

The engine is built by the application's recipe (``recipe.json``), which the page may run:

    GET  /setup                          its one step, and whether it is done
    POST /setup/run[?again=1]            build the engine (again)

Only 127.0.0.1 is bound unless ``--host`` says otherwise. On a public address, building and
checking need the server's token (``--token``), and requests are held within smaller limits.

Usage::

    python applications/history-complex/serve.py [--host 127.0.0.1] [--port 8870] [--open]

then open http://127.0.0.1:8870/ in a browser.
"""

from __future__ import annotations

import argparse
import hashlib
import http.server
import json
import subprocess
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
    ".jsonl": "text/plain; charset=utf-8",
    ".md": "text/markdown; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
}
CANONS = ROOT / "domains" / "history" / "canons"
CHECKS = ROOT / "data" / "cache" / "products" / "history"
CHECKER = ROOT / "data" / "products" / "history_check.py"

sys.path.insert(0, str(APP))
import engine  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
# The limits of a request: the lines of a canon and the wall time of each kind of request.
LIMITS = {"lines": 200_000, "seconds": {"build": 120.0, "moment": 30.0, "power": 180.0, "field": 120.0}}
RUNNING: set = set()
SLOTS = threading.BoundedSemaphore(3)
OPS = {"/api/build": "build", "/api/moment": "moment", "/api/power": "power", "/api/field": "field"}
# Builds are kept, keyed by the engine's build and the canons' text: fitting the couplings exactly
# takes up to half a minute on the largest canons, and a canon's build never changes otherwise.
BUILDS = ROOT / "data" / "cache" / "history-complex" / "builds"


def build_key(paths: list) -> str:
    h = hashlib.sha256(json.dumps(engine.read_build() or {}, sort_keys=True).encode("utf-8"))
    for path in paths:
        h.update(path.name.encode("utf-8"))
        h.update(path.read_bytes())
    return h.hexdigest()[:32]
# The canons joined into one record (the engine's join_records): a canon of its own in the list.
INTEGRATED = "integrated"
# The checks against Wikidata in progress or ended, by canon: {"status", "started", "message"}.
CHECKING: dict[str, dict] = {}
CHECK_LOCK = threading.Lock()


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


def write_refusal(handler) -> str | None:
    if not WRITES["allowed"]:
        return "this server does not build or check: start it with --token to allow it with a token"
    token = WRITES["token"]
    if token and handler.headers.get("Authorization", "") != f"Bearer {token}":
        return "building and checking on this server need its token"
    return None


def canon_header(path: Path) -> dict | None:
    """A canon's first line, which names it, or None if it is not one."""
    try:
        with path.open(encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    head = json.loads(line)
                    break
            else:
                return None
    except (OSError, ValueError):
        return None
    if not isinstance(head, dict) or head.get("kind") != "canon":
        return None
    return {**head, "id": head.get("id") or path.stem, "file": path.name}


def canons() -> list[dict]:
    if not CANONS.is_dir():
        return []
    found = [h for h in (canon_header(p) for p in sorted(CANONS.glob("*.jsonl"))) if h]
    if len(found) > 1:
        found.append({"kind": "canon", "id": INTEGRATED, "title": "All canons, joined", "integrated": True,
                      "about": "Every canon of the domain in one record, glued along the beings, occasions and ties they share.",
                      "period": "/".join([min((h.get("period") or "/").split("/")[0] for h in found), max((h.get("period") or "/").split("/")[-1] for h in found)])})
    return found


def canon_paths(canon_id: str) -> list[Path] | None:
    """The files of a canon, or of every canon when it is the integrated one."""
    if canon_id == INTEGRATED:
        paths = sorted(p.resolve() for p in CANONS.glob("*.jsonl") if canon_header(p)) if CANONS.is_dir() else []
        return paths if len(paths) > 1 else None
    path = canon_path(canon_id)
    return [path] if path else None


def canon_path(canon_id: str) -> Path | None:
    """The file of a canon, by its id (its first line's, or its file's name)."""
    if not isinstance(canon_id, str) or not canon_id:
        return None
    for p in sorted(CANONS.glob("*.jsonl")) if CANONS.is_dir() else []:
        h = canon_header(p)
        if h and (h["id"] == canon_id or p.stem == canon_id):
            return p.resolve()
    return None


def check_file(canon_id: str) -> Path:
    return CHECKS / f"{canon_id}.check.json"


def check_state(canon_id: str) -> dict:
    with CHECK_LOCK:
        job = dict(CHECKING.get(canon_id, {}))
    if job.get("status") == "running":
        return {"status": "running", "started": job.get("started")}
    try:
        result = json.loads(check_file(canon_id).read_text(encoding="utf-8"))
        out = {"status": "done", **result}
    except (OSError, ValueError):
        out = {"status": "none"}
    if job.get("status") == "failed":
        out = {**out, "status": "failed" if out["status"] == "none" else out["status"], "message": job.get("message")}
    return out


def start_check(canon_id: str, path: Path) -> dict:
    """Checks a canon's beings against Wikidata in the background (it fetches, politely, for minutes)."""
    with CHECK_LOCK:
        if CHECKING.get(canon_id, {}).get("status") == "running":
            return {"status": "running", "started": CHECKING[canon_id]["started"]}
        CHECKING[canon_id] = {"status": "running", "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    def work() -> None:
        try:
            r = subprocess.run([sys.executable, str(CHECKER), str(path)], cwd=ROOT, capture_output=True, timeout=3600,
                               text=True, encoding="utf-8", errors="replace")
            ok = r.returncode == 0
            message = (r.stderr or r.stdout).strip().splitlines()[-1] if not ok and (r.stderr or r.stdout).strip() else None
        except (OSError, subprocess.TimeoutExpired) as error:
            ok, message = False, str(error)
        with CHECK_LOCK:
            CHECKING[canon_id] = {"status": "done" if ok else "failed", "message": message}

    threading.Thread(target=work, daemon=True).start()
    return {"status": "running", "started": CHECKING[canon_id]["started"]}


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
        if url.path == "/api/canons":
            return self.reply_json({"canons": canons()}, 200)
        if url.path == "/api/check":
            canon_id = dict(urllib.parse.parse_qsl(url.query)).get("canon", "")
            if canon_path(canon_id) is None:
                return self.reply_json({"message": f"no canon {canon_id!r}"}, 404)
            return self.reply_json(check_state(canon_id), 200)
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
        if url.path == "/api/check/run":
            body = self.read_body()
            if body is None:
                return
            refusal = write_refusal(self)
            if refusal:
                return self.reply_json({"message": refusal}, 403)
            path = canon_path(body.get("canon", ""))
            if path is None:
                return self.reply_json({"message": f"no canon {body.get('canon')!r}"}, 404)
            return self.reply_json(start_check(body["canon"], path), 200)
        if url.path in OPS:
            return self.engine_request(OPS[url.path])
        self.send_error(404)

    def read_body(self) -> dict | None:
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > 1024 * 1024:
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
        paths = canon_paths(body.get("canon", ""))
        if paths is None:
            return self.reply_json({"message": f"no canon {body.get('canon')!r} in domains/history/canons/"}, 404)
        try:
            lines = 0
            for path in paths:
                with path.open(encoding="utf-8") as f:
                    lines += sum(1 for _ in f)
        except OSError as error:
            return self.reply_json({"message": f"the canon cannot be read: {error}"}, 500)
        if lines > LIMITS["lines"]:
            return self.reply_json({"message": f"{lines:,} lines: this server builds canons of at most {LIMITS['lines']:,}"}, 400)
        if op == "moment" and not isinstance(body.get("date"), str):
            return self.reply_json({"message": "a moment needs a date (EDTF, such as 1814-11-15)"}, 400)
        state = engine.status()
        if not state["current"]:
            return self.reply_json({"message": "the engine is not built from the present sources: build it first", "engine": state}, 503)
        # A page may ask for a build and several moments at once: each waits for one of the engine's slots.
        if not SLOTS.acquire(timeout=120):
            return self.reply_json({"message": "the engine is busy with other requests; try again in a moment"}, 409)
        key = build_key(paths) if op == "build" else None
        kept = BUILDS / f"{key}.json" if key else None
        if kept is not None and kept.exists():
            SLOTS.release()
            return self.reply_text(kept.read_text(encoding="utf-8"), 200)
        toolchain = (engine.read_build() or {}).get("toolchain") or "native"
        request = {"op": op}
        if len(paths) > 1 or body.get("canon") == INTEGRATED:
            request["canons"] = [engine.path_for(toolchain, p) for p in paths]
        else:
            request["canon"] = engine.path_for(toolchain, paths[0])
        if op == "moment":
            request["date"] = body["date"]
        # The couplings a build fitted, so that a moment, the power or the field need not fit again.
        couplings = body.get("couplings")
        if isinstance(couplings, dict) and all(isinstance(k, str) and isinstance(v, (int, float)) for k, v in couplings.items()):
            request["couplings"] = couplings
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
            message = err.strip().splitlines()[-1] if err.strip() else "the engine failed"
            return self.reply_json({"message": message}, 499 if code < 0 else 500)
        if kept is not None:
            try:
                BUILDS.mkdir(parents=True, exist_ok=True)
                kept.write_text(out, encoding="utf-8")
            except OSError:
                pass
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
    parser = argparse.ArgumentParser(description="Serve the History Complex locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8870)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows building the engine and checking canons")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
        LIMITS.update(lines=20_000, seconds={"build": 30.0, "moment": 10.0, "power": 60.0, "field": 30.0})
    if not engine.status()["current"]:
        print("The engine is not built from the present sources. Build it from the page, or run:\n"
              "  python data/acquisition/recipe.py applications/history-complex/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if not public else f"{url} and on {args.host}"
    print(f"History Complex at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
