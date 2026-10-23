"""Serve the Spiking Simulator on this computer.

The page is served from this directory, and the code it shares from the repository's
``domains/``, ``foundations/`` and ``visualization/`` under the same paths. The networks it
composes are run by the engine (``engine.py``: the neuroscience domain's C++ simulator):

    GET  /api/engine                     whether the engine is built and current, and the limits of a run
    GET  /api/models                     the catalogue of models, from the engine
    POST /api/run                        run the specification in the body; its results, as JSON
    POST /api/cancel                     stop the run in progress

The engine is built by the application's recipe (``recipe.json``), which the page may run:

    GET  /setup                          its one step, and whether it is done
    POST /setup/run[?again=1]            build the engine (again)

Only 127.0.0.1 is bound unless ``--host`` says otherwise. On a public address, building
needs the server's token (``--token``), and runs are held within smaller limits.

Usage::

    python applications/spiking-simulator/serve.py [--host 127.0.0.1] [--port 8830] [--open]

then open http://127.0.0.1:8830/ in a browser.
"""

from __future__ import annotations

import argparse
import http.server
import json
import os
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
# The limits of a run: the neurons of a network, its simulated duration, the wall time.
LIMITS = {"neurons": 250_000, "duration": 60.0, "seconds": 900.0, "threads": os.cpu_count() or 1}
RUNNING = {"run": None}
RUN_LOCK = threading.Lock()


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


def seconds_of(q) -> float:
    """A duration written as the specification writes it ("2 s", "500 ms", "1.5"); NaN if not."""
    if isinstance(q, (int, float)):
        return float(q)
    try:
        number, unit = str(q).strip().split()
        return float(number) * {"s": 1, "ms": 1e-3, "us": 1e-6, "min": 60, "h": 3600}[unit]
    except (ValueError, KeyError):
        return float("nan")


def refuse_beyond_limits(spec: dict) -> str | None:
    """The reason a specification exceeds this server's limits, or None. The engine
    checks everything else."""
    try:
        neurons = sum(int(p.get("size", 0)) for p in spec.get("populations", []))
    except (TypeError, ValueError):
        return "every population needs a whole number as its size"
    if neurons > LIMITS["neurons"]:
        return f"{neurons:,} neurons: this server runs at most {LIMITS['neurons']:,}"
    duration = seconds_of(spec.get("duration", 0))
    if duration == duration and duration > LIMITS["duration"]:
        return f"{duration:g} s of activity: this server simulates at most {LIMITS['duration']:g} s in a run"
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
            return self.reply_json({**engine.status(), "limits": LIMITS, "running": RUNNING["run"] is not None}, 200)
        if url.path == "/api/models":
            m = engine.models()
            return self.reply_json({"models": m} if m is not None else {"message": "the engine is not built"}, 200 if m is not None else 503)
        if url.path == "/setup":
            return self.reply_json(*self.recipe_request(url))
        return super().do_GET()

    def do_POST(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/setup/run":
            return self.reply_json(*self.recipe_request(url))
        if url.path == "/api/cancel":
            run = RUNNING["run"]
            if run:
                run.cancel()
            return self.reply_json({"cancelled": run is not None}, 200)
        if url.path == "/api/run":
            return self.run_request()
        self.send_error(404)

    def run_request(self) -> None:
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > 32 * 1024 * 1024:
            return self.reply_json({"message": "the specification is too large"}, 413)
        try:
            spec = json.loads(self.rfile.read(length).decode("utf-8"))
            if not isinstance(spec, dict):
                raise ValueError("a specification is a JSON object")
        except ValueError as error:
            return self.reply_json({"message": f"not a specification: {error}"}, 400)
        reason = refuse_beyond_limits(spec)
        if reason:
            return self.reply_json({"message": reason}, 400)
        state = engine.status()
        if not state["current"]:
            return self.reply_json({"message": "the engine is not built from the present sources: build it first", "engine": state}, 503)
        if not RUN_LOCK.acquire(blocking=False):
            return self.reply_json({"message": "another run is in progress"}, 409)
        try:
            threads = min(int(spec.pop("threads", LIMITS["threads"]) or 1), LIMITS["threads"])
            run = engine.Run(json.dumps(spec), threads)
            RUNNING["run"] = run
            code, out, err = run.execute(timeout=LIMITS["seconds"])
        finally:
            RUNNING["run"] = None
            RUN_LOCK.release()
        if code != 0:
            message = err.strip().splitlines()[-1] if err.strip() else "the run failed"
            return self.reply_text(json.dumps({"message": message.removeprefix("snn: ")}), 400 if code > 0 else 499)
        # The results are passed through as the engine wrote them.
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
    parser = argparse.ArgumentParser(description="Serve the Spiking Simulator locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8830)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows building the engine")
    parser.add_argument("--max-neurons", type=int, help="the most neurons in a run")
    parser.add_argument("--max-duration", type=float, help="the most simulated seconds in a run")
    parser.add_argument("--max-seconds", type=float, help="the most wall-clock seconds a run may take")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
        LIMITS.update(neurons=50_000, duration=10.0, seconds=60.0, threads=min(LIMITS["threads"], 4))
    for key, value in (("neurons", args.max_neurons), ("duration", args.max_duration), ("seconds", args.max_seconds)):
        if value is not None:
            LIMITS[key] = value
    if not engine.status()["current"]:
        print("The engine is not built from the present sources. Build it from the page, or run:\n"
              "  python data/acquisition/recipe.py applications/spiking-simulator/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if not public else f"{url} and on {args.host}"
    print(f"Spiking Simulator at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
