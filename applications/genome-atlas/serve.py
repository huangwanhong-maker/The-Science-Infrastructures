"""Serve the Genome Atlas on this computer.

The page is served from this directory, and the code it shares from the repository's
``domains/``, ``foundations/`` and ``visualization/`` under the same paths, with the data cache
at ``/data/``. The genomes it shows are NCBI's, fetched by the data library into
``data/cache/ncbi-genomes/`` (``data/registry/ncbi-genomes.json``), and genomes simulated from a
seed. Everything is computed by the engine (``engine.py``: the genomics and evolution domains'
C++), kept alive as workers by the cluster (``cluster.py``). The protocol is in the README, §3:

    GET  /api/engine            the engine's build and the cluster's workers, with their targets
    GET  /api/genomes           the genomes there are: NCBI's (fetched or not) and the simulated presets
    POST /api/run               {"op", "genome", ...}: one operation on one worker
    POST /api/job               {"op", "genome", ...}: a job cut into shards across the workers
    POST /api/restart           stop the workers (they start again, from the present build, when next needed)

The engine and the genomes are made by the application's recipe (``recipe.json``):

    GET  /setup                 its steps, and whether each is done
    POST /setup/run[?again=1]   run them

Only 127.0.0.1 is bound unless ``--host`` says otherwise. On a public address, running the recipe
and restarting need the server's token (``--token``), and requests are held within smaller limits.

Usage::

    python applications/genome-atlas/serve.py [--host 127.0.0.1] [--port 8890] [--workers 2] [--open]

then open http://127.0.0.1:8890/ in a browser.
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
    ".md": "text/markdown; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".wgsl": "text/plain; charset=utf-8",
}
MANIFEST = ROOT / "data" / "registry" / "ncbi-genomes.json"
GENOMES = ROOT / "data" / "cache" / "ncbi-genomes"

sys.path.insert(0, str(APP))
import cluster  # noqa: E402
import engine  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
LIMITS = {"seconds": 300.0, "job_seconds": 1800.0, "body": 16 * 1024 * 1024}
CLUSTER: cluster.Cluster | None = None

# Genomes simulated from a seed: what the atlas shows without the network, with known truth.
SIMULATED = {
    "sim:asm-30k": {"title": "A small genome to assemble, 30 kb, with two repeats of 600 bp and a tandem array",
                    "spec": {"length": 30_000, "gc": 0.5, "seed": 21,
                             "repeats": {"elements": 2, "copies": 3, "length": 600, "divergence": 0.0},
                             "tandem": {"arrays": 1, "unit": 40, "copies": 20, "divergence": 0.0}}},
    "sim:plain-100k": {"title": "A random genome, 100 kb, GC 50%", "spec": {"length": 100_000, "gc": 0.5, "seed": 1}},
    "sim:repeats-200k": {"title": "A genome with repeats, 200 kb: interspersed copies of 3 kb and tandem arrays",
                         "spec": {"length": 200_000, "gc": 0.52, "seed": 7,
                                  "repeats": {"elements": 3, "copies": 6, "length": 3000, "divergence": 0.01},
                                  "tandem": {"arrays": 3, "unit": 37, "copies": 40, "divergence": 0.03}}},
    "sim:gc-rich-1m": {"title": "A GC-rich genome, 1 Mb, GC 66%, with repeats",
                       "spec": {"length": 1_000_000, "gc": 0.66, "seed": 11,
                                "repeats": {"elements": 8, "copies": 4, "length": 5000, "divergence": 0.005}}},
    "sim:bacterium-5m": {"title": "A bacterium-sized genome, 5 Mb, GC 51%",
                         "spec": {"length": 5_000_000, "gc": 0.51, "seed": 13,
                                  "repeats": {"elements": 7, "copies": 7, "length": 5500, "divergence": 0.002}}},
}


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


def write_refusal(handler) -> str | None:
    if not WRITES["allowed"]:
        return "this server does not build or restart: start it with --token to allow it with a token"
    token = WRITES["token"]
    if token and handler.headers.get("Authorization", "") != f"Bearer {token}":
        return "building and restarting on this server need its token"
    return None


def genome_catalogue() -> list[dict]:
    """NCBI's genomes of the manifest (with whether each is fetched) and the simulated presets."""
    out = []
    try:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        for f in manifest.get("files", []):
            acc = f["file"].removesuffix(".gb")
            title = f.get("title", acc)
            group = "mitochondria" if "mitochondri" in title.lower() else "viruses" if any(w in title.lower() for w in ("phage", "virus", "sars")) else "bacteria"
            out.append({"id": f"ncbi:{acc}", "title": title, "accession": acc, "group": group, "fetched": (GENOMES / f["file"]).exists(),
                        "size": f.get("size")})
    except (OSError, ValueError):
        pass
    for gid, s in SIMULATED.items():
        out.append({"id": gid, "title": s["title"], "group": "simulated", "fetched": True, "spec": s["spec"]})
    return out


def genome_spec(value, toolchain: str) -> dict:
    """The engine's spec of the genome a request names."""
    if isinstance(value, dict):
        if "simulated" in value and isinstance(value["simulated"], dict):
            return {"simulated": value["simulated"]}
        if "text" in value and isinstance(value["text"], str):
            return {"text": value["text"], "name": str(value.get("name") or "pasted")}
        raise ValueError("a genome is an id, {simulated: spec} or {text, name}")
    if not isinstance(value, str):
        raise ValueError("the request names no genome")
    if value.startswith("ncbi:"):
        acc = value[5:]
        path = GENOMES / f"{acc}.gb"
        known = {g["id"] for g in genome_catalogue()}
        if value not in known:
            raise ValueError(f"no genome {acc} in the manifest")
        if not path.exists():
            raise ValueError(f"{acc} is not fetched yet: run the recipe (Setup)")
        return {"file": engine.path_for(toolchain, path)}
    if value in SIMULATED:
        return {"simulated": SIMULATED[value]["spec"]}
    raise ValueError(f"no genome {value!r}")


def translate(body: dict) -> dict:
    """The request as the engine takes it: genomes named by id become specs; so do lists of them."""
    toolchain = (engine.read_build() or {}).get("toolchain") or "native"
    out = dict(body)
    if "genome" in body:
        out["genome"] = genome_spec(body["genome"], toolchain)
    if isinstance(body.get("genomes"), list):
        out["genomes"] = [{"id": g if isinstance(g, str) else None, **genome_spec(g, toolchain)} for g in body["genomes"]]
    return out


def local_only(request: dict) -> bool:
    """Whether a request reads this computer's files (and so stays with the local workers)."""
    specs = [request.get("genome")] + list(request.get("genomes") or [])
    return any(isinstance(s, dict) and "file" in s for s in specs)


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
            return self.reply_json({**engine.status(), "cluster": CLUSTER.describe(), "limits": LIMITS}, 200)
        if url.path == "/api/genomes":
            return self.reply_json({"genomes": genome_catalogue()}, 200)
        if url.path == "/setup":
            return self.reply_json(*self.recipe_request(url))
        return super().do_GET()

    def do_POST(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/setup/run":
            return self.reply_json(*self.recipe_request(url))
        if url.path == "/api/restart":
            refusal = write_refusal(self)
            if refusal:
                return self.reply_json({"message": refusal}, 403)
            CLUSTER.restart()
            return self.reply_json({"restarted": True}, 200)
        if url.path in ("/api/run", "/api/job"):
            return self.engine_request(job=url.path == "/api/job")
        self.send_error(404)

    def read_body(self) -> dict | None:
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > LIMITS["body"]:
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

    def engine_request(self, job: bool) -> None:
        body = self.read_body()
        if body is None:
            return
        op = body.get("op")
        if not isinstance(op, str) or not op:
            return self.reply_json({"message": "a request names its op"}, 400)
        state = engine.status()
        if not state["built"]:
            return self.reply_json({"message": "the engine is not built: build it first (Setup)", "engine": state}, 503)
        # A stale engine (built from older sources) still answers; the page says so.
        try:
            request = translate(body)
        except ValueError as error:
            return self.reply_json({"message": str(error)}, 400)
        try:
            if job:
                result = CLUSTER.job(op, request, timeout=LIMITS["job_seconds"], local_only=local_only(request))
            else:
                result = CLUSTER.request(request, timeout=LIMITS["seconds"])
        except cluster.WorkerError as error:
            return self.reply_json({"message": str(error)}, 400)
        self.reply_json(result if isinstance(result, dict) else {"result": result}, 200)

    def recipe_request(self, url) -> tuple[dict, int]:
        refusal = write_refusal(self)
        policy = {"writable": refusal is None, "refusal": refusal}
        try:
            if url.path == "/setup":
                return {**recipe_job().view(), **policy}, 200
            if refusal:
                return {"message": refusal, **policy}, 403
            query = dict(urllib.parse.parse_qsl(url.query))
            CLUSTER.restart()  # a new build replaces the program the workers run
            return {**recipe_job().start(None, again=query.get("again") == "1"), **policy}, 200
        except RuntimeError as error:
            return {"message": str(error)}, 409

    def reply_json(self, body: dict, code: int) -> None:
        data = json.dumps(body, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
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
    global CLUSTER
    parser = argparse.ArgumentParser(description="Serve the Genome Atlas locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8890)
    parser.add_argument("--workers", type=int, default=2, help="local engine workers (each uses every hardware thread)")
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows running the recipe and restarting")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
        LIMITS.update(seconds=60.0, job_seconds=300.0, body=2 * 1024 * 1024)
    CLUSTER = cluster.Cluster(local=args.workers)
    if not engine.status()["current"]:
        print("The engine is not built from the present sources. Build it from the page (Setup), or run:\n"
              "  python data/acquisition/recipe.py applications/genome-atlas/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    print(f"Genome Atlas at {url if not public else f'{url} and on {args.host}'}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        CLUSTER.stop()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
