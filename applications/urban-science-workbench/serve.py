"""Serve the Urban Science Workbench on this computer.

The page is served from this directory, with the repository's ``domains/``, ``foundations/`` and
``visualization/`` under the same paths and the data cache at ``/data/``. Its plug-ins live in
``plugins/<name>/``: a manifest (``plugin.json``) naming its views, the views' modules, and, for
those with C++, an engine built by ``engine.py`` and kept running here, one process each.

    GET  /api/plugins                    the plug-ins: their manifests and whether their engines are current
    POST /api/<plugin>/<op>              {…}: one request to a plug-in's engine
    GET  /tiles/terrain/earth/<z>/<x>/<y>  a terrain tile (the data library's tiles), built on first request
    GET  /world/city/<z>/<x>/<y>         any place's buildings and ways (Overture Maps; data/acquisition/world_city.py), kept once asked
    GET  /world/base/<z>/<x>/<y>         any place's water and open land
    GET  /world/heights/<z>/<x>/<y>      the ground and the buildings' average height on a grid over a tile
    GET  /world/raster?kind=&lat=&lon=[&rp=]  a raster about a place (settlement, urbanisation, flood): {meta, data's address}
    GET  /world/space2stats?lat=&lon=    a tile of Space2Stats' cells
    GET  /world/thinkhazard?name= | code= the hazards of an administrative unit
    GET  /world/data360?q= | db=&id=     indicators by meaning, or an indicator's values by economy
    GET  /studio/scripts                 the studio's scripts in plugins/studio/scripts
    GET  /studio/script?file=            one's text          GET /studio/describe?file=   its description
    POST /studio/script                  {file, code}: keep it there
    POST /studio/run                     {file, records, params, context}: run it; returns {job}
    GET  /studio/job/<id>                its state, progress and results
    POST /studio/job/<id>/cancel         stop it
    GET  /setup                          the recipe (recipe.json): building the engines, and whether it is done
    POST /setup/run[?again=1]            build them (again)

Answers to requests without parameters are kept until the engine is rebuilt. Only 127.0.0.1 is
bound unless ``--host`` says otherwise; on a public address, building needs ``--token``.

Usage::

    python applications/urban-science-workbench/serve.py [--host 127.0.0.1] [--port 8970] [--open]

then open http://127.0.0.1:8970/ in a browser.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import http.server
import json
import re
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
sys.path.insert(0, str(ROOT / "data"))
import engine  # noqa: E402
from tiles.pyramid import TileKey  # noqa: E402
from tiles.terrain import TerrainService, decode as decode_terrain  # noqa: E402
sys.path.insert(0, str(ROOT / "data" / "acquisition"))
sys.path.insert(0, str(ROOT / "api" / "python"))
import services  # noqa: E402
import world_city  # noqa: E402
from science import spacetime  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
WORKERS: dict[str, engine.Worker] = {}
CACHE: dict[tuple[str, str], str] = {}
LOCK = threading.Lock()
TERRAIN = TerrainService(ROOT / "data" / "cache")
TILE = re.compile(r"^/tiles/terrain/earth/(\d+)/(\d+)/(\d+)$")
# Terrain tiles are built a few at a time, so that a page asking for many does not starve the rest.
BUILDS = threading.BoundedSemaphore(8)
MAX_TERRAIN_LEVEL = 18
WORLD = re.compile(r"^/world/(city|base|heights)/(\d+)/(\d+)/(\d+)$")
WORLD_TILES = threading.BoundedSemaphore(6)
# A city tile's vector data is read in pure Python, tens of thousands of footprints: apart from the server's own thread, so that the
# page goes on being answered while a place it has not seen is prepared.
READERS: concurrent.futures.ProcessPoolExecutor | None = None


def read_city_tile(z: int, x: int, y: int) -> str:
    global READERS
    with LOCK:
        if READERS is None:
            READERS = concurrent.futures.ProcessPoolExecutor(max_workers=3)
    return READERS.submit(world_city.city_tile, z, x, y).result(timeout=240)


# ---- the studio: scripts over data in space and time (ADR 0079) ---------------------------------------------------------------

SCRIPTS = APP / "plugins" / "studio" / "scripts"
SCRIPT_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}\.(py|cpp|js)$")
STUDIO_SECONDS = 600
_studio: dict = {"host": None, "jobs": {}}


def studio_host() -> spacetime.Host:
    if _studio["host"] is None:
        _studio["host"] = spacetime.Host(SCRIPTS, ROOT / "data" / "cache" / "urban-workbench" / "scripts", include_dirs=engine.include_dirs("studio"))
    return _studio["host"]


def studio_job(job: dict, body: dict) -> None:
    def progress(fraction: float, text: str) -> None:
        job["progress"], job["text"] = fraction, text
    try:
        job["results"] = studio_host().run(job["file"], {"records": body.get("records") or [], "context": body.get("context") or {}}, body.get("params") or {},
                                           progress=progress, cancel=job["cancel"], seconds=STUDIO_SECONDS)
        job["state"], job["progress"] = "done", 1.0
    except Exception as error:  # a script's failure is reported to the page, never raised in the server
        job["state"] = "cancelled" if job["cancel"].is_set() else "failed"
        job["message"] = str(error)


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
        tile = TILE.match(url.path)
        if tile:
            return self.tile_request(*map(int, tile.groups()))
        world = WORLD.match(url.path)
        if world:
            return self.world_tile(world.group(1), *map(int, world.groups()[1:]))
        if url.path.startswith("/world/"):
            return self.world_service(url.path.split("/", 2)[2], dict(urllib.parse.parse_qsl(url.query)))
        if url.path.startswith("/studio/"):
            return self.studio_get(url.path[len("/studio/"):], dict(urllib.parse.parse_qsl(url.query)))
        return super().do_GET()

    def world_tile(self, kind: str, z: int, x: int, y: int) -> None:
        if x >= 2 ** z or y >= 2 ** z or z > 14 or (kind != "heights" and z < 13):
            return self.reply_json({"message": "no such tile"}, 404)
        try:
            with WORLD_TILES:
                if kind == "city":
                    text = read_city_tile(z, x, y)
                elif kind == "base":
                    text = world_city.base_tile(z, x, y)
                else:
                    text = json.dumps(world_city.heights_of(z, x, y, TERRAIN, TileKey, decode_terrain))
        except Exception as error:  # a source that fails is reported for this tile; the place stays as it was
            return self.reply_json({"message": f"{kind} {z}/{x}/{y}: {error}"}, 502)
        self.reply_text(text, 200)

    def world_service(self, name: str, q: dict) -> None:
        try:
            if name == "space2stats":
                return self.reply_json(services.space2stats_tile(float(q["lat"]), float(q["lon"])), 200)
            if name == "thinkhazard":
                return self.reply_json({"report": services.thinkhazard_report(int(q["code"]))} if "code" in q else {"units": services.thinkhazard_search(q["name"])}, 200)
            if name == "data360":
                return self.reply_json({"values": services.data360_values(q["db"], q["id"])} if "id" in q else {"found": services.data360_search(q["q"])}, 200)
            if name == "raster":
                return self.reply_json(services.raster_window(q["kind"], float(q["lat"]), float(q["lon"]), int(q["rp"]) if "rp" in q else None), 200)
        except KeyError as error:
            return self.reply_json({"message": f"the request names no {error}"}, 400)
        except ValueError as error:
            return self.reply_json({"message": str(error)}, 400)
        except Exception as error:  # the world's services go down; the page says so where it was asked
            return self.reply_json({"message": f"the service could not answer: {error}"}, 502)
        self.reply_json({"message": f"no service {name}"}, 404)

    def tile_request(self, z: int, x: int, y: int) -> None:
        if z > MAX_TERRAIN_LEVEL:
            return self.reply_json({"message": f"no terrain below level {MAX_TERRAIN_LEVEL}"}, 404)
        key = TileKey(z, x, y)
        try:
            key.check()
            with BUILDS:
                data = TERRAIN.tile("earth", key)
        except ValueError as error:
            return self.reply_json({"message": str(error)}, 404)
        except FileNotFoundError as error:
            return self.reply_json({"message": f"the terrain's sources are not in the cache ({error})"}, 503)
        except Exception as error:  # a source that fails is reported for this tile; the page draws its parent
            return self.reply_json({"message": f"terrain {z}/{x}/{y}: {error}"}, 502)
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def studio_get(self, path: str, q: dict) -> None:
        try:
            if path == "scripts":
                files = sorted(p.name for p in SCRIPTS.iterdir() if p.is_file() and SCRIPT_NAME.match(p.name))
                return self.reply_json({"scripts": [{"file": f} for f in files], "refusal": self.studio_refusal(False)}, 200)
            if path == "script":
                return self.reply_json({"file": q["file"], "code": self.script_path(q["file"]).read_text(encoding="utf-8")}, 200)
            if path == "describe":
                return self.reply_json(studio_host().describe(q["file"]), 200)
            if path.startswith("job/"):
                job = _studio["jobs"].get(path[4:])
                if job is None:
                    return self.reply_json({"message": "no such run"}, 404)
                return self.reply_json({k: v for k, v in job.items() if k != "cancel"}, 200)
        except KeyError as error:
            return self.reply_json({"message": f"the request names no {error}"}, 400)
        except (ValueError, FileNotFoundError) as error:
            return self.reply_json({"message": str(error)}, 400)
        except Exception as error:
            return self.reply_json({"message": str(error)}, 500)
        self.reply_json({"message": "no such request"}, 404)

    @staticmethod
    def script_path(name: str) -> Path:
        if not SCRIPT_NAME.match(str(name)):
            raise ValueError(f"{name}: a script's name is letters, digits, - and _, ending .py, .cpp or .js")
        return SCRIPTS / name

    def studio_refusal(self, check_request: bool = True) -> str | None:
        """Why this server would not write or run a script for this request, or None. A script runs with the authority of whoever started the server, so that
        only the server's own page, at its own address, may ask: any other page in the same browser could otherwise post here."""
        if not WRITES["allowed"]:
            return "this server runs no scripts: start it with --token to allow it"
        if not check_request:
            return None
        if WRITES["token"] and self.headers.get("Authorization", "") != f"Bearer {WRITES['token']}":
            return "running a script on this server needs its token"
        host = self.headers.get("Host", "")
        if self.headers.get("X-Workbench") != "studio" or "application/json" not in self.headers.get("Content-Type", ""):
            return "a script is run only from the workbench's own page"
        origin = self.headers.get("Origin")
        if origin is None or urllib.parse.urlparse(origin).netloc != host:
            return "a script is run only from the workbench's own page"
        name = host.rsplit(":", 1)[0]
        if not WRITES["token"] and name not in ("127.0.0.1", "localhost", "[::1]"):
            return "a script is run only at this computer's own address"
        return None

    def studio_post(self, path: str) -> None:
        refusal = self.studio_refusal()
        if refusal:
            return self.reply_json({"message": refusal}, 403)
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > 512 * 1024 * 1024:
            return self.reply_json({"message": "the request is too large"}, 413)
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            if path == "script":
                target = self.script_path(body["file"])
                code = str(body["code"])
                SCRIPTS.mkdir(parents=True, exist_ok=True)
                if not target.exists() or target.read_text(encoding="utf-8") != code:
                    target.write_text(code, encoding="utf-8")
                return self.reply_json({"file": target.name}, 200)
            if path == "run":
                target = self.script_path(body["file"])
                if target.suffix == ".js":
                    raise ValueError("a script in JavaScript runs in the page")
                ident = f"{int(time.time() * 1000):x}"
                job = {"id": ident, "file": target.name, "state": "running", "progress": 0.0, "text": "", "cancel": threading.Event()}
                _studio["jobs"][ident] = job
                threading.Thread(target=studio_job, args=(job, body), daemon=True).start()
                return self.reply_json({"job": ident}, 200)
            if path.startswith("job/") and path.endswith("/cancel"):
                job = _studio["jobs"].get(path[4:-7])
                if job is None:
                    return self.reply_json({"message": "no such run"}, 404)
                job["cancel"].set()
                return self.reply_json({"cancelled": job["id"]}, 200)
        except KeyError as error:
            return self.reply_json({"message": f"the request names no {error}"}, 400)
        except ValueError as error:
            return self.reply_json({"message": str(error)}, 400)
        self.reply_json({"message": "no such request"}, 404)

    def do_POST(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path.startswith("/studio/"):
            return self.studio_post(url.path[len("/studio/"):])
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
        if length > 96 * 1024 * 1024:
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
        if not TILE.match(urllib.parse.urlparse(self.path).path):
            self.send_header("Cache-Control", "no-cache")
        else:
            self.send_header("Cache-Control", "public, max-age=86400")
        super().end_headers()

    def log_message(self, format: str, *args) -> None:
        if args and str(args[1]).startswith(("4", "5")):
            super().log_message(format, *args)


class Server(http.server.ThreadingHTTPServer):
    request_queue_size = 64
    daemon_threads = True


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Serve the Urban Science Workbench locally.")
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
              "  python data/acquisition/recipe.py applications/urban-science-workbench/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    threading.Thread(target=warm, daemon=True).start()
    url = f"http://127.0.0.1:{args.port}/"
    print(f"Urban Science Workbench at {url if not public else url + ' and on ' + args.host}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        for w in WORKERS.values():
            w.stop()
        if READERS is not None:
            READERS.shutdown(wait=False, cancel_futures=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
