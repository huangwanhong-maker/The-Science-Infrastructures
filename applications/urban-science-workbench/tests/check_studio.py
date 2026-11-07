"""Checks of the studio's server: who may write and run a script, and that scripts in Python and C++ answer what they should. Each line is "ok  name" or "FAIL  name: why",
read by run-tests.mjs. A script runs with the authority of whoever started the server, so that the refusals here are the point of the file.
"""

from __future__ import annotations

import json
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))
import serve  # noqa: E402

failures = 0


def check(name: str, ok: bool, why: str = "") -> None:
    global failures
    if not ok:
        failures += 1
    print(f"{'ok  ' if ok else 'FAIL'}  {name}{': ' + str(why)[:300] if why else ''}", flush=True)


server = serve.Server(("127.0.0.1", 0), serve.Handler)
threading.Thread(target=server.serve_forever, daemon=True).start()
PORT = server.server_address[1]
BASE = f"http://127.0.0.1:{PORT}"
OWN = {"Content-Type": "application/json", "X-Workbench": "studio", "Origin": BASE}
NAME = "check-studio-script.py"


def post(path: str, body: dict, headers: dict) -> tuple[int, dict]:
    request = urllib.request.Request(BASE + path, data=json.dumps(body).encode(), headers=headers, method="POST")
    try:
        with urllib.request.urlopen(request, timeout=300) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read() or b"{}")


def get(path: str) -> tuple[int, dict]:
    try:
        with urllib.request.urlopen(BASE + path, timeout=300) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read() or b"{}")


def wait(job: str) -> dict:
    for _ in range(1200):
        _, a = get(f"/studio/job/{job}")
        if a.get("state") != "running":
            return a
        time.sleep(0.25)
    return {"state": "timed out"}


CODE = '''
import numpy as np
def create(api):
    class Twice(api.SpaceTimeAnalysis):
        id = "twice"; title = "Twice"; description = "Doubles a value."
        params = [api.param("by", "By", 2.0)]
        def run(self, records, params, context):
            r = records[0]
            context.progress(0.5, "halfway")
            return [api.events("Doubled", r.time, r.lat, r.lon, values={"v": r.values["v"] * params["by"]}), api.note(f"looking at {context.lat}")]
    return Twice()
'''
target = serve.SCRIPTS / NAME
try:
    # ---- who may ask -----------------------------------------------------------------------------------------------------------------------------
    body = {"file": NAME, "code": CODE}
    code, _ = post("/studio/script", body, {"Content-Type": "application/json"})
    check("a page that does not say it is the workbench is refused", code == 403 and not target.exists(), code)
    code, _ = post("/studio/script", body, {**OWN, "Origin": "https://example.org"})
    check("a request from another site is refused, whatever it says", code == 403 and not target.exists(), code)
    code, _ = post("/studio/script", body, {"Content-Type": "application/json", "X-Workbench": "studio"})
    check("a request that names no origin is refused", code == 403 and not target.exists(), code)
    code, _ = post("/studio/script", body, {**OWN, "Content-Type": "text/plain"})
    check("a form's content type, which needs no permission to send, is refused", code == 403 and not target.exists(), code)
    code, _ = post("/studio/run", {"file": "density.py", "records": []}, {**OWN, "Origin": "https://example.org"})
    check("nor can another site run a script already kept", code == 403, code)
    # ---- what may be named -----------------------------------------------------------------------------------------------------------------------
    for bad in ("../../evil.py", "..\\evil.py", "evil.exe", ".hidden.py", "a/b.py"):
        code, _ = post("/studio/script", {"file": bad, "code": "x"}, OWN)
        check(f"a script cannot be named {bad}", code == 400, code)
    code, a = get("/studio/script?file=../../serve.py")
    check("nor can a file outside be read through it", code == 400, code)
    # ---- what it does ------------------------------------------------------------------------------------------------------------------------------
    code, a = post("/studio/script", body, OWN)
    check("the workbench's own page may keep a script", code == 200 and target.exists(), a)
    code, a = get(f"/studio/describe?file={NAME}")
    check("it is described", code == 200 and a.get("id") == "twice" and a.get("kind") == "spacetime-analysis" and a.get("language") == "python", a)
    record = {"kind": "events", "title": "t", "time": [1e12, 1e12 + 1000], "lat": [1, 2], "lon": [3, 4], "values": {"v": [10, None]}}
    code, a = post("/studio/run", {"file": NAME, "records": [record], "params": {"by": 3}, "context": {"lat": 5.5}}, OWN)
    done = wait(a.get("job", ""))
    results = done.get("results") or []
    check("a Python script runs, with its parameter and the view", done.get("state") == "done" and results[0]["values"]["v"][0] == 30 and results[1]["text"] == "looking at 5.5", done.get("message") or results)
    check("a missing value stays missing through it", done.get("state") == "done" and results[0]["values"]["v"][1] is None, results)
    (serve.SCRIPTS / NAME).write_text("def create(api):\n    raise ValueError('it failed on purpose')\n", encoding="utf-8")
    code, a = post("/studio/run", {"file": NAME, "records": [record]}, OWN)
    done = wait(a.get("job", ""))
    check("a script's failure comes back with its reason, and the server stands", done.get("state") == "failed" and "on purpose" in (done.get("message") or ""), done.get("message"))
    # The examples, on events of a known count.
    n = 400
    events = {"kind": "events", "title": "ring", "time": [1.7e12 + k * 3.6e6 for k in range(n)], "lat": [10 + (k % 20) * 0.5 for k in range(n)], "lon": [100 + (k % 30) * 0.5 for k in range(n)], "values": {}}
    code, a = post("/studio/run", {"file": "hex_bins.cpp", "records": [events], "params": {"resolution": 3, "steps": 8}}, OWN)
    done = wait(a.get("job", ""))
    cells = (done.get("results") or [{}])[0]
    total = sum(sum(f) for f in cells.get("values", {}).get("events", [])) if cells else 0
    check("the C++ example compiles and runs, and every event is counted once", done.get("state") == "done" and total == n, done.get("message") or total)
    check("and its cells have rings, and a step for each time", done.get("state") == "done" and all(len(c["ring"]) >= 5 for c in cells["cells"]) and len(cells["time"]) == 8, "")
    code, a = post("/studio/run", {"file": "density.py", "records": [events], "params": {"bandwidth_km": 200, "steps": 6, "window": 1.0, "width": 90}}, OWN)
    done = wait(a.get("job", ""))
    grid = (done.get("results") or [{}])[0]
    check("the Python example answers a grid of the size it says", done.get("state") == "done" and len(grid["values"]["events per km²"][0]) == grid["width"] * grid["height"], done.get("message"))
    # Its last frame's integral over the grid is the count of events, up to what its kernels spread beyond the grid and what it cut below a thousandth.
    if done.get("state") == "done":
        import math
        south, west, north, east = grid["box"]
        w, h = grid["width"], grid["height"]
        last = grid["values"]["events per km²"][-1]
        area = 0.0
        for j in range(h):
            lat = south + (j + 0.5) * (north - south) / h
            dy, dx = (north - south) / h * 111.0, (east - west) / w * 111.0 * math.cos(math.radians(lat))
            area += sum(v for v in last[j * w:(j + 1) * w] if v is not None) * dx * dy
        check("and the integral of its surface is near the count of events in its window", abs(area / n - 1) < 0.1, f"{area:.1f} of {n}")
finally:
    target.unlink(missing_ok=True)
    server.shutdown()
print(f"\n{failures} failed", flush=True)
sys.exit(1 if failures else 0)
