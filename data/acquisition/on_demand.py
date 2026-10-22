"""Requests made on demand: a dataset fetched with chosen parameters and built
into a product while an application waits.

Most data are fetched ahead of time and listed in the index of products. Some
cannot be: the ocean beneath a place a user has just chosen is one column out
of millions, and no one fetches them all. A manifest offers itself for such
requests with an `on_demand` block:

    "on_demand": {
      "parameters": ["lat", "lon"],     the parameters a request may set
      "builder": "ocean_column",        a module of data/products with build(folder, manifest) -> path
      "seconds": 240                    how long a request usually takes, to tell the user
    }

A request names the dataset and values for those parameters, which are checked
as any parameter is (type and bounds, from the manifest). The retrieval and the
product are cached like any other: the same request is answered from the cache
while it is current. Requests run in threads, two at a time, so that the
servers asked are not crowded; identical requests share one job.

The laboratory's local servers expose this to their pages (for example
applications/geography/earth-atlas/serve.py, at /request); nothing here is
reachable from outside the computer.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib
import io
import json
import sys
import threading
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(DATA / "products"))
import fetch  # noqa: E402

FRESH_SECONDS = 3600          # a finished job is reused for this long, then asked again (often from the cache)


class RequestError(ValueError):
    pass


class _PerThread(io.TextIOBase):
    """Standard output that a thread may divert to its own buffer (a job's log),
    while every other thread keeps writing to the console."""

    def __init__(self, base):
        self.base, self.local = base, threading.local()

    def write(self, s: str) -> int:
        target = getattr(self.local, "buffer", None) or self.base
        return target.write(s)

    def flush(self) -> None:
        self.base.flush()

    @contextlib.contextmanager
    def capture(self, buffer: io.StringIO):
        self.local.buffer = buffer
        try:
            yield buffer
        finally:
            self.local.buffer = None


def offered() -> dict:
    """The datasets that accept requests on demand: id -> manifest."""
    out = {}
    for path in sorted(fetch.REGISTRY.glob("*.json")):
        manifest = fetch.load(path)
        if manifest.get("on_demand"):
            out[manifest["id"]] = manifest
    return out


class Jobs:
    """Jobs by id; each is a dict that the server returns as JSON."""

    def __init__(self, at_once: int = 2):
        self.jobs: dict[str, dict] = {}
        self.lock = threading.Lock()
        self.slots = threading.Semaphore(at_once)
        if not isinstance(sys.stdout, _PerThread):
            sys.stdout = _PerThread(sys.stdout)
        self.output = sys.stdout

    def request(self, dataset: str, parameters: dict) -> dict:
        manifests = offered()
        if dataset not in manifests:
            raise RequestError(f"no dataset {dataset!r} accepts requests on demand")
        manifest = manifests[dataset]
        spec = manifest["on_demand"]
        unknown = set(parameters) - set(spec["parameters"])
        if unknown:
            raise RequestError(f"{dataset}: parameters {', '.join(sorted(unknown))} cannot be set on demand")
        try:
            values = {k: fetch.coerce(k, manifest["parameters"][k], v) for k, v in parameters.items()}
        except SystemExit as error:          # coerce() reports bad values as the command line does
            raise RequestError(str(error)) from None
        key = json.dumps([dataset, sorted(values.items())])
        ident = hashlib.sha256(key.encode()).hexdigest()[:16]
        with self.lock:
            job = self.jobs.get(ident)
            stale = job and job["state"] == "done" and time.time() - job["finished"] > FRESH_SECONDS
            if job and job["state"] != "failed" and not stale:
                return self.view(job)
            job = {"id": ident, "dataset": dataset, "parameters": values, "state": "queued",
                   "seconds": spec.get("seconds"), "created": time.time(), "started": None, "finished": None,
                   "product": None, "message": spec.get("description", "")}
            self.jobs[ident] = job
        threading.Thread(target=self._run, args=(job, manifest), daemon=True).start()
        return self.view(job)

    def status(self, ident: str) -> dict:
        with self.lock:
            job = self.jobs.get(ident)
        if not job:
            raise RequestError(f"no request {ident!r}")
        return self.view(job)

    @staticmethod
    def view(job: dict) -> dict:
        out = dict(job)
        now = time.time()
        if job["started"]:
            out["elapsed"] = round((job["finished"] or now) - job["started"], 1)
        return out

    def _run(self, job: dict, manifest: dict) -> None:
        with self.slots:
            job["state"], job["started"] = "running", time.time()
            log = io.StringIO()
            try:
                with self.output.capture(log):
                    kind = manifest.get("kind", "files")
                    overrides = {k: str(v) for k, v in job["parameters"].items()}
                    if not fetch.ACQUIRE[kind](manifest, None, overrides, False):
                        raise RuntimeError(log.getvalue().strip().splitlines()[-1] if log.getvalue().strip() else "the retrieval failed")
                    _, tag = fetch.resolve(manifest, None, overrides)
                    builder = importlib.import_module(manifest["on_demand"]["builder"])
                    path = builder.build(fetch.CACHE / manifest["id"] / tag, manifest)
                job["product"] = Path(path).resolve().relative_to(fetch.CACHE.resolve()).as_posix()
                job["state"] = "done"
            except Exception as error:      # reported to the page, which shows it
                job["state"] = "failed"
                job["message"] = f"{type(error).__name__}: {error}".strip()
                traceback.print_exc()
            finally:
                job["finished"] = time.time()
                job["log"] = log.getvalue()[-2000:]
