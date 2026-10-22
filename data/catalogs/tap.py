"""A minimal client for the IVOA Table Access Protocol (TAP), standard library only.

TAP is the protocol that astronomical archives (ESA Gaia, CDS VizieR, NASA
Exoplanet Archive, MAST, ...) use to answer ADQL queries. This client runs
asynchronous (UWS) jobs, which suit large results, and can also run
synchronous queries for small ones.

Reference: IVOA, "Table Access Protocol, Version 1.1" (2019); "Universal Worker
Service Pattern, Version 1.1" (2016).
"""

from __future__ import annotations

import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

USER_AGENT = "laboratory-data-tap/1"
CHUNK = 1 << 20


class TapError(RuntimeError):
    pass


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


_opener = urllib.request.build_opener(_NoRedirect)


def _request(url: str, data: dict | None = None) -> urllib.request.Request:
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    return urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT})


def _read_text(url: str) -> str:
    with urllib.request.urlopen(_request(url), timeout=60) as response:
        return response.read().decode("utf-8", "replace").strip()


def _stream(response, dest: Path, progress: bool) -> int:
    part = dest.with_name(dest.name + ".part")
    done = 0
    with part.open("wb") as out:
        for block in iter(lambda: response.read(CHUNK), b""):
            out.write(block)
            done += len(block)
            if progress:
                print(f"\r  received {done / 1e6:7.1f} MB", end="", flush=True)
    if progress:
        print()
    part.replace(dest)
    return done


def run_async(service: str, adql: str, dest: Path, fmt: str = "csv", maxrec: int | None = None,
              poll_seconds: float = 3.0, timeout_seconds: float = 3600.0, progress: bool = True) -> int:
    """Runs an ADQL query as an asynchronous job and writes the result to dest.

    maxrec asks the service for up to that many rows; without it, many services
    apply a default limit smaller than the query's TOP. Returns the number of
    bytes written. Raises TapError if the job fails.
    """
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt, "PHASE": "RUN", "QUERY": adql}
    if maxrec is not None:
        params["MAXREC"] = str(maxrec)
    try:
        _opener.open(_request(f"{service.rstrip('/')}/async", params), timeout=60)
        raise TapError("the service did not redirect to a job")
    except urllib.error.HTTPError as error:
        if error.code not in (301, 302, 303):
            raise TapError(f"job submission failed: HTTP {error.code}: {error.read()[:500]!r}") from None
        job = error.headers["Location"]
    if progress:
        print(f"  job {job}")

    started = time.monotonic()
    while True:
        phase = _read_text(f"{job}/phase")
        if phase == "COMPLETED":
            break
        if phase in ("ERROR", "ABORTED"):
            try:
                detail = _read_text(f"{job}/error")
            except urllib.error.URLError:
                detail = ""
            raise TapError(f"job {phase.lower()}: {detail[:1000]}")
        if time.monotonic() - started > timeout_seconds:
            raise TapError(f"job still {phase.lower()} after {timeout_seconds:.0f} s: {job}")
        if progress:
            print(f"\r  {phase.lower():12s} {time.monotonic() - started:6.0f} s", end="", flush=True)
        time.sleep(poll_seconds)
    if progress:
        print(f"\r  completed in {time.monotonic() - started:.0f} s" + " " * 20)

    with urllib.request.urlopen(_request(f"{job}/results/result"), timeout=600) as response:
        return _stream(response, dest, progress)


def run_sync(service: str, adql: str, dest: Path, fmt: str = "csv", progress: bool = True) -> int:
    """Runs a small ADQL query synchronously and writes the result to dest."""
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt, "QUERY": adql}
    try:
        with urllib.request.urlopen(_request(f"{service.rstrip('/')}/sync", params), timeout=600) as response:
            return _stream(response, dest, progress)
    except urllib.error.HTTPError as error:
        raise TapError(f"query failed: HTTP {error.code}: {error.read()[:1000]!r}") from None
