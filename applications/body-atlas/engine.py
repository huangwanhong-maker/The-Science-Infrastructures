"""The Body Atlas's engine, built and kept running.

The engine (``engine/main.cpp``) is header-only C++20 on the anatomy domain, the foundations (morphe, calx,
physis) and the API's JSON (``api/cpp/science/json.hpp``). It is compiled here by g++ 11 or later, natively
where there is one and otherwise in WSL, into ``data/cache/engines/body/``, beside a ``build.json`` holding
the hash of every file it was built from (its source and the headers it includes, followed transitively), so
that work elsewhere in the foundations does not make it stale. Under WSL the build takes the lock that the
sessions working on this repository share (``%TEMP%/science-wsl.lock``).

A built engine speaks JSON lines (``--serve``); ``Worker`` keeps one process and passes it one request at a
time, telling it where the data library's products are (``--data``).

Usage::

    python applications/body-atlas/engine.py build [--force]
    python applications/body-atlas/engine.py status
    python applications/body-atlas/engine.py ask '<request json>'
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
SOURCE = APP / "engine" / "main.cpp"
OUT = ROOT / "data" / "cache" / "engines" / "body"
PROGRAM = OUT / "body"
PRODUCTS = ROOT / "data" / "cache" / "products"
INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]', re.M)


def include_dirs() -> list[Path]:
    dirs = [ROOT / "foundations" / "morphe" / "include"]
    for pattern in ("foundations/morphe/*/include", "foundations/calx/numerics/*/include", "foundations/calx/compute/*/include",
                    "foundations/calx/compute/backends/*/include", "foundations/physis/*/include", "domains/anatomy/*/include"):
        dirs += sorted(p for p in ROOT.glob(pattern) if p.is_dir())
    dirs += [ROOT / "api" / "cpp", APP / "engine"]
    return [d for d in dirs if d.is_dir()]


def sources() -> list[Path]:
    """The engine's source and the project headers it includes, transitively."""
    dirs = include_dirs()
    seen: set[Path] = set()
    todo = [SOURCE]
    while todo:
        f = todo.pop()
        if f in seen:
            continue
        seen.add(f)
        for inc in INCLUDE.findall(f.read_text(encoding="utf-8", errors="replace")):
            for d in [f.parent, *dirs]:
                candidate = d / inc
                if candidate.is_file():
                    todo.append(candidate.resolve())
                    break
    return sorted(seen)


def digest() -> str:
    h = hashlib.sha256()
    for f in sources():
        h.update(f.relative_to(ROOT).as_posix().encode())
        h.update(b"\0")
        h.update(f.read_bytes())
    return h.hexdigest()


def wsl_path(p: Path) -> str:
    s = str(p.resolve())
    if len(s) > 1 and s[1] == ":":
        return f"/mnt/{s[0].lower()}{s[2:]}".replace("\\", "/")
    return s.replace("\\", "/")


def toolchain() -> str | None:
    """'native' where g++ is on the PATH (outside Windows), 'wsl' where it is in WSL, else None."""
    if shutil.which("g++") and os.name != "nt":
        return "native"
    if os.name == "nt" and shutil.which("wsl"):
        try:
            r = subprocess.run(["wsl", "--exec", "sh", "-c", "command -v g++"], capture_output=True, timeout=60)
            if r.returncode == 0 and r.stdout.strip():
                return "wsl"
        except (OSError, subprocess.TimeoutExpired):
            pass
    return "native" if shutil.which("g++") else None


def command(kind: str, args: list[str]) -> list[str]:
    return ["wsl", "--exec", *args] if kind == "wsl" else args


def path_for(kind: str, p: Path) -> str:
    return wsl_path(p) if kind == "wsl" else str(p)


@contextlib.contextmanager
def wsl_lock(kind: str, log=print, wait: float = 1800):
    """The sessions' shared lock on WSL, held while compiling there."""
    if kind != "wsl":
        yield
        return
    lock = Path(tempfile.gettempdir()) / "science-wsl.lock"
    t0, told = time.time(), False
    while True:
        try:
            lock.mkdir()
            break
        except FileExistsError:
            if not told:
                owner = (lock / "owner").read_text(encoding="utf-8", errors="replace").strip() if (lock / "owner").exists() else "another session"
                log(f"Waiting for WSL, held by {owner}.")
                told = True
            if time.time() - t0 > wait:
                raise RuntimeError("WSL stayed busy; try again later")
            time.sleep(5)
    try:
        (lock / "owner").write_text(f"body-atlas build {time.strftime('%Y-%m-%dT%H:%M:%S%z')}", encoding="utf-8")
        yield
    finally:
        shutil.rmtree(lock, ignore_errors=True)


def read_build() -> dict | None:
    try:
        return json.loads((OUT / "build.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def status() -> dict:
    b = read_build()
    built = b is not None and PROGRAM.exists()
    return {"built": built, "toolchain": b.get("toolchain") if b else None, "built_at": b.get("built_at") if b else None,
            "current": built and b.get("sources") == digest(), "sources": b.get("sources") if b else None}


def build(force: bool = False, log=print) -> bool:
    if not force and status()["current"]:
        log("The Body Atlas's engine is current.")
        return True
    kind = toolchain()
    if kind is None:
        log("No C++ compiler: install g++ 11 or later (on Windows, in WSL).")
        return False
    OUT.mkdir(parents=True, exist_ok=True)
    stamp = digest()
    args = ["g++", "-std=c++20", "-O2", "-DNDEBUG", "-pthread"]
    args += [f"-I{path_for(kind, d)}" for d in include_dirs()]
    args += [path_for(kind, SOURCE), "-o", path_for(kind, PROGRAM)]
    log(f"Compiling the Body Atlas's engine with g++ ({kind}); about a minute.")
    t0 = time.time()
    with wsl_lock(kind, log):
        r = subprocess.run(command(kind, args), capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        log(r.stderr[-6000:])
        log("The engine did not compile.")
        return False
    record = {"sources": stamp, "toolchain": kind, "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "seconds": round(time.time() - t0, 1)}
    (OUT / "build.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    log(f"Built the engine in {record['seconds']} s: {PROGRAM}")
    return True


class Worker:
    """One engine process (--serve), asked one request at a time; restarted if it ends."""

    def __init__(self):
        self.lock = threading.Lock()
        self.process: subprocess.Popen | None = None

    def _start(self) -> None:
        b = read_build() or {}
        kind = b.get("toolchain") or toolchain() or "native"
        args = [path_for(kind, PROGRAM), "--serve", "--data", path_for(kind, PRODUCTS)]
        self.process = subprocess.Popen(command(kind, args), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=subprocess.DEVNULL, text=True, encoding="utf-8", bufsize=1)

    def ask(self, request: dict) -> dict:
        with self.lock:
            for _attempt in range(2):
                if self.process is None or self.process.poll() is not None:
                    self._start()
                try:
                    self.process.stdin.write(json.dumps(request) + "\n")
                    self.process.stdin.flush()
                    line = self.process.stdout.readline()
                    if line:
                        return json.loads(line)
                except (OSError, ValueError):
                    pass
                self.stop()
            raise RuntimeError("the Body Atlas's engine did not answer")

    def stop(self) -> None:
        if self.process is not None:
            with contextlib.suppress(Exception):
                self.process.kill()
            self.process = None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--force", action="store_true")
    sub.add_parser("status")
    a = sub.add_parser("ask")
    a.add_argument("request")
    args = parser.parse_args()
    if args.cmd == "build":
        return 0 if build(args.force) else 1
    if args.cmd == "status":
        print(json.dumps(status(), indent=1))
        return 0
    w = Worker()
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        print(json.dumps(w.ask(json.loads(args.request)), ensure_ascii=False))
    finally:
        w.stop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
