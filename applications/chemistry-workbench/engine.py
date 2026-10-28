"""The engines of the Chemistry Workbench's plug-ins, built and kept running.

Each plug-in with C++ has ``plugins/<name>/engine/main.cpp``, header-only C++20 on the foundations
(morphe, calx, physis) and the API's JSON (``api/cpp/science/json.hpp``). It is compiled here by
g++ 11 or later, natively where there is one and otherwise in WSL, into
``data/cache/engines/chemistry/<name>/``, beside a ``build.json`` holding the hash of every file it
was built from: its sources and the headers they include, followed transitively, so that work
elsewhere in the foundations does not make it stale. Under WSL the build takes the lock that the
sessions working on this repository share (``%TEMP%/science-wsl.lock``, a directory created
atomically, with an ``owner`` file), so that two compilations do not run there at once.

A built engine speaks JSON lines (``--serve``): a request a line, an answer a line. ``Worker``
keeps one process for each plug-in and passes it one request at a time.

Usage::

    python applications/chemistry-workbench/engine.py build [<plugin>] [--force]
    python applications/chemistry-workbench/engine.py status
    python applications/chemistry-workbench/engine.py ask <plugin> '<request json>'
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
PLUGINS = APP / "plugins"
OUT = ROOT / "data" / "cache" / "engines" / "chemistry"
INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]', re.M)


def plugins() -> list[dict]:
    """Every plug-in's manifest (plugins/<name>/plugin.json), with its name and directory."""
    out = []
    for manifest in sorted(PLUGINS.glob("*/plugin.json")):
        try:
            p = json.loads(manifest.read_text(encoding="utf-8"))
        except ValueError:
            continue
        p["name"] = manifest.parent.name
        p["engine"] = (manifest.parent / "engine" / "main.cpp").is_file()
        out.append(p)
    return out


def source(name: str) -> Path:
    return PLUGINS / name / "engine" / "main.cpp"


def program(name: str) -> Path:
    return OUT / name / name


def include_dirs(name: str) -> list[Path]:
    dirs = [ROOT / "foundations" / "morphe" / "include"]
    for pattern in ("foundations/morphe/*/include", "foundations/calx/numerics/*/include",
                    "foundations/calx/compute/*/include", "foundations/physis/*/include"):
        dirs += sorted(p for p in ROOT.glob(pattern) if p.is_dir())
    dirs += [ROOT / "api" / "cpp", PLUGINS / name / "engine"]
    return [d for d in dirs if d.is_dir()]


def sources(name: str) -> list[Path]:
    """The engine's source and the project headers it includes, transitively."""
    dirs = include_dirs(name)
    seen: set[Path] = set()
    todo = [source(name)] if source(name).exists() else []
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


def digest(name: str) -> str:
    h = hashlib.sha256()
    for f in sources(name):
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
        (lock / "owner").write_text(f"chemistry-workbench build {time.strftime('%Y-%m-%dT%H:%M:%S%z')}", encoding="utf-8")
        yield
    finally:
        shutil.rmtree(lock, ignore_errors=True)


def read_build(name: str) -> dict | None:
    try:
        return json.loads((OUT / name / "build.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def status(name: str) -> dict:
    b = read_build(name)
    built = b is not None and program(name).exists()
    return {"plugin": name, "built": built, "toolchain": b.get("toolchain") if b else None,
            "built_at": b.get("built_at") if b else None, "current": built and b.get("sources") == digest(name)}


def build(name: str, force: bool = False, log=print) -> bool:
    if not source(name).exists():
        log(f"The plug-in {name} has no engine.")
        return True
    if not force and status(name)["current"]:
        log(f"The engine of {name} is current.")
        return True
    kind = toolchain()
    if kind is None:
        log("No C++ compiler: install g++ 11 or later (on Windows, in WSL).")
        return False
    (OUT / name).mkdir(parents=True, exist_ok=True)
    stamp = digest(name)
    args = ["g++", "-std=c++20", "-O2", "-DNDEBUG"]
    args += [f"-I{path_for(kind, d)}" for d in include_dirs(name)]
    args += [path_for(kind, source(name)), "-o", path_for(kind, program(name))]
    log(f"Compiling the engine of {name} with g++ ({kind}); a minute or two.")
    t0 = time.time()
    with wsl_lock(kind, log):
        r = subprocess.run(command(kind, args), capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        log(r.stderr[-6000:])
        log("The engine did not compile.")
        return False
    record = {"sources": stamp, "toolchain": kind, "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "seconds": round(time.time() - t0, 1)}
    (OUT / name / "build.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    log(f"Built {name} in {record['seconds']} s: {program(name)}")
    return True


class Worker:
    """One engine process (--serve) for a plug-in, asked one request at a time; restarted if it ends."""

    def __init__(self, name: str):
        self.name = name
        self.lock = threading.Lock()
        self.process: subprocess.Popen | None = None

    def _start(self) -> None:
        b = read_build(self.name) or {}
        kind = b.get("toolchain") or toolchain() or "native"
        self.process = subprocess.Popen(command(kind, [path_for(kind, program(self.name)), "--serve"]), stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, encoding="utf-8", bufsize=1)

    def ask(self, request: dict) -> dict:
        with self.lock:
            for attempt in range(2):
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
            raise RuntimeError(f"the engine of {self.name} did not answer")

    def stop(self) -> None:
        if self.process is not None:
            with contextlib.suppress(Exception):
                self.process.kill()
            self.process = None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("plugin", nargs="?")
    b.add_argument("--force", action="store_true")
    sub.add_parser("status")
    a = sub.add_parser("ask")
    a.add_argument("plugin")
    a.add_argument("request")
    args = parser.parse_args()
    if args.cmd == "build":
        names = [args.plugin] if args.plugin else [p["name"] for p in plugins() if p["engine"]]
        ok = all(build(n, args.force) for n in names)
        return 0 if ok else 1
    if args.cmd == "status":
        print(json.dumps([status(p["name"]) for p in plugins() if p["engine"]], indent=1))
        return 0
    w = Worker(args.plugin)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        print(json.dumps(w.ask(json.loads(args.request)), ensure_ascii=False))
    finally:
        w.stop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
