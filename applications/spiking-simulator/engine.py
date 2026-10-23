"""The engine of the Spiking Simulator: the neuroscience domain's simulator of spiking networks.

The engine is the domain's command-line tool (``domains/neuroscience/snn/tools/snn.cpp``),
header-only C++20 on calx and physis. It is compiled here by g++ 11 or later: natively
where there is one, and otherwise in WSL (Windows without a compiler), into
``data/cache/engines/snn/``. Beside the program it keeps ``build.json``, with the hash of
every source it was built from (a change to any header makes it stale), and
``models.json``, the catalogue of models that the page reads.

A run is the tool reading a specification (JSON, quantities with their units) on its
standard input and writing its results (JSON, SI units) on its standard output; the
library it runs is pure.

Usage::

    python applications/spiking-simulator/engine.py build [--force]
    python applications/spiking-simulator/engine.py status
    python applications/spiking-simulator/engine.py run <spec.json> [--threads N]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
OUT = ROOT / "data" / "cache" / "engines" / "snn"
PROGRAM = OUT / "snn"
SOURCE = ROOT / "domains" / "neuroscience" / "snn" / "tools" / "snn.cpp"


def include_dirs() -> list[Path]:
    """The include directories of the engine: every module of calx and physis, and snn."""
    dirs = []
    for pattern in ("foundations/calx/numerics/*/include", "foundations/calx/compute/*/include", "foundations/physis/*/include"):
        dirs += sorted(p for p in ROOT.glob(pattern) if p.is_dir())
    dirs.append(ROOT / "domains" / "neuroscience" / "snn" / "include")
    return dirs


def digest() -> str:
    """The hash of the engine's sources: every header it may include, and the tool itself."""
    h = hashlib.sha256()
    files = sorted({f for d in include_dirs() for f in d.rglob("*.hpp")} | {SOURCE})
    for f in files:
        h.update(f.relative_to(ROOT).as_posix().encode())
        h.update(b"\0")
        h.update(f.read_bytes())
    return h.hexdigest()


def wsl_path(p: Path) -> str:
    """C:\\Users\\x -> /mnt/c/Users/x."""
    s = str(p.resolve())
    if len(s) > 1 and s[1] == ":":
        return f"/mnt/{s[0].lower()}{s[2:]}".replace("\\", "/")
    return s.replace("\\", "/")


def toolchain() -> str | None:
    """'native' where g++ is on the PATH, 'wsl' where it is in WSL, else None."""
    if shutil.which("g++") and os.name != "nt":
        return "native"
    if os.name == "nt" and shutil.which("wsl"):
        try:
            r = subprocess.run(["wsl", "--exec", "sh", "-c", "command -v g++"], capture_output=True, timeout=60)
            if r.returncode == 0 and r.stdout.strip():
                return "wsl"
        except (OSError, subprocess.TimeoutExpired):
            pass
    if shutil.which("g++"):
        return "native"
    return None


def command(kind: str, args: list[str]) -> list[str]:
    """A command for the toolchain: paths in args are already in its form."""
    return ["wsl", "--exec", *args] if kind == "wsl" else args


def path_for(kind: str, p: Path) -> str:
    return wsl_path(p) if kind == "wsl" else str(p)


def read_build() -> dict | None:
    try:
        return json.loads((OUT / "build.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def status() -> dict:
    """Whether the engine is built, and from the present sources."""
    b = read_build()
    built = b is not None and PROGRAM.exists()
    out = {"built": built, "toolchain": b.get("toolchain") if b else None, "built_at": b.get("built_at") if b else None}
    out["current"] = built and b.get("sources") == digest()
    return out


def build(force: bool = False, log=print) -> bool:
    """Compiles the engine unless it is current; writes the catalogue of models."""
    if not force and status()["current"]:
        log("The engine is current.")
        return True
    kind = toolchain()
    if kind is None:
        log("No C++ compiler: install g++ 11 or later (on Windows, in WSL).")
        return False
    OUT.mkdir(parents=True, exist_ok=True)
    sources = digest()
    args = ["g++", "-std=c++20", "-O3", "-march=native", "-pthread", "-DNDEBUG"]
    args += [f"-I{path_for(kind, d)}" for d in include_dirs()]
    args += [path_for(kind, SOURCE), "-o", path_for(kind, PROGRAM)]
    log(f"Compiling the engine with g++ ({kind}); this takes a minute or two.")
    t0 = time.time()
    r = subprocess.run(command(kind, args), capture_output=True, text=True)
    if r.returncode != 0:
        log(r.stderr[-4000:])
        log("The engine did not compile.")
        return False
    m = subprocess.run(command(kind, [path_for(kind, PROGRAM), "models"]), capture_output=True, text=True)
    if m.returncode != 0:
        log(m.stderr)
        return False
    (OUT / "models.json").write_text(m.stdout, encoding="utf-8")
    v = subprocess.run(command(kind, [path_for(kind, PROGRAM), "version"]), capture_output=True, text=True)
    record = {"sources": sources, "toolchain": kind, "version": v.stdout.strip(), "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "seconds": round(time.time() - t0, 1), "flags": args[1:6]}
    (OUT / "build.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    log(f"Built in {record['seconds']} s: {PROGRAM}")
    return True


def models() -> list | None:
    try:
        return json.loads((OUT / "models.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


class Run:
    """One run of the engine, which another thread may cancel."""

    def __init__(self, spec: str, threads: int | None = None):
        b = read_build() or {}
        self.kind = b.get("toolchain") or toolchain() or "native"
        args = [path_for(self.kind, PROGRAM), "run", "-"]
        if threads:
            args += ["--threads", str(int(threads))]
        self.args = command(self.kind, args)
        self.spec = spec
        self.process: subprocess.Popen | None = None
        self.cancelled = False
        self.lock = threading.Lock()

    def execute(self, timeout: float) -> tuple[int, str, str]:
        """(return code, results, messages); the code is -1 on a timeout or a cancellation."""
        with self.lock:
            if self.cancelled:
                return -1, "", "cancelled"
            self.process = subprocess.Popen(self.args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                            text=True, encoding="utf-8")
        try:
            out, err = self.process.communicate(self.spec, timeout=timeout)
        except subprocess.TimeoutExpired:
            self.cancel()
            out, err = self.process.communicate()
            return -1, "", f"the run took longer than {timeout:g} s and was stopped"
        if self.cancelled:
            return -1, "", "cancelled"
        return self.process.returncode, out, err

    def cancel(self) -> None:
        with self.lock:
            self.cancelled = True
            if self.process and self.process.poll() is None:
                self.process.kill()
                if self.kind == "wsl":
                    # Killing wsl.exe does not always end the program inside it.
                    subprocess.run(["wsl", "--exec", "pkill", "-f", wsl_path(PROGRAM) + " run"], capture_output=True)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Build or run the Spiking Simulator's engine.")
    sub = parser.add_subparsers(dest="what", required=True)
    b = sub.add_parser("build", help="compile the engine unless it is current")
    b.add_argument("--force", action="store_true", help="compile it again")
    sub.add_parser("status", help="whether the engine is built and current")
    r = sub.add_parser("run", help="run a specification and print its results")
    r.add_argument("spec")
    r.add_argument("--threads", type=int)
    args = parser.parse_args(argv)
    if args.what == "build":
        return 0 if build(args.force) else 1
    if args.what == "status":
        print(json.dumps(status(), indent=2))
        return 0
    if not status()["current"] and not build():
        return 1
    code, out, err = Run(Path(args.spec).read_text(encoding="utf-8"), args.threads).execute(timeout=24 * 3600)
    sys.stderr.write(err)
    sys.stdout.write(out)
    return code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
