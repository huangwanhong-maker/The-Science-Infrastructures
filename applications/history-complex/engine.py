"""The engine of the History Complex: the history domain's construction and readings, as one program.

The engine is ``engine/history.cpp``, header-only C++20 on the history domain
(``domains/history/``: the record and its chronology, the construction of the complex, the
readings of its signed ties, its spin foam and field theory), morphe (homology, the fundamental
group, cycles, canonical forms), physis (the ℤ₂ gauge theory of the foam, exact and sampled;
recoupling) and calx, with the API's JSON (``api/cpp/science/json.hpp``). It is compiled here by g++ 11 or later: natively where
there is one, and otherwise in WSL (Windows without a compiler), into
``data/cache/engines/history/``. Beside the program it keeps ``build.json``, with the hash of
every source it was built from (a change to any header makes it stale).

A run is the program reading one request (JSON) on its standard input and writing one answer
(JSON) on its standard output: a canon built into its complex, with its invariants and readings,
or the network of the moment at a date (the protocol is in this application's README, §3). The
canon is named by its path; the library it runs reads nothing else.

Usage::

    python applications/history-complex/engine.py build [--force]
    python applications/history-complex/engine.py status
    python applications/history-complex/engine.py run <request.json | ->
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
OUT = ROOT / "data" / "cache" / "engines" / "history"
PROGRAM = OUT / "history"
SOURCE = APP / "engine" / "history.cpp"
CANONS = ROOT / "domains" / "history" / "canons"


def include_dirs() -> list[Path]:
    """The include directories of the engine: the history domain, morphe, calx, physis, the API's C++ and the engine's own."""
    dirs = sorted(p for p in ROOT.glob("domains/history/*/include") if p.is_dir())
    dirs += [ROOT / "foundations" / "morphe" / "include"]
    for pattern in ("foundations/morphe/*/include", "foundations/calx/numerics/*/include", "foundations/calx/compute/*/include",
                    "foundations/physis/*/include"):
        dirs += sorted(p for p in ROOT.glob(pattern) if p.is_dir())
    dirs += [ROOT / "api" / "cpp", APP / "engine"]
    return [d for d in dirs if d.is_dir()]


def sources() -> list[Path]:
    """Every file the engine may be built from: the headers it may include, and its own sources."""
    files = {f for d in include_dirs() for f in d.rglob("*.hpp")}
    if (APP / "engine").is_dir():
        files |= {f for f in (APP / "engine").rglob("*") if f.suffix in (".cpp", ".hpp")}
    return sorted(files)


def digest() -> str:
    h = hashlib.sha256()
    for f in sources():
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
    out = {"built": built, "toolchain": b.get("toolchain") if b else None, "built_at": b.get("built_at") if b else None,
           "source": SOURCE.exists()}
    out["current"] = built and SOURCE.exists() and b.get("sources") == digest()
    return out


def ask(kind: str, request: dict, timeout: float = 120) -> dict:
    """One request to the built program; its answer."""
    r = subprocess.run(command(kind, [path_for(kind, PROGRAM)]), input=json.dumps(request).encode("utf-8"),
                       capture_output=True, timeout=timeout)
    out = r.stdout.decode("utf-8", errors="replace")
    if r.returncode != 0 and not out.strip():
        raise RuntimeError(r.stderr.decode("utf-8", errors="replace").strip() or f"the engine ended with status {r.returncode}")
    return json.loads(out)


def canon_files() -> list[Path]:
    return sorted(CANONS.glob("*.jsonl")) if CANONS.is_dir() else []


def build(force: bool = False, log=print) -> bool:
    """Compiles the engine unless it is current; then builds the first canon with it, as a proof."""
    if not SOURCE.exists():
        log(f"The engine's source, {SOURCE.relative_to(ROOT).as_posix()}, is not there yet.")
        return False
    if not force and status()["current"]:
        log("The engine is current.")
        return True
    kind = toolchain()
    if kind is None:
        log("No C++ compiler: install g++ 11 or later (on Windows, in WSL).")
        return False
    OUT.mkdir(parents=True, exist_ok=True)
    stamp = digest()
    args = ["g++", "-std=c++20", "-O2", "-march=native", "-pthread", "-DNDEBUG"]
    args += [f"-I{path_for(kind, d)}" for d in include_dirs()]
    args += [path_for(kind, SOURCE), "-o", path_for(kind, PROGRAM)]
    log(f"Compiling the engine with g++ ({kind}); this takes a minute or two.")
    t0 = time.time()
    r = subprocess.run(command(kind, args), capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        log(r.stderr[-4000:])
        log("The engine did not compile.")
        return False
    canons = canon_files()
    built = 0
    for canon in canons[:1]:
        try:
            answer = ask(kind, {"op": "build", "canon": path_for(kind, canon)})
        except (RuntimeError, ValueError, subprocess.TimeoutExpired) as error:
            log(f"The engine compiled but did not build {canon.name}: {error}")
            return False
        if "complex" not in answer:
            log(f"The engine compiled but did not build {canon.name}: {answer.get('error', answer)}")
            return False
        built += 1
    record = {"sources": stamp, "toolchain": kind, "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "seconds": round(time.time() - t0, 1), "flags": args[1:6], "canons_tried": built}
    (OUT / "build.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    log(f"Built in {record['seconds']} s{', and the first canon built with it' if built else ''}: {PROGRAM}")
    return True


class Run:
    """One request to the engine, which another thread may cancel."""

    def __init__(self, request: str):
        b = read_build() or {}
        self.kind = b.get("toolchain") or toolchain() or "native"
        self.args = command(self.kind, [path_for(self.kind, PROGRAM)])
        self.request = request
        self.process: subprocess.Popen | None = None
        self.cancelled = False
        self.lock = threading.Lock()

    def execute(self, timeout: float) -> tuple[int, str, str]:
        """(return code, answer, messages); the code is -1 on a timeout or a cancellation."""
        with self.lock:
            if self.cancelled:
                return -1, "", "cancelled"
            self.process = subprocess.Popen(self.args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            out, err = self.process.communicate(self.request.encode("utf-8"), timeout=timeout)
        except subprocess.TimeoutExpired:
            self.cancel()
            self.process.communicate()
            return -1, "", f"the computation took longer than {timeout:g} s and was stopped"
        if self.cancelled:
            return -1, "", "cancelled"
        return self.process.returncode, out.decode("utf-8", errors="replace"), err.decode("utf-8", errors="replace")

    def cancel(self) -> None:
        with self.lock:
            self.cancelled = True
            if self.process and self.process.poll() is None:
                self.process.kill()
                if self.kind == "wsl":
                    # Killing wsl.exe does not always end the program inside it.
                    subprocess.run(["wsl", "--exec", "pkill", "-f", wsl_path(PROGRAM)], capture_output=True)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Build or run the History Complex's engine.")
    sub = parser.add_subparsers(dest="what", required=True)
    b = sub.add_parser("build", help="compile the engine unless it is current")
    b.add_argument("--force", action="store_true", help="compile it again")
    sub.add_parser("status", help="whether the engine is built and current")
    r = sub.add_parser("run", help="answer one request (a JSON file, or - for standard input)")
    r.add_argument("request")
    args = parser.parse_args(argv)
    if args.what == "build":
        return 0 if build(args.force) else 1
    if args.what == "status":
        print(json.dumps(status(), indent=2))
        return 0
    if not status()["current"] and not build():
        return 1
    # The engine speaks UTF-8 whatever the console's code page.
    text = sys.stdin.buffer.read().decode("utf-8") if args.request == "-" else Path(args.request).read_text(encoding="utf-8")
    code, out, err = Run(text).execute(timeout=24 * 3600)
    sys.stderr.write(err)
    sys.stdout.buffer.write(out.encode("utf-8"))
    sys.stdout.flush()
    return code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
