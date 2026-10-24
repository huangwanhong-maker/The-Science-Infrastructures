"""The engine of the Genome Atlas: the genomics and evolution domains' C++, as one program.

The engine is the units of ``engine/`` (``genome.cpp``, the loop that answers, and one
``ops_*.cpp`` for each group of operations), header-only C++20 on the genomics domain
(``domains/genomics/``), the evolution domain (``domains/evolution/``), morphe and calx, with the
API's JSON (``api/cpp/science/json.hpp``). It is compiled here by g++ 11 or later, natively where
there is one and otherwise in WSL, into ``data/cache/engines/genome/``, with -O3 and
-march=native, so that calx's lanes use the widest vectors of this machine. Each unit is compiled
on its own into ``obj/``, one after another, and only when a file it includes has changed; then
the objects are linked. One unit holding every module of both domains needed more memory than a
small WSL machine gives without failing. Beside the program it keeps ``build.json``: the hash of
every source it was built from (a change to any header makes it stale), each unit's, the flags,
and the target the compiler was allowed.

The program speaks JSON lines (engine/protocol.hpp): with --serve it is a worker that answers
requests until its input ends, which is how the server and its cluster keep it (cluster.py).

Usage::

    python applications/genome-atlas/engine.py build [--force]
    python applications/genome-atlas/engine.py status
    python applications/genome-atlas/engine.py run <request.json | ->
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
OUT = ROOT / "data" / "cache" / "engines" / "genome"
PROGRAM = OUT / "genome"
OBJECTS = OUT / "obj"
SOURCE = APP / "engine" / "genome.cpp"


def units() -> list[Path]:
    """The translation units: the loop (genome.cpp) and one for each group of operations."""
    return sorted((APP / "engine").glob("*.cpp"))


def include_dirs() -> list[Path]:
    """The include directories: the two domains' modules, morphe's and calx's, the API's C++, the engine's own."""
    dirs: list[Path] = []
    for pattern in ("domains/genomics/*/include", "domains/evolution/*/include", "foundations/morphe/*/include",
                    "foundations/calx/numerics/*/include", "foundations/calx/compute/*/include",
                    "foundations/calx/compute/backends/*/include"):
        dirs += sorted(p for p in ROOT.glob(pattern) if p.is_dir())
    dirs += [ROOT / "foundations" / "morphe" / "include", ROOT / "api" / "cpp", APP / "engine"]
    return [d for d in dirs if d.is_dir()]


def sources(listed: list[str] | None = None) -> list[Path]:
    """The files the engine is built from: those the compiler said it read (build.json), or, before
    a build, every header it may include."""
    files = {ROOT / f for f in listed} if listed else {f for d in include_dirs() for f in d.rglob("*.hpp")}
    files |= {f for f in (APP / "engine").rglob("*") if f.suffix in (".cpp", ".hpp")}
    return sorted(files)


def digest(listed: list[str] | None = None) -> str:
    h = hashlib.sha256()
    for f in sources(listed):
        h.update(f.relative_to(ROOT).as_posix().encode())
        h.update(b"\0")
        h.update(f.read_bytes() if f.exists() else b"(missing)")
    return h.hexdigest()


def run_tool(kind: str, args: list[str], **kw) -> subprocess.CompletedProcess:
    """A command of the toolchain. Under WSL, a launch the service refuses (no output from the
    tool at all, and a status of −1) is tried again a few times; the tool's own failures are not."""
    for attempt in range(4):
        r = subprocess.run(command(kind, args), capture_output=True, text=True, encoding="utf-8", errors="replace", **kw)
        refused = kind == "wsl" and r.returncode != 0 and not r.stderr.strip() and r.returncode in (-1, 4294967295, 1)
        if not refused:
            return r
        time.sleep(5 * (attempt + 1))
    return r


def quote(a: str) -> str:
    return "'" + a.replace("'", "'\\''") + "'"


def run_script(kind: str, lines: list[str]) -> subprocess.CompletedProcess:
    """A bash script of the toolchain's commands, run by one launch (under WSL, a launch the
    service may refuse; run_tool tries it again)."""
    OUT.mkdir(parents=True, exist_ok=True)
    script = OUT / "step.sh"
    script.write_text("set -u\n" + "\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return run_tool(kind, ["bash", path_for(kind, script)])


def parse_make(text: str, kind: str) -> list[str] | None:
    """The repository's files of a rule in make's syntax ("x.o: a.cpp b.hpp \\<newline> c.hpp", a
    space within a path escaped by a backslash)."""
    if ":" not in text:
        return None
    text = text.split(":", 1)[1].replace("\\\n", " ")
    paths, cur, i = [], "", 0
    while i < len(text):
        if text[i] == "\\" and i + 1 < len(text) and text[i + 1] == " ":
            cur += " "
            i += 2
            continue
        if text[i].isspace():
            if cur:
                paths.append(cur)
            cur = ""
        else:
            cur += text[i]
        i += 1
    if cur:
        paths.append(cur)
    root = path_for(kind, ROOT).rstrip("/") + "/"
    out = sorted({p[len(root):] for p in paths if p.startswith(root)})
    return out or None


def dependencies_of(kind: str, flags: list[str], todo: list[Path]) -> dict[str, list[str] | None]:
    """Each unit's files, as the compiler lists them (g++ -MM), in one launch."""
    inc = " ".join(quote(f"-I{path_for(kind, d)}") for d in include_dirs())
    lines = []
    for unit in todo:
        lines.append(f"echo '@@@ {unit.name}'")
        lines.append(f"g++ {' '.join(flags)} {inc} -MM {quote(path_for(kind, unit))} || echo '@@@ failed'")
    r = run_script(kind, lines)
    out: dict[str, list[str] | None] = {}
    last = None
    for block in r.stdout.split("@@@ ")[1:]:
        name, _, body = block.partition("\n")
        name = name.strip()
        if name == "failed":  # the unit before it could not be read
            if last:
                out[last] = None
            continue
        out[name] = parse_make(body, kind)
        last = name
    return out


def wsl_path(p: Path) -> str:
    """C:\\Users\\x -> /mnt/c/Users/x."""
    s = str(Path(p).resolve())
    if len(s) > 1 and s[1] == ":":
        return f"/mnt/{s[0].lower()}{s[2:]}".replace("\\", "/")
    return s.replace("\\", "/")


def toolchain() -> str | None:
    """'native' where g++ is on the PATH, 'wsl' where it is in WSL, else None."""
    if shutil.which("g++") and os.name != "nt":
        return "native"
    if os.name == "nt" and shutil.which("wsl"):
        # WSL's service can refuse a launch while many are starting at once: try a few times.
        for attempt in range(4):
            try:
                r = subprocess.run(["wsl", "--exec", "sh", "-c", "command -v g++"], capture_output=True, timeout=60)
                if r.returncode == 0 and r.stdout.strip():
                    return "wsl"
            except (OSError, subprocess.TimeoutExpired):
                pass
            time.sleep(3 * (attempt + 1))
    if shutil.which("g++"):
        return "native"
    return None


def command(kind: str, args: list[str]) -> list[str]:
    return ["wsl", "--exec", *args] if kind == "wsl" else args


def path_for(kind: str, p: Path) -> str:
    return wsl_path(p) if kind == "wsl" else str(p)


def read_build() -> dict | None:
    try:
        return json.loads((OUT / "build.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def status() -> dict:
    b = read_build()
    built = b is not None and PROGRAM.exists()
    out = {"built": built, "toolchain": b.get("toolchain") if b else None, "built_at": b.get("built_at") if b else None,
           "target": b.get("target") if b else None, "source": SOURCE.exists()}
    out["current"] = built and SOURCE.exists() and b.get("sources") == digest(b.get("files"))
    return out


def worker_command(kind: str | None = None) -> list[str]:
    """The command that starts a local worker."""
    kind = kind or (read_build() or {}).get("toolchain") or toolchain() or "native"
    return command(kind, [path_for(kind, PROGRAM), "--serve"])


def ask(request: dict, timeout: float = 120) -> dict:
    """One request to a fresh engine, its answer (for builds and the command line)."""
    kind = (read_build() or {}).get("toolchain") or toolchain() or "native"
    r = subprocess.run(command(kind, [path_for(kind, PROGRAM)]), input=json.dumps(request).encode("utf-8"),
                       capture_output=True, timeout=timeout)
    out = r.stdout.decode("utf-8", errors="replace")
    if not out.strip():
        raise RuntimeError(r.stderr.decode("utf-8", errors="replace").strip() or f"the engine ended with status {r.returncode}")
    return json.loads(out)


def build(force: bool = False, log=print) -> bool:
    """Compiles the engine's units that changed, links them, and asks the engine for its target,
    as a proof; nothing when it is current."""
    if not SOURCE.exists():
        log(f"The engine's source, {SOURCE.relative_to(ROOT).as_posix()}, is not there.")
        return False
    if not force and status()["current"]:
        log("The engine is current.")
        return True
    kind = toolchain()
    if kind is None:
        log("No C++ compiler: install g++ 11 or later (on Windows, in WSL).")
        return False
    OBJECTS.mkdir(parents=True, exist_ok=True)
    flags = ["-std=c++20", "-O3", "-march=native", "-pthread", "-DNDEBUG"]
    previous = (read_build() or {}).get("units", {}) if not force else {}
    t0 = time.time()
    todo = units()
    deps = dependencies_of(kind, flags, todo)
    if not deps:
        log("The toolchain did not answer (under WSL, its service refused the launch); try again in a minute.")
        return False
    all_files: set[str] = set()
    stamps: dict[str, str] = {}
    objects: list[Path] = []
    changed: list[tuple[Path, Path]] = []
    for unit in todo:
        files = deps.get(unit.name)
        if files is None:
            log(f"The compiler could not read {unit.name}'s includes.")
            files = [unit.relative_to(ROOT).as_posix()]
        all_files |= set(files)
        stamp = digest(files)
        stamps[unit.name] = stamp
        obj = OBJECTS / (unit.stem + ".o")
        objects.append(obj)
        if not (previous.get(unit.name) == stamp and obj.exists()):
            changed.append((unit, obj))
    # The changed units one after another, then the link, in one launch.
    inc = " ".join(quote(f"-I{path_for(kind, d)}") for d in include_dirs())
    lines = []
    for unit, obj in changed:
        lines.append(f"echo '@@@ {unit.name}'; t=$(date +%s)")
        lines.append(f"g++ {' '.join(flags)} {inc} -c {quote(path_for(kind, unit))} -o {quote(path_for(kind, obj))} 2>&1 || {{ echo '@@@ failed'; exit 1; }}")
        lines.append("echo \"@@@ seconds $(( $(date +%s) - t ))\"")
    lines.append("echo '@@@ link'")
    lines.append(f"g++ {' '.join(flags)} {' '.join(quote(path_for(kind, o)) for o in objects)} -o {quote(path_for(kind, PROGRAM))} 2>&1 || {{ echo '@@@ failed'; exit 1; }}")
    log(f"Compiling {len(changed)} of the engine's {len(todo)} units with g++ ({kind}), each on its own" + ("; a few minutes the first time." if changed else "."))
    r = run_script(kind, lines)
    if r.returncode != 0 or "@@@ failed" in r.stdout:
        log(r.stdout[-6000:] or r.stderr[-2000:] or f"the build ended with status {r.returncode}")
        log("The engine did not compile.")
        return False
    files = sorted(all_files)
    stamp = digest(files)
    try:
        answer = ask({"id": 0, "op": "target"})
    except (RuntimeError, ValueError, subprocess.TimeoutExpired) as error:
        log(f"The engine compiled but did not answer: {error}")
        return False
    if not answer.get("ok"):
        log(f"The engine compiled but did not answer: {answer.get('error')}")
        return False
    target = answer["result"]
    record = {"sources": stamp, "files": files, "units": stamps, "toolchain": kind, "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "seconds": round(time.time() - t0, 1), "flags": flags,
              "target": {k: target.get(k) for k in ("brand", "threads", "built_for", "lane_bytes", "best")},
              "operations": target.get("operations", [])}
    (OUT / "build.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    log(f"Built in {record['seconds']} s for {target.get('built_for')} ({target.get('brand')}): {PROGRAM}")
    return True


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Build or run the Genome Atlas's engine.")
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
    text = sys.stdin.buffer.read().decode("utf-8") if args.request == "-" else Path(args.request).read_text(encoding="utf-8")
    print(json.dumps(ask(json.loads(text), timeout=24 * 3600)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
