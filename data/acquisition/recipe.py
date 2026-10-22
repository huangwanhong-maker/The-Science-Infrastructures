"""A recipe: the data an application needs, fetched and built in order.

An application's README lists the commands that fetch its datasets and build
its products. A recipe is that list as data, kept beside the application
(``applications/<name>/recipe.json``), so that one list serves the command
line, the application's server and its page::

    {
      "title": "The Earth Atlas's data",
      "steps": [
        {"id": "natural-earth-50m", "title": "Natural Earth 1:50 million",
         "run": ["data/acquisition/fetch.py", "natural-earth-50m"],
         "done": "data/cache/natural-earth-50m/**/retrieval.json"},
        {"id": "geography", "title": "The map",
         "run": ["data/products/geography.py"],
         "done": "data/cache/products/geography-ne50m.json"},
        ...
      ]
    }

Each step runs one of the library's scripts from the repository root, with the
interpreter running this one. A step is done when each pattern in ``done`` (one
or a list) matches a file. Steps run in order, one at a time, because the
product builders all rewrite the index of products. A step that fails does not
stop the rest; the failures are reported at the end.

Usage::

    python data/acquisition/recipe.py applications/earth-atlas/recipe.json             # the steps not yet done
    python data/acquisition/recipe.py applications/earth-atlas/recipe.json --status    # which steps are done
    python data/acquisition/recipe.py applications/earth-atlas/recipe.json --all       # every step, again
    python data/acquisition/recipe.py applications/earth-atlas/recipe.json gfs oisst   # chosen steps, again

A server runs a recipe for its page through ``Job``, one run at a time, as the
atlases' servers do at ``/setup``.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]            # data/acquisition -> the repository
LOG_LINES = 400                   # a job keeps this many lines of output


def load(path: Path) -> dict:
    recipe = json.loads(Path(path).read_text(encoding="utf-8"))
    for step in recipe["steps"]:
        if isinstance(step.get("done"), str):
            step["done"] = [step["done"]]
    return recipe


def is_done(step: dict, root: Path = ROOT) -> bool:
    return all(any(root.glob(pattern)) for pattern in step.get("done", []))


def status(recipe: dict, root: Path = ROOT) -> list[dict]:
    """Each step, whether it is done, and when its newest output was written (Unix time)."""
    out = []
    for step in recipe["steps"]:
        files = [p for pattern in step.get("done", []) for p in root.glob(pattern)]
        newest = max((p.stat().st_mtime for p in files), default=None)
        out.append({"id": step["id"], "title": step["title"], "run": step["run"],
                    "done": is_done(step, root), "updated": newest})
    return out


def run_step(step: dict, root: Path = ROOT, log=print) -> bool:
    """Runs one step from the repository root, passing its output to ``log`` line by line."""
    log(f"$ python {' '.join(step['run'])}")
    process = subprocess.Popen(
        [sys.executable, *step["run"]], cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace", bufsize=1,
        env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
    for line in process.stdout:
        log(line.rstrip("\r\n").split("\r")[-1])      # a progress line keeps its last state
    return process.wait() == 0


def run(recipe: dict, chosen: list[str] | None = None, again: bool = False,
        root: Path = ROOT, log=print) -> list[str]:
    """Runs steps in order and returns the ids of those that failed. Without
    ``chosen``, runs the steps not yet done, or every step with ``again``."""
    failed = []
    for step in recipe["steps"]:
        if chosen is not None:
            if step["id"] not in chosen:
                continue
        elif is_done(step, root) and not again:
            continue
        log(f"\n### {step['title']}")
        if not run_step(step, root, log):
            failed.append(step["id"])
            log(f"### failed: {step['id']}")
    return failed


class Job:
    """A recipe run in a thread for a server, one run at a time. The server
    returns ``view()`` as JSON; the page shows the steps' state and the log."""

    def __init__(self, recipe_path: Path, root: Path = ROOT):
        self.path, self.root = Path(recipe_path), root
        self.lock = threading.RLock()
        self.thread: threading.Thread | None = None
        self.lines: list[str] = []
        self.state = {"state": "idle", "started": None, "finished": None, "failed": [], "chosen": None}

    def view(self) -> dict:
        recipe = load(self.path)
        with self.lock:
            return {"title": recipe["title"], "steps": status(recipe, self.root), **self.state,
                    "running": bool(self.thread and self.thread.is_alive()),
                    "log": "\n".join(self.lines[-60:])}

    def start(self, chosen: list[str] | None = None, again: bool = False) -> dict:
        with self.lock:
            if self.thread and self.thread.is_alive():
                raise RuntimeError("a run is already in progress")
            # The schedule's lock (schedule.py): a scheduled job and a run from the page
            # would both rewrite the index of products.
            import schedule
            self.schedule_lock = schedule.Lock()
            if not self.schedule_lock.acquire(wait=False):
                raise RuntimeError("a scheduled job is running; try again when it ends")
            self.lines = []
            self.state = {"state": "running", "started": time.time(), "finished": None, "failed": [],
                          "chosen": chosen}
            self.thread = threading.Thread(target=self._run, args=(chosen, again), daemon=True)
            self.thread.start()
        return self.view()

    def _log(self, line: str) -> None:
        with self.lock:
            self.lines.append(line)
            del self.lines[:-LOG_LINES]

    def _run(self, chosen: list[str] | None, again: bool) -> None:
        try:
            failed = run(load(self.path), chosen, again, self.root, self._log)
            with self.lock:
                self.state.update(state="failed" if failed else "done", failed=failed)
        except Exception as error:          # reported to the page, which shows it
            self._log(f"{type(error).__name__}: {error}")
            with self.lock:
                self.state.update(state="failed", failed=["recipe"])
        finally:
            with self.lock:
                self.state["finished"] = time.time()
            self.schedule_lock.release()


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("recipe", help="the recipe, such as applications/earth-atlas/recipe.json")
    parser.add_argument("steps", nargs="*", help="ids of steps to run again (default: the steps not yet done)")
    parser.add_argument("--status", action="store_true", help="show which steps are done, and run nothing")
    parser.add_argument("--all", action="store_true", help="run every step, done or not")
    args = parser.parse_args(argv)
    # Consoles that cannot show ′, ° or CJK; and lines as they come when piped to a log.
    sys.stdout.reconfigure(errors="replace", line_buffering=True)

    recipe = load(Path(args.recipe))
    ids = {step["id"] for step in recipe["steps"]}
    unknown = [s for s in args.steps if s not in ids]
    if unknown:
        print(f"unknown step(s): {', '.join(unknown)}; see --status", file=sys.stderr)
        return 2
    if args.status:
        print(recipe["title"])
        for s in status(recipe):
            when = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime(s["updated"])) if s["updated"] else ""
            print(f"  {'done' if s['done'] else '   —'}  {s['id']:36s} {when}")
        return 0
    failed = run(recipe, args.steps or None, args.all)
    if failed:
        print(f"\n{len(failed)} step(s) failed: {', '.join(failed)}")
        return 1
    print(f"\n{recipe['title']}: every step is done")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
