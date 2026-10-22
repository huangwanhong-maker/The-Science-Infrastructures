"""The compute bridge from Python to the domains' JavaScript: a pool of persistent
Node workers (api/service/compute.mjs), each answering one call at a time over
JSON lines. A call that takes longer than its limit has its worker killed and a
new one started in its place, so that a runaway computation costs one worker for
a few seconds and never the service.

The bridge stands in for the C++ bindings to come (ADR 0003): when a domain's
C++ exists, the same calls may be answered by it, and the API does not change.

    pool = Pool(root)                                  # Node from PATH, two workers, 10 s a call
    pool.call("astronomy", "sunPosition", [2460000.5]) # {"ra": ..., "dec": ..., ...}
    pool.units()                                       # the units, their functions and constants
"""

from __future__ import annotations

import collections
import itertools
import json
import os
import queue
import shutil
import subprocess
import threading
from pathlib import Path

from .local import BadRequest, NotFound, Unavailable


class Timeout(RuntimeError):
    """A call that ran past its limit (the service answers 504)."""


class Failed(RuntimeError):
    """A function that raised an error, or a result that is not JSON (the service answers 422)."""


class Worker:
    """One Node process running compute.mjs, with a thread reading its answers."""

    def __init__(self, node: str, script: Path, root: Path):
        self.process = subprocess.Popen(
            [node, str(script)], cwd=str(root), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, encoding="utf-8", bufsize=1,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        self.answers: queue.Queue = queue.Queue()
        self.errors = collections.deque(maxlen=40)
        threading.Thread(target=self._read, args=(self.process.stdout, self.answers.put, True), daemon=True).start()
        threading.Thread(target=self._read, args=(self.process.stderr, self.errors.append, False), daemon=True).start()

    @staticmethod
    def _read(stream, sink, closing: bool) -> None:
        """Passes each line to `sink`, and None at the end when `closing`."""
        try:
            for line in stream:
                sink(line)
        except (OSError, ValueError):
            pass
        if closing:
            sink(None)

    def alive(self) -> bool:
        return self.process.poll() is None

    def ask(self, request: dict, seconds: float) -> dict:
        try:
            self.process.stdin.write(json.dumps(request, allow_nan=False) + "\n")
            self.process.stdin.flush()
        except ValueError:
            raise BadRequest("the arguments are not JSON (NaN and infinities are not; send null)") from None
        except OSError:
            raise Unavailable(f"the compute worker has stopped: {''.join(self.errors)[-600:]}") from None
        while True:
            try:
                line = self.answers.get(timeout=seconds)
            except queue.Empty:
                self.kill()
                raise Timeout(f"the call took longer than {seconds:g} s and was stopped") from None
            if line is None:
                raise Unavailable(f"the compute worker stopped: {''.join(self.errors)[-600:].strip()}")
            answer = json.loads(line)
            if answer.get("id") == request["id"]:
                return answer

    def kill(self) -> None:
        try:
            self.process.kill()
        except OSError:
            pass


class Pool:
    """Persistent Node workers for calls into the domains and foundations."""

    def __init__(self, root: Path, *, size: int = 2, seconds: float = 10.0, node: str | None = None):
        self.root = Path(root)
        self.script = self.root / "api" / "service" / "compute.mjs"
        self.node = node or os.environ.get("SCIENCE_NODE") or shutil.which("node")
        self.seconds = seconds
        self.size = max(1, size)
        self.idle: queue.LifoQueue = queue.LifoQueue()
        self.started = 0
        self.lock = threading.Lock()
        self.ids = itertools.count(1)
        self._units = None

    def available(self) -> str | None:
        """None when calls can be made; else why not."""
        if not self.node:
            return "calling the domains' functions needs Node (18 or later) on this computer"
        if not self.script.is_file():
            return f"the compute worker is missing: {self.script}"
        return None

    def _take(self) -> Worker:
        reason = self.available()
        if reason:
            raise Unavailable(reason)
        with self.lock:
            if self.idle.empty() and self.started < self.size:
                self.started += 1
                return Worker(self.node, self.script, self.root)
        try:
            worker = self.idle.get(timeout=max(60.0, 3 * self.seconds))
        except queue.Empty:
            raise Unavailable("every compute worker is busy: try again in a while") from None
        return worker if worker.alive() else self._replace()

    def _replace(self) -> Worker:
        return Worker(self.node, self.script, self.root)

    def _ask(self, request: dict, seconds: float | None = None) -> dict:
        worker = self._take()
        request = {"id": next(self.ids), **request}
        try:
            answer = worker.ask(request, seconds or self.seconds)
        except (Timeout, Unavailable):
            worker.kill()
            self.idle.put(self._replace())
            raise
        except BaseException:
            self.idle.put(worker)
            raise
        self.idle.put(worker)
        if answer.get("ok"):
            return answer["result"]
        if answer.get("kind") == "not-found":
            raise NotFound(answer.get("error", "not found"))
        raise Failed(answer.get("error", "the call failed"))

    def units(self) -> list[dict]:
        """Every unit with its functions (name, parameters, the comment above it), constants and classes."""
        if self._units is None:
            self._units = self._ask({"op": "units"}, max(self.seconds, 30))
        return self._units

    def unit(self, name: str) -> dict:
        found = next((u for u in self.units() if u["unit"] == name), None)
        if not found:
            raise NotFound(f"no unit {name!r}: the units are {', '.join(u['unit'] for u in self.units())}")
        return found

    def call(self, unit: str, function: str, args: list, seconds: float | None = None):
        """A function's result, as JSON."""
        if not isinstance(args, list):
            raise BadRequest("args is a list")
        return self._ask({"op": "call", "unit": unit, "function": function, "args": args}, seconds)

    def close(self) -> None:
        while not self.idle.empty():
            worker = self.idle.get()
            try:
                worker.process.stdin.close()
            except OSError:
                pass
            worker.kill()
