"""The Genome Workbench's cluster: engine workers kept alive, here and on other machines, and jobs
cut into shards among them.

A worker is the engine started with --serve (engine/protocol.hpp): one JSON request a line on
its input, one answer a line on its output, its genomes kept in memory between requests. The
cluster keeps several:

- **local workers**, engine processes on this computer (in WSL on Windows), each with every
  hardware thread in its own pool;
- **remote workers**, started by a command the user writes in
  ``data/cache/engines/genome/cluster.json`` (not part of the repository), such as
  ``["ssh", "eu_region_server", "~/The-Science-Infrastructures/data/cache/engines/genome/genome --serve"]``:
  the same pipes, carried by ssh. A remote engine must be built there from the same sources; the
  cluster compares the build each worker reports with this one's.

A request goes to the least busy worker. A **job** is an operation with three parts in the
engine: ``<op>.plan`` cuts it into shards (a shard is a few numbers: which reads, which
replicates, which block of a matrix; simulated inputs are regenerated from their seed on the
worker, so nothing large travels), ``<op>.shard`` computes one, and ``<op>.merge`` combines the
partial results in the order of the shards, so that the answer is the same whichever worker
computed which shard. Free workers pull shards from a queue (dynamic placement: a fast machine
takes more), and every shard's worker, start, end and the engine's own seconds are kept: the
trace the Run view draws.

Only simulated or inline inputs go to remote workers; a genome read from this computer's files
stays with the local ones.
"""

from __future__ import annotations

import json
import queue
import subprocess
import threading
import time
from pathlib import Path

import engine

CONFIG = engine.OUT / "cluster.json"


class WorkerError(RuntimeError):
    pass


class Worker:
    """One engine process and the requests waiting for its answers."""

    def __init__(self, name: str, command: list[str], remote: bool = False):
        self.name, self.command, self.remote = name, command, remote
        self.process: subprocess.Popen | None = None
        self.pending: dict[int, dict] = {}
        self.lock = threading.Lock()
        self.next_id = 1
        self.busy = 0
        self.served = 0
        self.seconds = 0.0
        self.target: dict | None = None
        self.error: str | None = None
        self.started_at: float | None = None

    def start(self) -> None:
        self.process = subprocess.Popen(self.command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.started_at = time.time()
        threading.Thread(target=self._read, daemon=True).start()
        threading.Thread(target=self._drain, daemon=True).start()
        try:
            self.target = self.request({"op": "target"}, timeout=60)
            self.error = None
        except WorkerError as error:
            self.error = str(error)
            raise

    def alive(self) -> bool:
        return self.process is not None and self.process.poll() is None

    def _drain(self) -> None:
        """Keeps the worker's messages, the last of which says why it ended."""
        proc = self.process
        for line in iter(proc.stderr.readline, b""):
            self.error = line.decode("utf-8", errors="replace").strip() or self.error

    def _read(self) -> None:
        proc = self.process
        for raw in iter(proc.stdout.readline, b""):
            try:
                answer = json.loads(raw.decode("utf-8"))
            except ValueError:
                continue
            slot = self.pending.get(answer.get("id"))
            if slot is not None:
                slot["answer"] = answer
                slot["event"].set()
        # The process ended: every waiting request fails.
        for slot in list(self.pending.values()):
            slot["event"].set()

    def request(self, request: dict, timeout: float = 120) -> dict:
        """Sends one request and waits for its answer; returns the result, or raises with the engine's reason."""
        if not self.alive():
            raise WorkerError(f"worker {self.name} is not running")
        with self.lock:
            rid = self.next_id
            self.next_id += 1
            slot = {"event": threading.Event(), "answer": None}
            self.pending[rid] = slot
            self.busy += 1
            line = json.dumps({**request, "id": rid}, ensure_ascii=False) + "\n"
            try:
                self.process.stdin.write(line.encode("utf-8"))
                self.process.stdin.flush()
            except OSError as error:
                self.pending.pop(rid, None)
                self.busy -= 1
                raise WorkerError(f"worker {self.name} cannot be written to: {error}") from error
        t0 = time.time()
        try:
            if not slot["event"].wait(timeout):
                self.stop()
                raise WorkerError(f"the computation took longer than {timeout:g} s and was stopped")
            answer = slot["answer"]
            if answer is None:
                raise WorkerError(f"worker {self.name} ended: {self.error or 'no reason given'}")
            self.served += 1
            self.seconds += time.time() - t0
            if not answer.get("ok"):
                raise WorkerError(answer.get("error") or "the engine failed")
            result = answer.get("result")
            if isinstance(result, dict):
                result.setdefault("_seconds", answer.get("seconds"))
            return result
        finally:
            with self.lock:
                self.pending.pop(rid, None)
                self.busy -= 1

    def stop(self) -> None:
        proc = self.process
        if proc and proc.poll() is None:
            proc.kill()
            if self.command[:2] == ["wsl", "--exec"]:
                # Killing wsl.exe does not always end the program inside it.
                subprocess.run(["wsl", "--exec", "pkill", "-f", self.command[2]], capture_output=True)

    def describe(self) -> dict:
        t = self.target or {}
        return {"name": self.name, "remote": self.remote, "alive": self.alive(), "busy": self.busy, "served": self.served,
                "seconds": round(self.seconds, 3), "brand": t.get("brand"), "threads": t.get("threads"),
                "built_for": t.get("built_for"), "lane_bytes": t.get("lane_bytes"), "best": t.get("best"), "error": self.error}


class Cluster:
    """The workers, started when first needed; requests to the least busy; jobs cut into shards."""

    def __init__(self, local: int = 2):
        self.local = max(1, local)
        self.workers: list[Worker] = []
        self.lock = threading.Lock()
        self.last_job: dict | None = None
        self.remote_config: list[dict] = []
        self.holders: dict[str, list[str]] = {}  # the genomes each worker was last given

    def configure(self) -> None:
        try:
            config = json.loads(CONFIG.read_text(encoding="utf-8"))
            self.remote_config = [w for w in config.get("workers", []) if isinstance(w.get("command"), list)]
        except (OSError, ValueError):
            self.remote_config = []

    def ensure(self) -> list[Worker]:
        """Starts the workers that are not running (the remote ones best effort)."""
        with self.lock:
            if not self.workers:
                self.configure()
                for k in range(self.local):
                    self.workers.append(Worker(f"local-{k + 1}", engine.worker_command()))
                for w in self.remote_config:
                    self.workers.append(Worker(w.get("name") or "remote", w["command"], remote=True))
            for w in self.workers:
                if not w.alive():
                    try:
                        w.start()
                    except (OSError, WorkerError) as error:
                        w.error = str(error)
            alive = [w for w in self.workers if w.alive()]
        if not alive:
            raise WorkerError("no engine worker could be started: " + "; ".join(f"{w.name}: {w.error}" for w in self.workers))
        return alive

    def restart(self) -> None:
        """Stops every worker, so that the next request starts them from the present build."""
        with self.lock:
            for w in self.workers:
                w.stop()
            self.workers = []
            self.holders = {}

    def request(self, request: dict, timeout: float = 120, local_only: bool = True) -> dict:
        workers = [w for w in self.ensure() if not (local_only and w.remote)]
        if not workers:
            raise WorkerError("no local worker is running")
        # A worker that has the genome in memory is preferred while it is not much busier.
        key = json.dumps(request.get("genome"), sort_keys=True) if request.get("genome") is not None else None
        holder = next((w for w in workers if key and key in self.holders.get(w.name, ())), None)
        least = min(workers, key=lambda w: (w.busy, w.served))
        worker = holder if holder is not None and holder.busy <= least.busy + 1 else least
        if key:
            held = self.holders.setdefault(worker.name, [])
            if key not in held:
                held.append(key)
                del held[:-6]  # a worker keeps six genomes (genomes.hpp, Worker::KEEP)
        result = worker.request(request, timeout)
        if isinstance(result, dict):
            result["_worker"] = worker.name
        return result

    def job(self, op: str, params: dict, timeout: float = 600, local_only: bool = True) -> dict:
        """A sharded job: plan, shards pulled by free workers, merge; with its trace."""
        t0 = time.time()
        plan = self.request({**params, "op": f"{op}.plan"}, timeout=timeout, local_only=True)
        shards = plan.get("shards", [])
        workers = [w for w in self.ensure() if not (local_only and w.remote)]
        todo: queue.Queue = queue.Queue()
        for k, s in enumerate(shards):
            todo.put((k, s))
        parts: list = [None] * len(shards)
        trace: list[dict] = []
        failures: list[str] = []
        trace_lock = threading.Lock()

        def pull(worker: Worker) -> None:
            while not failures:
                try:
                    k, shard = todo.get_nowait()
                except queue.Empty:
                    return
                start = time.time() - t0
                try:
                    part = worker.request({**params, "op": f"{op}.shard", "shard": shard}, timeout=timeout)
                except WorkerError as error:
                    if worker.alive():
                        failures.append(f"shard {k} on {worker.name}: {error}")
                        return
                    todo.put((k, shard))  # the worker died: another takes the shard
                    return
                end = time.time() - t0
                inner = part.pop("_trace", None) if isinstance(part, dict) else None
                seconds = part.pop("_seconds", None) if isinstance(part, dict) else None
                parts[k] = part
                with trace_lock:
                    trace.append({"shard": k, "worker": worker.name, "start": round(start, 4), "end": round(end, 4),
                                  "engine_seconds": seconds, "size": shard.get("size") if isinstance(shard, dict) else None, "threads": inner})

        # Rounds of pulling: a shard put back by a worker that ended is taken by one that lives.
        while not todo.empty() and not failures:
            live = [w for w in workers if w.alive()]
            if not live:
                break
            threads = [threading.Thread(target=pull, args=(w,), daemon=True) for w in live]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
        if failures:
            raise WorkerError(failures[0])
        if any(p is None for p in parts):
            raise WorkerError("some shards were not computed: every worker that took them ended")
        merge_start = time.time() - t0
        result = self.request({**params, "op": f"{op}.merge", "parts": parts, "plan": plan}, timeout=timeout)
        end = time.time() - t0
        job = {"op": op, "shards": len(shards), "seconds": round(end, 4), "plan_seconds": plan.get("_seconds"),
               "merge": {"start": round(merge_start, 4), "end": round(end, 4), "worker": result.get("_worker")},
               "workers": [w.describe() for w in workers], "trace": sorted(trace, key=lambda r: r["start"])}
        self.last_job = job
        if isinstance(result, dict):
            result["_job"] = job
        return result

    def describe(self) -> dict:
        return {"local": self.local, "configured_remote": len(self.remote_config), "config": str(CONFIG),
                "workers": [w.describe() for w in self.workers], "last_job": self.last_job}

    def stop(self) -> None:
        for w in self.workers:
            w.stop()
