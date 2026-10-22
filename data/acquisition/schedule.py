"""Run an application's data on a schedule: observations every hour, forecasts every six hours, and so on.

A recipe (recipe.py) says what an application needs; a schedule says when to
fetch it again. A schedule file, kept beside the application
(``applications/<name>/schedule.json``), lists jobs, each a list of commands,
scripts of the data library run in order from the repository root with the
interpreter running this one, as a recipe's steps are::

    {"about": "...",
     "jobs": [{"name": "observations", "every": "1h", "offset_minutes": 10, "enabled": true,
               "timeout_minutes": 30,
               "run": [["data/acquisition/recipe.py", "applications/earth-atlas/recipe.json", "metar", "observations-metar"],
                       ["data/products/archive.py"]]}]}

**Due times.** A job of period P (``every``: minutes, hours or days, as "30m",
"6h", "1d") and offset O (``offset_minutes``) is due at the instants k P + O
from 1970-01-01T00:00:00Z, in UTC: "6h" with an offset of 300 minutes is due at
05:00, 11:00, 17:00 and 23:00 UTC. A job is due when the latest of its instants
at or before the present is later than its last start; a job that has never
run is due at once. Instants missed (a computer asleep, a long job before it)
give one run, not one for each.

**One job at a time.** The product builders all rewrite the library's index of
products, so jobs never run together: a lock file (data/cache/schedule.lock,
locked through the operating system, so that a runner that dies releases it)
keeps two runners, or a runner and ``--once``, apart, and a job due while
another runs waits for it. A job checks again that it is due once it holds
the lock, so that two runners do not run it twice.

**Failures.** A command that fails (an exit status other than 0, a script that
is missing, or a run longer than the job's ``timeout_minutes``, 120 unless
set, after which it is stopped with the processes it started) is logged, and
the job goes on with its next command; the job's status is then "failed", with
the commands that failed. recipe.py refuses a step it does not know and then
runs nothing, so a recipe.py command naming steps that its recipe does not (yet)
have runs without them, with a note, and is skipped if none is left.

**State**, in data/cache/schedule-state.json: for each schedule file and job,
its last start and end (UTC), its status (running, ok, failed, or interrupted
when the runner was stopped, which stops the command in progress too) and the
commands that failed. Every line of the log begins with the time in UTC.

Usage::

    python data/acquisition/schedule.py applications/earth-atlas/schedule.json                  # wait for each job's time
    python data/acquisition/schedule.py applications/earth-atlas/schedule.json --list           # jobs, last runs, next times
    python data/acquisition/schedule.py applications/earth-atlas/schedule.json --once forecasts # one job now, then exit
    python data/acquisition/schedule.py applications/earth-atlas/schedule.json --dry-run        # what is due, and its commands
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import re
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]            # data/acquisition -> the repository
STATE = ROOT / "data" / "cache" / "schedule-state.json"
LOCK = ROOT / "data" / "cache" / "schedule.lock"
RECIPE = "data/acquisition/recipe.py"
TIMEOUT_MINUTES = 120
POLL_SECONDS = 60                 # the longest sleep between looks at the schedule
UTC = dt.timezone.utc
PERIOD = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([mhd])\s*$")


def utcnow() -> dt.datetime:
    return dt.datetime.now(UTC)


def stamp(t: dt.datetime) -> str:
    return t.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(message: str) -> None:
    for line in str(message).split("\n"):
        print(f"{stamp(utcnow())}  {line}", flush=True)


# ---- due times --------------------------------------------------------------------

def period_seconds(every: str) -> float:
    """The period of a job, "30m", "1h", "6h" or "1d", in seconds."""
    m = PERIOD.match(every)
    if not m:
        raise ValueError(f"a period is a number and m, h or d, as '6h'; got {every!r}")
    return float(m.group(1)) * {"m": 60, "h": 3600, "d": 86400}[m.group(2)]


def latest_slot(every: str, offset_minutes: float, now: dt.datetime) -> dt.datetime:
    """The latest instant k P + O (from 1970-01-01T00:00Z) at or before now."""
    p, o = period_seconds(every), offset_minutes * 60
    k = math.floor((now.timestamp() - o) / p)
    return dt.datetime.fromtimestamp(k * p + o, UTC)


def next_slot(every: str, offset_minutes: float, now: dt.datetime) -> dt.datetime:
    """The first instant k P + O after now."""
    return latest_slot(every, offset_minutes, now) + dt.timedelta(seconds=period_seconds(every))


def is_due(job: dict, last_start: str | None, now: dt.datetime) -> bool:
    """Whether a job is due at now, given its last start (ISO 8601, or None if it never ran)."""
    if not job.get("enabled", True):
        return False
    if not last_start:
        return True
    last = dt.datetime.fromisoformat(last_start.replace("Z", "+00:00"))
    return latest_slot(job["every"], job.get("offset_minutes", 0), now) > last


# ---- state and the lock -------------------------------------------------------------------

def load_state(path: Path = STATE) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def save_state(state: dict, path: Path = STATE) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(state, indent=1) + "\n", encoding="utf-8")
    os.replace(part, path)


class Lock:
    """A lock on a file, held through the operating system (msvcrt on Windows, flock elsewhere), so that
    it is released when its process ends, however it ends."""

    def __init__(self, path: Path = LOCK):
        self.path, self.handle = path, None

    def _try(self) -> bool:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = open(self.path, "a+")
        try:
            if os.name == "nt":
                import msvcrt
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            handle.close()
            return False
        self.handle = handle
        return True

    def acquire(self, wait: bool = True, say=log) -> bool:
        if self._try():
            return True
        if not wait:
            return False
        say("another job is running; waiting for it to end")
        while not self._try():
            time.sleep(5)
        return True

    def release(self) -> None:
        if not self.handle:
            return
        try:
            if os.name == "nt":
                import msvcrt
                self.handle.seek(0)
                msvcrt.locking(self.handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        finally:
            self.handle.close()
            self.handle = None


# ---- commands ------------------------------------------------------------------------

def known_steps(command: list[str], root: Path = ROOT) -> tuple[list[str] | None, list[str]]:
    """A recipe.py command with only the steps its recipe has, and the steps left out; None if none is left.
    Any other command is returned as it is."""
    if len(command) < 3 or Path(command[0]).as_posix() != RECIPE:
        return command, []
    try:
        ids = {s["id"] for s in json.loads((root / command[1]).read_text(encoding="utf-8"))["steps"]}
    except (OSError, ValueError, KeyError):
        return command, []
    named = [a for a in command[2:] if not a.startswith("-")]
    unknown = [a for a in named if a not in ids]
    if not unknown:
        return command, []
    kept = [a for a in named if a in ids]
    if not kept:
        return None, unknown
    return [command[0], command[1], *kept, *(a for a in command[2:] if a.startswith("-"))], unknown


def stop_tree(process: subprocess.Popen) -> None:
    """Stops a process and those it started."""
    if os.name == "nt":
        subprocess.run(["taskkill", "/T", "/F", "/PID", str(process.pid)], capture_output=True)
    else:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except OSError:
            process.kill()


def run_command(command: list[str], root: Path = ROOT, timeout_s: float | None = None, say=log) -> int | None:
    """Runs a script from the repository root with this interpreter, passing its output to `say`;
    returns its exit status, or None if it was stopped for running too long."""
    if not (root / command[0]).exists():
        say(f"no such script: {command[0]}")
        return 127
    say(f"$ python {' '.join(command)}")
    options = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP} if os.name == "nt" else {"start_new_session": True}
    process = subprocess.Popen(
        [sys.executable, *command], cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        encoding="utf-8", errors="replace", bufsize=1,
        env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"}, **options)

    def relay():
        for line in process.stdout:
            say("  " + line.rstrip("\r\n").split("\r")[-1])      # a progress line keeps its last state
    reader = threading.Thread(target=relay, daemon=True)
    reader.start()
    try:
        status = process.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        stop_tree(process)
        process.wait()
        status = None
    except KeyboardInterrupt:          # the command is in a group of its own, and would outlive the runner
        stop_tree(process)
        raise
    reader.join(timeout=10)
    return status


def run_job(schedule_path: Path, job: dict, root: Path = ROOT, state_path: Path = STATE, lock: Lock | None = None,
            force: bool = False, say=log) -> str | None:
    """Runs a job's commands in order, holding the lock; returns its status ("ok" or "failed"),
    or None if, once the lock was held, the job was no longer due (another runner ran it)."""
    lock = lock or Lock(state_path.with_name(LOCK.name))
    key = schedule_key(schedule_path, root)
    lock.acquire(say=say)
    try:
        state = load_state(state_path)
        entry = state.setdefault(key, {}).get(job["name"], {})
        if not force and not is_due(job, entry.get("last_start"), utcnow()):
            return None
        started = utcnow()
        entry = {**entry, "last_start": stamp(started), "status": "running", "failed": []}
        state[key][job["name"]] = entry
        save_state(state, state_path)
        say(f"### {job['name']}: {len(job['run'])} command(s)")
        timeout = float(job.get("timeout_minutes", TIMEOUT_MINUTES)) * 60
        failed = []
        for command in job["run"]:
            chosen, unknown = known_steps(command, root)
            if unknown:
                say(f"not in the recipe (yet): {', '.join(unknown)}" + ("; skipped" if chosen is None else ""))
            if chosen is None:
                continue
            try:
                status = run_command(chosen, root, timeout, say)
            except OSError as error:
                say(f"could not run it: {error}")
                status = 126
            except KeyboardInterrupt:
                state = load_state(state_path)
                state.setdefault(key, {})[job["name"]] = {
                    "last_start": stamp(started), "last_end": stamp(utcnow()), "status": "interrupted",
                    "failed": failed + [" ".join(chosen)]}
                save_state(state, state_path)
                raise
            if status != 0:
                say(f"failed: {' '.join(chosen)} ({'stopped after ' + str(round(timeout / 60)) + ' min' if status is None else 'status ' + str(status)})")
                failed.append(" ".join(chosen))
        ended = utcnow()
        state = load_state(state_path)
        state.setdefault(key, {})[job["name"]] = {
            "last_start": stamp(started), "last_end": stamp(ended), "seconds": round((ended - started).total_seconds()),
            "status": "failed" if failed else "ok", "failed": failed}
        save_state(state, state_path)
        say(f"### {job['name']}: {'failed: ' + str(len(failed)) + ' command(s)' if failed else 'ok'}, "
            f"{(ended - started).total_seconds() / 60:.1f} min")
        return "failed" if failed else "ok"
    finally:
        lock.release()


# ---- the schedule -------------------------------------------------------------------------

def schedule_key(path: Path, root: Path = ROOT) -> str:
    path = Path(path).resolve()
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def load(path: Path) -> dict:
    schedule = json.loads(Path(path).read_text(encoding="utf-8"))
    names = [job["name"] for job in schedule["jobs"]]
    if len(set(names)) != len(names):
        raise ValueError("two jobs have the same name")
    for job in schedule["jobs"]:
        period_seconds(job["every"])
    return schedule


def due_jobs(schedule: dict, entries: dict, now: dt.datetime) -> list[dict]:
    """The jobs due at now, the one due the longest first."""
    due = [job for job in schedule["jobs"] if is_due(job, (entries.get(job["name"]) or {}).get("last_start"), now)]
    return sorted(due, key=lambda job: latest_slot(job["every"], job.get("offset_minutes", 0), now))


def show(schedule_path: Path, schedule: dict, state: dict, now: dt.datetime, commands: bool = False) -> None:
    entries = state.get(schedule_key(schedule_path), {})
    print(schedule.get("about", ""))
    for job in schedule["jobs"]:
        e = entries.get(job["name"]) or {}
        due = is_due(job, e.get("last_start"), now)
        when = "disabled" if not job.get("enabled", True) else ("due now" if due else
                                                                  f"next {stamp(next_slot(job['every'], job.get('offset_minutes', 0), now))}")
        last = f"last {e['last_start']} {e.get('status', '')}" if e.get("last_start") else "never run"
        print(f"  {job['name']:14s} every {job['every']:4s} +{job.get('offset_minutes', 0):>3} min  {when:28s} {last}")
        if commands and (due or not job.get("enabled", True)):
            for command in job["run"]:
                chosen, unknown = known_steps(command)
                note = f"   (not in the recipe yet: {', '.join(unknown)})" if unknown else ""
                print(f"      {'(skipped)' if chosen is None else '$ python ' + ' '.join(chosen)}{note}")


def serve(schedule_path: Path, root: Path = ROOT, state_path: Path = STATE) -> None:
    """Runs the jobs as they fall due, one at a time, until interrupted; the schedule is read again each time,
    so that a change to it takes effect without a restart."""
    schedule, key = load(schedule_path), schedule_key(schedule_path, root)
    log(f"schedule {key}: {', '.join(j['name'] + ' every ' + j['every'] for j in schedule['jobs'] if j.get('enabled', True))}")
    while True:
        try:
            schedule = load(schedule_path)
        except (OSError, ValueError, KeyError) as error:
            log(f"the schedule could not be read ({error}); keeping the last one")
        now = utcnow()
        due = due_jobs(schedule, load_state(state_path).get(key, {}), now)
        if due:
            run_job(schedule_path, due[0], root, state_path)
            continue
        enabled = [j for j in schedule["jobs"] if j.get("enabled", True)]
        wake = min((next_slot(j["every"], j.get("offset_minutes", 0), now) for j in enabled), default=None)
        time.sleep(POLL_SECONDS if wake is None else min(POLL_SECONDS, max(1.0, (wake - utcnow()).total_seconds())))


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("schedule", help="the schedule, such as applications/earth-atlas/schedule.json")
    parser.add_argument("--list", action="store_true", help="the jobs, their last runs and next times; run nothing")
    parser.add_argument("--once", metavar="JOB", help="run one job now (waiting for one in progress), then exit")
    parser.add_argument("--dry-run", action="store_true", help="show what is due now and its commands; run nothing")
    args = parser.parse_args(argv)
    # Consoles that cannot show ′, ° or CJK; and lines as they come when piped to a log.
    sys.stdout.reconfigure(errors="replace", line_buffering=True)

    path = Path(args.schedule)
    schedule = load(path)
    now = utcnow()
    if args.list or args.dry_run:
        if args.once:
            job = next((j for j in schedule["jobs"] if j["name"] == args.once), None)
            if not job:
                print(f"no job {args.once}", file=sys.stderr)
                return 2
            schedule = {**schedule, "jobs": [{**job, "enabled": True}]}
            show(path, schedule, {}, now, commands=args.dry_run)
        else:
            show(path, schedule, load_state(), now, commands=args.dry_run)
        return 0
    if args.once:
        job = next((j for j in schedule["jobs"] if j["name"] == args.once), None)
        if not job:
            print(f"no job {args.once}; see --list", file=sys.stderr)
            return 2
        return 0 if run_job(path, job, force=True) == "ok" else 1
    try:
        serve(path)
    except KeyboardInterrupt:
        log("stopped")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
