"""Serve the API on its own: the data library's products, the sensors'
observations, the domains' algorithms and the plug-ins, for programs in any
language (api/README.md). Only 127.0.0.1 is bound unless ``--host`` says otherwise.

Usage::

    python api/service/server.py [--host 127.0.0.1] [--port 8800] [--token T] [--read-only]
                                 [--compute-seconds 10] [--plugin-seconds 600] [--product-seconds 60]

then ask http://127.0.0.1:8800/api/v1/ for the endpoints. The token may come
from the environment variable SCIENCE_API_TOKEN instead of the command line. The
limits bound a call into the domains and a plug-in's run; posted readings are
built into their source's product at most every --product-seconds (0: never).
The checks also set --heartbeat-seconds and --sensor-store, and the environment
variables SCIENCE_OBSERVATIONS (the store's directory) and SCIENCE_PLUGIN_OUTPUT
(where runs write their products), so as to touch nothing of a person's own.

Writes (readings posted by sensors) are open on this computer; with a token they
need ``Authorization: Bearer <token>``; with ``--read-only`` they are refused. A
service on another address than the loopback, without a token, refuses writes:
secure by default. Reads and computations are open to every caller.
"""

from __future__ import annotations

import argparse
import http.server
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from handler import LOOPBACK, PREFIX, Api  # noqa: E402

API: Api | None = None


def write_policy(host: str, token: str | None, read_only: bool) -> tuple[bool, str | None]:
    """Whether the service accepts writes, and why not when it does not."""
    if read_only:
        return False, "this service is read-only (it was started with --read-only)"
    if token:
        return True, None
    if host not in LOOPBACK:
        return False, ("this service serves other computers without a token, so it accepts no writes: "
                       "start it with --token (or $SCIENCE_API_TOKEN)")
    return True, None


class Handler(http.server.BaseHTTPRequestHandler):
    server_version = "ScienceAPI/0.1"
    # A client that stops sending, or stops reading, is let go after this many
    # seconds; a stream's heartbeats keep a reading client well within it.
    timeout = 120

    def do_GET(self) -> None:
        if self.path in ("/", ""):
            self.send_response(302)
            self.send_header("Location", PREFIX + "/")
            self.send_header("Content-Length", "0")
            self.end_headers()
        elif not API.handle(self):
            self.not_found()

    def do_POST(self) -> None:
        if not API.handle(self):
            self.not_found()

    def do_OPTIONS(self) -> None:
        if not API.handle(self):
            self.not_found()

    def not_found(self) -> None:
        data = json.dumps({"error": f"no endpoint {self.path}: the API is at {PREFIX}/"}).encode("utf-8")
        self.send_response(404)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def end_headers(self) -> None:
        # Products are rebuilt in place when data are fetched again, and readings
        # arrive all the time: answers are revalidated rather than cached.
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, format: str, *args) -> None:
        if args and str(args[1]).startswith(("4", "5")):
            super().log_message(format, *args)


class Server(http.server.ThreadingHTTPServer):
    # Streams hold a thread each; clients may ask for much at once.
    request_queue_size = 128
    daemon_threads = True


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Serve the API of The Science Infrastructures.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8800)
    parser.add_argument("--token", default=os.environ.get("SCIENCE_API_TOKEN"),
                        help="the token that writes need (default: $SCIENCE_API_TOKEN)")
    parser.add_argument("--read-only", action="store_true", help="refuse every write")
    parser.add_argument("--compute-seconds", type=float, default=10.0, help="the time limit of a call into the domains")
    parser.add_argument("--plugin-seconds", type=float, default=600.0, help="the time limit of a plug-in's run")
    parser.add_argument("--product-seconds", type=float, default=60.0,
                        help="rebuild a sensor source's product at most this often after posts (0: never)")
    parser.add_argument("--heartbeat-seconds", type=float, default=15.0, help="seconds between a stream's heartbeats")
    parser.add_argument("--sensor-store", metavar="PATH",
                        help="a module with the interface of observation/sensors/store.py, used in its place (for tests)")
    args = parser.parse_args(argv)
    token = args.token or None
    allowed, reason = write_policy(args.host, token, args.read_only)
    global API
    API = Api(ROOT, allow_writes=allowed, token=token, read_only_reason=reason,
              store=Path(args.sensor_store).resolve() if args.sensor_store else None,
              compute_seconds=args.compute_seconds, plugin_seconds=args.plugin_seconds, product_seconds=args.product_seconds,
              heartbeat=args.heartbeat_seconds)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{server.server_address[1]}{PREFIX}/"
    where = url if args.host == "127.0.0.1" else f"{url} and on {args.host}"
    writes = "writes need the token" if token else ("writes are open on this computer" if allowed else reason)
    print(f"The Science API at {where}  ({writes}; Ctrl+C to stop)", flush=True)
    if not allowed and not args.read_only:
        print(reason, file=sys.stderr, flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        API.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
