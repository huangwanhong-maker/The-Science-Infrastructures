"""Readings as they arrive: a client of the service's stream (Server-Sent Events),
printing each new reading of a source, or of every source. It reconnects by
itself when the connection drops, and asks for what it missed.

    python api/examples/python/sensor_stream.py [example-garden] [--url ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import science  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("source", nargs="?", help="a sensor source (default: every source)")
parser.add_argument("--url")
args = parser.parse_args()

api = science.Client(args.url)
print(f"Waiting for readings of {args.source or 'every source'} (Ctrl+C to stop)", flush=True)
try:
    for event in api.stream(args.source):
        for r in event["readings"]:
            values = ", ".join(f"{k} {v:g}" for k, v in r["values"].items())
            print(f"{r['time']}  {event['source']}/{r['station']}  {values}", flush=True)
except KeyboardInterrupt:
    pass
