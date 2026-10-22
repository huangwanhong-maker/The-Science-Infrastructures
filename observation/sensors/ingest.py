"""Ingest readings from a file into the sensor store, and build the product.

A file is CSV or JSON lines:

- CSV: a header, then one reading a row. The columns are station, lon, lat,
  time (ISO 8601, UTC unless an offset is given, or milliseconds since 1970),
  then one column per variable, named as the source's variables are; optional
  columns name (the station's) and elevation (m). An empty cell is a missing
  value.
- JSON lines (.jsonl, .ndjson): one reading a line, as store.append takes it:
  {"station", "lon", "lat", "elevation"?, "time", "values": {variable: value}}.

Values are in the units the source declares for its variables (degree_C, hPa
and others are converted to SI on arrival; see store.CONVERSIONS). A source
not yet defined is defined from --define and --variable; failing those, from
the CSV header, when every variable column is a CF standard name the store
knows, taken in its SI unit.

With --follow the file is watched as it grows, as a logger or a weather
station's software writes it, and each complete new line is ingested; a file
that shrinks (rotated or rewritten) is read again from its start. The product
is built after the file is read, and while following at most every
--build-every seconds.

Usage::

    python observation/sensors/ingest.py garden --define title="Garden station"         --variable air_temperature:K --variable relative_humidity:% --variable air_pressure:Pa   # define only
    python observation/sensors/ingest.py garden readings.csv
    python observation/sensors/ingest.py garden logger.csv --follow \\
        --define title="Garden station" --variable air_temperature:degree_C --variable relative_humidity:%
    python observation/sensors/ingest.py tide-pole tide.jsonl --define title="Tide pole, Kamakura" \\
        --variable water_surface_height_above_reference_datum:cm:"Water level above the pole's zero"
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import store  # noqa: E402

FIXED = {"station", "lon", "lat", "time", "name", "elevation"}


def number(text: str):
    text = (text or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return text            # kept, so that validation names the bad value


def csv_readings(header: list[str], rows) -> list[dict]:
    out = []
    for row in rows:
        if not row or not any(cell.strip() for cell in row):
            continue
        cell = dict(zip(header, row))
        reading = {"station": cell.get("station", "").strip(), "lon": number(cell.get("lon")),
                   "lat": number(cell.get("lat")), "time": cell.get("time", "").strip(),
                   "values": {k: number(v) for k, v in cell.items() if k not in FIXED}}
        if cell.get("elevation", "").strip():
            reading["elevation"] = number(cell["elevation"])
        if cell.get("name", "").strip():
            reading["name"] = cell["name"].strip()
        if reading["time"] and reading["time"].lstrip("-").replace(".", "", 1).isdigit():
            reading["time"] = float(reading["time"])          # milliseconds since 1970
        out.append(reading)
    return out


def json_readings(lines: list[str]) -> tuple[list[dict], int]:
    out, bad = [], 0
    for line in lines:
        if not line.strip():
            continue
        try:
            out.append(json.loads(line))
        except ValueError:
            bad += 1
    return out, bad


def is_csv(path: Path) -> bool:
    return path.suffix.lower() not in (".jsonl", ".ndjson", ".json")


def ensure_source(name: str, args, header: list[str] | None) -> None:
    """Defines the source from the options, or from the CSV header if it is new and none are given."""
    about = dict(item.split("=", 1) for item in args.define) if args.define else {}
    if args.variable or about:
        current = {s["name"]: s for s in store.sources()}.get(name)
        variables = args.variable or ([] if current else None)
        if variables is None:
            raise SystemExit(f"{name}: a new source needs its variables (--variable name:units)")
        store.define(name, about.get("title") or (current or {}).get("title") or name, variables,
                     credit=about.get("credit"), licence=about.get("licence"), domain=about.get("domain"))
        return
    if any(s["name"] == name for s in store.sources()):
        return
    columns = [c for c in (header or []) if c not in FIXED]
    unknown = [c for c in columns if c not in store.KNOWN]
    if not columns or unknown:
        raise SystemExit(f"{name} is not defined; give --define title=... and --variable name:units "
                         f"({'the columns ' + ', '.join(unknown) + ' are not CF names the store knows' if unknown else 'no variable columns'})")
    store.define(name, name, columns)
    print(f"{name}: defined from the header, in SI units: {', '.join(columns)}")


def ingest(name: str, readings: list[dict], where: str) -> int:
    rejected = []
    accepted = store.append(name, readings, rejected)
    print(f"{where}: {accepted} of {len(readings)} reading(s) accepted")
    for i, why in rejected[:10]:
        print(f"  refused reading {i + 1}: {why}")
    if len(rejected) > 10:
        print(f"  and {len(rejected) - 10} more refused")
    return accepted


def build(name: str, hours: float) -> None:
    path = store.build_product(name, hours)
    product = json.loads(path.read_text(encoding="utf-8"))
    print(f"{path.name}: {len(product['rows'])} rows, {product['window']['start']} to {product['window']['end']}")


def follow(name: str, path: Path, args) -> None:
    """Tails the file: each complete new line is ingested as it is written."""
    position, header, pending, last_build, dirty = 0, None, "", 0.0, False
    print(f"following {path} (Ctrl-C to stop)")
    try:
        while True:
            size = path.stat().st_size if path.exists() else 0
            if size < position:                                   # rotated or rewritten
                print(f"{path.name} shrank; reading it again from the start")
                position, header, pending = 0, None, ""
            if size > position:
                with path.open("rb") as handle:
                    handle.seek(position)
                    data = handle.read(size - position)
                position = size
                text = pending + data.decode("utf-8", errors="replace")
                lines = text.split("\n")
                pending = lines.pop()                             # a line not yet ended waits
                lines = [line.rstrip("\r") for line in lines]
                if is_csv(path):
                    if header is None and lines:
                        header = [h.strip() for h in next(csv.reader([lines.pop(0).lstrip("﻿")]))]
                        ensure_source(name, args, header)
                    readings = csv_readings(header or [], csv.reader(lines)) if lines else []
                else:
                    readings, bad = json_readings(lines)
                    if bad:
                        print(f"  {bad} line(s) are not JSON")
                if readings and ingest(name, readings, path.name):
                    dirty = True
            if dirty and time.monotonic() - last_build >= args.build_every:
                build(name, args.hours)
                last_build, dirty = time.monotonic(), False
            time.sleep(args.poll)
    except KeyboardInterrupt:
        if dirty:
            build(name, args.hours)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", help="the source's name: lowercase, digits and hyphens")
    parser.add_argument("file", type=Path, nargs="?", help="a .csv file, or JSON lines (.jsonl); none: only define the source")
    parser.add_argument("--follow", action="store_true", help="watch the file as it grows")
    parser.add_argument("--define", action="append", default=[], metavar="KEY=VALUE",
                        help="describe the source: title=..., and optionally credit=..., licence=..., domain=...")
    parser.add_argument("--variable", action="append", default=[], metavar="NAME:UNITS[:LONG NAME]",
                        help="a variable and the units its values arrive in, as air_temperature:degree_C (repeatable)")
    parser.add_argument("--hours", type=float, default=48, help="the product's window, in hours (default 48)")
    parser.add_argument("--poll", type=float, default=2.0, help="seconds between looks at a followed file")
    parser.add_argument("--build-every", type=float, default=60.0, help="seconds between builds while following")
    args = parser.parse_args(argv)
    for item in args.define:
        if "=" not in item:
            raise SystemExit(f"--define expects KEY=VALUE, got {item!r}")
    try:
        store._folder(args.source)
    except ValueError as error:
        raise SystemExit(str(error)) from None
    if args.file is None:
        if not (args.define or args.variable):
            raise SystemExit("give a file to ingest, or --define and --variable to define the source")
        ensure_source(args.source, args, None)
        desc = {s["name"]: s for s in store.sources()}[args.source]
        print(f"{args.source}: {desc['title']}; " + ", ".join(
            f"{v['name']} ({v.get('input_units', v['units'])})" for v in desc["variables"]))
        return 0
    if args.follow:
        if not is_csv(args.file):
            ensure_source(args.source, args, None)
        follow(args.source, args.file, args)
        return 0
    if not args.file.exists():
        raise SystemExit(f"no file {args.file}")
    text = args.file.read_text(encoding="utf-8-sig")
    if is_csv(args.file):
        rows = csv.reader(io.StringIO(text))
        header = [h.strip() for h in next(rows, [])]
        ensure_source(args.source, args, header)
        readings = csv_readings(header, rows)
    else:
        ensure_source(args.source, args, None)
        readings, bad = json_readings(text.splitlines())
        if bad:
            print(f"{bad} line(s) are not JSON")
    accepted = ingest(args.source, readings, args.file.name)
    if accepted or any(s["name"] == args.source for s in store.sources()):
        build(args.source, args.hours)
    return 0 if accepted or not readings else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
