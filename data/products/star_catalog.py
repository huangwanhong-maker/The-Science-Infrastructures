"""Derive compact star catalogues for the applications from registered datasets.

Sources:
    hyg-v41   data/cache/hyg-v41/hygdata_v41.csv          -> products/stars-v1.*
    gaia-dr3  data/cache/gaia-dr3/<profile>/result.csv     -> products/gaia-dr3-<profile>.*

Each product is a pair (or triple) of files in data/cache/products/:
    <name>.bin       little-endian float32 records of eight values:
                     ra, dec (radians, ICRS); mag; colour index; dist (pc, NaN if unknown);
                     x, y, z (equatorial, pc, NaN if the distance is unknown)
    <name>.json      layout, photometric system, names, designations, spectra, provenance
    <name>.ids.bin   Gaia only: little-endian uint64 source_id per record

The photometric system is 'johnson' (mag = V, colour = B-V) for HYG and 'gaia'
(mag = G, colour = BP-RP) for Gaia. Gaia stars receive proper names by
cross-matching HYG's named stars, allowing for proper motion between J2000 and
J2016. data/cache/products/index.json lists the star products for applications.

Usage::

    python data/products/star_catalog.py                                   # HYG
    python data/products/star_catalog.py --source gaia-dr3 --profile small
    python data/products/star_catalog.py --source gaia-dr3 --profile custom-1a2b3c4d
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import struct
import sys
from array import array
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import register  # noqa: E402

DATA = Path(__file__).resolve().parents[1]
CACHE = DATA / "cache"
OUT = CACHE / "products"
HYG = CACHE / "hyg-v41" / "hygdata_v41.csv"
FIELDS = ["ra", "dec", "mag", "ci", "dist", "x", "y", "z"]
NAN = math.nan

UNKNOWN_DISTANCE = 100000.0          # HYG's marker for an unknown distance
SPECTRA_WITHIN_PC = 100.0
DESIGNATE_BRIGHTER_THAN = 6.5
DESIGNATE_WITHIN_PC = 25.0
GREEK = {
    "Alp": "α", "Bet": "β", "Gam": "γ", "Del": "δ", "Eps": "ε", "Zet": "ζ", "Eta": "η",
    "The": "θ", "Iot": "ι", "Kap": "κ", "Lam": "λ", "Mu": "μ", "Nu": "ν", "Xi": "ξ",
    "Omi": "ο", "Pi": "π", "Rho": "ρ", "Sig": "σ", "Tau": "τ", "Ups": "υ", "Phi": "φ",
    "Chi": "χ", "Psi": "ψ", "Ome": "ω",
}
SUPERSCRIPT = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def number(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        return NAN


def record(ra_deg: float, dec_deg: float, mag: float, colour: float, dist: float) -> bytes:
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    if math.isfinite(dist):
        c = math.cos(dec)
        x, y, z = dist * c * math.cos(ra), dist * c * math.sin(ra), dist * math.sin(dec)
    else:
        x = y = z = NAN
    return struct.pack("<8f", ra, dec, mag, colour, dist, x, y, z)


def designation(row: dict) -> str:
    """Bayer, Flamsteed, Gliese or Hipparcos designation, in that order of preference."""
    con = row["con"].strip()
    bayer = row["bayer"].strip()
    if bayer:
        letter, _, index = bayer.partition("-")
        return f"{GREEK.get(letter, letter)}{index.translate(SUPERSCRIPT)} {con}".strip()
    if row["flam"].strip():
        return f"{row['flam'].strip()} {con}".strip()
    if row["gl"].strip():
        return row["gl"].strip()
    if row["hip"].strip():
        return f"HIP {row['hip'].strip()}"
    return ""


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---- HYG --------------------------------------------------------------------

def build_hyg() -> tuple[str, bytes, dict]:
    if not HYG.exists():
        raise SystemExit(f"missing {HYG}; run: python data/acquisition/fetch.py hyg-v41")
    records = bytearray()
    names, designations, spectra = [], [], []
    count = 0
    with HYG.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["id"] == "0":
                continue
            mag, dist = number(row["mag"]), number(row["dist"])
            known = math.isfinite(dist) and 0.0 < dist < UNKNOWN_DISTANCE
            records += record(number(row["ra"]) * 15.0, number(row["dec"]), mag, number(row["ci"]),
                              dist if known else NAN)
            if row["proper"].strip():
                names.append([count, row["proper"].strip()])
            if mag <= DESIGNATE_BRIGHTER_THAN or (known and dist <= DESIGNATE_WITHIN_PC):
                label = designation(row)
                if label:
                    designations.append([count, label])
            if known and dist <= SPECTRA_WITHIN_PC and row["spect"].strip():
                spectra.append([count, row["spect"].strip()])
            count += 1
    meta = {
        "title": "HYG v4.1 (Hipparcos, Yale, Gliese)",
        "description": "About 120 000 stars, complete to about magnitude 7.5; Hipparcos distances",
        "count": count, "photometry": "johnson", "epoch": "J2000",
        "names": names, "designations": designations, "spectra": spectra,
        "source": {"dataset": "hyg-v41", "file": HYG.name, "sha256": sha256_of(HYG), "licence": "CC BY-SA 4.0"},
    }
    return "stars-v1", bytes(records), meta


# ---- Gaia -------------------------------------------------------------------

def hyg_named_grid() -> dict:
    """HYG's named stars in 1-degree cells, for cross-matching."""
    grid = {}
    if not HYG.exists():
        return grid
    with HYG.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = row["proper"].strip()
            if not name or row["id"] == "0":
                continue
            ra, dec, mag = number(row["ra"]) * 15.0, number(row["dec"]), number(row["mag"])
            grid.setdefault((math.floor(dec), math.floor(ra)), []).append((ra, dec, mag, name))
    return grid


def separation_arcsec(ra1, dec1, ra2, dec2) -> float:
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    s = math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) * math.sin((r2 - r1) / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(s)))) * 3600


def build_gaia(profile: str) -> tuple[str, bytes, dict, bytes]:
    folder = CACHE / "gaia-dr3" / profile
    source, provenance = folder / "result.csv", folder / "retrieval.json"
    if not source.exists():
        raise SystemExit(f"missing {source}; run: python data/acquisition/fetch.py gaia-dr3 --profile {profile}")
    retrieval = json.loads(provenance.read_text(encoding="utf-8")) if provenance.exists() else {}
    grid = hyg_named_grid()
    records, ids = bytearray(), array("Q")
    names, taken = [], set()
    count = 0
    with source.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            ra, dec = float(row["ra"]), float(row["dec"])
            g, bp_rp = number(row["phot_g_mean_mag"]), number(row["bp_rp"])
            dist = number(row["r_med_geo"])
            if not (math.isfinite(dist) and dist > 0):
                parallax = number(row["parallax"])
                dist = 1000.0 / parallax if parallax > 0 else NAN
            records += record(ra, dec, g, bp_rp, dist)
            ids.append(int(row["source_id"]))
            if g < 9.5 and grid:
                pm = math.hypot(number(row["pmra"]) or 0.0, number(row["pmdec"]) or 0.0) / 1000.0
                tolerance = 10.0 + 16.0 * (pm if math.isfinite(pm) else 0.0)   # J2000 to J2016
                best = None
                for ddec in (-1, 0, 1):
                    for dra in (-1, 0, 1):
                        for hra, hdec, hmag, name in grid.get((math.floor(dec) + ddec, (math.floor(ra) + dra) % 360), ()):
                            sep = separation_arcsec(ra, dec, hra, hdec)
                            if sep < tolerance and abs(g - hmag) < 1.5 and (best is None or sep < best[0]):
                                best = (sep, name)
                if best and best[1] not in taken:
                    names.append([count, best[1]])
                    taken.add(best[1])
            count += 1
    if sys.byteorder != "little":
        ids.byteswap()
    meta = {
        "title": f"Gaia DR3 ({profile})",
        "description": ", ".join(f"{k}={v}" for k, v in retrieval.get("parameters", {}).items()),
        "count": count, "photometry": "gaia", "epoch": "J2016.0",
        "names": names, "designations": [], "spectra": [],
        "source": {"dataset": "gaia-dr3", "profile": profile, "parameters": retrieval.get("parameters"),
                   "file": source.name, "sha256": sha256_of(source), "licence": "CC BY-SA 3.0 IGO"},
    }
    return f"gaia-dr3-{profile}", bytes(records), meta, ids.tobytes()


# ---- output -----------------------------------------------------------------

def write(name: str, records: bytes, meta: dict, ids: bytes | None = None) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.bin").write_bytes(records)
    meta = {"name": name, **meta, "fields": FIELDS, "record_bytes": 4 * len(FIELDS),
            "built_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
    if ids is not None:
        (OUT / f"{name}.ids.bin").write_bytes(ids)
        meta["ids"] = f"{name}.ids.bin"
    (OUT / f"{name}.json").write_text(json.dumps(meta, ensure_ascii=False), encoding="utf-8")
    register("stars", {k: meta[k] for k in ("name", "title", "description", "count", "photometry", "epoch")},
             first="stars-v1")
    print(f"{name}: {meta['count']} stars, {len(meta['names'])} named -> {OUT}")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default="hyg-v41", choices=["hyg-v41", "gaia-dr3"])
    parser.add_argument("--profile", default="small", help="Gaia profile or custom tag (folder in data/cache/gaia-dr3)")
    args = parser.parse_args(argv)
    if args.source == "hyg-v41":
        write(*build_hyg())
    else:
        write(*build_gaia(args.profile))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
