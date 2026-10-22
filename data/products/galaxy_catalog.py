"""Derive galaxy catalogues for the applications from registered datasets.

Sources:
    ungc                 data/cache/ungc/<profile>/result.csv
                         -> products/galaxies-ungc.json
    tully-2mass-groups   data/cache/tully-2mass-groups/<profile>/result.csv
                         -> products/galaxies-2mrs-groups.{json,bin,groups.bin}

The Local Volume catalogue (UNGC) is small and is written as JSON columns; the
application computes positions, tidal indices and groups from it.

The 2MRS group catalogue is written as two little-endian float32 tables:
    .bin         per galaxy: SGX, SGY, SGZ (h^-1 Mpc), log L_K, V_CMB (km/s),
                 distance (h^-1 Mpc), morphological type, Ks, PGC, group index
    .groups.bin  per group:  centre SGX, SGY, SGZ, distance, log L_K (group),
                 M_lum, M_vir (10^12 solar masses), sigma_V (km/s), R_2t (Mpc),
                 members here, members in the catalogue, nest id, PGC of the
                 brightest member, index of the brightest member here (-1 if absent)
Integers are stored exactly: every one is below 2^24.

Names come from SIMBAD (CDS) when it can be reached: the brightest member of
each group with four or more members, and catalogued galaxy clusters ("Virgo
Cluster", ...) whose position lies within a rich group's turnaround radius and
whose velocity agrees with the group's. Results are cached beside the source
with their provenance, so later builds need no network.

Usage::

    python data/products/galaxy_catalog.py --source ungc
    python data/products/galaxy_catalog.py --source tully-2mass-groups [--profile near] [--no-names]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import re
import struct
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import register  # noqa: E402

DATA = Path(__file__).resolve().parents[1]
CACHE = DATA / "cache"
OUT = CACHE / "products"
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
NAN = math.nan
EXACT_INTEGER_LIMIT = 2 ** 24
CLUSTER_NAME = re.compile(r"^NAME ([A-Z][A-Za-z-]+) Cluster$")


def number(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return NAN


def json_number(x: float):
    return None if not math.isfinite(x) else x


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_rows(dataset: str, profile: str) -> tuple[Path, list[dict], dict]:
    folder = CACHE / dataset / profile
    path = folder / "result.csv"
    if not path.exists():
        raise SystemExit(f"missing {path}; run: python data/acquisition/fetch.py {dataset} --profile {profile}")
    retrieval = folder / "retrieval.json"
    meta = json.loads(retrieval.read_text(encoding="utf-8")) if retrieval.exists() else {}
    with path.open(newline="", encoding="utf-8") as handle:
        return path, list(csv.DictReader(handle)), meta


# ---- UNGC -------------------------------------------------------------------

def build_ungc(profile: str) -> None:
    path, rows, retrieval = source_rows("ungc", profile)
    columns = {
        "name": [r["Name"].strip() for r in rows],
        "ra": [number(r["RAJ2000"]) for r in rows],
        "dec": [number(r["DEJ2000"]) for r in rows],
        "dist": [number(r["Dist"]) for r in rows],
        "distMethod": [r["f_Dist"].strip() for r in rows],
        "logL": [number(r["KLum"]) for r in rows],
        "absB": [number(r["BMag"]) for r in rows],
        "type": [number(r["TT"]) for r in rows],
        "dwarfType": [r["Mcl"].strip() for r in rows],
        "vHelio": [number(r["HRV"]) for r in rows],
        "vLG": [number(r["Vlg"]) for r in rows],
        "theta1": [number(r["Ti1"]) for r in rows],
        "disturber": [r["MD"].strip() for r in rows],
        "theta5": [number(r["Ti5"]) for r in rows],
        "thetaJ": [number(r["Tij"]) for r in rows],
        "logMass26": [number(r["M26"]) for r in rows],
        "logMHI": [number(r["MHI"]) for r in rows],
        "diameter": [number(r["A26"]) for r in rows],
    }
    columns = {k: [json_number(v) if isinstance(v, float) else v for v in values] for k, values in columns.items()}
    name = "galaxies-ungc" if profile == "all" else f"galaxies-ungc-{profile}"
    meta = {
        "name": name,
        "title": "Local Volume (UNGC, Karachentsev et al. 2013)",
        "description": "Galaxies within about 11 Mpc with measured distances, tidal indices and main disturbers",
        "kind": "tidal", "distances": "measured", "units": "Mpc",
        "count": len(rows),
        "columns": columns,
        "source": {"dataset": "ungc", "profile": profile, "parameters": retrieval.get("parameters"),
                   "sha256": sha256_of(path), "citation": "Karachentsev, Makarov & Kaisina (2013), AJ 145, 101"},
        "built_at": now(),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(meta, ensure_ascii=False), encoding="utf-8")
    register("galaxies", {k: meta[k] for k in ("name", "title", "description", "kind", "distances", "units", "count")},
             first="galaxies-ungc")
    print(f"{name}: {len(rows)} galaxies -> {OUT}")


# ---- 2MRS groups --------------------------------------------------------------

def simbad(adql: str) -> list[dict]:
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql}).encode()
    request = urllib.request.Request(SIMBAD, data=data, headers={"User-Agent": "laboratory-products/1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        text = response.read().decode("utf-8", "replace")
    return list(csv.DictReader(text.splitlines()))


def tidy(identifier: str) -> str:
    return re.sub(r"\s+", " ", identifier.replace("NAME ", "", 1)).strip()


def resolve_names(rows: list[dict], groups: list[dict], cache_path: Path) -> dict:
    """Names of brightest members and clusters, from SIMBAD, cached with provenance."""
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))

    by_pgc = {r["PGC"]: r for r in rows}
    wanted = sorted({g["pgc1"] for g in groups if g["n_here"] >= 4 and g["pgc1"] in by_pgc}, key=int)
    members = {}
    for start in range(0, len(wanted), 300):
        chunk = wanted[start:start + 300]
        ids = ", ".join(f"'PGC {p}'" for p in chunk)
        for r in simbad(f"SELECT i.id, b.main_id FROM ident AS i JOIN basic AS b ON b.oid = i.oidref WHERE i.id IN ({ids})"):
            pgc = re.sub(r"\D", "", r["id"])
            members[pgc] = tidy(r["main_id"])
    print(f"  SIMBAD: {len(members)} brightest-member names")

    # Named clusters: ClG objects with a 'NAME <Word> Cluster' identifier, inside a
    # rich group's turnaround radius and within 1500 km/s of its velocity.
    rich = sorted((g for g in groups if g["n_here"] >= 15 and g["pgc1"] in by_pgc), key=lambda g: -g["n_here"])[:80]
    clusters = {}
    for start in range(0, len(rich), 20):
        cones = []
        for g in rich[start:start + 20]:
            r = by_pgc[g["pgc1"]]
            radius = max(1.0, min(8.0, math.degrees(math.atan2(g["r2t"], g["dist"]))) if g["dist"] > 0 else 1.0)
            cones.append(f"CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {r['_RA_icrs']}, {r['_DE_icrs']}, {radius:.3f})) = 1")
        query = ("SELECT b.main_id, b.ra, b.dec, b.rvz_radvel, i.id FROM basic AS b JOIN ident AS i ON i.oidref = b.oid "
                 f"WHERE b.otype = 'ClG' AND i.id LIKE 'NAME %' AND ({' OR '.join(cones)})")
        for r in simbad(query):
            m = CLUSTER_NAME.match(re.sub(r"\s+", " ", r["id"]).strip())
            if not m:
                continue
            entry = clusters.setdefault(r["main_id"], {"names": set(), "ra": number(r["ra"]), "dec": number(r["dec"]),
                                                       "v": number(r["rvz_radvel"])})
            entry["names"].add(m.group(1))
    named = {}
    for main_id, c in clusters.items():
        if not math.isfinite(c["v"]):
            continue
        best = None
        for g in rich:
            r = by_pgc[g["pgc1"]]
            ra, dec = number(r["_RA_icrs"]), number(r["_DE_icrs"])
            sep = math.degrees(math.acos(max(-1.0, min(1.0,
                math.sin(math.radians(dec)) * math.sin(math.radians(c["dec"]))
                + math.cos(math.radians(dec)) * math.cos(math.radians(c["dec"])) * math.cos(math.radians(ra - c["ra"]))))))
            radius = max(1.0, min(8.0, math.degrees(math.atan2(g["r2t"], g["dist"]))))
            if sep <= radius and abs(c["v"] - g["v"]) <= 1500 and (best is None or sep < best[0]):
                best = (sep, g["nest"])
        if best:
            word = max(c["names"], key=len)
            previous = named.get(best[1])
            if previous is None or best[0] < previous[0]:
                named[best[1]] = (best[0], f"{word} Cluster")
    result = {
        "members": members,
        "clusters": {str(nest): name for nest, (_, name) in named.items()},
        "service": SIMBAD,
        "retrieved_at": now(),
    }
    cache_path.write_text(json.dumps(result, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"  SIMBAD: clusters named: {', '.join(sorted(result['clusters'].values())) or 'none'}")
    return result


def build_groups(profile: str, names_wanted: bool) -> None:
    path, rows, retrieval = source_rows("tully-2mass-groups", profile)
    galaxies = bytearray()
    groups: dict[str, dict] = {}
    order: list[str] = []
    kept = []
    for r in rows:
        dist, lon, lat = number(r["Dist"]), number(r["SGLON"]), number(r["SGLAT"])
        if not (math.isfinite(dist) and math.isfinite(lon) and math.isfinite(lat)):
            continue
        kept.append(r)
        nest = r["Nest"].strip()
        if nest not in groups:
            glon, glat = number(r["GSGLON"]), number(r["GSGLAT"])
            if not (math.isfinite(glon) and math.isfinite(glat)):
                glon, glat = lon, lat
            cl, cb = math.radians(glon), math.radians(glat)
            groups[nest] = {
                "nest": nest, "index": len(order), "dist": dist,
                "centre": (dist * math.cos(cb) * math.cos(cl), dist * math.cos(cb) * math.sin(cl), dist * math.sin(cb)),
                "logLK": number(r["logLK"]), "mlum": number(r["Mlum"]), "mvir": number(r["Mvir"]),
                "sigV": number(r["sigV"]), "r2t": number(r["R2t"]), "n_catalog": int(number(r["Nmb"]) or 1),
                "pgc1": r["PGC1"].strip(), "n_here": 0, "brightest": -1, "v_sum": 0.0,
            }
            order.append(nest)
        g = groups[nest]
        g["n_here"] += 1
        g["v_sum"] += number(r["Vcmba"])
        index = len(kept) - 1
        if r["PGC"].strip() == g["pgc1"]:
            g["brightest"] = index
        l, b = math.radians(lon), math.radians(lat)
        pgc = int(r["PGC"])
        for value in (pgc, g["index"], len(order)):
            assert value < EXACT_INTEGER_LIMIT, value
        galaxies += struct.pack("<10f",
            dist * math.cos(b) * math.cos(l), dist * math.cos(b) * math.sin(l), dist * math.sin(b),
            number(r["logLKi"]), number(r["Vcmba"]), dist, number(r["MType"]), number(r["Kmag"]),
            pgc, g["index"])
    group_list = [groups[n] for n in order]
    for g in group_list:
        g["v"] = g["v_sum"] / g["n_here"]

    names = {"members": {}, "clusters": {}}
    if names_wanted:
        try:
            names = resolve_names(kept, group_list, path.parent / "simbad-names.json")
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            print(f"  SIMBAD unavailable ({error}); building without names")

    table = bytearray()
    named_groups = []
    for g in group_list:
        table += struct.pack("<14f", *g["centre"], g["dist"], g["logLK"], g["mlum"], g["mvir"], g["sigV"], g["r2t"],
                             g["n_here"], g["n_catalog"], int(g["nest"]), int(g["pgc1"] or 0), g["brightest"])
        label = names["clusters"].get(g["nest"]) or (
            f"{names['members'][g['pgc1']]} group" if g["pgc1"] in names["members"] and g["n_here"] >= 4 else None)
        if label:
            named_groups.append([g["index"], label])

    name = "galaxies-2mrs-groups" if profile == "all" else f"galaxies-2mrs-groups-{profile}"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.bin").write_bytes(bytes(galaxies))
    (OUT / f"{name}.groups.bin").write_bytes(bytes(table))
    member_names = {}
    for i, r in enumerate(kept):
        if r["PGC"] in names["members"]:
            member_names[str(i)] = names["members"][r["PGC"]]
    meta = {
        "name": name,
        "title": "Local Universe (2MRS groups, Tully 2015)",
        "description": "Galaxies brighter than Ks = 11.75 with their group membership; distances from velocities",
        "kind": "groups", "distances": "kinematic", "units": "h^-1 Mpc",
        "count": len(kept), "groups": len(group_list),
        "fields": ["sgx", "sgy", "sgz", "logL", "vCMB", "dist", "type", "kmag", "pgc", "group"],
        "group_fields": ["cx", "cy", "cz", "dist", "logLK", "mLum", "mVir", "sigmaV", "r2t",
                         "nHere", "nCatalog", "nest", "pgc1", "brightest"],
        "names": member_names,
        "groupNames": named_groups,
        "source": {"dataset": "tully-2mass-groups", "profile": profile, "parameters": retrieval.get("parameters"),
                   "sha256": sha256_of(path),
                   "citation": "Tully (2015), AJ 149, 171; Huchra et al. (2012), ApJS 199, 26",
                   "names": "SIMBAD (CDS)" if names["members"] else None},
        "built_at": now(),
    }
    (OUT / f"{name}.json").write_text(json.dumps(meta, ensure_ascii=False), encoding="utf-8")
    register("galaxies", {k: meta[k] for k in ("name", "title", "description", "kind", "distances", "units", "count")},
             first="galaxies-ungc")
    print(f"{name}: {len(kept)} galaxies in {len(group_list)} groups, {len(named_groups)} groups named -> {OUT}")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True, choices=["ungc", "tully-2mass-groups"])
    parser.add_argument("--profile", default="all")
    parser.add_argument("--no-names", action="store_true", help="do not query SIMBAD for names")
    args = parser.parse_args(argv)
    if args.source == "ungc":
        build_ungc(args.profile)
    else:
        build_groups(args.profile, not args.no_names)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
