"""The households and dwellings of the City of Melbourne, for the urban-science domain's housing (ADR 0075).

    python data/acquisition/fetch.py melbourne-clue-dwellings melbourne-clue-blocks homes-victoria-rents abs-census-2021-sa1-vic abs-asgs-2021-sa1
    python data/acquisition/fetch.py melbourne-development-activity --set south=-37.86 --set west=144.88 --set north=-37.77 --set east=145.0
    python data/products/housing_records.py

Writes data/cache/products/housing-melbourne.json, listed in the index under the kind `housing`, in the city model's
local frame (x east, y north, metres about its origin, by street_records.enu). Tables are kept as columns and rows.

1. **Blocks.** CLUE's 603 blocks: each block's id, small area, outer rings and point.
2. **Buildings.** CLUE's residential dwellings by property (a building) and type (House/Townhouse, Residential
   Apartments, Student Apartments) at every census year 2002-2024: the stock, never its occupants. A property's place
   is the mean of its rows' points; a property never given a point is placed at its block's point (counted).
3. **Rents.** Homes Victoria's moving annual medians of new lettings, quarterly from March 2000, for the suburb groups
   that hold the municipality's small areas and for each region's total (Inner Melbourne, the metropolis's other
   regions, the regional cities), every type of dwelling of the report (count and median a quarter; null where the
   report publishes none), and the map from small areas to groups (declared here: the groups cut across the
   municipality's boundary).
4. **Census.** The 2021 Census's SA1s whose point lies in a CLUE block: each one's ring and small area, and its
   households by income band (G33), by size (G35), its dwellings by structure and its unoccupied dwellings (G36), its
   tenure by landlord and structure (G37), its rents in bands (G40), its medians (G02) and its persons at a different
   address a year before (G44). Taken on 10 August 2021, in a lockdown.
5. **Developments.** The Development Activity Monitor's developments in the municipality, by status: completed (with
   the year), under construction, approved, applied for; their dwellings by bedrooms and student apartments and beds,
   their property ids and place.

The names of places are the sources' own names of the city's areas and suburbs; none is in dispute (naming.json).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402
from street_records import contains, enu, retrieved, write  # noqa: E402

DATA = Path(__file__).resolve().parents[1]
CACHE = DATA / "cache"
REGISTRY = DATA / "registry"
CITY_MODEL = OUT / "city-model-melbourne.json"
DWELLINGS = CACHE / "melbourne-clue-dwellings" / "municipality"
BLOCKS = CACHE / "melbourne-clue-blocks" / "all"
RENTS = CACHE / "homes-victoria-rents"
CENSUS = CACHE / "abs-census-2021-sa1-vic"
SA1 = CACHE / "abs-asgs-2021-sa1" / "city-of-melbourne"
DEVELOPMENTS = CACHE / "melbourne-development-activity" / "custom-1a4dca7e"
WORKBOOK = "moving-annual-rent-suburb-september-quarter-2025.xlsx"
ORIGIN = {"lat": -37.8085, "lon": 144.9525}  # the city model's, when it is not built

YEARS = list(range(2002, 2025))
TYPES = ["House/Townhouse", "Residential Apartments", "Student Apartments"]

# CLUE's small areas and the report's suburb groups that hold them (the groups reach beyond the municipality:
# Flemington, most of South Yarra and of Port Melbourne lie outside it).
GROUP_OF = {
    "Melbourne (CBD)": "CBD-St Kilda Rd",
    "Melbourne (Remainder)": "CBD-St Kilda Rd",
    "Southbank": "Southbank",
    "Docklands": "Docklands",
    "Carlton": "Carlton-Parkville",
    "Parkville": "Carlton-Parkville",
    "North Melbourne": "North Melbourne-West Melbourne",
    "West Melbourne (Residential)": "North Melbourne-West Melbourne",
    "West Melbourne (Industrial)": "North Melbourne-West Melbourne",
    "East Melbourne": "East Melbourne",
    "Kensington": "Flemington-Kensington",
    "South Yarra": "South Yarra",
    "Port Melbourne": "Port Melbourne",
}
SHEETS = ["1 bedroom flat", "2 bedroom flat", "3 bedroom flat", "2 bedroom house", "3 bedroom house", "4 bedroom house", "All properties"]

G33_BANDS = [("Negative_Nil_income_Tot", -1, 0), ("HI_1_149_Tot", 1, 149), ("HI_150_299_Tot", 150, 299), ("HI_300_399_Tot", 300, 399),
             ("HI_400_499_Tot", 400, 499), ("HI_500_649_Tot", 500, 649), ("HI_650_799_Tot", 650, 799), ("HI_800_999_Tot", 800, 999),
             ("HI_1000_1249_Tot", 1000, 1249), ("HI_1250_1499_Tot", 1250, 1499), ("HI_1500_1749_Tot", 1500, 1749),
             ("HI_1750_1999_Tot", 1750, 1999), ("HI_2000_2499_Tot", 2000, 2499), ("HI_2500_2999_Tot", 2500, 2999),
             ("HI_3000_3499_Tot", 3000, 3499), ("HI_3500_3999_Tot", 3500, 3999), ("HI_4000_more_Tot", 4000, None)]
G40_BANDS = [("R_1_74_Tot", 1, 74), ("R_75_99_Tot", 75, 99), ("R_100_149_Tot", 100, 149), ("R_150_199_Tot", 150, 199),
             ("R_200_224_Tot", 200, 224), ("R_225_274_Tot", 225, 274), ("R_275_349_Tot", 275, 349), ("R_350_449_Tot", 350, 449),
             ("R_450_549_Tot", 450, 549), ("R_550_649_Tot", 550, 649), ("R_650_749_Tot", 650, 749), ("R_750_849_Tot", 750, 849),
             ("R_850_949_Tot", 850, 949), ("R_950_over_Tot", 950, None)]
G35_SIZES = ["Num_Psns_UR_1_Total", "Num_Psns_UR_2_Total", "Num_Psns_UR_3_Total", "Num_Psns_UR_4_Total", "Num_Psns_UR_5_Total", "Num_Psns_UR_6mo_Total"]
G36_COLUMNS = ["OPDs_Separate_house_Dwellings", "OPDs_SD_r_t_h_th_Tot_Dwgs", "OPDs_Flt_apart_Tot_Dwgs", "OPDs_Other_dwelling_Tot_Dwgs",
               "OPDs_Dwlling_structur_NS_Dwgs", "OPDs_Tot_OPDs_Dwellings", "Unoccupied_PDs_Dwgs", "Total_PDs_Dwellings"]
G37_COLUMNS = ["O_OR_DS_Sep_house", "O_OR_DS_SemiD_ro_or_tce_h_th", "O_OR_DS_Flat_apart", "O_OR_Total",
               "O_MTG_DS_Sep_house", "O_MTG_DS_SemiD_ro_or_tce_h_th", "O_MTG_DS_Flat_apart", "O_MTG_Total",
               "R_Tot_DS_Sep_house", "R_Tot_DS_SemiD_ro_or_tce_h_th", "R_Tot_DS_Flat_apart", "R_Tot_Total",
               "R_RE_Agt_Total", "R_ST_h_auth_Total", "R_Com_Hp_Total", "R_Psn_not_in_s_hh_Total", "R_Ot_landld_typ_Total", "R_Landld_typ_NS_Total",
               "Oth_ten_type_Total", "Ten_type_NS_Total", "Total_Total"]
G02_COLUMNS = ["Median_rent_weekly", "Median_tot_hhd_inc_weekly", "Average_household_size", "Median_age_persons"]
G44_COLUMNS = ["Sme_Usl_ad_1_yr_ago_as_2021_P", "Dif_Usl_ad_1_yr_ago_Sme_SA2_P", "Dif_Us_ad_1_ago_Dif_SA2_Tot_P", "Difnt_Usl_add_1_yr_ago_OS_P",
               "Difnt_Usl_add_1_yr_ago_Ns_P", "Difnt_Usl_add_1_yr_ago_Tot_P", "N_stated_P", "Tot_P"]


# ---- reading ----------------------------------------------------------------------------------------------------------

def read_csv(path: Path, delimiter=";") -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def number(text):
    text = (text or "").strip()
    if not text or text == "-":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def manifest(ident: str) -> dict:
    return json.loads((REGISTRY / f"{ident}.json").read_text(encoding="utf-8"))


def source_entry(ident: str, folder: Path) -> dict:
    m = manifest(ident)
    return {"id": ident, "title": m["title"], "credit": m["credit"], "licence": m["licence"], "retrieved": retrieved(folder)}


NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def column_index(ref: str) -> int:
    n = 0
    for ch in ref:
        if not ch.isalpha():
            break
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_workbook(path: Path) -> dict[str, list[list]]:
    """Every sheet of an .xlsx workbook as rows of cell values (shared strings resolved; numbers as text)."""
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    book = ET.fromstring(z.read("xl/workbook.xml"))
    target = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    sheets = {}
    for s in book.find("m:sheets", NS):
        t = target[s.get("{%s}id" % NS["r"])].lstrip("/")
        t = t if t.startswith("xl/") else "xl/" + t
        rows = []
        for row in ET.fromstring(z.read(t)).iter("{%s}row" % NS["m"]):
            cells = {}
            for c in row.findall("m:c", NS):
                v, kind = c.find("m:v", NS), c.get("t")
                if kind == "s" and v is not None:
                    value = shared[int(v.text)]
                elif kind == "inlineStr":
                    value = "".join(x.text or "" for x in c.iter("{%s}t" % NS["m"]))
                else:
                    value = v.text if v is not None else None
                cells[column_index(c.get("r"))] = value
            rows.append([cells.get(i) for i in range(max(cells) + 1)] if cells else [])
        sheets[s.get("name")] = rows
    return sheets


def quarter_label(text: str) -> str:
    month, year = text.split()
    return f"{year}Q{ {'Mar': 1, 'Jun': 2, 'Sep': 3, 'Sept': 3, 'Dec': 4}[month[:4] if month.startswith('Sept') else month[:3]] }"


def area_of_ring(ring) -> float:
    a = 0.0
    for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
        a += x0 * y1 - x1 * y0
    return a / 2


def centroid_of_ring(ring) -> tuple[float, float]:
    a = cx = cy = 0.0
    for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
        w = x0 * y1 - x1 * y0
        a += w
        cx += (x0 + x1) * w
        cy += (y0 + y1) * w
    if abs(a) < 1e-9:
        return float(np.mean([p[0] for p in ring])), float(np.mean([p[1] for p in ring]))
    return cx / (3 * a), cy / (3 * a)


class Frame:
    def __init__(self, origin: dict):
        self.lat0, self.lon0 = origin["lat"], origin["lon"]

    def xy(self, lat, lon):
        e, n, _ = enu(np.asarray(lat, float), np.asarray(lon, float), self.lat0, self.lon0)
        return e, n

    def ring(self, coords) -> list[list[float]]:
        lon = np.array([c[0] for c in coords], float)
        lat = np.array([c[1] for c in coords], float)
        e, n = self.xy(lat, lon)
        pts = [[round(float(x), 1), round(float(y), 1)] for x, y in zip(np.atleast_1d(e), np.atleast_1d(n))]
        if len(pts) > 1 and pts[0] == pts[-1]:
            pts = pts[:-1]
        return pts


def outer_rings(geometry) -> list:
    if geometry["type"] == "Polygon":
        return [geometry["coordinates"][0]]
    if geometry["type"] == "MultiPolygon":
        return [p[0] for p in geometry["coordinates"]]
    return []


def table(columns, rows):
    return {"columns": columns, "rows": rows}


# ---- building ---------------------------------------------------------------------------------------------------------

def build() -> dict:
    origin = json.loads(CITY_MODEL.read_text(encoding="utf-8"))["origin"] if CITY_MODEL.exists() else ORIGIN
    frame = Frame(origin)
    notes = []

    # 1. Blocks.
    blocks_raw = json.loads((BLOCKS / "blocks.geojson").read_text(encoding="utf-8"))
    blocks, block_rings, block_point, block_area = [], {}, {}, {}
    for f in blocks_raw["features"]:
        p = f["properties"]
        bid = int(p["block_id"])
        rings = [frame.ring(r) for r in outer_rings(f["geometry"])]
        rings = [r for r in rings if len(r) >= 3]
        if not rings:
            continue
        big = max(rings, key=lambda r: abs(area_of_ring(r)))
        cx, cy = centroid_of_ring(big)
        block_rings[bid], block_point[bid], block_area[bid] = rings, (cx, cy), p.get("region_name") or p.get("area_name")
        blocks.append([bid, block_area[bid], round(cx, 1), round(cy, 1), round(sum(abs(area_of_ring(r)) for r in rings), 0), rings])
    blocks.sort(key=lambda r: r[0])

    # 2. Buildings and their dwellings by type and year.
    rows = read_csv(DWELLINGS / "dwellings.csv")
    stock = defaultdict(lambda: np.zeros((len(TYPES), len(YEARS)), int))
    points = defaultdict(list)
    blocks_of, areas_of, address_of, base_of = defaultdict(Counter), defaultdict(Counter), {}, {}
    dropped = Counter()
    for r in rows:
        y = int(r["census_year"][:4])
        t = r["dwelling_type"]
        if t not in TYPES or y not in YEARS:
            dropped[f"type or year ({t}, {y})"] += 1
            continue
        pid = r["property_id"].strip()
        n = int(float(r["dwelling_number"] or 0))
        stock[pid][TYPES.index(t), YEARS.index(y)] += n
        lat, lon = number(r["latitude"]), number(r["longitude"])
        if lat is not None and lon is not None:
            points[pid].append((lat, lon))
        if r["block_id"].strip():
            blocks_of[pid][int(float(r["block_id"]))] += 1
        areas_of[pid][r["clue_small_area"]] += 1
        address_of[pid] = r["building_address"]
        base_of[pid] = r["base_property_id"]
    buildings, placed = [], Counter()
    for pid in sorted(stock, key=lambda s: int(s) if s.isdigit() else 10**12):
        block = blocks_of[pid].most_common(1)[0][0] if blocks_of[pid] else None
        if points[pid]:
            lat = float(np.mean([p[0] for p in points[pid]]))
            lon = float(np.mean([p[1] for p in points[pid]]))
            e, n = frame.xy(lat, lon)
            x, y = float(e), float(n)
            placed["its own point"] += 1
        elif block in block_point:
            x, y = block_point[block]
            placed["its block's point"] += 1
        else:
            dropped["no point and no block"] += 1
            continue
        area = areas_of[pid].most_common(1)[0][0]
        if block not in block_area:
            # The block the point lies in.
            hit = next((b for b, rings in block_rings.items() if any(contains(rg, (x, y)) for rg in rings)), None)
            if hit is not None:
                block = hit
        buildings.append([int(pid) if pid.isdigit() else pid, int(base_of[pid]) if base_of[pid].isdigit() else base_of[pid], block, area,
                          round(x, 1), round(y, 1), [stock[pid][k].tolist() for k in range(len(TYPES))]])
    totals = {str(y): [int(sum(b[6][k][i] for b in buildings)) for k in range(len(TYPES))] for i, y in enumerate(YEARS)}

    # 3. Rents.
    book = read_workbook(RENTS / WORKBOOK)
    groups = sorted(set(GROUP_OF.values()))
    rents, quarters = {}, None
    for sheet in SHEETS:
        rs = book[sheet]
        header = rs[1]
        qs = [quarter_label(header[i]) for i in range(2, len(header), 2) if header[i]]
        quarters = quarters or qs
        if qs != quarters:
            raise SystemExit(f"housing: the sheets' quarters differ ({sheet})")
        region = None
        for r in rs[3:]:
            if not r or len(r) < 3:
                continue
            region = r[0] or region
            name = r[1]
            key = None
            if name in groups:
                key = name
            elif name == "Group Total":
                key = region  # each region's total: Inner Melbourne's, and the metropolis's others (the outside option)
            if key is None:
                continue
            counts = [number(r[2 + 2 * i]) if 2 + 2 * i < len(r) else None for i in range(len(qs))]
            medians = [number(r[3 + 2 * i]) if 3 + 2 * i < len(r) else None for i in range(len(qs))]
            rents.setdefault(key, {})[sheet] = {"count": counts, "median": medians}
    missing = [g for g in groups if g not in rents]
    if missing:
        notes.append(f"suburb groups not found in the workbook: {missing}")

    # 4. Census.
    sa1_raw = json.loads((SA1 / "sa1.geojson").read_text(encoding="utf-8"))
    tables = {}
    for code in ["G02", "G33", "G35", "G36", "G37", "G40", "G44"]:
        tables[code] = {r["SA1_CODE_2021"]: r for r in read_csv(CENSUS / f"2021Census_{code}_VIC_SA1.csv", ",")}
    census = []
    for f in sa1_raw["features"]:
        p = f["properties"]
        code = p["sa1_code_2021"]
        rings = [frame.ring(r) for r in outer_rings(f["geometry"])]
        rings = [r for r in rings if len(r) >= 3]
        if not rings:
            continue
        big = max(rings, key=lambda r: abs(area_of_ring(r)))
        cx, cy = centroid_of_ring(big)
        block = next((b for b, brs in block_rings.items() if any(contains(rg, (cx, cy)) for rg in brs)), None)
        if block is None:
            continue
        row = [code, p["sa2_name_2021"], block, block_area[block], round(cx, 1), round(cy, 1), rings]
        for t, cols in [("G02", G02_COLUMNS), ("G33", [c for c, _, _ in G33_BANDS] + ["Partial_income_stated_Tot", "All_incomes_not_stated_Tot", "Tot_Tot"]),
                        ("G35", G35_SIZES + ["Total_Total"]), ("G36", G36_COLUMNS), ("G37", G37_COLUMNS),
                        ("G40", [c for c, _, _ in G40_BANDS] + ["Rent_ns_Tot", "Tot_Tot"]), ("G44", G44_COLUMNS)]:
            rec = tables[t].get(code)
            row.append([number(rec.get(c)) if rec else None for c in cols])
        census.append(row)
    census.sort(key=lambda r: r[0])

    # 5. Developments.
    devs = []
    for r in read_csv(DEVELOPMENTS / "developments.csv"):
        lat, lon = number(r["latitude"]), number(r["longitude"])
        if lat is None or lon is None:
            dropped["development without a point"] += 1
            continue
        e, n = frame.xy(lat, lon)
        ids = [r[k].strip() for k in ("property_id", "property_id_2", "property_id_3", "property_id_4", "property_id_5") if (r.get(k) or "").strip()]
        devs.append([r["status"], int(r["year_completed"]) if (r["year_completed"] or "").strip().isdigit() else None, r["clue_small_area"],
                     int(float(r["clue_block"])) if number(r["clue_block"]) is not None else None, ids, round(float(e), 1), round(float(n), 1),
                     number(r["floors_above"]), number(r["resi_dwellings"]), number(r["studio_dwe"]), number(r["one_bdrm_dwe"]),
                     number(r["two_bdrm_dwe"]), number(r["three_bdrm_dwe"]), number(r["student_apartments"]), number(r["student_beds"])])

    product = {
        "title": "The City of Melbourne's dwellings by building and type 2002-2024, its blocks, the rents of new lettings by suburb group, the 2021 Census's households at SA1, and its developments",
        "frame": {"origin": origin, "units": "metres", "axes": "x east, y north (local ENU about the origin, as the city model's)"},
        "years": YEARS,
        "types": TYPES,
        "blocks": table(["block", "small_area", "x", "y", "area_m2", "rings"], blocks),
        "buildings": table(["property", "base_property", "block", "small_area", "x", "y", "dwellings_by_type_by_year"], buildings),
        "totals": totals,
        "rents": {"quarters": quarters, "sheets": SHEETS, "groups": rents, "group_of_small_area": GROUP_OF,
                  "meaning": "moving annual median weekly rent (A$) of new lettings over the year ending in the quarter, and their count; null where the report publishes none"},
        "census": {
            "date": "2021-08-10",
            "columns": ["sa1", "sa2_name", "block", "small_area", "x", "y", "rings", "G02", "G33", "G35", "G36", "G37", "G40", "G44"],
            "G02": G02_COLUMNS,
            "G33": [c for c, _, _ in G33_BANDS] + ["Partial_income_stated_Tot", "All_incomes_not_stated_Tot", "Tot_Tot"],
            "G33_bands": [[lo, hi] for _, lo, hi in G33_BANDS],
            "G35": G35_SIZES + ["Total_Total"],
            "G36": G36_COLUMNS,
            "G37": G37_COLUMNS,
            "G40": [c for c, _, _ in G40_BANDS] + ["Rent_ns_Tot", "Tot_Tot"],
            "G40_bands": [[lo, hi] for _, lo, hi in G40_BANDS],
            "G44": G44_COLUMNS,
            "rows": census,
            "note": "SA1s whose point lies in a CLUE block; the ABS perturbs small cells, so that tables' sums differ slightly",
        },
        "developments": table(["status", "year_completed", "small_area", "block", "properties", "x", "y", "floors", "dwellings", "studio", "one_bedroom",
                               "two_bedroom", "three_bedroom", "student_apartments", "student_beds"], devs),
        "sources": [source_entry("melbourne-clue-dwellings", DWELLINGS), source_entry("melbourne-clue-blocks", BLOCKS),
                    source_entry("homes-victoria-rents", RENTS), source_entry("abs-census-2021-sa1-vic", CENSUS),
                    source_entry("abs-asgs-2021-sa1", SA1), source_entry("melbourne-development-activity", DEVELOPMENTS)],
        "credit": "CLUE and the Development Activity Monitor: City of Melbourne (CC BY 4.0); Rental Report: Homes Victoria (CC BY 4.0); "
                  "Source: ABS, 2021 Census and ASGS Edition 3 (CC BY 4.0), ABS data used with permission from the Australian Bureau of Statistics",
        "report": {"rows": len(rows), "buildings": len(buildings), "placed": dict(placed), "dropped": dict(dropped), "blocks": len(blocks),
                   "sa1": len(census), "developments": Counter(d[0] for d in devs), "quarters": [quarters[0], quarters[-1]], "notes": notes},
    }
    path = OUT / "housing-melbourne.json"
    write(path, product)
    register("housing", {"name": "housing-melbourne", "file": path.name, "title": product["title"], "domain": "urban-science",
                         "years": [YEARS[0], YEARS[-1]], "buildings": len(buildings), "sa1": len(census),
                         "source": "City of Melbourne (CLUE, Development Activity Monitor); Homes Victoria; Australian Bureau of Statistics"})
    print(f"  housing-melbourne: {len(buildings)} buildings ({dict(placed)}), {len(blocks)} blocks, {len(census)} SA1s, "
          f"{len(devs)} developments, rents {quarters[0]}-{quarters[-1]} for {len(rents)} groups; {path.stat().st_size / 1e6:.2f} MB")
    print(f"  dwellings 2002: {totals['2002']}, 2024: {totals['2024']} (house/townhouse, apartments, student apartments); dropped {dict(dropped)}")
    return product


NEEDS = [DWELLINGS / "dwellings.csv", BLOCKS / "blocks.geojson", RENTS / WORKBOOK, CENSUS / "2021Census_G33_VIC_SA1.csv", SA1 / "sa1.geojson",
         DEVELOPMENTS / "developments.csv"]


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace")
    lacking = [p for p in NEEDS if not p.exists()]
    if lacking:
        print("housing: missing " + ", ".join(str(p.relative_to(DATA.parent)) for p in lacking) + "; run: python data/acquisition/fetch.py "
              "melbourne-clue-dwellings melbourne-clue-blocks homes-victoria-rents abs-census-2021-sa1-vic abs-asgs-2021-sa1, and "
              "melbourne-development-activity --set south=-37.86 --set west=144.88 --set north=-37.77 --set east=145.0", file=sys.stderr)
        return 1
    build()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
