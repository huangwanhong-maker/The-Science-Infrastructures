"""Build the anatomy product of the whole human body from BodyParts3D 4.0, for the anatomy domain and the
Body Atlas.

BodyParts3D (Mitsuhashi et al. 2009) gives one adult man as 2,234 meshes ("elements"), each the
representation of a concept of the Foundational Model of Anatomy (FMA), with FMA's is-a and part-of trees
over the concepts and, for each compound concept, the elements that compose it. Its coordinates are
millimetres in the LPS patient frame of DICOM (+x the subject's left, +y posterior, +z superior; the
database's coordinate_system.png), which the product keeps.

What the product adds to the files, and nothing else (the science is the domain's):

- Each element's concept. Seventeen headers name none; such an element takes the most specific compound
  concept of the is-a tree that contains it (the one composed of the fewest elements), and says so.
- Duplicate copies. Some concepts are drawn twice by meshes of the same shape and place, some remeshed (the
  pelvic floor's muscles, the external anal sphincter three times, the hyoid bone). An element is a copy of
  an earlier one when it has the same concept, its bounds within 3 mm and its published volume within 3 %;
  it is kept, marked `duplicate_of`, so that a drawing shows one. Looser bounds would take the two halves of
  a midline set (the interspinales on either side, each called "set of interspinales lumborum") for copies.
- The meshes quantised: positions as unsigned 16-bit integers over the element's own bounding box (a step
  of at most 1/65535 of its extent: 0.03 mm for the skin, finer than the files' 0.01 mm digits for parts
  under 650 mm), normals as signed 8-bit integers, triangles as 16- or 32-bit indices, each block aligned to
  4 bytes; packed into four files of about equal size, in element order.

Products (data/cache/products/):

- body-bp3d.json: the catalogue the page reads (schema body/1): the elements with their concept, name,
  pack, byte offsets, counts, quantisation box, published volume and bounds, and copies; the packs; the
  sources.
- body-bp3d-elements.tsv, -concepts.tsv, -isa.tsv, -partof.tsv, -compounds-isa.tsv, -compounds-partof.tsv:
  the same as tab-separated tables, which the domain's C++ reads without a JSON parser.
- body-bp3d-0.bin ... -3.bin and their .bin.gz: the packs, raw for the engine and gzipped for the page.

Usage::

    python data/acquisition/fetch.py bodyparts3d
    python data/products/body_anatomy.py
"""

from __future__ import annotations

import collections
import datetime as dt
import gzip
import json
import re
import sys
import time
import zipfile
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

CACHE = DATA / "cache"
DATASET = "bodyparts3d"
ZIP = "isa_BP3D_4.0_obj_99.zip"
NAME = "body-bp3d"
PACKS = 4
HEADER = re.compile(r"^# ([^:\n]+?) ?: (.*)$", re.M)


def table(name: str) -> list[list[str]]:
    rows = (CACHE / DATASET / name).read_text(encoding="utf-8").splitlines()
    return [r.split("\t") for r in rows[1:] if r.strip()]


def read_obj(text: str) -> tuple[dict, np.ndarray, np.ndarray | None, np.ndarray]:
    """The header, positions (n × 3), normals (n × 3, or None) and triangles (m × 3, 0-based) of an OBJ file."""
    head = dict(HEADER.findall(text))
    v, vn, f = [], [], []
    for line in text.splitlines():
        if line.startswith("v "):
            v.append(line[2:])
        elif line.startswith("vn "):
            vn.append(line[3:])
        elif line.startswith("f "):
            f.append(line[2:])
    P = np.array(" ".join(v).split(), dtype=np.float64).reshape(-1, 3)
    N = np.array(" ".join(vn).split(), dtype=np.float64).reshape(-1, 3) if vn else None
    corners = " ".join(f).replace("//", "/").split()
    pairs = [c.split("/") for c in corners]
    vi = np.array([int(p[0]) for p in pairs], dtype=np.int64).reshape(-1, 3) - 1
    if N is not None:
        ni = np.array([int(p[-1]) for p in pairs], dtype=np.int64).reshape(-1, 3) - 1
        if len(N) != len(P) or not np.array_equal(vi, ni):
            N = None  # normals indexed apart from positions: computed by the reader instead
    return head, P, N, vi


def bounds_of(text: str) -> list[float] | None:
    m = re.match(r"\(([^)]*)\)-\(([^)]*)\)", text or "")
    if not m:
        return None
    return [float(x) for x in m.group(1).split(",")] + [float(x) for x in m.group(2).split(",")]


def provenance() -> list:
    r = CACHE / DATASET / "retrieval.json"
    rec = json.loads(r.read_text(encoding="utf-8")) if r.exists() else {}
    return [{"dataset": DATASET, "file": x["file"], "sha256": x.get("sha256"), "retrieved_at": rec.get("retrieved_at")} for x in rec.get("files", [])]


def build() -> Path:
    t0 = time.time()
    path = CACHE / DATASET / ZIP
    if not path.exists():
        raise SystemExit(f"{path} is missing; run: python data/acquisition/fetch.py {DATASET}")

    # The trees and the compounds.
    concept_name: dict[str, str] = {}
    in_tree = collections.defaultdict(set)
    for tree in ("isa", "partof"):
        for cid, _rep, name in table(f"{tree}_parts_list_e.txt"):
            concept_name.setdefault(cid, name)
            in_tree[cid].add(tree)
    relations = {}
    for tree in ("isa", "partof"):
        rel = []
        for p, pn, c, cn in table(f"{tree}_inclusion_relation_list.txt"):
            concept_name.setdefault(p, pn)
            concept_name.setdefault(c, cn)
            in_tree[p].add(tree)
            in_tree[c].add(tree)
            rel.append((p, c))
        relations[tree] = rel
    compounds = {tree: [(c, f) for c, _n, f in table(f"{tree}_element_parts.txt")] for tree in ("isa", "partof")}
    size_of = collections.Counter(c for c, _f in compounds["isa"])
    containing = collections.defaultdict(list)
    for c, f in compounds["isa"]:
        containing[f].append(c)

    # The meshes.
    z = zipfile.ZipFile(path)
    names = sorted(n for n in z.namelist() if n.endswith(".obj"))
    elements, blocks = [], []
    for n in names:
        fid = Path(n).stem
        head, P, N, T = read_obj(z.read(n).decode("utf-8", "replace"))
        concept = head.get("Concept ID", "").strip()
        source = "header"
        if not concept:
            options = sorted(containing.get(fid, []), key=lambda c: (size_of[c], c))
            if not options:
                raise SystemExit(f"{fid}: no concept in its header and no compound containing it")
            concept, source = options[0], "compound"
        name = head.get("English name", "").strip() or concept_name.get(concept, concept)
        lo, hi = P.min(axis=0), P.max(axis=0)
        ext = np.where(hi - lo > 0, hi - lo, 1e-3)
        q = np.rint((P - lo) / ext * 65535).astype(np.uint16)
        if N is None:
            N = np.zeros_like(P)
            a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
            fn = np.cross(b - a, c - a)
            for k in range(3):
                np.add.at(N, T[:, k], fn)
        N = N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30)
        nq = np.zeros((len(P), 4), dtype=np.int8)
        nq[:, :3] = np.clip(np.rint(N * 127), -127, 127).astype(np.int8)
        wide = len(P) > 65535
        idx = T.astype(np.uint32 if wide else np.uint16).reshape(-1)
        published = head.get("Volume(cm3)")
        elements.append({
            "id": fid, "concept": concept, "name": name, "conceptFrom": source,
            "vertices": int(len(P)), "triangles": int(len(T)), "indexBytes": 4 if wide else 2,
            "box": [float(x) for x in lo] + [float(x) for x in hi],
            "publishedVolume": float(published) if published else None,
            "publishedBounds": bounds_of(head.get("Bounds(mm)")),
        })
        blocks.append((q.tobytes(), nq.tobytes(), idx.tobytes()))
    print(f"read {len(elements)} meshes, {sum(e['vertices'] for e in elements)} vertices, {sum(e['triangles'] for e in elements)} triangles ({time.time() - t0:.1f} s)")

    # Copies: same concept, bounds within 3 mm and published volume within 3 %.
    seen: dict[str, list[int]] = collections.defaultdict(list)
    copies = 0
    for i, e in enumerate(elements):
        e["duplicateOf"] = None
        for j in seen[e["concept"]]:
            o = elements[j]
            same_box = max(abs(a - b) for a, b in zip(e["box"], o["box"])) <= 3.0
            va, vb = e["publishedVolume"], o["publishedVolume"]
            same_volume = va is not None and vb is not None and abs(va - vb) <= 0.03 * max(va, vb)
            if same_box and same_volume:
                e["duplicateOf"] = o["id"]
                copies += 1
                break
        seen[e["concept"]].append(i)
    print(f"copies of earlier elements: {copies}")

    # Packs of about equal size, in element order.
    def pad(b: bytes) -> bytes:
        return b + b"\0" * (-len(b) % 4)
    total = sum(len(pad(a)) + len(pad(b)) + len(pad(c)) for a, b, c in blocks)
    target = total / PACKS
    packs = [bytearray() for _ in range(PACKS)]
    k = 0
    for e, (pos, nrm, idx) in zip(elements, blocks):
        if len(packs[k]) >= target and k < PACKS - 1:
            k += 1
        buf = packs[k]
        e["pack"] = k
        e["positionOffset"] = len(buf)
        buf += pad(pos)
        e["normalOffset"] = len(buf)
        buf += pad(nrm)
        e["indexOffset"] = len(buf)
        buf += pad(idx)

    OUT.mkdir(parents=True, exist_ok=True)
    for k, buf in enumerate(packs):
        (OUT / f"{NAME}-{k}.bin").write_bytes(bytes(buf))
        (OUT / f"{NAME}-{k}.bin.gz").write_bytes(gzip.compress(bytes(buf), compresslevel=6))

    def write_tsv(name: str, header: list[str], rows) -> None:
        lines = ["\t".join(header)] + ["\t".join("" if x is None else str(x) for x in r) for r in rows]
        (OUT / f"{NAME}-{name}.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    write_tsv("elements", ["index", "id", "concept", "name", "concept_from", "pack", "position_offset", "normal_offset", "index_offset",
                           "vertices", "triangles", "index_bytes", "box_lo_x", "box_lo_y", "box_lo_z", "box_hi_x", "box_hi_y", "box_hi_z",
                           "published_volume_cm3", "duplicate_of"],
              [[i, e["id"], e["concept"], e["name"], e["conceptFrom"], e["pack"], e["positionOffset"], e["normalOffset"], e["indexOffset"],
                e["vertices"], e["triangles"], e["indexBytes"], *(repr(x) for x in e["box"]), e["publishedVolume"], e["duplicateOf"]]
               for i, e in enumerate(elements)])
    write_tsv("concepts", ["concept", "name", "in_isa", "in_partof"],
              [[c, concept_name[c], int("isa" in in_tree[c]), int("partof" in in_tree[c])] for c in sorted(concept_name)])
    for tree in ("isa", "partof"):
        write_tsv(tree, ["parent", "child"], relations[tree])
        write_tsv(f"compounds-{tree}", ["concept", "element"], compounds[tree])

    meta = {
        "name": NAME, "kind": "body", "schema": "body/1",
        "title": "The whole body of an adult man from BodyParts3D 4.0: 2,234 meshes of its parts, named by the Foundational Model of Anatomy",
        "frame": {"name": "LPS", "unit": "mm", "about": "DICOM's patient frame: +x the subject's left, +y posterior, +z superior"},
        "packs": [{"file": f"{NAME}-{k}.bin.gz", "raw": f"{NAME}-{k}.bin", "bytes": len(buf)} for k, buf in enumerate(packs)],
        "tables": {t: f"{NAME}-{t}.tsv" for t in ("elements", "concepts", "isa", "partof", "compounds-isa", "compounds-partof")},
        "quantisation": "positions uint16 over each element's box: p = lo + q (hi − lo) / 65535; normals int8 (x, y, z, 0) / 127; indices uint16 or uint32 by indexBytes; blocks 4-aligned",
        "elements": elements,
        "counts": {"elements": len(elements), "vertices": int(sum(e["vertices"] for e in elements)), "triangles": int(sum(e["triangles"] for e in elements)),
                   "concepts": len(concept_name), "copies": copies},
        "credit": "BodyParts3D, © The Database Center for Life Science, CC BY 4.0",
        "sources": provenance(),
        "built_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }
    out = OUT / f"{NAME}.json"
    out.write_text(json.dumps(meta, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("body", {"name": NAME, "file": f"{NAME}.json", "title": meta["title"], "frame": "LPS"})
    sizes = ", ".join(f"{(OUT / f'{NAME}-{k}.bin.gz').stat().st_size / 1e6:.1f} MB" for k in range(PACKS))
    print(f"wrote {out.name}, the tables and the packs (gzipped {sizes}) in {time.time() - t0:.1f} s")
    return out


if __name__ == "__main__":
    build()
