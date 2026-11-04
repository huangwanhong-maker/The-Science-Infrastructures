"""Build the anatomy product of a woman's whole body from the Human Reference Atlas (HuBMAP), for the anatomy
domain and Fabrica.

The Human Reference Atlas's 3D Reference Organ Set for United, Female (v1.8) places every organ modelled from the
Visible Human Female in one frame: a binary glTF file whose scene nodes are the anatomical structures, grouped
in system layers (VH_F_integumentary_system, VH_F_circulatory_system, ...), and a graph giving each structure's
UBERON identifier and label. Its coordinates are metres with y up, z anterior and +x the subject's left (the
heart and spleen at +x, the liver at −x), one node (the skin) carrying its own transform.

What the product does, and nothing else (the science is the domain's):

- **The frame.** Each mesh is carried by its nodes' transforms into the model's frame, then into the LPS patient
  frame of DICOM in millimetres, as BodyParts3D's body is: x_LPS = 1000 x, y_LPS = −1000 z, z_LPS = 1000 (y −
  y_soles), a proper rotation (triangles keep their outward orientation), with the soles at z = 0.
- **The parts.** Every node with a mesh is an element; every node is a concept, named by its structure's label
  in the graph (its side added from the node's _L or _R where the label lacks it), so that the node tree is
  the part-of tree over the concepts, a tree as the domain needs. The UBERON identifier of each concept is
  kept in its own table: the Atlas names structures by UBERON, not by FMA, and no is-a tree is given.
- **The meshes quantised** exactly as BodyParts3D's product quantises them (positions uint16 over each
  element's box, normals int8, indices uint16 or uint32, blocks aligned to 4 bytes, four packs).

Products (data/cache/products/), in the layout of body-bp3d (body_anatomy.py), named body-hra-f:

- body-hra-f.json: the catalogue the page reads (schema body/1), with the sex, the source and each element's UBERON.
- body-hra-f-elements.tsv, -concepts.tsv, -isa.tsv (empty), -partof.tsv, -compounds-isa.tsv (empty),
  -compounds-partof.tsv (each concept's elements, the tree's closure), -ontology.tsv (concept, UBERON, label).
- body-hra-f-0.bin ... -3.bin and their .bin.gz.

Usage::

    python data/acquisition/fetch.py hra-female
    python data/products/body_female.py
"""

from __future__ import annotations

import datetime as dt
import gzip
import json
import struct
import sys
import time
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

CACHE = DATA / "cache" / "hra-female"
NAME = "body-hra-f"
PACKS = 4
COMPONENTS = {5120: np.int8, 5121: np.uint8, 5122: np.int16, 5123: np.uint16, 5125: np.uint32, 5126: np.float32}
WIDTH = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4}


def read_glb(path: Path) -> tuple[dict, bytes]:
    b = path.read_bytes()
    magic, version, length = struct.unpack_from("<4sII", b, 0)
    if magic != b"glTF" or version != 2:
        raise SystemExit(f"{path}: not a glTF 2.0 binary")
    at, chunks = 12, {}
    while at < length:
        n, kind = struct.unpack_from("<I4s", b, at)
        chunks[kind] = b[at + 8:at + 8 + n]
        at += 8 + n
    return json.loads(chunks[b"JSON"]), chunks[b"BIN\x00"]


def accessor(g: dict, binary: bytes, k: int) -> np.ndarray:
    a = g["accessors"][k]
    bv = g["bufferViews"][a["bufferView"]]
    dtype = np.dtype(COMPONENTS[a["componentType"]])
    width = WIDTH[a["type"]]
    off = bv.get("byteOffset", 0) + a.get("byteOffset", 0)
    stride = bv.get("byteStride", 0) or dtype.itemsize * width
    rows = np.frombuffer(binary, np.uint8, a["count"] * stride, off).reshape(-1, stride)[:, :dtype.itemsize * width]
    out = rows.copy().view(dtype).reshape(-1, width).astype(np.float64 if dtype.kind == "f" else np.int64)
    if a.get("normalized") and dtype.kind != "f":
        out = out / np.iinfo(dtype).max
    return out if width > 1 else out.reshape(-1)


def local_matrix(node: dict) -> np.ndarray:
    M = np.eye(4)
    if "matrix" in node:
        return np.array(node["matrix"], dtype=float).reshape(4, 4).T
    if "rotation" in node:
        x, y, z, w = node["rotation"]
        M[:3, :3] = [[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]]
    if "scale" in node:
        M[:3, :3] = M[:3, :3] @ np.diag(node["scale"])
    if "translation" in node:
        M[:3, 3] = node["translation"]
    return M


def humane(node_name: str) -> str:
    s = node_name
    for p in ("VH_F_", "Allen_", "Yao_"):
        s = s.removeprefix(p)
    return s.replace("_", " ").strip()


def label_for(node_name: str, entity: dict | None) -> str:
    label = (entity or {}).get("pref_label") or humane(node_name)
    label = label[0].lower() + label[1:] if label[:1].isupper() and not label[:2].isupper() else label
    side = "left" if node_name.endswith(("_L", "_l")) else "right" if node_name.endswith(("_R", "_r")) else ""
    if side and side not in label.lower():
        label = f"{side} {label}"
    return label


def build() -> Path:
    t0 = time.time()
    g, binary = read_glb(CACHE / "3d-vh-f-united.glb")
    graph = json.loads((CACHE / "graph.json").read_text(encoding="utf-8"))
    entity = {e["object_reference"]["file_subpath"]: e for e in graph["data"] if "object_reference" in e}
    nodes = g["nodes"]
    root = g["scenes"][g.get("scene", 0)]["nodes"][0]

    # Node names made unique concept ids; the tree, depth first; world transforms.
    seen: dict[str, int] = {}
    cid: dict[int, str] = {}
    world: dict[int, np.ndarray] = {}
    parent: dict[int, int] = {}
    order: list[int] = []
    stack = [(root, np.eye(4), -1)]
    while stack:
        i, M, p = stack.pop()
        nm = nodes[i].get("name") or f"node{i}"
        seen[nm] = seen.get(nm, 0) + 1
        cid[i] = nm if seen[nm] == 1 else f"{nm}#{seen[nm]}"
        world[i] = M @ local_matrix(nodes[i])
        if p >= 0:
            parent[i] = p
        order.append(i)
        for c in reversed(nodes[i].get("children", [])):
            stack.append((c, world[i], i))

    # The meshes in the model's frame (metres), then the soles.
    raw = []
    for i in order:
        nd = nodes[i]
        if "mesh" not in nd:
            continue
        Ps, Ns, Ts, base = [], [], [], 0
        R = world[i][:3, :3]
        for prim in g["meshes"][nd["mesh"]]["primitives"]:
            if prim.get("mode", 4) != 4 or "indices" not in prim:
                continue
            P = accessor(g, binary, prim["attributes"]["POSITION"])
            P = P @ R.T + world[i][:3, 3]
            N = accessor(g, binary, prim["attributes"]["NORMAL"]) @ np.linalg.inv(R).T if "NORMAL" in prim["attributes"] else None
            T = accessor(g, binary, prim["indices"]).reshape(-1, 3) + base
            Ps.append(P), Ns.append(N), Ts.append(T)
            base += len(P)
        if not Ps:
            continue
        N = None if any(n is None for n in Ns) else np.vstack(Ns)
        raw.append((i, np.vstack(Ps), N, np.vstack(Ts)))
    y_soles = min(P[:, 1].min() for _, P, _, _ in raw)

    def to_lps(P: np.ndarray) -> np.ndarray:
        return np.stack([1000 * P[:, 0], -1000 * P[:, 2], 1000 * (P[:, 1] - y_soles)], axis=1)

    def normals_lps(N: np.ndarray) -> np.ndarray:
        return np.stack([N[:, 0], -N[:, 2], N[:, 1]], axis=1)

    elements, blocks = [], []
    for i, P, N, T in raw:
        P = to_lps(P)
        if N is None:
            N = np.zeros_like(P)
            a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
            fn = np.cross(b - a, c - a)
            for k in range(3):
                np.add.at(N, T[:, k], fn)
        else:
            N = normals_lps(N)
        N = N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30)
        lo, hi = P.min(axis=0), P.max(axis=0)
        ext = np.where(hi - lo > 0, hi - lo, 1e-3)
        q = np.rint((P - lo) / ext * 65535).astype(np.uint16)
        nq = np.zeros((len(P), 4), dtype=np.int8)
        nq[:, :3] = np.clip(np.rint(N * 127), -127, 127).astype(np.int8)
        wide = len(P) > 65535
        idx = T.astype(np.uint32 if wide else np.uint16).reshape(-1)
        ent = entity.get(nodes[i].get("name", ""))
        elements.append({
            "id": cid[i], "concept": cid[i], "name": label_for(nodes[i].get("name", cid[i]), ent), "conceptFrom": "node",
            "uberon": (ent or {}).get("representation_of"),
            "vertices": int(len(P)), "triangles": int(len(T)), "indexBytes": 4 if wide else 2,
            "box": [float(x) for x in lo] + [float(x) for x in hi],
            "publishedVolume": None, "publishedBounds": None, "duplicateOf": None,
        })
        blocks.append((q.tobytes(), nq.tobytes(), idx.tobytes()))
    print(f"read {len(elements)} meshes, {sum(e['vertices'] for e in elements)} vertices, {sum(e['triangles'] for e in elements)} triangles ({time.time() - t0:.1f} s)")

    # Concepts, the part-of tree and each concept's elements (the tree's closure).
    concept_name = {cid[i]: label_for(nodes[i].get("name", cid[i]), entity.get(nodes[i].get("name", ""))) for i in order}
    uberon = {cid[i]: (entity.get(nodes[i].get("name", "")) or {}).get("representation_of", "") for i in order}
    partof = [[cid[parent[i]], cid[i]] for i in order if i in parent]
    element_of = {e["id"]: e["id"] for e in elements}
    below: dict[str, list[str]] = {cid[i]: [] for i in order}
    for i in order:
        if cid[i] in element_of:
            j = i
            while True:
                below[cid[j]].append(cid[i])
                if j not in parent:
                    break
                j = parent[j]
    compounds = [[c, e] for c in below for e in below[c] if not (len(below[c]) == 1 and below[c][0] == c)]

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
                e["vertices"], e["triangles"], e["indexBytes"], *(repr(x) for x in e["box"]), None, None]
               for i, e in enumerate(elements)])
    write_tsv("concepts", ["concept", "name", "in_isa", "in_partof"], [[c, concept_name[c], 0, 1] for c in sorted(concept_name)])
    write_tsv("isa", ["parent", "child"], [])
    write_tsv("partof", ["parent", "child"], partof)
    write_tsv("compounds-isa", ["concept", "element"], [])
    write_tsv("compounds-partof", ["concept", "element"], compounds)
    write_tsv("ontology", ["concept", "uberon", "label"], [[c, uberon[c], concept_name[c]] for c in sorted(concept_name)])

    meta = {
        "name": NAME, "kind": "body", "schema": "body/1", "sex": "female", "source": "hra",
        "title": "The whole body of an adult woman from the Human Reference Atlas (HuBMAP), 3D Reference Organ Set for United, Female v1.8: the Visible Human Female's organs in one frame",
        "frame": {"name": "LPS", "unit": "mm", "about": "DICOM's patient frame: +x the subject's left, +y posterior, +z superior; the soles at z = 0"},
        "packs": [{"file": f"{NAME}-{k}.bin.gz", "raw": f"{NAME}-{k}.bin", "bytes": len(buf)} for k, buf in enumerate(packs)],
        "tables": {t: f"{NAME}-{t}.tsv" for t in ("elements", "concepts", "isa", "partof", "compounds-isa", "compounds-partof", "ontology")},
        "quantisation": "positions uint16 over each element's box: p = lo + q (hi − lo) / 65535; normals int8 (x, y, z, 0) / 127; indices uint16 or uint32 by indexBytes; blocks aligned to 4 bytes",
        "elements": elements,
        "counts": {"elements": len(elements), "vertices": int(sum(e["vertices"] for e in elements)), "triangles": int(sum(e["triangles"] for e in elements)),
                   "concepts": len(concept_name), "copies": 0},
        "credit": "Female body: the Human Reference Atlas, 3D Reference Organ Set for United, Female v1.8 (HuBMAP), CC BY 4.0",
        "sources": [{"dataset": "hra-female", "files": ["3d-vh-f-united.glb", "graph.json"], "version": graph["metadata"].get("version"),
                     "licence": graph["metadata"].get("license")}],
        "built_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }
    out = OUT / f"{NAME}.json"
    out.write_text(json.dumps(meta, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("body", {"name": NAME, "file": f"{NAME}.json", "title": meta["title"], "frame": "LPS", "sex": "female"})
    sizes = ", ".join(f"{(OUT / f'{NAME}-{k}.bin.gz').stat().st_size / 1e6:.1f} MB" for k in range(PACKS))
    print(f"wrote {out.name}, the tables and the packs (gzipped {sizes}) in {time.time() - t0:.1f} s")
    return out


if __name__ == "__main__":
    build()
