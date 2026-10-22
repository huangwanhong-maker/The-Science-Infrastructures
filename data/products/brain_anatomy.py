"""Build the anatomy product of the human brain: smooth surfaces of its parts,
in the MNI 152 ICBM 2009c asymmetric space, for the Brain Atlas to draw.

Every source lies on the same 1 mm grid of the 2009c template (193 x 229 x 193
voxels from (-96, -132, -78) mm, RAS), so that their structures fit together:

- the tissue probability maps of the template (grey matter, white matter), for
  the pial and white surfaces of each hemisphere and the outer surface of each
  half of the cerebellum;
- FreeSurfer's segmentation of the template (aseg), for the deep nuclei, the
  ventricles, the brainstem, the cerebellar white matter, the corpus callosum
  and the optic chiasm, and to tell the hemispheres apart;
- MASSP 2.0 (Bazin et al. 2020, 2025), for the finer nuclei and tracts;
- the thalamic nuclei of Najdenovska et al. (2018);
- the deep white-matter bundles of Guevara et al. (2012), from EBRAINS;
- the cerebellar nuclei of Julich-Brain 3.1 (Amunts et al. 2020);
- the Harvard-Oxford cortical atlas, lateralised, for the gyrus beneath every
  vertex of the pial surfaces.

A surface is extracted by surface nets (Gibson 1998), the method of the shared
visualization layer's isosurface.js, from a continuous field: for the cortex
the tissue probabilities themselves, for a labelled structure its mask blurred
by a Gaussian and upsampled twice by trilinear interpolation, so that the
surfaces have no voxel steps. The level is 0.5, lowered for thin structures
until the enclosed volume reaches 85 % of the voxels'. The surfaces are then
relaxed by Taubin's (1995) smoothing, which does not shrink them, and each
vertex gets its normal and its openness, a volumetric obscurance (Loos and
Sloan 2010): one minus the mean of the blurred tissue sampled at 1.5 to 7 mm
along the normal, so that sulci and creases are darker when drawn.

The product also holds one label volume joining the sources (the finest
structure at each voxel: a gyrus in the cortex, a nucleus in the thalamus, a
subfield in the hippocampus), the volume, centroid and bounds of every
structure, and their contacts: for each pair of structures, the number of
voxel faces they share.

Products (data/cache/products/):

- anatomy-mni152.json: the catalogue of meshes and regions (schema anatomy/1);
- anatomy-mni152-core.bin.gz, -white.bin.gz and -bundles.bin.gz: the meshes in
  three packs (the white-matter surfaces and the fibre bundles apart, to be read
  only when they are shown), gzipped: per mesh, positions as int16 in
  hundredths of a millimetre, normals as int8 (x, y, z, 0), openness as uint8,
  the gyrus label of each vertex as uint8 (pial surfaces only), and triangles
  as uint16 or uint32 indices, and for the large surfaces a sulcal depth as uint8,
  each block aligned to 4 bytes;
- anatomy-mni152-labels.nii.gz: the joined label volume, uint16, with its table
  in the catalogue.

Usage::

    python data/acquisition/fetch.py templateflow-mni152 massp-2 harvard-oxford mial-thalamic-nuclei ebrains-deep-white-matter ebrains-julich-brain
    python data/products/brain_anatomy.py
"""

from __future__ import annotations

import datetime as dt
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
from nifti import read_nifti, write_nifti  # noqa: E402
from product_index import OUT, register  # noqa: E402

CACHE = DATA / "cache"
TF = "tpl-MNI152NLin2009cAsym_res-01_"
FILES = {
    "gm": ("templateflow-mni152", TF + "label-GM_probseg.nii.gz"),
    "wm": ("templateflow-mni152", TF + "label-WM_probseg.nii.gz"),
    "aseg": ("templateflow-mni152", TF + "seg-aseg_dseg.nii.gz"),
    "massp": ("massp-2", TF + "atlas-MASSP20_dseg.nii.gz"),
    "ho": ("harvard-oxford", TF + "atlas-HOCPAL_desc-th25_dseg.nii.gz"),
    "mial": ("mial-thalamic-nuclei", TF + "atlas-MIAL67ThalamicNuclei_dseg.nii.gz"),
    "dwm": ("ebrains-deep-white-matter", "long-bundles_maxprob.nii.gz"),
    "jlh": ("ebrains-julich-brain", "JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152.nii.gz"),
    "jrh": ("ebrains-julich-brain", "JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152.nii.gz"),
}

# FreeSurfer's aseg labels (FreeSurferColorLUT.txt), by the part of the brain they belong to.
CEREBRUM = {"lh": [2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 26, 28, 30, 31],
            "rh": [41, 42, 43, 44, 49, 50, 51, 52, 53, 54, 58, 60, 62, 63]}
CEREBELLUM = {"lh": [7, 8], "rh": [46, 47]}
BRAINSTEM = [16, 15]
MIDLINE = [14, 24, 72, 77, 85, 251, 252, 253, 254, 255]
FILL = {"lh": [2, 4, 5, 10, 11, 12, 13, 17, 18, 26, 28, 31], "rh": [41, 43, 44, 49, 50, 51, 52, 53, 54, 58, 60, 63]}
DEEP = {"lh": [4, 5, 10, 11, 12, 13, 26, 28, 31], "rh": [43, 44, 49, 50, 51, 52, 58, 60, 63]}
ASEG_MESHES = [4, 5, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 26, 28, 31,
               43, 44, 46, 49, 50, 51, 52, 53, 54, 58, 60, 63, 85, 251, 252, 253, 254, 255]
# MASSP 2.0 structures that the aseg already gives (caudate, thalamus, ventricles,
# amygdala, putamen, accumbens) are left to the aseg, which shares the cortex's segmentation.
MASSP_SKIP = {1, 2, 13, 14, 15, 16, 17, 18, 19, 20, 38, 39, 40, 41}
# Julich-Brain 3.1's cerebellar nuclei, by their labels in each hemisphere's maximum
# probability map (siibra-configurations, maps/mni152-jba31-labelled.json).
JULICH_NUCLEI = {"dentate": [56, 115], "interposed": [87], "fastigial": [187]}
# Structures of more voxels than this are meshed at 1 mm; smaller ones at 0.5 mm.
SMALL = 4000
# Blur (voxels) and upsampling for each source.
PARAMS = {"aseg": (0.8, 2), "massp": (0.65, 2), "mial": (0.8, 2), "dwm": (1.0, 2), "jlh": (0.6, 2), "jrh": (0.6, 2)}
# The blur and the smoothing of large structures, by source: the bundles are territories of
# a maximum-probability map, ragged at the voxel scale, and are drawn smoother than nuclei.
LARGE = {"dwm": (1.7, 14)}

# The 12 edges of a cell, its corners numbered by bits (i, j, k), as isosurface.js has them.
EDGES = [(0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7)]


# ---- fields -------------------------------------------------------------------------------

def blur(a: np.ndarray, sigma: float) -> np.ndarray:
    """A separable Gaussian blur of standard deviation sigma (voxels), zero outside."""
    if sigma <= 0:
        return a.astype(np.float32)
    r = max(1, int(np.ceil(3 * sigma)))
    x = np.arange(-r, r + 1, dtype=np.float64)
    w = np.exp(-0.5 * (x / sigma) ** 2)
    w = (w / w.sum()).astype(np.float32)
    out = a.astype(np.float32)
    for axis in range(3):
        src = np.moveaxis(out, axis, 0)
        pad = np.zeros((src.shape[0] + 2 * r,) + src.shape[1:], np.float32)
        pad[r:r + src.shape[0]] = src
        acc = np.zeros_like(src)
        for t in range(2 * r + 1):
            acc += w[t] * pad[t:t + src.shape[0]]
        out = np.moveaxis(acc, 0, axis)
    return np.ascontiguousarray(out)


def upsample2(a: np.ndarray) -> np.ndarray:
    """Trilinear upsampling by two: n points become 2n - 1, the old ones kept at even indices."""
    out = a.astype(np.float32)
    for axis in range(3):
        src = np.moveaxis(out, axis, 0)
        n = src.shape[0]
        new = np.empty((2 * n - 1,) + src.shape[1:], np.float32)
        new[0::2] = src
        new[1::2] = 0.5 * (src[:-1] + src[1:])
        out = np.moveaxis(new, 0, axis)
    return np.ascontiguousarray(out)


def sample(field: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Trilinear interpolation of a field at points p (n x 3, in its index coordinates), 0 outside."""
    nx, ny, nz = field.shape
    i0 = np.floor(p).astype(np.int64)
    f = (p - i0).astype(np.float32)
    out = np.zeros(len(p), np.float32)
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                i, j, k = i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz
                ok = (i >= 0) & (j >= 0) & (k >= 0) & (i < nx) & (j < ny) & (k < nz)
                w = (f[:, 0] if dx else 1 - f[:, 0]) * (f[:, 1] if dy else 1 - f[:, 1]) * (f[:, 2] if dz else 1 - f[:, 2])
                v = np.zeros(len(p), np.float32)
                v[ok] = field[i[ok], j[ok], k[ok]]
                out += w * v
    return out


# ---- surface nets -------------------------------------------------------------------------

def surface_nets(field: np.ndarray, level: float) -> tuple[np.ndarray, np.ndarray]:
    """Vertices (n x 3, index coordinates) and triangles (m x 3) of the surface at level,
    wound so that normals point toward lower values (outward for a mask)."""
    f = field.astype(np.float32) - level
    ins = f > 0
    nx, ny, nz = f.shape
    corner = lambda a, n: a[n & 1:nx - 1 + (n & 1), (n >> 1) & 1:ny - 1 + ((n >> 1) & 1), (n >> 2) & 1:nz - 1 + ((n >> 2) & 1)]  # noqa: E731
    mask = np.zeros((nx - 1, ny - 1, nz - 1), np.uint8)
    for n in range(8):
        mask |= corner(ins, n).astype(np.uint8) << n
    ci, cj, ck = np.nonzero((mask != 0) & (mask != 255))
    m = mask[ci, cj, ck].astype(np.int32)
    vals = [f[ci + (n & 1), cj + ((n >> 1) & 1), ck + ((n >> 2) & 1)] for n in range(8)]
    s = np.zeros((len(ci), 3), np.float32)
    cnt = np.zeros(len(ci), np.float32)
    for a, b in EDGES:
        cross = ((m >> a) & 1) != ((m >> b) & 1)
        t = np.zeros(len(ci), np.float32)
        t[cross] = vals[a][cross] / (vals[a][cross] - vals[b][cross])
        for axis in range(3):
            pa, pb = (a >> axis) & 1, (b >> axis) & 1
            s[:, axis] += np.where(cross, pa + (pb - pa) * t, 0)
        cnt += cross
    verts = np.stack([ci, cj, ck], 1).astype(np.float32) + s / cnt[:, None]
    vid = np.full((nx - 1, ny - 1, nz - 1), -1, np.int64)
    vid[ci, cj, ck] = np.arange(len(ci))
    tris = []

    def quads(a, b, c, e, flip):
        t1 = np.where(flip[:, None], np.stack([a, b, c], 1), np.stack([a, c, b], 1))
        t2 = np.where(flip[:, None], np.stack([a, c, e], 1), np.stack([a, e, c], 1))
        tris.append(t1)
        tris.append(t2)

    # Edges along x from (i, j, k): the four cells (i, j-1..j, k-1..k).
    i, j, k = np.nonzero(ins[:-1, 1:-1, 1:-1] != ins[1:, 1:-1, 1:-1])
    j, k = j + 1, k + 1
    quads(vid[i, j - 1, k - 1], vid[i, j, k - 1], vid[i, j, k], vid[i, j - 1, k], ins[i, j, k])
    i, j, k = np.nonzero(ins[1:-1, :-1, 1:-1] != ins[1:-1, 1:, 1:-1])
    i, k = i + 1, k + 1
    quads(vid[i - 1, j, k - 1], vid[i - 1, j, k], vid[i, j, k], vid[i, j, k - 1], ins[i, j, k])
    i, j, k = np.nonzero(ins[1:-1, 1:-1, :-1] != ins[1:-1, 1:-1, 1:])
    i, j = i + 1, j + 1
    quads(vid[i - 1, j - 1, k], vid[i, j - 1, k], vid[i, j, k], vid[i - 1, j, k], ins[i, j, k])
    tri = np.concatenate(tris) if tris else np.zeros((0, 3), np.int64)
    return verts, tri


# ---- meshes ---------------------------------------------------------------------------------

def edges_of(tri: np.ndarray) -> np.ndarray:
    e = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    e.sort(axis=1)
    return np.unique(e, axis=0)


def components(nv: int, tri: np.ndarray) -> np.ndarray:
    """The component of every vertex, by union-find with pointer jumping."""
    parent = np.arange(nv)
    e = edges_of(tri)
    a, b = e[:, 0], e[:, 1]
    while True:
        while True:
            p2 = parent[parent]
            if np.array_equal(p2, parent):
                break
            parent = p2
        ra, rb = parent[a], parent[b]
        diff = ra != rb
        if not diff.any():
            return parent
        lo, hi = np.minimum(ra[diff], rb[diff]), np.maximum(ra[diff], rb[diff])
        np.minimum.at(parent, hi, lo)


def keep_components(v: np.ndarray, tri: np.ndarray, largest_only: bool, min_fraction: float = 0.02) -> tuple[np.ndarray, np.ndarray]:
    if len(tri) == 0:
        return v, tri
    comp = components(len(v), tri)
    roots, counts = np.unique(comp[tri[:, 0]], return_counts=True)
    keep_roots = roots[[np.argmax(counts)]] if largest_only else roots[counts >= max(8, min_fraction * counts.max())]
    keep_t = np.isin(comp[tri[:, 0]], keep_roots)
    tri = tri[keep_t]
    used = np.unique(tri)
    remap = np.full(len(v), -1, np.int64)
    remap[used] = np.arange(len(used))
    return v[used], remap[tri]


def taubin(v: np.ndarray, tri: np.ndarray, iterations: int, lam: float = 0.5, mu: float = -0.53) -> np.ndarray:
    """Taubin's smoothing: a Laplacian step lam, then an inflating step mu, repeated."""
    e = edges_of(tri)
    a = np.concatenate([e[:, 0], e[:, 1]])
    b = np.concatenate([e[:, 1], e[:, 0]])
    deg = np.bincount(a, minlength=len(v)).astype(np.float64)
    deg[deg == 0] = 1
    v = v.astype(np.float64)
    for _ in range(iterations):
        for factor in (lam, mu):
            mean = np.stack([np.bincount(a, weights=v[b, d], minlength=len(v)) for d in range(3)], 1) / deg[:, None]
            v = v + factor * (mean - v)
    return v


def normals(v: np.ndarray, tri: np.ndarray) -> np.ndarray:
    fn = np.cross(v[tri[:, 1]] - v[tri[:, 0]], v[tri[:, 2]] - v[tri[:, 0]])
    n = np.zeros_like(v)
    for c in range(3):
        np.add.at(n, tri[:, c], fn)
    length = np.linalg.norm(n, axis=1)
    length[length == 0] = 1
    return n / length[:, None]


def enclosed_volume(v: np.ndarray, tri: np.ndarray) -> float:
    a, b, c = v[tri[:, 0]], v[tri[:, 1]], v[tri[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6)


class Grid:
    """The common 1 mm grid: index coordinates are millimetres offset by the origin."""

    def __init__(self, affine: np.ndarray, shape):
        if not np.allclose(affine[:3, :3], np.eye(3)):
            raise ValueError("the sources are expected on an axis-aligned 1 mm grid")
        self.affine, self.shape, self.origin = affine, tuple(shape), affine[:3, 3].copy()

    def world(self, idx: np.ndarray) -> np.ndarray:
        return idx + self.origin


def bbox(mask: np.ndarray, margin: int, shape) -> tuple[slice, slice, slice]:
    idx = np.argwhere(mask)
    lo = np.maximum(idx.min(0) - margin, 0)
    hi = np.minimum(idx.max(0) + margin + 1, shape)
    return tuple(slice(int(l), int(h)) for l, h in zip(lo, hi))


def mesh_field(field: np.ndarray, offset: np.ndarray, up: int, level: float, grid: Grid, smooth: int,
               largest_only: bool, occluder: np.ndarray | None, occluder_offset: np.ndarray | None,
               reach=(1.5, 3.0, 5.0), envelope: float = 0.0) -> dict:
    """A mesh of field at level; field's index u is grid index offset + u / up."""
    pad = np.zeros(tuple(s + 2 for s in field.shape), np.float32)
    pad[1:-1, 1:-1, 1:-1] = field
    verts, tri = surface_nets(pad, level)
    if len(tri) == 0:
        return None
    verts = (verts - 1) / up + offset
    verts, tri = keep_components(verts, tri, largest_only)
    verts = taubin(verts, tri, smooth)
    n = normals(verts, tri)
    openness = np.ones(len(verts), np.float32)
    if occluder is not None:
        occ = np.zeros(len(verts), np.float32)
        for d in reach:
            occ += sample(occluder, verts + d * n - occluder_offset)
        openness = np.clip(1 - occ / len(reach), 0, 1)
    out = {"index": verts, "tri": tri, "normal": n, "openness": openness, "level": level}
    if envelope > 0:
        # Sulcal depth: the shape's solid blurred over envelope mm, sampled at each vertex; a
        # vertex on a gyral crown sits near the envelope's edge (about one half), one at the
        # fundus of a sulcus deep inside it (toward one).
        solid = blur((field > level).astype(np.float32), envelope * up)
        out["depth"] = sample(solid, (verts - offset) * up)
    return out


def mesh_mask(mask: np.ndarray, grid: Grid, sigma: float, up: int, smooth: int = 8, large: tuple = (0.9, 8)) -> dict | None:
    voxels = int(mask.sum())
    if voxels == 0:
        return None
    if voxels > SMALL:
        # Large structures are smooth enough at the grid's own resolution, a little more blurred.
        up, sigma, smooth = 1, max(sigma, large[0]), max(smooth, large[1])
    margin = int(np.ceil(3 * sigma)) + 3
    sl = bbox(mask, margin, grid.shape)
    sub = mask[sl].astype(np.float32)
    field = blur(sub, sigma)
    fine = upsample2(field) if up == 2 else field
    level = 0.5
    for level in (0.5, 0.45, 0.4, 0.35, 0.3, 0.25):
        if (fine > level).sum() / up ** 3 >= 0.85 * voxels:
            break
    offset = np.array([s.start for s in sl], np.float64)
    occluder = blur(sub, 2.0)
    m = mesh_field(fine, offset, up, level, grid, smooth, largest_only=False, occluder=occluder, occluder_offset=offset)
    if m is not None:
        m["sigma"], m["upsample"] = sigma, up
    return m


# ---- the hemispheres ---------------------------------------------------------------------------

def assign_parts(aseg: np.ndarray, grid: Grid, iterations: int = 14) -> np.ndarray:
    """Every voxel near the brain to one part: 1 left cerebrum, 2 right cerebrum, 3 left
    cerebellum, 4 right cerebellum, 5 brainstem; midline structures by the side of x;
    the rest by growing the parts into their neighbours (6-connectivity)."""
    part = np.zeros(aseg.shape, np.uint8)
    part[np.isin(aseg, CEREBRUM["lh"])] = 1
    part[np.isin(aseg, CEREBRUM["rh"])] = 2
    part[np.isin(aseg, CEREBELLUM["lh"])] = 3
    part[np.isin(aseg, CEREBELLUM["rh"])] = 4
    part[np.isin(aseg, BRAINSTEM)] = 5
    x = (np.arange(aseg.shape[0]) + grid.origin[0])[:, None, None]
    mid = np.isin(aseg, MIDLINE)
    part[mid & (x < 0)] = 1
    part[mid & (x >= 0)] = 2
    for _ in range(iterations):
        empty = part == 0
        if not empty.any():
            break
        grow = np.zeros_like(part)
        for axis in range(3):
            for shift in (1, -1):
                src = np.roll(part, shift, axis=axis)
                edge = [slice(None)] * 3
                edge[axis] = slice(0, 1) if shift == 1 else slice(-1, None)
                src[tuple(edge)] = 0
                take = empty & (grow == 0) & (src > 0)
                grow[take] = src[take]
        part[empty] = grow[empty]
    return part


def cortex_labels(verts: np.ndarray, n: np.ndarray, ho: np.ndarray, side: str, grid: Grid) -> np.ndarray:
    """The Harvard-Oxford label beneath each vertex: at 1.5 mm inward along the normal,
    else the nearest labelled voxel of the same hemisphere within 4 mm of it."""
    parity = 1 if side == "lh" else 0
    good = (ho > 0) & ((ho % 2) == parity)
    lab = np.where(good, ho, 0).astype(np.uint8)
    shape = np.array(grid.shape)
    out = np.zeros(len(verts), np.uint8)
    offs = [(dx, dy, dz) for dx in range(-4, 5) for dy in range(-4, 5) for dz in range(-4, 5) if dx * dx + dy * dy + dz * dz <= 16]
    offs.sort(key=lambda o: o[0] ** 2 + o[1] ** 2 + o[2] ** 2)
    for depth in (1.5, 0.5, 2.5):
        todo = np.nonzero(out == 0)[0]
        if not len(todo):
            break
        p = np.rint(verts[todo] - depth * n[todo]).astype(np.int64)
        p = np.clip(p, 0, shape - 1)
        out[todo] = lab[p[:, 0], p[:, 1], p[:, 2]]
    base = np.rint(verts - 1.5 * n).astype(np.int64)
    for o in offs[1:]:
        todo = np.nonzero(out == 0)[0]
        if not len(todo):
            break
        p = np.clip(base[todo] + np.array(o), 0, shape - 1)
        out[todo] = lab[p[:, 0], p[:, 1], p[:, 2]]
    return out


def smooth_labels(labels: np.ndarray, tri: np.ndarray, iterations: int) -> np.ndarray:
    """A majority filter along the surface: each labelled vertex takes the label most of
    its neighbours and itself have, so that borders between gyri run smoothly."""
    e = edges_of(tri)
    a = np.concatenate([e[:, 0], e[:, 1], np.arange(len(labels))])
    b = np.concatenate([e[:, 1], e[:, 0], np.arange(len(labels))])
    out = labels.astype(np.int64)
    for _ in range(iterations):
        votes = np.zeros((len(labels), 256), np.int32)
        np.add.at(votes, (a, out[b]), 1)
        votes[:, 0] = 0
        best = votes.argmax(1)
        keep = (out > 0) & (votes[np.arange(len(out)), best] > 0)
        out = np.where(keep, best, out)
    return out.astype(np.uint8)


# ---- the joined label volume -------------------------------------------------------------------

def joined_labels(src: dict, part: np.ndarray) -> tuple[np.ndarray, list]:
    """One label at each voxel: the finest structure there, from the sources in order
    of precedence (aseg; the Harvard-Oxford gyri in the aseg's cortex, grown to fill it;
    the thalamic nuclei in the thalamus; MASSP's finer structures; the cerebellar nuclei)."""
    keys, out = [], np.zeros(part.shape, np.uint16)

    def put(mask, key):
        if not mask.any():
            return
        if key not in keys:
            keys.append(key)
        out[mask] = keys.index(key) + 1

    aseg, ho = src["aseg"], src["ho"]
    for v in np.unique(aseg):
        if v and v not in (3, 42):
            put(aseg == v, f"aseg:{int(v)}")
    # The cortex by gyrus: Harvard-Oxford labels of the right side, grown through the cortex.
    for side, ctx, parity in (("lh", 3, 1), ("rh", 42, 0)):
        cortex = aseg == ctx
        lab = np.where(cortex & (ho > 0) & ((ho % 2) == parity), ho, 0).astype(np.int32)
        for _ in range(12):
            empty = cortex & (lab == 0)
            if not empty.any():
                break
            for axis in range(3):
                for shift in (1, -1):
                    nb = np.roll(lab, shift, axis=axis)
                    take = empty & (lab == 0) & (nb > 0)
                    lab[take] = nb[take]
        for v in np.unique(lab):
            if v:
                put(lab == v, f"ho:{int(v)}")
        rest = cortex & (lab == 0)
        put(rest, f"aseg:{ctx}")
    thal = np.isin(aseg, [10, 49])
    for v in np.unique(src["mial"]):
        if v:
            put((src["mial"] == v) & thal, f"mial:{int(v)}")
    massp = np.rint(src["massp"]).astype(np.int32)
    for v in np.unique(massp):
        if v and v not in MASSP_SKIP:
            put(massp == v, f"massp:{int(v)}")
    for side in ("jlh", "jrh"):
        for name, labels in JULICH_NUCLEI.items():
            put(np.isin(src[side], labels) & np.isin(aseg, [7, 8, 46, 47]), f"{side}:{name}")
    return out, keys


def contacts(labels: np.ndarray, n: int) -> np.ndarray:
    """Faces shared by every pair of labels (6-adjacency), as an (n+1) x (n+1) matrix."""
    m = np.zeros((n + 1, n + 1), np.int64)
    for axis in range(3):
        a = np.moveaxis(labels, axis, 0)
        x, y = a[:-1].ravel().astype(np.int64), a[1:].ravel().astype(np.int64)
        d = x != y
        code = x[d] * (n + 1) + y[d]
        cnt = np.bincount(code, minlength=(n + 1) ** 2).reshape(n + 1, n + 1)
        m += cnt + cnt.T
    np.fill_diagonal(m, 0)
    return m


def stats(mask: np.ndarray, grid: Grid) -> dict:
    idx = np.argwhere(mask)
    c = idx.mean(0) + grid.origin
    lo, hi = idx.min(0) + grid.origin, idx.max(0) + grid.origin
    return {"voxels": int(len(idx)), "volume": float(len(idx)), "centroid": [round(float(x), 2) for x in c],
            "bounds": {"min": [float(x) for x in lo], "max": [float(x) for x in hi]}}


# ---- writing ----------------------------------------------------------------------------------

class Writer:
    """The bytes of one pack of meshes, each block aligned to 4 bytes."""

    def __init__(self):
        self.parts, self.size = [], 0

    def add(self, array: np.ndarray) -> int:
        b = np.ascontiguousarray(array).tobytes()
        offset = self.size
        self.parts.append(b)
        self.size += len(b)
        pad = (-self.size) % 4
        if pad:
            self.parts.append(b"\0" * pad)
            self.size += pad
        return offset

    def bytes(self) -> bytes:
        return b"".join(self.parts)


PACKS = ("core", "white", "bundles")


def pack_of(key: str) -> str:
    return "white" if key.startswith("white:") else "bundles" if key.startswith("dwm:") else "core"


def pack(key: str, m: dict, grid: Grid, writers: dict, labels: np.ndarray | None = None) -> dict:
    pack_name = pack_of(key)
    w = writers[pack_name]
    world = grid.world(m["index"])
    q = np.rint(world * 100).astype(np.int16)
    nrm = np.zeros((len(world), 4), np.int8)
    nrm[:, :3] = np.clip(np.rint(m["normal"] * 127), -127, 127)
    ao = np.rint(np.clip(m["openness"], 0, 1) * 255).astype(np.uint8)
    nv = len(world)
    index_type = "uint16" if nv < 65536 else "uint32"
    tri = m["tri"].astype(np.uint16 if index_type == "uint16" else np.uint32)
    entry = {"key": key, "pack": pack_name, "vertices": nv, "triangles": int(len(tri)), "position": w.add(q), "normal": w.add(nrm),
             "openness": w.add(ao), "index": w.add(tri), "indexType": index_type,
             "level": round(float(m["level"]), 3), "sigma": m.get("sigma"), "upsample": m.get("upsample", 1),
             "enclosedVolume": round(enclosed_volume(world, m["tri"]), 1)}
    if labels is not None:
        entry["label"] = w.add(labels.astype(np.uint8))
    if "depth" in m:
        entry["depth"] = w.add(np.rint(np.clip(m["depth"], 0, 1) * 255).astype(np.uint8))
    return entry


def provenance() -> list:
    out = []
    for key, (dataset, name) in FILES.items():
        r = CACHE / dataset / "retrieval.json"
        rec = json.loads(r.read_text(encoding="utf-8")) if r.exists() else {}
        f = next((x for x in rec.get("files", []) if x["file"] == name), {})
        out.append({"source": key, "dataset": dataset, "file": name, "sha256": f.get("sha256"), "retrieved_at": rec.get("retrieved_at")})
    return out


def build() -> Path:
    t0 = time.time()
    src, grid = {}, None
    for key, (dataset, name) in FILES.items():
        path = CACHE / dataset / name
        if not path.exists():
            raise SystemExit(f"{path} is missing; run: python data/acquisition/fetch.py {dataset}")
        data, affine, _ = read_nifti(path)
        if grid is None:
            grid = Grid(affine, data.shape)
        elif data.shape != grid.shape or not np.allclose(affine, grid.affine):
            raise SystemExit(f"{name} is not on the grid of the template")
        src[key] = data
    print(f"read the sources ({time.time() - t0:.1f} s)")
    aseg = src["aseg"].astype(np.int32)
    part = assign_parts(aseg, grid)
    w = {name: Writer() for name in PACKS}
    meshes, regions = [], {}

    # The surfaces of the hemispheres and of the cerebellum, from the tissue probabilities.
    tissue = (src["gm"] + src["wm"]).astype(np.float32)
    for side, cls in (("lh", 1), ("rh", 2)):
        weight = blur(part == cls, 0.8)
        sl = bbox(part == cls, 4, grid.shape)
        offset = np.array([s.start for s in sl], np.float64)
        fill = blur(np.isin(aseg, FILL[side]), 1.0)
        pial = np.maximum(tissue, fill)[sl] * weight[sl]
        m = mesh_field(pial, offset, 1, 0.5, grid, 6, True, blur(pial, 2.5), offset, reach=(2.0, 4.0, 7.0), envelope=5.0)
        labels = cortex_labels(grid.world(m["index"]) - grid.origin, m["normal"], src["ho"], side, grid)
        labels = smooth_labels(labels, m["tri"], 3)
        meshes.append(pack(f"cortex:{side}", m, grid, w, labels))
        regions[f"cortex:{side}"] = stats(aseg == (3 if side == "lh" else 42), grid)
        print(f"  pial surface {side}: {len(m['index'])} vertices, {np.mean(labels > 0):.1%} labelled ({time.time() - t0:.1f} s)")
        deep = blur(np.isin(aseg, DEEP[side]), 1.0)
        white = np.maximum(src["wm"].astype(np.float32), deep)[sl] * weight[sl]
        m = mesh_field(white, offset, 1, 0.5, grid, 6, True, blur(white, 2.0), offset, reach=(1.5, 3.0, 5.0), envelope=4.0)
        meshes.append(pack(f"white:{side}", m, grid, w))
        regions[f"white:{side}"] = stats(aseg == (2 if side == "lh" else 41), grid)
        print(f"  white surface {side}: {len(m['index'])} vertices")
    for side, cls in (("lh", 3), ("rh", 4)):
        weight = blur(part == cls, 0.8)
        sl = bbox(part == cls, 4, grid.shape)
        offset = np.array([s.start for s in sl], np.float64)
        field = tissue[sl] * weight[sl]
        m = mesh_field(field, offset, 1, 0.5, grid, 6, True, blur(field, 2.0), offset, reach=(1.5, 3.0, 5.0), envelope=3.0)
        meshes.append(pack(f"cerebellum:{side}", m, grid, w))
        regions[f"cerebellum:{side}"] = stats(np.isin(aseg, CEREBELLUM[side]), grid)
        print(f"  cerebellum {side}: {len(m['index'])} vertices")

    # The labelled structures.
    jobs = [("aseg", v, aseg == v) for v in ASEG_MESHES]
    massp = np.rint(src["massp"]).astype(np.int32)
    jobs += [("massp", int(v), massp == v) for v in np.unique(massp) if v and v not in MASSP_SKIP]
    jobs += [("mial", int(v), src["mial"] == v) for v in np.unique(src["mial"]) if v]
    jobs += [("dwm", int(v), src["dwm"] == v) for v in np.unique(src["dwm"]) if v]
    for side in ("jlh", "jrh"):
        jobs += [(side, name, np.isin(src[side], labels) & np.isin(aseg, [7, 8, 46, 47])) for name, labels in JULICH_NUCLEI.items()]
    for source, label, mask in jobs:
        key = f"{source}:{label}"
        sigma, up = PARAMS[source]
        m = mesh_mask(mask, grid, sigma, up, large=LARGE.get(source, (0.9, 8)))
        if m is None:
            print(f"  {key}: empty")
            continue
        meshes.append(pack(key, m, grid, w))
        regions[key] = stats(mask, grid)
    print(f"{len(meshes)} meshes, {sum(m['vertices'] for m in meshes)} vertices, {sum(m['triangles'] for m in meshes)} triangles ({time.time() - t0:.1f} s)")

    # The joined labels, the regions' measures and their contacts.
    joined, keys = joined_labels(src, part)
    n = len(keys)
    touch = contacts(joined, n)
    for i, key in enumerate(keys, start=1):
        mask = joined == i
        r = regions.setdefault(key, stats(mask, grid))
        r["joined"] = stats(mask, grid) if key in regions else None
        row = touch[i]
        order = np.argsort(-row)
        r["contacts"] = [[keys[j - 1], int(row[j])] for j in order if j > 0 and row[j] > 0][:16]
        if r.get("joined") is None:
            r.pop("joined")
    # Bundles: the gyri within 4 mm of each.
    ho_index = {k: i + 1 for i, k in enumerate(keys) if k.startswith("ho:")}
    for v in np.unique(src["dwm"]):
        if not v:
            continue
        mask = src["dwm"] == v
        near = mask.copy()
        for _ in range(4):
            g = near.copy()
            for axis in range(3):
                g |= np.roll(near, 1, axis) | np.roll(near, -1, axis)
            near = g
        found = joined[near & ~mask]
        counts = np.bincount(found, minlength=n + 1)
        hits = [[k, int(counts[i])] for k, i in ho_index.items() if counts[i] > 0]
        hits.sort(key=lambda x: -x[1])
        regions[f"dwm:{int(v)}"]["near"] = hits[:12]
    print(f"joined labels: {n} structures ({time.time() - t0:.1f} s)")

    OUT.mkdir(parents=True, exist_ok=True)
    name = "anatomy-mni152"
    for p in PACKS:
        (OUT / f"{name}-{p}.bin.gz").write_bytes(gzip.compress(w[p].bytes(), compresslevel=6))
    write_nifti(OUT / f"{name}-labels.nii.gz", joined, grid.affine, code=4, descrip="anatomy-mni152 joined labels")
    meta = {
        "name": name, "kind": "anatomy", "schema": "anatomy/1",
        "title": "The anatomy of the human brain in the MNI 152 ICBM 2009c asymmetric space: surfaces of the cortex, the white matter, the cerebellum, the deep nuclei, the ventricles and the fibre bundles",
        "space": "mni152-nlin-2009c-asym",
        "grid": {"dims": list(grid.shape), "affine": [float(x) for x in grid.affine.reshape(-1)]},
        "packs": {p: f"{name}-{p}.bin.gz" for p in PACKS}, "labels": f"{name}-labels.nii.gz", "positionScale": 0.01,
        "meshes": meshes, "regions": regions, "labelTable": [[i + 1, k] for i, k in enumerate(keys)],
        "sources": provenance(),
        "method": "surface nets on blurred, upsampled masks (tissue probabilities for the cortex and cerebellum), Taubin smoothing, volumetric obscurance; see data/products/brain_anatomy.py",
        "built_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }
    out = OUT / f"{name}.json"
    out.write_text(json.dumps(meta, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("anatomy", {"name": name, "file": f"{name}.json", "title": meta["title"], "space": meta["space"]})
    sizes = ", ".join(f"{p} {(OUT / f'{name}-{p}.bin.gz').stat().st_size / 1e6:.1f} MB" for p in PACKS)
    print(f"wrote {out.name} and the packs ({sizes}) in {time.time() - t0:.1f} s")
    return out


if __name__ == "__main__":
    build()
