"""NIfTI-1 and NIfTI-2 volumes, read with the standard library and numpy.

A NIfTI file is a binary header (348 bytes for version 1, 540 for version 2)
followed in a single .nii file by the voxel data at vox_offset; .nii.gz is the
same, gzipped. The header places the voxels in the world by an affine: the
sform when its code is set, else the qform (a quaternion rotation with the
voxel sizes and an offset), else the voxel sizes alone, as nibabel chooses.
World coordinates are RAS millimetres. The field offsets follow nifti1.h and
nifti2.h (NIFTI-Imaging/nifti_clib).

The reader returns the voxel array in index order, a[i, j, k], so that voxel
(i, j, k) is centred at affine @ (i, j, k, 1); NIfTI stores the first index
fastest, which numpy reads with order="F". The domain's JavaScript reader
(domains/neuroscience/js/nifti.js) is the reference for the same format.

Usage::

    from nifti import read_nifti
    data, affine, header = read_nifti(path)
"""

from __future__ import annotations

import gzip
import struct
from pathlib import Path

import numpy as np

DTYPES = {2: np.uint8, 4: np.int16, 8: np.int32, 16: np.float32, 64: np.float64, 256: np.int8,
          512: np.uint16, 768: np.uint32, 1024: np.int64, 1280: np.uint64}


def _quaternion_affine(b: float, c: float, d: float, sizes, qfac: float, offset) -> np.ndarray:
    a = np.sqrt(max(0.0, 1.0 - (b * b + c * c + d * d)))
    r = np.array([
        [a * a + b * b - c * c - d * d, 2 * b * c - 2 * a * d, 2 * b * d + 2 * a * c],
        [2 * b * c + 2 * a * d, a * a + c * c - b * b - d * d, 2 * c * d - 2 * a * b],
        [2 * b * d - 2 * a * c, 2 * c * d + 2 * a * b, a * a + d * d - c * c - b * b]])
    m = np.eye(4)
    m[:3, :3] = r * np.array([sizes[0], sizes[1], qfac * sizes[2]])
    m[:3, 3] = offset
    return m


def read_header(raw: bytes) -> dict:
    """The header's fields that place and decode the data."""
    for endian in "<>":
        size = struct.unpack(endian + "i", raw[:4])[0]
        if size in (348, 540):
            break
    else:
        raise ValueError("not a NIfTI header")
    e = endian
    if size == 348:
        dim = struct.unpack(e + "8h", raw[40:56])
        datatype, = struct.unpack(e + "h", raw[70:72])
        pixdim = struct.unpack(e + "8f", raw[76:108])
        vox_offset, slope, inter = struct.unpack(e + "3f", raw[108:120])
        qform_code, sform_code = struct.unpack(e + "2h", raw[252:256])
        quat = struct.unpack(e + "6f", raw[256:280])
        srow = struct.unpack(e + "12f", raw[280:328])
        magic = raw[344:347]
    else:
        datatype, = struct.unpack(e + "h", raw[12:14])
        dim = struct.unpack(e + "8q", raw[16:80])
        pixdim = struct.unpack(e + "8d", raw[104:168])
        vox_offset, = struct.unpack(e + "q", raw[168:176])
        slope, inter = struct.unpack(e + "2d", raw[176:192])
        qform_code, sform_code = struct.unpack(e + "2i", raw[344:352])
        quat = struct.unpack(e + "6d", raw[352:400])
        srow = struct.unpack(e + "12d", raw[400:496])
        magic = raw[4:7]
    if magic not in (b"n+1", b"ni1", b"n+2", b"ni2"):
        raise ValueError(f"not a NIfTI file: magic {magic!r}")
    nd = max(1, min(7, dim[0]))
    shape = tuple(int(x) for x in dim[1:1 + nd])
    if sform_code > 0:
        affine = np.vstack([np.array(srow, dtype=np.float64).reshape(3, 4), [0, 0, 0, 1]])
    elif qform_code > 0:
        qfac = -1.0 if pixdim[0] == -1 else 1.0
        affine = _quaternion_affine(quat[0], quat[1], quat[2], pixdim[1:4], qfac, quat[3:6])
    else:
        affine = np.diag([pixdim[1], pixdim[2], pixdim[3], 1.0])
    return {"version": 1 if size == 348 else 2, "endian": e, "shape": shape, "datatype": datatype,
            "vox_offset": int(vox_offset), "slope": float(slope), "inter": float(inter), "affine": affine}


def read_nifti(path: str | Path, scale: bool = True) -> tuple[np.ndarray, np.ndarray, dict]:
    """The voxel array (index order), the voxel-to-world affine, and the header."""
    path = Path(path)
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    h = read_header(raw)
    dtype = DTYPES.get(h["datatype"])
    if dtype is None:
        raise ValueError(f"unsupported NIfTI datatype {h['datatype']}")
    count = int(np.prod(h["shape"]))
    data = np.frombuffer(raw, dtype=np.dtype(dtype).newbyteorder(h["endian"]), count=count, offset=h["vox_offset"])
    data = data.reshape(h["shape"], order="F")
    if scale and h["slope"] not in (0.0, 1.0) or (scale and h["inter"] != 0.0 and h["slope"] != 0.0):
        data = data.astype(np.float32) * h["slope"] + h["inter"]
    return np.ascontiguousarray(data), h["affine"], h


def write_nifti(path: str | Path, data: np.ndarray, affine: np.ndarray, code: int = 4, descrip: str = "") -> None:
    """A gzipped single-file NIfTI-1 of a 3D array, placed by its sform (and the
    same affine as qform when it is a scaled rotation, which the products'
    axis-aligned grids are). code 4 is MNI 152."""
    dtypes = {np.dtype(v): k for k, v in DTYPES.items()}
    dt = np.dtype(data.dtype)
    if dt not in dtypes:
        raise ValueError(f"unsupported dtype {dt}")
    header = bytearray(352)
    shape = data.shape + (1,) * (7 - data.ndim)
    struct.pack_into("<i", header, 0, 348)
    struct.pack_into("<8h", header, 40, data.ndim, *shape)
    struct.pack_into("<2h", header, 70, dtypes[dt], dt.itemsize * 8)
    sizes = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
    struct.pack_into("<8f", header, 76, 1.0, *sizes, 1.0, 1.0, 1.0, 1.0)
    struct.pack_into("<3f", header, 108, 352.0, 1.0, 0.0)
    header[123] = 2                                         # xyzt_units: mm
    header[148:148 + min(79, len(descrip))] = descrip.encode("ascii", "replace")[:79]
    struct.pack_into("<2h", header, 252, code, code)
    r = affine[:3, :3] / sizes
    if np.allclose(r, np.diag(np.diag(r))) and np.all(np.diag(r) > 0):
        struct.pack_into("<6f", header, 256, 0.0, 0.0, 0.0, *affine[:3, 3])
    else:
        struct.pack_into("<2h", header, 252, 0, code)       # sform only
    struct.pack_into("<12f", header, 280, *affine[:3, :4].reshape(-1))
    header[344:348] = b"n+1\x00"
    body = np.asfortranarray(data).tobytes(order="F")
    Path(path).write_bytes(gzip.compress(bytes(header) + body, compresslevel=6))
