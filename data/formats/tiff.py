"""A reader and a writer for GeoTIFF, including Cloud-Optimised GeoTIFFs read
over the web by ranges.

A TIFF file (Adobe, "TIFF Revision 6.0", 1992, section 2) is a header and a
chain of image file directories (IFDs), each a table of tagged fields that
describe one image: its size and samples, how it is cut into strips or tiles,
how each is compressed and where each lies in the file. A Cloud-Optimised
GeoTIFF (COG; OGC 21-026, "OGC Cloud Optimized GeoTIFF Standard", 2023) is a
tiled GeoTIFF whose directories come first, the full image followed by its
reduced-resolution versions (overviews), so that a reader that has the first
few kilobytes knows where every tile lies, and reads a window of any image by
range requests for its tiles alone. A 20 km box of a 190 MB Sentinel-2 band, at
20 m (1000 x 1000 pixels from nine tiles of the first overview), takes five
requests: the address, the directories, and three runs of tiles, 3.2 MB in all.

Reading supports classic TIFF and BigTIFF, in either byte order; strips and
tiles, interleaved or in planes; compression none (1), LZW (5; TIFF 6.0,
section 13) and Deflate (8, and the older 32946: zlib streams, RFC 1950; Adobe,
"PageMaker 6.0 TIFF Technical Notes", 1995); the horizontal differencing
predictor (2; TIFF 6.0, section 14) and the floating-point predictor (3; Adobe,
"Photoshop TIFF Technical Note 3", 2005); unsigned and signed integers of 8, 16,
32 and 64 bits and IEEE floats of 32 and 64 bits (TIFF 6.0, section 19). JPEG,
PackBits and the CCITT codings are not supported, nor the LZW of TIFF 5.0, whose
bits run the other way.

Georeferencing follows the OGC GeoTIFF Standard 1.1 (OGC 19-008r4, 2019): the
GeoKey directory, with the EPSG code of the projected or geographic CRS and the
raster type (a pixel is an area, or a point at its centre); and the model tie
point with the pixel scale, or a model transformation without rotation. Here
the georeferencing of an image is its *transform* (x0, dx, y0, dy): the outer
corner of its first pixel, west and north, and the size of a pixel, so that
pixel (column c, row r) covers [x0 + c dx, x0 + (c + 1) dx] in x and
[y0 - (r + 1) dy, y0 - r dy] in y, whatever the raster type of the file.

Coordinates of the UTM zones on WGS 84 (EPSG:32601-32660 north, 32701-32760
south) are converted to and from longitude and latitude by the transverse
Mercator projection in Krüger's series to sixth order in the third flattening n,
as Karney (2011) gives it; its error is below 5 nm within 3900 km of the central
meridian (Karney 2011, section 4). (The projection belongs to geodesy; it is
here until that domain has its Python or C++.)

Writing stores a float32 grid of the library (a regular longitude-latitude grid
on WGS 84, rows north to south) as a GeoTIFF in EPSG:4326 that GDAL, and so
QGIS, read: strips compressed by Deflate with the floating-point predictor, the
GeoKeys of GeoTIFF 1.1, NaN declared as no value by GDAL's GDAL_NODATA tag
(42113), and the variable's name and units in its GDAL_METADATA tag (42112).
The grid's values are at the centres of pixels that are areas, which is how
GDAL reads a global grid such as the GFS's.

Header parsing and LZW decoding need only the standard library; pixels are
numpy arrays.

References: Adobe Developers Association, "TIFF Revision 6.0" (1992); OGC 19-008r4,
"OGC GeoTIFF Standard", version 1.1 (2019); Karney, C. F. F. (2011), "Transverse
Mercator with an accuracy of a few nanometers", J. Geodesy 85, 475-485,
doi:10.1007/s00190-011-0445-3.
"""

from __future__ import annotations

import io
import math
import struct
import sys
import zlib
from pathlib import Path

# TIFF field types (TIFF 6.0, section 2; BigTIFF adds 16-18) -> (struct code, size)
TYPES = {1: ("B", 1), 2: ("s", 1), 3: ("H", 2), 4: ("I", 4), 5: ("I", 8), 6: ("b", 1), 7: ("B", 1),
         8: ("h", 2), 9: ("i", 4), 10: ("i", 8), 11: ("f", 4), 12: ("d", 8), 16: ("Q", 8), 17: ("q", 8), 18: ("Q", 8)}

NEW_SUBFILE_TYPE, IMAGE_WIDTH, IMAGE_LENGTH, BITS_PER_SAMPLE, COMPRESSION = 254, 256, 257, 258, 259
PHOTOMETRIC, IMAGE_DESCRIPTION, STRIP_OFFSETS, SAMPLES_PER_PIXEL, ROWS_PER_STRIP = 262, 270, 273, 277, 278
STRIP_BYTE_COUNTS, PLANAR_CONFIGURATION, SOFTWARE, PREDICTOR = 279, 284, 305, 317
TILE_WIDTH, TILE_LENGTH, TILE_OFFSETS, TILE_BYTE_COUNTS, SAMPLE_FORMAT = 322, 323, 324, 325, 339
MODEL_PIXEL_SCALE, MODEL_TIEPOINT, MODEL_TRANSFORMATION = 33550, 33922, 34264
GEO_KEY_DIRECTORY, GEO_DOUBLE_PARAMS, GEO_ASCII_PARAMS = 34735, 34736, 34737
GDAL_METADATA, GDAL_NODATA = 42112, 42113

# GeoKeys (OGC 19-008r4, section 7 and annex C)
GT_MODEL_TYPE, GT_RASTER_TYPE, GT_CITATION = 1024, 1025, 1026
GEODETIC_CRS, GEODETIC_CITATION, GEOG_ANGULAR_UNITS = 2048, 2049, 2054
PROJECTED_CRS, PROJ_LINEAR_UNITS = 3072, 3076
MODEL_PROJECTED, MODEL_GEOGRAPHIC = 1, 2
PIXEL_IS_AREA, PIXEL_IS_POINT = 1, 2
USER_DEFINED = 32767

NONE, LZW, DEFLATE, DEFLATE_OLD = 1, 5, 8, 32946
HORIZONTAL, FLOATING_POINT = 2, 3
GAP = 1 << 16            # ranges closer than this are read in one request


class TIFFError(ValueError):
    pass


# ---- bytes by ranges --------------------------------------------------------


class _Bytes:
    """The bytes of a seekable file, read by ranges. A TIFF's directories are
    read in many small pieces, mostly near its start: small reads are kept,
    and served from what was read before; the pixels of a window are read as
    few ranges as their positions allow, and not kept."""

    CHUNK = 1 << 16

    def __init__(self, file):
        self.file = file
        self.kept: list[tuple[int, bytes]] = []
        self.requests = 0

    def _read(self, start: int, n: int) -> bytes:
        self.file.seek(start)
        self.requests += 1
        out = bytearray()
        while len(out) < n:
            data = self.file.read(n - len(out))
            if not data:
                break
            out += data
        return bytes(out)

    def get(self, start: int, n: int) -> bytes:
        for s, block in self.kept:
            if s <= start and start + n <= s + len(block):
                return block[start - s:start - s + n]
        block = self._read(start, max(n, self.CHUNK))
        if len(block) < n:
            raise TIFFError(f"the file ends at {start + len(block)}, before byte {start + n}")
        self.kept.append((start, block))
        return block[:n]

    def get_many(self, ranges: list[tuple[int, int]]) -> list[bytes]:
        """The bytes of each (start, length), reading neighbouring ranges together."""
        order = sorted(range(len(ranges)), key=lambda i: ranges[i][0])
        out: list[bytes] = [b""] * len(ranges)
        group: list[int] = []

        def flush():
            start = ranges[group[0]][0]
            end = max(ranges[i][0] + ranges[i][1] for i in group)
            data = self._read(start, end - start)
            for i in group:
                s, n = ranges[i]
                out[i] = data[s - start:s - start + n]

        for i in order:
            if ranges[i][1] == 0:
                continue
            if group and ranges[i][0] > max(ranges[j][0] + ranges[j][1] for j in group) + GAP:
                flush()
                group = []
            group.append(i)
        if group:
            flush()
        return out


# ---- decompression and predictors -------------------------------------------


def lzw_decode(data: bytes) -> bytes:
    """Decodes TIFF LZW (TIFF 6.0, section 13). Codes are 9 to 12 bits wide,
    most significant bit first; 256 clears the table and 257 ends the strip or
    tile; the first new entry is 258. The width grows by one bit when the next
    free entry would need it one code early, at 511, 1023 and 2047, as the
    encoders of libtiff and GDAL write it."""
    if len(data) >= 2 and data[0] == 0 and data[1] & 1:
        raise TIFFError("LZW of TIFF 5.0 (bits least significant first) is not supported")
    table = [bytes([i]) for i in range(256)] + [b"", b""]
    out = bytearray()
    bits, acc, held, i, n = 9, 0, 0, 0, len(data)
    old = None
    while True:
        while held < bits:
            if i >= n:
                return bytes(out)           # no end code: the data end with the strip
            acc = (acc << 8) | data[i]
            i += 1
            held += 8
        held -= bits
        code = (acc >> held) & ((1 << bits) - 1)
        acc &= (1 << held) - 1
        if code == 257:
            break
        if code == 256:
            del table[258:]
            bits, old = 9, None
            continue
        if old is None:
            if code > 255:
                raise TIFFError("corrupt LZW data: a code after a clear is not a byte")
            entry = table[code]
        elif code < len(table):
            entry = table[code]
            if len(table) < 4096:
                table.append(old + entry[:1])
        elif code == len(table):
            entry = old + old[:1]       # the code being defined: KwKwK
            table.append(entry)
        else:
            raise TIFFError(f"corrupt LZW data: code {code} with {len(table)} entries")
        out += entry
        old = entry
        if len(table) >= (1 << bits) - 1 and bits < 12:
            bits += 1
    return bytes(out)


def inflate(data: bytes) -> bytes:
    """Decodes Deflate (compression 8 and 32946): a zlib stream (RFC 1950)."""
    return zlib.decompress(data)


def undo_horizontal(a):
    """Undoes the horizontal differencing predictor (TIFF 6.0, section 14) on an
    array (rows, columns, samples) in native byte order: each sample was stored
    as its difference from the same sample of the pixel to its left, modulo
    2^bits. Floats are summed as the integers of their bits, as libtiff does."""
    import numpy as np
    kind = a.dtype
    as_int = a.view(f"u{kind.itemsize}") if kind.kind == "f" else a
    summed = np.cumsum(as_int, axis=1, dtype=as_int.dtype)
    return summed.view(kind) if kind.kind == "f" else summed


def undo_floating_point(data: bytes, rows: int, width: int, samples: int, size: int):
    """Undoes the floating-point predictor (Adobe Photoshop TIFF Technical Note 3,
    2005) on decompressed bytes: each row holds the bytes of its values in planes,
    most significant first, and each byte was stored as its difference from the
    byte `samples` places before it in the row, modulo 256. Returns the values
    as big-endian floats (rows, width, samples), whatever the file's byte order."""
    import numpy as np
    count = width * samples
    raw = np.frombuffer(data, np.uint8, rows * count * size).reshape(rows, count * size // samples, samples)
    raw = np.cumsum(raw, axis=1, dtype=np.uint8).reshape(rows, size, count)
    return np.ascontiguousarray(raw.transpose(0, 2, 1)).view(f">f{size}").reshape(rows, width, samples)


def horizontal_difference(a):
    """The horizontal differencing of an array (rows, columns, samples) in native
    byte order, the inverse of undo_horizontal()."""
    kind = a.dtype
    as_int = a.view(f"u{kind.itemsize}") if kind.kind == "f" else a
    out = as_int.copy()
    out[:, 1:] = as_int[:, 1:] - as_int[:, :-1]
    return out.view(kind) if kind.kind == "f" else out


def floating_point_difference(values) -> bytes:
    """The floating-point predictor's bytes of float rows (rows, width, samples),
    the inverse of undo_floating_point()."""
    import numpy as np
    rows, width, samples = values.shape
    size = values.dtype.itemsize
    big = np.ascontiguousarray(values.astype(f">f{size}")).view(np.uint8).reshape(rows, width * samples, size)
    planes = np.ascontiguousarray(big.transpose(0, 2, 1)).reshape(rows, -1, samples)
    out = planes.copy()
    out[:, 1:] = planes[:, 1:] - planes[:, :-1]
    return out.tobytes()


# ---- directories ------------------------------------------------------------


class Image:
    """One image of the file (an IFD): the full image, an overview or a mask."""

    def __init__(self, tags: dict, order: str, source: _Bytes):
        self.tags, self.order, self.source = tags, order, source
        one = lambda tag, default=None: tags[tag][0] if tag in tags else default  # noqa: E731
        self.width, self.height = one(IMAGE_WIDTH), one(IMAGE_LENGTH)
        self.samples = one(SAMPLES_PER_PIXEL, 1)
        bits = tags.get(BITS_PER_SAMPLE, (1,))
        if len(set(bits)) != 1:
            raise TIFFError(f"samples of different sizes {bits} are not supported")
        self.bits = bits[0]
        self.format = one(SAMPLE_FORMAT, 1)            # 1 unsigned, 2 signed, 3 IEEE float
        self.compression = one(COMPRESSION, NONE)
        self.predictor = one(PREDICTOR, 1)
        self.planar = one(PLANAR_CONFIGURATION, 1)
        self.subfile = one(NEW_SUBFILE_TYPE, 0)         # bit 0: reduced resolution; bit 2: a mask
        self.tiled = TILE_WIDTH in tags
        if self.tiled:
            self.block_width, self.block_height = one(TILE_WIDTH), one(TILE_LENGTH)
            self.offsets, self.counts = tags[TILE_OFFSETS], tags[TILE_BYTE_COUNTS]
        else:
            self.block_width = self.width
            self.block_height = min(one(ROWS_PER_STRIP, self.height), self.height)
            self.offsets, self.counts = tags[STRIP_OFFSETS], tags[STRIP_BYTE_COUNTS]
        nodata = tags.get(GDAL_NODATA)
        self.nodata = float(nodata.strip()) if isinstance(nodata, str) and nodata.strip() else None
        codes = {(1, 8): "u1", (1, 16): "u2", (1, 32): "u4", (1, 64): "u8", (2, 8): "i1", (2, 16): "i2",
                 (2, 32): "i4", (2, 64): "i8", (3, 32): "f4", (3, 64): "f8"}
        if (self.format, self.bits) not in codes:
            raise TIFFError(f"samples of format {self.format} and {self.bits} bits are not supported")
        self.code = codes[(self.format, self.bits)]
        if self.compression not in (NONE, LZW, DEFLATE, DEFLATE_OLD):
            raise TIFFError(f"compression {self.compression} is not supported (only none, LZW and Deflate)")
        if self.predictor not in (1, HORIZONTAL, FLOATING_POINT):
            raise TIFFError(f"predictor {self.predictor} is not supported")

    def __repr__(self):
        return f"Image({self.width} x {self.height}, {self.samples} x {self.code}, compression {self.compression})"

    @property
    def across(self) -> int:
        return -(-self.width // self.block_width)

    @property
    def down(self) -> int:
        return -(-self.height // self.block_height)

    def decode(self, raw: bytes, rows: int, samples: int):
        """One strip or tile, (rows, block width, samples), in native byte order."""
        import numpy as np
        if self.compression == LZW:
            data = lzw_decode(raw)
        elif self.compression in (DEFLATE, DEFLATE_OLD):
            data = inflate(raw)
        else:
            data = raw
        width, size = self.block_width, self.bits // 8
        need = rows * width * samples * size
        if len(data) < need:
            raise TIFFError(f"a block decodes to {len(data)} bytes, not {need}")
        if self.predictor == FLOATING_POINT:
            if self.format != 3:
                raise TIFFError("the floating-point predictor applies to floats only")
            block = undo_floating_point(data, rows, width, samples, size)
        else:
            block = np.frombuffer(data, self.order + self.code, rows * width * samples).reshape(rows, width, samples)
            block = block.astype(self.code)              # native byte order, before the predictor sums
            if self.predictor == HORIZONTAL:
                block = undo_horizontal(block)
        return block.astype(self.code, copy=False)

    def window(self, col0: int, row0: int, ncols: int, nrows: int):
        """Pixels [row0, row0 + nrows) x [col0, col0 + ncols): an array (nrows, ncols),
        or (nrows, ncols, samples) for several samples. Sparse blocks (no bytes)
        hold the no-data value, or zero."""
        import numpy as np
        if col0 < 0 or row0 < 0 or ncols <= 0 or nrows <= 0 or col0 + ncols > self.width or row0 + nrows > self.height:
            raise TIFFError(f"window {col0}, {row0}, {ncols} x {nrows} is not within {self.width} x {self.height}")
        bw, bh = self.block_width, self.block_height
        planes = self.samples if self.planar == 2 else 1
        per_block = 1 if self.planar == 2 else self.samples
        wanted = []
        for plane in range(planes):
            for by in range(row0 // bh, (row0 + nrows - 1) // bh + 1):
                for bx in range(col0 // bw, (col0 + ncols - 1) // bw + 1):
                    wanted.append((plane, by, bx, plane * self.across * self.down + by * self.across + bx))
        raws = self.source.get_many([(self.offsets[k], self.counts[k]) for *_, k in wanted])
        fill = self.nodata if self.nodata is not None and (self.format == 3 or math.isfinite(self.nodata)) else 0
        out = np.full((nrows, ncols, self.samples), fill, dtype=np.dtype(self.code))
        for (plane, by, bx, k), raw in zip(wanted, raws):
            if not raw:
                continue
            rows = bh if self.tiled else min(bh, self.height - by * bh)
            block = self.decode(raw, rows, per_block)
            r0, r1 = max(row0, by * bh), min(row0 + nrows, by * bh + rows)
            c0, c1 = max(col0, bx * bw), min(col0 + ncols, bx * bw + bw)
            target = slice(plane, plane + 1) if self.planar == 2 else slice(None)
            out[r0 - row0:r1 - row0, c0 - col0:c1 - col0, target] = block[r0 - by * bh:r1 - by * bh, c0 - bx * bw:c1 - bx * bw]
        return out[:, :, 0] if self.samples == 1 else out

    def read(self):
        """The whole image."""
        return self.window(0, 0, self.width, self.height)


class GeoTIFF:
    """A GeoTIFF, from a seekable binary file: a local file, bytes, or a file on a
    web server whose reads are range requests (GeoTIFF.url)."""

    def __init__(self, file):
        self.source = _Bytes(file)
        head = self.source.get(0, 16)
        self.order = {b"II": "<", b"MM": ">"}.get(head[:2])
        if not self.order:
            raise TIFFError("not a TIFF file")
        version = struct.unpack_from(self.order + "H", head, 2)[0]
        if version == 42:
            self.big, offset = False, struct.unpack_from(self.order + "I", head, 4)[0]
        elif version == 43:
            if struct.unpack_from(self.order + "H", head, 4)[0] != 8:
                raise TIFFError("BigTIFF with offsets other than 8 bytes")
            self.big, offset = True, struct.unpack_from(self.order + "Q", head, 8)[0]
        else:
            raise TIFFError(f"not a TIFF file (version {version})")
        images, seen = [], set()
        while offset and offset not in seen:
            seen.add(offset)
            tags, offset = self._directory(offset)
            images.append(Image(tags, self.order, self.source))
        if not images:
            raise TIFFError("no image in the file")
        self.images = [im for im in images if not im.subfile & 4]     # transparency masks set aside
        self.masks = [im for im in images if im.subfile & 4]
        self.main = self.images[0]
        self.geokeys = geokeys(self.main.tags)
        self.transform = self._transform()

    @classmethod
    def open(cls, path: str | Path) -> "GeoTIFF":
        return cls(io.BytesIO(Path(path).read_bytes()))

    @classmethod
    def url(cls, url: str) -> "GeoTIFF":
        """A GeoTIFF on a web server that answers range requests (data/acquisition/remotezip.py)."""
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "acquisition"))
        import remotezip
        return cls(remotezip.RemoteFile(url))

    # ---- the directory ------------------------------------------------------

    def _directory(self, offset: int) -> tuple[dict, int]:
        e = self.order
        if self.big:
            n = struct.unpack(e + "Q", self.source.get(offset, 8))[0]
            size, head, link = 20, 8, "Q"
        else:
            n = struct.unpack(e + "H", self.source.get(offset, 2))[0]
            size, head, link = 12, 2, "I"
        raw = self.source.get(offset + head, n * size + (8 if self.big else 4))
        tags = {}
        for k in range(n):
            base = k * size
            tag, typ = struct.unpack_from(e + "HH", raw, base)
            count = struct.unpack_from(e + ("Q" if self.big else "I"), raw, base + 4)[0]
            inline = raw[base + (12 if self.big else 8):base + size]
            if typ not in TYPES:
                continue                     # an unknown type is skipped (TIFF 6.0, section 2)
            code, width = TYPES[typ]
            total = count * width
            data = inline[:total] if total <= len(inline) else self.source.get(struct.unpack(e + link, inline)[0], total)
            if typ == 2:
                # GeoAsciiParams is kept whole: the GeoKeys index into it byte by byte.
                tags[tag] = data.decode("latin-1") if tag == GEO_ASCII_PARAMS else data.split(b"\x00")[0].decode("latin-1")
            elif typ in (5, 10):
                pairs = struct.unpack(f"{e}{2 * count}{code}", data)
                tags[tag] = tuple(pairs[i] / pairs[i + 1] if pairs[i + 1] else math.nan for i in range(0, len(pairs), 2))
            else:
                tags[tag] = struct.unpack(f"{e}{count}{code}", data)
        following = struct.unpack_from(e + link, raw, n * size)[0]
        return tags, following

    def _transform(self) -> tuple[float, float, float, float] | None:
        tags = self.main.tags
        if MODEL_PIXEL_SCALE in tags and MODEL_TIEPOINT in tags:
            sx, sy = tags[MODEL_PIXEL_SCALE][:2]
            i, j, _, x, y, _ = tags[MODEL_TIEPOINT][:6]
            x0, y0 = x - i * sx, y + j * sy
        elif MODEL_TRANSFORMATION in tags:
            m = tags[MODEL_TRANSFORMATION]
            if m[1] or m[4]:
                raise TIFFError("a rotated model transformation is not supported")
            sx, sy, x0, y0 = m[0], -m[5], m[3], m[7]
        else:
            return None
        if self.raster_type == "point":       # the tie point is the centre of its pixel
            x0, y0 = x0 - sx / 2, y0 + sy / 2
        return (x0, sx, y0, sy)

    # ---- the coordinate reference system ---------------------------------------

    @property
    def raster_type(self) -> str:
        return "point" if self.geokeys.get(GT_RASTER_TYPE) == PIXEL_IS_POINT else "area"

    @property
    def epsg(self) -> int | None:
        """The EPSG code of the CRS, projected or geographic; None if it is user-defined or not stated."""
        model = self.geokeys.get(GT_MODEL_TYPE)
        code = self.geokeys.get(PROJECTED_CRS if model == MODEL_PROJECTED else GEODETIC_CRS)
        return code if isinstance(code, int) and code not in (0, USER_DEFINED) else None

    # ---- levels and windows -----------------------------------------------------

    def levels(self) -> list[tuple[Image, float, float]]:
        """Each image, full resolution first, with its pixel size (dx, dy) in the CRS's units."""
        x0, dx, y0, dy = self.transform
        w, h = self.main.width, self.main.height
        return [(im, dx * w / im.width, dy * h / im.height) for im in self.images]

    def level_for(self, resolution: float) -> int:
        """The coarsest level whose pixels are no larger than `resolution`, or the full image."""
        best = 0
        for k, (_, dx, dy) in enumerate(self.levels()):
            if max(dx, dy) <= resolution * (1 + 1e-9):
                best = k
        return best

    def read_box(self, x_min: float, y_min: float, x_max: float, y_max: float, level: int = 0):
        """The pixels of an image covering a box in the CRS, clipped to the image:
        (array, transform of the window). Raises TIFFError if the box misses the image."""
        image, dx, dy = self.levels()[level]
        x0, _, y0, _ = self.transform
        c0 = max(0, math.floor((x_min - x0) / dx))
        c1 = min(image.width, math.ceil((x_max - x0) / dx))
        r0 = max(0, math.floor((y0 - y_max) / dy))
        r1 = min(image.height, math.ceil((y0 - y_min) / dy))
        if c1 <= c0 or r1 <= r0:
            raise TIFFError("the box does not meet the image")
        return image.window(c0, r0, c1 - c0, r1 - r0), (x0 + c0 * dx, dx, y0 - r0 * dy, dy)


def geokeys(tags: dict) -> dict:
    """The GeoKeys of a directory (OGC 19-008r4, section 7.1): key -> value, a
    short, the doubles of GeoDoubleParams, or the text of GeoAsciiParams."""
    d = tags.get(GEO_KEY_DIRECTORY)
    if not d:
        return {}
    doubles, text = tags.get(GEO_DOUBLE_PARAMS, ()), tags.get(GEO_ASCII_PARAMS, "")
    keys = {}
    for k in range(d[3]):
        key, location, count, value = d[4 + 4 * k:8 + 4 * k]
        if location == 0:
            keys[key] = value
        elif location == GEO_DOUBLE_PARAMS:
            keys[key] = doubles[value] if count == 1 else tuple(doubles[value:value + count])
        elif location == GEO_ASCII_PARAMS:
            keys[key] = text[value:value + count].rstrip("|\x00")
        elif location == GEO_KEY_DIRECTORY:
            keys[key] = tuple(d[value:value + count])
    return keys


# ---- writing ----------------------------------------------------------------


def _entry(tag: int, typ: int, values) -> tuple[int, int, int, bytes]:
    """A directory entry as (tag, type, count, little-endian bytes of the values)."""
    if typ == 2:
        raw = values.encode("latin-1") + b"\x00"
        return tag, typ, len(raw), raw
    values = list(values) if isinstance(values, (list, tuple)) else [values]
    code = TYPES[typ][0]
    return tag, typ, len(values), struct.pack(f"<{len(values)}{code}", *values)


def _xml(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def encode_grid(values, grid: dict, *, name: str = "", units: str = "", description: str = "",
                metadata: dict | None = None) -> bytes:
    """A GeoTIFF of one float32 grid (nlat, nlon) on the library's grid
    {lon0, dlon, nlon, lat0, dlat, nlat}: values at the centres of the pixels,
    rows north to south, NaN where there is no value. `metadata` adds items
    (name -> text) to GDAL_METADATA, such as the time."""
    import numpy as np
    a = np.ascontiguousarray(np.asarray(values, dtype=np.float32))
    nlat, nlon = a.shape
    if (nlat, nlon) != (grid["nlat"], grid["nlon"]):
        raise TIFFError(f"values {a.shape} do not match the grid {grid['nlat']} x {grid['nlon']}")
    rows_per_strip = max(1, min(nlat, (1 << 16) // (4 * nlon)))
    strips = []
    for r in range(0, nlat, rows_per_strip):
        block = a[r:r + rows_per_strip][:, :, None]
        strips.append(zlib.compress(floating_point_difference(block), 6))

    # The model space: the outer corner of the first pixel, half a cell west and north of its centre.
    corner_lon = grid["lon0"] - grid["dlon"] / 2
    corner_lat = grid["lat0"] + grid["dlat"] / 2
    citation = "WGS 84|"
    keys = [(GT_MODEL_TYPE, 0, 1, MODEL_GEOGRAPHIC), (GT_RASTER_TYPE, 0, 1, PIXEL_IS_AREA),
            (GEODETIC_CRS, 0, 1, 4326), (GEODETIC_CITATION, GEO_ASCII_PARAMS, len(citation), 0),
            (GEOG_ANGULAR_UNITS, 0, 1, 9102)]
    directory = [1, 1, 1, len(keys)] + [x for k in keys for x in k]       # GeoTIFF 1.1: revision 1.1
    items = [f'  <Item name="DESCRIPTION" sample="0" role="description">{_xml(name)}</Item>'] if name else []
    if units:
        items.append(f'  <Item name="UNITTYPE" sample="0" role="unittype">{_xml(units)}</Item>')
    for key, text in (metadata or {}).items():
        items.append(f'  <Item name="{_xml(key)}">{_xml(str(text))}</Item>')
    entries = [
        _entry(IMAGE_WIDTH, 4, nlon), _entry(IMAGE_LENGTH, 4, nlat), _entry(BITS_PER_SAMPLE, 3, 32),
        _entry(COMPRESSION, 3, DEFLATE), _entry(PHOTOMETRIC, 3, 1),
        _entry(STRIP_OFFSETS, 4, [0] * len(strips)), _entry(SAMPLES_PER_PIXEL, 3, 1),
        _entry(ROWS_PER_STRIP, 4, rows_per_strip), _entry(STRIP_BYTE_COUNTS, 4, [len(s) for s in strips]),
        _entry(PLANAR_CONFIGURATION, 3, 1), _entry(SOFTWARE, 2, "The Science Infrastructures, data/formats/tiff.py"),
        _entry(PREDICTOR, 3, FLOATING_POINT), _entry(SAMPLE_FORMAT, 3, 3),
        _entry(MODEL_PIXEL_SCALE, 12, [grid["dlon"], grid["dlat"], 0.0]),
        _entry(MODEL_TIEPOINT, 12, [0.0, 0.0, 0.0, corner_lon, corner_lat, 0.0]),
        _entry(GEO_KEY_DIRECTORY, 3, directory), _entry(GEO_ASCII_PARAMS, 2, citation),
        _entry(GDAL_NODATA, 2, "nan"),
    ]
    if description:
        entries.append(_entry(IMAGE_DESCRIPTION, 2, description))
    if items:
        entries.append(_entry(GDAL_METADATA, 2, "<GDALMetadata>\n" + "\n".join(items) + "\n</GDALMetadata>"))
    entries.sort(key=lambda t: t[0])             # entries in ascending order of tag (TIFF 6.0, section 2)

    # Layout: header, directory, values too long for their entries (on word boundaries), strips.
    extra_at = 8 + 2 + 12 * len(entries) + 4
    at = extra_at + sum(len(raw) + len(raw) % 2 for *_, raw in entries if len(raw) > 4)
    strip_offsets = []
    for s in strips:
        strip_offsets.append(at)
        at += len(s)
    entries = [_entry(STRIP_OFFSETS, 4, strip_offsets) if tag == STRIP_OFFSETS else (tag, typ, count, raw)
               for tag, typ, count, raw in entries]
    out = bytearray(b"II" + struct.pack("<HI", 42, 8))
    out += struct.pack("<H", len(entries))
    tail = bytearray()
    for tag, typ, count, raw in entries:
        if len(raw) > 4:
            out += struct.pack("<HHII", tag, typ, count, extra_at + len(tail))
            tail += raw + b"\x00" * (len(raw) % 2)
        else:
            out += struct.pack("<HHI", tag, typ, count) + raw + b"\x00" * (4 - len(raw))
    out += struct.pack("<I", 0)
    out += tail
    if len(out) != strip_offsets[0]:
        raise TIFFError("internal error: the strips are not where the directory says")
    for s in strips:
        out += s
    return bytes(out)


def write_grid(path: str | Path, values, grid: dict, **options) -> None:
    """Writes encode_grid(values, grid, **options) to a file, through a .part file."""
    path = Path(path)
    part = path.with_name(path.name + ".part")
    part.write_bytes(encode_grid(values, grid, **options))
    part.replace(path)


def grid_of(tif: GeoTIFF) -> dict:
    """The library's grid {lon0, dlon, nlon, lat0, dlat, nlat} of a GeoTIFF in EPSG:4326."""
    if tif.epsg != 4326:
        raise TIFFError(f"the file is in EPSG:{tif.epsg}, not EPSG:4326")
    x0, dx, y0, dy = tif.transform
    return {"lon0": x0 + dx / 2, "dlon": dx, "nlon": tif.main.width,
            "lat0": y0 - dy / 2, "dlat": dy, "nlat": tif.main.height}


# ---- the transverse Mercator projection: UTM ---------------------------------

WGS84_A = 6378137.0                      # m, semi-major axis (NIMA TR8350.2)
WGS84_F = 1 / 298.257223563              # flattening
UTM_K0 = 0.9996                          # scale on the central meridian
UTM_FALSE_EASTING = 500000.0             # m
UTM_FALSE_NORTHING_SOUTH = 10000000.0    # m, in the southern zones


class TransverseMercator:
    """The transverse Mercator projection of an ellipsoid by Krüger's series to
    order n^6 in the third flattening n = f / (2 - f) (Karney 2011, section 2,
    with the coefficients of eqs. 29, 35 and 36). Longitudes and latitudes in
    degrees, easting and northing in metres. Error below 5 nm within 3900 km of
    the central meridian (Karney 2011, section 4); the series diverge beyond
    about 8000 km, so use it within a few thousand km of the central meridian."""

    def __init__(self, central_meridian: float = 0.0, k0: float = UTM_K0, false_easting: float = 0.0,
                 false_northing: float = 0.0, a: float = WGS84_A, f: float = WGS84_F):
        self.lon0, self.k0, self.fe, self.fn = central_meridian, k0, false_easting, false_northing
        n = f / (2 - f)
        self.e2 = f * (2 - f)
        self.e = math.sqrt(self.e2)
        self.A = a / (1 + n) * (1 + n**2 / 4 + n**4 / 64 + n**6 / 256)                       # eq. (29)
        self.alpha = [                                                                         # eq. (35)
            n / 2 - 2 * n**2 / 3 + 5 * n**3 / 16 + 41 * n**4 / 180 - 127 * n**5 / 288 + 7891 * n**6 / 37800,
            13 * n**2 / 48 - 3 * n**3 / 5 + 557 * n**4 / 1440 + 281 * n**5 / 630 - 1983433 * n**6 / 1935360,
            61 * n**3 / 240 - 103 * n**4 / 140 + 15061 * n**5 / 26880 + 167603 * n**6 / 181440,
            49561 * n**4 / 161280 - 179 * n**5 / 168 + 6601661 * n**6 / 7257600,
            34729 * n**5 / 80640 - 3418889 * n**6 / 1995840,
            212378941 * n**6 / 319334400,
        ]
        self.beta = [                                                                          # eq. (36)
            n / 2 - 2 * n**2 / 3 + 37 * n**3 / 96 - n**4 / 360 - 81 * n**5 / 512 + 96199 * n**6 / 604800,
            n**2 / 48 + n**3 / 15 - 437 * n**4 / 1440 + 46 * n**5 / 105 - 1118711 * n**6 / 3870720,
            17 * n**3 / 480 - 37 * n**4 / 840 - 209 * n**5 / 4480 + 5569 * n**6 / 90720,
            4397 * n**4 / 161280 - 11 * n**5 / 504 - 830251 * n**6 / 7257600,
            4583 * n**5 / 161280 - 108847 * n**6 / 3991680,
            20648693 * n**6 / 638668800,
        ]

    def _conformal(self, tau):
        """tau' = tan of the conformal latitude, from tau = tan of the geodetic latitude: eqs. (7)-(9)."""
        import numpy as np
        sigma = np.sinh(self.e * np.arctanh(self.e * tau / np.sqrt(1 + tau**2)))
        return tau * np.sqrt(1 + sigma**2) - sigma * np.sqrt(1 + tau**2)

    def forward(self, lon, lat):
        """(easting, northing) in metres of geodetic longitude and latitude in degrees."""
        import numpy as np
        lam = np.radians(np.asarray(lon, dtype=np.float64) - self.lon0)
        lam = (lam + np.pi) % (2 * np.pi) - np.pi
        phi = np.radians(np.asarray(lat, dtype=np.float64))
        tau_p = self._conformal(np.tan(phi))
        xi_p = np.arctan2(tau_p, np.cos(lam))                                    # eq. (10)
        eta_p = np.arcsinh(np.sin(lam) / np.sqrt(tau_p**2 + np.cos(lam)**2))
        xi, eta = xi_p.copy(), eta_p.copy()
        for j, a in enumerate(self.alpha, start=1):                               # eq. (11)
            xi += a * np.sin(2 * j * xi_p) * np.cosh(2 * j * eta_p)
            eta += a * np.cos(2 * j * xi_p) * np.sinh(2 * j * eta_p)
        return self.fe + self.k0 * self.A * eta, self.fn + self.k0 * self.A * xi    # eq. (13)

    def inverse(self, x, y):
        """(longitude, latitude) in degrees of easting and northing in metres."""
        import numpy as np
        eta = (np.asarray(x, dtype=np.float64) - self.fe) / (self.k0 * self.A)      # eq. (15)
        xi = (np.asarray(y, dtype=np.float64) - self.fn) / (self.k0 * self.A)
        xi_p, eta_p = xi.copy(), eta.copy()
        for j, b in enumerate(self.beta, start=1):                                 # eq. (16)
            xi_p -= b * np.sin(2 * j * xi) * np.cosh(2 * j * eta)
            eta_p -= b * np.cos(2 * j * xi) * np.sinh(2 * j * eta)
        tau_p = np.sin(xi_p) / np.sqrt(np.sinh(eta_p)**2 + np.cos(xi_p)**2)           # eq. (18)
        lam = np.arctan2(np.sinh(eta_p), np.cos(xi_p))
        tau = tau_p.copy()                                                         # eqs. (19)-(21), Newton's method
        for _ in range(4):
            t_i = self._conformal(tau)
            tau = tau + (tau_p - t_i) / np.sqrt(1 + t_i**2) * (1 + (1 - self.e2) * tau**2) / ((1 - self.e2) * np.sqrt(1 + tau**2))
        return self.lon0 + np.degrees(lam), np.degrees(np.arctan(tau))              # eq. (22)


def utm(epsg: int) -> TransverseMercator:
    """The projection of a UTM zone on WGS 84 by its EPSG code: 326zz north, 327zz south."""
    zone, south = epsg % 100, epsg // 100 == 327
    if epsg // 100 not in (326, 327) or not 1 <= zone <= 60:
        raise TIFFError(f"EPSG:{epsg} is not a UTM zone on WGS 84")
    return TransverseMercator(central_meridian=6 * zone - 183, k0=UTM_K0, false_easting=UTM_FALSE_EASTING,
                              false_northing=UTM_FALSE_NORTHING_SOUTH if south else 0.0)
