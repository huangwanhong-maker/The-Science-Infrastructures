"""A reader and a writer for netCDF classic files (CDF-1 and CDF-2).

ERDDAP and many other services deliver gridded data in this format. The reader
parses the header (dimensions, attributes, variables) and reads variables as
Python arrays, unpacking scale_factor and add_offset and turning fill values
into NaN for floating-point results; it needs only the standard library.

The writer stores fixed-size variables (no record dimension) in the 64-bit
offset format, so that data fetched in other forms, such as OPeNDAP responses,
are cached in a standard, self-describing file. It needs numpy.

Reference: Unidata, "NetCDF Classic and 64-bit Offset File Formats"
(https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html).
All values are big-endian.
"""

from __future__ import annotations

import math
import struct
import sys
from array import array
from pathlib import Path

NC_DIMENSION, NC_VARIABLE, NC_ATTRIBUTE = 0x0A, 0x0B, 0x0C
# nc_type -> (struct code, size, array typecode)
TYPES = {1: ("b", 1, "b"), 2: ("c", 1, "B"), 3: ("h", 2, "h"), 4: ("i", 4, "i"), 5: ("f", 4, "f"), 6: ("d", 8, "d")}
STREAMING = 0xFFFFFFFF


class NetCDFError(ValueError):
    pass


def _pad(n: int) -> int:
    return (4 - n % 4) % 4


class Variable:
    def __init__(self, name, dimensions, attributes, nc_type, vsize, begin, record):
        self.name, self.dimensions, self.attributes = name, dimensions, attributes
        self.nc_type, self.vsize, self.begin, self.record = nc_type, vsize, begin, record

    def __repr__(self):
        return f"Variable({self.name}, {self.dimensions})"


class Dataset:
    """A netCDF classic dataset held in memory."""

    def __init__(self, data: bytes):
        self.data = data
        if data[:3] != b"CDF" or data[3] not in (1, 2):
            raise NetCDFError("not a netCDF classic (CDF-1 or CDF-2) file")
        self.offset_size = 4 if data[3] == 1 else 8
        self.pos = 4
        self.numrecs = self._u32()
        self.dimensions = {}            # name -> length (0 for the record dimension)
        self._dim_names = []
        self.attributes = {}
        self.variables = {}
        self._read_dimensions()
        self.attributes = self._read_attributes()
        self._read_variables()
        if self.numrecs == STREAMING:
            self.numrecs = self._infer_records()

    @classmethod
    def open(cls, path: str | Path) -> "Dataset":
        return cls(Path(path).read_bytes())

    # ---- header -------------------------------------------------------------

    def _u32(self) -> int:
        value = struct.unpack_from(">I", self.data, self.pos)[0]
        self.pos += 4
        return value

    def _name(self) -> str:
        n = self._u32()
        text = self.data[self.pos:self.pos + n].decode("utf-8")
        self.pos += n + _pad(n)
        return text

    def _values(self, nc_type: int, n: int):
        code, size, _ = TYPES[nc_type]
        raw = self.data[self.pos:self.pos + n * size]
        self.pos += n * size + _pad(n * size)
        if nc_type == 2:
            return raw.rstrip(b"\x00").decode("utf-8", "replace")
        values = struct.unpack(f">{n}{code}", raw)
        return values[0] if n == 1 else list(values)

    def _list_header(self, tag: int) -> int:
        kind, n = self._u32(), self._u32()
        if kind == 0 and n == 0:
            return 0
        if kind != tag:
            raise NetCDFError(f"expected tag {tag:#x}, found {kind:#x}")
        return n

    def _read_dimensions(self) -> None:
        for _ in range(self._list_header(NC_DIMENSION)):
            name = self._name()
            length = self._u32()
            self.dimensions[name] = length
            self._dim_names.append(name)

    def _read_attributes(self) -> dict:
        attributes = {}
        for _ in range(self._list_header(NC_ATTRIBUTE)):
            name = self._name()
            nc_type, n = self._u32(), self._u32()
            attributes[name] = self._values(nc_type, n)
        return attributes

    def _read_variables(self) -> None:
        for _ in range(self._list_header(NC_VARIABLE)):
            name = self._name()
            dims = [self._dim_names[self._u32()] for _ in range(self._u32())]
            attributes = self._read_attributes()
            nc_type, vsize = self._u32(), self._u32()
            begin = struct.unpack_from(">I" if self.offset_size == 4 else ">Q", self.data, self.pos)[0]
            self.pos += self.offset_size
            record = bool(dims) and self.dimensions[dims[0]] == 0
            self.variables[name] = Variable(name, dims, attributes, nc_type, vsize, begin, record)

    def _record_size(self) -> int:
        records = [v for v in self.variables.values() if v.record]
        if len(records) == 1:            # a single record variable is not padded
            v = records[0]
            return TYPES[v.nc_type][1] * self._slab(v)
        return sum(v.vsize for v in records)

    def _infer_records(self) -> int:
        records = [v for v in self.variables.values() if v.record]
        if not records:
            return 0
        first = min(v.begin for v in records)
        return (len(self.data) - first) // max(1, self._record_size())

    # ---- data ---------------------------------------------------------------

    def shape(self, name: str) -> tuple[int, ...]:
        v = self.variables[name]
        return tuple(self.numrecs if self.dimensions[d] == 0 else self.dimensions[d] for d in v.dimensions)

    def _slab(self, v: Variable) -> int:
        n = 1
        for d in v.dimensions[1 if v.record else 0:]:
            n *= self.dimensions[d]
        return n

    def raw(self, name: str) -> array:
        """The stored values of a variable, in C order, as a native-endian array."""
        v = self.variables[name]
        code, size, typecode = TYPES[v.nc_type]
        out = array(typecode)
        if not v.record:
            n = 1
            for s in self.shape(name):
                n *= s
            out.frombytes(self.data[v.begin:v.begin + n * size])
        else:
            slab = self._slab(v) * size
            stride = self._record_size()
            for r in range(self.numrecs):
                start = v.begin + r * stride
                out.frombytes(self.data[start:start + slab])
        if sys.byteorder == "little" and size > 1:
            out.byteswap()
        return out

    def array(self, name: str):
        """Values as a numpy float64 array of the variable's shape, unpacked as by values().

        Needs numpy; values() gives the same numbers without it, more slowly.
        """
        import numpy as np
        v = self.variables[name]
        code = TYPES[v.nc_type][0]
        dtype = np.dtype(">" + {"b": "i1", "c": "u1", "h": "i2", "i": "i4", "f": "f4", "d": "f8"}[code])
        shape = self.shape(name)
        if not v.record:
            stored = np.frombuffer(self.data, dtype, int(np.prod(shape)), v.begin)
        else:
            slab, stride = self._slab(v), self._record_size()
            stored = np.concatenate([np.frombuffer(self.data, dtype, slab, v.begin + r * stride)
                                     for r in range(self.numrecs)])
        stored = stored.reshape(shape)
        a = v.attributes
        out = stored.astype(np.float64)
        for key in ("_FillValue", "missing_value"):
            for fill in (a[key] if isinstance(a.get(key), list) else [a.get(key)]):
                if fill is not None and not (isinstance(fill, float) and math.isnan(fill)):
                    out[stored == fill] = np.nan
        return out * a.get("scale_factor", 1.0) + a.get("add_offset", 0.0)

    def values(self, name: str) -> array:
        """Values with scale_factor and add_offset applied and fill values set to NaN, as doubles."""
        v = self.variables[name]
        raw = self.raw(name)
        a = v.attributes
        scale = a.get("scale_factor", 1.0)
        offset = a.get("add_offset", 0.0)
        fills = {x for key in ("_FillValue", "missing_value") if key in a
                 for x in (a[key] if isinstance(a[key], list) else [a[key]])}
        out = array("d", bytes(8 * len(raw)))
        for i, x in enumerate(raw):
            out[i] = math.nan if (x in fills or (isinstance(x, float) and math.isnan(x))) else x * scale + offset
        return out


# ---- writing ----------------------------------------------------------------

NC_TYPE_OF = {"int8": 1, "int16": 3, "int32": 4, "float32": 5, "float64": 6}


def _name_bytes(name: str) -> bytes:
    raw = name.encode("utf-8")
    return struct.pack(">I", len(raw)) + raw + b"\x00" * _pad(len(raw))


def _typed(value):
    """An attribute value as (nc_type, big-endian bytes, count)."""
    import numpy as np
    if isinstance(value, str):
        raw = value.encode("utf-8")
        return 2, raw, len(raw)
    if isinstance(value, (bool, int)):
        a = np.asarray([value], ">i4")
    elif isinstance(value, float):
        a = np.asarray([value], ">f8")
    elif isinstance(value, (list, tuple)):
        a = np.asarray(value)
        a = a.astype(">i4" if a.dtype.kind in "iu" else ">f8")
    else:                                   # a numpy scalar or array keeps its type
        a = np.atleast_1d(np.asarray(value))
        a = a.astype(a.dtype.newbyteorder(">"))
    return NC_TYPE_OF[a.dtype.newbyteorder("=").name], a.tobytes(), a.size


def _attributes_bytes(attributes: dict) -> bytes:
    if not attributes:
        return struct.pack(">II", 0, 0)
    out = struct.pack(">II", NC_ATTRIBUTE, len(attributes))
    for name, value in attributes.items():
        nc_type, raw, count = _typed(value)
        out += _name_bytes(name) + struct.pack(">II", nc_type, count) + raw + b"\x00" * _pad(len(raw))
    return out


def write(path: str | Path, dimensions: dict, variables: dict, attributes: dict | None = None) -> None:
    """Writes a netCDF classic file (64-bit offset format) of fixed-size variables.

    dimensions: name -> length. variables: name -> (dimension names, numpy
    array of int8, int16, int32, float32 or float64, attributes). An attribute
    is a str (written as text), a Python int or float (as int32 or double), a
    list of numbers, or a numpy scalar or array (in its own type, as CF
    requires of _FillValue).
    """
    import numpy as np
    dim_names = list(dimensions)
    arrays = {}
    for name, (dims, data, _) in variables.items():
        data = np.ascontiguousarray(data)
        if data.dtype.name not in NC_TYPE_OF:
            raise NetCDFError(f"{name}: cannot write values of type {data.dtype}")
        if tuple(dimensions[d] for d in dims) != data.shape:
            raise NetCDFError(f"{name}: shape {data.shape} does not match dimensions {dims}")
        arrays[name] = data.astype(data.dtype.newbyteorder(">"))

    def header(begins: list[int]) -> bytes:
        out = b"CDF\x02" + struct.pack(">I", 0)
        out += struct.pack(">II", NC_DIMENSION, len(dim_names)) if dim_names else struct.pack(">II", 0, 0)
        for d in dim_names:
            out += _name_bytes(d) + struct.pack(">I", dimensions[d])
        out += _attributes_bytes(attributes or {})
        out += struct.pack(">II", NC_VARIABLE, len(variables)) if variables else struct.pack(">II", 0, 0)
        for (name, (dims, _, attrs)), begin in zip(variables.items(), begins):
            a = arrays[name]
            out += _name_bytes(name) + struct.pack(">I", len(dims))
            out += b"".join(struct.pack(">I", dim_names.index(d)) for d in dims)
            out += _attributes_bytes(attrs)
            vsize = a.nbytes + _pad(a.nbytes)
            out += struct.pack(">II", NC_TYPE_OF[a.dtype.newbyteorder("=").name], min(vsize, 0xFFFFFFFF))
            out += struct.pack(">Q", begin)
        return out

    size = len(header([0] * len(variables)))
    begins = []
    for name in variables:
        begins.append(size)
        size += arrays[name].nbytes + _pad(arrays[name].nbytes)
    path = Path(path)
    part = path.with_name(path.name + ".part")
    with part.open("wb") as out:
        out.write(header(begins))
        for name in variables:
            raw = arrays[name].tobytes()
            out.write(raw + b"\x00" * _pad(len(raw)))
    part.replace(path)
