"""Members of a tar archive on a web server, read without downloading the archive.

A tar archive is a sequence of members, each a 512-byte header (its name, its size in octal, its type)
followed by its bytes padded to a multiple of 512. There is no directory at the end, as a zip has; the
members are found by walking the headers, one range request each (RFC 9110, section 14), from the first
to the one sought. The walk is kept as an index beside the cache (`index_path`), so that it is made once:
a 9 GB archive of 900 scans yields one 30 MB scan in a few hundred small requests the first time, and in
one request afterwards.

Names: the ustar prefix field is joined to the name; GNU long names (type 'L') and pax extended headers
(type 'x', their 'path' record) are read. Members of other types than regular files are listed but not
read. The server must answer range requests (remotezip.RemoteFile).
"""

from __future__ import annotations

import json
from pathlib import Path

from remotezip import RemoteFile

BLOCK = 512


def _octal(field: bytes) -> int:
    text = field.rstrip(b"\0 ").strip()
    return int(text, 8) if text else 0


def _pax_path(data: bytes) -> str | None:
    """The 'path' record of a pax extended header ("<length> path=<value>\\n" records)."""
    for line in data.split(b"\n"):
        if b" path=" in line:
            return line.split(b" path=", 1)[1].decode("utf-8", "replace")
    return None


class RemoteTar:
    """The index of a remote tar archive, walked as far as needed and kept in a JSON file."""

    def __init__(self, url: str, index_path: Path | None = None):
        self.url = url
        self.file = RemoteFile(url)
        self.index_path = index_path
        self.index = {"url": url, "size": self.file.size, "next": 0, "complete": False, "members": {}}
        if index_path and index_path.exists():
            kept = json.loads(index_path.read_text(encoding="utf-8"))
            if kept.get("url") == url and kept.get("size") == self.file.size:
                self.index = kept

    def _save(self) -> None:
        if self.index_path:
            self.index_path.parent.mkdir(parents=True, exist_ok=True)
            self.index_path.write_text(json.dumps(self.index, indent=1) + "\n", encoding="utf-8")

    def _read(self, offset: int, n: int) -> bytes:
        self.file.seek(offset)
        return self.file.read(n)

    def walk(self, until: str | None = None, report=None) -> None:
        """Walks the headers from where the index stopped, until `until` is found or the archive ends."""
        offset, long_name = self.index["next"], None
        steps = 0
        while not self.index["complete"]:
            if until is not None and until in self.index["members"]:
                break
            header = self._read(offset, BLOCK)
            steps += 1
            if len(header) < BLOCK or header == b"\0" * BLOCK:
                self.index["complete"] = True
                break
            name = header[0:100].rstrip(b"\0").decode("utf-8", "replace")
            prefix = header[345:500].rstrip(b"\0").decode("utf-8", "replace") if header[257:262] == b"ustar" else ""
            if prefix:
                name = f"{prefix}/{name}"
            size = _octal(header[124:136])
            kind = header[156:157].decode("ascii", "replace") or "0"
            data_offset = offset + BLOCK
            offset = data_offset + ((size + BLOCK - 1) // BLOCK) * BLOCK
            if kind == "L":                               # GNU long name for the next member
                long_name = self._read(data_offset, size).rstrip(b"\0").decode("utf-8", "replace")
                continue
            if kind == "x":                               # pax extended header for the next member
                long_name = _pax_path(self._read(data_offset, size)) or long_name
                continue
            if long_name:
                name, long_name = long_name, None
            self.index["members"][name] = {"offset": data_offset, "size": size, "type": kind}
            self.index["next"] = offset
            if report and steps % 50 == 0:
                report(f"  indexing the archive: {len(self.index['members'])} members")
            if steps % 200 == 0:
                self._save()
        self._save()

    def read(self, name: str, report=None) -> bytes:
        """The bytes of a member, walking the headers to it if the index does not yet reach it."""
        if name not in self.index["members"]:
            self.walk(until=name, report=report)
        entry = self.index["members"].get(name)
        if entry is None:
            raise KeyError(f"{self.url}: no member {name}")
        if entry["type"] not in ("0", "\0", "7"):
            raise ValueError(f"{self.url}: {name} is not a regular file (type {entry['type']})")
        return self._read(entry["offset"], entry["size"])

    def names(self) -> list[str]:
        self.walk(report=print)
        return list(self.index["members"])
