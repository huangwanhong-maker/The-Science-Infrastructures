"""Members of a zip archive on a web server, read without downloading the archive.

A zip archive ends with a central directory that lists its members and where
each one starts. Python's zipfile reads an archive through any seekable file;
given one whose reads are HTTP range requests (RFC 9110, section 14), it reads
the directory and then only the bytes of the members asked for. A 130 MB
archive of maps thus yields one 200 kB map in three or four requests.

The server must answer range requests (status 206). Redirects are followed
once, and the final address is kept, since signed addresses (as from figshare
to Amazon S3) are valid for a few seconds only for the first request; the
address is renewed if the server refuses a later one.
"""

from __future__ import annotations

import io
import urllib.error
import urllib.request
import zipfile

AGENT = "laboratory-data-fetch/1 (+https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory)"


class RangeError(OSError):
    pass


class RemoteFile(io.RawIOBase):
    """A read-only, seekable file whose reads are HTTP range requests."""

    def __init__(self, url: str, timeout: int = 120):
        self.url, self.timeout = url, timeout
        self.requests = 0
        self.pos = 0
        self._locate()

    def _open(self, start: int, end: int):
        request = urllib.request.Request(self.final, headers={"Range": f"bytes={start}-{end}", "User-Agent": AGENT})
        self.requests += 1
        return urllib.request.urlopen(request, timeout=self.timeout)

    def _locate(self) -> None:
        """Follows redirects to the final address and learns the size."""
        request = urllib.request.Request(self.url, headers={"Range": "bytes=0-0", "User-Agent": AGENT})
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            if response.status != 206:
                raise RangeError(f"{self.url}: the server does not answer range requests")
            self.final = response.url
            self.size = int(response.headers["Content-Range"].rsplit("/", 1)[1])

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self.pos, io.SEEK_END: self.size}[whence]
        self.pos = max(0, base + offset)
        return self.pos

    def read(self, n: int = -1) -> bytes:
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        try:
            response = self._open(self.pos, self.pos + n - 1)
        except urllib.error.HTTPError as error:
            if error.code not in (400, 403):
                raise
            self._locate()                      # a signed address has expired
            response = self._open(self.pos, self.pos + n - 1)
        with response:
            data = response.read()
        self.pos += len(data)
        return data

    def readinto(self, buffer) -> int:
        data = self.read(len(buffer))
        buffer[:len(data)] = data
        return len(data)


def members(url: str) -> list[zipfile.ZipInfo]:
    """The members of a remote archive."""
    return zipfile.ZipFile(RemoteFile(url)).infolist()


def read_members(url: str, names: list[str]) -> dict[str, bytes]:
    """The contents of the named members of a remote archive (checked by their CRC)."""
    archive = zipfile.ZipFile(RemoteFile(url))
    return {name: archive.read(name) for name in names}
