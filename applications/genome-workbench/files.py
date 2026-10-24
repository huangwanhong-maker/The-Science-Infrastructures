"""The workbench's library of the user's files: uploaded through the page, kept in the data cache,
recognised by their content, and named in requests as ``file:<id>``.

A file is stored once under its content's hash, so the same file added twice is one record, in
``data/cache/genome-workbench/files/<id>/<name>``, beside ``index.json``, the records of all of
them. Gzipped files are kept as they come and also unpacked, since the engine reads plain text.
Each record says what the file is (its kind, recognised from its first lines rather than its name)
and what is in it, counted here cheaply: sequences and bases for FASTA, reads for FASTQ, records
for GenBank, samples and sites for VCF, leaves for Newick, records for SAM.

The library is an edge (CLAUDE.md §3.3): it touches the user's files; the engine reads only the
paths a request gives it.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STORE = ROOT / "data" / "cache" / "genome-workbench" / "files"
INDEX = STORE / "index.json"
LOCK = threading.Lock()

# What each kind is, and which kinds can stand for a genome.
KINDS = {
    "fasta": "sequences (FASTA)",
    "fastq": "reads (FASTQ)",
    "genbank": "an annotated genome (GenBank)",
    "vcf": "variants (VCF)",
    "newick": "a tree (Newick)",
    "sam": "alignments (SAM)",
    "text": "text",
}
GENOME_KINDS = {"fasta", "genbank"}


def _load() -> dict:
    try:
        return json.loads(INDEX.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"files": {}}


def _save(index: dict) -> None:
    STORE.mkdir(parents=True, exist_ok=True)
    tmp = INDEX.with_suffix(".tmp")
    tmp.write_text(json.dumps(index, indent=1, ensure_ascii=False), encoding="utf-8")
    tmp.replace(INDEX)


def safe_name(name: str) -> str:
    """A file name kept to letters, digits, dots, dashes and underscores."""
    base = Path(name or "file").name
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", base).strip("._") or "file"
    return base[:120]


def sniff(head: str) -> str:
    """The kind of a text from its first lines."""
    lines = [l for l in head.splitlines() if l.strip()][:8]
    if not lines:
        return "text"
    first = lines[0].lstrip()
    if first.startswith("##fileformat=VCF"):
        return "vcf"
    if first.startswith("LOCUS"):
        return "genbank"
    if first.startswith(">"):
        return "fasta"
    if first.startswith("@") and len(lines) >= 3 and lines[2].startswith("+"):
        return "fastq"
    if first.startswith(("@HD", "@SQ", "@PG", "@RG")) or (len(first.split("\t")) >= 11 and first.split("\t")[1].isdigit()):
        return "sam"
    if first.startswith("(") and ";" in head:
        return "newick"
    if re.fullmatch(r"[ACGTUNacgtun\s]+", "\n".join(lines)):
        return "fasta"  # a bare sequence
    return "text"


def summarise(path: Path, kind: str, limit: int = 400 * 1024 * 1024) -> dict:
    """What is in a file, counted in one pass (files beyond `limit` bytes are counted in part)."""
    out: dict = {}
    seen = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            if kind == "fasta":
                n = bases = gc = 0
                lengths: list[int] = []
                names: list[str] = []
                cur = 0
                for line in f:
                    seen += len(line)
                    if line.startswith(">"):
                        if n:
                            lengths.append(cur)
                        n += 1
                        cur = 0
                        if len(names) < 50:
                            names.append(line[1:].strip().split(" ")[0][:60])
                    else:
                        s = line.strip()
                        cur += len(s)
                        bases += len(s)
                        gc += sum(s.count(c) for c in "GCgc")
                    if seen > limit:
                        out["partial"] = True
                        break
                lengths.append(cur) if (n or cur) else None
                n = max(n, 1 if cur else 0)
                lengths.sort(reverse=True)
                half, acc, n50 = bases / 2, 0, 0
                for L in lengths:
                    acc += L
                    if acc >= half:
                        n50 = L
                        break
                out.update(sequences=n, bases=bases, gc=round(gc / bases, 4) if bases else None, longest=lengths[0] if lengths else 0,
                           n50=n50, names=names)
            elif kind == "fastq":
                reads = bases = 0
                for k, line in enumerate(f):
                    seen += len(line)
                    if k % 4 == 1:
                        reads += 1
                        bases += len(line.strip())
                    if seen > limit:
                        out["partial"] = True
                        break
                out.update(reads=reads, bases=bases, mean_length=round(bases / reads, 1) if reads else 0)
            elif kind == "genbank":
                records = 0
                organism = definition = None
                for line in f:
                    if line.startswith("LOCUS"):
                        records += 1
                    elif definition is None and line.startswith("DEFINITION"):
                        definition = line[10:].strip()
                    elif organism is None and line.strip().startswith("ORGANISM"):
                        organism = line.strip()[8:].strip()
                out.update(records=records, definition=definition, organism=organism)
            elif kind == "vcf":
                sites = 0
                samples: list[str] = []
                for line in f:
                    if line.startswith("#CHROM"):
                        samples = line.rstrip("\n").split("\t")[9:]
                    elif not line.startswith("#") and line.strip():
                        sites += 1
                out.update(samples=len(samples), sites=sites, sample_names=samples[:50])
            elif kind == "newick":
                text = f.read(4 * 1024 * 1024)
                out.update(trees=text.count(";"), leaves=len(re.findall(r"[(,]\s*([^(),:;\[\]]+)", text)))
            elif kind == "sam":
                records = 0
                for line in f:
                    if not line.startswith("@") and line.strip():
                        records += 1
                out.update(records=records)
    except OSError as error:
        out["error"] = str(error)
    return out


def add(name: str, stream, length: int, max_bytes: int) -> dict:
    """Stores an uploaded file (read from `stream`, `length` bytes) and returns its record."""
    if length > max_bytes:
        raise ValueError(f"a file of {length / 1e6:.0f} MB is larger than this server takes ({max_bytes / 1e6:.0f} MB)")
    STORE.mkdir(parents=True, exist_ok=True)
    tmp = STORE / f"upload-{threading.get_ident()}-{time.time_ns()}.part"
    h = hashlib.sha256()
    left = length
    with tmp.open("wb") as out:
        while left > 0:
            chunk = stream.read(min(1 << 20, left))
            if not chunk:
                break
            out.write(chunk)
            h.update(chunk)
            left -= len(chunk)
    if left:
        tmp.unlink(missing_ok=True)
        raise ValueError("the upload ended early")
    fid = h.hexdigest()[:16]
    name = safe_name(name)
    folder = STORE / fid
    folder.mkdir(parents=True, exist_ok=True)
    stored = folder / name
    tmp.replace(stored)
    # Gzipped: unpack beside it; the engine reads the plain text.
    with stored.open("rb") as fb:
        magic = fb.read(2)
    plain = stored
    if magic == b"\x1f\x8b":
        plain = folder / (name[:-3] if name.endswith(".gz") else name + ".txt")
        with gzip.open(stored, "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
    with plain.open("r", encoding="utf-8", errors="replace") as f:
        kind = sniff(f.read(64 * 1024))
    record = {"id": f"file:{fid}", "name": name, "path": str(plain.relative_to(ROOT)).replace("\\", "/"), "bytes": plain.stat().st_size,
              "kind": kind, "about": KINDS[kind], "genome": kind in GENOME_KINDS, "added": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    record["summary"] = summarise(plain, kind)
    with LOCK:
        index = _load()
        index["files"][fid] = record
        _save(index)
    return record


def records() -> list[dict]:
    with LOCK:
        return sorted(_load()["files"].values(), key=lambda r: r.get("added", ""), reverse=True)


def get(fid: str) -> dict | None:
    fid = fid.removeprefix("file:")
    with LOCK:
        return _load()["files"].get(fid)


def path_of(fid: str) -> Path:
    """The plain file of a record, or ValueError."""
    r = get(fid)
    if r is None:
        raise ValueError(f"no file {fid} in the library")
    p = ROOT / r["path"]
    if not p.exists():
        raise ValueError(f"the file {r['name']} is no longer in the cache")
    return p


def remove(fid: str) -> bool:
    fid = fid.removeprefix("file:")
    with LOCK:
        index = _load()
        r = index["files"].pop(fid, None)
        _save(index)
    if r is None:
        return False
    shutil.rmtree(STORE / fid, ignore_errors=True)
    return True


def text(fid: str, limit: int = 16 * 1024 * 1024) -> str:
    """The whole text of a small file; ValueError beyond `limit` bytes."""
    p = path_of(fid)
    if p.stat().st_size > limit:
        raise ValueError(f"{p.name} is larger than {limit // (1024 * 1024)} MB: give it to a tool as a file instead")
    return p.read_text(encoding="utf-8", errors="replace")


def head(fid: str, nbytes: int = 4096) -> str:
    with path_of(fid).open("r", encoding="utf-8", errors="replace") as f:
        return f.read(nbytes)
