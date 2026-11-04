"""The Medical Segmentation Decathlon's hepatic vessels (Task 08) unpacked for the C++ readers (ADR 0070).

The registry's dataset `msd-hepatic-vessels` holds four training cases as gzipped NIfTI-1, read from the archive by
range requests (data/acquisition/fetch.py). The C++ reads uncompressed NIfTI-1 (eikon's grid/nifti.hpp), so this
product writes each case's image and label, unpacked, to data/cache/products/msd-hepatic-vessels/, with an index:

    <case>-image.nii   CT in Hounsfield units (as stored, with its scaling)
    <case>-label.nii   0 background, 1 hepatic vessel, 2 tumour
    index.json         the cases, their files and the dataset's credit

Usage: python data/products/msd_vessels.py   (after python data/acquisition/fetch.py msd-hepatic-vessels)
"""
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache" / "msd-hepatic-vessels"
OUT = ROOT / "data" / "cache" / "products" / "msd-hepatic-vessels"
REGISTRY = ROOT / "data" / "registry" / "msd-hepatic-vessels.json"


def main() -> int:
    manifest = json.loads(REGISTRY.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    cases = {}
    for entry in manifest["files"]:
        src = CACHE / entry["file"]
        if not src.exists():
            print(f"missing {src}: run python data/acquisition/fetch.py msd-hepatic-vessels", file=sys.stderr)
            return 1
        case, kind = entry["file"].removesuffix(".nii.gz").rsplit("-", 1)
        dest = OUT / f"{case}-{kind}.nii"
        if not dest.exists() or dest.stat().st_mtime < src.stat().st_mtime:
            dest.write_bytes(gzip.decompress(src.read_bytes()))
            print(f"  {dest.name}: {dest.stat().st_size / 1e6:.1f} MB")
        cases.setdefault(case, {})[kind] = dest.name
    index = {"dataset": manifest["id"], "credit": manifest["credit"], "licence": manifest["licence"], "cases": cases}
    (OUT / "index.json").write_text(json.dumps(index, indent=1), encoding="utf-8")
    print(f"{len(cases)} cases in {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
