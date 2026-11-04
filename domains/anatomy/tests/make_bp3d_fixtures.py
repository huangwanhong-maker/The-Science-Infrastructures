"""Write the anatomy domain's fixture of BodyParts3D's published volumes, from the database's own files.

Each OBJ file of BodyParts3D 4.0 carries a header with its element's volume in cm³ as the database computed
it ("# Volume(cm3): ..."). This script reads those headers from the zip that the data library fetched
(python data/acquisition/fetch.py bodyparts3d) and writes them, with the archive's checksum and the date, to
tests/fixtures/bodyparts3d-volumes.tsv, against which check_measures.cpp holds the volumes of the solids
the meshes bound. BodyParts3D is © The Database Center for Life Science, CC BY 4.0.

Usage::

    python domains/anatomy/tests/make_bp3d_fixtures.py
"""

from __future__ import annotations

import datetime as dt
import hashlib
import re
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ZIP = ROOT / "data" / "cache" / "bodyparts3d" / "isa_BP3D_4.0_obj_99.zip"
OUT = HERE / "fixtures" / "bodyparts3d-volumes.tsv"
FIELD = re.compile(r"^# ([^:\n]+?) ?: (.*)$", re.M)


def main() -> None:
    if not ZIP.exists():
        raise SystemExit(f"{ZIP} is missing; run: python data/acquisition/fetch.py bodyparts3d")
    digest = hashlib.sha256(ZIP.read_bytes()).hexdigest()
    rows = []
    with zipfile.ZipFile(ZIP) as z:
        for name in sorted(n for n in z.namelist() if n.endswith(".obj")):
            head = z.read(name)[:2048].decode("utf-8", "replace")
            fields = dict(FIELD.findall(head))
            rows.append((Path(name).stem, fields.get("Volume(cm3)", "").strip(), fields.get("English name", "").strip()))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# BodyParts3D 4.0: each element's volume in cm³ as published in its OBJ file's header.",
        f"# Source: {ZIP.name} (sha256 {digest}), from https://dbarchive.biosciencedbc.jp/data/bodyparts3d/20130619/;",
        f"# read {dt.date.today().isoformat()} by domains/anatomy/tests/make_bp3d_fixtures.py.",
        "# BodyParts3D, © The Database Center for Life Science, licensed under CC Attribution 4.0 International.",
        "id\tvolume_cm3\tname",
    ] + [f"{i}\t{v}\t{n}" for i, v, n in rows]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(rows)} volumes to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
