"""Builds the self-contained replay page from viewer.html and the record's results.

    python build_viewer.py [output.html]

The page embeds results/trajectories.json and results/landscapes.json, so it opens from a file or as a
published page with no server. Default output: viewer/searches-on-a-landscape.html (not committed; built).
"""

import json
import pathlib
import sys

here = pathlib.Path(__file__).resolve().parent
results = here.parent / "results"
data = {
    "trajectories": json.loads((results / "trajectories.json").read_text(encoding="utf-8")),
    "landscapes": json.loads((results / "landscapes.json").read_text(encoding="utf-8")),
}
text = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
page = (here / "viewer.html").read_text(encoding="utf-8").replace("/*DATA*/", text)
out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else here / "searches-on-a-landscape.html"
out.write_text(page, encoding="utf-8")
print(f"{out} ({out.stat().st_size / 1e6:.2f} MB)")
