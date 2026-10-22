"""A server plug-in run through the service: the Earth Atlas's vorticity and
divergence in Python (applications/earth-atlas/plugins/vorticity.py), on the
source it would choose by default, or on another; then its result read like any
field. One run is made at a time; the call waits until it is done.

    python api/examples/python/run_plugin.py [--level "gfs-levels-forecast|500 hPa"] [--url ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import science  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--level", help="the wind to use, as product|level (see the description's candidates)")
parser.add_argument("--url")
args = parser.parse_args()

api = science.Client(args.url)
described = api.describe_plugin("earth-atlas", "plugins/vorticity.py")
print(f"{described['title']}: {described['description']}")
need = described["needs"][0]
print(f"  sources of '{need['label']}': {', '.join(c['key'] for c in need['candidates'])}")
chosen = {need["name"]: args.level} if args.level else None
job = api.run_plugin("earth-atlas", "plugins/vorticity.py", chosen=chosen)
print(f"  done in {job['finished'] - job['started']:.1f} s")
for result in job["results"]:
    if result["kind"] == "field":
        made = api.field(result["field"])
        zeta = made.sample("zeta", 139.69, 35.69, time="now")
        print(f"  {made.name}: relative vorticity over Tokyo now {zeta:.2e} s-1")
