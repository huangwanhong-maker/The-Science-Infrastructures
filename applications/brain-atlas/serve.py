"""Serve the Brain Atlas on this computer.

The page is served from this directory; the code it shares, from the
repository's ``domains/``, ``foundations/`` and ``visualization/`` under the
same paths; and the atlas data from the repository's cache (``data/cache/``)
under ``/data/``. Only 127.0.0.1 is bound unless ``--host`` says otherwise.

Usage::

    python applications/brain-atlas/serve.py [--host 127.0.0.1] [--port 8810] [--open]

then open http://127.0.0.1:8810/ in a browser.
"""

from __future__ import annotations

import argparse
import http.server
import sys
import urllib.parse
import webbrowser
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
DATA = ROOT / "data" / "cache"
MOUNTS = {
    "/data/": DATA,
    "/domains/": ROOT / "domains",
    "/foundations/": ROOT / "foundations",
    "/visualization/": ROOT / "visualization",
}

# Explicit types: on Windows the registry may otherwise map .js to text/plain,
# which browsers refuse for ES modules.
TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json",
    ".tsv": "text/tab-separated-values; charset=utf-8",
    ".gz": "application/gzip",
    ".nrrd": "application/octet-stream",
    ".bin": "application/octet-stream",
    ".svg": "image/svg+xml",
    ".png": "image/png",
}


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **TYPES}

    def translate_path(self, path: str) -> str:
        path = urllib.parse.unquote(path.split("?", 1)[0].split("#", 1)[0])
        prefix = next((p for p in MOUNTS if path.startswith(p)), None)
        base, rel = (MOUNTS[prefix], path[len(prefix):]) if prefix else (APP, path.lstrip("/"))
        target = (base / rel).resolve()
        if target != base and base not in target.parents:
            return str(base / "__outside__")
        return str(target)

    def guess_type(self, path: str) -> str:
        # A .nii.gz is served as gzip bytes, which the page inflates itself.
        if str(path).endswith(".nii.gz"):
            return "application/gzip"
        return super().guess_type(path)

    def end_headers(self) -> None:
        # The cache is refreshed in place when data are fetched again, so it is
        # revalidated rather than cached for long.
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, format: str, *args) -> None:
        if args and str(args[1]).startswith(("4", "5")):
            super().log_message(format, *args)


class Server(http.server.ThreadingHTTPServer):
    request_queue_size = 128
    daemon_threads = True


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Serve the Brain Atlas locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8810)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    args = parser.parse_args(argv)
    if not (DATA / "templateflow-mni152" / "tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz").exists():
        print("The MNI template is missing. From the repository root, run:\n"
              "  python data/acquisition/fetch.py templateflow-mni152 schaefer-2018 ebrains-julich-brain allen-mouse-ccf\n"
              "and see applications/brain-atlas/README.md.", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if args.host == "127.0.0.1" else f"{url} and on {args.host}"
    print(f"Brain Atlas at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
