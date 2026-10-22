"""The Science Infrastructures from Python: the API's client, and the library
read on this computer without a server.

    import science
    api = science.Client("http://127.0.0.1:8800")          # a service (api/service/server.py), or an atlas's server
    gfs = api.field("gfs-forecast")
    gfs.sample("tmp2m", lon=139.69, lat=35.69, time="now+1d")  # kelvin, at Tokyo, tomorrow

    here = science.Local()                                  # the same reads, from data/cache, with no server

`Client` needs only the standard library (numpy, when present, makes grids
arrays); `Local` needs numpy for fields. The modules `fields` (field products in
the page's format), `plugins` (the classes of plug-ins run in a server) and
`hosting` (their host) serve the servers and the plug-ins, and import numpy.
See api/README.md.
"""

from __future__ import annotations

from .client import ApiError, Client
from .local import Local

__version__ = "0.1.0"
__all__ = ["ApiError", "Client", "Local", "__version__"]
