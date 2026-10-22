"""Plug-ins written in Python (ADR 0029): the same classes as the Earth Atlas's
src/plugins/plugin.js, over numpy, run in a server (an application's, or the
API's) rather than in a page. The classes belong to no application: the hosting
of plug-ins is science/hosting.py, and an application binds it to its own
plugins/ directory (applications/earth-atlas/py/host.py). A script exports

    def create(api):
        class Mine(api.Diagnostic): ...
        return Mine()

where `api` is this module. The classes are named by what the algorithm is, and
each writes run() from smaller hooks that a subclass overrides:
  Plugin            identity, parameters, needs and run(): the general case;
  Diagnostic        a field computed from fields, at every time of the input;
  ForecastModel     a state advanced in time from an analysis, verified against later analyses;
  DataAssimilation  a background and observations to an analysis and its increment;
  Trajectories      parcels carried through a flow.
Results are data: a field product (drawn, contoured, compared and read like any
product of the library), a sketch (polylines, markers, labels and moving points
in longitude and latitude), tables and notes. Times are datetimes in UTC; fields
are float32 arrays (nlat, nlon) with NaN where there is no value, rows north to south.
"""

from __future__ import annotations

import datetime as dt
import math

import numpy as np

from .fields import Field, bilinear, iso, UTC

EARTH_RADIUS_M = 6371008.8
DEG = math.pi / 180
HOUR = dt.timedelta(hours=1)

# ---- results ---------------------------------------------------------------------------


def field(product: Field, show: dict | None = None) -> dict:
    """A field product to show; `show` names what to select, by variable or vector name:
    {"surface": "zeta", "flow": "wind@500hPa", "lines": "z@500hPa"}."""
    return {"kind": "field", "product": product, "show": show or {}}


def sketch(title: str = "", polylines=None, markers=None, labels=None, moving=None) -> dict:
    """A sketch over the globe: polylines [{coords: [lon, lat, ...], colour, width, alpha, dashed}],
    markers [{lon, lat, colour, size}], labels [{lon, lat, text}], moving [{times: [ms], coords, colour, size}]."""
    return {"kind": "sketch", "title": title, "polylines": polylines or [], "markers": markers or [],
            "labels": labels or [], "moving": moving or []}


def table(title: str, columns: list, rows: list, note: str = "") -> dict:
    return {"kind": "table", "title": title, "columns": columns, "rows": rows, "note": note}


def note(text: str) -> dict:
    return {"kind": "note", "text": text}


# ---- products and numerics -------------------------------------------------------------------


def make_product(name: str, title: str, grid: dict, times: list, variables: list, *, model: str | None = None,
                 domain: str = "meteorology", vectors: list | None = None, source: dict | None = None,
                 verification: dict | None = None) -> Field:
    """A field product made in memory, with the interface of one read from the library.
    variables: [{name, long_name, standard_name, units, level?, difference?, slices: [array per time]}]."""
    slices = {v["name"]: [np.asarray(s, dtype=np.float32) for s in v["slices"]] for v in variables}
    described = [{k: v[k] for k in v if k != "slices"} for v in variables]
    meta = {
        "name": name, "computed": True, "title": title, "model": model or title, "domain": domain, "grid": grid,
        "times": [iso(t) if isinstance(t, dt.datetime) else t for t in times],
        "variables": described, "vectors": vectors or [],
        "source": {"provider": "computed in the server", "credit": title, **(source or {})},
    }
    if verification:
        meta["verification"] = verification
    return Field(meta, slices=slices)


def sample_on_grid(product: Field, name: str, index: int, grid: dict) -> np.ndarray:
    """A variable's values at one time, sampled bilinearly onto another grid."""
    lons = grid["lon0"] + np.arange(grid["nlon"]) * grid["dlon"]
    lats = grid["lat0"] - np.arange(grid["nlat"]) * grid["dlat"]
    LON, LAT = np.meshgrid(lons, lats)
    return bilinear(product.slice(name, index), product.grid, LON, LAT).astype(np.float32)


def rms_difference(a: np.ndarray, b: np.ndarray, grid: dict, lat_min: float = -90, lat_max: float = 90) -> float:
    """The root-mean-square difference of two fields, weighted by the cosine of latitude, between two latitudes."""
    lats = grid["lat0"] - np.arange(grid["nlat"]) * grid["dlat"]
    rows = (lats >= lat_min) & (lats <= lat_max)
    w = np.cos(lats[rows] * DEG)[:, None] * np.ones((1, grid["nlon"]))
    d = (np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64))[rows]
    ok = np.isfinite(d)
    return float(math.sqrt(np.sum(w[ok] * d[ok] ** 2) / np.sum(w[ok]))) if ok.any() else float("nan")


def lattice(n: int) -> list[dict]:
    """Points spread evenly over the sphere (a Fibonacci lattice): [{lon, lat}]."""
    golden = math.pi * (3 - math.sqrt(5))
    out = []
    for i in range(n):
        z = 1 - (2 * i + 1) / n
        lon = ((i * golden / DEG + 180) % 360) - 180
        out.append({"lon": lon, "lat": math.asin(z) / DEG})
    return out


def wrap_longitude(lon: float) -> float:
    return ((lon + 180) % 360 + 360) % 360 - 180


def parcel_step(flow, lon: float, lat: float, when: dt.datetime, seconds: float):
    """One fourth-order Runge–Kutta step of a parcel on the sphere through
    flow(lon, lat, when) -> (u, v) in m/s or None. Returns (lon, lat), or None."""
    def rate(lo, la, t):
        uv = flow(wrap_longitude(lo), max(-89.99, min(89.99, la)), t)
        if uv is None:
            return None
        return (uv[0] / (EARTH_RADIUS_M * math.cos(la * DEG)) / DEG, uv[1] / EARTH_RADIUS_M / DEG)
    half = dt.timedelta(seconds=seconds / 2)
    k1 = rate(lon, lat, when)
    if k1 is None:
        return None
    k2 = rate(lon + 0.5 * seconds * k1[0], lat + 0.5 * seconds * k1[1], when + half)
    if k2 is None:
        return None
    k3 = rate(lon + 0.5 * seconds * k2[0], lat + 0.5 * seconds * k2[1], when + half)
    if k3 is None:
        return None
    k4 = rate(lon + seconds * k3[0], lat + seconds * k3[1], when + 2 * half)
    if k4 is None:
        return None
    nlon = lon + seconds / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
    nlat = lat + seconds / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    if nlat > 90:
        nlat, nlon = 180 - nlat, nlon + 180
    elif nlat < -90:
        nlat, nlon = -180 - nlat, nlon + 180
    return wrap_longitude(nlon), nlat


# ---- the classes ----------------------------------------------------------------------------


class Plugin:
    """Identity, parameters, needs and run(): the general case."""
    id = "plugin"
    title = "A plug-in"
    description = ""
    domain = "meteorology"
    # {name: {label, default, min, max, step, unit}} or {name: {label, default, choices: [[value, label]]}}
    parameters: dict = {}

    def needs(self) -> list[dict]:
        """[{name, label, standard_name | standard_names: [u, v], level?, times?, distinct?, optional?}]; run() receives
        inputs[name] = Input(product, names, variables, meta, model, level) from the source the user chose.
        A need {name, label, kind: 'observations', variable, optional?} is served by the observation sets
        carrying the variable, and received as Observed(observations, name, model, variable); an optional
        need left without a source is received as None."""
        return []

    def run(self, inputs: dict, params: dict, context) -> list[dict]:
        return []

    def breathe(self, context, fraction: float, text: str = "") -> None:
        context.progress(fraction, text)

    @property
    def product_name(self) -> str:
        return f"plugin-{self.id}"


class Input:
    def __init__(self, product: Field, names: list[str], variables: list[dict], meta: dict, model: str, level: str | None):
        self.product, self.names, self.variables, self.meta, self.model, self.level = product, names, variables, meta, model, level


class Observed:
    """The source chosen for a need of kind 'observations' ({name, label, kind: 'observations',
    variable, optional?}): the observation set as the library holds it (schema observations/1:
    {fields, rows, variables, window, source}), its name, its title as `model`, and the variable
    asked for; as the page's host passes { observations, name, model, variable }."""

    def __init__(self, observations: dict, name: str, model: str, variable: str | None):
        self.observations, self.name, self.model, self.variable = observations, name, model, variable


class Context:
    """What a plug-in sees: the instant shown, the chosen place, the library's
    products (fields, and the index's entries, events and observations), progress.
    The host fills the readers in from the products on disk; the page's context
    (src/plugins/context.js) has the same members."""

    def __init__(self, time: dt.datetime, pinned=None, progress=None, fields=None, field=None,
                 listing=None, events=None, observations=None):
        self.time = time
        self.pinned = pinned
        self._progress = progress or (lambda fraction, text="": None)
        self._fields = fields or (lambda: [])
        self._field = field or (lambda name: None)
        self._list = listing or (lambda kind: [])
        self._events = events or (lambda name: None)
        self._observations = observations or (lambda name: None)

    def progress(self, fraction: float, text: str = "") -> None:
        self._progress(fraction, text)

    def fields(self) -> list[dict]:
        return self._fields()

    def field(self, name: str) -> Field | None:
        return self._field(name)

    def list(self, kind: str) -> list[dict]:
        """The library's index entries of a kind: 'fields', 'events', 'observations', 'tracks', 'regions', ..."""
        return self._list(kind)

    def events(self, name: str) -> dict | None:
        """An event set by name (earthquakes, mechanisms, disasters): {fields, events: [[row]], window, source}."""
        return self._events(name)

    def observations(self, name: str) -> dict | None:
        """An observation set by name (stations, buoys, sensors): {variables, fields, rows, window, source},
        read again each time, since observations change while a server runs."""
        return self._observations(name)


class Diagnostic(Plugin):
    """A field computed from fields, at every time of the first need, on its grid.
    Override outputs(), compute(), and vectors() and show() if wanted."""

    def outputs(self) -> list[dict]:
        return []

    def vectors(self) -> list[dict]:
        return []

    def compute(self, inputs: dict, params: dict, index: int, context) -> dict:
        raise NotImplementedError(f"{self.id}: compute() is not written")

    def show(self) -> dict:
        outputs = self.outputs()
        return {"surface": outputs[0]["name"]} if outputs else {}

    def run(self, inputs, params, context):
        first = inputs[self.needs()[0]["name"]]
        product = first.product
        outputs = self.outputs()
        lists = {o["name"]: [] for o in outputs}
        n = max(1, len(product.times))
        for k in range(n):
            values = self.compute(inputs, params, k, context)
            for o in outputs:
                lists[o["name"]].append(np.asarray(values[o["name"]], dtype=np.float32))
            self.breathe(context, (k + 1) / n, f"time {k + 1} of {n}")
        made = make_product(
            self.product_name, self.title, product.grid, product.times or [context.time],
            [{"level": first.level, **o, "slices": lists[o["name"]]} for o in outputs],
            model=f"{self.title} (plug-in)", domain=self.domain, vectors=self.vectors(),
            source={"dataset": first.meta["name"], "provider": f"computed in the server from {first.model}",
                    "credit": f"{self.title} (plug-in), from {first.model}"})
        return [field(made, self.show())]


class ForecastModel(Plugin):
    """A state advanced in time from an analysis: the first need gives the start,
    at the input time at or before the instant shown. Override grid(), outputs(),
    initial(), advance() and record(); verify() names what to compare with the input."""
    dt = 1800          # the step, in seconds
    every = 6          # hours between recorded times
    parameters = {"hours": {"label": "Forecast length", "default": 72, "min": 6, "max": 240, "step": 6, "unit": "h"}}

    def grid(self) -> dict:
        raise NotImplementedError(f"{self.id}: grid() is not written")

    def outputs(self) -> list[dict]:
        return []

    def vectors(self) -> list[dict]:
        return []

    def verify(self) -> dict | None:
        """{output, need, label, unit}: the output variable and the need holding the truth, or None."""
        return None

    def initial(self, inputs, index, params, context):
        raise NotImplementedError(f"{self.id}: initial() is not written")

    def advance(self, state, seconds, params) -> None:
        raise NotImplementedError(f"{self.id}: advance() is not written")

    def record(self, state, params) -> dict:
        raise NotImplementedError(f"{self.id}: record() is not written")

    def show(self) -> dict:
        outputs = self.outputs()
        return {"surface": outputs[0]["name"]} if outputs else {}

    def run(self, inputs, params, context):
        first = inputs[self.needs()[0]["name"]]
        product = first.product
        i0 = product.index_at_or_before(context.time)
        t0 = product.times[i0]
        hours = float(params.get("hours", 72))
        steps = round(hours * 3600 / self.dt)
        record_every = max(1, round(self.every * 3600 / self.dt))
        outputs = self.outputs()
        lists = {o["name"]: [] for o in outputs}
        times = [t0]
        state = self.initial(inputs, i0, params, context)

        def take():
            r = self.record(state, params)
            for o in outputs:
                lists[o["name"]].append(np.asarray(r[o["name"]], dtype=np.float32))
        take()
        for s in range(1, steps + 1):
            self.advance(state, self.dt, params)
            if s % record_every == 0:
                take()
                times.append(t0 + dt.timedelta(seconds=s * self.dt))
            if s % 4 == 0:
                self.breathe(context, s / steps, f"+{round(s * self.dt / 3600)} h")
        grid = self.grid()
        v = self.verify()
        verification, rows, truth = None, [], None
        if v and v["need"] in inputs and t0 in inputs[v["need"]].product.times:
            truth = inputs[v["need"]]
            name = truth.names[0]
            start = sample_on_grid(truth.product, name, truth.product.times.index(t0), grid)
            now = dt.datetime.now(UTC)
            for k, t in enumerate(times):
                if k == 0 or t not in truth.product.times:
                    continue
                analysed = sample_on_grid(truth.product, name, truth.product.times.index(t), grid)
                rows.append({"hours": (t - t0) / HOUR, "model": rms_difference(lists[v["output"]][k], analysed, grid, 20, 80),
                             "persistence": rms_difference(start, analysed, grid, 20, 80),
                             "against": "analysis" if t <= now else "forecast"})
            verification = {"quantity": v["label"], "against": truth.model, "start": iso(t0), "leads": rows}
        made = make_product(
            self.product_name, self.title, grid, times, [{**o, "slices": lists[o["name"]]} for o in outputs],
            model=f"{self.title} (plug-in)", domain=self.domain, vectors=self.vectors(), verification=verification,
            source={"dataset": first.meta["name"], "provider": f"computed in the server from {first.model} at {iso(t0)}",
                    "credit": f"{self.title} (plug-in), from {first.model}"})
        results = [field(made, self.show())]
        if rows:
            unit = f" {v['unit']}" if v.get("unit") else ""
            shown = [r for r in rows if r["hours"] % 24 == 0] if len(rows) > 8 else rows
            results.append(table(
                f"Root-mean-square error of {v['label']}, 20–80° N, against {truth.model}", ["Lead", self.title, "Persistence"],
                [[f"+{int(r['hours'])} h", f"{round(r['model'])}{unit}", f"{round(r['persistence'])}{unit}"] for r in shown],
                "Beyond the present, against the forecast rather than analyses." if any(r["against"] == "forecast" for r in rows) else "Against analyses."))
        return results


class DataAssimilation(Plugin):
    """A background field and observations to an analysis, at the instant shown
    (clamped to the background's span), the background interpolated linearly in
    time to it, so that hourly reports and a sensor's latest reading meet a
    background of six-hourly times. Override observations() and analyse(); the
    result is the analysis and its increment at that one time, with the observations sketched."""

    def observations(self, inputs, when: dt.datetime, params, context) -> list[dict]:
        """[{lon, lat, value, error?}] observed at `when`, the analysis's time (UTC)."""
        return []

    def analyse(self, background: np.ndarray, grid: dict, observations: list[dict], params: dict, context) -> np.ndarray:
        return background

    def markers(self, observations: list[dict], innovations) -> list[dict]:
        return [{"lon": o["lon"], "lat": o["lat"], "size": 2.6, "colour": "#f6a07e" if d > 0 else "#8fb8e8"}
                for o, d in zip(observations, innovations)]

    def run(self, inputs, params, context):
        first = inputs[self.needs()[0]["name"]]
        product, name, v = first.product, first.names[0], first.variables[0]
        times = product.times
        t = max(times[0], min(times[-1], context.time)) if times else context.time
        i0, i1, w = product.bracket(t, name)
        a0 = product.slice(name, max(0, i0))
        if w > 0:      # in double precision, rounded once to float32, as the page interpolates
            a1 = product.slice(name, max(0, i1))
            background = (a0.astype(np.float64) + (a1.astype(np.float64) - a0) * w).astype(np.float32)
        else:
            background = a0
        observations = self.observations(inputs, t, params, context)
        lons = np.array([o["lon"] for o in observations]) if observations else np.zeros(0)
        lats = np.array([o["lat"] for o in observations]) if observations else np.zeros(0)
        innovations = np.array([o["value"] for o in observations]) - product.sample(name, t, lons, lats) if observations else np.zeros(0)
        analysis = np.asarray(self.analyse(background, product.grid, observations, params, context), dtype=np.float32)
        increment = analysis - background
        made = make_product(
            self.product_name, self.title, product.grid, [t],
            [{"name": v["name"], "long_name": f"{v.get('long_name') or v['name']}, analysed", "standard_name": v.get("standard_name"),
              "units": v.get("units"), "level": v.get("level"), "slices": [analysis]},
             {"name": f"{v['name']}_increment", "long_name": f"Analysis increment of {(v.get('long_name') or v['name']).lower()}",
              "units": v.get("units"), "level": v.get("level"), "difference": True, "slices": [increment]}],
            model=f"{self.title} (plug-in)", domain=self.domain,
            source={"dataset": first.meta["name"], "credit": f"{self.title} (plug-in), from {first.model}",
                    "provider": f"analysed in the server from {first.model} at {iso(t)} and {len(observations)} observations"})
        return [field(made, {"surface": f"{v['name']}_increment"}),
                sketch(f"{len(observations)} observations (warm where above the background)", markers=self.markers(observations, innovations))]


class Trajectories(Plugin):
    """Parcels carried through a flow: the first need is a vector field. Override
    seeds(); velocity() samples the need by default. The result is a sketch of the paths,
    with the parcels moving along them as the time shown changes."""
    parameters = {
        "hours": {"label": "Hours ahead", "default": 48, "min": 6, "max": 240, "step": 6, "unit": "h"},
        "step": {"label": "Time step", "default": 1, "min": 0.25, "max": 6, "step": 0.25, "unit": "h"},
    }
    colour = "#e8c3dc"

    def seeds(self, inputs, params, context) -> list[dict]:
        return []

    def velocity(self, inputs, params):
        w = inputs[self.needs()[0]["name"]]
        u, v = w.names

        def flow(lon, lat, when):
            a = float(w.product.sample(u, when, np.array([lon]), np.array([lat]))[0])
            b = float(w.product.sample(v, when, np.array([lon]), np.array([lat]))[0])
            return (a, b) if math.isfinite(a) and math.isfinite(b) else None
        return flow

    def run(self, inputs, params, context):
        w = inputs[self.needs()[0]["name"]]
        times = w.product.times
        flow = self.velocity(inputs, params)
        t_end = times[-1]
        t0 = max(times[0], min(context.time, t_end))
        seconds = float(params.get("step", 1)) * 3600
        steps = max(0, min(round(float(params.get("hours", 48)) * 3600 / seconds), int((t_end - t0).total_seconds() // seconds)))
        seeds = self.seeds(inputs, params, context)
        paths = []
        for s, seed in enumerate(seeds):
            coords, when = [seed["lon"], seed["lat"]], [t0]
            lon, lat = seed["lon"], seed["lat"]
            for k in range(steps):
                nxt = parcel_step(flow, lon, lat, t0 + dt.timedelta(seconds=k * seconds), seconds)
                if nxt is None:
                    break
                lon, lat = nxt
                coords += [lon, lat]
                when.append(t0 + dt.timedelta(seconds=(k + 1) * seconds))
            paths.append({"coords": coords, "times": [int(t.timestamp() * 1000) for t in when], "colour": seed.get("colour", self.colour)})
            if s % 8 == 7:
                self.breathe(context, (s + 1) / len(seeds), f"parcel {s + 1} of {len(seeds)}")
        return [sketch(
            f"{len(paths)} parcels through the wind of {w.model}, from {iso(t0)} for {round(steps * seconds / 3600)} h",
            polylines=[{"coords": p["coords"], "colour": p["colour"], "width": 1.3, "alpha": 0.85} for p in paths],
            moving=[{"times": p["times"], "coords": p["coords"], "colour": "#f6c47e", "size": 3} for p in paths])]
