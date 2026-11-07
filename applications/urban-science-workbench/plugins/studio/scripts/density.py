"""How thick events lie, through time: a Gaussian kernel density of the events shown, on a grid over them, at equal steps of time, each step counting a window of events
behind it. A surface one can watch move. An example of a script in Python: it is given the records shown and answers a record of its own, here a grid, which the workbench
draws on the ground through time.

The kernel is an isotropic Gaussian on the sphere, near enough for bandwidths far below a radian: a cell's distance from an event is measured on the chord of its
latitude. Its integral over the grid is the count of events in the window, up to what falls outside the grid.

Limits: a window of events is a boxcar in time, so that the surface jumps as a large event enters; the edges of the grid lose what their kernels spread beyond them.
"""

import numpy as np


def create(api):
    class Density(api.SpaceTimeAnalysis):
        id = "density"
        title = "How thick events lie, through time"
        description = "A kernel density of the events shown, on a grid over them, at equal steps of time."
        params = [
            api.param("bandwidth_km", "Bandwidth (km)", 250.0, minimum=1, maximum=2000),
            api.param("steps", "Steps of time", 40, kind="integer", minimum=2, maximum=200),
            api.param("window", "Window (share of the span)", 0.1, minimum=0.01, maximum=1, step=0.01),
            api.param("width", "Grid across", 240, kind="integer", minimum=20, maximum=720),
        ]

        def run(self, records, params, context):
            found = [r for r in records if r.kind == "events"]
            if not found:
                raise ValueError("This needs events: show a record of events first.")
            r = found[0]
            ok = np.isfinite(r.lat) & np.isfinite(r.lon) & np.isfinite(r.time)
            lat, lon, time = r.lat[ok], r.lon[ok], r.time[ok]
            pad = params["bandwidth_km"] * 2 / 111.0
            south, north = max(-90.0, lat.min() - pad), min(90.0, lat.max() + pad)
            west, east = max(-180.0, lon.min() - pad), min(180.0, lon.max() + pad)
            width = int(params["width"])
            height = max(2, int(round(width * (north - south) / max(1e-6, east - west))))
            ys = south + (np.arange(height) + 0.5) * (north - south) / height
            xs = west + (np.arange(width) + 0.5) * (east - west) / width
            h2 = (params["bandwidth_km"] / 111.0) ** 2
            steps = int(params["steps"])
            t0, t1 = time.min(), time.max()
            at = np.linspace(t0, t1, steps)
            window = params["window"] * (t1 - t0)
            frames = np.zeros((steps, height, width), dtype=float)
            cos_y = np.cos(np.radians(ys))[:, None]
            reach_y = int(np.ceil(4 * np.sqrt(h2) / ((north - south) / height)))
            reach_x = int(np.ceil(4 * np.sqrt(h2) / ((east - west) / width)))
            for k, t in enumerate(at):
                keep = (time <= t) & (time > t - window)
                context.progress(k / steps, f"step {k + 1} of {steps}: {int(keep.sum())} events")
                if not keep.any():
                    continue
                # Binned onto the grid first, so that the cost grows with the grid and not with the catalogue; each bin reaches four bandwidths, beyond which
                # its kernel is below a part in three thousand.
                counts, _, _ = np.histogram2d(lat[keep], lon[keep], bins=[height, width], range=[[south, north], [west, east]])
                row, col = np.nonzero(counts)
                for i, j, c in zip(row, col, counts[row, col]):
                    i0, i1 = max(0, i - reach_y), min(height, i + reach_y + 1)
                    x_reach = int(np.ceil(reach_x / max(0.05, cos_y[i, 0])))
                    j0, j1 = max(0, j - x_reach), min(width, j + x_reach + 1)
                    d2 = (ys[i0:i1, None] - ys[i]) ** 2 + ((xs[None, j0:j1] - xs[j]) * cos_y[i0:i1]) ** 2
                    frames[k, i0:i1, j0:j1] += c * np.exp(-0.5 * d2 / h2)
            frames /= 2 * np.pi * (params["bandwidth_km"] ** 2)
            frames[frames < frames.max() * 1e-3] = np.nan
            return [
                api.grid(f"How thick: {r.title}", at, [south, west, north, east], {"events per km²": frames}, units={"events per km²": ""},
                         show={"map": "plum", "window_ms": float((t1 - t0) / steps)}, note=f"A Gaussian kernel of {params['bandwidth_km']:g} km over a window of {window / 86400000:.0f} days."),
                api.note(f"{int(ok.sum()):,} events of {r.title}, on a grid of {width} × {height}, at {steps} steps."),
            ]

    return Density()
