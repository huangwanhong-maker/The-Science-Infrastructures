"""Relative vorticity and divergence of the wind, in Python: the algorithm of
vorticity.js, run by the server on numpy (applications/earth-atlas/py). Centred
differences on the sphere,

    ζ = (1 / (a cos φ)) [ ∂v/∂λ − ∂(u cos φ)/∂φ ],    δ = (1 / (a cos φ)) [ ∂u/∂λ + ∂(v cos φ)/∂φ ],

rows running north to south; the poles and, on a grid that does not wrap, the
edge columns are left without a value.
"""

import numpy as np

RADIUS = 6371008.8
DEG = np.pi / 180


def create(api):
    class Vorticity(api.Diagnostic):
        id = "vorticity-py"
        title = "Vorticity and divergence (Python)"
        description = "Relative vorticity and divergence of the wind, by centred differences on the sphere, computed in Python by the server."

        def needs(self):
            return [{"name": "wind", "label": "Wind", "standard_names": ["eastward_wind", "northward_wind"]}]

        def outputs(self):
            return [
                {"name": "zeta", "long_name": "Relative vorticity", "standard_name": "atmosphere_relative_vorticity", "units": "s-1"},
                {"name": "delta", "long_name": "Divergence of the wind", "standard_name": "divergence_of_wind", "units": "s-1"},
            ]

        def compute(self, inputs, params, index, context):
            wind = inputs["wind"]
            u = wind.product.slice(wind.names[0], index).astype(np.float64)
            v = wind.product.slice(wind.names[1], index).astype(np.float64)
            g = wind.product.grid
            lats = g["lat0"] - np.arange(g["nlat"]) * g["dlat"]
            cos = np.cos(lats * DEG)[:, None]
            dx, dy = 2 * g["dlon"] * DEG * RADIUS, 2 * g["dlat"] * DEG * RADIUS
            # Neighbours east and west, wrapping round the Earth on a periodic grid.
            if g.get("periodic"):
                ue, uw, ve, vw = np.roll(u, -1, axis=1), np.roll(u, 1, axis=1), np.roll(v, -1, axis=1), np.roll(v, 1, axis=1)
            else:
                ue, uw, ve, vw = (np.full_like(a, np.nan) for a in (u, u, v, v))
                ue[:, :-1], uw[:, 1:], ve[:, :-1], vw[:, 1:] = u[:, 1:], u[:, :-1], v[:, 1:], v[:, :-1]
            ucos, vcos = u * cos, v * cos
            zeta, delta = np.full_like(u, np.nan), np.full_like(u, np.nan)
            # Interior rows: the row above is north (φ + Δφ), the row below south.
            inner = slice(1, -1)
            with np.errstate(invalid="ignore", divide="ignore"):
                zeta[inner] = ((ve - vw)[inner] / dx - (ucos[:-2] - ucos[2:]) / dy) / cos[inner]
                delta[inner] = ((ue - uw)[inner] / dx + (vcos[:-2] - vcos[2:]) / dy) / cos[inner]
            polar = np.abs(cos[:, 0]) < 1e-6
            zeta[polar], delta[polar] = np.nan, np.nan
            return {"zeta": zeta.astype(np.float32), "delta": delta.astype(np.float32)}

    return Vorticity()
