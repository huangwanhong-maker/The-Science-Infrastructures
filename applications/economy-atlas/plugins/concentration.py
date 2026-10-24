"""The concentration of each economy's exports over its partners: the Herfindahl–Hirschman
index (Hirschman 1945; Herfindahl 1950),

    H_i = Σ_j (F_ij / X_i)²,

1 when everything goes to one partner, 1/k when it is spread evenly over k. An economy whose
exports depend on few partners is exposed to them; a high H with a large partner share names
the dependence.
"""

import numpy as np


def create(api):
    class Concentration(api.NetworkAnalysis):
        id = "concentration"
        title = "Concentration of exports"
        description = ("The Herfindahl–Hirschman index of each economy's exports over its partners: 1 when one partner "
                       "takes everything, 1/k when k partners share alike.")
        params = [
            api.param("threshold", "Least exports listed", 10000.0, minimum=0,
                      help="in the network's unit; the table lists only economies exporting more"),
        ]

        def run(self, network, params, context):
            F, X = network.matrix(), network.exports()
            exporting = X > 0
            shares = np.divide(F, X[:, None], out=np.zeros_like(F), where=exporting[:, None])
            H = (shares ** 2).sum(axis=1)
            main = shares.argmax(axis=1)
            index = {network.ids[k]: float(H[k]) for k in range(network.n) if exporting[k]}
            listed = [k for k in np.argsort(-H) if exporting[k] and X[k] > params["threshold"]][:15]
            rows = [[network.places[k]["name"], float(H[k]), network.places[main[k]]["name"], float(shares[k, main[k]])]
                    for k in listed]
            return [
                api.values("Concentration of exports (Herfindahl–Hirschman)", index, scale="sequential",
                           note="1: one partner takes everything; 1/k: k partners share alike."),
                api.table("Most concentrated", ["economy", "H", "main partner", "its share"], rows,
                          note=f"Among the {int((X > params['threshold']).sum())} economies exporting more than "
                               f"{params['threshold']:g} {network.unit}.".rstrip()),
            ]

    return Concentration()
