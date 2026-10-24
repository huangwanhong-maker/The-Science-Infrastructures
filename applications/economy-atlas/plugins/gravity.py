"""The gravity equation of trade (Tinbergen 1962; Anderson and van Wincoop 2003 for its theory):
the flow from i to j grows with the exporter's and the importer's size and falls with the
distance between them,

    ln F_ij = c + a ln X_i + b ln M_j − d ln D_ij + ε_ij,

X_i the exports of i, M_j the imports of j, D_ij the great-circle distance (floored at 100 km,
so that neighbours sharing a border are not infinitely near). Fitted by ordinary least squares
over the positive flows of the network. Its residuals are the trade that size and distance do
not explain: an economy whose exports are mostly above the line trades with its partners more
than the equation predicts (a common language, a trade agreement, a supply chain). The arcs are
the flows whose value most exceeds the prediction, actual − exp(fitted), so that the large flows
above the line show, not the smallest, whose ratios are the largest.

The run is one least-squares solve on the flows: a few hundred milliseconds for the world's.
"""

import numpy as np


def create(api):
    class Gravity(api.NetworkAnalysis):
        id = "gravity"
        title = "Gravity equation"
        description = ("Tinbergen's gravity equation fitted to the flows: ln F = c + a ln X + b ln M − d ln D. "
                       "Colours each economy by how far its exports lie above or below the line; arcs are the flows "
                       "furthest above it.")
        params = [
            api.param("minimum", "Least flow fitted", 1.0, minimum=0, help="in the network's unit; smaller flows are left out"),
            api.param("arcs", "Arcs drawn", 40, kind="integer", minimum=0, maximum=500),
        ]

        def run(self, network, params, context):
            n = network.n
            F, X, M = network.matrix(), network.exports(), network.imports()
            D = np.maximum(network.distances_km(), 100.0)
            context.progress(0.2, "choosing the flows")
            mask = (F >= params["minimum"]) & (F > 0) & ~np.eye(n, dtype=bool)
            mask &= (X[:, None] > 0) & (M[None, :] > 0)
            i, j = np.nonzero(mask)
            if len(i) < 5:
                raise ValueError(f"{len(i)} flows above the least, too few to fit four coefficients")
            y = np.log(F[i, j])
            A = np.column_stack([np.ones(len(i)), np.log(X[i]), np.log(M[j]), -np.log(D[i, j])])
            context.progress(0.5, f"fitting {len(i)} flows")
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            residual = y - A @ coef
            dof = max(len(i) - A.shape[1], 1)
            sigma2 = float(residual @ residual) / dof
            se = np.sqrt(np.clip(np.diag(sigma2 * np.linalg.pinv(A.T @ A)), 0, None))
            total = float(((y - y.mean()) ** 2).sum())
            r2 = 1 - float(residual @ residual) / total if total > 0 else 1.0
            # Each exporter's mean residual over its flows.
            count = np.bincount(i, minlength=n)
            mean = np.bincount(i, weights=residual, minlength=n) / np.maximum(count, 1)
            by_exporter = {network.ids[k]: float(mean[k]) for k in range(n) if count[k] > 0}
            excess = F[i, j] - np.exp(A @ coef)
            order = [k for k in np.argsort(-excess) if excess[k] > 0][: params["arcs"]]
            above = [{"from": network.ids[i[k]], "to": network.ids[j[k]], "value": float(F[i[k], j[k]])} for k in order]
            unit = network.unit
            rows = [[name, float(c), float(s)] for name, c, s in
                    zip(["constant c", "exporter's size a", "importer's size b", "distance d"], coef, se)]
            top = [[f"{network.place(network.ids[i[k]])['name']} → {network.place(network.ids[j[k]])['name']}",
                    float(F[i[k], j[k]]), float(excess[k]), float(np.exp(residual[k]))] for k in order[:12]]
            return [
                api.values("Exports above the gravity line (mean log residual)", by_exporter, unit="ln(actual / predicted)",
                           scale="diverging",
                           note="Above 0: an economy's exports exceed what its size, its partners' and the distances predict."),
                api.arcs("Flows furthest above the line", above, unit=unit,
                         note="The flows whose value most exceeds the gravity equation's prediction, actual − predicted."),
                api.table("The gravity equation", ["coefficient", "estimate", "standard error"], rows,
                          note=f"Least squares over {len(i)} flows; R² = {r2:.3f}. A distance elasticity d near 1 is the "
                               f"literature's usual finding (Disdier and Head 2008)."),
                api.table("Most above the line", ["flow", "value", "above the prediction", "actual / predicted"], top,
                          note=f"Values in {unit}." if unit else ""),
            ]

    return Gravity()
