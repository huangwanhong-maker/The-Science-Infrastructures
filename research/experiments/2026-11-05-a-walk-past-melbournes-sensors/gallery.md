### A Walk past Melbourne's Sensors

The shortest closed walk through the City of Melbourne's 134 pedestrian counters, on geodesics of the WGS 84 ellipsoid, searched by every permutation search of heuresis at a budget of a million evaluations ([the record](../research/experiments/2026-11-05-a-walk-past-melbournes-sensors/README.md)).

<table>
<tr>
<td width="50%" valign="top"><a href="../assets/gallery/search-30-melbourne-best-tour.jpg"><img src="../assets/gallery/search-30-melbourne-best-tour.jpg" alt="The best known closed tour through 134 pedestrian sensors drawn over Melbourne's streets, with the results of 32 restarts"></a><br><b>The best known walk.</b> 26.31 km through all 134 counters of the City of Melbourne's Pedestrian Counting System (positions CC BY 4.0), by iterated local search with white-box 2-opt and double-bridge kicks; 27 of 32 restarts end on it. Distances are Karney's geodesics from the geodesy domain; streets © OpenStreetMap contributors (ODbL), from the data library's city model.</td>
<td width="50%" valign="top"><a href="../assets/gallery/search-31-melbourne-tours-by-algorithm.jpg"><img src="../assets/gallery/search-31-melbourne-tours-by-algorithm.jpg" alt="Eight small maps, each the best tour of one search algorithm after a million evaluations"></a><br><b>Each search's best walk.</b> The shortest tour among each algorithm's runs at 10⁶ evaluations: random tours, the Plackett–Luce EDA, the GA with cycle, partially mapped and order crossover, the MAX–MIN ant system, simulated annealing and iterated local search, each a black box on heuresis's ask–tell protocol.</td>
</tr>
<tr>
<td width="50%" valign="top"><a href="../assets/gallery/search-32-melbourne-convergence.jpg"><img src="../assets/gallery/search-32-melbourne-convergence.jpg" alt="Excess tour length over the best known against evaluations for eight algorithms, log–log"></a><br><b>How fast each search closes in.</b> Median excess over the best known walk, with quartiles, against evaluations (10 seeds; the ant system 5). The ant system, given the distances as its heuristic, leads early and stagnates at 1.7 %; annealing ends best of the black boxes at 2.3 %; the positional Plackett–Luce law barely beats random tours.</td>
<td width="50%" valign="top"></td>
</tr>
</table>
