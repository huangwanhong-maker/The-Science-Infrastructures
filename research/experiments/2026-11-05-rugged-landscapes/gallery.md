# Rugged Landscapes: the Pictures

Kauffman's NK landscapes of twenty bits, from smooth (k = 0) to uncorrelated (k = 19), each counted whole, all 2²⁰ points, by heuresis's `landscapes` module, and seven of heuresis's searches against them ([the record](README.md)). The pictures are drawn by `figures.py` from `results/`.

<table>
<tr>
<td width="50%" valign="top"><a href="../assets/gallery/search-27-ruggedness.jpg"><img src="../assets/gallery/search-27-ruggedness.jpg" alt="Six panels of ruggedness against k for NK landscapes of twenty bits"></a><br><b>From smooth to rugged, counted whole.</b> The figure shows, against k, five measures of sixty NK landscapes (five per k): the local optima under one-bit flips, the one-step correlation ρ(1), the correlation length, the fitness–distance correlation and the share of the space in the global optimum's basin. A sixth panel gives the autocorrelation along a random walk. At k = 19 the landscape is uncorrelated, and its 49 830–50 119 local optima match 2²⁰/21 = 49 932. Weinberger's ρ(1) = 1 − (k+1)/n, corrected for tables of only 2^(k+1) entries, meets the exact values within 0.004.</td>
<td width="50%" valign="top"><a href="../assets/gallery/search-28-performance-against-k.jpg"><img src="../assets/gallery/search-28-performance-against-k.jpg" alt="Gap, rank and success of seven searches against k at two budgets"></a><br><b>Seven searches against ruggedness.</b> Random search, the (1+1) EA, a GA, simulated annealing, tabu search, the compact GA and UMDA each had 5 000 or 50 000 evaluations. Each panel is over fifty runs per k: the gap to the enumerated optimum, the rank of the best point among the 2²⁰, and the share of runs that found the optimum. Tabu search leads at 5 000 evaluations from k = 4 to 10. Annealing on the Lam schedule finds the optimum in every run up to k = 4 at 50 000. At k = 19 nothing does clearly better than random search.</td>
</tr>
<tr>
<td width="50%" valign="top"><a href="../assets/gallery/search-29-local-optima-and-basins.jpg"><img src="../assets/gallery/search-29-local-optima-and-basins.jpg" alt="Scatter plots of local optima by fitness and basin size for k = 2, 6 and 10"></a><br><b>Local optima and their basins.</b> Every local optimum of three landscapes (instance 0 at k = 2, 6 and 10: 56, 1 577 and 7 351 optima) is placed by its fitness and the size of its basin under steepest descent, and coloured by its Hamming distance to the global optimum (the star). Fitter optima have larger basins (correlation of F with the logarithm of the basin size 0.86–0.89), but from k = 4 on the optima lie at random distances from the best.</td>
<td width="50%" valign="top"></td>
</tr>
</table>
