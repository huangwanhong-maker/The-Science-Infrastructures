### Quantum Search

Quantum search as heuresis states it: a device is a sampler of a law, emulated exactly, never simulated as a quantum computer, and its speed-up is counted in oracle calls; a quantum-inspired algorithm is a product law ([the record](../research/experiments/2026-11-05-quantum-search/README.md)).

<table>
<tr>
<td width="50%" valign="top"><a href="../assets/gallery/search-37-grover-law.jpg"><img src="../assets/gallery/search-37-grover-law.jpg" alt="The probability of measuring a marked point after r Grover iterations, against r and against the rotation angle, with emulated frequencies"></a><br><b>Grover's law and its emulator.</b> sin²((2r + 1)θ), sin²θ = t/N, on N = 4096 points, against the iterations and against the rotation angle, where every t/N falls on one curve. Rings are the frequencies of heuresis's emulator, 4000 measurements each: Σz² = 157.9 on 176 cells, p = 0.83.</td>
<td width="50%" valign="top"><a href="../assets/gallery/search-38-grover-adaptive-search.jpg"><img src="../assets/gallery/search-38-grover-adaptive-search.jpg" alt="Oracle calls of Grover adaptive search against the size of the space, beside Dürr and Høyer's bound and random search; threshold descents on a QUBO"></a><br><b>Grover adaptive search in oracle calls.</b> Mean oracle calls to the unique minimum over 200 runs at each size, from 7 at 2⁶ points to 2685 at 2²⁰, below Dürr and Høyer's bound and with a slope tending to ½, against random search's N. Right: the threshold's rank falling on a QUBO of 20 bits, the minimum in a median of 2758 oracle calls.</td>
</tr>
<tr>
<td width="50%" valign="top"><a href="../assets/gallery/search-39-quantum-inspired.jpg"><img src="../assets/gallery/search-39-quantum-inspired.jpg" alt="Heat maps of the probability of reading each bit as the optimum's for the quantum-inspired EA and PBIL, with convergence curves"></a><br><b>A quantum-inspired search is a product law.</b> Han and Kim's QIEA beside PBIL on the same QUBO: each bit's probability of reading as the optimum's, generation by generation. Its 1 600 000 observed bits equal the product-Bernoulli family's reference draws at sin²φ on the same keyed streams, without one exception.</td>
<td width="50%" valign="top"></td>
</tr>
</table>
