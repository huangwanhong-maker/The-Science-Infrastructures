# signal: The Spectra of DNA

*Status: draft (2026-10-24).* A DNA sequence read as four signals, one for each base: their discrete Fourier transforms, the spectrum they sum to, and the period-3 signal by which protein-coding DNA shows itself. It serves learning, as the abelian case of the Fourier transform of a group, and use, as a reference-free measure of coding potential along a genome. The theory is [theory/signals-classification-and-comparison.md](../theory/signals-classification-and-comparison.md), §1.

## 1. Concepts

- **Four indicator sequences.** u_b(n) = 1 where the base at n is b, 0 elsewhere (Voss 1992). Their transforms U_b(k) over ℤ_N are computed by calx's FFT (`calx/transforms/fft.hpp`), and S(k) = Σ_b |U_b(k)|² is the spectrum. Parseval makes its mean over the frequencies exactly the number of bases; it is the same for a sequence and its reverse complement.
- **The DFT is the Fourier transform of the cyclic group.** The frequencies k are the characters n ↦ e^{−2πi nk/N} of ℤ_N. When 3 divides N, the character of k = N/3 factors through ℤ_N → ℤ_3, so U_b(N/3) is a transform on ℤ_3 of the counts of b at the three positions of codons: an Eisenstein integer, whose norm is an integer. The period-3 power is therefore computed exactly from twelve counts, with no transform, and slides along a genome in O(1) a base.
- **Coding DNA has period 3.** Codons use the bases unequally at their three positions, and S has a peak at N/3 (Tiwari et al. 1997; Anastassiou 2001). Its strength P = S(N/3)/S̄ is near 1 − Σ_b p_b² in random sequence (3/4 at equal frequencies) and larger in genes. On *Escherichia coli* K-12, windows of 120 bases lying in CDSs have a mean strength of 1.69, against 0.74 between genes, the value for random sequence of their composition.

## 2. Content

| Header | Content |
|---|---|
| [`spectrum.hpp`](include/genomics/signal/spectrum.hpp) | `indicator(x, b)`; `indicator_transforms(x)` → U_b(k) for the four bases (calx's FFT); `spectrum(x)` → S(k); `CodonCounts{c}` with `bases()` and `power()` (S(N/3) as Σ_b of Eisenstein norms); `codon_counts(x, begin, length)`, `period3_power(x, begin, length)`; `period3_profile(x, window, step)` → `Period3{begin, power, bases}` with `strength()` |

## 3. Dependencies

`genomics::sequence` (codes); calx `transforms` (the FFT). The checks use `genomics::simulate` and calx `random`.

## 4. Checks

| Check | What it verifies |
|---|---|
| `spectrum` (21) | U_b(k) by the FFT against the direct sum over ℤ_N, at N = 1024 (radix 2) and 999 (Bluestein), with bases N, to 10⁻¹⁶ of N; S(k) against the direct sum; Parseval (the mean of S the number of bases); S(N − k) = S(k); the reverse complement's spectrum the same. (ACG)^400 peaking at N/3 exactly with S(N/3) = N²/3 and nothing else in (0, N/2]; random words of period 5, 7, 10, 21 supported on the multiples of N/p. The exact period-3 power against S(N/3) by the FFT on 50 random stretches; the sliding profile against counts made afresh, for windows of 3 to 351 and steps 1, 3, 7. The mean strength of 25 000 disjoint windows of random sequence against 1 − Σ_b p_b² at GC 0.3, 0.5 and 0.7, within 4 standard errors |
| `ecoli` (4) | When the cache has NC_000913.3: windows of 120 every 3 bases, 1 188 315 in CDSs and 69 170 touching no gene: the mean strength in CDSs more than twice that between genes (measured 2.26 times), a CDS window above an intergenic one in more than 0.83 of pairs (AUC measured 0.857), intergenic windows within 0.03 of 1 − Σ p² of their composition (0.744 against 0.741); windows of 351: ratio above 4 and AUC above 0.96 (measured 4.57 and 0.972). Skipped without the cache |

## 5. References

- Voss, R. F. (1992), Evolution of long-range fractal correlations and 1/f noise in DNA base sequences, *Physical Review Letters* 68, 3805–3808.
- Tiwari, S., Ramachandran, S., Bhattacharya, A., Bhattacharya, S. and Ramaswamy, R. (1997), Prediction of probable genes by Fourier analysis of genomic sequences, *CABIOS* 13, 263–270.
- Anastassiou, D. (2001), Genomic signal processing, *IEEE Signal Processing Magazine* 18(4), 8–20.
- Cooley, J. W. and Tukey, J. W. (1965) and Bluestein, L. I. (1970), as calx's FFT cites them.

## 6. Open Questions

1. **The non-abelian transform.** The DFT is the Fourier transform of ℤ_N; a group acting on the bases or on codons (the Klein group of complementarity and transition, or the symmetries of the genetic code) has a Fourier transform over its irreducible representations, of which this is the abelian case. Its general and exact forms belong to morphe's `algebra`; the spectra of a genome under them would join here.
2. **A gene finder from the signal.** The strength separates long coding windows well and short ones poorly (AUC 0.74 at 60 bases); a detector would combine it with the codon statistics of `annotation/genes.hpp`, or take the frame from the phase of U_b(N/3).
3. **Other mappings.** A complex or numerical encoding of the bases, x̂(n) = Σ_b w_b u_b(n) (Anastassiou's), has the transform Σ_b w_b U_b(k), a linear function of the four computed here; such spectra can be added from them when a use calls for one.
