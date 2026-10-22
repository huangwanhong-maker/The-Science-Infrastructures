# 0019. A Second Laboratory Model: the Shallow-Water Equations

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review. Extends [ADR 0015](0015-earth-atlas-models-imagery-detail.md) §2.3.
- **Date:** 2026-10-22
- **Decided by:** the user asked for weather prediction models. Claude chose the design below.

## 1. Context

ADR 0015 gave the laboratory one model of its own, the barotropic vorticity equation. It named the shallow-water equations as the natural next step, with Williamson et al.'s test suite to verify them.

The barotropic model's main fault is known: it lacks divergence, so its longest waves run westward too fast. A model of our own should grow step by step, as numerical weather prediction did, with each step verified.

## 2. Decision

1. **The shallow-water model.**
   - The shallow-water equations on the sphere, in vorticity–divergence form, by the spectral transform method at T42.
   - Fourth-order Runge–Kutta steps of 15 minutes, and ∇⁴ damping (`src/core/shallowwater.js`).
   - It starts from the analysed wind and height at 500 hPa.
2. **Initialisation by a digital filter.** Before the forecast, a digital filter over six hours removes the gravity waves of an unbalanced start (Lynch and Huang 1992).
3. **Two models, compared.**
   - Each run becomes a product of its own, `lab-barotropic` or `lab-shallow-water`. The two can be drawn, differenced and read like any other models.
   - The panel runs either model and reports the error of each beside persistence, against the same analyses.
   - The worker that runs them is renamed `labmodel.worker.js`.
4. **Fast Fourier transforms.**
   - The spectral transforms use a radix-2 FFT when the number of longitudes is a power of two.
   - This halves the cost of both models and changes no result beyond rounding.
5. **Tests.**
   - Williamson et al.'s case 2, about the pole and across the poles, exact to 10⁻¹⁴.
   - Conservation of mass and energy.
   - The period of gravity waves on a non-rotating sphere (Lamb 1932).
   - The filter's removal of fast waves and keeping of slow ones.
   - The FFT against the direct sum.
6. **Theory.** The theory is in [theory/meteorology/shallow-water-model.md](../../theory/meteorology/shallow-water-model.md).

## 3. Alternatives Considered

- **A semi-implicit time step**, the standard for operational spectral models. It allows steps several times longer, but it couples the equations through a Helmholtz problem and is harder to read. At T42 the explicit step costs about twenty seconds for three days. Clarity is worth more here than speed.
- **Normal-mode initialisation** (Machenhauer 1977). Exact for the linear modes, but it needs the model's normal modes. The digital filter needs only the model itself, and it can be tested on single waves.
- **A multi-level primitive-equation model.** The real next step toward baroclinic dynamics. It needs vertical structure, thermodynamics and many more data. It is recorded as the next idea, not attempted yet.

## 4. Consequences

- **A better model of our own.** From the analysis of 2026-10-19 00:00 UTC, the shallow-water model's 500 hPa error over 20–80° N was:
  - 43 m at a day, against 58 m for the barotropic model;
  - 74 m against 85 m at two days;
  - 100 m against 108 m at three days.

  Persistence was 81, 108 and 127 m.
- **Tests:** the Earth Atlas has 105 checks.
- **Running time:** a shallow-water forecast of three days takes about twenty seconds in the browser, shown by a progress bar.
