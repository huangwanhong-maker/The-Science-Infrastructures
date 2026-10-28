# 0052. Fields Known in Part: Potential Fields, the Moon's Magnetism, the Planetary-Science Domain, the Kernel Asked, and Instruments Declared

- **Status:** Accepted, provisionally. The unresolved part of a field is a stationary prior; the carriage of the solar wind is the planetary-science domain's until a space-weather bridge is founded.
- **Date:** 2026-10-28
- **Decided by:** the user asked Cosmopoiesis whether its scheduler was "RT kernel aware". They asked for the compass to be a sensor "only available near Earth's ground, or simply remove compass and only visualize the magnet fields". They asked for "a generative model … from historical data of moon magnet field etc. And input to virtual sensor instances". They asked for classes for the encyclopedia's instruments, which produce no data until written. And they asked whether phenomena and their visual expressions are well categorised, "if you think better design exist, just use your idea". Asked where the Moon's magnetism should live, since geophysics is the solid Earth, the user answered "planetary can founded". Claude designed it and wrote it, with this record.

  The work's record is Cosmopoiesis's [decision 0005](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0005-fields-known-in-part-and-declared-instruments.md). The foundations' records are calx's [decision 0010](https://github.com/huangwanhong-maker/calx/blob/main/development/decisions/0010-kernel-awareness-and-extended-range-legendre.md) and physis's [decision 0010](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0010-potential-fields.md). ADR 0051 is another work's.

## 1. Context

ADR 0050 gave Cosmopoiesis instruments as kernels of observation, read on the domains' models and scheduled with deadlines. The user's five requests needed general methods the foundations did not have:

- a question to the operating system's kernel;
- Legendre functions that survive high degree at the poles;
- potential fields evaluated at points and on grids, with their spectra;
- random fields of a given spectrum;
- a dataset of the Moon's crustal field;
- a fault fixed in the shared globe.

## 2. Decision

### 2.1 What Went Where

| Where | What |
|---|---|
| calx `scheduler` | `kernel.hpp`: the kernel's real-time classes probed on short-lived threads (SCHED_FIFO, SCHED_DEADLINE, RLIMIT_RTPRIO, PREEMPT_RT, WSL), reported with the reason for a refusal; the calling thread raised to SCHED_FIFO and lowered again |
| calx `special_functions` | Schmidt's Legendre functions in extended range (Fukushima 2012): to degree 2700 at $2\times10^{-6}$ rad of the poles, the sum rule to $n^2\varepsilon$ |
| physis `analysis` (new) | internal potential fields outside a sphere at points and on grids (by rows, $O(LN^2)+O(LMN)$), the Lowes–Mauersberger spectrum, SHTOOLS's format; spectral laws $A\,n^\beta\rho^{2n}$ fitted by least squares; isotropic Gaussian fields of a given spectrum; the part an estimate leaves unresolved, by the law of total variance. The convention for spherical harmonics, the semi-normalisation without the Condon–Shortley phase, is fixed (physis's open decision 2, in part) |
| the data library | `moon-magnetic-t2015`: T2015_449, the Moon's crustal field to degree 449 (Tsunakawa et al. 2015, expanded by Wieczorek, Zenodo, CC BY 4.0); `fetch.py` decompresses a file that arrives gzipped (`"gunzip": true`) |
| the visualization layer | `globe/lod.js`: a layer of one level (a whole map) seen close is asked of the drawn tiles' ancestor at that level; it was never loaded before |
| the planetary-science domain (new) | `magnetism`: the Moon's crustal field (T2015_449 read, the law of its spectrum on degrees 30–250 and its effective depth, the unresolved part drawn from a seed, the reading on low ground below the sphere) and the interplanetary field carried from L1 to a body near the Earth (the travel, the magnetotail, GSE's axes); its theory note; 18 checks, against pyshtools on T2015_449 and closed forms |
| Cosmopoiesis | the compass with its coverage; the Moon's magnetometer (the crust, its finer scales, the solar wind carried from L1); the field lines and the Moon's field map; twelve kernel classes and the encyclopedia's 101 families declared; the envelope of a reading with its quality; the kernel's report in each session; the visual expression of phenomena as restriction, reduction and mark (its formalism §3.6, §8.5, §8.7) |

### 2.2 Instruments: Kernel Classes and Declared Families

Each kernel of observation is now an abstract class, from which every family with a model derives. The encyclopedia's families are rows of a generated table, refused when mounted until a model is written. The user asked for classes that tell developers and users what can be implemented; a row says that, and is checked against the encyclopedia. The kernels stay Cosmopoiesis's until a second work needs them (ADR 0050).

### 2.3 The Planetary-Science Domain

The domain the README's §4 planned in tier 3, a bridge over astronomy, geophysics and meteorology for the planets and moons as worlds, is founded now at the user's request, at height 3 ([domains/planetary-science](../../domains/planetary-science/)). Its first module, `magnetism`, holds what names the Moon and its environment: T2015_449, its band and its lowest ground, and the carriage of the solar wind's field to it. Its methods are physis's and calx's. Cosmopoiesis's magnetometer is an instrument class of the work on the domain's model. The carriage of the solar wind and the Earth's magnetotail belong, strictly, to the space-weather bridge (tier 4): until it is founded they are this domain's, as the plasma environment of a body.

### 2.4 The Moon's Field Known in Part

The Moon's crust is known to degree 449 by an assimilation of orbital magnetometry. A spectral law fitted to it gives the power of the scales it does not resolve, and they are drawn from the session's seed, labelled synthetic. The solar wind's field, measured at L1, is carried to the Moon as physical forcing. No learned model is used: the crust's field is static, and a generator trained on one model would invent structure that is not a potential field's.

## 3. Alternatives Considered

- **Planetary magnetism in physis.** Rejected: T2015, the Moon's radius and the magnetotail's geometry name objects in nature. physis holds what any potential field needs.
- **The Moon's binding in the geophysics domain.** Rejected: geophysics is the solid Earth (ADR 0024); the user chose to found planetary science.
- **The binding kept in the work.** Rejected: by CLAUDE.md §3.7 what names an object in nature belongs to a domain, and a second work (the Celestial Atlas) may read the Moon's field.
- **The kernels in a foundation.** Kept in the work, as ADR 0050 decided, until a second work needs them.

## 4. Consequences

- **Fifteen domains have code.** Planetary science is C++ from the start, at height 3; it passes its checks alone with its down-set.
- **The checks.**
  - calx gains `check_calx_kernel`, and the extended range in `check_calx_legendre`.
  - physis gains `check_physis_analysis`, against pyshtools and 40-digit arithmetic (25 programs).
  - The visualization layer's 109, with the whole map seen close, a check that fails without the fix.
  - planetary science's 18, with T2015_449 against pyshtools at five mounts to $10^{-10}$ of $\lvert B\rvert$.
  - Cosmopoiesis's 175, with the Moon's magnetometer through the engine against the same fixture to $10^{-7}$, and on low ground.
- **Speed.** The Moon's field map on a half-degree grid crosses the WSL pipe as float32 in base64. It takes 2.1 s, down from 8.3 s as decimal digits, and the server keeps the last few.
