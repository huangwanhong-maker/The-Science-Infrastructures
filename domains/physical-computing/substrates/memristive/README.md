# memristive: Memristors, Crossbars, Mott Switches and Phase-Change Memories

*Status: planned (ADR 0076).* Modes: programmed, reservoir, learning. **Breadth:** networks of one kind of unit. Nearest to practice of all families. Ionic memristors are in [`soft-matter`](../soft-matter/), and mycelial ones in [`living`](../living/).

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Memristors (Chua 1971; Strukov et al. 2008) | pinched hysteresis; non-volatile memory of charge | TiO₂ (abundant); platinum electrodes (a critical raw material) |
| Crossbars (Ielmini and Wong 2018) | matrix–vector products by Ohm's and Kirchhoff's laws | as above; analogue precision of a few bits, with drift |
| Mott switches and oscillators: VO₂, NbO₂ (Pickett and Williams 2012; Kumar, Strachan and Williams 2017) | threshold switching, hysteresis, relaxation oscillation, chaos at the edge | vanadium and niobium, both critical raw materials |
| Phase-change memories: GeSbTe (Le Gallo and Sebastian 2020) | analogue weights by partial crystallisation; resistance drift | germanium and antimony, both critical raw materials |

## 2. The Simulators

1. **The HP memristor.** $v = (R_{on}w/D + R_{off}(1 - w/D))\,i$ and $\dot w = \mu_v R_{on} i/D$, with Joglekar's window $f(x) = 1 - (2x - 1)^{2p}$. Checks:
   - the loop is pinched: $v = 0$ exactly when $i = 0$;
   - without a window, the memristance is the closed form in the charge, $M(q) = R_{off} - (R_{off} - R_{on})\,\mu_v R_{on} q/D^2$;
   - the lobes' area tends to 0 as the frequency grows.
2. **The crossbar.** Conductances $G_{ij}$ are programmed with lognormal variation, read with noise and allowed to drift. The wires' resistance is solved as a Kirchhoff network (physis `discretisation`), and the outputs are read by ADCs of $b$ bits. Checks:
   - with no wire resistance, the outputs are $I_j = \sum_i G_{ij} V_i$ exactly;
   - the error from wire resistance is first order in that resistance;
   - quantisation errors are uniform within half a least significant bit;
   - the energy of a product comes from the declared device energies and the ADCs' figure of merit.
3. **The Mott oscillator.** A threshold switch $(V_{th}, V_{hold}, R_{on}, R_{off})$ is a hybrid system (physis `dynamics`, `hybrid`) in an RC circuit. Checks:
   - the period equals the closed form of the charging and discharging times;
   - two capacitively coupled oscillators lock in anti-phase.
4. **The Mott element at the edge of chaos** (Kumar, Strachan and Williams 2017). Joule heating, a conductance activated thermally and a parallel capacitance. Check: the sign of the largest Lyapunov exponent across the drive, against the paper's regimes.
5. **The phase-change memory.** The crystallised fraction follows $X = 1 - \exp(-(kt)^n)$ under SET pulses (Avrami 1939), the conductance follows from $X$, and the resistance drifts as $R(t) = R_0 (t/t_0)^\nu$. Checks: both closed forms.

## 3. Tasks and their Lowering

- **Functions** (simulator 2): linear algebra; a matrix is programmed into conductances, with differential pairs for signs.
- **Filters** (simulators 1, 3 and 4 as reservoirs).
- **Phenomena** (simulator 4).
- **Learning:** updates in place on the crossbar.

## 4. What the Foundations Must Gain

- physis `discretisation`: Kirchhoff networks on graphs;
- physis `dynamics`: Lyapunov spectra; its `hybrid` is already written.

## 5. References

- Chua (1971, 2014); Strukov et al. (2008); Kumar, Strachan and Williams (2017): in the theory note.
- Avrami, M. (1939), Kinetics of phase change. I. General theory, *Journal of Chemical Physics* 7, 1103–1112.
- Ielmini, D. and Wong, H.-S. P. (2018), In-memory computing with resistive switching devices, *Nature Electronics* 1, 333–343.
- Joglekar, Y. N. and Wolf, S. J. (2009), The elusive memristor: properties of basic electrical circuits, *European Journal of Physics* 30, 661–675.
- Le Gallo, M. and Sebastian, A. (2020), An overview of phase-change memory device physics, *Journal of Physics D* 53, 213002.
- Pickett, M. D. and Williams, R. S. (2012), Sub-100 fJ and sub-nanosecond thermally driven threshold switching in niobium oxide crosspoint nanodevices, *Nanotechnology* 23, 215202.
