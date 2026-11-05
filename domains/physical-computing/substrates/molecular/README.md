# molecular: DNA Strand Displacement and Enzymatic DNA Networks

*Status: planned (ADR 0076).* Mode: programmed. **Breadth:** DNA strand displacement realises any reaction network up to its rate constants (Soloveichik, Seelig and Winfree 2010), and so, by [`chemical`](../chemical/)'s dual rails, any polynomial ODE. Provisional until molecular-science is founded (ADR 0076 §2.5).

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| DNA strand displacement: compiled networks (Soloveichik, Seelig and Winfree 2010), seesaw circuits (Qian and Winfree 2011) | any formal reaction network; digital circuits | DNA (C, H, N, O, P), biodegradable; its fuel consumed; synthesis energy-intensive; room temperature |
| Enzymatic DNA networks, the PEN toolbox (Montagne et al. 2011; Fujii and Rondelez 2013; Zadorin et al. 2017) | oscillators, bistability, predator–prey cycles, travelling waves, a French flag | DNA templates, a polymerase, a nickase, an exonuclease; fed by dNTPs |

## 2. The Simulators

1. **Compiled strand displacement.** A formal network is compiled by Soloveichik's scheme into species, fuel complexes and intermediates, with rates $q_i$ and fuel at $C_{max}$. It runs by mass action and as a Markov chain (physis `chemistry`). Checks:
   - as $C_{max} \to \infty$ with $q_i C_{max}$ held fixed, the trajectories converge to the formal network's;
   - the fuel consumed equals the reactions fired;
   - leaks are declared and switchable.
2. **Seesaw circuits.** Gates, thresholds and amplification. Check: the truth tables of AND and OR.
3. **PEN networks.** Templates; production by polymerase and nickase with saturable enzymes; degradation by an exonuclease with Michaelis–Menten kinetics, shared among species, so that load couples them; patterns by reaction–diffusion. Checks:
   - the oscillator's Hopf point by continuation agrees with the Routh–Hurwitz condition;
   - an autocatalytic front moves at Fisher's speed $2\sqrt{Dk}$.

## 3. Tasks and their Lowering

- **Functions and flows:** a formal network, from [`chemical`](../chemical/)'s programs or dual rails, is lowered to strands.
- **Phenomena:** oscillation and patterns (simulator 3).

## 4. What the Foundations Must Gain

physis `chemistry`: networks with Michaelis–Menten terms. Mass action and Markov chains are already written.

## 5. References

- Soloveichik, Seelig and Winfree (2010); Qian and Winfree (2011); Montagne et al. (2011); Zadorin et al. (2017): in the theory note.
- Fisher, R. A. (1937), The wave of advance of advantageous genes, *Annals of Eugenics* 7, 355–369.
- Fujii, T. and Rondelez, Y. (2013), Predator–prey molecular ecosystems, *ACS Nano* 7, 27–34.
