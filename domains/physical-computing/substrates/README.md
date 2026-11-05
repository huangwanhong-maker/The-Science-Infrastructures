# substrates: The Backends, One Family per Directory

*Status: planned (ADR 0076 §2.4).* Each directory is a family of substrates that share their laws and their lowering, and each family is a target of the pipeline. Its first backend is a **simulator**; later backends are a **device**, through `observation/`, and a learned **surrogate** ([`runtime/`](../runtime/)).

## 1. The Families

| Family | Substrates | First simulator to write |
|---|---|---|
| [`analogue-electronics`](analogue-electronics/) | analogue computers, programmable analogue arrays, CMOS oscillators and oscillator-based Ising machines | a netlist of integrators, summers and multipliers |
| [`spintronics`](spintronics/) | spin-torque oscillators; magnetic networks, artificial spin ices; stochastic tunnel junctions | Slavin and Tiberkevich's auto-oscillator |
| [`memristive`](memristive/) | memristors and crossbars; Mott switches and oscillators; phase-change memories | the HP memristor; the crossbar |
| [`probabilistic`](probabilistic/) | p-bits and probabilistic computers; thermodynamic sampling and linear algebra | a network of p-bits |
| [`photonics`](photonics/) | interferometer meshes; delay-loop reservoirs; coherent Ising machines | the delay loop |
| [`chemical`](chemical/) | the BZ medium and its arrays; reaction networks as programs; electrochemical oscillators | the Oregonator |
| [`molecular`](molecular/) | DNA strand displacement; enzymatic DNA networks | a formal network compiled to strand displacement |
| [`fluids`](fluids/) | water in a tank; convection cells; rotating, stratified tanks; Faraday waves, walking droplets; capillary and foam networks; bubble logic; Phillips's machine | the convection cell and its Lorenz map |
| [`soft-matter`](soft-matter/) | ionic membranes, hydrogels; liquid crystals, active nematics; active colloids | the liquid-crystal cell |
| [`mechanical`](mechanical/) | metamaterials and hysterons; metronomes; soft bodies; granular media and landscapes | interacting hysterons |
| [`learning-networks`](learning-networks/) | resistor, flow and spring networks that learn | coupled learning in a resistor network |
| [`living`](living/) | mycelium; *Physarum*; ecological microcosms; gene circuits; biofilms; neural cultures and organoids | the repressilator and the toggle switch |
| [`quantum`](quantum/) | quantum reservoirs | a transverse-field Ising reservoir |
| [`swarms`](swarms/) | robot swarms | gradients on the swarm's graph |

## 2. The Contract of a Simulator

Each simulator, in every family, has these parts:

- **A model** with its source and equation numbers (CLAUDE.md §2.4), its state, its parameters in SI units with their sources, and its range of validity.
- **Noise and non-idealities**, declared one by one, each switchable, each drawn from calx's keyed streams.
- **Its own error**: the integrator and step it uses, and its error, kept apart from any task's tolerance.
- **Checks** against closed forms or published results, each named in its family's README before the code is written.
- **A target description** that the [compiler](../compiler/) reads: the kinds of task the substrate accepts, its dimensions and ranges, its speed, its energy per operation and its materials.
- **Laws from the foundations.** The directory holds only what names its substrates. A substrate whose material's domain is only planned holds the least model it needs, marked provisional, with the domain it will move to (ADR 0076 §2.5).

## 3. The README of a Family

Each family's README states:

1. its substrates, their phenomena and their materials;
2. its simulators, each with its model, its noise and non-idealities, what it stands on and its checks;
3. the tasks it accepts and how they are lowered;
4. what the foundations must gain for it;
5. the references it adds to the theory note's, all written from memory and to be checked.
