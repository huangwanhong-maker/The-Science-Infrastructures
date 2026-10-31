# physiology

*Status: planned. Field, height 2, on cell biology and neuroscience.* Organisms as regulated systems: excitable tissues and muscle, the heart and the circulation, breath and the exchange of gases, the feedback that holds the body's variables in range, immunity, the energetics of the whole body, and the physiology of plants. Immunology, endocrinology, nutrition and plant physiology are its modules. Placed in the order with the other domains of life at the user's request (2026-10-31, [ADR 0061](../../development/decisions/0061-the-domains-of-life.md)); its work is in tier 4 of the [README's order](../../README.md#46-the-order-of-implementation).

## 1. Purpose

- **Homeostasis is feedback.** A regulated variable, such as glucose in the blood, the body's temperature or arterial pressure, is held by a loop of sensor, controller and effector. Its models are control systems, whose gains, delays and nonlinearities set how closely and how stably the variable is held: Bergman's minimal model of glucose and insulin, the baroreflex.
- **Excitable tissue is channels on a cable.** The heart's cells fire as neurons do, by voltage-gated currents, but with plateaus and pacemakers of their own (Noble 1962). The wave they carry through the myocardium is reaction–diffusion on a tissue (the monodomain equation). The brain is neuroscience's; the body's excitable tissues are here, on the same equations of physis.
- **Circulation and breath are flows on networks.** The vessels form a network of conductances (Poiseuille's law; Murray's law of branching), and the arteries' compliance makes a Windkessel; this is physis's flow with conductances on a graph. Oxygen binds haemoglobin cooperatively (Hill, Adair) and crosses the alveolar membrane by diffusion.
- **Size sets rate.** Metabolic rate scales as a power of body mass near 3/4 (Kleiber 1932); networks of supply give one account of the exponent (West, Brown and Enquist 1997).
- **Plants regulate too.** Photosynthesis is limited by Rubisco or by electron transport (Farquhar, von Caemmerer and Berry 1980); stomata trade carbon for water; water rises through a plant down a gradient of water potential, the soil–plant–atmosphere continuum.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `excitable` | Cardiac and smooth-muscle cells: Noble's model and its successors, pacemaking, the monodomain equation on tissue, conduction and its block; muscle contraction by Huxley's cross-bridges | planned |
| `circulation` | The heart as a pump (time-varying elastance), Windkessel models of the arteries, vascular networks as conductances, Murray's law, the delivery of oxygen | planned |
| `respiration` | Ventilation and gas exchange: the alveolar gas equation, oxygen–haemoglobin binding (Hill, Adair), diffusion across the membrane | planned |
| `regulation` | Homeostasis as control: glucose and insulin (Bergman's minimal model), endocrine axes with their delays, thermoregulation, the baroreflex | planned |
| `immunity` | The immune response as dynamics: clonal expansion and contraction, affinity maturation; the host's side of an infection, whose pathogen is microbiology's | planned |
| `energetics` | Metabolic rate and its scaling, energy balance and nutrition, body temperature | planned |
| `plants` | Photosynthesis (Farquhar–von Caemmerer–Berry), stomatal conductance, transpiration, water transport along the soil–plant–atmosphere continuum | planned |

## 3. Dependencies

| Source | What physiology uses from it |
|---|---|
| cell-biology | Membranes and their potentials, signalling, metabolism |
| neuroscience | The nervous control of organs, and the models of excitable cells it shares, through its public interface |
| physis | `dynamics` (hybrid systems; delays, planned), `chemistry` (reaction–diffusion), `discretisation` (flows with conductances on graphs), `continuum` (diffusion; later the flow of fluids) |
| calx | Integrators (the θ-method for stiff systems), linear algebra, optimisation |
| tyche | Fitting models to physiological records |

## 4. Used By

No work uses it yet. Ecology will declare it when its `ecosystems` module needs photosynthesis, or the metabolic theory of ecology.

## 5. Open Questions

- Where neuroscience ends and physiology begins: the brain, its maps and its signals there; the heart, muscle, gut and the autonomic control of organs here. Whether the peripheral nervous system is neuroscience's.
- The host and the pathogen: immunity here, infection within a host in microbiology, spread between hosts in epidemiology; which bridge couples them.
- Pharmacokinetics is drug discovery's (pharmacology, README §4.5) until it graduates onto physiology.
- Anatomy, the body's structures as images and models, is a domain of its own, on which physiology may come to stand.

## 6. References

- Keener, J. and Sneyd, J. (2009), *Mathematical Physiology*, 2nd ed., 2 vols, Springer.
- Hall, J. E. and Hall, M. E. (2021), *Guyton and Hall Textbook of Medical Physiology*, 14th ed., Elsevier.
- Taiz, L., Zeiger, E., Møller, I. M. and Murphy, A. (2015), *Plant Physiology and Development*, 6th ed., Sinauer.
- Noble, D. (1962), A modification of the Hodgkin–Huxley equations applicable to Purkinje fibre action and pace-maker potentials, *Journal of Physiology* 160, 317–352.
- Huxley, A. F. (1957), Muscle structure and theories of contraction, *Progress in Biophysics and Biophysical Chemistry* 7, 255–318.
- Bergman, R. N., Ider, Y. Z., Bowden, C. R. and Cobelli, C. (1979), Quantitative estimation of insulin sensitivity, *American Journal of Physiology* 236, E667–E677.
- Westerhof, N., Lankhaar, J.-W. and Westerhof, B. E. (2009), The arterial Windkessel, *Medical & Biological Engineering & Computing* 47, 131–141.
- Murray, C. D. (1926), The physiological principle of minimum work: I. The vascular system and the cost of blood volume, *PNAS* 12, 207–214.
- Kleiber, M. (1932), Body size and metabolism, *Hilgardia* 6, 315–353.
- West, G. B., Brown, J. H. and Enquist, B. J. (1997), A general model for the origin of allometric scaling laws in biology, *Science* 276, 122–126.
- Farquhar, G. D., von Caemmerer, S. and Berry, J. A. (1980), A biochemical model of photosynthetic CO₂ assimilation in leaves of C₃ species, *Planta* 149, 78–90.
