# schedule: When a Program Runs

*Status: draft (2026-11-06): sample-and-hold, washout and repetitions, `include/physical_computing/schedule/schedule.hpp`.* A schedule says when inputs are applied and states are read. It is kept apart from the program, as Halide keeps a schedule apart from its algorithm (Ragan-Kelley et al. 2013), so that one program can be run under many schedules and each compared on its cost.

## 1. What a Schedule Holds

- **Injection.** When each input is applied, and for how long it is held.
- **Washout.** The time allowed for the initial state to be forgotten before readings count (the echo state property, theory note §1.2).
- **Virtual nodes.** Their number and spacing in a delay-loop reservoir.
- **Readings.** Their times and rates, and the conversion each needs.
- **Annealing.** The inverse temperature of p-bits, or the pump of a coherent Ising machine, as a function of time.
- **Data.** The split into training, validation and test sets.
- **Repetitions.** Repeated runs and their seeds, for the statistics of chaos and noise.
- **Clocks.** The ratio of the substrate's clock to the task's, from the compiler's scaling.

## 2. The Cost of a Schedule

Every schedule reports its duration and, with the backend's power, its energy. Together with the program, these give the cost of a verified result (theory note §1.5).

## 3. References

- Ragan-Kelley, J. et al. (2013), Halide: a language and compiler for optimizing parallelism, locality, and recomputation in image processing pipelines, *PLDI 2013*, 519–530.
