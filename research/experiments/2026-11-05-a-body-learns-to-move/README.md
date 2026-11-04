# A Body Learns to Move: the Rabbit's Controller Tuned by Search

*Status: draft (2026-11-05).* Research: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) on the robotics domain's rabbit ([ADR 0058](../../../development/decisions/0058-a-body-with-joints.md)).

## 1. The Question

The robotics domain's rabbit hops on the Moon under a controller of Raibert's three parts, whose gains were set by hand for hopping at 0.3 m/s. Two questions follow. How well do the general searches of heuresis tune those gains for a task the defaults cannot do, and how do they compare? And what is the trade between speed and the power drawn when the gains and the speed asked are searched together?

This is a black-box problem of the kind heuresis is for. Each evaluation is an 8.5-s simulation of a body with compliant contact, and it has no gradient. It fails outright in a third of the box: the rabbit falls or tips over. Its objectives are measured, not modelled. The body's physics, its refusals and its failures are the domain's. The search sees only what a trial returns.

## 2. Method

### 2.1 The Body, Its Controller, and the Loop That Moves Them

The body is `robotics::legged::rabbit()`: 1.69 kg, six links on a floating base, five geared motors, an Achilles tendon and compliant feet on regolith. It is advanced by semi-implicit Euler at 0.5 ms, the cheapest of the domain's three methods. The controller, `robotics::legged::Controller`, runs at 500 Hz. Gravity is the Moon's 1.62 m/s², the domain's default.

The domain's `Simulation` builds its controller with the default `Gains`: its `SimulationConfig` has no field for them, so a search cannot tune the controller through it (§8). The experiment therefore moves the body with its own loop, [hopper.hpp](hopper.hpp). It is the `Simulation`'s semi-implicit path step for step, built from the domain's public parts (`Rabbit`, `RabbitForces`, `potential_force`, `Controller`, `GroundBelief`) and physis's `semi_implicit_euler`. The controller is sampled every 2 ms with its torques held, and the motors' electrical energy is counted over each held period.

**The loop reproduces the domain.** With the default gains, after a second of standing and nine of hopping at 0.3 m/s, its state equals the domain's `Simulation` bit for bit, with the same electrical energy, 41.526742079 J (`rabbit check`, `results/check.json`).

**The feet are not read as sensors.** The trials leave out the feet's readings and the ground's belief that they feed. The controller asks at most 0.5 of friction, and its belief of the ground is presumed to be μ ≈ 0.69, so three quarters of that belief (0.52) never binds. The motion without the readings is again bit for bit the same, and the trials are cheaper.

### 2.2 A Trial

The rabbit stands for 0.5 s and is then commanded to hop at a speed and an apex. It is simulated for 8 s more. It is measured over the last 4 s: the mean speed, the electrical energy and its copper part, the cost of transport E/(m g d), and the mean apex at each landing. The largest pitch is taken over the whole hop. A trial is a `heuresis::search::failure` with its reason when:

- **fell**: the torso's corners or head touched the ground;
- **tipped over**: the pitch went beyond 60°;
- **refused**: the controller declined the command;
- **diverged**: the state was not finite;
- **stopped hopping**: fewer than three landings in the window.

Its cost is the simulated time.

### 2.3 The Gains Searched

Ten fields of `robotics::legged::Gains`, each mapped from [0, 1]:

| Gain | Range | Default | Best found (CMA-ES, seed 2005) | Median over the 32 runs of the four searches |
|---|---|---|---|---|
| leg_stiffness (N/m) | 200 – 800 | 420 | 440 | 477 |
| leg_damping (N s/m) | 0 – 10 | 3 | 3.1 | 4.2 |
| leg_rest (m) | 0.17 – 0.215 | 0.20 | 0.215 | 0.211 |
| hop_ankle (°) | 55 – 85 | 70 | 78.4 | 80.3 |
| attitude_kp (N m/rad) | 0.2 – 3 | 1.0 | 0.80 | 0.97 |
| momentum_rate (1/s) | 5 – 80 | 30 | 26.8 | 36.5 |
| placement (s) | 0.02 – 0.2 | 0.09 | 0.077 | 0.067 |
| speed_integral | 0 – 1 | 0.4 | 0.36 | 0.52 |
| height_gain (N/m) | 5 – 60 | 20 | 12.4 | 24.1 |
| thrust_max (N) | 4 – 20 | 12 | 18.9 | 11.1 |

The medians are those of each run's best feasible gains, over the runs of the (1+1)-ES, CMA-ES, DE and PSO that found one.

### 2.4 One Task, Five Searches

The task is to hop at 0.6 m/s with an apex of 12 cm, at the least cost of transport. The constraints, each ≤ 0 when held, are:

- the speed within 3 cm/s of 0.6 m/s;
- the apex within 2 cm of 12 cm;
- the pitch within 20°.

Deb's order compares the trials: feasible before infeasible, then by cost or by total violation, failures last. The default gains reach 0.59 m/s at a cost of 1.89, but pitch to 55° as they gather speed, which is infeasible.

Five searches run on the unit box, each driven by `heuresis::search::run`, 8 seeds of 800 trials:

| Search | heuresis | Settings |
|---|---|---|
| random search | `search::RandomSearch` | uniform on the box |
| (1+1)-ES | `strategies::OnePlusOneEs` | a uniform start, σ₀ = 0.25, the one-fifth rule, reflection at the box |
| CMA-ES | `strategies::Cma` | a uniform start, σ₀ = 0.25, λ = 10, active, the box by its transformation |
| DE | `strategies::DifferentialEvolution` | rand/1/bin, 20 members, F = 0.5, CR = 0.9 |
| PSO | `swarms::ParticleSwarm` | the inertia rule (w = 0.7298, c₁ = c₂ = 1.49618), a ring, 20 particles |

### 2.5 Two Objectives

The second search adds the speed asked, 0.1 to 1.0 m/s, as an eleventh coordinate, with the apex asked held at 12 cm. It minimises the speed reached (negated) and the mean electrical power. Its constraints are the apex within 2 cm and the pitch within 20°. It runs `multiobjective::Nsga2` and `multiobjective::SmsEmoa`, μ = 40, with SBX (η = 15) and polynomial mutation (η = 20), 3 seeds of 2 000 trials. Every feasible trial is kept; the union's non-dominated set is the front. The hypervolume, by `multiobjective::hypervolume`, is taken against the reference (0 m/s, 6 W).

The cost of transport is not an objective. Copper loss, (τ/k_τ)² R, is 78 to 94 % of the power drawn in every gait measured here, and it changes little with speed. So the cost per metre falls with speed up to the gait's limit: with the default gains it is 10.4 at 0.09 m/s and 2.24 at 0.51 m/s. A front of speed against it collapses to its fastest point. Speed against power is the real trade: the battery's drain against the distance it buys.

## 3. Code Version

- The infrastructure at 4153cb9 (heuresis is a submodule there).
- heuresis 6b3cc88, calx 90aa710, morphe 14b3fb8, physis e2b98d8, tyche b4034b8.
- The robotics domain as committed at 4153cb9.

To build and run (in WSL, g++ 11.5 and CMake 3.28), from this directory:

```
cmake -S . -B ~/.cache/science-build/exp-rabbit -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-rabbit -j4
B=~/.cache/science-build/exp-rabbit/rabbit
$B check results; $B explore results; $B tune results 8 800; $B front results 3 2000
.venv/Scripts/python.exe figures.py      # from Windows, the repository's environment
```

The machine was shared with other work while these ran (a load average of about 50 on eight cores). The 40 tuning runs took 510 s of wall time and 3 944 s of their threads, 0.12 s a trial on average (trials that fall end early); the six two-objective runs took 120 s. A full trial of 8.5 s took 0.2 s of a thread under that load.

## 4. Parameters

| | |
|---|---|
| Body | `robotics::legged::rabbit()`, the build of ADR 0058 |
| Site | the Moon, g = 1.62 m/s²; the rabbit's own ground law (k = 3 × 10⁵ N/m^{3/2}, μ = 0.7) |
| Method | semi-implicit Euler, step 0.5 ms; controller 2 ms |
| Trial | 0.5 s standing, 8 s hopping, measured over the last 4 s |
| Task | 0.6 m/s, apex 0.12 m; tolerances 3 cm/s, 2 cm, 20° |
| Tuning | 8 seeds (2000 to 2007) × 800 trials × 5 searches |
| Front | 3 seeds (3000 to 3002) × 2 000 trials × 2 searches; μ = 40 |

## 5. Data

No external data. The body, its motors and its ground are the robotics domain's ([theory: legged](../../../domains/robotics/theory/legged.md)); the Moon's gravity is its default. The tables in `results/` are this study's own:

| File | Content |
|---|---|
| `check.json` | the loop against the domain's `Simulation` |
| `explore.json` | the default gains over a grid of commands |
| `tune.json` | the 40 tuning runs, their checkpoints, failures, best gains and the measures of the best |
| `tune_motions.json` | the poses of the default and tuned gains |
| `front.json` | the six two-objective runs, their hypervolumes and fronts, and the union's front with its gains |
| `front_motions.json` | the poses of four points of the front |
| `*.log` | the program's reports |

## 6. Results

### 6.1 What the Defaults Do

With the default gains and an apex of 12 cm asked, the rabbit tracks the speed asked up to 0.5 m/s, reaching 0.51 m/s at 3.14 W. At 0.6 m/s it reaches 0.59 m/s but pitches nose-down to 55° in its first hops as it gathers speed, then recovers. At 0.7 and 0.8 m/s it tips over. Hopping in place it draws 2.50 W, of which 2.28 W is copper loss.

### 6.2 Five Searches on One Task

| Search | Final best cost of transport: median | best | worst | Feasible trials | Fell | Tipped over |
|---|---|---|---|---|---|---|
| random | 1.655 | 1.494 | 1.742 | 3.6 % | 32.4 % | 36.1 % |
| (1+1)-ES | 1.466 | 1.396 | 1.569 | 78.5 % | 1.3 % | 3.0 % |
| CMA-ES | 1.441 | 1.362 | none on 1 seed | 51.1 % | 12.8 % | 10.0 % |
| DE | 1.480 | 1.454 | 1.511 | 44.9 % | 3.2 % | 11.8 % |
| PSO | 1.494 | 1.446 | 1.551 | 41.1 % | 4.7 % | 15.2 % |

- **Every search beats random search**, and every seed of the four but one of CMA-ES's reaches a cost of 1.57 or less, against random search's median of 1.66. Most of the box is a body on the ground: 69 % of the random trials fell or tipped over, and 3.6 % were feasible.
- **The (1+1)-ES is the fastest to a good gait.** Its median is 1.47 by 200 trials, and it then stalls. It stays near what works: 79 % of its trials are feasible.
- **CMA-ES ends with the best median, 1.441, and the best gait of all, 1.362.** It is slowest to start, its first generations spread over σ₀ = 0.25 of the box. On one seed (2006) it began in a region where all 800 trials fell or tipped over. Its ranking there saw only failures, tied, so it carried no direction, and the search never left.
- **DE and PSO** are the most consistent (DE's eight seeds span 1.454 to 1.511) but end higher than CMA-ES and the (1+1)-ES.
- **No hop was refused**, none diverged, and none stopped.

The best tuned gains, measured again, reach 0.615 m/s at a cost of transport of 1.362, 28 % below the default's 1.89 at 0.59 m/s. They draw 2.29 W against 3.05 W, with an apex of 11.0 cm and a pitch within 7.4°, against the default's 55°.

**The gains the searches agree on.** Across the 32 runs of the four searches:

- **a longer virtual leg**: 0.211 m in the median against 0.20, and the best at the bound, 0.215;
- **the ankle stretched further in a hop**: 78.6 to 83.3° in the middle half, against 70°. The tendon is slack at 75°;
- **a smaller placement gain**: 0.059 to 0.071 s, against 0.09.

The other seven gains spread over their ranges. The cost is flat in them near the optimum, or they trade against one another.

### 6.3 Speed Against Power

The six runs' union front holds 149 gaits, 91 of them from NSGA-II and 58 from SMS-EMOA. They run from 0.124 m/s at 1.42 W to 0.988 m/s at 3.80 W. At the speeds where the defaults hop, the front draws about a third less power:

| Speed (m/s) | Default gains (W) | Front (W) | |
|---|---|---|---|
| 0.19 | 2.45 | 1.48 | −39 % |
| 0.30 | 2.50 | 1.67 | −33 % |
| 0.40 | 2.76 | 1.81 | −34 % |
| 0.51 | 3.14 | 2.08 | −34 % |

The front reaches 0.99 m/s, which the default gains tip over at. Its cost of transport falls from 4.18 at 0.12 m/s to 1.33 at 0.71 m/s, then rises to 1.40 at 0.99 m/s, where the pitch reaches the 20° bound.

The final hypervolumes are 2.99, 3.58 and 3.62 for NSGA-II, and 3.26, 3.39 and 3.53 for SMS-EMOA. The six runs end between 2.99 and 3.62; three seeds cannot rank the two searches. Of each run's trials, 14 to 16 % tipped over or fell.

The four gaits drawn (search-18):

| Asked (m/s) | Reached (m/s) | Power (W) | Cost of transport | Apex (cm) | Pitch (°) |
|---|---|---|---|---|---|
| 0.13 | 0.124 | 1.42 | 4.18 | 10.9 | 7.7 |
| 0.42 | 0.413 | 1.82 | 1.61 | 10.5 | 7.0 |
| 1.00 | 0.706 | 2.58 | 1.33 | 10.5 | 6.7 |
| 0.99 | 0.988 | 3.80 | 1.40 | 13.1 | 19.6 |

The speed asked is not the speed reached: at 0.71 m/s the gait asked for 1.0 m/s, and the gains, not the command, set where it ends.

## 7. Figures

See [gallery.md](gallery.md):

- `search-14-rabbit-convergence`: the best feasible cost over trials.
- `search-15-rabbit-budget`: what became of the trials, and the final costs.
- `search-16-rabbit-gains`: the gains settled on.
- `search-17-rabbit-front`: speed against power, and the hypervolumes.
- `search-18-rabbit-strobes`: four gaits of the front.
- `search-19-rabbit-before-after`: the default and the tuned gains at 0.6 m/s.

## 8. Caveats

- **The loop is the experiment's, not the domain's.** `robotics::legged::Simulation` (and likewise `swimming::Simulation` and `quadruped::Simulation`) cannot be given controller gains. The loop is checked bit for bit against it, with the default gains and the feet's readings, but it is a copy of the `Simulation`'s semi-implicit path, and it will not follow changes to that path. A `gains` field in `SimulationConfig` would remove it.
- **The feet's readings are left out.** The belief they feed could bind only if three quarters of it fell below the 0.5 of friction asked, below 0.667. On this uniform ground of μ = 0.7 the default motion is unchanged without them, but a tuned gait that slides more could, in the domain's `Simulation`, feel a little differently.
- **Semi-implicit Euler is of order 1.** Its energy residual is the largest of the three methods (theory §5.2). Costs of transport here carry that error; the variational integrator would cost about a hundred times as much a trial.
- **The ranges bound the answer.** The best gains put the leg's rest length at its upper bound, 0.215 m. A wider range might go further, if the inverse kinematics and the launch's refusal allow it.
- **Failures carry no degree.** A trial that falls at 0.3 s and one that falls at 8 s are the same failure. That cost CMA-ES one seed. A graded measure of the failure as a constraint (the time survived), in place of a failure, would give the searches a direction in the region of falls; the choice is the experimenter's, and here the failures are as the domain reports them.
- **Budgets are small.** 800 trials for ten gains is about 80 per dimension, and 8 seeds. The ranking of DE, PSO and CMA-ES at the end is within their spread, and three seeds do not rank NSGA-II against SMS-EMOA.
- **One body, one ground, the Moon only.** On the Earth the domain refuses the hop (the ankle's motor is too weak for the launch), so no Earth task was posed. The planar model has no roll or yaw, no spine and no motor dynamics; braking energy is not recovered.
- **A window of 4 s** measures a gait after 4.5 s of hopping. A gait still changing then is measured as it was.

## 9. References

- Beume, N., Naujoks, B. and Emmerich, M. (2007), SMS-EMOA: multiobjective selection based on dominated hypervolume, *Eur. J. Oper. Res.* 181, 1653–1669.
- Deb, K. (2000), An efficient constraint handling method for genetic algorithms, *Comput. Methods Appl. Mech. Engrg.* 186, 311–338.
- Deb, K., Pratap, A., Agarwal, S. and Meyarivan, T. (2002), A fast and elitist multiobjective genetic algorithm: NSGA-II, *IEEE Trans. Evol. Comput.* 6, 182–197.
- Hansen, N. (2016), The CMA evolution strategy: a tutorial, arXiv:1604.00772.
- Kennedy, J. and Eberhart, R. (1995), Particle swarm optimization, *Proc. IEEE ICNN* 4, 1942–1948.
- Raibert, M. H. (1986), *Legged Robots That Balance*, MIT Press.
- Rechenberg, I. (1973), *Evolutionsstrategie*, Frommann-Holzboog.
- Storn, R. and Price, K. (1997), Differential evolution, *J. Global Optim.* 11, 341–359.
- Tucker, V. A. (1975), The energetic cost of moving about, *American Scientist* 63, 413–419 (the cost of transport).
