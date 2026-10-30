# A Dog: Four Legs on the Earth's Ground

*Status: draft (2026-10-31). Serves use, learning and research.* The theory of the module `quadruped` ([ADR 0059](../../../development/decisions/0059-robots-in-their-environments.md); `dog.hpp`). Its mechanics is physis's planar trees and compliant contact (decision 0014). Its ground and its belief of it are the rabbit's (`legged/ground.hpp`, [theory: legged §7](legged.md)), on the Earth's ground presumed by the geophysics domain (`ground`).

## 1. The Body as Built

A torso of 7.0 kg, 0.40 m long, and four legs, each a thigh and a shank of 0.20 m with a point foot. The fore legs hang at the shoulders, 0.18 m ahead of the torso's centre, and the hind legs at the hips, 0.18 m behind. The dog weighs 10.2 kg. Its knees bend backward, as a robot dog's do.

Its configuration space is $SE(2)\times T^8$ in its sagittal plane. The left and right legs of a pair are separate links in that one plane, drawn at ±0.10 m across it. This is what lets it trot: its diagonal pairs (right fore with left hind, left fore with right hind) stand in turn. At every instant a fore foot and a hind foot are down, so the support spans the body in the plane, as it does in three dimensions.

The motors are geared (the domain's `actuation`): 20 N m at the hips and 30 N m at the knees, with their torque–speed lines, armature, friction and copper loss.

## 2. The Feet on the Ground

The contact is physis's compliant contact. On ice and rock its stiffness and damping are those of the foot's rubber pad, in series with a ground far stiffer than it: a load of 30–40 N presses the foot 1 mm in, with a damping of 1.5 s/m. On sediment they are the softer ground's own.

**Remark (an integrable friction).** Below the slip speed $v_s$, Coulomb's regularised friction acts on a foot like a damper of $\mu f_n/v_s$. Under a dog's foot, 30 N on 0.2 kg, at $v_s = 0.01$ m/s and a step of 0.5 ms, $h\mu f_n/(v_sm) = 5$. An explicit step is unstable above 2, and the feet chattered. The dog's feet slip at $v_s \ge 0.03$ m/s and its step is 0.25 ms, which makes that number 0.8. The rabbit's lighter feet did not need it.

## 3. The Controller

At 1 kHz.
- **Standing.** Each joint follows its posture by a PD law, and the weight is carried by $-J^{\mathsf T}f$ for each foot's share, the fore and the hind pairs by the lever rule.
- **Trotting.** A gait clock of period 0.36 s. The pair in stance gives the body the force and the moment asked:
  - vertically, the weight and a servo of the torso's height to 0.30 m;
  - forward, $m\,(k_v(v^* - v) + k_i\!\int(v^* - v))$, the integral removing the steady error the placement leaves;
  - the moment that holds the pitch level.

  The force is distributed between the stance pair's fore foot and hind foot. Each takes half the forward force, and their vertical forces solve the balance of the moment. No foot is asked for more friction than three quarters of what the dog believes is there, less one deviation. The swing pair follows an arc of 6 cm's lift to Raibert's placement ahead of its mount, $x_f = vT/4 + k_v(v - v^*)$, by each leg's inverse kinematics.
- **Stopping.** The speed asked is zero, and when the dog is slow it stands.

## 4. What the Checks Hold

`check_quadruped` (18):
- **The body:** its mass, and the legs' inverse kinematics.
- **Standing on rock:** still ($6\times10^{-5}$ m in a second), level, its weight carried by its four feet ($2\times10^{-6}$), its energy accounted for (22 mJ).
- **Trotting:** in place within 12 cm of where it stood. Forward at the 0.3 m/s asked, its mean over the last 4 s within $10^{-4}$ m/s, level within 1°, its torso within 1.1 cm of its height.
- **Stopping:** from 0.3 m/s it stops and stands, still and level.
- **The three methods** agree on the distance trotted in 4 s: semi-implicit Euler within $4\times10^{-4}$ of Runge–Kutta, the variational integrator within $10^{-4}$. Their costs are 0.03, 0.10 and 6.3 s a second. Their energy residuals, over 4.5 s with 25 impacts of feet, are 23, 4.9 and 9.5 J: the contacts' stiffness at a 0.25 ms step is where the methods' errors of energy lie.
- **On ice** (μ 0.15, told so), it gathers speed more gently than on rock: 0.296 m against 0.303 m in its first 1.5 s. A steady trot asks little of the ground; gathering speed asks more, within its belief. It trots on ice level.

Through Cosmopoiesis, on the Greenland ice sheet, it was told μ 0.15. The ice under its feet was 0.317, and after 400 readings of sliding feet it believed 0.296 ± 17 %; its readings' whitened mean square was 1.09.

## 5. Open Questions

- **Three dimensions**, where a trot is a trot of four legs and a turn is possible: the spatial tree.
- **Slipping on ice.** Its controller asks less of icy ground, but a gait that keeps its feet from sliding (shorter strides, a slower clock) is not chosen by it.
- **Its energy.** The cost of transport, and the copper loss of holding a stance.

## References

- Raibert, M. H. (1986), *Legged Robots That Balance*, MIT Press.
- Di Carlo, J., Wensing, P. M., Katz, B., Bledt, G. and Kim, S. (2018), Dynamic locomotion in the MIT Cheetah 3 through convex model-predictive control, *IEEE/RSJ International Conference on Intelligent Robots and Systems* (the force distribution of a trot, which this planar dog does by its moment balance alone).
