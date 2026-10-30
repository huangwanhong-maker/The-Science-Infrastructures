# The Ground as a Foot Meets It: the Moon's Regolith Presumed

*Status: draft (2026-10-30). Serves use and research.* The theory of `surfaces/regolith.hpp` ([ADR 0059](../../../development/decisions/0059-robots-in-their-environments.md)).

## 1. The Law of a Contact

A foot pressed into the ground by $\delta$ and sliding at $v_t$ feels physis's compliant contact (its theory §7):
$$f_n = k\,\delta^n(1 + c\,\dot\delta), \qquad f_t = -\mu f_n\tanh(v_t/v_s).$$
Three parameters name the ground: the friction $\mu$, the stiffness $k$ (with Hertz's exponent $n = 3/2$) and the damping $c$. The stiffness is given as the load that presses a foot 1 mm in, $k\,(1\text{ mm})^n$, in newtons. $k$ itself has the unit N/m$^{3/2}$, whose half-integer exponent of length physis's units do not carry.

## 2. The Presumption

No instrument has measured these parameters for a small foot at most places of the Moon. Each is presumed as a log-normal law, with a median and the standard deviation of its logarithm:

| Parameter | Median | Spread of the logarithm | Basis |
|---|---|---|---|
| $\mu$ | $\sqrt{\tan 30°\,\tan 50°} = 0.830$ | $\ln(\tan 50°/\tan 30°)/3.92 = 0.185$ | $\tan\varphi$, with $\varphi$ the regolith's angle of internal friction, 30° to 50° (Carrier, Olhoeft and Mendell 1991): the law's 95 % interval |
| $k\,(1\text{ mm})^n$ | 9.5 N | 0.5 | assumed: the rabbit's contact as built |
| $c$ | 1.5 s/m | 0.3 | assumed: a restitution near 0.9 at 0.1 m/s |

*Remark.* The friction of a foot on a granular ground is bounded by the ground's own shear strength. $\tan\varphi$ is therefore the friction of a rough foot, the case in which the regolith shears before the interface slips. A smooth foot slips at less. The assumed parameters are labelled as assumptions in their fields' sources, so that a reading shows what it rests on.

## 3. Its Variation, and a Realisation

Each parameter varies over the ground as $p(x) = p_{50}\,e^{s\,z(x)}$, where $s$ is its spread and $z$ a stationary Gaussian field of unit variance with correlation length $\ell = 1$ m (assumed: a robot's step). $z$ is read at the point's place in the Moon's fixed frame by physis's spectral representation (its decision 0015), independently for the three parameters. A world's seed fixes $z$: that realisation is the ground of the simulated world (provenance synthetic), and the presumption is its law (provenance model).

**Proposition 3.1.** The realised parameter is log-normal at each point, with the presumed median and spread. Its logarithms at two points correlate as $e^{-|x - y|^2/2\ell^2}$.

*Proof.* $\ln p(x) = \ln p_{50} + s\,z(x)$ is affine in $z$, and $z$ has unit variance and the squared-exponential covariance. $\square$

The checks see over 2000 places a kilometre apart each parameter's $z$ of mean 0 and deviation 1, within four standard errors. Two feet 5 cm apart stand on ground of correlation 0.999, and two places 4 m apart on independent ground.

## 4. What Is Not Yet Here

- **Plastic sinkage.** Regolith compacts under a foot and keeps the print. A plastic law, such as Bekker's pressure–sinkage relation or a critical-state model, would let the ground remember a robot's steps: a field the robot and the ground make together.
- **Measured mechanics.** Apollo's penetrometer and the depth of its boot prints would constrain the stiffness in place of an assumption.
- **The ground's variety.** Rocks, crater rims and the ice of the permanently shadowed regions would each be a class with its own law.

## References

- Carrier, W. D., Olhoeft, G. R. and Mendell, W. (1991), Physical properties of the lunar surface, in Heiken, G. H., Vaniman, D. T. and French, B. M. (eds.), *Lunar Sourcebook*, Cambridge University Press, ch. 9.
- Shinozuka, M. and Deodatis, G. (1991), Simulation of stochastic processes by spectral representation, *Appl. Mech. Rev.* 44, 191–204.
