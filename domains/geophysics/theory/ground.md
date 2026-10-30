# The Earth's Ground as a Foot Meets It

*Status: draft (2026-10-31). Serves use and research.* The theory of `ground/presumed.hpp` ([ADR 0059](../../../development/decisions/0059-robots-in-their-environments.md)).

## 1. What the World Knows of the Ground at a Place

Nothing has measured how a robot's foot meets the ground at most places of the Earth. What the world does hold everywhere is CRUST1.0 (Laske et al. 2013), which describes each cell of 1° by its layers: water, ice, three of sediments, three of crust. From it the ground's class is presumed:
- **sea:** where the cell has water; not ground;
- **ice:** where it has ice (Greenland, Antarctica, the ice caps);
- **sediment:** where it has sediments (soil, sand, gravel);
- **rock:** elsewhere.

A cell of a degree is the whole of what is known: a rock outcrop in a basin, or snow on rock, is not seen. The checks see ice at Summit on the Greenland ice sheet and on Antarctica's plateau, sediment in the Paris Basin, and the sea in the equatorial Pacific.

## 2. Each Class's Law

physis's compliant contact, with its three parameters each log-normal, a median and the spread of its logarithm:

| Class | Friction μ | Load at 1 mm | Damping |
|---|---|---|---|
| ice | 0.15, spread 0.45 | 30 N, 0.4 | 1.5 s/m, 0.3 |
| sediment | 0.60, spread 0.25 | 15 N, 0.6 | 2.0 s/m, 0.3 |
| rock | 0.80, spread 0.15 | 40 N, 0.4 | 1.5 s/m, 0.3 |

All are assumptions, labelled so in the fields' sources:
- **Friction.** Rubber on ice slides at a few hundredths near the melting point and at several tenths in deep cold, so ice's spread is wide; rubber grips dry soil and rock at the order of 0.5 to 1.
- **Stiffness and damping.** On ice and rock they are those of the foot's rubber pad, in series with a ground far stiffer than it. On sediment they are the softer ground's own.

A measured law, with the temperature for ice, is the next step.

## 3. Its Realisation

Each parameter varies as its class's median times $e^{s\,z(x)}$, with $z$ a stationary Gaussian field of unit variance and a correlation length of 1 m (physis's spectral representation, its decision 0015). Each is read at the point's place in the Earth's fixed frame, on its own stream of the world's seed. The presumption is a field of provenance model; the realisation, synthetic, is the ground a simulated world's robots stand on.

## References

- Laske, G., Masters, G., Ma, Z. and Pasyanos, M. (2013), Update on CRUST1.0: a 1-degree global model of Earth's crust, *Geophys. Res. Abstracts* 15, EGU2013-2658.
