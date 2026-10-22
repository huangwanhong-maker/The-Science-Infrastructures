# 0013. Galaxies and Their Relations in the Celestial Atlas

- **Status:** Accepted
- **Date:** 2026-10-22
- **Decided by:** the user asked how to investigate relations between galaxies, and asked to visualise them. Claude chose the catalogues, the methods and the views.

## 1. Context

Until now the atlas reached only the stellar neighbourhood. Galaxies are related in several ways:

- by gravity: tidal influence, group and cluster membership;
- by motion: the Hubble flow and departures from it;
- by structure: groups, filaments and voids.

Investigating these relations means computing them, not only drawing published values.

## 2. Decision

1. **Two catalogues**, both queried from CDS VizieR through the data registry ([ADR 0012](0012-selectable-parameterised-datasets.md)):
   - **The Updated Nearby Galaxy Catalog** (Karachentsev, Makarov and Kaisina 2013). It covers 869 galaxies within about 11 Mpc, with measured distances, luminosities, tidal indices and main disturbers.
   - **The 2MASS groups of Tully (2015).** It covers 43 037 galaxies of the 2MASS Redshift Survey, with their group membership. Distances follow from velocities, in h⁻¹ Mpc.
2. **Recomputed relations.** The tidal index and the main disturber are recomputed from positions and luminosities in `src/core/galaxies.js`, and tested against the published values:
   - the main disturbers agree for 99.2 % of galaxies;
   - Θ₁ agrees to a median of 0.027 dex.

   The Milky Way's luminosity, which the catalogue omits, is inferred from the published indices of its satellites (log L = 10.55).
3. **Tidal groups.** Groups are defined by joining each galaxy with Θ₁ > 0 to its main disturber.
4. **Names from SIMBAD.** Group names come from SIMBAD at build time:
   - the name of each group's brightest member;
   - the names of catalogued clusters lying within a group's turnaround radius whose velocity agrees with the group's to within 1500 km/s.
5. **The Galaxies view.** It has three parts:
   - **a 3D map** in supergalactic coordinates, with relations drawn as lines. Each galaxy links to its main disturber, or to its group centre;
   - **a Hubble-flow diagram**, distance against Local Group velocity, for the Local Volume;
   - **a mass–light diagram**, virial mass against K-band luminosity, for groups.
6. **Distances labelled by kind.** Kinematic distances are labelled as such, with units of h⁻¹ Mpc. They are never presented as measured.

## 3. Alternatives Considered

- **Show only published relations.** Rejected: that would not be investigation, and it would leave the method untested.
- **Cosmicflows-4 first.** Deferred. Its peculiar velocities support the study of flows, which is the natural next step (see the [idea](../ideas/2026-10-22-cosmic-flows.md)).
- **Hand-written names for famous clusters.** Rejected in favour of names resolved from SIMBAD and checked by position and velocity.

## 4. Consequences

- The supergalactic transformation is tested against Tully's own coordinates: the median difference is 0.00012° over 43 036 galaxies. One catalogue anomaly is recorded in the tests.
- The tests depend on the data cache, as the DE440 tests already did.
- The building blocks are reusable in physis and in `core/domains/astronomy`: the tidal index, the inference of missing masses, and groups from relations.
