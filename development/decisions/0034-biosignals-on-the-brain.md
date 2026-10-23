# 0034. Biosignals on the Brain: EEG, MEG and fMRI Carried onto the Template

- **Status:** Accepted
- **Date:** 2026-10-23
- **Amends:** [0028](0028-neuroscience-and-the-brain-atlas.md), [0032](0032-the-anatomy-of-the-brain-atlas.md) §2.4 (a sixth mode)
- **Decided by:** the user asked that biosignals, such as EEG, MEG and fMRI, be mapped onto the brain in three dimensions and visualised with high quality and creativity. Claude chose the design below and wrote this record.

## 1. Context

The Brain Atlas's Head mode (0032) placed the 10–20 and 10–10 electrodes on the template's scalp and carried a set of channel values onto the cortex by interpolation. That is a picture of the scalp, not of the brain. A potential on the scalp is made by currents in the cortex, smeared by the skull. Carrying it to the cortex needs a *forward model* (from currents to sensors) and an *inverse* (from sensors back to currents). fMRI is different again: its maps are already volumes in the template's space, but they must be read through the cortex onto its surface, and one person's unsmoothed map is mostly noise at single voxels.

Open data cover all three: PhysioNet's EEG Motor Movement/Imagery Dataset (Schalk et al. 2004; 64 channels of the 10-10 system, ODC-By 1.0), and the Individual Brain Charting project's z maps on NeuroVault (Pinho et al. 2018; CC0).

## 2. Decision

1. **Three modules in the neuroscience domain,** pure and tested (theory in `theory/forward-and-inverse.md`):
   - `signals.js`: the readers (EDF and EDF+ with their annotations, BDF, BrainVision, tables) and the analysis of recordings (FFT, zero-phase band-pass, the analytic signal, Welch's spectrum, band power in the IFCN bands, the average reference, global field power, covariance, intervals from annotations).
   - `forward.js`: the concentric-shell EEG model solved in Legendre series (three shells after Rush and Driscoll 1968, or any number), Sarvas's MEG field (1987), sphere fitting and the radial normalisation of a head to its sphere, lead fields, an idealised helmet, and source spaces decimated from the cortex.
   - `inverse.js`: the regularised minimum-norm operator with depth weighting (Lin et al. 2006), dSPM (Dale et al. 2000) and sLORETA (Pascual-Marqui 2002), with the power of an estimate from a sensor covariance.
   - `maps.js` and `volume.js` gain the reading of a volume through a surface, the peaks of a statistic, Gaussian smoothing confined to a mask, and robust restandardisation.
2. **A Signals mode** in the Brain Atlas, in three parts:
   - *EEG:* a recording on the template's head; its lead field computed in a worker; traces and an overview with the annotations; band power, the potential at an instant and the contrast of two conditions as estimates on the cortex; the scalp map and the spectrum on the card.
   - *EEG and MEG:* one current placed by a click, with its potential and its field, to show what each modality sees; the estimate from each.
   - *fMRI:* a map read through the cortex, smoothed (6 mm by default for the examples) and restandardised; peaks named by the anatomy; sections with the map over the template's image.
   Files of the same kinds can be opened or dropped on the brain.
3. **Two datasets join the data library:** `physionet-eegmmidb` (four runs of subject 1) and `neurovault-ibc` (six z maps of subject 1), with the Brain Atlas's recipe fetching them.
4. **The physics and the science are checked:** the shell model against its closed forms and numerical derivatives; Sarvas's identities; sLORETA's exact localisation without noise. On the real data: Berger's alpha rhythm (13 times stronger at O1 with the eyes closed, peaking at 10 Hz, estimated at the occipital poles); the mu rhythm subsiding over the central gyri as a fist moves; the right hand's fMRI peak in the left central gyri.

## 3. Alternatives Considered

- **A boundary-element model** of the template's head. It is more faithful, but it needs meshes of the skull's surfaces and a solver of a few thousand unknowns. The three-shell model after radial normalisation is the standard first approximation, exact for spheres; a BEM can replace it later behind the same lead-field interface.
- **Running MNE-Python in the server.** It is mature, but it would put the domain's science outside the domain, and the page would need the server for what is a few seconds of computation in a worker.
- **Showing only scalp maps.** They are honest but do not answer where in the brain the signal comes from. The estimate is shown as an estimate: the legend and the card say so.
- **The HCP group maps** for fMRI. They carry the HCP's data-use terms; the IBC maps are CC0, and one person's maps are more instructive once smoothed.
- **The IBC faces-against-shapes contrast.** For this subject it is negative in the fusiform gyrus (the shapes drive early vision harder), so the faces example uses faces against bodies, places and tools.

## 4. Consequences

- The neuroscience domain now models the physics of measurement, not only the anatomy; its checks rise from 120 to 162.
- The forward model is spherical; estimates near the temporal base and the frontal pole, where the head is least spherical, are the least reliable. The cards say that the template is not the subject.
- The recipe fetches about 36 MB more data.
- Spiking network models, which the user asked for next, belong in the domain as a module of their own, with an application apart from the Brain Atlas.
