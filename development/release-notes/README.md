# release-notes: What Each Commit Brought

## 1. Purpose

A running record, for readers, of what each commit made possible or better, in the infrastructure's own repository; the foundations and the domains in repositories of their own keep theirs ([ADR 0037](../decisions/0037-release-notes.md)). It is written at each commit, and not between commits.

## 2. Convention

- One file per day, `YYYY-MM-DD.md`; within it, one short section per commit (or per group of commits made together), newest first, headed by what it brought.
- Each section is a clear, concise summary in plain prose: what is new, what one can now do, what works better ("In this update the atlas can read EEG in EDF, EDF+, BDF and BrainVision …"). It is not a list of changed files or functions; the git log keeps that. The decision or note behind it is linked.

## 3. Index

| Date | Commits |
|---|---|
| [2026-10-24](2026-10-24.md) | The history as a field theory; the history domain and the History Complex |
| [2026-10-23](2026-10-23.md) | The generation of spacetime; spin foams, vertex by vertex; complexes, symmetry and fields; the Complex Atlas; the simulator in the Quick Start; spiking networks in the gallery; the simulator's visualizer and profiler; the Spiking Simulator; a quantum neuron; seven new neuron models; the mathematical foundations and the first C++; spiking networks; biosignals on the brain; the Brain Atlas rebuilt; the API; observations, scores and hazards |
