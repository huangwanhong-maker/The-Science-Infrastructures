# development: Decisions, Ideas, Release Notes, the Roadmap and the Tools

*Status: draft.* The process of the infrastructure, kept apart from its content.

| Path | Purpose |
|---|---|
| [`decisions/`](decisions/) | Records of the decisions that shape the infrastructure, with their reasons. Each foundation library keeps its own. |
| [`ideas/`](ideas/) | Seeds: ideas written down early, and promoted or abandoned later. An idea about one foundation goes in that library's own. |
| [`release-notes/`](release-notes/) | What each commit brought: code and features developed, faults fixed; one file per day, written at each commit (ADR 0037). Each foundation and each domain in a repository of its own keeps its own. |
| [`roadmap.md`](roadmap.md) | The phases of the infrastructure: ordered intentions, not commitments. |
| [`tools/`](tools/) | The checks' harness (`checks.mjs`, which also builds and runs C++ checks with CMake); the checker of the dependency rules (`check-dependencies.mjs`); each domain run alone with its down-set (`check-alone.mjs`); one library's C++ as a suite (`run-cpp.mjs`); and the runner of every check (`run-tests.mjs`) |
