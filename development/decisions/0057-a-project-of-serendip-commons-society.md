# 0057. A Project of Serendip Commons Society: Stewardship without Ownership, and the Record of Contributors

- **Status:** Accepted
- **Date:** 2026-10-30
- **Decided by:** the user, Wanhong Huang, who is also the Founder and sole Director of Serendip Commons Society: "We want make the Science Infrastructure a project of it, but license still open source, and SCS not get the IP. Contributors.txt need serious made", and "consider a reasonable, and anyone free to continue develop and use this infrastructure ways". The form of the instruments and of the contributor terms is Claude's, after the Society's Bylaws and policies.

## 1. Context

Serendip Commons Society is a Delaware nonstock nonprofit corporation, founded to develop infrastructures for relational being. It brings a project under its governance by a Project Approval (the `SCS/PA` series) and publishes what it decides with an OpenTimestamps proof (its guidance `SCS/PGD/2026/01`). Four constraints bear on how The Science Infrastructures can join it.

- **Its Intellectual Property Policy** (`SCS/POL/2026/01`) claims no ownership of what a project produces: intellectual property belongs to its creators unless the project's terms provide otherwise.
- **Its Initial Period** (Bylaws, Article XIV), while the Founder is the only Director, forbids accepting donations or other transfers of value from natural persons other than the Founder (§14.10(a)) and contracts longer than twelve months (§14.11(a)). An assignment of copyright, or a contributor licence agreement, from a contributor to the Society would be one or the other.
- **The infrastructure was the user's own work before it joined.** It began on 2026-10-21 (ADR 0001); every repository is MIT-licensed, with the copyright line "Wanhong Huang" (CLAUDE.md §7).
- **Most commits were written with AI assistants** and say so in a `Co-Authored-By:` trailer. None carries a sign-off, and no file recorded who the contributors are or on what terms they contribute.

## 2. Decision

### 2.1 The Project and Its Instruments

The infrastructure is a project of the Society, approved as `SCS/PA/2026/05` and announced by the notice `SCS/NTC/2026/02`, both of 30 October 2026, both published with their OpenTimestamps proofs. The approval classifies it under the Society's Schedule A (Rev.1), Domain I, Theoretical and Ontological, with Domains VI, Educational, and II, Epistemic, as secondary; finds that its present scope needs no ethics examination under the Bylaws' §3.5, with stated conditions; and keeps the project's technical decisions with its maintainers, in this log. The Society presents it at `https://serendip.ngo/science-infrastructures/`.

The approval covers the seven repositories: this one, calx, morphe, physis, tyche, neuroscience and Cosmopoiesis-Lab. A foundation or domain that later becomes a repository of its own joins by having a `CONTRIBUTORS.txt` on these terms.

### 2.2 What the Society Holds

Nothing. The Society holds no copyright, patent, database right or other intellectual property in the infrastructure, takes no licence beyond the MIT licence everyone receives, and cannot relicense the work. It is a steward: it keeps the project's record and presents it. A contribution whose author cannot be identified stays under the MIT licence; the Society records it and holds it only as custodian (its IP Policy, §6).

The public instances of the three atlases on the user's server stay the user's own, not an activity of the Society. Bringing them under it would first need the documentation of the jurisdiction the Bylaws require (§14.12(e)) and a privacy notice, since a server keeps visitors' addresses in its logs.

### 2.3 The Licence and the Copyright

The licence is unchanged: MIT, in every repository. The copyright line becomes "Wanhong Huang and the contributors listed in CONTRIBUTORS.txt", which is true now (there is one contributor) and stays true as others join, without anyone assigning anything.

### 2.4 The Record of Contributors

Every repository carries `CONTRIBUTORS.txt`, with the same terms in each:

1. **Who holds what:** each contributor keeps their copyright; the Society holds nothing; the MIT licence cannot be withdrawn from a copy received, so anyone may fork and continue the work with or without the project or the Society.
2. **How a contribution is made:** inbound = outbound, under MIT to everyone; no contributor licence agreement and no agreement with the Society; each commit signed off (`git commit -s`), by which its author makes the Developer Certificate of Origin 1.1, reproduced in the file. The commits made before the file existed are certified by their one author, through the file.
3. **AI assistance:** a `Co-Authored-By:` trailer records the tool used and names no contributor. An AI system holds no rights and is not listed; the person who submits AI-assisted material is its contributor and answers for it.
4. **The steward:** what the Society does and does not do, how the stewardship ends (nothing is transferred, because nothing is held), and the Society's name kept as its own.
5. **Material of others:** the files that carry third-party material under its own licence, with their notices; the facts in the fixtures; the datasets, never committed; and, in this repository, the three screenshots of Pangu-Weather's forecasts, which are CC BY-NC-SA 4.0 and not MIT.
6. **The register:** each contributor, added in the change of their first contribution.
7. **The certificate,** verbatim.
8. **Changes:** the terms change only in every repository together, by a decision, and never for contributions already made.

## 3. Alternatives Considered

- **The Society holds the copyright,** as in its approvals of G³ (`SCS/PA/2026/03`) and of Fluid Dynamics in Sociological Systems (`SCS/PA/2026/04`). Rejected: the user asked that the Society not take the intellectual property; the Society's later IP Policy makes non-retention its default; and contributors other than the Founder could not assign to it during the Initial Period (§14.10(a)).
- **A contributor licence agreement with the Society.** Rejected: it is a contract with every contributor, of indefinite duration (§14.11(a)), and it would give the Society a power to relicense that no one asked for. Inbound = outbound under MIT with the certificate gives every party what it needs.
- **A copyleft or a non-commercial licence.** Rejected: MIT is this infrastructure's decision (CLAUDE.md §7), and the user wants anyone free to continue and use it.
- **Joining without a record of contributors.** Rejected: the Society's IP Policy settles affiliation from the project's record (§5), which must therefore exist and say who made what.

## 4. Consequences

- From now on every commit carries the human author's `Signed-off-by:` line (`git commit -s`), beside any `Co-Authored-By:` trailer.
- A new contributor is added to the `CONTRIBUTORS.txt` of the repository they contribute to, in the same change.
- The Society's record of the project (its approval, the notice and its reports) lives in the Society's repository, not here; the release notes here remain the project's own account of its progress.
- Nothing in the code, the checks or the licence of any dataset changes.
