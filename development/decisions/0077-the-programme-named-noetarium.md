# 0077. The Programme Named Noetarium: the Science and Humanity Infrastructure

- **Status:** Accepted
- **Date:** 2026-11-06
- **Decided by:** the user, who renamed the repository on GitHub to `Noetarium-The-Science-and-Humanity-Infrastructure` and asked: "Remote rename to git@github.com:huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure.git. And also revise the README title and related name." Claude carried the name through the repository.

## 1. Context

The programme was named The Science Infrastructures when it was founded (ADR 0025). Since then:
- it has gained the humanities (the history domain, ADR 0041) and the social sciences (economics, nomos);
- its foundations have grown to seven.

Its repository is now `huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure` on GitHub, which redirects the old address.

## 2. Decision

1. **The name.** The programme is *Noetarium*, in full *Noetarium: The Science and Humanity Infrastructure*. "The infrastructure" stays its short name in prose.
2. **The remote.** The root repository fetches from `https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure.git` and pushes over the SSH alias `serendip_wanhong` to the same repository, as every repository of the programme does (CLAUDE.md §7).
3. **What is renamed in the root repository:**
   - the README's title and opening line, which now names the humanities beside the sciences;
   - CLAUDE.md;
   - CONTRIBUTORS.txt;
   - the API's name and its OpenAPI title;
   - the CMake messages of the domains;
   - the data library's User-Agent (`noetarium-data/1`, with the new address);
   - the applications' READMEs and pages, and the documents.
4. **What is left as written:**
   - the decision records, which are never rewritten;
   - past release notes;
   - the experiment records, whose cited versions belong to the repository as it was then named;
   - the Serendip Commons Society's page for the project (`serendip.ngo/science-infrastructures/`), which is the Society's to change;
   - Project Approval SCS/PA/2026/05, which approved the project under its former name. CONTRIBUTORS.txt says so.
5. **The repositories attached as submodules** (the seven foundations, the neuroscience domain, Cosmopoiesis) still name and link the programme by its former name; GitHub's redirect keeps their links working. Renaming them needs a commit in each repository, and waits for the user's word.

## 3. Alternatives Considered

- **Renaming every mention, records included.** Rejected: records are never rewritten, and their names belong to their dates.
- **Keeping the old name in prose and changing only the remote.** Rejected: the user asked for the README's title and the related name to change.

## 4. Consequences

- The Python package keeps its import name, `science`, so that no caller breaks; only its docstrings name the programme.
- The development machine's working directory (`The Science`) and the public server's clone are unaffected.
