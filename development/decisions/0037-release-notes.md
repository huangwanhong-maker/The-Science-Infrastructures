# 0037. Release Notes at Each Commit

- **Status:** Accepted
- **Date:** 2026-10-23
- **Decided by:** the user asked for a release-notes folder recording what was developed, improved or fixed, updated at each commit and not continuously, written as a concise summary in natural language rather than as commit messages. Claude chose its place and form and wrote this record.

## 1. Context

The decision log records why the structure is as it is, the ideas what may come, and the git history every change; none says, for a reader, what each commit brought.

## 2. Decision

1. **`development/release-notes/`** in every repository of the project (the infrastructure, its foundations and its domains in repositories of their own), beside decisions and ideas.
2. **One file per day,** `YYYY-MM-DD.md`, with one short section per commit, newest first: a clear, concise summary in plain prose of what is new, what one can now do and what works better, linking the decision or note behind it; not a list of changed files, which the git log keeps (the user's correction of a first draft written as such lists). The folder's README indexes the days.
3. **Written at each commit,** in each repository committed (submodules first, as their commits come first), and not between commits.
4. A shared convention, with the writing style, the status vocabulary, decisions and ideas: changed in all repositories together.

## 3. Alternatives Considered

- **One CHANGELOG file at the root.** It grows without bound in one file and mixes repositories; a folder by day reads as the decision log does.
- **Notes generated from commit messages, or written like them.** Commit messages are written for the history, not for readers of what was built; the notes say what a change means for its reader.

## 4. Consequences

- Every CLAUDE.md states the convention; the first notes are of 23 October 2026.
