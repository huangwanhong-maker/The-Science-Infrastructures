# 0001. Record Architecture Decisions

- **Status:** Accepted; location of the log amended by [0009](0009-development-folder.md)
- **Date:** 2026-10-21
- **Decided by:** the user and Claude, as part of the revisability principle

## 1. Context

The project is long-lived, and its direction will change as ideas emerge. Decisions made early will be questioned later, both by the user and by Claude in sessions that do not remember the original discussion. A decision whose reasons have been lost can only be kept out of habit, or changed blindly.

## 2. Decision

- Record each decision that shapes structure, conventions or tooling as a numbered record in `docs/decisions/`, using the template in the [index](README.md).
- Never rewrite a record. A changed decision gets a new record that supersedes the old one.

## 3. Alternatives Considered

- **Commit messages only.** Rejected: the reasons become hard to find, and the folder is not yet under version control.
- **A single design document, edited in place.** Rejected: it loses the history of why things changed.

## 4. Consequences

- Each decision costs a small amount of writing.
- CLAUDE.md and `docs/architecture.md` describe the current state. The records explain how that state came about.
