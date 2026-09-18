# 1. Record architecture decisions

Status: Accepted

## Context

celda is a long-lived Bioconductor package maintained across multiple
contributors and AI coding agents. Decisions about dependencies, module
structure, and class design get re-litigated when the reasoning behind them
isn't written down anywhere.

## Decision

We will use Architecture Decision Records, as described by Michael Nygard,
for decisions that are hard to reverse or that future contributors will
question: dependency changes, file/module splits, S4 class redesign, and
build/deploy machinery changes. Routine choices (bug fixes, new arguments,
small helpers) do not need one.

Proposals go through a GitHub issue first; once approved, the ADR is added
here, numbered sequentially, and linked from `dev/adr/README.md`'s index.
ADRs are append-only — a later decision that reverses an earlier one is a new
ADR marking the old one "Superseded."

## Consequences

Contributors and agents have one place to check "why is it built this way"
before proposing a structural change, instead of guessing from git blame.
