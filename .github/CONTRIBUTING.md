# Contributing to celda

Thanks for your interest in contributing! This applies equally to human
contributors and AI coding agents working in this repo.

## Ground rules

- Read [`AGENTS.md`](../AGENTS.md) first — it's the canonical instruction file
  for this repo (conventions, sanctioned commands, safety rules).
- Branch from `devel`. Never push directly to `devel` or `master`.
- All changes land via pull request, reviewed before merge.

## Workflow

1. Fork the repo (external contributors) or branch from `devel` (lab members).
2. Make your change. Run `make test` after every change; run `make check`
   before opening a PR.
3. Update `NEWS.md` for any user-facing change.
4. Open a PR against `devel` using the pull request template.
5. Architectural changes (file splits, dependency changes, S4 redesign)
   require an approved ADR first — see [`dev/adr/`](../dev/adr/README.md).

## Pull request checklist

See [`PULL_REQUEST_TEMPLATE.md`](PULL_REQUEST_TEMPLATE.md) — `make test` and
`make check` must pass, `NEWS.md` updated, and results reviewed for scientific
correctness by a human before merge.

## Code of conduct

This project follows the guidelines in [`CONDUCT.md`](../CONDUCT.md).
