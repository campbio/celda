# Contributing to celda

Thanks for your interest in contributing! This applies equally to human
contributors and AI coding agents working in this repo.

## Ground rules

- celda follows the shared
  [r-bioc-dev-standards](https://github.com/campbio/r-bioc-dev-standards);
  [`AGENTS.md`](../AGENTS.md) adds what's specific to celda.
- Branch from `devel` (`fix/<topic>` or `feature/<topic>`). Never push
  directly to `devel`, a `RELEASE_*` branch, or `master`.
- `master` is an automatic copy of the current Bioconductor release
  (kept in sync by `.github/workflows/sync-stable.yaml`). Don't commit to
  it or open PRs against it; a check fails any PR aimed at it.
- All changes land via pull request against `devel`, reviewed before merge.
  The one exception is the maintainer's Bioconductor sync on release day,
  which pushes `devel` and `RELEASE_X_Y` directly (see
  [`dev/RELEASE.md`](../dev/RELEASE.md)).

## Workflow

1. Fork the repo (external contributors) or branch from `devel` (lab
   members).
2. Make your change. While developing, run `make test-one FILTER=<pattern>`;
   before handing off, run `make test` and `make coverage`.
3. Update `NEWS.md` for any user-facing change.
4. Before opening a PR, run `make check-full` and `make bioccheck`.
5. Open a PR against `devel` using the pull request template.
6. Architectural changes (file splits, dependency changes, S4 redesign)
   require an approved ADR first; see [`dev/adr/`](../dev/adr/README.md).

The standard `make` targets come from the shared standards, so the first
`make` on a new machine needs internet access to download them. Run
`make help` to list them.

## Pull request checklist

See [`PULL_REQUEST_TEMPLATE.md`](PULL_REQUEST_TEMPLATE.md). Checks must
pass, `NEWS.md` must be updated, and a person must review changed results
for scientific correctness before merge.

## Code of conduct

This project follows the guidelines in [`CONDUCT.md`](../CONDUCT.md).
