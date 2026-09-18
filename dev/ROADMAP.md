# Roadmap

This file gives agents and contributors context on which direction proposals
should lean. It is a living document — update it as priorities shift, and
record any resulting structural decision as an ADR (`dev/adr/`).

## Current focus

- Keep pace with the Bioconductor release cadence (see `dev/RELEASE.md`).
- Maintain the Bayesian clustering (`celda_C`/`celda_G`/`celda_CG`) and
  DecontX contamination-removal models as the two core feature areas.

## Not currently planned

- A Shiny app front-end (none exists today; would need an ADR to add one).
- Migrating away from the Rcpp/RcppEigen-backed Gibbs sampling implementation.

## Open questions worth an ADR if pursued

- Migrating the pkgdown site off a committed `docs/` directory to a
  `gh-pages` branch (see `dev/RELEASE.md`).
