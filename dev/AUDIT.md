# Dependency / deprecation audit

Run periodically (not on a fixed schedule yet — trigger manually or from the
weekly BioC-check cron once it's promoted to required). Prompt to run:

> Run `BiocCheck::BiocCheck()` and `BiocManager::valid()` on this package.
> Use r-lib lifecycle practices to find deprecated functions or S4 methods
> and report replacements compatible with the current Bioconductor release.
> Write findings to `dev/agent-log.md`. Do not change code — findings become
> GitHub issues (and ADRs where structural).

## Log

No formal audits recorded yet. Findings get appended to `dev/agent-log.md`
(created on first run) and summarized here with a date and links to the
resulting issues.

- **2026-09-18** (found during AI-agent tooling setup, not a formal audit):
  `pkgdown::check_pkgdown()` fails with `In _pkgdown.yml, url is missing.` —
  pre-existing on `devel`/`master`, unrelated to the setup PR. The new
  `pkgdown-check` CI job is non-blocking until this is fixed (add a `url:`
  field to `_pkgdown.yml`) and the job is flipped to required.
- **2026-09-18**: `make lint` (new `.lintr` config) surfaces pre-existing
  style issues in `vignettes/` (line length, indentation, object naming) —
  not touched by the setup PR. The `lintr` CI job should stay non-required
  until these are triaged; consider filing a GitHub issue to track cleanup.
- **2026-09-19** (formal audit — first run of this file's own prompt; run
  inside the official `bioconductor/bioconductor_docker:RELEASE_3_23`
  container, R 4.6.1, current Bioconductor release, celda's full dependency
  tree rebuilt from source, BiocCheck 1.48.1): `BiocCheck::BiocCheck()` on
  the built tarball came back 0 ERRORS | 1 WARNING | 18 NOTES;
  `BiocCheckGitClone()` came back 1 real NOTE (no `CITATION` file — its "2
  ERRORS" were local container-copy artifacts, not real, see log). `make
  test` passed with 0 failures but **21 warnings** — 20 from
  `scuttle::librarySizeFactors()`/`normalizeCounts()` being deprecated
  upstream (hit via `scater::logNormCounts()` in every DecontX call,
  `R/decon.R`), 1 from the lifecycle sweep's big finding:
  `ggplot2::aes_string()` (deprecated since ggplot2 3.0.0), used **36
  times** across `R/perplexity.R` (22), `R/plot_dr.R` (11), and
  `R/plot_decontx.R` (3). `BiocManager::valid()` on a clean release image:
  8 out-of-date, 0 too new, no archived dependencies among celda's own
  Imports/Depends/Suggests. Full detail and suggested next steps in
  `dev/agent-log.md`. Nothing fixed here (scope guard) — recommend filing
  GitHub issues for the `aes_string()` migration, the `scater`/`scuttle`
  deprecation, and the quick `DESCRIPTION` fixes (URL field, R version
  bump to 4.6.0, LazyData setting, fnd role, vignette chunk labels).
