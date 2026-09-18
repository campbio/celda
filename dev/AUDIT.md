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
