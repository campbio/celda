# Agent audit log

Findings from `dev/AUDIT.md`'s audit prompt. No code was changed to produce
these — this is a read-only investigation. Findings below should become
GitHub issues (and ADRs where structural) rather than being fixed silently.

## 2026-09-19 — BiocCheck / BiocManager::valid() / lifecycle sweep (Bioconductor 3.23 release, R 4.6.1)

Run inside the official `bioconductor/bioconductor_docker:RELEASE_3_23`
container (current Bioconductor release; devel is 3.24) so results match
what a real Bioconductor build machine would report, rather than whatever
Bioconductor version a contributor's local R happens to resolve to. celda's
full dependency tree (~309 packages) was rebuilt from source in that
environment first. BiocCheck version: 1.48.1.

### `BiocCheck::BiocCheck()` on the built tarball: 0 ERRORS | 1 WARNING | 18 NOTES

- **WARNING**: avoid `T`/`F` as logical literals — use `TRUE`/`FALSE` (16
  occurrences, mostly in `R/findMarkersTree.R`).
- **NOTE**: `LazyData: true` in `DESCRIPTION` should be `false` or removed.
- **NOTE**: bump `R (>= 4.0)` dependency to `R (>= 4.6.0)` to match current
  Bioconductor requirements.
- **NOTE**: add a `URL` field to `DESCRIPTION`.
- **NOTE**: add an `'fnd'` (funder) role to `Authors@R`, if applicable.
- **NOTE**: maintainer ORCID iD suggested in `Authors@R`.
- **NOTE**: vignettes with missing chunk labels — `celda.Rmd`, `decontX.Rmd`.
- **NOTE**: coding-practice, with exact locations —
  `sapply()` at `R/findMarkersTree.R:2414` (use `vapply()`);
  `1:...` idioms in `R/findMarkersTree.R` (use `seq_len()`/`seq_along()`);
  `cat()` outside a `show` method at `R/findMarkersTree.R:1895`;
  `paste` in condition signals at `R/celda_functions.R:390` and
  `R/findMarkersTree.R:281`;
  `suppressWarnings()` at `R/moduleHeatmap.R:352` (4 occurrences total).
- **NOTE**: 108 functions exceed the recommended 50-line length. The 5
  longest are all in `R/findMarkersTree.R` (726 lines) and
  `R/recursiveSplit.R` (401 lines) — `findMarkersTree.R` is the single
  densest file for coding-practice cleanup.
- **NOTE**: consider runnable examples on exported man pages; prefer
  `donttest{}` over `dontrun{}` (`celdaGridSearch.Rd`, `reportceldaCG.Rd`).
- **NOTE**: style — 102 lines >80 chars, 6081 lines (26%) not indented in
  multiples of 4 spaces (matches what `make lint` already surfaces).
- Remaining NOTE is informational (bioc-devel mailing list subscription
  can't be verified without admin credentials).

### `BiocCheck::BiocCheckGitClone(".")`: 1 real NOTE

- No `CITATION` file. Optional — only worth adding if there's a
  preprint/publication to cite for celda.
- (Its "2 ERRORS" — a stray `.DS_Store` picked up from the host's `dev/`
  folder and BiocCheck's own output folder — were artifacts of running the
  tool inside a temporary container copy, confirmed via `git ls-files` that
  neither path is actually tracked. Not a real repo problem; cleaned up.)

### `make test`: FAIL 0 | WARN 21 | SKIP 2 | PASS 101

No test failures — celda's suite passes cleanly. All 21 warnings are
deprecation notices, not correctness problems:

- **20 warnings** trace to **`scuttle`/`scater`** (transitive dependencies),
  triggered via `scater::logNormCounts()` in `R/decon.R:1067`, called from
  `.decontxInitializeZ()`:
  - `scuttle::librarySizeFactors()` is deprecated → use
    `scrapper::centerSizeFactors()` (10 occurrences)
  - `scuttle::normalizeCounts()` is deprecated → use
    `scrapper::normalizeCounts()` (10 occurrences)

  All 20 fire during `tests/testthat/test-decon.R` (DecontX tests) — every
  DecontX call currently exercises two newly-deprecated Bioconductor
  functions under current dependency versions. Since celda doesn't call
  `librarySizeFactors`/`normalizeCounts` directly, the fix is
  upstream-dependent: `scater::logNormCounts()` needs a newer call pattern,
  or celda needs to route around it (e.g. call
  `scrapper::centerSizeFactors()`/`scrapper::normalizeCounts()` directly)
  once `scater` itself updates its internals.
- **1 warning**: `ggplot2::aes_string()` (`test-celda_C.R:98`) — see below.

### Lifecycle / deprecated-function sweep

- **`ggplot2::aes_string()`** — deprecated since ggplot2 3.0.0 (2018),
  superseded by tidy `aes()` evaluation. Used **36 times** across 3 files:
  `R/perplexity.R` (22), `R/plot_dr.R` (11), `R/plot_decontx.R` (3). The
  installed ggplot2 (4.0.3) still only warns, but this is old enough that a
  future ggplot2 release could remove it outright, which would break most of
  celda's plotting functions at once. This is the single highest-value,
  celda-code-only finding — worth its own issue (and possibly an ADR if the
  fix touches a lot of `perplexity.R`/`plot_dr.R` at once, since the
  playbook's scope guard treats file-wide edits as needing one).
- No archived or CRAN-removed packages found among celda's direct
  `Imports`/`Depends`/`Suggests`. `plyr` and `reshape2` are both in
  tidyverse "superseded" status (long-term maintenance only, not
  recommended for new code) but still on CRAN and not throwing deprecation
  warnings — low-priority modernization candidate, not urgent.
- `BiocManager::valid()` (clean `RELEASE_3_23` image, no celda installed):
  8 packages out-of-date, 0 too new — normal drift in the base Docker image
  between rebuilds, not a celda-specific dependency-pinning issue.

### `pkgdown::check_pkgdown()` and `lintr::lint_package()`

- `pkgdown::check_pkgdown()` still errors: `In _pkgdown.yml, url is missing.`
  A real, pre-existing config gap (add a `url:` field to `_pkgdown.yml`).
- `lintr::lint_package()`: 5339 total lints (line length, indentation,
  naming) across `R/` and `vignettes/` — same categories `make lint`
  already surfaces, at full-package scale.

## Suggested next steps

1. File a GitHub issue for the `aes_string()` migration (36 occurrences,
   `R/perplexity.R`/`R/plot_dr.R`/`R/plot_decontx.R`).
2. File a **separate** GitHub issue for the `scater`/`scuttle` deprecation
   in the DecontX path (`R/decon.R`) — this one may need to wait on
   upstream `scater` before celda's own code can change; worth flagging to
   the scater/scuttle maintainers too if it's not already tracked there.
3. File a GitHub issue (or a short PR) for the quick `DESCRIPTION` fixes:
   `URL` field, `R (>= 4.6.0)`, `LazyData` setting, `fnd` role, vignette
   chunk labels.
4. Leave style/line-length/function-length NOTEs for a later, dedicated
   cleanup pass — `R/findMarkersTree.R` is the single densest file for this.
5. Re-run this same current-release audit periodically (e.g. from the
   weekly `BioC-check` cron once it's promoted to required) rather than
   relying on whatever Bioconductor version a given contributor's local R
   happens to resolve to — that's exactly what hid the `scater`/`scuttle`
   finding until this run.
