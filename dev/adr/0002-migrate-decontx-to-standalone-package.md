# 2. Migrate DecontX implementation to the standalone decontX package

Status: Accepted

Proposed in campbio/celda#414. Approved by the maintainer (2026-09-22).

## Context

celda has historically bundled the full **DecontX** ambient-RNA
decontamination method: the R orchestration and EM helpers in `R/decon.R`, a
C++ backend in `src/DecontX.cpp`, plotting in `R/plot_decontx.R`, and
simulation utilities. The method has since been extracted into a dedicated
Bioconductor package, `decontX`, which is now the primary home for the
algorithm. Maintaining two copies of the same numerical code invites drift and
double the maintenance/review burden, and celda's copy is the source of ~20
`scater`/`scuttle` deprecation warnings surfaced during the release audit
(`dev/agent-log.md`).

A complication constrains the dependency direction: the `decontX` package
**Imports celda** because its EM routines call `celda::normalizeCounts()`.
Bioconductor forbids cycles in the strong dependency graph
(Depends/Imports/LinkingTo), so celda **cannot** add `decontX` to `Imports` —
that would create `celda -> decontX -> celda`.

## Decision

Make celda's DecontX surface a set of **thin, backward-compatible wrappers**
around the `decontX` package, so the algorithm lives in exactly one place:

- Keep exporting `decontX`, `decontXcounts`/`decontXcounts<-`,
  `plotDecontXContamination`, `plotDecontXMarkerExpression`,
  `plotDecontXMarkerPercentage`, and `simulateContamination`. Each (except the
  `decontXcounts` accessors) delegates to the corresponding `decontX::` function
  via `...` passthrough, behind a `requireNamespace("decontX")` guard.
- The `decontXcounts`/`decontXcounts<-` accessors stay self-contained in celda
  (they are plain `assay(object, "decontXcounts")` accessors and need no
  dependency).
- List `decontX` in **Suggests** (not Imports) to respect the circular-
  dependency rule. When `decontX` is not installed, the wrappers error with an
  actionable install message.
- Delete celda's DecontX algorithm internals and `src/DecontX.cpp`; regenerate
  `RcppExports` with `Rcpp::compileAttributes()`.
- Drop DecontX-only `Imports` that nothing else in celda uses
  (`MCMCprecision`; `DelayedArray` if unused). Retain `uwot`, `dbscan`,
  `Rtsne`, and `RcppEigen`/`LinkingTo` (used by non-DecontX celda code).

Once the `decontX` package removes its `celda::normalizeCounts` dependency, a
follow-up ADR may promote celda's dependency from Suggests to Imports and
re-export the `decontX` symbols directly instead of wrapping.

## Consequences

- **Easier:** one implementation to maintain and review; smaller celda C++
  footprint and build; removal of the DecontX-originated deprecation warnings.
- **Harder / follow-up:** celda's DecontX functions now require the `decontX`
  package at runtime (Suggests-guarded). Two packages temporarily export the
  same generics (`decontX`, `decontXcounts`), producing harmless masking
  messages when both are attached; this resolves when celda later re-exports
  instead of wrapping. The vignettes now require `decontX` to be installed to
  knit.
