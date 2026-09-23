# 3. Declare Seurat as a Suggests dependency

Status: Accepted

Approved by the maintainer (2026-09-23).

## Context

`findMarkersTree()` accepts a `seurat` argument (a documented, supported input
type) and calls `Seurat::Idents()` and `Seurat::RunUMAP()` on it. Seurat was
never declared in `DESCRIPTION`, so these were bare calls to an undeclared
package. `R CMD check` flagged them as "no visible global function definition
for 'Idents' / 'RunUMAP'", and the code would fail at runtime if a user passed
a Seurat object without Seurat installed via some other path.

## Decision

Declare **Seurat** in `Suggests` (not Imports — it is only used for the
optional Seurat-object input path) and call it namespaced as `Seurat::Idents`
and `Seurat::RunUMAP`. The calls are already guarded by
`methods::hasArg(seurat)`, so they only run when the user supplies a Seurat
object (which implies Seurat is available).

## Consequences

- Clears the remaining R CMD check "no visible global function" NOTE.
- The Seurat-object input path to `findMarkersTree()` now has a declared,
  namespaced dependency and fails with a clear error if Seurat is absent.
- Adds Seurat (a large package) to the check/CI install surface, but only as a
  Suggests, so it does not affect a standard `install("celda")`.
