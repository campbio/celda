# findMarkersTree: top-level markers and informational forced markers (#450)

## Problem

With metaclusters, many cells are predicted as "MISSING". On the #410
data (40,349 cells, 48 clusters, 96 modules) accuracy is 0.460 with five
automatic metaclusters against 0.585 without.

Two mechanisms, both traced on that data:

1. **Top level.** At 0.9 the one-off split marks M1, M3, M4 and M5. M2 has
   no marker. `altSplit` forces one (L19 >= 73, stat 0.76, no threshold
   check), `.allTopMarked()` counts it, and the tree is accepted. In
   `.predictClass()` a class with an up rule matches on its top-stat up
   rule alone and its down rules are ignored, so M2 becomes "L19 >= 73":
   1,017 M4 cells also pass (two matches, so MISSING) and 22% of M2 cells
   do not (no match). Without the forced marker, M2 defined by the absence
   of the other markers, top-level routing of M2 goes from 77% to 94% and
   of M4 from 66% to 96%. With `altSplit = FALSE` the threshold loop lowers
   to 0.729 and picks L19 as a real marker, so the same weak cut returns.
2. **Subtrees.** Forced markers inside metacluster subtrees do the same to
   down-only clusters: with only the top level fixed, defaults give 0.495
   (M1 2,177 and M2 2,275 MISSING); with subtree forced markers off, 0.557
   (M2 5,851 correct / 203 MISSING, M4 2,717 / 32).

## Goal

The top level finds the metaclusters that can be predicted with one set of
mutually exclusive markers at the user's threshold. Metaclusters without
such a marker are defined by being low on all of them and are separated by
the rest of the procedure. Forced markers stay visible but no longer decide
predictions.

## A. Top level of a metacluster tree

In the metacluster branch of `findMarkersTree()`:

1. Build the top-level tree with `.generateTreeList()` at `threshold`. No
   alternative split at the top level, whatever `altSplit` is. The
   threshold is never lowered.
2. If the top-level split is not a one-off split (no metacluster has a
   marker), stop building with metaclusters: `message()` that the tree is
   built without metaclusters, and return the same result as
   `autoMetaclusters = FALSE` with no `metaclusters`. The no-metacluster
   code moves into a helper, `.findMarkersTreeFlat()`, called from both
   places.
3. Otherwise get the unmarked metaclusters with `.unmarkedMetaclusters()`.
   With two or more, merge them with `.mergeMetaclusters()` (named by
   joining with "+", as now) and go back to step 1 with the merged labels.
   Each merge lowers the number of metaclusters, so the loop ends. A merge
   that would leave a single metacluster is never needed: a one-off split
   marks at least one metacluster.
4. With at most one unmarked metacluster, keep the tree. Its single
   level is the one-off split; the unmarked metacluster has only down
   rules there.
5. Messages, not warnings, as this is now a normal outcome:
   - one unmarked: "Metacluster M2 has no positive marker at the top
     level. It is identified by the absence of the other metaclusters'
     markers."
   - merged: "Metaclusters M2, M4 have no positive marker at the top
     level. They were merged into M2+M4 and are separated in its own
     subtree."
   The names use `.metaclusterDisplayNames()`, as now.
6. Removed, as they only served the old loop: threshold lowering,
   `tmpThreshold`, `markerThreshold` and the `markersToRemove` step,
   `buildTopTree`, the balanced-split merge, and the functions
   `.canAddTopMarker()`, `.allTopMarked()`, `.balancedMergeClasses()` and
   `.mergeUnmarkedMetaclusters()`. `.unmarkedMetaclusters()` and
   `.mergeMetaclusters()` stay.

7. With `featureLabels` and `counts`, the gene-level rules drop the
   top-level down rows only for metaclusters that have a top-level up row.
   A metacluster without a marker keeps its down rows, which describe it.
   Before, its rules became empty and the run failed with "missing value
   where TRUE/FALSE needed"; this was rare on devel, where every
   metacluster was given a marker. (Added after the spec review.)

Subtrees are built as now, with `altSplit` applied there.

## B. Forced markers are informational

Applies wherever `.addAlternativeSplit()` runs: metacluster subtrees and
trees without metaclusters.

1. `.addAlternativeSplitClass()` sets `alternative = TRUE` on the split it
   adds.
2. `.mapClass2features()` adds a logical column `alternative` to each
   rules table: TRUE for rows from a forced split (the forced class's up
   row and the other classes' down rows), FALSE otherwise.
3. `.predictClass()` drops rows with `alternative == TRUE` before checking
   a level, when the column exists. Rules without the column (older saved
   trees) are handled as before.
4. So predictions equal those with `altSplit = FALSE`; the rules and the
   dendrogram still show the forced marker. `statUsed` is unchanged, so
   branch points and plots are unaffected.

## Docs

- `@param threshold`: drop the lowering text; say metaclusters without a
  top-level marker at the threshold are identified by the absence of the
  others' markers, and two or more are merged ("+") and separated in their
  own subtree; with none marked, the tree is built without metaclusters.
- `@param altSplit`: the forced marker is shown in the rules and the
  dendrogram but not used for prediction; the top level of a metacluster
  tree never gets one.
- `@return`: document `alternative`; fix "five columns".
- NEWS entry under the development version; bump to 1.29.14.

## Tests (`tests/testthat/test-findMarkersTree.R`)

New, each failing on its assertion on devel:

- #450 fixture with metaclusters: two strongly marked metaclusters and a
  third with no marker, plus a weak feature high in the third and in part
  of a marked one. For `altSplit` TRUE and FALSE: the third has no
  top-level up rule; its cells and the overlapping metacluster's cells are
  predicted correctly, not MISSING; no "Decreasing classifier threshold"
  message.
- Forced marker informational: a tree that gets a forced marker
  (the existing `f3Class3Only` fixture) has the same `$prediction` with
  `altSplit` TRUE and FALSE, and the forced row has `alternative = TRUE`.
- `getDecisions()` on rules without the `alternative` column gives the
  same predictions as before.
- Fall-back: metaclusters where no feature marks any one of them gives the
  message and the same rules as the tree without metaclusters.

Changed existing tests:

- Expected warning becomes an expected message: "merges metaclusters left
  without markers", "merges a multi-cluster metacluster", "keeps a single
  metacluster without markers", "splits a metacluster with one feature
  left", "merged metaclusters keep names with parentheses", "stops
  lowering the threshold with reuseFeatures" (renamed: it checks the
  merge with `reuseFeatures = TRUE`).
- Rewritten: the two #443 tests "altSplit = FALSE applies to the top-level
  tree" and "altSplit = FALSE does not try a top-level alternative split"
  become "the top level never gets a forced marker" for both settings.
- Replaced: "metaclusters are not all merged into one" by the fall-back
  test.
- Deleted with their functions: "a balanced split is merged only when that
  leaves two classes" (`.balancedMergeClasses()`) and "a balanced
  top-level split's feature can still become a marker"
  (`.canAddTopMarker()`).

## Checks

`make test`, `make check`. Then the #410 data with defaults: expect about
0.557 accuracy (0.460 on devel) and few MISSING in M2 and M4; report the
numbers in the PR.

## Effect on users

Trees, rules and predictions change for every metacluster call, and
predictions change for any tree with a forced marker. Rules tables gain a
column. singleCellTK does not call `findMarkersTree()`, `getDecisions()`,
`plotDendro()` or `plotMarkerHeatmap()` (checked in a local clone at
04110f3c, 2026-10-09).
