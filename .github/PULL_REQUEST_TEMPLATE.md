<!-- Base branch: devel. Use RELEASE_X_Y only for an approved release fix.
     Never target main/master, which is updated automatically. -->

## What changed and why

<!-- Link the issue if there is one. -->

## How it was tested

## Checklist

- [ ] Tests added or updated, `make test` passes, and `make coverage`
      didn't drop
- [ ] `make check-full` and `make bioccheck` pass with no new errors or warnings
- [ ] `make docs` run, if roxygen comments changed
- [ ] New exports added to `_pkgdown.yml`, and `make site-check` passes
- [ ] NEWS.md updated for user-facing changes
- [ ] Version bumped (z) if this will be pushed to Bioconductor
- [ ] Plan review and `/code-review` run; findings fixed or answered
- [ ] Related issue linked

## ADR

<!-- Link to dev/adr/NNNN-*.md for a structural change (file splits,
     dependencies, class redesign), or write "N/A". -->

## Scientific correctness

<!-- REQUIRES HUMAN JUDGMENT. An agent must never fill this in.
     If this PR changes the celda models (celda_C, celda_G, celda_CG), the
     Gibbs sampling or EM steps, clustering or module results, or what a
     plot or reported metric shows, say who verified the output is still
     scientifically correct and how. Passing tests is not enough.
     Otherwise write "No change to results". -->

## Generated content

<!-- If an AI agent wrote part of this PR, say which parts, so reviewers
     know where to focus. -->
