# AGENTS.md

## Campbell Lab Playbook (common across lab packages — v2.0, do not edit per-repo)

### Common commands
make test / make check / make bioccheck / make docs / make lint / make site
(See Makefile for definitions. These are the ONLY sanctioned entry points.
Run `make test` after every change; `make check` before opening a PR.)

### Git and PR workflow
- Branch from devel; all work lands via PR. Never push to devel or master.
- Use plan mode for any non-trivial change.
- Run /code-review before requesting human review.
- Every user-facing change gets a NEWS.md entry.

### Coding conventions
- Style enforced by lintr/styler (config in repo); <= 80-char lines (BiocCheck).
- roxygen2 owns man/ and NAMESPACE — NEVER hand-edit them.
- Use accessor functions, not @ slot access, outside class definition files.

### Documentation (pkgdown)
- The website is GENERATED. Improve docs by editing roxygen comments, vignettes,
  and _pkgdown.yml — never files under docs/ or the gh-pages branch.
- New exported functions MUST be added to the _pkgdown.yml reference index;
  verify with pkgdown::check_pkgdown().
- To preview one changed page: pkgdown::build_article("<name>") or
  build_reference_index(). NEVER run a full build_site() as verification —
  full site builds/deploys are a local maintainer action (make site-deploy).
- Files under vignettes/articles/ are pkgdown-only and NOT checked by
  R CMD check — knit locally when you edit them.

### Shiny app rules (packages with inst/shiny only)
- The app contains NO analysis logic. Server code only wires inputs to
  exported package functions and renders results. New app features are
  implemented as tested, exported functions first.
- Reactive logic is tested with shiny::testServer(); the golden path is
  covered by a small shinytest2 smoke suite (make test-app).
- UI changes are verified with a screenshot of the RUNNING app
  (make app + browser), not just passing tests.
- inst/ code is invisible to R CMD check — tests and lintr are the only
  guards; inst/shiny is included in the lint paths.

### Versioning and releases
- Bioconductor even/odd x.y.z scheme; releases ~April and ~October.
- Follow dev/RELEASE.md for the release checklist.

### Safety rules
- No structural refactors (file splits, DESCRIPTION dependency changes,
  class redesign) without an approved ADR — propose via a GitHub issue.
- Never commit secrets, tokens, or absolute local paths.
- Architectural decisions are recorded in dev/adr/ (see template and index
  there). Never store anything in docs/ — that is pkgdown build output.
- Maintainer docs (release, roadmap, audits) live in dev/, not the root.

## This package: celda

### Project overview
celda is a suite of Bayesian hierarchical models for clustering single-cell
RNA-seq data, able to bi-cluster genes into modules and cells into
subpopulations simultaneously. It also includes DecontX, a Bayesian method to
estimate and remove ambient RNA contamination from droplet-based scRNA-seq
without requiring empty-droplet data.

### Repository map
- Class/model definitions: `R/aaa.R` (S4 classes), `R/celda_C.R`, `R/celda_G.R`,
  `R/celda_CG.R`, `R/celdaGridSearch.R`, `R/decon.R`, `R/accessors.R`
- Plotting: `R/celda_heatmap.R`, `R/plotHeatmap.R`, `R/plot_dr.R`,
  `R/plot_decontx.R`, `R/moduleHeatmap.R`, `R/celdaProbabilityMap.R`,
  `R/semi_pheatmap.R`, `R/elbow.R`
- Dimensionality reduction: `R/celdatSNE.R`, `R/celdaUMAP.R`
- Data import/simulation: `R/simulateCells.R`, `R/celdatosce.R`, `R/data.R`,
  precomputed objects in `data/*.rda`
- Vignettes: `vignettes/celda.Rmd`, `vignettes/decontX.Rmd`,
  `vignettes/articles/`; report templates in `inst/rmarkdown/`

### Object model
S4 classes rooted in `celdaModel` (slots: `params`, `names`, `completeLogLik`,
`finalLogLik`, `clusters`), extended by `celda_C` (cell clustering, adds
`sampleLabel`), `celda_G` (feature/gene-module clustering), and `celda_CG`
(bi-clustering, `contains = c("celda_C", "celda_G")`). `celdaList` holds grid
search results (`runParams`, `resList`, `countChecksum`, `perplexity`). Most
user-facing workflows wrap results into a `SingleCellExperiment` rather than
manipulating the S4 model objects directly. `R/accessors.R` defines get/set
generics dispatching on either `SingleCellExperiment` or `celdaModel`
(`celdaClusters`, `celdaModules`, `sampleLabel`, `params`, `matrixNames`,
`runParams`, `resList`, `celdaModel`, `celdaPerplexity`, `countChecksum`).
DecontX (`R/decon.R`) defines its own generics (`decontX`, `decontXcounts`)
operating on `SingleCellExperiment` or matrix-like input. Use these accessors,
not `@` slot access, outside the class definition files.

### Environment setup
R >= 4.0. `BiocManager::install("celda", dependencies = TRUE)`. Uses Rcpp/
RcppEigen (`LinkingTo: Rcpp, RcppEigen` in DESCRIPTION) — a C++ toolchain is
required to build from source.

### Package-specific notes
- No Shiny app in this package (`inst/` only holds `rmarkdown/` report
  templates) — the Shiny app rules above do not apply here.
- Rcpp/C++ (`src/*.cpp`, `src/*.c`) backs the Gibbs sampling helpers
  (`cG_calcGibbsProbY.cpp`), DecontX's EM/log-likelihood routines
  (`DecontX.cpp`), and matrix utilities (`eigenMatMultInt.cpp`,
  `matrixNorm.cpp`, `matrixSums.c`, `matrixSumsSparse.cpp`, `perplexity.c`).
  `R/RcppExports.R` is auto-generated (`Rcpp::compileAttributes()`) — never
  hand-edit it.
- Slow tests: `tests/testthat/test-celda_G.R` and `test-celda_C.R` run real
  Gibbs sampling via `simulateCells()` + `celdaGridSearch()` and are the
  slowest part of `make test`.
- pkgdown site (`docs/`) is currently committed directly to the main branch
  rather than served from `gh-pages`; see `dev/RELEASE.md` for the planned
  migration.
