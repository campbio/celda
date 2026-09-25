# Architecture Decision Records (ADRs)

This directory records decisions for celda that are hard to reverse or that
future contributors will reasonably question: dependency changes, module/file
structure, class design, and build/deploy machinery. It is **not** for
routine choices (a bug fix, a new argument, a small helper function).

Decisions are never edited after acceptance. A reversal is a new ADR that
marks the old one as superseded.

Start from [`template.md`](template.md) when proposing one. Number files
sequentially: `000N-short-title.md`.

## Index

| # | Title | Status | Date |
|---|-------|--------|------|
| 0001 | [Record architecture decisions](0001-record-architecture-decisions.md) | Accepted | 2026-09-18 |
| 0002 | [Migrate DecontX implementation to the standalone decontX package](0002-migrate-decontx-to-standalone-package.md) | Accepted | 2026-09-22 |
| 0003 | [Declare Seurat as a Suggests dependency](0003-declare-seurat-as-suggests.md) | Accepted | 2026-09-23 |
