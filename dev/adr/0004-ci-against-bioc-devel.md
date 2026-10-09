# 4. Run GitHub CI against Bioconductor devel

Status: Proposed

Date: 2026-10-09. Needs the maintainer's approval.

## Context

The `devel` branch is the Bioconductor devel version of celda, but GitHub
CI installed the Bioconductor *release* packages (Bioc 3.23 for R 4.6).
pak, which `r-lib/actions/setup-r-dependencies` uses, picks the
Bioconductor version that matches R unless told otherwise, so this was true
even of the job running in the `bioconductor/bioconductor_docker:devel`
container. CI therefore tested devel code against older dependencies than
the Bioconductor builders use. It failed on warnings the builders don't
see: decontX 1.10 calls scuttle functions deprecated in scuttle 1.22
(`normalizeCounts`, `librarySizeFactors`), which makes R CMD check fail on
the `decontX` example. decontX 1.11 in devel uses scrapper instead.

The macOS jobs never reached the check: `alabaster.base` (pulled in through
singleCellTK in Suggests) has no macOS binary, and its source build fails
to link OpenSSL. The BiocCheck job ran on macOS, so it failed the same way.

singleCellTK hit the same problems and made the same changes (its ADR
0007).

## Decision

- Set `R_BIOC_VERSION` to the current Bioconductor devel version ("3.24")
  in `check-standard.yaml` (both jobs) and `BioC-check.yaml`, so pak
  installs Bioconductor devel packages.
- Before installing dependencies on macOS, install Homebrew's `openssl@3`
  and add its `lib` and `include` directories to `LDFLAGS` and `CPPFLAGS`
  in `~/.R/Makevars`.
- Run BiocCheck on `ubuntu-latest` and install it with the other
  dependencies (`bioc::BiocCheck`), so it also comes from devel.

## Consequences

- CI tests what the Bioconductor devel builders build, so a green CI
  should predict a clean build report.
- `R_BIOC_VERSION` must be updated at each Bioconductor release (April and
  October); `dev/RELEASE.md` has the step. If it is forgotten, CI silently
  tests against an old devel.
- The first CI runs build more packages from source until the cache is
  warm.
- The macOS OpenSSL step is a workaround. Remove it once `alabaster.base`
  has a macOS binary.
- The coverage, lint, and pkgdown workflows still use the release
  packages. They pass today; move them to devel too if they start to
  differ from the builders.

## Alternatives considered

- **Keep testing against release.** Rejected: it tests the wrong
  dependencies and fails on warnings that are already fixed in devel.
- **Muffle decontX's deprecation warnings in celda's wrapper.** Rejected:
  it would hide warnings that devel no longer raises, and it changes
  package code to work around a CI setting.
- **Drop macOS from CI.** Rejected: one of the Bioconductor builders is
  macOS, and the OpenSSL fix is small.
