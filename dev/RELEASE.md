# Release checklist

celda follows the Bioconductor release cadence: releases in about April and
October, each with a freeze shortly before. Versions are x.y.z: `y` is odd
on `devel` and even on the release branch. Bioconductor makes the `y` bumps
at release time; you bump `z` by 1 for every change that goes to the
Bioconductor git server, and never change `x` or `y` yourself (see the
shared standards).

Remotes: `bioc` is `git@git.bioconductor.org:packages/celda.git`; the
shared GitHub repo is `campbio` in the maintainer's clone. Check with
`git remote -v` before any push.

## Before the freeze

1. Sync with Bioconductor: `git fetch bioc`, then on `devel`
   `git merge bioc/devel`, and push `devel` to GitHub.
2. Confirm the `Version:` field in `DESCRIPTION` follows the even/odd rule
   for the branch.
3. Run `make check-full` and `make bioccheck`. Triage every WARNING and NOTE
   into GitHub issues rather than silent fixes, then fix, re-run, and open
   PRs for the fixes. `make bioccheck` prints the tarball size; the limit is
   10 MB (5 MB per file).
4. Make sure `NEWS.md` on `devel` has a section for the version going into
   the release, and that it's pushed to Bioconductor before the release
   branch is cut, so the release branch has it too.

## On release day

Bioconductor creates `RELEASE_X_Y` and bumps devel. Bring both to GitHub,
so the stable-branch sync (`.github/workflows/sync-stable.yaml`) can update
`master` and tag the release:

```bash
git fetch bioc
git checkout devel && git merge bioc/devel && git push campbio devel
git push campbio bioc/RELEASE_X_Y:refs/heads/RELEASE_X_Y
```

The sync runs daily, or start it from the Actions tab. It picks up the new
release once https://bioconductor.org/config.yaml names it.

Then set `R_BIOC_VERSION` in `.github/workflows/check-standard.yaml` and
`BioC-check.yaml` to the new devel version (the `devel_version` in
config.yaml), so CI keeps testing against Bioconductor devel (ADR 0004).
Open the change as a PR against `devel`.

## After release

- Check the Bioconductor build report for celda
  (https://bioconductor.org/checkResults/) across all platforms. It is the
  authoritative status; CI is only an early warning.
- File issues for any platform-specific failures that CI didn't catch.

## pkgdown site

The public site is https://www.camplab.net/celda. `docs/` is committed on
`master`, but this repo has no GitHub Pages site of its own (checked
2026-09-27), so find out where camplab.net serves the site from before
changing how `docs/` is published. Because `master` is now an automatic
copy of the current release branch, `master`'s `docs/` is the release
branch's `docs/`. The sync's first run replaces `master`'s old `docs/`
(from 1.18.2) with the release branch's, so confirm where the site is
served from before that first run.

Planned migration to a `gh-pages` deployment (a maintainer action, with a
person at the keyboard):

1. Confirm where the site is served from today.
2. `make site` locally and review the output.
3. `make site-deploy` to push the rendered site to `gh-pages`.
4. Point the site at `gh-pages` and confirm the live URL still serves
   (reversible by pointing it back).
5. Remove the committed `docs/` from the branches and add `docs/` to
   `.gitignore` (leave prior history alone).

CI only runs `pkgdown::check_pkgdown()` (`make site-check` locally); it
never builds or deploys the site.
