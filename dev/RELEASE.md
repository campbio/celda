# Release checklist

celda follows the Bioconductor release cadence: devel version bumps happen
~April and ~October, with each Bioconductor release freeze shortly before.
Bioconductor uses an even/odd `x.y.z` minor-version scheme — `y` is odd on
`devel`, even on `release`; `z` resets to `0`/`1` at each bump (see
`NEWS.md` for the pattern, e.g. the `1.22.0` "match Bioconductor release
version" entries).

## Before the freeze

1. Sync with the upstream Bioconductor git mirror
   (`git.bioconductor.org/packages/celda`) if not already configured as a
   remote — celda currently only has the GitHub `origin` remote; add
   `git.bioconductor.org` here when first needed.
2. Confirm the `Version:` field in `DESCRIPTION` follows the even/odd rule for
   the branch you're releasing.
3. Run `make check` and `make bioccheck`; triage any WARNING/NOTE into a fix
   plan (as GitHub issues, not silent fixes) — fix, re-run, PR the fixes.
4. Update `NEWS.md` with a summary of user-facing changes since the last
   release.

## After release

- Check the Bioconductor build report for celda
  (bioconductor.org/checkResults) across all platforms.
- File issues for any platform-specific failures that CI didn't catch.
- Bump the devel branch version per the even/odd scheme for the next cycle.

## pkgdown site

celda's `docs/` is currently committed directly to the branch that GitHub
Pages serves, rather than deployed to a dedicated `gh-pages` branch. Planned
migration (do this deliberately, with a human at the keyboard — see
`AGENTS.md` / the setup playbook for the full procedure):

1. Check the repo's Settings → Pages source.
2. `make site` locally and review the output.
3. `pkgdown::deploy_to_branch()` to push the rendered site to `gh-pages`.
4. Flip Settings → Pages source to `gh-pages` root; confirm the live URL
   still serves (reversible by flipping back).
5. Remove the committed `docs/` from the main branch and add `docs/` to
   `.gitignore` (leave prior history alone).

CI only runs `pkgdown::check_pkgdown()` (structure check) — it never builds
or deploys the site. Site deploys are a maintainer-only local action
(`make site-deploy`).
