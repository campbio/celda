# Updating your package to playbook v2 + security hardening

For packages that already have the original (v1) AI-agent scaffolding.
Work on a branch off `devel`. All steps are file moves/config — no code
changes. Reference implementation: decontX `ai_dev` branch, commits
`b5ebb1b` (v2 reorg) and `14e4958` (hardening).

## 1. Move files (use `git mv` so history follows)

```bash
git mv CONTRIBUTING.md SECURITY.md .github/
mkdir -p dev
git mv RELEASE.md ROADMAP.md AUDIT.md dev/
git mv docs/adr dev/adr        # docs/ is now reserved for pkgdown output only
```

## 2. Add `dev/adr/README.md`

Purpose + index table (number / title / status / date) + the granularity
rule: ADRs are for decisions that are hard to reverse or that future
contributors will question — copy it from decontX
(`dev/adr/README.md` on the `ai_dev` branch).

## 3. Fix cross-references

Grep your scaffolding for old paths and update: `RELEASE.md` →
`dev/RELEASE.md`, `docs/adr` → `dev/adr`, `agent-log.md` →
`dev/agent-log.md`, `SECURITY.md` → `.github/SECURITY.md`. Also bump the
AGENTS.md common section header to **v2.0** and replace its Safety-rules
bullet with the v2 wording ("Architectural decisions are recorded in
dev/adr/… Never store anything in docs/…").

```bash
grep -rn "docs/adr\|RELEASE\.md\|ROADMAP\.md\|AUDIT\.md\|agent-log" \
  AGENTS.md .github dev .claude
```

## 4. Consolidate `.Rbuildignore`

Delete the per-file entries (`^CONTRIBUTING\.md$`, `^SECURITY\.md$`,
`^RELEASE\.md$`, `^ROADMAP\.md$`, `^AUDIT\.md$`, `^agent-log\.md$`) and
make sure these exist:

```
^Makefile$  ^AGENTS\.md$  ^CLAUDE\.md$  ^GEMINI\.md$
^dev$  ^docs$  ^\.claude$  ^\.lintr$  ^\.github$
```

(One entry per line in the actual file.)

## 5. Harden `.claude/settings.json`

- **allow**: remove `Bash(Rscript:*)` and `Bash(make:*)`; list your
  Makefile targets individually instead — `Bash(make test)`,
  `Bash(make check)`, `Bash(make bioccheck)`, `Bash(make docs)`,
  `Bash(make lint)`, `Bash(make build)` — plus read-only git
  (`git status/diff/log/show/branch/blame`). Do **not** pre-approve
  `make clean`, `make site`, or `make site-deploy`.
  (Why: `Rscript *` allows arbitrary R code, which defeats every deny
  rule; site deploys are a maintainer action.)
- **deny**: add `Edit(docs/**)`, `Write(docs/**)`, `Bash(rm -rf:*)`,
  `Bash(rm -r:*)`, `Bash(git clean:*)` alongside the existing `man/**`,
  `NAMESPACE`, force-push, and generated-file rules.
- **hook**: replace the styler auto-rewrite hook with report-only lint:
  copy `dev/hooks/lint-changed.sh` from decontX, `chmod +x` it, and
  point the PostToolUse command at `bash dev/hooks/lint-changed.sh`.
  (Why: auto-styling touched files buries real changes in formatting
  diffs when a lint backlog exists.)

## 6. Verify, then PR

```bash
python3 -m json.tool .claude/settings.json    # settings parse
make lint && make docs && make test           # still green; revert any
                                              # RoxygenNote / stanExports churn
R CMD build --no-build-vignettes . \
  && tar -tzf *_*.tar.gz | grep -iE "AGENTS|Makefile|dev/|\.claude|\.lintr" \
  || echo "tarball clean"
git status -- R src tests DESCRIPTION NAMESPACE   # must be empty
```

Commit and open a PR against `devel`.
