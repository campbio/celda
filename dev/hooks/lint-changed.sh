#!/usr/bin/env bash
# Report-only lint on a single touched file. Does NOT rewrite anything —
# auto-styling buries real changes in formatting diffs while there's a lint
# backlog (see dev/AUDIT.md). Run manually with `make lint` for the full
# package, or `styler::style_file()` yourself if you want a file reformatted.
set -euo pipefail

file="${1:-${CLAUDE_TOOL_INPUT_FILE_PATH:-}}"

case "$file" in
  *.R)
    Rscript -e "lintr::lint('$file')"
    ;;
esac
