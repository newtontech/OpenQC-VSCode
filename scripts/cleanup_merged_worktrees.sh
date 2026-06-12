#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

git fetch --all --prune

if [ ! -d .worktrees ]; then
  echo "No .worktrees directory; pruning git metadata only."
  git worktree prune
  git remote prune origin || true
  exit 0
fi

for wt in .worktrees/*; do
  [ -d "$wt" ] || continue
  branch="$(git -C "$wt" branch --show-current 2>/dev/null || true)"
  [ -n "$branch" ] || continue

  pr_json="$(gh pr view "$branch" --json state,mergedAt,headRefName 2>/dev/null || true)"
  state="$(printf '%s\n' "$pr_json" | python -c 'import json,sys; data=sys.stdin.read().strip(); print(json.loads(data).get("state","") if data else "")' 2>/dev/null || true)"

  if [ "$state" = "MERGED" ]; then
    echo "Cleaning merged worktree: $wt ($branch)"
    git worktree remove --force "$wt" || true
    git branch -d "$branch" || git branch -D "$branch" || true
  fi
done

git worktree prune
git remote prune origin || true
