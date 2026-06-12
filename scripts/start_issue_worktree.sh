#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <issue-number> <short-slug> [base-branch]" >&2
  exit 2
fi

issue="$1"
slug="$2"
base="${3:-}"

root="$(git rev-parse --show-toplevel)"
cd "$root"

if [ -z "$base" ]; then
  base="$(git symbolic-ref --quiet --short refs/remotes/origin/HEAD 2>/dev/null | sed 's#^origin/##' || true)"
  base="${base:-$(git branch --show-current)}"
fi

branch="agent/issue-${issue}-${slug}"
worktree=".worktrees/issue-${issue}-${slug}"

git fetch origin "$base"
mkdir -p .worktrees

if git show-ref --verify --quiet "refs/heads/${branch}"; then
  echo "Branch exists: $branch"
else
  git branch "$branch" "origin/$base"
fi

if [ -d "$worktree" ]; then
  echo "Worktree exists: $worktree"
else
  git worktree add "$worktree" "$branch"
fi

echo "$worktree"
