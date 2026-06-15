#!/usr/bin/env bash
set -euo pipefail

PR_BODY="${PR_BODY:-}"
PR_HEAD="${PR_HEAD:-}"
PR_AUTHOR="${PR_AUTHOR:-}"
BASE_REF="${BASE_REF:-}"

changed_files() {
  if [[ -n "${CHANGED_FILES:-}" ]]; then
    printf '%s\n' "$CHANGED_FILES"
    return
  fi

  if [[ -z "$BASE_REF" ]]; then
    echo "BASE_REF is required when CHANGED_FILES is not provided" >&2
    return 2
  fi

  git fetch origin "$BASE_REF" --depth=1 >/dev/null
  git diff --name-only "origin/${BASE_REF}"...HEAD
}

is_dependabot_pr() {
  [[ "$PR_AUTHOR" == "dependabot[bot]" || "$PR_HEAD" == dependabot/* ]]
}

is_dependency_maintenance_file() {
  case "$1" in
    .github/dependabot.yml|.github/dependabot.yaml) return 0 ;;
    .github/workflows/*.yml|.github/workflows/*.yaml) return 0 ;;
    package.json|package-lock.json|npm-shrinkwrap.json|pnpm-lock.yaml|yarn.lock) return 0 ;;
    pyproject.toml|poetry.lock|uv.lock|requirements*.txt|setup.cfg|setup.py) return 0 ;;
    Cargo.toml|Cargo.lock|go.mod|go.sum) return 0 ;;
    Gemfile|Gemfile.lock|composer.json|composer.lock) return 0 ;;
    Dockerfile|docker-compose.yml|docker-compose.yaml) return 0 ;;
    *) return 1 ;;
  esac
}

all_changes_are_dependency_maintenance() {
  local path
  local saw_file=0
  while IFS= read -r path; do
    [[ -z "$path" ]] && continue
    saw_file=1
    if ! is_dependency_maintenance_file "$path"; then
      echo "Non-dependency-maintenance change: $path"
      return 1
    fi
  done
  [[ "$saw_file" -eq 1 ]]
}

changed="$(changed_files)"
printf '%s\n' "$changed"

dependabot_dependency_only=0
if is_dependabot_pr && all_changes_are_dependency_maintenance <<<"$changed"; then
  dependabot_dependency_only=1
  echo "Dependabot dependency-maintenance PR detected; body issue/test prose contract is exempt."
fi

if [[ "$dependabot_dependency_only" -ne 1 ]]; then
  printf '%s\n' "$PR_BODY" | grep -E "(Fixes|Closes|Resolves) #[0-9]+" >/dev/null \
    || { echo "PR body must include Fixes/Closes/Resolves #<issue-number>"; exit 1; }

  printf '%s\n' "$PR_BODY" | grep -E "(TDD Evidence|Tests|Test Plan|No-test justification)" >/dev/null \
    || { echo "PR body must include test evidence or a no-test justification"; exit 1; }
fi

if printf '%s\n' "$PR_HEAD" | grep -E '^agent/issue-[0-9]+-' >/dev/null; then
  echo "Agent branch naming is valid: $PR_HEAD"
else
  echo "Non-agent or legacy branch name: $PR_HEAD"
fi

code_changed="$(printf '%s\n' "$changed" | grep -E '\.(py|ts|tsx|js|jsx|rs)$' || true)"
if [[ -z "$code_changed" ]]; then
  echo "No source-code changes detected."
  exit 0
fi

test_changed="$(printf '%s\n' "$changed" | grep -E '(^tests/|/tests/|test_|_test\.|\.spec\.|\.test\.)' || true)"
if [[ -n "$test_changed" ]]; then
  echo "Test changes detected."
  exit 0
fi

printf '%s\n' "$PR_BODY" | grep -Ei "no-test justification|docs-only|config-only|refactor-only" >/dev/null \
  || { echo "Code PRs require test changes or explicit no-test justification."; exit 1; }
