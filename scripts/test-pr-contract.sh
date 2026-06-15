#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CHECK="$ROOT_DIR/scripts/check-pr-contract.sh"

expect_pass() {
  local name="$1"
  shift
  if "$@"; then
    echo "PASS expected: $name"
  else
    echo "Expected pass but failed: $name" >&2
    exit 1
  fi
}

expect_fail() {
  local name="$1"
  shift
  if "$@"; then
    echo "Expected failure but passed: $name" >&2
    exit 1
  else
    echo "FAIL expected: $name"
  fi
}

expect_fail "human PR without issue link" env \
  PR_AUTHOR="octocat" \
  PR_HEAD="feature/no-issue" \
  PR_BODY="Tests: npm test" \
  CHANGED_FILES="src/extension.ts" \
  "$CHECK"

expect_pass "dependabot workflow-only update" env \
  PR_AUTHOR="dependabot[bot]" \
  PR_HEAD="dependabot/github_actions/actions/checkout-6" \
  PR_BODY="Bumps actions/checkout from 4 to 6." \
  CHANGED_FILES=".github/workflows/pr-contract.yml" \
  "$CHECK"

expect_fail "dependabot source change still needs tests or justification" env \
  PR_AUTHOR="dependabot[bot]" \
  PR_HEAD="dependabot/npm/typescript-6" \
  PR_BODY="Bumps typescript." \
  CHANGED_FILES="src/extension.ts" \
  "$CHECK"

expect_pass "agent issue branch with tests" env \
  PR_AUTHOR="newtontech" \
  PR_HEAD="agent/issue-231-dependabot-contract" \
  PR_BODY=$'Fixes #231\n\nTests: scripts/test-pr-contract.sh' \
  CHANGED_FILES=$'scripts/check-pr-contract.sh\nscripts/test-pr-contract.sh' \
  "$CHECK"
