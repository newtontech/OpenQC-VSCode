#!/usr/bin/env bash
#
# OpenQC Smoke Test Script
#
# Runs all verification steps for the OpenQC VSCode extension and
# its sibling LSP repositories. This script implements the final
# verification steps from issues #168-#177.
#
# Usage:
#   ./scripts/smoke-test.sh [--code-root DIR] [--run-id ID]
#
# Options:
#   --code-root DIR   Parent directory containing sibling LSP repos (default: ..)
#   --run-id ID       Unique run identifier (default: auto-generated)
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CODE_ROOT="$(cd "$REPO_ROOT/.." && pwd)"
RUN_ID="run-$(date +%Y%m%d-%H%M%S)"
RUN_DIR="/tmp/agent-runs/${RUN_ID}"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --code-root)
      CODE_ROOT="$2"
      shift 2
      ;;
    --run-id)
      RUN_ID="$2"
      RUN_DIR="/tmp/agent-runs/${RUN_ID}"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

mkdir -p "$RUN_DIR"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
RESET='\033[0m'

log_step() {
  echo -e "\n${BLUE}=== $1 ===${RESET}"
}

log_pass() {
  echo -e "  ${GREEN}[PASS]${RESET} $1"
}

log_fail() {
  echo -e "  ${RED}[FAIL]${RESET} $1"
}

log_warn() {
  echo -e "  ${YELLOW}[WARN]${RESET} $1"
}

log_info() {
  echo -e "  ${BLUE}[INFO]${RESET} $1"
}

TOTAL_FAILURES=0

# ---------------------------------------------------------------------------
# Step 1: Compile TypeScript
# ---------------------------------------------------------------------------

log_step "Compiling TypeScript"

if npm run compile --prefix "$REPO_ROOT" > "$RUN_DIR/compile.log" 2>&1; then
  log_pass "TypeScript compilation succeeded"
else
  log_fail "TypeScript compilation failed (see $RUN_DIR/compile.log)"
  TOTAL_FAILURES=$((TOTAL_FAILURES + 1))
fi

# ---------------------------------------------------------------------------
# Step 2: Run unit and integration tests
# ---------------------------------------------------------------------------

log_step "Running unit and integration tests"

if npm test --prefix "$REPO_ROOT" > "$RUN_DIR/tests.log" 2>&1; then
  log_pass "All tests passed"
else
  log_fail "Some tests failed (see $RUN_DIR/tests.log)"
  TOTAL_FAILURES=$((TOTAL_FAILURES + 1))
fi

# ---------------------------------------------------------------------------
# Step 3: Run lint and typecheck
# ---------------------------------------------------------------------------

log_step "Running lint and typecheck"

if npm run lint --prefix "$REPO_ROOT" > "$RUN_DIR/lint.log" 2>&1; then
  log_pass "Lint passed"
else
  log_warn "Lint had issues (see $RUN_DIR/lint.log)"
fi

if npm run typecheck --prefix "$REPO_ROOT" > "$RUN_DIR/typecheck.log" 2>&1; then
  log_pass "Typecheck passed"
else
  log_fail "Typecheck failed (see $RUN_DIR/typecheck.log)"
  TOTAL_FAILURES=$((TOTAL_FAILURES + 1))
fi

# ---------------------------------------------------------------------------
# Step 4: Check sibling LSP repos exist
# ---------------------------------------------------------------------------

log_step "Checking sibling LSP repository checkouts"

CHECKOUTS_FOUND=0
CHECKOUTS_MISSING=0

while IFS= read -r repo_name; do
  if [ -d "$CODE_ROOT/$repo_name/.git" ]; then
    log_pass "$repo_name checkout exists"
    CHECKOUTS_FOUND=$((CHECKOUTS_FOUND + 1))
  else
    log_warn "$repo_name not found at $CODE_ROOT/$repo_name"
    CHECKOUTS_MISSING=$((CHECKOUTS_MISSING + 1))
  fi
done <<EOF
abacus-lsp
abinit-lsp
cif-lsp
cp2k-lsp-enhanced
VASP-LSP
gaussian-lsp
orca-lsp
qe-lsp
gamess-lsp
nwchem-lsp
gpumd-lsp
gromacs-lsp
lammps-lsp
mlip-lsp
pyatb-lsp
pyscf-lsp
EOF

log_info "Found $CHECKOUTS_FOUND / $((CHECKOUTS_FOUND + CHECKOUTS_MISSING)) sibling repos"

# ---------------------------------------------------------------------------
# Step 5: Check repo cleanliness
# ---------------------------------------------------------------------------

log_step "Checking sibling repo cleanliness"

for repo_dir in "$CODE_ROOT"/*/; do
  if [ -d "$repo_dir/.git" ]; then
    repo_name="$(basename "$repo_dir")"
    status="$(git -C "$repo_dir" status --short 2>/dev/null || echo "git-failed")"
    if [ -z "$status" ]; then
      log_pass "$repo_name is clean"
    else
      log_warn "$repo_name has uncommitted changes"
      echo "$status" | head -5
    fi
  fi
done

# ---------------------------------------------------------------------------
# Step 6: Verify VSIX can be built
# ---------------------------------------------------------------------------

log_step "Verifying VSIX package dependencies"

if [ -f "$REPO_ROOT/package.json" ]; then
  HAS_LC=$(node -e "
    const pkg = require('$REPO_ROOT/package.json');
    const deps = Object.keys(pkg.dependencies || {});
    const devDeps = Object.keys(pkg.devDependencies || {});
    const all = [...deps, ...devDeps];
    console.log(all.includes('vscode-languageclient') ? 'yes' : 'no');
  ")
  if [ "$HAS_LC" = "yes" ]; then
    log_pass "vscode-languageclient found in dependencies"
  else
    log_fail "vscode-languageclient missing from dependencies"
    TOTAL_FAILURES=$((TOTAL_FAILURES + 1))
  fi
fi

# ---------------------------------------------------------------------------
# Step 7: Run compatibility report
# ---------------------------------------------------------------------------

log_step "Generating compatibility report"

node -e "
  const path = require('path');
  const { generateCompatibilityReport, formatReportAsMarkdown } = require('$REPO_ROOT/out/smoke/compatibilityReport');
  const report = generateCompatibilityReport('$REPO_ROOT');
  const md = formatReportAsMarkdown(report);
  require('fs').writeFileSync('$RUN_DIR/compatibility-report.md', md);
  console.log('Total: ' + report.totalServers + ', Passed: ' + report.passedServers + ', Failed: ' + report.failedServers);
  process.exit(report.failedServers > 0 ? 1 : 0);
" > "$RUN_DIR/compat.log" 2>&1 && log_pass "Compatibility report generated" || log_warn "Some compatibility checks failed (see $RUN_DIR/compat.log)"

# ---------------------------------------------------------------------------
# Step 8: Final summary
# ---------------------------------------------------------------------------

log_step "Final Summary"

if [ $TOTAL_FAILURES -eq 0 ]; then
  log_pass "All smoke tests passed!"
  echo -e "\n${GREEN}Campaign report saved to: $RUN_DIR${RESET}"
  exit 0
else
  log_fail "$TOTAL_FAILURES step(s) failed"
  echo -e "\n${RED}Campaign report saved to: $RUN_DIR${RESET}"
  exit 1
fi
