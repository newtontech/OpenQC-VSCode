#!/usr/bin/env bash
set -euo pipefail

ran=0

has_npm_script() {
  local script="$1"
  [ -f package.json ] || return 1
  node -e "const p=require('./package.json'); process.exit(p.scripts && p.scripts[process.argv[1]] ? 0 : 1)" "$script"
}

if has_npm_script format:write; then
  npm run format:write
  ran=1
elif has_npm_script format; then
  npm run format
  ran=1
fi

if [ -f Cargo.toml ]; then
  cargo fmt
  ran=1
fi

if [ -f pyproject.toml ] || [ -f setup.py ]; then
  py_targets="$(python_format_targets)"
  if python -m black --version >/dev/null 2>&1; then
    python -m black $py_targets
    ran=1
  elif python -m ruff --version >/dev/null 2>&1; then
    python -m ruff format $py_targets
    ran=1
  fi
fi

if [ "$ran" -eq 0 ]; then
  echo "No formatter configured; skipping."
fi
