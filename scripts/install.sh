#!/usr/bin/env bash
set -euo pipefail

has_npm_script() {
  local script="$1"
  [ -f package.json ] || return 1
  node -e "const p=require('./package.json'); process.exit(p.scripts && p.scripts[process.argv[1]] ? 0 : 1)" "$script"
}

if [ -f package.json ]; then
  if [ -f package-lock.json ]; then
    npm ci
  else
    npm install
  fi
fi

if [ -f Cargo.toml ]; then
  cargo fetch
fi

if [ -f pyproject.toml ] || [ -f setup.py ]; then
  if [ -f pyproject.toml ] && grep -q '^\[tool.poetry\]' pyproject.toml && command -v poetry >/dev/null 2>&1; then
    poetry install --with dev || poetry install
  else
    python -m pip install --upgrade pip
    python -m pip install -e ".[dev]" || python -m pip install -e .
  fi
fi
