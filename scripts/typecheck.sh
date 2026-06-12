#!/usr/bin/env bash
set -euo pipefail

ran=0

has_npm_script() {
  local script="$1"
  [ -f package.json ] || return 1
  node -e "const p=require('./package.json'); process.exit(p.scripts && p.scripts[process.argv[1]] ? 0 : 1)" "$script"
}

python_typecheck_targets() {
  local targets=""
  [ -d src ] && targets="$targets src"
  for d in *_lsp cp2k_input_tools mdparser gromacs_lsp; do
    [ -d "$d" ] && targets="$targets $d"
  done
  if [ -z "${targets# }" ]; then
    echo "."
  else
    echo "$targets"
  fi
}

if has_npm_script typecheck; then
  npm run typecheck
  ran=1
elif [ -f tsconfig.json ] && [ -d node_modules ]; then
  npx tsc --noEmit
  ran=1
fi

if [ -f Cargo.toml ]; then
  cargo check --all-targets --all-features
  ran=1
fi

if [ -f pyproject.toml ] || [ -f setup.py ] || [ -f mypy.ini ]; then
  if python -m mypy --version >/dev/null 2>&1; then
    py_targets="$(python_typecheck_targets)"
    python -m mypy $py_targets
    ran=1
  fi
fi

if [ "$ran" -eq 0 ]; then
  echo "No typechecker configured; skipping."
fi
