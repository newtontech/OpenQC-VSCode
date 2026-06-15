# LSP Family Maturity Audit - 2026-06-15

This page is the current wiki-style summary for the OpenQC bundled LSP family. It captures what the latest release gates prove, what rules have generalized across the fleet, and what remains before the scientific LSPs reach top-tier Python/C++ LSP expectations such as Pyright/Pylance and clangd.

## Current Verdict

OpenQC can currently block a VSIX release on the latest managed LSP checkouts without drift or hard failures.

Verified on `2026-06-15` from the OpenQC masterline after PR `#190`:

| Gate | Result | Evidence |
|------|--------|----------|
| Latest checkout drift | Pass | `npm run lsp:check-latest -- --fail-on-drift` reported 17/17 registry entries latest, clean, and agent help passing. |
| Family strict release gate | Pass with warnings | `npm run lsp:check-family -- --strict` reported 17/17 passing, 0 blocking gaps, 38 warning gaps. |
| Bohrium skill probe | Pass | `lsp_skill_probe.py --root /Users/yhm/Desktop/code/.worktrees-lsp-latest --json` found 16/16 registered skill backends with available agent CLI help. |
| OpenQC CI for gate fix | Pass | PR `#190` passed TypeScript tests, extension tests, build, pre-commit, security, contract, and governance checks before merge. |

The remaining work is no longer basic wiring. It is maturity work: official-document provenance, full version-aware grammar/rule coverage, closed-loop repair/output diagnostics, and parity with the editor ergonomics expected from mature language servers.

## Fleet Rules That Have Generalized

The same release model now works across all bundled LSPs:

1. **Managed latest checkouts are authoritative.** Release gates must prefer `.worktrees-lsp-latest/<repo>` over dirty or stale sibling repositories.
2. **Every LSP needs an agent-facing CLI.** At minimum, `*-lsp-tool --help` and `check <file-or-dir> --format json` must work without editor activation.
3. **Every LSP needs a capability manifest.** `lsp-capabilities.json` should identify the registry entry, repository, schema or capabilities version, supported operations, file patterns, provenance, and OpenQC integration metadata.
4. **Diagnostic output should use a stable envelope.** The target contract is DiagnosticEnvelope/v1: stable severity, code, range, source, rule id, fixability, and machine-readable context.
5. **Fixtures are release evidence.** Each repo needs valid and invalid examples that exercise parser coverage, diagnostic coverage, and OpenQC smoke behavior.
6. **Provenance must be queryable.** Release quality requires source docs, version mapping, changelog or version files, and commit/tag traceability.
7. **Closed-loop features must be tested.** Repair previews, output/log diagnostics, and OpenQC smoke checks should be backed by fixtures and CLI evidence, not just declared in docs.
8. **OpenQC is the product gate.** A standalone LSP can be healthy but still not release-ready unless OpenQC can discover it, route files to it, and validate its agent CLI from the latest checkout.

## What "Top Python/C++ LSP" Means Here

The target is not only syntax highlighting or keyword completion. The scientific LSP family should eventually provide:

| Capability | Required behavior |
|------------|-------------------|
| Complete grammar | Parse every documented input construct for every supported software family. |
| Version awareness | Know which keywords, cards, namelists, blocks, units, and defaults exist for each software version. |
| Semantic diagnostics | Catch invalid combinations, missing required fields, wrong units, deprecated keywords, invalid references, and cross-file inconsistencies. |
| Fast incremental feedback | Return diagnostics quickly enough for edit-time use, with predictable CLI behavior for agents. |
| Completion and hover | Offer context-aware completions and explanations backed by versioned source documentation. |
| Navigation and symbols | Expose sections, blocks, groups, variables, atoms, species, files, and references as document/workspace symbols. |
| Formatting and code actions | Provide stable formatting and safe fix previews where a correction can be made mechanically. |
| Multi-file context | Validate related inputs together: structure files, topology files, potentials, pseudopotentials, logs, and generated outputs. |
| Output/log diagnostics | Parse runtime logs and map failures back to likely input causes when possible. |
| Provenance and wiki | Preserve doc source, version, extraction date, and rule derivation for every rule exposed to agents. |
| Regression harness | Keep valid/invalid fixtures, output fixtures, and OpenQC smoke tests in CI for every LSP. |

## Current Warning Matrix

The OpenQC family gate has no blocking errors, but it still reports these maturity warnings.

| Repo | Remaining warnings |
|------|--------------------|
| `abacus-lsp` | Missing canonical capabilities section; no git tags; no `CHANGELOG.md` or `VERSION`; closed-loop support still planned. |
| `abinit-lsp` | No git tags; no `CHANGELOG.md` or `VERSION`; closed-loop support still planned. |
| `cif-lsp` | Closed-loop support still planned. |
| `cp2k-lsp-enhanced` | No current family-gate warnings. |
| `VASP-LSP` | No git tags. |
| `gaussian-lsp` | Missing canonical capabilities section; no git tags; closed-loop support still planned. |
| `orca-lsp` | No git tags; closed-loop support still planned. |
| `qe-lsp` | No git tags. |
| `gamess-lsp` | No git tags; closed-loop support still planned. |
| `nwchem-lsp` | No invalid fixture files detected by the family gate; closed-loop support still planned. |
| `gpumd-lsp` | No git tags; no `CHANGELOG.md` or `VERSION`; closed-loop support still planned. |
| `gromacs-lsp` | No `CHANGELOG.md` or `VERSION`; closed-loop support still planned. |
| `lammps-lsp` | No current family-gate warnings. |
| `mlip-lsp` | No git tags; no `CHANGELOG.md` or `VERSION`; closed-loop support still planned. |
| `pyatb-lsp` | No git tags; no `CHANGELOG.md` or `VERSION`; closed-loop support still planned. |
| `pyscf-lsp` | No git tags; no `CHANGELOG.md` or `VERSION`; closed-loop support still planned. |
| `dpgen-lsp` | No git tags; no `CHANGELOG.md` or `VERSION`; no valid/invalid fixtures detected by the family gate; closed-loop support still planned. |

## Open Issue Summary

The first implementation waves closed the immediate smoke/fixture/runtime PRs and left no open PRs. The remaining open issues are intentionally larger capability tracks.

| Area | Open work |
|------|-----------|
| OpenQC product gate | `newtontech/OpenQC-VSCode#187` tracks editor-grade parity across the family. |
| Version-aware grammar/rules | Still open for most active LSP repos. These issues cover exhaustive parser/schema/rule coverage across software versions. |
| Official-docs to wiki/provenance | Still open for most LSP repos. These issues cover manual ingestion, source provenance, rule extraction, and version mapping. |
| Closed-loop evidence | Still open where repair previews, output diagnostics, or OpenQC smoke evidence are not complete. |

The important distinction is that the release gate is now green, but the "perfect LSP" goal is not complete. The gate proves integration readiness; it does not prove full grammar coverage for every historical and current version of every scientific package.

## Next Iteration Order

1. **Remove low-risk family warnings.**
   Add missing `CHANGELOG.md` or `VERSION` files, create initial release tags where appropriate, normalize `capabilities` sections, and add fixture names that the gate can detect.

2. **Graduate closed-loop status from planned to implemented.**
   For each repo still marked planned, add CLI-backed repair previews or output/log diagnostics, then update OpenQC `diagnosticReadiness` only after fixtures and local evidence exist.

3. **Build the official-docs to wiki pipeline.**
   Each repo needs a repeatable extraction path from official manuals to versioned rule tables, with source URL, software version, extraction commit/date, and generated rule ids.

4. **Complete version-aware grammar coverage.**
   Convert wiki rule tables into parser/schema diagnostics. The coverage target is every keyword/block/card/section in the supported version range, including deprecated and version-gated constructs.

5. **Add cross-file and runtime diagnostics.**
   Prioritize failures that cost agents the most iterations: missing potential files, inconsistent atom/species names, topology/index mismatches, unsupported methods, unit mistakes, and common runtime-log failures.

6. **Expand editor parity beyond diagnostics.**
   Add context-aware completion, hover, symbols, formatting, and code actions with the same version provenance used by diagnostics.

7. **Keep OpenQC as the final acceptance gate.**
   A repo-level PR should not be considered done until `.worktrees-lsp-latest` is refreshed and both OpenQC gates pass.

## Acceptance Checklist

A repo can be treated as mature only when all of the following are true:

- `lsp-capabilities.json` has canonical identity, version, repository, capabilities, file patterns, and OpenQC metadata.
- `*-lsp-tool --help` and `*-lsp-tool check ... --format json` run from a clean latest checkout.
- Valid, invalid, and runtime/log fixtures are present and exercised in tests.
- DiagnosticEnvelope/v1 includes stable codes, ranges, severities, rule ids, and fix metadata.
- Official-source wiki/provenance exists for every implemented rule.
- Version gating is explicit in parser, completion, hover, and diagnostics.
- OpenQC `lsp:check-latest -- --fail-on-drift` passes.
- OpenQC `lsp:check-family -- --strict` passes with zero blocking gaps.
- Remaining warnings are either zero or documented as non-blocking product decisions.

