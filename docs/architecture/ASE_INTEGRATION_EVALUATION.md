# ASE Integration Architecture Evaluation

Status: proposed simplification for issue #17

## Context

`PLAN.md` currently describes ASE as a universal intermediate layer for Phase 3
cross-code migration. The implemented code path is narrower:

- `src/converters/FormatConverter.ts` shells out to `python/format_converter.py`.
- `python/format_converter.py` imports `dpdata` on demand and does not import ASE.
- `python/requirements.txt` previously installed ASE even though the default
  converter does not require it.
- The architecture overview mentioned a future `ASEAdapter`, but no
  `src/ase/` or adapter implementation exists in this checkout.

That means the current risk is not an overbuilt ASE runtime yet. The risk is
architectural drift: docs and dependencies imply that ASE is a core layer before
there is evidence that users need calculator orchestration or full workflow
migration inside the VS Code extension.

## Current Usage Scenarios

### Core scenarios

1. Open a supported input file and get syntax, parsing, validation, and
   visualization in VS Code.
2. Convert structure files between practical editor-facing formats such as
   POSCAR, XYZ, PDB, CIF, Gaussian, and ORCA.
3. Batch-convert files from the command palette or the Python CLI.

These scenarios are already served by TypeScript parsers, visualization code,
and the dpdata subprocess converter. They do not require ASE as a mandatory
dependency.

### Advanced scenarios

1. Convert structures through an `ase.Atoms` representation when dpdata cannot
   preserve enough geometry, cell, constraint, or trajectory information.
2. Generate inputs for multiple simulation engines from a common structure.
3. Use ASE calculators to launch jobs and read results back into OpenQC.
4. Validate migration quality with round-trip tests and reference fixtures.

These scenarios are useful, but they are not essential to the extension's
current editing and visualization loop. Calculator execution also overlaps with
job lifecycle and remote execution concerns, so it should not be part of a KISS
format-conversion baseline.

## Alternatives

| Alternative | What it means | Advantages | Costs and risks |
|-------------|---------------|------------|-----------------|
| Keep ASE as mandatory Python backend | Install ASE with the default converter and route conversions through it | Maximum future flexibility, broad chemistry ecosystem | Higher install cost, dependency conflicts, dual-language maintenance before clear demand |
| Optional ASE plugin path | Keep dpdata as the default converter and enable ASE only for advanced commands or explicitly configured workflows | Small default install, preserves future escape hatch, lets tests prove value before expansion | Requires capability detection and clear UX for optional features |
| Pure TypeScript baseline | Implement key structure conversions directly in the extension | Simplest deployment, no Python for core flows, easiest packaging | Limited format coverage, duplicated chemistry logic, higher risk for edge cases |
| Remote conversion service | Send conversion/migration requests to a hosted backend | No local Python dependency, centrally maintained dependencies | Network/privacy concerns, operating cost, credential and availability burden |
| WebAssembly chemistry backend | Bundle compiled conversion logic into the extension | Local execution without Python, predictable packaging | Significant build complexity, limited ASE reuse, unclear library maturity for target formats |

## Recommendation

Use a two-tier architecture:

1. **Default tier: TypeScript + dpdata subprocess.** Keep the current converter
   as the only required conversion backend. It is already implemented, easy to
   reason about, and supports the current user-facing conversion commands.
2. **Optional tier: ASE capability module.** Add ASE later only for commands
   that explicitly need `ase.Atoms`, structure manipulation, or migration
   validation. Gate it behind backend capability checks and optional dependency
   installation.
3. **Defer calculator execution.** Do not put ASE calculators into the
   extension baseline until job execution/lifecycle boundaries are separately
   designed and tested.

This keeps the current product simple while preserving a path for scientifically
meaningful ASE features.

## Acceptance Criteria for Future ASE Work

Before introducing an ASE adapter or `src/ase/` module, require:

- A named command or workflow that cannot be implemented well with dpdata or a
  small TypeScript parser.
- Fixture-based round-trip tests for each supported source and target format.
- Capability detection that reports whether Python, dpdata, and ASE are present
  separately.
- Documentation that labels ASE as optional and lists the exact features it
  enables.
- Telemetry or user research showing demand for migration features beyond basic
  format conversion.

## Concrete Simplification

- Keep ASE out of the default `python/requirements.txt`.
- Document ASE as optional future functionality, not an implemented core
  adapter.
- Treat Phase 3 ASE work as an evaluation gate: first prove one high-value
  migration workflow, then add the smallest adapter surface needed for it.

