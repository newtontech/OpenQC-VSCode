# Agent Memory

## 2026-07-16: Viewer and calculation workflow closure

- Issue #238 is implemented on `agent/issue-238-viewer-workflows`; the branch also merges the completed issue #236 release pipeline through `d14151c`.
- The production 3Dmol viewer now supports bounded editing, selective-dynamics constraints, undo/redo, supercells, trajectories, PNG/native export, and limited source writeback. Structures above 10,000 atoms use a deterministic, explicit read-only 10,000-atom LOD.
- Viewer loading now consumes the compute worker, structure LRU cache, and incremental parser. Calculator, Jobs, analyzer, output-result, cancellation, restart, trajectory, and export paths are wired to real handlers.
- 3Dmol, NGL, and Plotly runtime bundles are vendored under `media/vendor`; they remain development dependencies and must not be re-included from `node_modules` in `.vscodeignore`.
- The release verifier must require the `media/vendor` bundle paths. Packaged visual smoke extracts `vsix/openqc-0.0.1.vsix` and uses those exact runtime files.
- Final local gates use `make format`, `make lint`, `make typecheck`, `make test`, and `make check`, plus source and packaged visual smoke. The only lint output is the existing non-error curly warning in `src/agent/OpenQCToolRegistry.ts`.
- The installed Node 22 binary was killed by the host before startup, so local commands ran on Node 26 with an engine warning; CI and release workflows enforce Node 22.
- Ubuntu CI must install Playwright Chromium after `npm ci`. On software WebGL, `waitForFunction` can time out after the viewer marker is already `ready`; accept only a directly re-read `ready` marker, then retain screenshot pixel validation so genuine boot/error states still fail.

## 2026-07-16 - Issue #237

- The production AI bridge must launch `openqc.ai.client` with `python -m`; eager imports from
  `openqc.ai.__init__` trigger a `runpy` warning, so client exports are lazy.
- OpenAI uses the Responses API and receives its credential only from VS Code SecretStorage.
  Fake local HTTP servers cover provider behavior; tests must never send live provider requests.
- The dependency-free MCP server declares five tools. Keep the tool schemas, handler routing,
  JSON-RPC validation, request IDs, and normalized optional-dependency errors in one module.
- The PR contract recognizes test evidence only when the body contains `Tests`, `Test Plan`, or
  `TDD Evidence`; a generic `Validation` heading does not satisfy the check.
- After issue #236 retired the duplicate `core/` tree, the combined AI/MCP and release baseline
  passes the maintained format, lint, typecheck, Jest, Python bridge, and release-candidate gates.
