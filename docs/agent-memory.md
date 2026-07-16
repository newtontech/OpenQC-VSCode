# Agent Memory

## 2026-07-16: Viewer and calculation workflow closure

- Issue #238 is implemented on `agent/issue-238-viewer-workflows`; the branch also merges the completed issue #236 release pipeline through `d14151c`.
- The production 3Dmol viewer now supports bounded editing, selective-dynamics constraints, undo/redo, supercells, trajectories, PNG/native export, and limited source writeback. Structures above 10,000 atoms use a deterministic, explicit read-only 10,000-atom LOD.
- Viewer loading now consumes the compute worker, structure LRU cache, and incremental parser. Calculator, Jobs, analyzer, output-result, cancellation, restart, trajectory, and export paths are wired to real handlers.
- 3Dmol, NGL, and Plotly runtime bundles are vendored under `media/vendor`; they remain development dependencies and must not be re-included from `node_modules` in `.vscodeignore`.
- The release verifier must require the `media/vendor` bundle paths. Packaged visual smoke extracts `vsix/openqc-0.0.1.vsix` and uses those exact runtime files.
- Final local gates use `make format`, `make lint`, `make typecheck`, `make test`, and `make check`, plus source and packaged visual smoke. The only lint output is the existing non-error curly warning in `src/agent/OpenQCToolRegistry.ts`.
- The installed Node 22 binary was killed by the host before startup, so local commands ran on Node 26 with an engine warning; CI and release workflows enforce Node 22.
