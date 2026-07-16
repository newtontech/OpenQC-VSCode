# Agent Memory

## 2026-07-16 - Issue #237

- The production AI bridge must launch `openqc.ai.client` with `python -m`; eager imports from
  `openqc.ai.__init__` trigger a `runpy` warning, so client exports are lazy.
- OpenAI uses the Responses API and receives its credential only from VS Code SecretStorage.
  Fake local HTTP servers cover provider behavior; tests must never send live provider requests.
- The dependency-free MCP server declares five tools. Keep the tool schemas, handler routing,
  JSON-RPC validation, request IDs, and normalized optional-dependency errors in one module.
- `make lint` is currently blocked by 544 pre-existing Ruff errors under `core/`. `make test`
  reaches a missing `pytest` executable and a nonexistent `core/tests/unit` path after the Jest
  unit suite passes; use the maintained root Python tests plus npm/Jest gates for this slice.
