# LSP Traceability Report Contract

OpenQC consumes backend reports that prove scientific docstrings link to LLM
Wiki nodes and that wiki source claims link to raw evidence assets. Backend
repositories own language-specific scanning; OpenQC owns the shared report
contract and family gate.

Schema file: `docs/schemas/lsp-traceability-report.schema.json`

## Required Report

Backends should write one JSON report at one of the paths consumed by
`scripts/check-lsp-family.mjs`, preferably:

```text
reports/docstring-wiki-raw-traceability.json
```

The report must use:

```json
{
  "schemaVersion": "openqc.lsp.traceability.v1",
  "serverId": "cp2k-lsp-enhanced",
  "repository": "newtontech/cp2k-lsp-enhanced",
  "languageId": "cp2k",
  "generatedAt": "2026-06-15T00:00:00.000Z",
  "summary": {
    "docstringsTotal": 2,
    "docstringsLinked": 2,
    "brokenWikiLinks": 0,
    "wikiSourcesWithoutRaw": 0,
    "rawManifestFailures": 0
  },
  "docstrings": [],
  "wikiSources": [],
  "ruleIds": [],
  "sourceUrls": [],
  "rawManifest": {
    "path": "raw/assets/manifest.json",
    "ok": true
  }
}
```

## Rule IDs

Rule IDs must use:

```text
<BACKEND>-<FILE_ROLE>-<CATEGORY>-NNN
```

Examples:

- `CP2K-INPUT-SYNTAX-001`
- `DPGEN-MACHINE-SCHEMA-001`

`BACKEND`, `FILE_ROLE`, and `CATEGORY` are uppercase alphanumeric tokens. Use
underscores inside a token when needed.

## Path Rules

Reports must not require hidden local absolute paths. Use repository-relative
paths for code, wiki, and raw evidence:

- `src/...` for code docstrings and rule sources
- `wiki/...` or `docs/...` for LLM Wiki nodes
- `raw/...` or `raw/assets/...` for raw assets and manifests

Remote source claims must also include a raw evidence path, not only a URL.

## OpenQC Gate

`npm run lsp:check-family -- --strict` validates the report shape when a
backend report is present. Missing reports remain warnings during rollout;
invalid schema, server ID mismatch, unlinked docstrings, broken wiki links,
missing raw evidence links, and raw manifest failures are reported as
traceability gaps with explicit remediation actions.
