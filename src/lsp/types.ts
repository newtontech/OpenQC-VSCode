/**
 * LSP Registry Types
 *
 * Defines the shape of a bundled LSP server registry entry and the
 * stability level for each server.
 *
 * @module lsp/types
 */

/** Stability classification for an LSP server entry. */
export type LSPStability = 'stable' | 'experimental';

/**
 * A single entry in the bundled LSP server registry.
 *
 * This is the version-controlled single source of truth for known LSP servers.
 * It is intentionally a subset of the full `LSPServerDefinition` used at
 * runtime — fields that only make sense after network discovery (e.g.
 * `lastUpdated`, `description`) are omitted.
 */
export interface LSPServerRegistryEntry {
  /** Unique registry key, e.g. "gaussian-lsp". */
  readonly id: string;

  /** Human-readable software name, e.g. "Gaussian". */
  readonly name: string;

  /** Full GitHub repository path, e.g. "newtontech/gaussian-lsp". */
  readonly repository: string;

  /** Executable name on PATH, e.g. "gaussian-lsp". */
  readonly executable: string;

  /** VS Code language ID, e.g. "gaussian". */
  readonly languageId: string;

  /** File extensions this LSP handles (without leading dot). */
  readonly fileExtensions: readonly string[];

  /** Exact file names this LSP handles (e.g. "INCAR", "POSCAR"). */
  readonly fileNames: readonly string[];

  /** Whether the LSP is enabled by default. */
  readonly enabled: boolean;

  /** Repository URL for reference and linking. */
  readonly repositoryUrl: string;

  /** Stability classification. */
  readonly stability: LSPStability;

  /**
   * Default branch to clone or reference, when it differs from "main".
   * Omit this field when the default branch is the repo's default.
   */
  readonly defaultBranch?: string;
}
