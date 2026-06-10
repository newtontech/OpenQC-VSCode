/**
 * LSP Registry Types
 *
 * Defines the shape of a bundled LSP server registry entry and the
 * stability level for each server.
 *
 * Also defines types for domain language description APIs consumed
 * from standalone LSP servers.
 *
 * @module lsp/types
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/145
 */

/** Stability classification for an LSP server entry. */
export type LSPStability = 'stable' | 'experimental';

/** Local sibling-repository launch strategy for a bundled LSP server. */
export type LocalLspLaunch =
  | {
      readonly kind: 'pythonFunction';
      readonly repoName: string;
      readonly importPath: string;
      readonly functionName: string;
      readonly sourcePath?: string;
    }
  | {
      readonly kind: 'nodeScript';
      readonly repoName: string;
      readonly scriptPath: string;
    }
  | {
      readonly kind: 'cargoBinary';
      readonly repoName: string;
      readonly binaryName: string;
      readonly cargoBin?: string;
    };

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

  /** Default arguments passed to the executable. Defaults to ["--stdio"]. */
  readonly args?: readonly string[];

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

  /** How to launch this server directly from a sibling local repository. */
  readonly localLaunch?: LocalLspLaunch;

  /**
   * Default branch to clone or reference, when it differs from "main".
   * Omit this field when the default branch is the repo's default.
   */
  readonly defaultBranch?: string;
}

// ---------------------------------------------------------------------------
// Domain language description capability types
// ---------------------------------------------------------------------------

/**
 * Flags indicating which domain-specific capabilities a standalone LSP
 * server advertises during initialization.
 *
 * These are discovered from custom extensions in the server's
 * `ServerCapabilities` during the `initialize` handshake.
 */
export interface DomainLSPCapabilities {
  /** Server can provide a high-level language description. */
  readonly languageDescription: boolean;
  /** Server can look up sections and keywords. */
  readonly sectionKeywordLookup: boolean;
  /** Server can provide minimal usage examples. */
  readonly minimalExamples: boolean;
  /** Server can suggest the next token given a context. */
  readonly nextTokenGuidance: boolean;
}

/** Parameters for the `openqc/languageDescription` request. */
export interface LanguageDescriptionParams {
  /** The language ID to describe (e.g. "gaussian"). */
  readonly languageId: string;
}

/** Response from the `openqc/languageDescription` request. */
export interface LanguageDescriptionResponse {
  /** Human-readable name of the language. */
  readonly name: string;
  /** Brief summary of the language's purpose. */
  readonly description: string;
  /** Top-level sections or categories in the language. */
  readonly sections: readonly string[];
}

/** Parameters for the `openqc/sectionKeywordLookup` request. */
export interface SectionKeywordParams {
  /** The language ID to search within. */
  readonly languageId: string;
  /** Section name or keyword prefix to look up. */
  readonly query: string;
}

/** A single keyword entry returned by section/keyword lookup. */
export interface SectionKeywordEntry {
  /** Keyword or section name. */
  readonly name: string;
  /** Brief documentation string. */
  readonly documentation: string;
  /** Section this keyword belongs to, if applicable. */
  readonly section?: string;
}

/** Response from the `openqc/sectionKeywordLookup` request. */
export interface SectionKeywordResponse {
  /** Matching entries. */
  readonly entries: readonly SectionKeywordEntry[];
}

/** Parameters for the `openqc/minimalExample` request. */
export interface MinimalExampleParams {
  /** The language ID for the example. */
  readonly languageId: string;
  /** Optional keyword to focus the example around. */
  readonly keyword?: string;
}

/** Response from the `openqc/minimalExample` request. */
export interface MinimalExampleResponse {
  /** Example input file content. */
  readonly example: string;
  /** Short description of what the example demonstrates. */
  readonly description: string;
}

/** Parameters for the `openqc/nextTokenGuidance` request. */
export interface NextTokenGuidanceParams {
  /** The language ID. */
  readonly languageId: string;
  /** Document URI for context. */
  readonly textDocument: { readonly uri: string };
  /** Cursor position (line, character) within the document. */
  readonly position: { readonly line: number; readonly character: number };
}

/** A single token suggestion with its documentation. */
export interface NextTokenCandidate {
  /** Suggested token text. */
  readonly token: string;
  /** Brief explanation of what this token does. */
  readonly documentation?: string;
}

/** Response from the `openqc/nextTokenGuidance` request. */
export interface NextTokenGuidanceResponse {
  /** Ranked candidate tokens. */
  readonly candidates: readonly NextTokenCandidate[];
}

/**
 * Custom server capabilities extension shape that standalone LSP servers
 * may include in their `ServerCapabilities` during initialization.
 *
 * Negotiated via the `capabilities` field in the `InitializeResult`.
 */
export interface OpenQCServerCapabilities {
  /** Domain-specific capability flags, if the server supports them. */
  readonly openqc?: {
    readonly domainCapabilities?: {
      readonly languageDescription?: boolean;
      readonly sectionKeywordLookup?: boolean;
      readonly minimalExamples?: boolean;
      readonly nextTokenGuidance?: boolean;
    };
  };
}
