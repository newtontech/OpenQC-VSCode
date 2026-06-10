/**
 * LSP Registry Types
 *
 * Defines the shape of a bundled LSP server registry entry and the
 * stability level for each server.
 *
 * Also defines types for the DSL authoring context aggregation API,
 * which compacts standalone LSP features into a single agent-facing
 * context bundle for OpenQC coding-agent workflows.
 *
 * @module lsp/types
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/146
 */

/** Stability classification for an LSP server entry. */
export type LSPStability = 'stable' | 'experimental';

// ---------------------------------------------------------------------------
// DSL Authoring Context types
// ---------------------------------------------------------------------------

/**
 * Whether a particular capability is present in the active LSP server.
 *
 * - `available`   -- the server reported support and returned data.
 * - `unavailable` -- the server explicitly does not support the capability.
 * - `unknown`     -- the server could not be reached, has not started, or
 *                    did not respond to a capability probe.
 */
export type CapabilityStatus = 'available' | 'unavailable' | 'unknown';

/**
 * High-level description of the domain-specific language handled by an LSP
 * server (e.g. "VASP INCAR", "Gaussian route section").
 */
export interface LanguageDescription {
  /** Human-readable name of the language (e.g. "VASP INCAR"). */
  readonly name: string;
  /** Short prose summary of the language's purpose. */
  readonly summary: string;
  /** URI or path to reference documentation, if known. */
  readonly documentationUri?: string;
}

/**
 * Schema metadata for the section or keyword surrounding the user's cursor.
 */
export interface SectionKeywordSchema {
  /** The section or keyword name (e.g. "FORCE_EVAL", "ENCUT"). */
  readonly name: string;
  /** Short description of what this keyword controls. */
  readonly description: string;
  /** Allowed values (if enumerable). */
  readonly allowedValues?: readonly string[];
  /** Default value. */
  readonly defaultValue?: string;
  /** Value type hint (e.g. "integer", "real", "string"). */
  readonly type?: string;
}

/**
 * A minimal, self-contained example for the detected language or keyword.
 */
export interface MinimalExample {
  /** Short label for the example (e.g. "Geometry optimization"). */
  readonly title: string;
  /** The example input snippet. */
  readonly code: string;
  /** Optional description of what the example demonstrates. */
  readonly description?: string;
}

/**
 * Next-token or completion guidance relevant at the cursor position.
 */
export interface NextTokenGuidance {
  /** Ordered list of candidate tokens the user might type next. */
  readonly candidates: readonly string[];
  /** Prose hint about what kind of token is expected. */
  readonly hint?: string;
}

/**
 * A single diagnostic entry surfaced through the context bundle.
 */
export interface ContextDiagnostic {
  /** 0-based line number. */
  readonly line: number;
  /** 0-based character offset. */
  readonly character?: number;
  /** Diagnostic message. */
  readonly message: string;
  /** Severity level. */
  readonly severity: 'error' | 'warning' | 'information' | 'hint';
}

/**
 * A code-action entry surfaced through the context bundle.
 */
export interface ContextCodeAction {
  /** Human-readable title of the action. */
  readonly title: string;
  /** Kind of the action (e.g. "quickfix", "refactor"). */
  readonly kind?: string;
  /** Whether the action is preferred. */
  readonly isPreferred?: boolean;
}

/**
 * The complete DSL authoring context bundle returned to calling agents.
 *
 * Every capability field is always present; optional capabilities that are
 * not supported degrade gracefully with an explicit status marker and empty
 * data.
 */
export interface DSLAuthoringContext {
  /** URI of the document this context describes. */
  readonly documentUri: string;
  /** VS Code language ID (e.g. "gaussian", "vasp"). */
  readonly languageId: string;
  /** Registry ID of the detected LSP server (e.g. "gaussian-lsp"). */
  readonly serverId: string;
  /** Human-readable software name (e.g. "Gaussian"). */
  readonly serverName: string;
  /** Stability of the server. */
  readonly stability: LSPStability;
  /** Whether the LSP server is currently running. */
  readonly serverRunning: boolean;

  /** Domain language description. */
  readonly languageDescription: {
    readonly status: CapabilityStatus;
    readonly data?: LanguageDescription;
  };

  /** Section / keyword schema at the cursor position. */
  readonly sectionKeywordSchema: {
    readonly status: CapabilityStatus;
    readonly data?: SectionKeywordSchema;
  };

  /** Minimal usage examples. */
  readonly examples: {
    readonly status: CapabilityStatus;
    readonly items: readonly MinimalExample[];
  };

  /** Next-token guidance at the cursor position. */
  readonly nextTokenGuidance: {
    readonly status: CapabilityStatus;
    readonly data?: NextTokenGuidance;
  };

  /** Current diagnostics for the document. */
  readonly diagnostics: {
    readonly status: CapabilityStatus;
    readonly items: readonly ContextDiagnostic[];
  };

  /** Available code actions at the cursor position. */
  readonly codeActions: {
    readonly status: CapabilityStatus;
    readonly items: readonly ContextCodeAction[];
  };

  /** ISO-8601 timestamp when the context was assembled. */
  readonly assembledAt: string;
}

/**
 * Output format requested by the caller.
 */
export type DSLContextOutputFormat = 'json' | 'markdown';

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
