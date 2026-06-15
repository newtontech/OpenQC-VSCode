/**
 * Smoke Test and Verification Types
 *
 * Defines types for the OpenQC smoke test infrastructure including rule manifests,
 * compatibility reports, verification results, and campaign reports.
 *
 * @module smoke/types
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/159
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/160
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/161
 */

// ---------------------------------------------------------------------------
// Rule Manifest types (Issues #160, #162)
// ---------------------------------------------------------------------------

/** Severity level for a diagnostic rule. */
export type RuleSeverity = 'error' | 'warning' | 'information' | 'hint';

/** A single diagnostic rule exported by an LSP repository. */
export interface RuleManifestEntry {
  /** Unique rule code, e.g. "E001", "W010". */
  readonly code: string;
  /** Human-readable description of what this rule checks. */
  readonly message: string;
  /** Severity level. */
  readonly severity: RuleSeverity;
  /** Category grouping, e.g. "syntax", "semantics", "best-practice". */
  readonly category?: string;
  /** Whether the rule blocks downstream processing. */
  readonly blocking?: boolean;
  /** Confidence level 0-1 for automated application. */
  readonly confidence?: number;
}

/**
 * The full rule manifest exported by an LSP repository via
 * `lsp-rules export --json` or equivalent.
 */
export interface RuleManifest {
  /** LSP server registry ID, e.g. "gaussian-lsp". */
  readonly serverId: string;
  /** Human-readable software name, e.g. "Gaussian". */
  readonly serverName: string;
  /** Version of the manifest schema. */
  readonly schemaVersion: string;
  /** ISO-8601 timestamp when the manifest was generated. */
  readonly generatedAt: string;
  /** Diagnostic rules in this manifest. */
  readonly rules: readonly RuleManifestEntry[];
  /** Optional metadata about the export source. */
  readonly source?: {
    readonly repository?: string;
    readonly commit?: string;
    readonly branch?: string;
  };
}

// ---------------------------------------------------------------------------
// Compatibility Report types (Issue #163)
// ---------------------------------------------------------------------------

/** Summary counts exported by backend docstring/wiki/raw traceability reports. */
export interface TraceabilityReportSummary {
  readonly docstringsTotal: number;
  readonly docstringsLinked: number;
  readonly brokenWikiLinks: number;
  readonly wikiSourcesWithoutRaw: number;
  readonly rawManifestFailures: number;
}

/** Full backend provenance traceability report consumed by OpenQC. */
export interface TraceabilityReport {
  readonly schemaVersion: 'openqc.lsp.traceability.v1';
  readonly serverId: string;
  readonly repository: string;
  readonly languageId: string;
  readonly generatedAt: string;
  readonly summary: TraceabilityReportSummary;
  readonly docstrings: readonly {
    readonly path: string;
    readonly symbol: string;
    readonly wikiPath: string;
  }[];
  readonly wikiSources: readonly {
    readonly wikiPath: string;
    readonly sourceUrl: string;
    readonly rawPath: string;
  }[];
  readonly ruleIds: readonly {
    readonly code: string;
    readonly sourcePath: string;
  }[];
  readonly sourceUrls: readonly {
    readonly url: string;
    readonly rawPath: string;
    readonly kind?: string;
  }[];
  readonly rawManifest: {
    readonly path: string;
    readonly ok: boolean;
  };
}

/** Status of an individual compatibility check. */
export type CheckStatus = 'pass' | 'fail' | 'skip' | 'warn';

/** Result of a single compatibility check. */
export interface CompatibilityCheck {
  /** Name of the check, e.g. "registry-entry-exists". */
  readonly name: string;
  /** Human-readable description of what was checked. */
  readonly description: string;
  /** Result status. */
  readonly status: CheckStatus;
  /** Optional detail message (e.g. failure reason). */
  readonly detail?: string;
}

/** Per-LSP compatibility report entry. */
export interface LspCompatibilityEntry {
  /** Registry ID. */
  readonly serverId: string;
  /** Human-readable name. */
  readonly serverName: string;
  /** Language ID. */
  readonly languageId: string;
  /** Stability classification. */
  readonly stability: string;
  /** Individual check results. */
  readonly checks: readonly CompatibilityCheck[];
  /** Overall pass/fail. */
  readonly passed: boolean;
}

/** Full compatibility report. */
export interface CompatibilityReport {
  /** ISO-8601 timestamp. */
  readonly generatedAt: string;
  /** Total number of LSP servers checked. */
  readonly totalServers: number;
  /** Number that passed all checks. */
  readonly passedServers: number;
  /** Number that failed or had warnings. */
  readonly failedServers: number;
  /** Per-LSP details. */
  readonly entries: readonly LspCompatibilityEntry[];
}

// ---------------------------------------------------------------------------
// Output/Log Document Detection types (Issue #159)
// ---------------------------------------------------------------------------

/** Classification of a document type for diagnostics routing. */
export type DocumentKind = 'input' | 'output' | 'log' | 'unknown';

/** Result of output/log document detection. */
export interface DocumentDetectionResult {
  /** The detected document kind. */
  readonly kind: DocumentKind;
  /** The registry ID of the associated LSP server, if any. */
  readonly serverId?: string;
  /** The language ID of the associated input language, if any. */
  readonly languageId?: string;
  /** Confidence score 0-1. */
  readonly confidence: number;
}

// ---------------------------------------------------------------------------
// Smoke Test types (Issue #161)
// ---------------------------------------------------------------------------

/** Result of a smoke test run for a single LSP. */
export interface SmokeTestResult {
  /** Registry ID. */
  readonly serverId: string;
  /** Whether valid input fixture passed. */
  readonly validInputPass: boolean;
  /** Whether invalid input fixture was correctly detected. */
  readonly invalidInputPass: boolean;
  /** Whether runtime log fixture was processed (or skipped if N/A). */
  readonly runtimeLogPass: boolean | null;
  /** Detail message if any check failed. */
  readonly detail?: string;
}

/** Overall smoke test summary. */
export interface SmokeTestSummary {
  /** ISO-8601 timestamp. */
  readonly generatedAt: string;
  /** Total LSP servers tested. */
  readonly totalServers: number;
  /** Number that passed all applicable smoke tests. */
  readonly passedServers: number;
  /** Per-LSP results. */
  readonly results: readonly SmokeTestResult[];
}

// ---------------------------------------------------------------------------
// VSIX Verification types (Issue #166)
// ---------------------------------------------------------------------------

/** Result of VSIX package verification. */
export interface VsixVerificationResult {
  /** Whether the VSIX was built successfully. */
  readonly buildSuccess: boolean;
  /** Whether vscode-languageclient is listed as a dependency. */
  readonly hasLanguageClient: boolean;
  /** List of required dependencies found. */
  readonly foundDependencies: readonly string[];
  /** List of missing dependencies. */
  readonly missingDependencies: readonly string[];
  /** VSIX file size in bytes, if built. */
  readonly fileSize?: number;
  /** Detail message. */
  readonly detail?: string;
}

// ---------------------------------------------------------------------------
// Final Verification types (Issues #168-#177)
// ---------------------------------------------------------------------------

/** GitHub repository status check result. */
export interface RepoGitHubStatus {
  /** Repository full name, e.g. "newtontech/gaussian-lsp". */
  readonly repository: string;
  /** Number of open issues. */
  readonly openIssues: number;
  /** Number of open PRs. */
  readonly openPullRequests: number;
  /** Whether the repo passed (zero open issues and PRs). */
  readonly passed: boolean;
}

/** Result of running gates for an LSP repo. */
export interface RepoGateResult {
  /** Repository full name. */
  readonly repository: string;
  /** Whether the local checkout exists. */
  readonly checkoutExists: boolean;
  /** Whether gate tests passed. */
  readonly gatePassed: boolean;
  /** Detail message. */
  readonly detail?: string;
}

/** Result of git cleanliness check. */
export interface RepoCleanStatus {
  /** Repository full name. */
  readonly repository: string;
  /** Whether the checkout exists. */
  readonly exists: boolean;
  /** Whether the working tree is clean. */
  readonly isClean: boolean;
  /** Output of `git status --short` if dirty. */
  readonly statusOutput?: string;
}

/** Final campaign report. */
export interface CampaignReport {
  /** ISO-8601 timestamp. */
  readonly generatedAt: string;
  /** Run identifier. */
  readonly runId: string;
  /** Whether all GitHub checks passed. */
  readonly githubChecksPassed: boolean;
  /** Whether all local gates passed. */
  readonly localGatesPassed: boolean;
  /** Whether OpenQC CI passed. */
  readonly openqcCiPassed: boolean;
  /** Whether VSIX build succeeded. */
  readonly vsixBuildPassed: boolean;
  /** Whether all repos are clean. */
  readonly reposClean: boolean;
  /** Per-repo GitHub status. */
  readonly repoStatuses: readonly RepoGitHubStatus[];
  /** Per-repo gate results. */
  readonly gateResults: readonly RepoGateResult[];
  /** Per-repo cleanliness. */
  readonly cleanStatuses: readonly RepoCleanStatus[];
  /** Compatibility report summary. */
  readonly compatibilityReport?: CompatibilityReport;
  /** Smoke test summary. */
  readonly smokeTestSummary?: SmokeTestSummary;
  /** VSIX verification result. */
  readonly vsixVerification?: VsixVerificationResult;
  /** Overall pass/fail. */
  readonly overallPassed: boolean;
}
