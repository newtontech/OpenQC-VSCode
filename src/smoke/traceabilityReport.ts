/**
 * LSP traceability report contract.
 *
 * Backend repositories emit this report after checking that scientific
 * docstrings link to LLM Wiki nodes and wiki source claims link to raw assets.
 *
 * @module smoke/traceabilityReport
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/233
 */

export const TRACEABILITY_REPORT_SCHEMA_VERSION = 'openqc.lsp.traceability.v1';

const RULE_ID_PATTERN = /^[A-Z0-9]+-[A-Z0-9_]+-[A-Z0-9_]+-\d{3}$/;

export interface TraceabilityReportSummary {
  readonly docstringsTotal: number;
  readonly docstringsLinked: number;
  readonly brokenWikiLinks: number;
  readonly wikiSourcesWithoutRaw: number;
  readonly rawManifestFailures: number;
}

export interface TraceabilityReportDocstring {
  readonly path: string;
  readonly symbol: string;
  readonly wikiPath: string;
}

export interface TraceabilityReportWikiSource {
  readonly wikiPath: string;
  readonly sourceUrl: string;
  readonly rawPath: string;
}

export interface TraceabilityReportRuleId {
  readonly code: string;
  readonly sourcePath: string;
}

export interface TraceabilityReportSourceUrl {
  readonly url: string;
  readonly rawPath: string;
  readonly kind?: string;
}

export interface TraceabilityReport {
  readonly schemaVersion: typeof TRACEABILITY_REPORT_SCHEMA_VERSION;
  readonly serverId: string;
  readonly repository: string;
  readonly languageId: string;
  readonly generatedAt: string;
  readonly summary: TraceabilityReportSummary;
  readonly docstrings: readonly TraceabilityReportDocstring[];
  readonly wikiSources: readonly TraceabilityReportWikiSource[];
  readonly ruleIds: readonly TraceabilityReportRuleId[];
  readonly sourceUrls: readonly TraceabilityReportSourceUrl[];
  readonly rawManifest: {
    readonly path: string;
    readonly ok: boolean;
  };
}

export interface TraceabilityValidationOptions {
  readonly expectedServerId?: string;
}

export function validateTraceabilityReport(
  data: unknown,
  options: TraceabilityValidationOptions = {}
): string[] {
  const errors: string[] = [];
  if (!isRecord(data)) {
    return ['Traceability report must be a non-null object'];
  }

  const schemaVersion = data.schemaVersion;
  if (schemaVersion !== TRACEABILITY_REPORT_SCHEMA_VERSION) {
    errors.push(
      `schemaVersion must be "${TRACEABILITY_REPORT_SCHEMA_VERSION}", got ${JSON.stringify(schemaVersion)}`
    );
  }

  const serverId = data.serverId;
  if (!nonEmptyString(serverId)) {
    errors.push('serverId must be a non-empty string');
  } else if (options.expectedServerId && serverId !== options.expectedServerId) {
    errors.push(`serverId mismatch: expected "${options.expectedServerId}", got "${serverId}"`);
  }

  for (const field of ['repository', 'languageId', 'generatedAt'] as const) {
    if (!nonEmptyString(data[field])) {
      errors.push(`${field} must be a non-empty string`);
    }
  }
  if (nonEmptyString(data.generatedAt) && Number.isNaN(Date.parse(data.generatedAt))) {
    errors.push('generatedAt must be an ISO-8601 timestamp');
  }

  validateSummary(data.summary, errors);
  validateDocstrings(data.docstrings, errors);
  validateWikiSources(data.wikiSources, errors);
  validateRuleIds(data.ruleIds, errors);
  validateSourceUrls(data.sourceUrls, errors);
  validateRawManifest(data.rawManifest, errors);

  return errors;
}

function validateSummary(value: unknown, errors: string[]): void {
  if (!isRecord(value)) {
    errors.push('summary must be an object');
    return;
  }
  for (const field of [
    'docstringsTotal',
    'docstringsLinked',
    'brokenWikiLinks',
    'wikiSourcesWithoutRaw',
    'rawManifestFailures',
  ] as const) {
    if (!nonNegativeInteger(value[field])) {
      errors.push(`summary.${field} must be a non-negative integer`);
    }
  }
  if (
    nonNegativeInteger(value.docstringsTotal) &&
    nonNegativeInteger(value.docstringsLinked) &&
    value.docstringsLinked > value.docstringsTotal
  ) {
    errors.push('summary.docstringsLinked cannot exceed summary.docstringsTotal');
  }
}

function validateDocstrings(value: unknown, errors: string[]): void {
  if (!Array.isArray(value)) {
    errors.push('docstrings must be an array');
    return;
  }
  value.forEach((entry, index) => {
    if (!isRecord(entry)) {
      errors.push(`docstrings[${index}] must be an object`);
      return;
    }
    requireNonEmpty(entry, 'path', `docstrings[${index}]`, errors);
    requireNonEmpty(entry, 'symbol', `docstrings[${index}]`, errors);
    requireNonEmpty(entry, 'wikiPath', `docstrings[${index}]`, errors);
    requireRepoRelativePath(entry, 'path', `docstrings[${index}]`, errors);
    requireRepoRelativePath(entry, 'wikiPath', `docstrings[${index}]`, errors);
  });
}

function validateWikiSources(value: unknown, errors: string[]): void {
  if (!Array.isArray(value)) {
    errors.push('wikiSources must be an array');
    return;
  }
  value.forEach((entry, index) => {
    if (!isRecord(entry)) {
      errors.push(`wikiSources[${index}] must be an object`);
      return;
    }
    requireNonEmpty(entry, 'wikiPath', `wikiSources[${index}]`, errors);
    requireNonEmpty(entry, 'sourceUrl', `wikiSources[${index}]`, errors);
    requireNonEmpty(entry, 'rawPath', `wikiSources[${index}]`, errors);
    requireRepoRelativePath(entry, 'wikiPath', `wikiSources[${index}]`, errors);
    requireRepoRelativePath(entry, 'rawPath', `wikiSources[${index}]`, errors);
  });
}

function validateRuleIds(value: unknown, errors: string[]): void {
  if (!Array.isArray(value)) {
    errors.push('ruleIds must be an array');
    return;
  }
  value.forEach((entry, index) => {
    if (!isRecord(entry)) {
      errors.push(`ruleIds[${index}] must be an object`);
      return;
    }
    if (!nonEmptyString(entry.code)) {
      errors.push(`ruleIds[${index}].code must be a non-empty string`);
    } else if (!RULE_ID_PATTERN.test(entry.code)) {
      errors.push(`ruleIds[${index}].code must match <BACKEND>-<FILE_ROLE>-<CATEGORY>-NNN`);
    }
    requireNonEmpty(entry, 'sourcePath', `ruleIds[${index}]`, errors);
    requireRepoRelativePath(entry, 'sourcePath', `ruleIds[${index}]`, errors);
  });
}

function validateSourceUrls(value: unknown, errors: string[]): void {
  if (!Array.isArray(value)) {
    errors.push('sourceUrls must be an array');
    return;
  }
  value.forEach((entry, index) => {
    if (!isRecord(entry)) {
      errors.push(`sourceUrls[${index}] must be an object`);
      return;
    }
    requireNonEmpty(entry, 'url', `sourceUrls[${index}]`, errors);
    requireNonEmpty(entry, 'rawPath', `sourceUrls[${index}]`, errors);
    requireRepoRelativePath(entry, 'rawPath', `sourceUrls[${index}]`, errors);
  });
}

function validateRawManifest(value: unknown, errors: string[]): void {
  if (!isRecord(value)) {
    errors.push('rawManifest must be an object');
    return;
  }
  requireNonEmpty(value, 'path', 'rawManifest', errors);
  requireRepoRelativePath(value, 'path', 'rawManifest', errors);
  if (typeof value.ok !== 'boolean') {
    errors.push('rawManifest.ok must be a boolean');
  }
}

function requireNonEmpty(
  object: Record<string, unknown>,
  field: string,
  prefix: string,
  errors: string[]
): void {
  if (!nonEmptyString(object[field])) {
    errors.push(`${prefix}.${field} must be a non-empty string`);
  }
}

function requireRepoRelativePath(
  object: Record<string, unknown>,
  field: string,
  prefix: string,
  errors: string[]
): void {
  const value = object[field];
  if (nonEmptyString(value) && !isRepoRelativePath(value)) {
    errors.push(`${prefix}.${field} must be a repository-relative path`);
  }
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return !!value && typeof value === 'object' && !Array.isArray(value);
}

function nonEmptyString(value: unknown): value is string {
  return typeof value === 'string' && value.trim().length > 0;
}

function nonNegativeInteger(value: unknown): value is number {
  return typeof value === 'number' && Number.isInteger(value) && value >= 0;
}

function isRepoRelativePath(value: string): boolean {
  return (
    !value.startsWith('/') &&
    !/^[A-Za-z]:[\\/]/.test(value) &&
    !value.startsWith('file:') &&
    !value.split(/[\\/]+/).includes('..')
  );
}
