/**
 * Rule Manifest Reader
 *
 * Reads and validates rule manifests exported by each LSP repository.
 * Accepts `lsp-rules export --json` output from any LSP repo without
 * embedding rule logic in OpenQC itself.
 *
 * @module smoke/ruleManifestReader
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/160
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/162
 */

import * as fs from 'fs';
import * as path from 'path';
import type { RuleManifest, RuleManifestEntry, RuleSeverity } from './types';

// ---------------------------------------------------------------------------
// Schema validation
// ---------------------------------------------------------------------------

const VALID_SEVERITIES = new Set<string>(['error', 'warning', 'information', 'hint']);
const SCHEMA_VERSIONS = new Set<string>(['1.0.0', '1.1.0']);

/**
 * Validate a parsed rule manifest object against the expected schema.
 *
 * Returns an array of validation error strings. An empty array means the
 * manifest is valid.
 *
 * @param data - The parsed JSON object to validate.
 * @returns Array of validation error messages.
 */
export function validateManifest(data: unknown): string[] {
  const errors: string[] = [];

  if (!data || typeof data !== 'object') {
    errors.push('Manifest must be a non-null object');
    return errors;
  }

  const manifest = data as Record<string, unknown>;

  if (typeof manifest.serverId !== 'string' || manifest.serverId.length === 0) {
    errors.push('Manifest must have a non-empty string "serverId"');
  }

  if (typeof manifest.serverName !== 'string' || manifest.serverName.length === 0) {
    errors.push('Manifest must have a non-empty string "serverName"');
  }

  if (typeof manifest.schemaVersion !== 'string') {
    errors.push('Manifest must have a string "schemaVersion"');
  } else if (!SCHEMA_VERSIONS.has(manifest.schemaVersion)) {
    errors.push(
      `Unsupported schema version "${manifest.schemaVersion}". ` +
        `Supported: ${Array.from(SCHEMA_VERSIONS).join(', ')}`
    );
  }

  if (typeof manifest.generatedAt !== 'string') {
    errors.push('Manifest must have a string "generatedAt" timestamp');
  }

  if (!Array.isArray(manifest.rules)) {
    errors.push('Manifest must have a "rules" array');
  } else {
    const ruleErrors = validateRules(manifest.rules as unknown[]);
    errors.push(...ruleErrors);
  }

  return errors;
}

/**
 * Validate the rules array entries.
 */
function validateRules(rules: unknown[]): string[] {
  const errors: string[] = [];

  for (let i = 0; i < rules.length; i++) {
    const rule = rules[i];
    if (!rule || typeof rule !== 'object') {
      errors.push(`Rule at index ${i} must be a non-null object`);
      continue;
    }

    const r = rule as Record<string, unknown>;

    if (typeof r.code !== 'string' || r.code.length === 0) {
      errors.push(`Rule at index ${i} must have a non-empty string "code"`);
    }

    if (typeof r.message !== 'string' || r.message.length === 0) {
      errors.push(`Rule at index ${i} must have a non-empty string "message"`);
    }

    if (typeof r.severity !== 'string' || !VALID_SEVERITIES.has(r.severity)) {
      errors.push(
        `Rule at index ${i} has invalid "severity": "${r.severity}". ` +
          `Must be one of: ${Array.from(VALID_SEVERITIES).join(', ')}`
      );
    }

    if (r.category !== undefined && typeof r.category !== 'string') {
      errors.push(`Rule at index ${i} "category" must be a string if present`);
    }

    if (r.blocking !== undefined && typeof r.blocking !== 'boolean') {
      errors.push(`Rule at index ${i} "blocking" must be a boolean if present`);
    }

    if (r.confidence !== undefined) {
      if (typeof r.confidence !== 'number' || r.confidence < 0 || r.confidence > 1) {
        errors.push(`Rule at index ${i} "confidence" must be a number between 0 and 1`);
      }
    }
  }

  return errors;
}

// ---------------------------------------------------------------------------
// Reader functions
// ---------------------------------------------------------------------------

/**
 * Read and parse a rule manifest from a JSON file path.
 *
 * @param filePath - Absolute or relative path to the manifest JSON file.
 * @returns The parsed and validated rule manifest.
 * @throws Error if the file cannot be read or validation fails.
 */
export function readManifestFromFile(filePath: string): RuleManifest {
  const absolutePath = path.resolve(filePath);

  if (!fs.existsSync(absolutePath)) {
    throw new Error(`Rule manifest file not found: ${absolutePath}`);
  }

  const raw = fs.readFileSync(absolutePath, 'utf8');
  return parseManifestString(raw);
}

/**
 * Parse and validate a rule manifest from a JSON string.
 *
 * @param json - The raw JSON string to parse.
 * @returns The parsed and validated rule manifest.
 * @throws Error if parsing or validation fails.
 */
export function parseManifestString(json: string): RuleManifest {
  let data: unknown;
  try {
    data = JSON.parse(json);
  } catch (parseError) {
    const msg = parseError instanceof Error ? parseError.message : String(parseError);
    throw new Error(`Failed to parse rule manifest JSON: ${msg}`);
  }

  const errors = validateManifest(data);
  if (errors.length > 0) {
    throw new Error(`Rule manifest validation failed:\n  ${errors.join('\n  ')}`);
  }

  return data as RuleManifest;
}

/**
 * Attempt to read a manifest, returning a result object instead of throwing.
 *
 * @param filePath - Path to the manifest JSON file.
 * @returns An object with either `manifest` on success or `error` on failure.
 */
export function tryReadManifest(filePath: string):
  | { readonly manifest: RuleManifest; readonly error?: undefined }
  | {
      readonly manifest?: undefined;
      readonly error: string;
    } {
  try {
    const manifest = readManifestFromFile(filePath);
    return { manifest };
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    return { error: msg };
  }
}

// ---------------------------------------------------------------------------
// Query helpers
// ---------------------------------------------------------------------------

/**
 * Get all rules from a manifest filtered by severity.
 *
 * @param manifest - The rule manifest to query.
 * @param severity - The severity level to filter by.
 * @returns Matching rule entries.
 */
export function getRulesBySeverity(
  manifest: RuleManifest,
  severity: RuleSeverity
): readonly RuleManifestEntry[] {
  return manifest.rules.filter(rule => rule.severity === severity);
}

/**
 * Get a specific rule by its code.
 *
 * @param manifest - The rule manifest to query.
 * @param code - The rule code to look up.
 * @returns The matching rule entry, or undefined if not found.
 */
export function getRuleByCode(manifest: RuleManifest, code: string): RuleManifestEntry | undefined {
  return manifest.rules.find(rule => rule.code === code);
}

/**
 * Get all unique categories present in a manifest.
 *
 * @param manifest - The rule manifest to query.
 * @returns Sorted array of unique category strings.
 */
export function getManifestCategories(manifest: RuleManifest): readonly string[] {
  const categories = new Set<string>();
  for (const rule of manifest.rules) {
    if (rule.category) {
      categories.add(rule.category);
    }
  }
  return Array.from(categories).sort();
}

/**
 * Count rules by severity level.
 *
 * @param manifest - The rule manifest to query.
 * @returns Record mapping severity to count.
 */
export function countRulesBySeverity(manifest: RuleManifest): Readonly<Record<string, number>> {
  const counts: Record<string, number> = {};
  for (const severity of VALID_SEVERITIES) {
    counts[severity] = 0;
  }
  for (const rule of manifest.rules) {
    counts[rule.severity] = (counts[rule.severity] || 0) + 1;
  }
  return counts;
}
