/**
 * Rule Manifest Reader Tests
 *
 * Tests for manifest parsing, validation, and query helpers.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/160
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/162
 */

import * as path from 'path';
import {
  validateManifest,
  readManifestFromFile,
  parseManifestString,
  tryReadManifest,
  getRulesBySeverity,
  getRuleByCode,
  getManifestCategories,
  countRulesBySeverity,
} from '../../../src/smoke/ruleManifestReader';
import type { RuleManifest } from '../../../src/smoke/types';

const FIXTURES_DIR = path.join(__dirname, '../../fixtures/smoke/manifests');

// ---------------------------------------------------------------------------
// validateManifest
// ---------------------------------------------------------------------------

describe('validateManifest', () => {
  it('returns no errors for a valid manifest', () => {
    const valid = {
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '1.0.0',
      generatedAt: '2026-06-12T00:00:00Z',
      rules: [{ code: 'E001', message: 'Test error', severity: 'error' }],
    };
    expect(validateManifest(valid)).toEqual([]);
  });

  it('rejects null input', () => {
    const errors = validateManifest(null);
    expect(errors).toHaveLength(1);
    expect(errors[0]).toContain('non-null object');
  });

  it('rejects non-object input', () => {
    const errors = validateManifest('string');
    expect(errors).toHaveLength(1);
    expect(errors[0]).toContain('non-null object');
  });

  it('requires serverId', () => {
    const errors = validateManifest({
      serverName: 'Test',
      schemaVersion: '1.0.0',
      generatedAt: '2026-01-01',
      rules: [],
    });
    expect(errors.some(e => e.includes('serverId'))).toBe(true);
  });

  it('requires serverName', () => {
    const errors = validateManifest({
      serverId: 'test-lsp',
      schemaVersion: '1.0.0',
      generatedAt: '2026-01-01',
      rules: [],
    });
    expect(errors.some(e => e.includes('serverName'))).toBe(true);
  });

  it('requires supported schemaVersion', () => {
    const errors = validateManifest({
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '99.0.0',
      generatedAt: '2026-01-01',
      rules: [],
    });
    expect(errors.some(e => e.includes('schema version'))).toBe(true);
  });

  it('requires rules to be an array', () => {
    const errors = validateManifest({
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '1.0.0',
      generatedAt: '2026-01-01',
      rules: 'not-array',
    });
    expect(errors.some(e => e.includes('"rules" array'))).toBe(true);
  });

  it('validates individual rule entries', () => {
    const errors = validateManifest({
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '1.0.0',
      generatedAt: '2026-01-01',
      rules: [{ code: '', message: '', severity: 'bad' }],
    });
    expect(errors.some(e => e.includes('code'))).toBe(true);
    expect(errors.some(e => e.includes('message'))).toBe(true);
    expect(errors.some(e => e.includes('severity'))).toBe(true);
  });

  it('accepts schema version 1.1.0', () => {
    const errors = validateManifest({
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '1.1.0',
      generatedAt: '2026-01-01',
      rules: [],
    });
    expect(errors).toEqual([]);
  });

  it('validates optional confidence range', () => {
    const errors = validateManifest({
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '1.0.0',
      generatedAt: '2026-01-01',
      rules: [{ code: 'E001', message: 'Test', severity: 'error', confidence: 1.5 }],
    });
    expect(errors.some(e => e.includes('confidence'))).toBe(true);
  });

  it('validates optional blocking is boolean', () => {
    const errors = validateManifest({
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '1.0.0',
      generatedAt: '2026-01-01',
      rules: [{ code: 'E001', message: 'Test', severity: 'error', blocking: 'yes' }],
    });
    expect(errors.some(e => e.includes('blocking'))).toBe(true);
  });
});

// ---------------------------------------------------------------------------
// parseManifestString
// ---------------------------------------------------------------------------

describe('parseManifestString', () => {
  it('parses a valid JSON manifest', () => {
    const json = JSON.stringify({
      serverId: 'test-lsp',
      serverName: 'Test',
      schemaVersion: '1.0.0',
      generatedAt: '2026-06-12T00:00:00Z',
      rules: [{ code: 'E001', message: 'Error', severity: 'error' }],
    });

    const manifest = parseManifestString(json);
    expect(manifest.serverId).toBe('test-lsp');
    expect(manifest.rules).toHaveLength(1);
  });

  it('throws on invalid JSON', () => {
    expect(() => parseManifestString('not-json')).toThrow('Failed to parse');
  });

  it('throws on invalid manifest structure', () => {
    expect(() => parseManifestString('{}')).toThrow('validation failed');
  });
});

// ---------------------------------------------------------------------------
// readManifestFromFile
// ---------------------------------------------------------------------------

describe('readManifestFromFile', () => {
  it('reads the Gaussian fixture manifest', () => {
    const manifest = readManifestFromFile(path.join(FIXTURES_DIR, 'gaussian-lsp-rules.json'));
    expect(manifest.serverId).toBe('gaussian-lsp');
    expect(manifest.serverName).toBe('Gaussian');
    expect(manifest.rules).toHaveLength(4);
  });

  it('reads the VASP fixture manifest', () => {
    const manifest = readManifestFromFile(path.join(FIXTURES_DIR, 'vasp-lsp-rules.json'));
    expect(manifest.serverId).toBe('vasp-lsp');
    expect(manifest.rules).toHaveLength(2);
  });

  it('throws for a missing file', () => {
    expect(() => readManifestFromFile(path.join(FIXTURES_DIR, 'nonexistent.json'))).toThrow(
      'not found'
    );
  });

  it('throws for an invalid manifest file', () => {
    expect(() => readManifestFromFile(path.join(FIXTURES_DIR, 'invalid-rules.json'))).toThrow(
      'validation failed'
    );
  });
});

// ---------------------------------------------------------------------------
// tryReadManifest
// ---------------------------------------------------------------------------

describe('tryReadManifest', () => {
  it('returns manifest on success', () => {
    const result = tryReadManifest(path.join(FIXTURES_DIR, 'gaussian-lsp-rules.json'));
    expect(result.manifest).toBeDefined();
    expect(result.error).toBeUndefined();
    expect(result.manifest!.serverId).toBe('gaussian-lsp');
  });

  it('returns error on failure', () => {
    const result = tryReadManifest(path.join(FIXTURES_DIR, 'nonexistent.json'));
    expect(result.manifest).toBeUndefined();
    expect(result.error).toBeDefined();
    expect(result.error).toContain('not found');
  });
});

// ---------------------------------------------------------------------------
// Query helpers
// ---------------------------------------------------------------------------

describe('getRulesBySeverity', () => {
  const manifest: RuleManifest = {
    serverId: 'test-lsp',
    serverName: 'Test',
    schemaVersion: '1.0.0',
    generatedAt: '2026-06-12T00:00:00Z',
    rules: [
      { code: 'E001', message: 'Error 1', severity: 'error' },
      { code: 'E002', message: 'Error 2', severity: 'error' },
      { code: 'W001', message: 'Warning 1', severity: 'warning' },
      { code: 'I001', message: 'Info 1', severity: 'information' },
    ],
  };

  it('filters rules by error severity', () => {
    const errors = getRulesBySeverity(manifest, 'error');
    expect(errors).toHaveLength(2);
    expect(errors.map(r => r.code)).toEqual(['E001', 'E002']);
  });

  it('returns empty for severity with no rules', () => {
    const hints = getRulesBySeverity(manifest, 'hint');
    expect(hints).toHaveLength(0);
  });
});

describe('getRuleByCode', () => {
  const manifest: RuleManifest = {
    serverId: 'test-lsp',
    serverName: 'Test',
    schemaVersion: '1.0.0',
    generatedAt: '2026-06-12T00:00:00Z',
    rules: [
      { code: 'E001', message: 'Error', severity: 'error' },
      { code: 'W001', message: 'Warning', severity: 'warning' },
    ],
  };

  it('finds a rule by code', () => {
    const rule = getRuleByCode(manifest, 'E001');
    expect(rule).toBeDefined();
    expect(rule?.message).toBe('Error');
  });

  it('returns undefined for unknown code', () => {
    expect(getRuleByCode(manifest, 'X999')).toBeUndefined();
  });
});

describe('getManifestCategories', () => {
  const manifest: RuleManifest = {
    serverId: 'test-lsp',
    serverName: 'Test',
    schemaVersion: '1.0.0',
    generatedAt: '2026-06-12T00:00:00Z',
    rules: [
      { code: 'E001', message: 'A', severity: 'error', category: 'syntax' },
      { code: 'W001', message: 'B', severity: 'warning', category: 'best-practice' },
      { code: 'E002', message: 'C', severity: 'error', category: 'syntax' },
      { code: 'E003', message: 'D', severity: 'error' },
    ],
  };

  it('returns unique sorted categories', () => {
    const categories = getManifestCategories(manifest);
    expect(categories).toEqual(['best-practice', 'syntax']);
  });
});

describe('countRulesBySeverity', () => {
  const manifest: RuleManifest = {
    serverId: 'test-lsp',
    serverName: 'Test',
    schemaVersion: '1.0.0',
    generatedAt: '2026-06-12T00:00:00Z',
    rules: [
      { code: 'E001', message: 'A', severity: 'error' },
      { code: 'E002', message: 'B', severity: 'error' },
      { code: 'W001', message: 'C', severity: 'warning' },
      { code: 'I001', message: 'D', severity: 'information' },
    ],
  };

  it('counts rules by severity', () => {
    const counts = countRulesBySeverity(manifest);
    expect(counts['error']).toBe(2);
    expect(counts['warning']).toBe(1);
    expect(counts['information']).toBe(1);
    expect(counts['hint']).toBe(0);
  });
});
