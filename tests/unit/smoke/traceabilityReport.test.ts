/**
 * Traceability report contract tests.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/233
 */

import * as fs from 'fs';
import * as path from 'path';
import {
  TRACEABILITY_REPORT_SCHEMA_VERSION,
  validateTraceabilityReport,
} from '../../../src/smoke/traceabilityReport';

const FIXTURES_DIR = path.join(__dirname, '../../fixtures/smoke/traceability');

function loadFixture(name: string): any {
  return JSON.parse(fs.readFileSync(path.join(FIXTURES_DIR, name), 'utf8'));
}

function clone<T>(value: T): T {
  return JSON.parse(JSON.stringify(value));
}

describe('validateTraceabilityReport', () => {
  it.each([
    ['cp2k-lsp-enhanced-report.json', 'cp2k-lsp-enhanced'],
    ['dpgen-lsp-report.json', 'dpgen-lsp'],
  ])('validates %s fixture', (fixtureName, serverId) => {
    const report = loadFixture(fixtureName);

    expect(report.schemaVersion).toBe(TRACEABILITY_REPORT_SCHEMA_VERSION);
    expect(validateTraceabilityReport(report, { expectedServerId: serverId })).toEqual([]);
  });

  it('reports server ID mismatches', () => {
    const report = loadFixture('cp2k-lsp-enhanced-report.json');

    const errors = validateTraceabilityReport(report, { expectedServerId: 'dpgen-lsp' });

    expect(errors).toContain('serverId mismatch: expected "dpgen-lsp", got "cp2k-lsp-enhanced"');
  });

  it('rejects malformed traceability rule IDs', () => {
    const report = clone(loadFixture('dpgen-lsp-report.json'));
    report.ruleIds[0].code = 'DPGEN_SCHEMA_001';

    const errors = validateTraceabilityReport(report, { expectedServerId: 'dpgen-lsp' });

    expect(errors).toContain('ruleIds[0].code must match <BACKEND>-<FILE_ROLE>-<CATEGORY>-NNN');
  });

  it('rejects hidden local absolute paths', () => {
    const report = clone(loadFixture('cp2k-lsp-enhanced-report.json'));
    report.docstrings[0].path = '/Users/example/private/cp2k_input_tools/diagnostics.py';

    const errors = validateTraceabilityReport(report, { expectedServerId: 'cp2k-lsp-enhanced' });

    expect(errors).toContain('docstrings[0].path must be a repository-relative path');
  });
});
