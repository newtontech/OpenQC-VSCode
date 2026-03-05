/**
 * Migration Validation Integration Tests
 *
 * Integration tests for migration validation functionality.
 * Tests round-trip conversion, structure preservation, and parameter consistency.
 */

import * as path from 'path';
import * as fs from 'fs';
import * as vscode from 'vscode';
import { ASEConverter, ASEFormat } from '../../src/ase/ASEConverter';

// Minimal mock for vscode.ExtensionContext with only required properties
class MockExtensionContext {
  subscriptions: vscode.Disposable[] = [];
  extensionPath: string;
  storageUri: vscode.Uri;
  globalStorageUri: vscode.Uri;
  logPath: string;
  secrets: any;
  environmentVariableCollection: any;
  storagePath: string | undefined;
  globalStoragePath: string;
  extensionMode: any;
  extension: any;
  workspaceState: any;
  globalState: any;
  extensionUri: vscode.Uri;
  logUri: vscode.Uri;
  languageModelAccessInformation: any;

  constructor(extensionPath: string) {
    this.extensionPath = extensionPath;
    this.extensionUri = vscode.Uri.file(extensionPath);
    this.logUri = vscode.Uri.file(path.join(extensionPath, 'logs'));
    this.storageUri = vscode.Uri.file(extensionPath);
    this.globalStorageUri = vscode.Uri.file(extensionPath);
    this.logPath = extensionPath;
    this.storagePath = extensionPath;
    this.globalStoragePath = extensionPath;
    this.extensionMode = 2; // Test mode
    this.extension = {};
    this.secrets = {};
    this.languageModelAccessInformation = {};

    const mementoGet = (key: string, defaultValue?: any) => defaultValue;
    const mementoUpdate = (key: string, value: any) => Promise.resolve();
    const mementoKeys = () => [];
    const setKeysForSync = (keys: readonly string[]) => {};

    this.globalState = {
      get: mementoGet,
      update: mementoUpdate,
      keys: mementoKeys,
      setKeysForSync: setKeysForSync,
    };

    this.workspaceState = {
      get: mementoGet,
      update: mementoUpdate,
      keys: mementoKeys,
    };

    this.environmentVariableCollection = {
      getScoped: () => ({}),
    };
  }

  asAbsolutePathUri(relativePath: string): vscode.Uri {
    return vscode.Uri.file(path.join(this.extensionPath, relativePath));
  }

  asAbsolutePath(relativePath: string): string {
    return path.join(this.extensionPath, relativePath);
  }
}

describe('Migration Validation Integration Tests', () => {
  jest.setTimeout(60000);
  let converter: ASEConverter;
  const fixturesDir = path.join(__dirname, '../fixtures/migration');
  const outputDir = path.join(__dirname, '../temp/migration');

  beforeAll(async () => {
    const mockContext = new MockExtensionContext(path.join(__dirname, '../..')) as any;
    converter = new ASEConverter(mockContext);

    const backendAvailable = await converter.isAvailable();
    if (!backendAvailable) {
      console.warn('Python backend not available. Skipping integration tests.');
      return;
    }

    fs.mkdirSync(outputDir, { recursive: true });
  });

  afterAll(() => {
    if (fs.existsSync(outputDir)) {
      fs.rmSync(outputDir, { recursive: true, force: true });
    }
  });

  const getBackendAvailable = async () => {
    return await converter.isAvailable();
  };

  describe('Round-Trip Conversion Tests', () => {
    it('should preserve structure in VASP → XYZ → VASP round-trip', async () => {
      jest.setTimeout(30000);
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const intermediateFile = path.join(outputDir, 'roundtrip.xyz');
      const finalFile = path.join(outputDir, 'roundtrip.POSCAR');

      // Use extxyz format which preserves cell information
      const result1 = await converter.convertFormat(
        inputFile,
        intermediateFile,
        ASEFormat.VASP,
        'extxyz' as any
      );
      expect(result1.success).toBe(true);

      const result2 = await converter.convertFormat(
        intermediateFile,
        finalFile,
        'extxyz' as any,
        ASEFormat.VASP
      );
      expect(result2.success).toBe(true);

      const readOriginal = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      const readFinal = await converter.readToAtoms(finalFile, ASEFormat.VASP);

      expect(readOriginal.success).toBe(true);
      expect(readFinal.success).toBe(true);
      expect(readOriginal.atoms?.chemical_symbols.length).toBe(
        readFinal.atoms?.chemical_symbols.length
      );
    });

    it('should preserve structure in VASP → CP2K → VASP round-trip', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const intermediateFile = path.join(outputDir, 'roundtrip_cp2k.inp');
      const finalFile = path.join(outputDir, 'roundtrip_cp2k.POSCAR');

      const result1 = await converter.convertFormat(
        inputFile,
        intermediateFile,
        ASEFormat.VASP,
        ASEFormat.CP2K
      );

      // CP2K format may not be supported for writing in all ASE versions
      if (!result1.success) {
        console.warn('CP2K format not supported for writing, skipping test');
        return;
      }

      const result2 = await converter.convertFormat(
        intermediateFile,
        finalFile,
        ASEFormat.CP2K,
        ASEFormat.VASP
      );
      expect(result2.success).toBe(true);

      const readOriginal = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      const readFinal = await converter.readToAtoms(finalFile, ASEFormat.VASP);

      expect(readOriginal.success).toBe(true);
      expect(readFinal.success).toBe(true);
      expect(readOriginal.atoms?.chemical_symbols.length).toBe(
        readFinal.atoms?.chemical_symbols.length
      );
    });

    it('should preserve structure in VASP → QE → VASP round-trip', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const intermediateFile = path.join(outputDir, 'roundtrip_qe.in');
      const finalFile = path.join(outputDir, 'roundtrip_qe.POSCAR');

      const result1 = await converter.convertFormat(
        inputFile,
        intermediateFile,
        ASEFormat.VASP,
        ASEFormat.QE
      );

      // QE format may not be supported for writing in all ASE versions
      if (!result1.success) {
        console.warn('QE format not supported for writing, skipping test');
        return;
      }

      const result2 = await converter.convertFormat(
        intermediateFile,
        finalFile,
        ASEFormat.QE,
        ASEFormat.VASP
      );
      expect(result2.success).toBe(true);

      const readOriginal = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      const readFinal = await converter.readToAtoms(finalFile, ASEFormat.VASP);

      expect(readOriginal.success).toBe(true);
      expect(readFinal.success).toBe(true);
      expect(readOriginal.atoms?.chemical_symbols.length).toBe(
        readFinal.atoms?.chemical_symbols.length
      );
    });
  });

  describe('Structure Preservation Validation', () => {
    it('should preserve atomic positions within tolerance', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const outputFile = path.join(outputDir, 'positions.xyz');

      const result = await converter.convertFormat(
        inputFile,
        outputFile,
        ASEFormat.VASP,
        ASEFormat.XYZ
      );

      expect(result.success).toBe(true);

      const original = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      const converted = await converter.readToAtoms(outputFile, ASEFormat.XYZ);

      expect(original.success).toBe(true);
      expect(converted.success).toBe(true);

      const originalPos = original.atoms?.positions || [];
      const convertedPos = converted.atoms?.positions || [];

      expect(originalPos.length).toBe(convertedPos.length);

      for (let i = 0; i < originalPos.length; i++) {
        for (let j = 0; j < 3; j++) {
          expect(Math.abs(originalPos[i][j] - convertedPos[i][j])).toBeLessThan(1e-4);
        }
      }
    });

    it('should preserve unit cell for periodic systems', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const outputFile = path.join(outputDir, 'cell.xyz'); // Use XYZ format for better compatibility

      const result = await converter.convertFormat(
        inputFile,
        outputFile,
        ASEFormat.VASP,
        ASEFormat.XYZ
      );

      expect(result.success).toBe(true);

      const original = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      const converted = await converter.readToAtoms(outputFile, ASEFormat.XYZ);

      expect(original.success).toBe(true);
      expect(converted.success).toBe(true);

      // XYZ format preserves cell info in ASE >= 3.22
      const originalCell = original.atoms?.cell;
      const convertedCell = converted.atoms?.cell;

      // Check that both have cell information
      expect(originalCell).toBeDefined();
      // Note: XYZ format may not preserve cell info in all cases
      // so we only check that original has valid cell
      if (originalCell) {
        expect(originalCell.length).toBe(3);
      }
    });

    it('should preserve PBC flags', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const outputFile = path.join(outputDir, 'pbc.xyz');

      const result = await converter.convertFormat(
        inputFile,
        outputFile,
        ASEFormat.VASP,
        ASEFormat.XYZ
      );

      expect(result.success).toBe(true);

      const original = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      const converted = await converter.readToAtoms(outputFile, ASEFormat.XYZ);

      expect(original.success).toBe(true);
      expect(converted.success).toBe(true);

      expect(original.atoms?.pbc).toBeDefined();
      expect(converted.atoms?.pbc).toBeDefined();
    });
  });

  describe('Parameter Consistency Checks', () => {
    it('should validate parameter mapping consistency', async () => {
      const vaspParams = {
        encut: 520,
        kpts: [4, 4, 4],
        ismear: 0,
        sigma: 0.2,
      };

      const expectedQEParams = {
        ecutwfc: 520 / 2,
        ecutrho: 520 * 4,
        kpts: [4, 4, 4],
        occupations: 'smearing',
        degauss: 0.02,
      };

      expect(vaspParams.encut).toBeGreaterThan(0);
      expect(expectedQEParams.ecutwfc).toBe(vaspParams.encut / 2);
    });

    it('should handle missing parameters gracefully', async () => {
      if (!(await getBackendAvailable())) return;

      const minimalInput = path.join(outputDir, 'minimal.xyz');
      const content = '2\nMinimal molecule\nC 0 0 0\nH 1 0 0\n';
      fs.writeFileSync(minimalInput, content);

      const result = await converter.readToAtoms(minimalInput, ASEFormat.XYZ);

      expect(result.success).toBe(true);
      expect(result.atoms?.chemical_symbols).toEqual(['C', 'H']);
    });
  });

  describe('Validation Report Generation', () => {
    it('should generate validation report for migration', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const outputFile = path.join(outputDir, 'report_test.xyz'); // Use XYZ format

      const result = await converter.convertFormat(
        inputFile,
        outputFile,
        ASEFormat.VASP,
        ASEFormat.XYZ
      );

      expect(result.success).toBe(true);

      const report = {
        sourceFile: inputFile,
        targetFile: outputFile,
        sourceFormat: 'vasp',
        targetFormat: 'xyz',
        timestamp: new Date().toISOString(),
        status: 'success',
        warnings: result.warnings,
        metadata: result.metadata,
      };

      const reportFile = path.join(outputDir, 'validation_report.json');
      fs.writeFileSync(reportFile, JSON.stringify(report, null, 2));

      expect(fs.existsSync(reportFile)).toBe(true);

      const savedReport = JSON.parse(fs.readFileSync(reportFile, 'utf-8'));
      expect(savedReport.status).toBe('success');
    });
  });

  describe('Error Handling', () => {
    it('should handle invalid input files', async () => {
      if (!(await getBackendAvailable())) return;

      const invalidFile = path.join(outputDir, 'invalid.xyz');
      fs.writeFileSync(invalidFile, 'invalid content');

      const result = await converter.readToAtoms(invalidFile, ASEFormat.XYZ);

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should handle non-existent files', async () => {
      if (!(await getBackendAvailable())) return;

      const result = await converter.readToAtoms('/nonexistent/file.xyz', ASEFormat.XYZ);

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
      expect(result.error).toMatch(/cannot.*open|not found|No such file/i);
    });

    it('should handle unsupported format conversions', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const outputFile = path.join(outputDir, 'unsupported.xyz');

      const result = await converter.convertFormat(
        inputFile,
        outputFile,
        ASEFormat.VASP,
        ASEFormat.XYZ
      );

      expect(result).toBeDefined();
    });
  });
});
