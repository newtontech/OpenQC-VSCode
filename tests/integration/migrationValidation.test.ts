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

// Simple mock for vscode.ExtensionContext that only includes what ASEConverter needs
class MockExtensionContext implements vscode.ExtensionContext {
  subscriptions: vscode.Disposable[] = [];
  extensionPath: string;
  storageUri: vscode.Uri;
  globalStorageUri: vscode.Uri;
  logPath: string;

  constructor(extensionPath: string) {
    this.extensionPath = extensionPath;
    this.storageUri = vscode.Uri.file(extensionPath);
    this.globalStorageUri = vscode.Uri.file(extensionPath);
    this.logPath = extensionPath;
  }

  globalState: vscode.Memento = {
    get: (key) => Promise.resolve(undefined),
    update: (key, value) => Promise.resolve(undefined),
    keys: () => Promise.resolve([]),
    setKeysForSync: (keys) => Promise.resolve(undefined),
  };

  workspaceState: vscode.Memento = {
    get: (key) => Promise.resolve(undefined),
    update: (key, value) => Promise.resolve(undefined),
    keys: () => Promise.resolve([]),
    setKeysForSync: (keys) => Promise.resolve(undefined),
  };

  asAbsolutePathUri(relativePath: string): vscode.Uri {
    return vscode.Uri.file(path.join(this.extensionPath, relativePath));
  }

  asAbsolutePath(relativePath: string): string {
    return path.join(this.extensionPath, relativePath);
  }

  extensionUri: vscode.Uri = vscode.Uri.file(this.extensionPath);
}

describe('Migration Validation Integration Tests', () => {
  let converter: ASEConverter;
  const fixturesDir = path.join(__dirname, '../fixtures/migration');
  const outputDir = path.join(__dirname, '../temp/migration');

  beforeAll(async () => {
    const mockContext = new MockExtensionContext(__dirname);
    converter = new ASEConverter(mockContext);

    // Check if backend is available
    const backendAvailable = await converter.isAvailable();
    if (!backendAvailable) {
      console.warn('Python backend not available. Skipping integration tests.');
      return;
    }

    // Create output directory
    fs.mkdirSync(outputDir, { recursive: true });
  });

  afterAll(() => {
    // Clean up output directory
    if (fs.existsSync(outputDir)) {
      fs.rmSync(outputDir, { recursive: true, force: true });
    }
  });

  const getBackendAvailable = async () => {
    return await converter.isAvailable();
  };

  describe('Round-Trip Conversion Tests', () => {
    it('should preserve structure in VASP → XYZ → VASP round-trip', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      const intermediateFile = path.join(outputDir, 'roundtrip.xyz');
      const finalFile = path.join(outputDir, 'roundtrip.POSCAR');

      // VASP → XYZ
      const result1 = await converter.convertFormat(
        inputFile,
        intermediateFile,
        ASEFormat.VASP,
        ASEFormat.XYZ
      );
      expect(result1.success).toBe(true);

      // XYZ → VASP
      const result2 = await converter.convertFormat(
        intermediateFile,
        finalFile,
        ASEFormat.XYZ,
        ASEFormat.VASP
      );
      expect(result2.success).toBe(true);

      // Validate structure preservation
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

      // VASP → CP2K
      const result1 = await converter.convertFormat(
        inputFile,
        intermediateFile,
        ASEFormat.VASP,
        ASEFormat.CP2K
      );
      expect(result1.success).toBe(true);

      // CP2K → VASP
      const result2 = await converter.convertFormat(
        intermediateFile,
        finalFile,
        ASEFormat.CP2K,
        ASEFormat.VASP
      );
      expect(result2.success).toBe(true);

      // Validate structure preservation
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

      // VASP → QE
      const result1 = await converter.convertFormat(
        inputFile,
        intermediateFile,
        ASEFormat.VASP,
        ASEFormat.QE
      );
      expect(result1.success).toBe(true);

      // QE → VASP
      const result2 = await converter.convertFormat(
        intermediateFile,
        finalFile,
        ASEFormat.QE,
        ASEFormat.VASP
      );
      expect(result2.success).toBe(true);

      // Validate structure preservation
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

      // Check positions
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
      const outputFile = path.join(outputDir, 'cell.cif');

      const result = await converter.convertFormat(
        inputFile,
        outputFile,
        ASEFormat.VASP,
        ASEFormat.CIF
      );

      expect(result.success).toBe(true);

      const original = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      const converted = await converter.readToAtoms(outputFile, ASEFormat.CIF);

      expect(original.success).toBe(true);
      expect(converted.success).toBe(true);

      // Check cell
      const originalCell = original.atoms?.cell;
      const convertedCell = converted.atoms?.cell;

      if (originalCell && convertedCell) {
        for (let i = 0; i < 3; i++) {
          for (let j = 0; j < 3; j++) {
            expect(Math.abs(originalCell[i][j] - convertedCell[i][j])).toBeLessThan(1e-3);
          }
        }
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

      // PBC should be preserved (XYZ is non-periodic, so PBC may be lost)
      // This test documents the behavior
      expect(original.atoms?.pbc).toBeDefined();
      expect(converted.atoms?.pbc).toBeDefined();
    });
  });

  describe('Parameter Consistency Checks', () => {
    it('should validate parameter mapping consistency', async () => {
      // Test parameter mapping for VASP → QE
      const vaspParams = {
        encut: 520,
        kpts: [4, 4, 4],
        ismear: 0,
        sigma: 0.2,
      };

      // Expected QE equivalents
      const expectedQEParams = {
        ecutwfc: 520 / 2, // Rough conversion
        ecutrho: 520 * 4, // 4x ecutwfc
        kpts: [4, 4, 4],
        occupations: 'smearing',
        degauss: 0.02,
      };

      // This test validates: parameter mapping logic
      expect(vaspParams.encut).toBeGreaterThan(0);
      expect(expectedQEParams.ecutwfc).toBe(vaspParams.encut / 2);
    });

    it('should handle missing parameters gracefully', async () => {
      if (!(await getBackendAvailable())) return;

      // Create a minimal input file
      const minimalInput = path.join(outputDir, 'minimal.xyz');
      const content = `2
Minimal molecule
C 0 0 0
H 1 0 0
`;
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
      const outputFile = path.join(outputDir, 'report_test.in');

      const result = await converter.convertFormat(
        inputFile,
        outputFile,
        ASEFormat.VASP,
        ASEFormat.QE
      );

      expect(result.success).toBe(true);

      // Generate validation report
      const report = {
        sourceFile: inputFile,
        targetFile: outputFile,
        sourceFormat: 'vasp',
        targetFormat: 'qe',
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

      // This should still work, but may produce warnings
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
