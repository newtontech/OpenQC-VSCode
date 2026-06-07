/**
 * Parameter Converter Tests
 */

import * as path from 'path';
import * as fs from 'fs';
import * as os from 'os';
import {
  ParameterConverter,
  ExtractionResult,
  ConversionResult,
} from '../../../src/utils/migration/parameterConverter';

describe('ParameterConverter', () => {
  let converter: ParameterConverter;
  const fixturesDir = path.join(__dirname, '../../fixtures');
  let tempDir: string;

  beforeEach(() => {
    converter = new ParameterConverter();
    // Create temp directory for test files
    tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-test-'));
  });

  afterEach(() => {
    // Clean up temp directory
    if (fs.existsSync(tempDir)) {
      fs.rmSync(tempDir, { recursive: true, force: true });
    }
  });

  describe('VASP Parameter Extraction', () => {
    it('should extract parameters from VASP INCAR', () => {
      const incarPath = path.join(fixturesDir, 'vasp', 'INCAR');
      const result = converter.extractParameters(incarPath);

      expect(result.success).toBe(true);
      expect(result.data).toBeDefined();
      expect(result.data?.format).toBe('vasp');
      expect(result.data?.parameters).toHaveProperty('ENCUT');
      expect(result.data?.parameters.ENCUT).toBe(520);
      expect(result.data?.parameters).toHaveProperty('ISMEAR');
      expect(result.data?.parameters.ISMEAR).toBe(0);
    });

    it('should parse boolean values correctly', () => {
      const incarPath = path.join(fixturesDir, 'vasp', 'INCAR');
      const result = converter.extractParameters(incarPath);

      expect(result.success).toBe(true);
      expect(result.data?.parameters.LCHARG).toBe(true);
      expect(result.data?.parameters.LWAVE).toBe(false);
    });

    it('should handle comments in INCAR', () => {
      const incarPath = path.join(fixturesDir, 'vasp', 'INCAR-Si');
      if (fs.existsSync(incarPath)) {
        const result = converter.extractParameters(incarPath);
        expect(result.success).toBe(true);
        expect(result.data?.parameters).toHaveProperty('ENCUT');
      }
    });
  });

  describe('Quantum ESPRESSO Parameter Extraction', () => {
    it('should extract parameters from QE input', () => {
      const qePath = path.join(fixturesDir, 'qe', 'Si.scf.in');
      const result = converter.extractParameters(qePath);

      expect(result.success).toBe(true);
      expect(result.data).toBeDefined();
      expect(result.data?.format).toBe('qe');
      expect(Object.keys(result.data?.parameters || {})).toContain('SYSTEM.ecutwfc');
      expect(result.data?.parameters['SYSTEM.ecutwfc']).toBe(40.0);
    });

    it('should parse namelist sections correctly', () => {
      const qePath = path.join(fixturesDir, 'qe', 'Si.scf.in');
      const result = converter.extractParameters(qePath);

      expect(result.success).toBe(true);
      expect(Object.keys(result.data?.parameters || {})).toContain('CONTROL.calculation');
      expect(Object.keys(result.data?.parameters || {})).toContain('ELECTRONS.conv_thr');
    });
  });

  describe('CP2K Parameter Extraction', () => {
    it('should extract parameters from CP2K input', () => {
      const cp2kPath = path.join(fixturesDir, 'cp2k', 'H2O.inp');
      const result = converter.extractParameters(cp2kPath);

      expect(result.success).toBe(true);
      expect(result.data).toBeDefined();
      expect(result.data?.format).toBe('cp2k');
    });
  });

  describe('Gaussian Parameter Extraction', () => {
    it('should extract parameters from Gaussian input', () => {
      const gaussianPath = path.join(fixturesDir, 'gaussian', 'H2O.gjf');
      if (fs.existsSync(gaussianPath)) {
        const result = converter.extractParameters(gaussianPath);

        expect(result.success).toBe(true);
        expect(result.data).toBeDefined();
        expect(result.data?.format).toBe('gaussian');
      }
    });
  });

  describe('Parameter Conversion', () => {
    it('should convert VASP parameters to QE format', async () => {
      const incarPath = path.join(fixturesDir, 'vasp', 'INCAR');
      const result = await converter.convertFile(incarPath, 'qe');

      expect(result.success).toBe(true);
      expect(result.source).toBeDefined();
      expect(result.target).toBeDefined();
      expect(result.target?.format).toBe('qe');
      expect(Object.keys(result.target?.parameters || {})).toContain('ecutwfc');
    });

    it('should convert VASP parameters to CP2K format', async () => {
      const incarPath = path.join(fixturesDir, 'vasp', 'INCAR');
      const result = await converter.convertFile(incarPath, 'cp2k');

      expect(result.success).toBe(true);
      expect(result.target?.format).toBe('cp2k');
      expect(Object.keys(result.target?.parameters || {})).toContain('CUTOFF');
    });

    it('should convert QE parameters to VASP format', async () => {
      const qePath = path.join(fixturesDir, 'qe', 'Si.scf.in');
      const result = await converter.convertFile(qePath, 'vasp');

      expect(result.success).toBe(true);
      expect(result.target?.format).toBe('vasp');
      expect(Object.keys(result.target?.parameters || {})).toContain('ENCUT');
    });

    it('should handle functional conversion from VASP to QE', async () => {
      // Create test INCAR in temp directory with proper name
      const testIncar = path.join(tempDir, 'INCAR');
      fs.writeFileSync(testIncar, 'GGA = PE\nENCUT = 500\n');

      const result = await converter.convertFile(testIncar, 'qe');

      expect(result.success).toBe(true);
      expect(Object.keys(result.target?.parameters || {})).toContain('input_dft');
      expect(result.target?.parameters.input_dft).toBe('PBE');
    });

    it('should warn about unmapped parameters', async () => {
      const incarPath = path.join(fixturesDir, 'vasp', 'INCAR');
      const result = await converter.convertFile(incarPath, 'qe');

      expect(result.success).toBe(true);
      expect(result.target?.unmapped.length).toBeGreaterThanOrEqual(0);
      expect(result.warnings.length).toBeGreaterThan(0);
    });
  });

  describe('Error Handling', () => {
    it('should handle missing files', () => {
      const result = converter.extractParameters('/nonexistent/file.incar');

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should handle unsupported formats', () => {
      const result = converter.extractParameters('/some/file.xyz');

      expect(result.success).toBe(false);
      expect(result.error).toContain('Unsupported format');
    });
  });

  describe('Mapping Suggestions', () => {
    it('should provide parameter mapping suggestions', () => {
      const mappings = converter.getMappingSuggestions('vasp', 'qe');

      expect(mappings.length).toBeGreaterThan(0);
      expect(mappings.some(m => m.sourceParam === 'ENCUT')).toBe(true);
    });

    it('should return empty array for unsupported pairs', () => {
      const mappings = converter.getMappingSuggestions('unknown', 'also_unknown');

      expect(mappings.length).toBe(0);
    });
  });
});
