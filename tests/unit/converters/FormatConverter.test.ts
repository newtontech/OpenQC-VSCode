/**
 * FormatConverter Unit Tests
 *
 * Tests for the format conversion functionality.
 * Following TDD approach with tests written before implementation.
 */

jest.setTimeout(15000);

// Mock util.promisify - must be before imports
const mockExecAsync = jest.fn();
jest.mock('util', () => {
  const actualUtil = jest.requireActual('util');
  return {
    ...actualUtil,
    promisify: jest.fn(() => mockExecAsync),
  };
});

import { FormatConverter, SupportedFormat, type ConversionResult } from '../../../src/converters';

// Mock vscode module
jest.mock(
  'vscode',
  () => ({
    window: {
      createOutputChannel: jest.fn(() => ({
        appendLine: jest.fn(),
        dispose: jest.fn(),
      })),
      showInformationMessage: jest.fn(),
      showErrorMessage: jest.fn(),
      activeTextEditor: {
        document: {
          uri: { fsPath: '/test/input.xyz' },
        },
      },
    },
    workspace: {
      openTextDocument: jest.fn(),
      fs: {
        readFile: jest.fn(),
      },
    },
    env: {
      clipboard: {
        writeText: jest.fn(),
      },
    },
    Uri: {
      file: jest.fn(path => ({ fsPath: path })),
    },
  }),
  { virtual: true }
);

// Mock child_process
jest.mock('child_process', () => ({
  exec: jest.fn(),
  spawn: jest.fn(),
}));

describe('FormatConverter', () => {
  let converter: FormatConverter;

  beforeEach(() => {
    converter = new FormatConverter();
    jest.clearAllMocks();
  });

  afterEach(() => {
    converter.dispose();
  });

  describe('Format Detection', () => {
    describe('detectFormat', () => {
      it('should detect VASP POSCAR files by name', () => {
        const result = FormatConverter.detectFormat('/path/to/POSCAR');
        expect(result.format).toBe(SupportedFormat.POSCAR);
        expect(result.confidence).toBe(1.0);
      });

      it('should detect VASP CONTCAR files by name', () => {
        const result = FormatConverter.detectFormat('/path/to/CONTCAR');
        expect(result.format).toBe(SupportedFormat.CONTCAR);
        expect(result.confidence).toBe(1.0);
      });

      it('should detect VASP INCAR files by name', () => {
        const result = FormatConverter.detectFormat('/path/to/INCAR');
        expect(result.format).toBe(SupportedFormat.VASP);
        expect(result.confidence).toBe(1.0);
      });

      it('should detect VASP KPOINTS files by name', () => {
        const result = FormatConverter.detectFormat('/path/to/KPOINTS');
        expect(result.format).toBe(SupportedFormat.VASP);
        expect(result.confidence).toBe(1.0);
      });

      it('should detect VASP POTCAR files by name', () => {
        const result = FormatConverter.detectFormat('/path/to/POTCAR');
        expect(result.format).toBe(SupportedFormat.VASP);
        expect(result.confidence).toBe(1.0);
      });

      it('should detect XYZ files by extension', () => {
        const result = FormatConverter.detectFormat('/path/to/molecule.xyz');
        expect(result.format).toBe(SupportedFormat.XYZ);
        expect(result.confidence).toBe(0.9);
      });

      it('should detect PDB files by extension', () => {
        const result = FormatConverter.detectFormat('/path/to/structure.pdb');
        expect(result.format).toBe(SupportedFormat.PDB);
        expect(result.confidence).toBe(0.9);
      });

      it('should detect CIF files by extension', () => {
        const result = FormatConverter.detectFormat('/path/to/crystal.cif');
        expect(result.format).toBe(SupportedFormat.CIF);
        expect(result.confidence).toBe(0.9);
      });

      it('should detect Gaussian files by .gjf extension', () => {
        const result = FormatConverter.detectFormat('/path/to/input.gjf');
        expect(result.format).toBe(SupportedFormat.GJF);
        expect(result.confidence).toBe(0.9);
      });

      it('should detect Gaussian files by .com extension', () => {
        const result = FormatConverter.detectFormat('/path/to/input.com');
        expect(result.format).toBe(SupportedFormat.COM);
        expect(result.confidence).toBe(0.9);
      });

      it('should detect ORCA files by .inp extension', () => {
        const result = FormatConverter.detectFormat('/path/to/orca.inp');
        expect(result.format).toBe(SupportedFormat.INP);
        expect(result.confidence).toBe(0.9);
      });

      it('should default to POSCAR for unknown formats', () => {
        const result = FormatConverter.detectFormat('/path/to/unknown.ext');
        expect(result.format).toBe(SupportedFormat.POSCAR);
        expect(result.confidence).toBe(0.5);
      });
    });
  });

  describe('Format Extensions', () => {
    describe('getExtensionForFormat', () => {
      it('should return correct extension for XYZ format', () => {
        expect(converter.getExtensionForFormat(SupportedFormat.XYZ)).toBe('xyz');
      });

      it('should return correct extension for PDB format', () => {
        expect(converter.getExtensionForFormat(SupportedFormat.PDB)).toBe('pdb');
      });

      it('should return correct extension for CIF format', () => {
        expect(converter.getExtensionForFormat(SupportedFormat.CIF)).toBe('cif');
      });

      it('should return POSCAR for VASP format', () => {
        expect(converter.getExtensionForFormat(SupportedFormat.VASP)).toBe('POSCAR');
      });

      it('should return POSCAR for POSCAR format', () => {
        expect(converter.getExtensionForFormat(SupportedFormat.POSCAR)).toBe('POSCAR');
      });

      it('should return gjf for Gaussian format', () => {
        expect(converter.getExtensionForFormat(SupportedFormat.Gaussian)).toBe('gjf');
      });

      it('should return inp for ORCA format', () => {
        expect(converter.getExtensionForFormat(SupportedFormat.ORCA)).toBe('inp');
      });
    });
  });

  describe('XYZ Conversion Utility', () => {
    describe('convertToXYZ', () => {
      it('should convert atoms array to XYZ format string', () => {
        const atoms = [
          { elem: 'C', x: 0.0, y: 0.0, z: 0.0 },
          { elem: 'H', x: 1.0, y: 0.0, z: 0.0 },
          { elem: 'H', x: -0.5, y: 0.866, z: 0.0 },
        ];

        const result = FormatConverter.convertToXYZ(atoms, 'Test Molecule');

        const lines = result.split('\n');
        expect(lines[0]).toBe('3');
        expect(lines[1]).toBe('Test Molecule');
        expect(lines[2]).toMatch(/C\s+0\.000000\s+0\.000000\s+0\.000000/);
        expect(lines[3]).toMatch(/H\s+1\.000000\s+0\.000000\s+0\.000000/);
        expect(lines[4]).toMatch(/H\s+-0\.500000\s+0\.866000\s+0\.000000/);
      });

      it('should use default comment when none provided', () => {
        const atoms = [{ elem: 'He', x: 0, y: 0, z: 0 }];
        const result = FormatConverter.convertToXYZ(atoms);
        const lines = result.split('\n');
        expect(lines[1]).toBe('molecule');
      });

      it('should format coordinates with 6 decimal places', () => {
        const atoms = [{ elem: 'O', x: 1.23456789, y: -9.87654321, z: 0.123456789 }];
        const result = FormatConverter.convertToXYZ(atoms);
        const lines = result.split('\n');
        expect(lines[2]).toMatch(/O\s+1\.234568\s+-9\.876543\s+0\.123457/);
      });

      it('should handle empty atom array', () => {
        const result = FormatConverter.convertToXYZ([], 'Empty');
        const lines = result.split('\n');
        expect(lines[0]).toBe('0');
        expect(lines[1]).toBe('Empty');
        expect(lines.length).toBe(2);
      });
    });
  });

  describe('Backend Integration', () => {
    describe('checkBackend', () => {
      beforeEach(() => {
        mockExecAsync.mockClear();
      });

      it('should return true when Python and dpdata are available', async () => {
        mockExecAsync.mockResolvedValueOnce({ stdout: 'Python 3.9.0', stderr: '' });
        mockExecAsync.mockResolvedValueOnce({ stdout: '', stderr: '' });

        const result = await converter.checkBackend();
        expect(result).toBe(true);
      });

      it('should return false when Python is not available', async () => {
        mockExecAsync.mockRejectedValueOnce(new Error('Command not found'));

        const result = await converter.checkBackend();
        expect(result).toBe(false);
      });

      it('should return false when dpdata is not installed', async () => {
        mockExecAsync.mockResolvedValueOnce({ stdout: 'Python 3.9.0', stderr: '' });
        mockExecAsync.mockRejectedValueOnce(
          new Error('ModuleNotFoundError: No module named dpdata')
        );

        const result = await converter.checkBackend();
        expect(result).toBe(false);
      });
    });
  });

  describe('Single File Conversion', () => {
    describe('convert', () => {
      beforeEach(() => {
        mockExecAsync.mockClear();
      });

      it('should successfully convert VASP to XYZ', async () => {
        const mockResult: ConversionResult = {
          success: true,
          input_format: 'vasp',
          output_format: 'xyz',
          atoms_count: 10,
          frames_count: 1,
          output_file: '/test/output.xyz',
        };

        mockExecAsync.mockResolvedValueOnce({
          stdout: JSON.stringify(mockResult),
          stderr: '',
        });

        const result = await converter.convert('/test/POSCAR', '/test/output.xyz');

        expect(result).toEqual(mockResult);
      });

      it('should handle conversion errors gracefully', async () => {
        mockExecAsync.mockRejectedValueOnce(new Error('Conversion failed'));

        const result = await converter.convert('/test/POSCAR', '/test/output.xyz');

        expect(result.success).toBe(false);
        expect(result.error).toBeDefined();
        expect(result.error_type).toBe('ExecutionError');
      });

      it('should pass from_format parameter when specified', async () => {
        const mockResult: ConversionResult = { success: true };
        mockExecAsync.mockResolvedValueOnce({ stdout: JSON.stringify(mockResult), stderr: '' });

        const result = await converter.convert(
          '/test/input',
          '/test/output',
          SupportedFormat.VASP,
          SupportedFormat.XYZ
        );

        expect(result.success).toBe(true);
      });

      it('should pass no-metadata flag when preserveMetadata is false', async () => {
        const converterNoMeta = new FormatConverter({ preserveMetadata: false });
        const mockResult: ConversionResult = { success: true };
        mockExecAsync.mockResolvedValueOnce({ stdout: JSON.stringify(mockResult), stderr: '' });

        const result = await converterNoMeta.convert('/test/input', '/test/output');

        expect(result.success).toBe(true);
      });
    });
  });

  describe('Batch Conversion', () => {
    describe('batchConvert', () => {
      beforeEach(() => {
        mockExecAsync.mockClear();
      });

      it('should successfully convert multiple files', async () => {
        const mockResult = {
          success: true,
          total: 3,
          successful: 3,
          failed: 0,
          results: [],
        };

        mockExecAsync.mockResolvedValueOnce({
          stdout: JSON.stringify(mockResult),
          stderr: '',
        });

        const result = await converter.batchConvert(
          ['/test/file1.xyz', '/test/file2.xyz', '/test/file3.xyz'],
          '/test/output',
          SupportedFormat.VASP
        );

        expect(result.success).toBe(true);
        expect(result.total).toBe(3);
        expect(result.successful).toBe(3);
      });

      it('should handle partial success in batch conversion', async () => {
        const mockResult = {
          success: false,
          total: 3,
          successful: 2,
          failed: 1,
          results: [],
        };

        mockExecAsync.mockResolvedValueOnce({
          stdout: JSON.stringify(mockResult),
          stderr: '',
        });

        const result = await converter.batchConvert(
          ['/test/file1.xyz', '/test/file2.xyz', '/test/file3.xyz'],
          '/test/output',
          SupportedFormat.VASP
        );

        expect(result.success).toBe(false);
        expect(result.successful).toBe(2);
        expect(result.failed).toBe(1);
      });

      it('should pass from_format parameter when specified', async () => {
        const mockResult = { success: true, total: 1, successful: 1, failed: 0, results: [] };
        mockExecAsync.mockResolvedValueOnce({ stdout: JSON.stringify(mockResult), stderr: '' });

        const result = await converter.batchConvert(
          ['/test/input.xyz'],
          '/test/output',
          SupportedFormat.VASP,
          SupportedFormat.XYZ
        );

        expect(result.total).toBe(1);
      });
    });
  });

  describe('Current Document Conversion', () => {
    describe('convertCurrentDocument', () => {
      beforeEach(() => {
        mockExecAsync.mockClear();
      });

      it('should convert active document when available', async () => {
        const mockResult: ConversionResult = {
          success: true,
          output_file: '/test/input.xyz',
        };
        mockExecAsync.mockResolvedValueOnce({ stdout: JSON.stringify(mockResult), stderr: '' });

        const result = await converter.convertCurrentDocument(SupportedFormat.XYZ);

        expect(result).toBeDefined();
        expect(result?.success).toBe(true);
      });
    });
  });
});
