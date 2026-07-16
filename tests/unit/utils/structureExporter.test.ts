import {
  StructureExporter,
  ExportFormat,
  StructureData,
  ExportOptions,
} from '../../../src/utils/structureExporter';
import * as vscode from 'vscode';

describe('StructureExporter', () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  describe('constructor', () => {
    it('should use default pythonPath from config', () => {
      const exporter = new StructureExporter();
      expect(exporter).toBeDefined();
    });

    it('should use context.extensionPath when context is provided', () => {
      const mockContext = {
        extensionPath: '/fake/extension/path',
      } as any;
      const exporter = new StructureExporter(mockContext);
      expect(exporter).toBeDefined();
    });
  });

  describe('exportStructure', () => {
    it('should reject empty structures', async () => {
      const exporter = new StructureExporter();

      const result = await exporter.exportStructure({ atoms: [] }, { format: 'xyz' });

      expect(result).toEqual({ success: false, error: 'No atoms to export' });
      expect(vscode.workspace.fs.writeFile).not.toHaveBeenCalled();
    });

    it('should export XYZ without requiring the Python backend', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [
          { element: 'H', x: 0, y: 0, z: 0 },
          { element: 'O', x: 0.96, y: 0, z: 0 },
        ],
      };
      const options: ExportOptions = { format: 'xyz', outputPath: '/tmp/water.xyz' };

      const result = await exporter.exportStructure(structure, options);
      expect(result.success).toBe(true);
      expect(result.outputPath).toBe('/tmp/water.xyz');
      expect(result.content).toContain('2\nOpenQC export');
      expect(result.content).toContain('H 0.000000 0.000000 0.000000');
      expect(result.content).toContain('O 0.960000 0.000000 0.000000');
      expect(vscode.workspace.fs.writeFile).toHaveBeenCalledWith(
        expect.objectContaining({ fsPath: '/tmp/water.xyz' }),
        expect.any(Buffer)
      );
    });

    it('should prompt for an output path when none is provided', async () => {
      const exporter = new StructureExporter();
      (vscode.window.showSaveDialog as jest.Mock).mockResolvedValue({
        fsPath: '/tmp/prompted.xyz',
      });

      const result = await exporter.exportStructure(
        { atoms: [{ element: 'He', x: 0, y: 0, z: 0 }] },
        { format: 'xyz' }
      );

      expect(result.success).toBe(true);
      expect(result.outputPath).toBe('/tmp/prompted.xyz');
      expect(vscode.window.showSaveDialog).toHaveBeenCalledWith(
        expect.objectContaining({
          saveLabel: 'Export Structure',
          filters: { XYZ: ['xyz'] },
        })
      );
    });

    it('should return a cancelled result when the save dialog is dismissed', async () => {
      const exporter = new StructureExporter();
      (vscode.window.showSaveDialog as jest.Mock).mockResolvedValue(undefined);

      const result = await exporter.exportStructure(
        { atoms: [{ element: 'He', x: 0, y: 0, z: 0 }] },
        { format: 'xyz' }
      );

      expect(result).toEqual({ success: false, error: 'Export cancelled' });
      expect(vscode.workspace.fs.writeFile).not.toHaveBeenCalled();
    });

    it('should export extended XYZ with lattice metadata from lattice parameters', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [{ element: 'C', x: Number.POSITIVE_INFINITY, y: 0.5, z: 1 }],
        cell: { a: 10, b: 20, c: 30, alpha: 90, beta: 90, gamma: 90 },
      };

      const result = await exporter.exportStructure(structure, {
        format: 'extxyz',
        outputPath: '/tmp/carbon.extxyz',
      });

      expect(result.success).toBe(true);
      expect(result.content).toContain(
        'Lattice="10.000000 0.000000 0.000000 0.000000 20.000000 0.000000 0.000000 0.000000 30.000000"'
      );
      expect(result.content).toContain('C 0.000000 0.500000 1.000000');
    });

    it('should export PDB records', async () => {
      const exporter = new StructureExporter();

      const result = await exporter.exportStructure(
        { atoms: [{ element: 'H', x: 1.2, y: 3.4, z: 5.6 }] },
        { format: 'pdb', outputPath: '/tmp/water.pdb' }
      );

      expect(result.success).toBe(true);
      expect(result.content).toContain('HETATM');
      expect(result.content).toContain('   1.200   3.400   5.600');
      expect(result.content).toContain('END');
    });

    it('should export PDB CONECT records for valid bonds', async () => {
      const exporter = new StructureExporter();

      const result = await exporter.exportStructure(
        {
          atoms: [
            { element: 'C', x: 0, y: 0, z: 0 },
            { element: 'O', x: 1.2, y: 0, z: 0 },
            { element: 'H', x: 2.2, y: 0, z: 0 },
          ],
          bonds: [
            { from: 0, to: 1, order: 2 },
            { from: 1, to: 2, order: 1 },
          ],
        },
        { format: 'pdb', outputPath: '/tmp/bonded.pdb' }
      );

      expect(result.success).toBe(true);
      expect(result.content).toContain('CONECT    1    2');
      expect(result.content).toContain('CONECT    2    3');
      expect(result.content?.trim().endsWith('END')).toBe(true);
    });

    it('should warn when non-PDB native exports drop edited bond topology', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [
          { element: 'C', x: 0, y: 0, z: 0 },
          { element: 'O', x: 1.2, y: 0, z: 0 },
        ],
        bonds: [{ from: 0, to: 1, order: 2 }],
      };

      const result = await exporter.exportStructure(structure, {
        format: 'xyz',
        outputPath: '/tmp/bonded.xyz',
      });

      expect(result.success).toBe(true);
      expect(result.warnings).toEqual([
        'XYZ export does not preserve edited bond topology or bond order; use PDB for basic CONECT connectivity.',
      ]);
    });

    it('should skip duplicate and out-of-range PDB CONECT bonds', async () => {
      const exporter = new StructureExporter();

      const result = await exporter.exportStructure(
        {
          atoms: [
            { element: 'C', x: 0, y: 0, z: 0 },
            { element: 'O', x: 1.2, y: 0, z: 0 },
          ],
          bonds: [
            { from: 0, to: 1, order: 1 },
            { from: 1, to: 0, order: 3 },
            { from: 0, to: 99, order: 1 },
          ],
        },
        { format: 'pdb', outputPath: '/tmp/deduped.pdb' }
      );

      expect(result.success).toBe(true);
      const conectLines = result.content?.split('\n').filter(line => line.startsWith('CONECT'));
      expect(conectLines).toEqual(['CONECT    1    2']);
      expect(result.warnings).toEqual([
        'Skipped PDB CONECT bond 2: atom indices 0-99 are out of range for 2 atoms',
      ]);
    });

    it('should export POSCAR when unit-cell vectors are available', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [
          { element: 'Si', x: 0, y: 0, z: 0 },
          { element: 'Si', x: 1.3575, y: 1.3575, z: 1.3575 },
        ],
        cell: {
          a: [2.715, 2.715, 0],
          b: [0, 2.715, 2.715],
          c: [2.715, 0, 2.715],
          pbc: [true, true, true],
          coordinateMode: 'cartesian',
        },
      };

      const result = await exporter.exportStructure(structure, {
        format: 'vasp',
        outputPath: '/tmp/POSCAR',
      });

      expect(result.success).toBe(true);
      expect(result.content).toContain('Si');
      expect(result.content).toContain('2');
      expect(result.content).toContain('Cartesian');
      expect(result.content).toContain('2.715000 2.715000 0.000000');
    });

    it('should preserve fractional coordinates as Direct in POSCAR', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [
          { element: 'Si', x: 0, y: 0, z: 0 },
          { element: 'Si', x: 0.25, y: 0.25, z: 0.25 },
        ],
        cell: {
          a: [2.715, 2.715, 0],
          b: [0, 2.715, 2.715],
          c: [2.715, 0, 2.715],
          pbc: [true, true, true],
          coordinateMode: 'fractional',
        },
      };

      const result = await exporter.exportStructure(structure, {
        format: 'vasp',
        outputPath: '/tmp/POSCAR',
      });

      expect(result.success).toBe(true);
      expect(result.content).toContain('Direct');
      expect(result.content).toContain('0.250000 0.250000 0.250000');
      expect(result.content).not.toContain('1.357500 1.357500 1.357500');
    });

    it('should preserve selective dynamics flags in POSCAR export and default missing flags to T', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [
          { element: 'Fe', x: 0, y: 0, z: 0, selectiveDynamics: [true, true, false] },
          { element: 'Fe', x: 0.5, y: 0.5, z: 0.5 },
        ] as any,
        cell: {
          a: [2.87, 0, 0],
          b: [0, 2.87, 0],
          c: [0, 0, 2.87],
          pbc: [true, true, true],
          coordinateMode: 'fractional',
        },
      };

      const result = await exporter.exportStructure(structure, {
        format: 'vasp',
        outputPath: '/tmp/POSCAR',
      });

      expect(result.success).toBe(true);
      const lines = result.content?.trim().split('\n') ?? [];
      expect(lines.slice(5, 8)).toEqual(['Fe', '2', 'Selective dynamics']);
      expect(lines[8]).toBe('Direct');
      expect(lines[9]).toBe('0.000000 0.000000 0.000000 T T F');
      expect(lines[10]).toBe('0.500000 0.500000 0.500000 T T T');
    });

    it('should export CIF cell parameters and atom sites', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [
          { element: 'Na', x: 0, y: 0, z: 0 },
          { element: 'Cl', x: 2.5, y: 2.5, z: 2.5 },
        ],
        cell: {
          a: [5, 0, 0],
          b: [0, 5, 0],
          c: [0, 0, 5],
          pbc: [true, true, true],
          coordinateMode: 'cartesian',
        },
      };

      const result = await exporter.exportStructure(structure, {
        format: 'cif',
        outputPath: '/tmp/nacl.cif',
      });

      expect(result.success).toBe(true);
      expect(result.content).toContain('_cell_length_a 5.000000');
      expect(result.content).toContain('_cell_angle_gamma 90.000000');
      expect(result.content).toContain('Na1 Na 0.000000 0.000000 0.000000');
      expect(result.content).toContain('Cl2 Cl 2.500000 2.500000 2.500000');
    });

    it('should convert fractional coordinates to Cartesian atom sites in CIF', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [{ element: 'Si', x: 0.25, y: 0.25, z: 0.25 }],
        cell: {
          a: [2.715, 2.715, 0],
          b: [0, 2.715, 2.715],
          c: [2.715, 0, 2.715],
          pbc: [true, true, true],
          coordinateMode: 'fractional',
        },
      };

      const result = await exporter.exportStructure(structure, {
        format: 'cif',
        outputPath: '/tmp/si.cif',
      });

      expect(result.success).toBe(true);
      expect(result.content).toContain('Si1 Si 1.357500 1.357500 1.357500');
    });

    it('should reject POSCAR and CIF exports without cell data', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      };

      const poscar = await exporter.exportStructure(structure, {
        format: 'vasp',
        outputPath: '/tmp/POSCAR',
      });
      const cif = await exporter.exportStructure(structure, {
        format: 'cif',
        outputPath: '/tmp/water.cif',
      });

      expect(poscar.success).toBe(false);
      expect(poscar.error).toContain('requires unit-cell vectors');
      expect(cif.success).toBe(false);
      expect(cif.error).toContain('requires unit-cell vectors');
    });

    it('should report unsupported native export formats clearly', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      };

      const result = await exporter.exportStructure(structure, {
        format: 'cp2k',
        outputPath: '/tmp/structure.inp',
      });

      expect(result.success).toBe(false);
      expect(result.error).toContain('requires the ASE backend');
      expect(vscode.workspace.fs.writeFile).not.toHaveBeenCalled();
    });

    it('should infer native writable formats from source filenames', () => {
      expect(StructureExporter.inferNativeFormatFromPath('/tmp/water.xyz')).toBe('xyz');
      expect(StructureExporter.inferNativeFormatFromPath('/tmp/water.extxyz')).toBe('extxyz');
      expect(StructureExporter.inferNativeFormatFromPath('/tmp/water.pdb')).toBe('pdb');
      expect(StructureExporter.inferNativeFormatFromPath('/tmp/water.cif')).toBe('cif');
      expect(StructureExporter.inferNativeFormatFromPath('/tmp/POSCAR')).toBe('vasp');
      expect(StructureExporter.inferNativeFormatFromPath('C:\\calc\\CONTCAR')).toBe('vasp');
      expect(StructureExporter.inferNativeFormatFromPath('/tmp/input.inp')).toBeUndefined();
    });

    it('should overwrite inferred native source formats without prompting', async () => {
      const exporter = new StructureExporter();
      const result = await exporter.overwriteStructureFile(
        { atoms: [{ element: 'Ne', x: 1, y: 2, z: 3 }] },
        '/tmp/neon.xyz'
      );

      expect(result.success).toBe(true);
      expect(result.outputPath).toBe('/tmp/neon.xyz');
      expect(result.content).toContain('Ne 1.000000 2.000000 3.000000');
      expect(vscode.window.showSaveDialog).not.toHaveBeenCalled();
      expect(vscode.workspace.fs.writeFile).toHaveBeenCalledWith(
        expect.objectContaining({ fsPath: '/tmp/neon.xyz' }),
        expect.any(Buffer)
      );
    });

    it('should reject source overwrite when the source format is not natively writable', async () => {
      const exporter = new StructureExporter();
      const result = await exporter.overwriteStructureFile(
        { atoms: [{ element: 'H', x: 0, y: 0, z: 0 }] },
        '/tmp/input.inp'
      );

      expect(result.success).toBe(false);
      expect(result.error).toContain('Cannot infer');
      expect(vscode.workspace.fs.writeFile).not.toHaveBeenCalled();
    });
  });

  describe('getFormatDisplayName', () => {
    it('should return correct display names for all formats', () => {
      const expectations: Record<ExportFormat, string> = {
        vasp: 'VASP POSCAR',
        cp2k: 'CP2K Input',
        qe: 'Quantum ESPRESSO',
        gaussian: 'Gaussian Input',
        orca: 'ORCA Input',
        nwchem: 'NWChem Input',
        gamess: 'GAMESS Input',
        lammps: 'LAMMPS Data',
        xyz: 'XYZ',
        extxyz: 'Extended XYZ',
        pdb: 'PDB',
        cif: 'CIF',
      };

      for (const [format, expectedName] of Object.entries(expectations)) {
        expect(StructureExporter.getFormatDisplayName(format as ExportFormat)).toBe(expectedName);
      }
    });
  });

  describe('getSupportedFormats', () => {
    it('should return all 12 supported export formats', () => {
      const formats = StructureExporter.getSupportedFormats();
      expect(formats).toHaveLength(12);
      expect(formats).toContain('vasp');
      expect(formats).toContain('cp2k');
      expect(formats).toContain('qe');
      expect(formats).toContain('gaussian');
      expect(formats).toContain('orca');
      expect(formats).toContain('nwchem');
      expect(formats).toContain('gamess');
      expect(formats).toContain('lammps');
      expect(formats).toContain('xyz');
      expect(formats).toContain('extxyz');
      expect(formats).toContain('pdb');
      expect(formats).toContain('cif');
    });

    it('should return formats as an array', () => {
      const formats = StructureExporter.getSupportedFormats();
      expect(Array.isArray(formats)).toBe(true);
    });

    it('should expose only directly implemented native export formats for command pickers', () => {
      expect(StructureExporter.getNativeFormats()).toEqual(['xyz', 'extxyz', 'pdb', 'vasp', 'cif']);
    });
  });
});
