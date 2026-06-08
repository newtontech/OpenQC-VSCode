import {
  StructureExporter,
  ExportFormat,
  StructureData,
  ExportOptions,
} from '../../../src/utils/structureExporter';

describe('StructureExporter', () => {
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
    it('should return not-implemented result', async () => {
      const exporter = new StructureExporter();
      const structure: StructureData = {
        atoms: [
          { element: 'H', x: 0, y: 0, z: 0 },
          { element: 'O', x: 0.96, y: 0, z: 0 },
        ],
      };
      const options: ExportOptions = { format: 'xyz' };

      const result = await exporter.exportStructure(structure, options);
      expect(result.success).toBe(false);
      expect(result.error).toBe('Not implemented yet');
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
  });
});
