import * as vscode from 'vscode';
import { StructureConverter, Atom } from '../../../src/visualizers/StructureConverter';
import { Molecule3D } from '../../../src/visualizers/Molecule3D';
import { VASPParser } from '../../../src/parsers/VASPParser';

// Mock vscode module
jest.mock('vscode', () => ({
  Uri: {
    file: (path: string) => ({ fsPath: path, scheme: 'file' }),
  },
  window: {
    createWebviewPanel: jest.fn(),
    activeTextEditor: undefined,
    showInformationMessage: jest.fn(),
    showErrorMessage: jest.fn(),
  },
  ViewColumn: {
    One: 1,
    Two: 2,
  },
}));

describe('Visualization Flow Integration', () => {
  describe('VASP POSCAR to 3D Viewer', () => {
    const poscarContent = `Si Diamond Structure
1.0
3.84 0.00 0.00
0.00 3.84 0.00
0.00 0.00 3.84
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
`;

    it('parses POSCAR and extracts atoms', () => {
      const parser = new VASPParser(poscarContent, 'POSCAR');
      const atoms = parser.getCoordinates();
      const atomTypes = parser.getAtomTypes();

      expect(atoms).toHaveLength(2);
      expect(atomTypes).toEqual(['Si']);
      expect(atoms[0]).toEqual([0, 0, 0]);
      expect(atoms[1]).toEqual([0.25, 0.25, 0.25]);
    });

    it('converts VASP structure to atoms array', () => {
      const molecule3D = new Molecule3D();
      const atoms = molecule3D.parseAtoms(poscarContent, 'VASP');

      expect(atoms).toHaveLength(2);
      expect(atoms[0]).toEqual({ elem: 'Si', x: 0, y: 0, z: 0 });
      expect(atoms[1]).toEqual({ elem: 'Si', x: 0.25, y: 0.25, z: 0.25 });
    });

    it('converts atoms to XYZ format', () => {
      const molecule3D = new Molecule3D();
      const atoms = molecule3D.parseAtoms(poscarContent, 'VASP');
      const xyz = StructureConverter.atomsToXYZ(atoms, 'Si Diamond');

      const lines = xyz.split('\n');
      expect(lines[0]).toBe('2'); // Atom count
      expect(lines[1]).toBe('Si Diamond'); // Comment
      expect(lines[2]).toContain('Si');
      expect(lines[3]).toContain('Si');
    });
  });

  describe('Gaussian to 3D Viewer', () => {
    const gaussianContent = `#P B3LYP/6-31G(d) Opt

Water molecule optimization

0 1
O  -0.464   0.177   0.0
H   0.441  -0.143   0.0
H  -0.441  -0.143   0.9
`;

    it('converts Gaussian structure to atoms array', () => {
      const molecule3D = new Molecule3D();
      const atoms = molecule3D.parseAtoms(gaussianContent, 'Gaussian');

      expect(atoms).toHaveLength(3);
      expect(atoms[0].elem).toBe('O');
      expect(atoms[1].elem).toBe('H');
      expect(atoms[2].elem).toBe('H');
    });

    it('converts to XYZ format', () => {
      const molecule3D = new Molecule3D();
      const atoms = molecule3D.parseAtoms(gaussianContent, 'Gaussian');
      const xyz = StructureConverter.atomsToXYZ(atoms, 'Water');

      const lines = xyz.split('\n');
      expect(lines[0]).toBe('3');
      expect(lines[1]).toBe('Water');
      expect(lines[2]).toContain('O');
      expect(lines[2]).toContain('-0.464');
    });
  });

  describe('ORCA to 3D Viewer', () => {
    const orcaContent = `! B3LYP def2-SVP Opt

* xyz 0 1
C   0.0000   0.0000   0.0000
H   0.0000   0.0000   1.0890
H   1.0267   0.0000  -0.3630
H  -0.5134  -0.8893  -0.3630
H  -0.5134   0.8893  -0.3630
*
`;

    it('converts ORCA structure to atoms array', () => {
      const molecule3D = new Molecule3D();
      const atoms = molecule3D.parseAtoms(orcaContent, 'ORCA');

      expect(atoms).toHaveLength(5);
      expect(atoms[0].elem).toBe('C');
      expect(atoms[1].elem).toBe('H');
    });
  });

  describe('Edge Cases', () => {
    it('handles empty structure gracefully', () => {
      const xyz = StructureConverter.atomsToXYZ([], 'empty');
      const lines = xyz.split('\n');
      expect(lines[0]).toBe('0');
      expect(lines[1]).toBe('empty');
    });

    it('handles special characters in molecule name', () => {
      const atoms: Atom[] = [{ elem: 'Fe', x: 0, y: 0, z: 0 }];
      const xyz = StructureConverter.atomsToXYZ(atoms, 'Iron(III) oxide [Fe2O3]');
      expect(xyz).toContain('Iron(III) oxide [Fe2O3]');
    });

    it('handles very small coordinates', () => {
      const atoms: Atom[] = [{ elem: 'H', x: 0.0000001, y: -0.0000001, z: 0.0 }];
      const xyz = StructureConverter.atomsToXYZ(atoms, 'test');
      expect(xyz).toContain('0.000000');
    });
  });

  describe('JSON format conversion', () => {
    it('converts atoms to NGL-compatible JSON', () => {
      const atoms: Atom[] = [
        { elem: 'H', x: 0, y: 0, z: 0 },
        { elem: 'O', x: 0, y: 0, z: 1 },
      ];

      const json = StructureConverter.atomsToJSON(atoms);
      const parsed = JSON.parse(json);

      expect(parsed.atoms).toHaveLength(2);
      expect(parsed.atoms[0]).toEqual({ elem: 'H', x: 0, y: 0, z: 0 });
    });

    it('includes unit cell data for periodic systems', () => {
      const atoms: Atom[] = [{ elem: 'Si', x: 0, y: 0, z: 0 }];
      const unitCell = {
        a: 5.43,
        b: 5.43,
        c: 5.43,
        alpha: 90,
        beta: 90,
        gamma: 90,
      };

      const json = StructureConverter.atomsToJSON(atoms, unitCell);
      const parsed = JSON.parse(json);

      expect(parsed.unitCell).toEqual(unitCell);
    });
  });
});
