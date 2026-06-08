/**
 * Unit tests for OpenQCStructure conversion utilities.
 */

import {
  openQCStructureToXYZ,
  xyzToOpenQCStructure,
  molecularStructureToOpenQCStructure,
  openQCStructureToMolecularStructure,
  poscarToOpenQCStructure,
  createOpenQCStructure,
} from '../../../src/structures/converters';
import { OPENQC_STRUCTURE_SCHEMA_VERSION } from '../../../src/structures/OpenQCStructure';
import {
  WATER_XYZ_CONTENT,
  METHANE_XYZ_CONTENT,
  SI_POSCAR_CONTENT,
  WATER_STRUCTURE,
  SILICON_STRUCTURE,
} from '../../../src/structures/fixtures';
import type { MolecularStructure } from '../../../src/visualizers/types';

// ============================================================================
// openQCStructureToXYZ
// ============================================================================

describe('openQCStructureToXYZ', () => {
  it('converts a water structure to XYZ format', () => {
    const xyz = openQCStructureToXYZ(WATER_STRUCTURE);

    const lines = xyz.trim().split('\n');
    expect(lines[0]).toBe('3');
    expect(lines[1]).toBe('Water');
    expect(lines.length).toBe(5); // 3 atoms + header + comment
    expect(lines[2].trim().split(/\s+/)[0]).toBe('O');
  });

  it('uses "OpenQCStructure" as default comment when name is undefined', () => {
    const noName = { ...WATER_STRUCTURE, name: undefined };
    const xyz = openQCStructureToXYZ(noName);
    expect(xyz.split('\n')[1]).toBe('OpenQCStructure');
  });

  it('formats coordinates to 6 decimal places', () => {
    const xyz = openQCStructureToXYZ(WATER_STRUCTURE);
    const atomLine = xyz.split('\n')[2].trim();
    const parts = atomLine.split(/\s+/);
    // x coordinate should have 6 decimal places
    expect(parts[1]).toMatch(/^\d+\.\d{6}$/);
  });
});

// ============================================================================
// xyzToOpenQCStructure
// ============================================================================

describe('xyzToOpenQCStructure', () => {
  it('parses a water XYZ string', () => {
    const structure = xyzToOpenQCStructure(WATER_XYZ_CONTENT);

    expect(structure.schemaVersion).toBe(OPENQC_STRUCTURE_SCHEMA_VERSION);
    expect(structure.kind).toBe('molecule');
    expect(structure.atoms).toHaveLength(3);
    expect(structure.atoms[0].element).toBe('O');
    expect(structure.atoms[1].element).toBe('H');
    expect(structure.atoms[2].element).toBe('H');
    expect(structure.name).toBe('Water molecule');
  });

  it('parses a methane XYZ string', () => {
    const structure = xyzToOpenQCStructure(METHANE_XYZ_CONTENT);

    expect(structure.atoms).toHaveLength(5);
    expect(structure.atoms[0].element).toBe('C');
  });

  it('throws on too-short content', () => {
    expect(() => xyzToOpenQCStructure('1\n')).toThrow('file too short');
  });

  it('throws on non-numeric atom count', () => {
    expect(() => xyzToOpenQCStructure('abc\ntest\nH 0 0 0\n')).toThrow('positive integer');
  });

  it('throws on zero atom count', () => {
    expect(() => xyzToOpenQCStructure('0\ntest\n')).toThrow('positive integer');
  });

  it('throws when no valid atoms found', () => {
    // Atom count says 1 but the line is malformed
    expect(() => xyzToOpenQCStructure('1\ntest\nbad line\n')).toThrow('no valid atom lines');
  });

  it('passes through source metadata', () => {
    const structure = xyzToOpenQCStructure(WATER_XYZ_CONTENT, {
      sourceSoftware: 'gaussian',
      name: 'Custom Name',
    });
    expect(structure.name).toBe('Custom Name');
    expect(structure.metadata?.source?.software).toBe('gaussian');
    expect(structure.metadata?.source?.parser).toBe('xyz');
  });

  it('round-trips through XYZ conversion', () => {
    const original = WATER_STRUCTURE;
    const xyz = openQCStructureToXYZ(original);
    const restored = xyzToOpenQCStructure(xyz);

    expect(restored.atoms).toHaveLength(original.atoms.length);
    for (let i = 0; i < original.atoms.length; i++) {
      expect(restored.atoms[i].element).toBe(original.atoms[i].element);
      expect(Math.abs(restored.atoms[i].x - original.atoms[i].x)).toBeLessThan(1e-5);
      expect(Math.abs(restored.atoms[i].y - original.atoms[i].y)).toBeLessThan(1e-5);
      expect(Math.abs(restored.atoms[i].z - original.atoms[i].z)).toBeLessThan(1e-5);
    }
  });
});

// ============================================================================
// molecularStructureToOpenQCStructure
// ============================================================================

describe('molecularStructureToOpenQCStructure', () => {
  it('converts a simple molecule-only MolecularStructure', () => {
    const ms: MolecularStructure = {
      atoms: [
        { element: 'H', x: 0, y: 0, z: 0 },
        { element: 'H', x: 0.74, y: 0, z: 0 },
      ],
      name: 'H2',
    };

    const structure = molecularStructureToOpenQCStructure(ms);

    expect(structure.kind).toBe('molecule');
    expect(structure.atoms).toHaveLength(2);
    expect(structure.atoms[0].element).toBe('H');
    expect(structure.name).toBe('H2');
    expect(structure.cell).toBeUndefined();
  });

  it('converts a periodic MolecularStructure with lattice', () => {
    const ms: MolecularStructure = {
      atoms: [
        { element: 'Si', x: 0, y: 0, z: 0 },
        { element: 'Si', x: 1.3575, y: 1.3575, z: 1.3575 },
      ],
      name: 'Si',
      lattice: {
        a: [5.43, 0, 0],
        b: [0, 5.43, 0],
        c: [0, 0, 5.43],
      },
    };

    const structure = molecularStructureToOpenQCStructure(ms);

    expect(structure.kind).toBe('periodic');
    expect(structure.cell).toBeDefined();
    expect(structure.cell?.a).toEqual([5.43, 0, 0]);
    expect(structure.cell?.pbc).toEqual([true, true, true]);
  });

  it('converts bonds correctly', () => {
    const ms: MolecularStructure = {
      atoms: [
        { element: 'C', x: 0, y: 0, z: 0 },
        { element: 'O', x: 1.13, y: 0, z: 0 },
      ],
      bonds: [{ atomIndex1: 0, atomIndex2: 1, length: 1.13, order: 2 }],
    };

    const structure = molecularStructureToOpenQCStructure(ms);

    expect(structure.bonds).toHaveLength(1);
    expect(structure.bonds![0].from).toBe(0);
    expect(structure.bonds![0].to).toBe(1);
    expect(structure.bonds![0].order).toBe(2);
  });

  it('passes source metadata through', () => {
    const ms: MolecularStructure = {
      atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
    };

    const structure = molecularStructureToOpenQCStructure(ms, {
      sourceSoftware: 'cp2k',
      sourceParser: 'cp2k-coord',
    });

    expect(structure.metadata?.source?.software).toBe('cp2k');
    expect(structure.metadata?.source?.parser).toBe('cp2k-coord');
  });
});

// ============================================================================
// openQCStructureToMolecularStructure
// ============================================================================

describe('openQCStructureToMolecularStructure', () => {
  it('converts a molecule OpenQCStructure back to MolecularStructure', () => {
    const ms = openQCStructureToMolecularStructure(WATER_STRUCTURE);

    expect(ms.atoms).toHaveLength(3);
    expect(ms.atoms[0].element).toBe('O');
    expect(ms.name).toBe('Water');
    expect(ms.lattice).toBeUndefined();
  });

  it('converts a periodic OpenQCStructure with cell', () => {
    const ms = openQCStructureToMolecularStructure(SILICON_STRUCTURE);

    expect(ms.lattice).toBeDefined();
    expect(ms.lattice?.a).toEqual([5.43, 0, 0]);
  });

  it('converts bonds back', () => {
    const structureWithBonds = {
      ...WATER_STRUCTURE,
      bonds: [
        { from: 0, to: 1, order: 1 },
        { from: 0, to: 2, order: 1 },
      ],
    };
    const ms = openQCStructureToMolecularStructure(structureWithBonds);

    expect(ms.bonds).toHaveLength(2);
    expect(ms.bonds![0].atomIndex1).toBe(0);
    expect(ms.bonds![0].atomIndex2).toBe(1);
  });

  it('round-trips molecule data through both converters', () => {
    const original: MolecularStructure = {
      atoms: [
        { element: 'N', x: 0, y: 0, z: 0 },
        { element: 'H', x: 0.94, y: 0, z: 0 },
        { element: 'H', x: -0.47, y: 0.81, z: 0 },
        { element: 'H', x: -0.47, y: -0.81, z: 0 },
      ],
      name: 'NH3',
    };

    const dto = molecularStructureToOpenQCStructure(original);
    const restored = openQCStructureToMolecularStructure(dto);

    expect(restored.atoms).toHaveLength(original.atoms.length);
    expect(restored.name).toBe('NH3');
    for (let i = 0; i < original.atoms.length; i++) {
      expect(restored.atoms[i].element).toBe(original.atoms[i].element);
    }
  });
});

// ============================================================================
// poscarToOpenQCStructure
// ============================================================================

describe('poscarToOpenQCStructure', () => {
  it('parses a silicon POSCAR', () => {
    const structure = poscarToOpenQCStructure(SI_POSCAR_CONTENT, 'POSCAR');

    expect(structure.kind).toBe('periodic');
    expect(structure.name).toBe('POSCAR');
    expect(structure.atoms).toHaveLength(2);
    expect(structure.atoms[0].element).toBe('Si');
    expect(structure.atoms[1].element).toBe('Si');
    expect(structure.cell).toBeDefined();
    expect(structure.cell?.a).toEqual([5.43, 0, 0]);
    expect(structure.cell?.coordinateMode).toBe('fractional');
    expect(structure.metadata?.source?.software).toBe('vasp');
  });

  it('throws on too-short POSCAR content', () => {
    expect(() => poscarToOpenQCStructure('short')).toThrow('too few lines');
  });

  it('handles selective dynamics flags', () => {
    const poscarWithSD = `Test SD
1.0
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
Si
1
Selective dynamics
0.0 0.0 0.0 T T F
`;

    const structure = poscarToOpenQCStructure(poscarWithSD);
    expect(structure.atoms[0].selectiveDynamics).toEqual([true, true, false]);
  });
});

// ============================================================================
// createOpenQCStructure
// ============================================================================

describe('createOpenQCStructure', () => {
  it('creates a minimal molecule structure', () => {
    const structure = createOpenQCStructure([
      { element: 'H', x: 0, y: 0, z: 0 },
      { element: 'H', x: 0.74, y: 0, z: 0 },
    ]);

    expect(structure.schemaVersion).toBe(OPENQC_STRUCTURE_SCHEMA_VERSION);
    expect(structure.kind).toBe('molecule');
    expect(structure.atoms).toHaveLength(2);
    expect(structure.cell).toBeUndefined();
  });

  it('creates a periodic structure when cell is provided', () => {
    const structure = createOpenQCStructure([{ element: 'Si', x: 0, y: 0, z: 0 }], {
      kind: 'periodic',
      cell: {
        a: [5.43, 0, 0],
        b: [0, 5.43, 0],
        c: [0, 0, 5.43],
        pbc: [true, true, true],
      },
    });

    expect(structure.kind).toBe('periodic');
    expect(structure.cell).toBeDefined();
  });

  it('infers periodic kind from cell presence when kind not specified', () => {
    const structure = createOpenQCStructure([{ element: 'Si', x: 0, y: 0, z: 0 }], {
      cell: {
        a: [5.43, 0, 0],
        b: [0, 5.43, 0],
        c: [0, 0, 5.43],
        pbc: [true, true, true],
      },
    });

    expect(structure.kind).toBe('periodic');
  });

  it('sets charge and multiplicity in metadata', () => {
    const structure = createOpenQCStructure([{ element: 'H', x: 0, y: 0, z: 0 }], {
      charge: 1,
      multiplicity: 2,
    });

    expect(structure.metadata?.charge).toBe(1);
    expect(structure.metadata?.multiplicity).toBe(2);
  });

  it('sets the name', () => {
    const structure = createOpenQCStructure([{ element: 'H', x: 0, y: 0, z: 0 }], {
      name: 'Test molecule',
    });

    expect(structure.name).toBe('Test molecule');
  });
});
