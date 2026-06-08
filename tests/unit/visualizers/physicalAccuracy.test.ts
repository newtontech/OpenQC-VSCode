/**
 * Physical accuracy validation for 3D molecular visualization
 *
 * Validates that the constants and calculations used in the visualization
 * are physically correct. These tests guard against regression of
 * scientific data (radii, colors, bond detection thresholds).
 *
 * Reference values are taken from established chemistry databases:
 * - Covalent radii: Cordero et al., Dalton Trans. 2008, 2832-2838
 * - Van der Waals radii: Bondi, J. Phys. Chem. 68, 441-452 (1964)
 * - CPK colors: Corey & Pauling / Koltun convention
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/19
 */

import { COVALENT_RADII, VDW_RADII, ELEMENT_COLORS } from '../../../src/visualizers/types';
import { StructureConverter } from '../../../src/visualizers/StructureConverter';
import { getElementSymbol } from '../../../src/visualizers/atomicData';

// ---------------------------------------------------------------------------
// 1. Covalent Radii
// ---------------------------------------------------------------------------

describe('Physical accuracy: covalent radii', () => {
  /**
   * Selected reference values from Cordero et al. (2008).
   * Tolerance: +/- 0.05 Angstroms to account for rounding.
   */
  const referenceCovalentRadii: Record<string, { value: number; tolerance: number }> = {
    H: { value: 0.31, tolerance: 0.05 },
    C: { value: 0.76, tolerance: 0.05 },
    N: { value: 0.71, tolerance: 0.05 },
    O: { value: 0.66, tolerance: 0.05 },
    F: { value: 0.57, tolerance: 0.05 },
    P: { value: 1.07, tolerance: 0.05 },
    S: { value: 1.05, tolerance: 0.05 },
    Cl: { value: 1.02, tolerance: 0.05 },
    Fe: { value: 1.32, tolerance: 0.08 },
    Cu: { value: 1.32, tolerance: 0.08 },
    Au: { value: 1.36, tolerance: 0.08 },
  };

  Object.entries(referenceCovalentRadii).forEach(([element, { value, tolerance }]) => {
    it(`covalent radius for ${element} should be ${value} +/- ${tolerance} A`, () => {
      const stored = COVALENT_RADII[element];
      expect(stored).toBeDefined();
      expect(Math.abs(stored - value)).toBeLessThanOrEqual(tolerance);
    });
  });

  it('should have a default fallback radius', () => {
    expect(COVALENT_RADII.default).toBeDefined();
    expect(COVALENT_RADII.default).toBeGreaterThan(0);
  });

  it('all covalent radii should be positive', () => {
    Object.entries(COVALENT_RADII).forEach(([key, radius]) => {
      expect(radius).toBeGreaterThan(0);
    });
  });
});

// ---------------------------------------------------------------------------
// 2. Van der Waals Radii
// ---------------------------------------------------------------------------

describe('Physical accuracy: van der Waals radii', () => {
  /**
   * Reference values from Bondi (1964) and Alvarez (2013).
   * Tolerance: +/- 0.1 Angstroms.
   */
  const referenceVdwRadii: Record<string, { value: number; tolerance: number }> = {
    H: { value: 1.2, tolerance: 0.1 },
    C: { value: 1.7, tolerance: 0.1 },
    N: { value: 1.55, tolerance: 0.1 },
    O: { value: 1.52, tolerance: 0.1 },
    F: { value: 1.47, tolerance: 0.1 },
    P: { value: 1.8, tolerance: 0.1 },
    S: { value: 1.8, tolerance: 0.1 },
    Cl: { value: 1.75, tolerance: 0.1 },
  };

  Object.entries(referenceVdwRadii).forEach(([element, { value, tolerance }]) => {
    it(`van der Waals radius for ${element} should be ${value} +/- ${tolerance} A`, () => {
      const stored = VDW_RADII[element];
      expect(stored).toBeDefined();
      expect(Math.abs(stored - value)).toBeLessThanOrEqual(tolerance);
    });
  });

  it('van der Waals radii should always be larger than covalent radii for the same element', () => {
    Object.keys(COVALENT_RADII).forEach(element => {
      if (element === 'default') {
        return;
      }
      const covalent = COVALENT_RADII[element];
      const vdw = VDW_RADII[element];
      if (vdw !== undefined) {
        expect(vdw).toBeGreaterThan(covalent);
      }
    });
  });
});

// ---------------------------------------------------------------------------
// 3. Element Colors (CPK convention)
// ---------------------------------------------------------------------------

describe('Physical accuracy: CPK element colors', () => {
  /**
   * CPK coloring convention reference values.
   * These are the standard colors used in molecular visualization.
   */
  const referenceColors: Record<string, string> = {
    H: '#FFFFFF', // White
    C: '#909090', // Dark gray
    N: '#3050F8', // Blue
    O: '#FF0D0D', // Red
    S: '#FFFF30', // Yellow
    Cl: '#1FF01F', // Green
    Fe: '#E06633', // Orange-brown
    Au: '#FFD123', // Gold
  };

  Object.entries(referenceColors).forEach(([element, expectedColor]) => {
    it(`${element} color should match CPK convention (${expectedColor})`, () => {
      expect(ELEMENT_COLORS[element]).toBe(expectedColor);
    });
  });

  it('should have a default fallback color', () => {
    expect(ELEMENT_COLORS.default).toBeDefined();
    expect(ELEMENT_COLORS.default).toMatch(/^#[0-9A-Fa-f]{6}$/);
  });

  it('all colors should be valid hex color codes', () => {
    Object.entries(ELEMENT_COLORS).forEach(([key, color]) => {
      expect(color).toMatch(/^#[0-9A-Fa-f]{6}$/);
    });
  });
});

// ---------------------------------------------------------------------------
// 4. Atomic Number Lookup
// ---------------------------------------------------------------------------

describe('Physical accuracy: atomic number lookup', () => {
  it('should return correct symbols for well-known elements', () => {
    expect(getElementSymbol(1)).toBe('H');
    expect(getElementSymbol(6)).toBe('C');
    expect(getElementSymbol(7)).toBe('N');
    expect(getElementSymbol(8)).toBe('O');
    expect(getElementSymbol(26)).toBe('Fe');
    expect(getElementSymbol(29)).toBe('Cu');
    expect(getElementSymbol(79)).toBe('Au');
  });

  it('should return unknown symbol for invalid atomic numbers', () => {
    expect(getElementSymbol(0)).toBe('X');
    expect(getElementSymbol(-1)).toBe('X');
    expect(getElementSymbol(1.5)).toBe('X');
  });
});

// ---------------------------------------------------------------------------
// 5. Bond Length Calculation
// ---------------------------------------------------------------------------

describe('Physical accuracy: bond length validation', () => {
  /**
   * Reference bond lengths from experimental data.
   * The visualization should detect bonds within 120% of the sum of covalent radii.
   */
  const referenceBondLengths: Array<{
    element1: string;
    element2: string;
    expectedMin: number;
    expectedMax: number;
    description: string;
  }> = [
    {
      element1: 'C',
      element2: 'C',
      expectedMin: 1.2,
      expectedMax: 1.8,
      description: 'C-C single bond (~1.54 A)',
    },
    {
      element1: 'C',
      element2: 'H',
      expectedMin: 0.9,
      expectedMax: 1.3,
      description: 'C-H bond (~1.09 A)',
    },
    {
      element1: 'O',
      element2: 'H',
      expectedMin: 0.8,
      expectedMax: 1.2,
      description: 'O-H bond (~0.96 A)',
    },
    {
      element1: 'C',
      element2: 'O',
      expectedMin: 1.0,
      expectedMax: 1.6,
      description: 'C-O bond (~1.43 A)',
    },
    {
      element1: 'N',
      element2: 'H',
      expectedMin: 0.8,
      expectedMax: 1.2,
      description: 'N-H bond (~1.01 A)',
    },
  ];

  referenceBondLengths.forEach(({ element1, element2, expectedMin, expectedMax, description }) => {
    it(`bond threshold for ${description} should be in reasonable range`, () => {
      const threshold =
        (COVALENT_RADII[element1] || COVALENT_RADII.default) +
        (COVALENT_RADII[element2] || COVALENT_RADII.default);

      // The threshold (sum of covalent radii) should be in a reasonable range
      expect(threshold).toBeGreaterThan(expectedMin * 0.8);
      expect(threshold).toBeLessThan(expectedMax * 1.5);
    });
  });
});

// ---------------------------------------------------------------------------
// 6. Structure Validation
// ---------------------------------------------------------------------------

describe('Physical accuracy: structure validation', () => {
  const converter = new StructureConverter();

  it('should reject structures with no atoms', () => {
    const result = converter.validateStructure({
      atoms: [],
      name: 'empty',
    });
    expect(result.valid).toBe(false);
    expect(result.errors).toContainEqual(expect.stringContaining('No atoms'));
  });

  it('should reject atoms with missing element', () => {
    const result = converter.validateStructure({
      atoms: [{ element: '', x: 0, y: 0, z: 0 }],
      name: 'no-element',
    });
    expect(result.valid).toBe(false);
  });

  it('should reject atoms with NaN coordinates', () => {
    const result = converter.validateStructure({
      atoms: [{ element: 'C', x: NaN, y: 0, z: 0 }],
      name: 'nan-coords',
    });
    expect(result.valid).toBe(false);
    expect(result.errors).toContainEqual(expect.stringContaining('invalid coordinates'));
  });

  it('should reject atoms with infinite coordinates', () => {
    const result = converter.validateStructure({
      atoms: [{ element: 'C', x: Infinity, y: 0, z: 0 }],
      name: 'inf-coords',
    });
    expect(result.valid).toBe(false);
  });

  it('should warn about overlapping atoms (bond length < 0.01 A)', () => {
    const result = converter.validateStructure({
      atoms: [
        { element: 'C', x: 0, y: 0, z: 0 },
        { element: 'C', x: 0.001, y: 0, z: 0 },
      ],
      name: 'overlapping',
    });
    expect(result.warnings.length).toBeGreaterThan(0);
    expect(result.warnings[0]).toContain('very close');
  });

  it('should accept a valid water molecule', () => {
    const result = converter.validateStructure({
      atoms: [
        { element: 'O', x: 0, y: 0, z: 0 },
        { element: 'H', x: 0.9572, y: 0, z: 0 },
        { element: 'H', x: -0.2399, y: 0.9267, z: 0 },
      ],
      name: 'water',
    });
    expect(result.valid).toBe(true);
    expect(result.errors).toHaveLength(0);
  });
});

// ---------------------------------------------------------------------------
// 7. XYZ Format Parsing Accuracy
// ---------------------------------------------------------------------------

describe('Physical accuracy: XYZ format parsing', () => {
  const converter = new StructureConverter();

  it('should parse a simple water molecule in XYZ format with correct coordinates', () => {
    const xyz = `3
water molecule
O  0.0000  0.0000  0.0000
H  0.9572  0.0000  0.0000
H -0.2399  0.9267  0.0000`;

    const structure = converter.fromXYZ(xyz, 'water.xyz');
    expect(structure.atoms).toHaveLength(3);

    // Check oxygen position
    expect(structure.atoms[0].element).toBe('O');
    expect(structure.atoms[0].x).toBeCloseTo(0.0, 4);
    expect(structure.atoms[0].y).toBeCloseTo(0.0, 4);
    expect(structure.atoms[0].z).toBeCloseTo(0.0, 4);

    // Check first hydrogen
    expect(structure.atoms[1].element).toBe('H');
    expect(structure.atoms[1].x).toBeCloseTo(0.9572, 4);

    // Check second hydrogen
    expect(structure.atoms[2].element).toBe('H');
    expect(structure.atoms[2].x).toBeCloseTo(-0.2399, 4);
    expect(structure.atoms[2].y).toBeCloseTo(0.9267, 4);
  });

  it('should reject XYZ files that are too short', () => {
    expect(() => converter.fromXYZ('1\nH', 'bad.xyz')).toThrow('file too short');
  });

  it('should reject XYZ files with invalid atom count', () => {
    expect(() => converter.fromXYZ('abc\ncomment\nH 0 0 0', 'bad.xyz')).toThrow(
      'first line must be number of atoms'
    );
  });
});

// ---------------------------------------------------------------------------
// 8. Physical Reasonableness Checks
// ---------------------------------------------------------------------------

describe('Physical accuracy: reasonableness checks', () => {
  it('all covalent radii should be between 0.2 and 2.5 A', () => {
    Object.entries(COVALENT_RADII).forEach(([key, radius]) => {
      expect(radius).toBeGreaterThanOrEqual(0.2);
      expect(radius).toBeLessThanOrEqual(2.5);
    });
  });

  it('all van der Waals radii should be between 1.0 and 3.0 A', () => {
    Object.entries(VDW_RADII).forEach(([key, radius]) => {
      expect(radius).toBeGreaterThanOrEqual(1.0);
      expect(radius).toBeLessThanOrEqual(3.0);
    });
  });

  it('the covalent radius sum for H-H should approximate experimental bond length (0.74 A)', () => {
    const hhSum = COVALENT_RADII['H'] * 2;
    // H2 bond length is 0.74 A; sum of covalent radii should be close
    expect(Math.abs(hhSum - 0.74)).toBeLessThan(0.15);
  });

  it('the covalent radius sum for C-C should approximate experimental bond length (1.54 A)', () => {
    const ccSum = COVALENT_RADII['C'] * 2;
    // Diamond/cc bond length is 1.54 A
    expect(Math.abs(ccSum - 1.54)).toBeLessThan(0.15);
  });
});
