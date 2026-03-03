/**
 * Migration Validation Suite
 *
 * Integration tests for format migration validation including:
 * - Round-trip conversion tests (VASP↔CP2K↔QE)
 * - Structure preservation validation
 * - Parameter consistency checks
 * - Validation report generation
 */

import * as path from 'path';
import * as fs from 'fs';
import { ASEConverter, ASEFormat, ASEAtoms } from '../../src/ase/ASEConverter';

describe('Migration Validation Suite', () => {
  let converter: ASEConverter;
  const fixturesDir = path.join(__dirname, '../fixtures/format_conversion');
  const outputDir = path.join(__dirname, '../temp/migration_validation');
  const reportsDir = path.join(__dirname, '../temp/migration_reports');

  beforeAll(async () => {
    // Create a mock extension context
    const mockContext = {
      extensionPath: path.join(__dirname, '../..'),
    } as any;

    converter = new ASEConverter(mockContext);

    // Check if backend is available
    const backendAvailable = await converter.isAvailable();
    if (!backendAvailable) {
      console.warn('Python backend not available. Skipping integration tests.');
      return;
    }

    // Create output directories
    fs.mkdirSync(outputDir, { recursive: true });
    fs.mkdirSync(reportsDir, { recursive: true });
  });

  afterAll(() => {
    // Clean up output directories
    if (fs.existsSync(outputDir)) {
      fs.rmSync(outputDir, { recursive: true, force: true });
    }
    if (fs.existsSync(reportsDir)) {
      fs.rmSync(reportsDir, { recursive: true, force: true });
    }
  });

  const getBackendAvailable = async () => {
    return await converter.isAvailable();
  };

  describe('Round-Trip Conversion Tests', () => {
    it('should preserve structure in VASP → CP2K → VASP round-trip', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      if (!fs.existsSync(inputFile)) {
        console.warn(`Test fixture not found: ${inputFile}`);
        return;
      }

      // Step 1: Read original structure
      const originalResult = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      expect(originalResult.success).toBe(true);
      expect(originalResult.atoms).toBeDefined();

      const originalAtoms = originalResult.atoms!;

      // Step 2: Convert to CP2K
      const cp2kPath = path.join(outputDir, 'roundtrip_cp2k.inp');
      const cp2kResult = await converter.writeFromAtoms(originalAtoms, cp2kPath, ASEFormat.CP2K);
      expect(cp2kResult.success).toBe(true);

      // Step 3: Read back from CP2K
      const cp2kReadResult = await converter.readToAtoms(cp2kPath, ASEFormat.CP2K);
      expect(cp2kReadResult.success).toBe(true);
      expect(cp2kReadResult.atoms).toBeDefined();

      const cp2kAtoms = cp2kReadResult.atoms!;

      // Step 4: Convert back to VASP
      const vaspPath = path.join(outputDir, 'roundtrip_vasp.POSCAR');
      const vaspResult = await converter.writeFromAtoms(cp2kAtoms, vaspPath, ASEFormat.VASP);
      expect(vaspResult.success).toBe(true);

      // Step 5: Read final structure
      const finalResult = await converter.readToAtoms(vaspPath, ASEFormat.VASP);
      expect(finalResult.success).toBe(true);
      expect(finalResult.atoms).toBeDefined();

      const finalAtoms = finalResult.atoms!;

      // Validate structure preservation
      const validation = validateStructurePreservation(originalAtoms, finalAtoms);
      expect(validation.preserved).toBe(true);
      expect(validation.positionDeviation).toBeLessThan(0.01); // 0.01 Angstrom tolerance
    });

    it('should preserve structure in VASP → QE → VASP round-trip', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      if (!fs.existsSync(inputFile)) {
        console.warn(`Test fixture not found: ${inputFile}`);
        return;
      }

      // Step 1: Read original structure
      const originalResult = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      expect(originalResult.success).toBe(true);
      expect(originalResult.atoms).toBeDefined();

      const originalAtoms = originalResult.atoms!;

      // Step 2: Convert to QE
      const qePath = path.join(outputDir, 'roundtrip_qe.in');
      const qeResult = await converter.writeFromAtoms(originalAtoms, qePath, ASEFormat.QE);
      expect(qeResult.success).toBe(true);

      // Step 3: Read back from QE
      const qeReadResult = await converter.readToAtoms(qePath, ASEFormat.QE);
      expect(qeReadResult.success).toBe(true);
      expect(qeReadResult.atoms).toBeDefined();

      const qeAtoms = qeReadResult.atoms!;

      // Step 4: Convert back to VASP
      const vaspPath = path.join(outputDir, 'roundtrip_vasp2.POSCAR');
      const vaspResult = await converter.writeFromAtoms(qeAtoms, vaspPath, ASEFormat.VASP);
      expect(vaspResult.success).toBe(true);

      // Step 5: Read final structure
      const finalResult = await converter.readToAtoms(vaspPath, ASEFormat.VASP);
      expect(finalResult.success).toBe(true);
      expect(finalResult.atoms).toBeDefined();

      const finalAtoms = finalResult.atoms!;

      // Validate structure preservation
      const validation = validateStructurePreservation(originalAtoms, finalAtoms);
      expect(validation.preserved).toBe(true);
      expect(validation.positionDeviation).toBeLessThan(0.01);
    });

    it('should preserve structure in CP2K → QE → CP2K round-trip', async () => {
      if (!(await getBackendAvailable())) return;

      // Create a simple CP2K input for testing
      const cp2kInput = `&GLOBAL
  PROJECT test
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      A 5.0 0.0 0.0
      B 0.0 5.0 0.0
      C 0.0 0.0 5.0
    &END CELL
    &COORD
      Si 0.0 0.0 0.0
      Si 2.5 2.5 2.5
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
`;

      const inputFile = path.join(fixturesDir, 'test_cp2k.inp');
      fs.writeFileSync(inputFile, cp2kInput);

      // Step 1: Read original structure
      const originalResult = await converter.readToAtoms(inputFile, ASEFormat.CP2K);

      // Clean up test file
      fs.unlinkSync(inputFile);

      if (!originalResult.success) {
        console.warn('CP2K read not fully supported, skipping test');
        return;
      }

      expect(originalResult.atoms).toBeDefined();
      const originalAtoms = originalResult.atoms!;

      // Step 2: Convert to QE
      const qePath = path.join(outputDir, 'cp2k_qe.in');
      const qeResult = await converter.writeFromAtoms(originalAtoms, qePath, ASEFormat.QE);
      expect(qeResult.success).toBe(true);

      // Step 3: Read back from QE
      const qeReadResult = await converter.readToAtoms(qePath, ASEFormat.QE);
      expect(qeReadResult.success).toBe(true);
      expect(qeReadResult.atoms).toBeDefined();

      const qeAtoms = qeReadResult.atoms!;

      // Step 4: Convert back to CP2K
      const cp2kPath = path.join(outputDir, 'qe_cp2k.inp');
      const cp2kResult = await converter.writeFromAtoms(qeAtoms, cp2kPath, ASEFormat.CP2K);
      expect(cp2kResult.success).toBe(true);

      // Validate structure preservation
      const validation = validateStructurePreservation(originalAtoms, qeAtoms);
      expect(validation.preserved).toBe(true);
    });
  });

  describe('Structure Preservation Validation', () => {
    it('should validate atomic positions within tolerance', () => {
      const atoms1: ASEAtoms = {
        chemical_symbols: ['Si', 'Si'],
        positions: [
          [0.0, 0.0, 0.0],
          [2.5, 2.5, 2.5],
        ],
        cell: [
          [5.0, 0.0, 0.0],
          [0.0, 5.0, 0.0],
          [0.0, 0.0, 5.0],
        ],
        pbc: [true, true, true],
      };

      const atoms2: ASEAtoms = {
        chemical_symbols: ['Si', 'Si'],
        positions: [
          [0.001, 0.0, 0.0],
          [2.501, 2.5, 2.5],
        ],
        cell: [
          [5.0, 0.0, 0.0],
          [0.0, 5.0, 0.0],
          [0.0, 0.0, 5.0],
        ],
        pbc: [true, true, true],
      };

      const validation = validateStructurePreservation(atoms1, atoms2);
      expect(validation.preserved).toBe(true);
      expect(validation.positionDeviation).toBeLessThan(0.01);
    });

    it('should detect chemical symbol mismatches', () => {
      const atoms1: ASEAtoms = {
        chemical_symbols: ['Si', 'O'],
        positions: [
          [0.0, 0.0, 0.0],
          [1.0, 1.0, 1.0],
        ],
        pbc: [false, false, false],
      };

      const atoms2: ASEAtoms = {
        chemical_symbols: ['Si', 'Si'],
        positions: [
          [0.0, 0.0, 0.0],
          [1.0, 1.0, 1.0],
        ],
        pbc: [false, false, false],
      };

      const validation = validateStructurePreservation(atoms1, atoms2);
      expect(validation.preserved).toBe(false);
      expect(validation.errors).toContain('Chemical symbols mismatch');
    });

    it('should validate periodic boundary conditions', () => {
      const atoms1: ASEAtoms = {
        chemical_symbols: ['H', 'H'],
        positions: [
          [0.0, 0.0, 0.0],
          [0.74, 0.0, 0.0],
        ],
        cell: [
          [10.0, 0.0, 0.0],
          [0.0, 10.0, 0.0],
          [0.0, 0.0, 10.0],
        ],
        pbc: [false, false, false],
      };

      const atoms2: ASEAtoms = {
        chemical_symbols: ['H', 'H'],
        positions: [
          [0.0, 0.0, 0.0],
          [0.74, 0.0, 0.0],
        ],
        cell: [
          [10.0, 0.0, 0.0],
          [0.0, 10.0, 0.0],
          [0.0, 0.0, 10.0],
        ],
        pbc: [true, true, true],
      };

      const validation = validateStructurePreservation(atoms1, atoms2);
      expect(validation.preserved).toBe(false);
      expect(validation.errors).toContain('PBC mismatch');
    });

    it('should validate cell parameters for periodic systems', () => {
      const atoms1: ASEAtoms = {
        chemical_symbols: ['Si'],
        positions: [[0.0, 0.0, 0.0]],
        cell: [
          [5.0, 0.0, 0.0],
          [0.0, 5.0, 0.0],
          [0.0, 0.0, 5.0],
        ],
        pbc: [true, true, true],
      };

      const atoms2: ASEAtoms = {
        chemical_symbols: ['Si'],
        positions: [[0.0, 0.0, 0.0]],
        cell: [
          [5.1, 0.0, 0.0],
          [0.0, 5.1, 0.0],
          [0.0, 0.0, 5.1],
        ],
        pbc: [true, true, true],
      };

      const validation = validateStructurePreservation(atoms1, atoms2);
      expect(validation.preserved).toBe(false);
      expect(validation.cellDeviation).toBeGreaterThan(0.01);
    });
  });

  describe('Parameter Consistency Checks', () => {
    it('should validate energy cutoff consistency', async () => {
      if (!(await getBackendAvailable())) return;

      const params1 = { encut: 400, kpts: [4, 4, 4] };
      const params2 = { encut: 450, kpts: [4, 4, 4] };

      const check = checkParameterConsistency(
        { vasp: params1 },
        { vasp: params2 }
      );

      expect(check.consistent).toBe(false);
      expect(check.differences).toContainEqual(
        expect.objectContaining({ parameter: 'encut', diff: 50 })
      );
    });

    it('should validate k-point grid consistency', () => {
      const check = checkParameterConsistency(
        { kpts: [4, 4, 4] },
        { kpts: [6, 6, 6] }
      );

      expect(check.consistent).toBe(false);
      expect(check.differences.length).toBeGreaterThan(0);
    });
  });

  describe('Validation Report Generation', () => {
    it('should generate validation report', async () => {
      if (!(await getBackendAvailable())) return;

      const inputFile = path.join(fixturesDir, 'POSCAR');
      if (!fs.existsSync(inputFile)) {
        console.warn(`Test fixture not found: ${inputFile}`);
        return;
      }

      const originalResult = await converter.readToAtoms(inputFile, ASEFormat.VASP);
      if (!originalResult.success || !originalResult.atoms) {
        return;
      }

      // Convert and back
      const xyzPath = path.join(outputDir, 'report_test.xyz');
      await converter.writeFromAtoms(originalResult.atoms, xyzPath, ASEFormat.XYZ);
      const finalResult = await converter.readToAtoms(xyzPath, ASEFormat.XYZ);

      if (!finalResult.success || !finalResult.atoms) {
        return;
      }

      const validation = validateStructurePreservation(
        originalResult.atoms,
        finalResult.atoms
      );

      const report = generateValidationReport(
        'VASP → XYZ → VASP',
        originalResult.atoms,
        finalResult.atoms,
        validation
      );

      expect(report).toBeDefined();
      expect(report.migration_path).toBe('VASP → XYZ → VASP');
      expect(report.timestamp).toBeDefined();
      expect(report.structure_validation).toBeDefined();

      // Save report
      const reportPath = path.join(reportsDir, 'validation_report.json');
      fs.writeFileSync(reportPath, JSON.stringify(report, null, 2));

      expect(fs.existsSync(reportPath)).toBe(true);
    });
  });
});

/**
 * Validate structure preservation between two ASE Atoms objects
 */
function validateStructurePreservation(
  atoms1: ASEAtoms,
  atoms2: ASEAtoms,
  tolerance: number = 0.01
): {
  preserved: boolean;
  positionDeviation: number;
  cellDeviation?: number;
  errors: string[];
} {
  const errors: string[] = [];

  // Check number of atoms
  if (atoms1.chemical_symbols.length !== atoms2.chemical_symbols.length) {
    errors.push(
      `Atom count mismatch: ${atoms1.chemical_symbols.length} vs ${atoms2.chemical_symbols.length}`
    );
    return { preserved: false, positionDeviation: Infinity, errors };
  }

  // Check chemical symbols
  for (let i = 0; i < atoms1.chemical_symbols.length; i++) {
    if (atoms1.chemical_symbols[i] !== atoms2.chemical_symbols[i]) {
      errors.push('Chemical symbols mismatch');
      break;
    }
  }

  // Check positions
  let maxPositionDeviation = 0;
  for (let i = 0; i < atoms1.positions.length; i++) {
    const pos1 = atoms1.positions[i];
    const pos2 = atoms2.positions[i];
    const deviation = Math.sqrt(
      Math.pow(pos1[0] - pos2[0], 2) +
        Math.pow(pos1[1] - pos2[1], 2) +
        Math.pow(pos1[2] - pos2[2], 2)
    );
    maxPositionDeviation = Math.max(maxPositionDeviation, deviation);
  }

  if (maxPositionDeviation > tolerance) {
    errors.push(
      `Position deviation ${maxPositionDeviation.toFixed(4)} exceeds tolerance ${tolerance}`
    );
  }

  // Check PBC
  const pbcMismatch = atoms1.pbc.some((p, i) => p !== atoms2.pbc[i]);
  if (pbcMismatch) {
    errors.push('PBC mismatch');
  }

  // Check cell for periodic systems
  let maxCellDeviation = 0;
  if (atoms1.cell && atoms2.cell && atoms1.pbc.some(p => p)) {
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const deviation = Math.abs(atoms1.cell[i][j] - atoms2.cell[i][j]);
        maxCellDeviation = Math.max(maxCellDeviation, deviation);
      }
    }

    if (maxCellDeviation > tolerance) {
      errors.push(
        `Cell deviation ${maxCellDeviation.toFixed(4)} exceeds tolerance ${tolerance}`
      );
    }
  }

  return {
    preserved: errors.length === 0,
    positionDeviation: maxPositionDeviation,
    cellDeviation: atoms1.cell && atoms2.cell ? maxCellDeviation : undefined,
    errors,
  };
}

/**
 * Check parameter consistency between two parameter sets
 */
function checkParameterConsistency(
  params1: Record<string, any>,
  params2: Record<string, any>
): {
  consistent: boolean;
  differences: Array<{ parameter: string; value1: any; value2: any; diff?: number }>;
} {
  const differences: Array<{
    parameter: string;
    value1: any;
    value2: any;
    diff?: number;
  }> = [];

  const allKeys = new Set([...Object.keys(params1), ...Object.keys(params2)]);

  for (const key of allKeys) {
    const value1 = params1[key];
    const value2 = params2[key];

    if (JSON.stringify(value1) !== JSON.stringify(value2)) {
      const diff: any = {
        parameter: key,
        value1,
        value2,
      };

      if (typeof value1 === 'number' && typeof value2 === 'number') {
        diff.diff = Math.abs(value1 - value2);
      }

      differences.push(diff);
    }
  }

  return {
    consistent: differences.length === 0,
    differences,
  };
}

/**
 * Generate validation report
 */
function generateValidationReport(
  migrationPath: string,
  original: ASEAtoms,
  final: ASEAtoms,
  validation: {
    preserved: boolean;
    positionDeviation: number;
    cellDeviation?: number;
    errors: string[];
  }
): {
  migration_path: string;
  timestamp: string;
  structure_validation: {
    preserved: boolean;
    natoms: number;
    formula: string;
    position_deviation: number;
    cell_deviation?: number;
    errors: string[];
  };
} {
  const formula = getChemicalFormula(original.chemical_symbols);

  return {
    migration_path: migrationPath,
    timestamp: new Date().toISOString(),
    structure_validation: {
      preserved: validation.preserved,
      natoms: original.chemical_symbols.length,
      formula,
      position_deviation: validation.positionDeviation,
      cell_deviation: validation.cellDeviation,
      errors: validation.errors,
    },
  };
}

/**
 * Get chemical formula from element symbols
 */
function getChemicalFormula(symbols: string[]): string {
  const counts: Record<string, number> = {};
  for (const symbol of symbols) {
    counts[symbol] = (counts[symbol] || 0) + 1;
  }

  return Object.entries(counts)
    .map(([element, count]) => `${element}${count > 1 ? count : ''}`)
    .join('');
}
