/**
 * Unit tests for Compute Worker
 */

import {
  ComputeWorker,
  WorkerMessageType,
  WorkerMessage,
} from '../../src/performance/computeWorker';
import { ASEAtoms } from '../../src/ase/ASEConverter';

describe('ComputeWorker', () => {
  let worker: ComputeWorker;

  const createMockAtoms = (count: number): ASEAtoms => ({
    chemical_symbols: Array(count).fill('C'),
    positions: Array(count)
      .fill(0)
      .map((_, i) => [i * 1.5, i * 1.5, i * 1.5]),
    pbc: [true, true, true],
    cell: [
      [10, 0, 0],
      [0, 10, 0],
      [0, 0, 10],
    ],
    info: {},
  });

  beforeEach(() => {
    worker = new ComputeWorker();
  });

  describe('PARSE_STRUCTURE', () => {
    it('should parse XYZ structure content into real atom symbols and coordinates', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.PARSE_STRUCTURE,
        id: 'test-1',
        payload: {
          content: '3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n',
          format: 'xyz',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result.chemical_symbols).toEqual(['O', 'H', 'H']);
      expect(response.result.positions).toEqual([
        [0, 0, 0],
        [0.96, 0, 0],
        [-0.24, 0.93, 0],
      ]);
      expect(response.result.pbc).toEqual([false, false, false]);
      expect(response.result.info).toMatchObject({ format: 'xyz', comment: 'water' });
      expect(response.duration).toBeGreaterThanOrEqual(0);
    });

    it('should parse VASP POSCAR Direct coordinates into Cartesian positions', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.PARSE_STRUCTURE,
        id: 'test-poscar',
        payload: {
          content: [
            'Si',
            '1.0',
            '5.43 0 0',
            '0 5.43 0',
            '0 0 5.43',
            'Si',
            '2',
            'Direct',
            '0 0 0',
            '0.25 0.25 0.25',
          ].join('\n'),
          format: 'poscar',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result.chemical_symbols).toEqual(['Si', 'Si']);
      expect(response.result.positions[1]).toEqual([1.3575, 1.3575, 1.3575]);
      expect(response.result.cell).toEqual([
        [5.43, 0, 0],
        [0, 5.43, 0],
        [0, 0, 5.43],
      ]);
      expect(response.result.pbc).toEqual([true, true, true]);
    });

    it('should reject unsupported structure formats instead of returning a fake structure', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.PARSE_STRUCTURE,
        id: 'test-unsupported-format',
        payload: {
          content: 'not a structure',
          format: 'unknown-format',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(false);
      expect(response.error).toContain('Unsupported structure format');
    });

    it('should reject XYZ content with mismatched declared atom count', async () => {
      const response = await worker.processMessage({
        type: WorkerMessageType.PARSE_STRUCTURE,
        id: 'test-xyz-mismatch',
        payload: {
          content: '2\nbad\nH 0 0 0\n',
          format: 'xyz',
        },
      });

      expect(response.success).toBe(false);
      expect(response.error).toContain('declared 2 atoms but parsed 1');
    });

    it('should reject XYZ content with non-finite or missing coordinates', async () => {
      const response = await worker.processMessage({
        type: WorkerMessageType.PARSE_STRUCTURE,
        id: 'test-xyz-bad-coords',
        payload: {
          content: '1\nbad\nH NaN 0 0\n',
          format: 'xyz',
        },
      });

      expect(response.success).toBe(false);
      expect(response.error).toContain('declared 1 atoms but parsed 0');
    });

    it('should fail with empty content', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.PARSE_STRUCTURE,
        id: 'test-2',
        payload: {
          content: '',
          format: 'xyz',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(false);
      expect(response.error).toBeDefined();
    });
  });

  describe('CONVERT_FORMAT', () => {
    it('should convert atoms to standard XYZ text', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['O', 'H', 'H'],
        positions: [
          [0, 0, 0],
          [0.9572, 0, 0],
          [-0.239, 0.927, 0],
        ],
        pbc: [false, false, false],
        info: { comment: 'water' },
      };
      const message: WorkerMessage = {
        type: WorkerMessageType.CONVERT_FORMAT,
        id: 'test-3',
        payload: {
          atoms,
          targetFormat: 'xyz',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result.split('\n')).toEqual([
        '3',
        'water',
        'O 0.000000 0.000000 0.000000',
        'H 0.957200 0.000000 0.000000',
        'H -0.239000 0.927000 0.000000',
      ]);
    });

    it('should convert atoms to grouped POSCAR text', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['Si', 'O', 'O', 'Si'],
        positions: [
          [0, 0, 0],
          [1.6, 0, 0],
          [0, 1.6, 0],
          [2.7, 2.7, 2.7],
        ],
        pbc: [true, true, true],
        cell: [
          [5.4, 0, 0],
          [0, 5.4, 0],
          [0, 0, 5.4],
        ],
        info: { title: 'SiO2' },
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.CONVERT_FORMAT,
        id: 'test-poscar-export',
        payload: { atoms, targetFormat: 'poscar' },
      });

      expect(response.success).toBe(true);
      expect(response.result.split('\n')).toEqual([
        'SiO2',
        '1.000000000000',
        '5.400000 0.000000 0.000000',
        '0.000000 5.400000 0.000000',
        '0.000000 0.000000 5.400000',
        'Si O',
        '2 2',
        'Cartesian',
        '0.000000 0.000000 0.000000',
        '2.700000 2.700000 2.700000',
        '1.600000 0.000000 0.000000',
        '0.000000 1.600000 0.000000',
      ]);
    });

    it('should reject POSCAR conversion without unit-cell vectors', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['H'],
        positions: [[0, 0, 0]],
        pbc: [false, false, false],
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.CONVERT_FORMAT,
        id: 'test-poscar-no-cell',
        payload: { atoms, targetFormat: 'vasp' },
      });

      expect(response.success).toBe(false);
      expect(response.error).toContain('without unit-cell vectors');
    });

    it('should reject unsupported worker conversion targets', async () => {
      const response = await worker.processMessage({
        type: WorkerMessageType.CONVERT_FORMAT,
        id: 'test-unsupported-convert',
        payload: {
          atoms: createMockAtoms(1),
          targetFormat: 'cube',
        },
      });

      expect(response.success).toBe(false);
      expect(response.error).toContain('Unsupported worker conversion target format');
    });

    it('should fail with invalid atoms', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.CONVERT_FORMAT,
        id: 'test-4',
        payload: {
          atoms: null,
          targetFormat: 'vasp',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(false);
      expect(response.error).toBeDefined();
    });
  });

  describe('VALIDATE_STRUCTURE', () => {
    it('should validate structure successfully', async () => {
      const atoms = createMockAtoms(10);
      const message: WorkerMessage = {
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-5',
        payload: {
          atoms,
          checks: ['bond_lengths', 'cell_consistency', 'atom_overlap'],
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result).toBeDefined();
      expect(response.result.bond_lengths).toBeDefined();
      expect(response.result.cell_consistency).toBeDefined();
      expect(response.result.atom_overlap).toBeDefined();
    });

    it('should detect atom overlap', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['C', 'C'],
        positions: [
          [0, 0, 0],
          [0.3, 0, 0],
        ], // Very close atoms
        pbc: [false, false, false],
      };

      const message: WorkerMessage = {
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-6',
        payload: {
          atoms,
          checks: ['atom_overlap'],
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result.atom_overlap.valid).toBe(false);
      expect(response.result.atom_overlap.overlaps.length).toBeGreaterThan(0);
    });

    it('should reject malformed atoms before running structure checks', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['C', 'O'],
        positions: [[0, 0, 0]],
        pbc: [false, false, false],
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-malformed-atoms',
        payload: {
          atoms,
          checks: ['bond_lengths'],
        },
      });

      expect(response.success).toBe(false);
      expect(response.error).toContain('chemical_symbols length must match positions length');
    });

    it('should reject periodic atoms with missing cell for cell consistency checks', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['C'],
        positions: [[0, 0, 0]],
        pbc: [true, true, true],
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-periodic-no-cell',
        payload: {
          atoms,
          checks: ['cell_consistency'],
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.cell_consistency.valid).toBe(false);
      expect(response.result.cell_consistency.warnings[0]).toContain('Periodic system is missing');
    });

    it('should validate charge neutrality from explicit charge arrays', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['Na', 'Cl'],
        positions: [
          [0, 0, 0],
          [2.8, 0, 0],
        ],
        pbc: [false, false, false],
        info: { charges: [1, -1] },
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-charge-neutral',
        payload: {
          atoms,
          checks: ['charge_neutrality'],
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.charge_neutrality).toMatchObject({
        valid: true,
        netCharge: 0,
        status: 'checked',
      });
      expect(response.result.charge_neutrality.warnings).toEqual([]);
    });

    it('should mark charge neutrality as skipped when no charge data is present', async () => {
      const atoms = createMockAtoms(2);

      const response = await worker.processMessage({
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-charge-skipped',
        payload: {
          atoms,
          checks: ['charge_neutrality'],
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.charge_neutrality).toMatchObject({
        valid: true,
        status: 'skipped',
        netCharge: null,
      });
    });

    it('should flag non-neutral explicit charge arrays', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['Na', 'Cl'],
        positions: [
          [0, 0, 0],
          [2.8, 0, 0],
        ],
        pbc: [false, false, false],
        info: { charges: [1, 0] },
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-charge-non-neutral',
        payload: {
          atoms,
          checks: ['charge_neutrality'],
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.charge_neutrality.valid).toBe(false);
      expect(response.result.charge_neutrality.status).toBe('checked');
      expect(response.result.charge_neutrality.netCharge).toBe(1);
      expect(response.result.charge_neutrality.warnings[0]).toContain('Net charge is 1.000000');
    });
  });

  describe('CALCULATE_PROPERTIES', () => {
    it('should calculate molecular properties', async () => {
      const atoms = createMockAtoms(10);
      const message: WorkerMessage = {
        type: WorkerMessageType.CALCULATE_PROPERTIES,
        id: 'test-7',
        payload: {
          atoms,
          properties: ['center_of_mass', 'bounding_box', 'atom_count'],
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result).toBeDefined();
      expect(response.result.center_of_mass).toBeDefined();
      expect(response.result.bounding_box).toBeDefined();
      expect(response.result.atom_count).toBe(10);
    });

    it('should calculate mass-weighted center of mass correctly', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['X', 'X'],
        positions: [
          [0, 0, 0],
          [2, 0, 0],
        ],
        pbc: [false, false, false],
        masses: [1, 3],
      };

      const message: WorkerMessage = {
        type: WorkerMessageType.CALCULATE_PROPERTIES,
        id: 'test-8',
        payload: {
          atoms,
          properties: ['center_of_mass'],
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result.center_of_mass).toEqual([1.5, 0, 0]);
    });

    it('should calculate a real moment of inertia tensor about the center of mass', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['X', 'X'],
        positions: [
          [-1, 0, 0],
          [1, 0, 0],
        ],
        pbc: [false, false, false],
        masses: [1, 1],
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.CALCULATE_PROPERTIES,
        id: 'test-inertia',
        payload: {
          atoms,
          properties: ['moment_of_inertia'],
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.moment_of_inertia).toEqual([
        [0, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
      ]);
    });

    it('should calculate exact finite bounding boxes for negative coordinates', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['C', 'O', 'H'],
        positions: [
          [-1, 2, 0],
          [3, -4, 5],
          [0, 1, -2],
        ],
        pbc: [false, false, false],
      };

      const response = await worker.processMessage({
        type: WorkerMessageType.CALCULATE_PROPERTIES,
        id: 'test-bounds',
        payload: {
          atoms,
          properties: ['bounding_box'],
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.bounding_box).toEqual({
        min: [-1, -4, -2],
        max: [3, 2, 5],
        size: [4, 6, 7],
      });
    });

    it('should reject empty atom payloads before calculating properties', async () => {
      const response = await worker.processMessage({
        type: WorkerMessageType.CALCULATE_PROPERTIES,
        id: 'test-empty-properties',
        payload: {
          atoms: { chemical_symbols: [], positions: [], pbc: [false, false, false] },
          properties: ['bounding_box'],
        },
      });

      expect(response.success).toBe(false);
      expect(response.error).toContain('atoms must contain at least one atom');
    });
  });

  describe('MIGRATE_PARAMETERS', () => {
    it('should migrate VASP to CP2K parameters', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.MIGRATE_PARAMETERS,
        id: 'test-9',
        payload: {
          sourceFormat: 'vasp',
          targetFormat: 'cp2k',
          parameters: {
            ENCUT: 520,
            EDIFF: 1e-6,
          },
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result.migrated).toBeDefined();
      expect(response.result.migrated.CUTOFF).toBeCloseTo(38.219295063, 9);
      expect(response.result.migrated.EPS_SCF).toBe(1e-6);
      expect(response.result.notes.join('\n')).toContain('ENCUT -> CUTOFF converted with factor');
    });

    it('should migrate QE to VASP parameters', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.MIGRATE_PARAMETERS,
        id: 'test-10',
        payload: {
          sourceFormat: 'qe',
          targetFormat: 'vasp',
          parameters: {
            ecutwfc: 60,
          },
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result.migrated.ENCUT).toBeCloseTo(816.34158738, 8);
    });

    it('should migrate VASP to QE parameters with unit conversion and convergence mapping', async () => {
      const response = await worker.processMessage({
        type: WorkerMessageType.MIGRATE_PARAMETERS,
        id: 'test-vasp-qe',
        payload: {
          sourceFormat: 'vasp',
          targetFormat: 'qe',
          parameters: {
            ENCUT: 408.17079368982,
            EDIFF: 1e-7,
            UNKNOWN: true,
          },
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.migrated.ecutwfc).toBeCloseTo(30, 10);
      expect(response.result.migrated.conv_thr).toBe(1e-7);
      expect(response.result.unmapped).toEqual(['UNKNOWN']);
      expect(response.result.warnings[0]).toContain('could not be mapped');
    });

    it('should report unsupported parameter migration pairs without pretending success', async () => {
      const response = await worker.processMessage({
        type: WorkerMessageType.MIGRATE_PARAMETERS,
        id: 'test-unsupported-migration',
        payload: {
          sourceFormat: 'orca',
          targetFormat: 'lammps',
          parameters: {
            MaxIter: 200,
          },
        },
      });

      expect(response.success).toBe(true);
      expect(response.result.supported).toBe(false);
      expect(response.result.migrated).toEqual({});
      expect(response.result.unmapped).toEqual(['MaxIter']);
      expect(response.result.warnings[0]).toContain('No worker-local parameter mappings');
    });

    it('should reject non-object migration parameter payloads', async () => {
      const response = await worker.processMessage({
        type: WorkerMessageType.MIGRATE_PARAMETERS,
        id: 'test-bad-params',
        payload: {
          sourceFormat: 'vasp',
          targetFormat: 'qe',
          parameters: null,
        },
      });

      expect(response.success).toBe(false);
      expect(response.error).toContain('parameters must be a key-value object');
    });
  });

  describe('Error Handling', () => {
    it('should handle unknown message type', async () => {
      const message: WorkerMessage = {
        type: 'UNKNOWN_TYPE' as WorkerMessageType,
        id: 'test-11',
        payload: {},
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(false);
      expect(response.error).toContain('Unknown message type');
    });
  });

  describe('Performance', () => {
    it('should track execution duration', async () => {
      const atoms = createMockAtoms(100);
      const message: WorkerMessage = {
        type: WorkerMessageType.VALIDATE_STRUCTURE,
        id: 'test-12',
        payload: {
          atoms,
          checks: ['bond_lengths'],
        },
      };

      const response = await worker.processMessage(message);

      expect(response.duration).toBeDefined();
      expect(response.duration).toBeGreaterThanOrEqual(0);
    });
  });
});
