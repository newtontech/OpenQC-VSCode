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
    it('should parse structure successfully', async () => {
      const message: WorkerMessage = {
        type: WorkerMessageType.PARSE_STRUCTURE,
        id: 'test-1',
        payload: {
          content: 'C 0 0 0\nC 1 0 0\nC 2 0 0',
          format: 'xyz',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result).toBeDefined();
      expect(response.duration).toBeGreaterThanOrEqual(0);
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
    it('should convert atoms to target format', async () => {
      const atoms = createMockAtoms(10);
      const message: WorkerMessage = {
        type: WorkerMessageType.CONVERT_FORMAT,
        id: 'test-3',
        payload: {
          atoms,
          targetFormat: 'vasp',
        },
      };

      const response = await worker.processMessage(message);

      expect(response.success).toBe(true);
      expect(response.result).toBeDefined();
      expect(response.result).toContain('Converted to vasp');
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

    it('should calculate center of mass correctly', async () => {
      const atoms: ASEAtoms = {
        chemical_symbols: ['C', 'C'],
        positions: [
          [0, 0, 0],
          [2, 0, 0],
        ],
        pbc: [false, false, false],
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
      expect(response.result.center_of_mass).toEqual([1, 0, 0]);
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
      expect(response.result.migrated.CUTOFF).toBe(572); // 520 * 1.1
      expect(response.result.migrated.EPSCF).toBe(1e-6);
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
      expect(response.result.migrated.ENCUT).toBe(60);
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
