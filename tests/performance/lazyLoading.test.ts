/**
 * Unit tests for Lazy Loading module
 */

import {
  StructureLazyLoader,
  createLazyLoadedStructure,
  DEFAULT_CHUNK_SIZE,
  LAZY_LOAD_THRESHOLD,
  DEFAULT_LAZY_CONFIG,
} from '../../src/performance/lazyLoading';
import { ASEAtoms } from '../../src/ase/ASEConverter';

describe('Lazy Loading Module', () => {
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

  describe('StructureLazyLoader', () => {
    let loader: StructureLazyLoader;

    beforeEach(() => {
      loader = new StructureLazyLoader();
    });

    describe('shouldUseLazyLoading', () => {
      it('should return false for small structures', () => {
        expect(loader.shouldUseLazyLoading(1000)).toBe(false);
        expect(loader.shouldUseLazyLoading(LAZY_LOAD_THRESHOLD - 1)).toBe(false);
      });

      it('should return true for large structures', () => {
        expect(loader.shouldUseLazyLoading(LAZY_LOAD_THRESHOLD + 1)).toBe(true);
        expect(loader.shouldUseLazyLoading(10000)).toBe(true);
      });
    });

    describe('chunk calculations', () => {
      it('should calculate correct chunk index', () => {
        expect(loader.getChunkIndex(0)).toBe(0);
        expect(loader.getChunkIndex(DEFAULT_CHUNK_SIZE - 1)).toBe(0);
        expect(loader.getChunkIndex(DEFAULT_CHUNK_SIZE)).toBe(1);
      });

      it('should calculate correct chunk count', () => {
        expect(loader.calculateChunkCount(500)).toBe(1);
        expect(loader.calculateChunkCount(DEFAULT_CHUNK_SIZE + 1)).toBe(2);
      });
    });
  });

  describe('createLazyLoadedStructure', () => {
    it('should create lazy loaded structure', () => {
      const atoms = createMockAtoms(2500);
      const lazy = createLazyLoadedStructure(atoms);

      expect(lazy.totalAtoms).toBe(2500);
      expect(lazy.chunkCount).toBe(3);
      expect(lazy.chunkSize).toBe(DEFAULT_CHUNK_SIZE);
    });

    it('should retrieve atoms in range', async () => {
      const atoms = createMockAtoms(2500);
      const lazy = createLazyLoadedStructure(atoms);

      const result = await lazy.getAtomsInRange(500, 1500);
      expect(result).not.toBeNull();
      expect(result!.chemical_symbols.length).toBe(1000);
    });

    it('should retrieve single chunk', async () => {
      const atoms = createMockAtoms(2500);
      const lazy = createLazyLoadedStructure(atoms);

      const chunk = await lazy.getChunk(0);
      expect(chunk).not.toBeNull();
      expect(chunk!.chemical_symbols.length).toBe(DEFAULT_CHUNK_SIZE);
    });
  });
});
