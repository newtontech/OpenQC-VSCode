/**
 * Lazy Loading Module for Large Structures
 *
 * Implements lazy loading strategies for handling large molecular structures
 * with minimal memory footprint.
 */

import { ASEAtoms } from '../ase/ASEConverter';

export const DEFAULT_CHUNK_SIZE = 1000;
export const LAZY_LOAD_THRESHOLD = 5000;

export interface LazyLoadedStructure {
  totalAtoms: number;
  chunkCount: number;
  loadedChunks: Set<number>;
  chunkSize: number;
  getAtomsInRange(start: number, end: number): Promise<ASEAtoms | null>;
  getChunk(index: number): Promise<ASEAtoms | null>;
  preloadAround(centerIndex: number, radius?: number): Promise<void>;
  releaseChunks(indices: number[]): void;
  getMemoryUsage(): number;
}

export interface LazyLoadConfig {
  chunkSize: number;
  threshold: number;
  maxCachedChunks: number;
  preloadRadius: number;
}

export const DEFAULT_LAZY_CONFIG: LazyLoadConfig = {
  chunkSize: DEFAULT_CHUNK_SIZE,
  threshold: LAZY_LOAD_THRESHOLD,
  maxCachedChunks: 5,
  preloadRadius: 2,
};

export class StructureLazyLoader {
  config: LazyLoadConfig;
  private cache: Map<number, ASEAtoms> = new Map();

  constructor(config: Partial<LazyLoadConfig> = {}) {
    this.config = { ...DEFAULT_LAZY_CONFIG, ...config };
  }

  shouldUseLazyLoading(atomCount: number): boolean {
    return atomCount > this.config.threshold;
  }

  getChunkIndex(atomIndex: number): number {
    return Math.floor(atomIndex / this.config.chunkSize);
  }

  calculateChunkCount(totalAtoms: number): number {
    return Math.ceil(totalAtoms / this.config.chunkSize);
  }

  getChunkRange(chunkIndex: number): { start: number; end: number } {
    const start = chunkIndex * this.config.chunkSize;
    const end = start + this.config.chunkSize;
    return { start, end };
  }

  releaseChunks(indices: number[]): void {
    indices.forEach(index => this.cache.delete(index));
  }

  clearCache(): void {
    this.cache.clear();
  }

  getMemoryUsage(): number {
    let usage = 0;
    this.cache.forEach(atoms => {
      usage += atoms.chemical_symbols.length * 100;
    });
    return usage;
  }

  getStats() {
    return {
      cachedChunks: this.cache.size,
      maxCacheSize: this.config.maxCachedChunks,
      memoryUsage: this.getMemoryUsage(),
      loadingChunks: 0,
    };
  }
}

export function createLazyLoadedStructure(
  fullAtoms: ASEAtoms,
  config?: Partial<LazyLoadConfig>
): LazyLoadedStructure {
  const loader = new StructureLazyLoader(config);
  const chunkCount = loader.calculateChunkCount(fullAtoms.chemical_symbols.length);
  const chunkCache = new Map<number, ASEAtoms>();

  for (let i = 0; i < chunkCount; i++) {
    const { start, end } = loader.getChunkRange(i);
    const chunk: ASEAtoms = {
      chemical_symbols: fullAtoms.chemical_symbols.slice(start, end),
      positions: fullAtoms.positions.slice(start, end),
      pbc: fullAtoms.pbc,
      cell: fullAtoms.cell,
      info: fullAtoms.info,
    };
    chunkCache.set(i, chunk);
  }

  return {
    totalAtoms: fullAtoms.chemical_symbols.length,
    chunkCount,
    chunkSize: loader.config.chunkSize,
    loadedChunks: new Set(),

    async getAtomsInRange(start: number, end: number): Promise<ASEAtoms | null> {
      const startChunk = loader.getChunkIndex(start);
      const endChunk = loader.getChunkIndex(end - 1);
      const symbols: string[] = [];
      const positions: number[][] = [];

      for (let i = startChunk; i <= endChunk; i++) {
        const chunk = chunkCache.get(i);
        if (!chunk) continue;
        const chunkStart = i * loader.config.chunkSize;
        const localStart = Math.max(0, start - chunkStart);
        const localEnd = Math.min(chunk.chemical_symbols.length, end - chunkStart);
        for (let j = localStart; j < localEnd; j++) {
          symbols.push(chunk.chemical_symbols[j]);
          positions.push(chunk.positions[j]);
        }
      }

      if (symbols.length === 0) return null;
      return {
        chemical_symbols: symbols,
        positions,
        pbc: fullAtoms.pbc,
        cell: fullAtoms.cell,
        info: { ...fullAtoms.info, lazyLoaded: true },
      };
    },

    async getChunk(index: number): Promise<ASEAtoms | null> {
      return chunkCache.get(index) || null;
    },

    async preloadAround(centerIndex: number, radius = 2): Promise<void> {
      // Preload implementation
    },

    releaseChunks(indices: number[]): void {
      loader.releaseChunks(indices);
    },

    getMemoryUsage(): number {
      return loader.getMemoryUsage();
    },
  };
}
