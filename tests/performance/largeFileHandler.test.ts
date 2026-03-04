/**
 * Large File Handler Tests
 */

import {
  ChunkedFileLoader,
  MemoryPool,
  VirtualScrollHandler,
  LargeFileHandler,
  LargeFileOptions,
  LoadProgress,
} from '../../src/performance/largeFileHandler';

describe('ChunkedFileLoader', () => {
  let loader: ChunkedFileLoader;

  beforeEach(() => {
    loader = new ChunkedFileLoader({
      chunkSize: 10,
      maxChunksInMemory: 3,
      preloadDistance: 1,
    });
  });

  describe('basic loading', () => {
    test('should initialize with correct configuration', () => {
      expect(loader).toBeDefined();
      const stats = loader.getStats();
      expect(stats.totalChunks).toBe(0);
    });

    test('should load file and emit events', async () => {
      const loadStartSpy = jest.fn();
      const loadCompleteSpy = jest.fn();
      
      loader.on('loadStart', loadStartSpy);
      loader.on('loadComplete', loadCompleteSpy);

      const mockParser = jest.fn().mockResolvedValue([
        { element: 'H', x: 0, y: 0, z: 0 },
      ]);

      const content = 'line1\\nline2\\nline3';
      await loader.loadFile('/test/file.xyz', mockParser, content);

      expect(loadStartSpy).toHaveBeenCalled();
      expect(loadCompleteSpy).toHaveBeenCalled();
    });

    test('should load atoms in range', async () => {
      const mockParser = jest.fn().mockImplementation((_, start, end) => {
        const atoms = [];
        for (let i = start; i < end; i++) {
          atoms.push({ element: 'H', index: i });
        }
        return Promise.resolve(atoms);
      });

      const content = Array(50).fill('H 0 0 0').join('\\n');
      await loader.loadFile('/test/file.xyz', mockParser, content);

      const atoms = await loader.getAtomsInRange(5, 15, mockParser, content);
      expect(atoms.length).toBe(10);
      expect(atoms[0].index).toBe(5);
      expect(atoms[9].index).toBe(14);
    });
  });

  describe('memory management', () => {
    test('should evict old chunks when max exceeded', async () => {
      const evictedSpy = jest.fn();
      loader.on('chunkEvicted', evictedSpy);

      const mockParser = jest.fn().mockImplementation((_, start, end) => {
        return Promise.resolve(Array(end - start).fill({ element: 'H' }));
      });

      const content = Array(100).fill('H 0 0 0').join('\\n');
      await loader.loadFile('/test/file.xyz', mockParser, content);

      // Load multiple chunks to trigger eviction
      await loader.getAtomsInRange(0, 10, mockParser, content);
      await loader.getAtomsInRange(10, 20, mockParser, content);
      await loader.getAtomsInRange(20, 30, mockParser, content);
      await loader.getAtomsInRange(30, 40, mockParser, content);

      expect(evictedSpy).toHaveBeenCalled();
    });

    test('should provide accurate statistics', async () => {
      const mockParser = jest.fn().mockResolvedValue([{ element: 'H' }]);
      const content = Array(20).fill('H 0 0 0').join('\\n');
      
      await loader.loadFile('/test/file.xyz', mockParser, content);
      
      const stats = loader.getStats();
      expect(stats.totalAtoms).toBeGreaterThan(0);
      expect(stats.totalChunks).toBeGreaterThan(0);
      expect(typeof stats.memoryUsage).toBe('number');
    });
  });

  describe('progress tracking', () => {
    test('should report progress during loading', async () => {
      const progressSpy = jest.fn();
      loader.on('progress', progressSpy);

      const mockParser = jest.fn().mockResolvedValue([{ element: 'H' }]);
      const content = Array(30).fill('H 0 0 0').join('\\n');
      
      await loader.loadFile('/test/file.xyz', mockParser, content);

      expect(progressSpy).toHaveBeenCalled();
      const progress: LoadProgress = progressSpy.mock.calls[0][0];
      expect(progress.percentage).toBeGreaterThanOrEqual(0);
      expect(progress.percentage).toBeLessThanOrEqual(100);
    });
  });

  describe('clear', () => {
    test('should clear all loaded chunks', async () => {
      const mockParser = jest.fn().mockResolvedValue([{ element: 'H' }]);
      const content = 'test content';
      
      await loader.loadFile('/test/file.xyz', mockParser, content);
      loader.clear();
      
      const stats = loader.getStats();
      expect(stats.chunksInMemory).toBe(0);
    });
  });
});

describe('MemoryPool', () => {
  let pool: MemoryPool;

  beforeEach(() => {
    pool = new MemoryPool({
      initialSize: 100,
      maxSize: 1000,
      blockSize: 10,
    });
  });

  test('should acquire and release blocks', () => {
    const block = (pool as any).acquire(5);
    expect(block).not.toBeNull();
    expect(block!.length).toBe(5);

    (pool as any).release(block!);
    
    const stats = (pool as any).getStats();
    expect(stats.usedBlocks).toBe(0);
    expect(stats.freeBlocks).toBeGreaterThan(0);
  });

  test('should return null when pool exhausted', () => {
    // Try to acquire more than max
    for (let i = 0; i < 200; i++) {
      (pool as any).acquire(10);
    }

    const block = (pool as any).acquire(10);
    expect(block).toBeNull();
  });

  test('should provide accurate statistics', () => {
    const block1 = (pool as any).acquire(5);
    const block2 = (pool as any).acquire(5);

    const stats = (pool as any).getStats();
    expect(stats.usedBlocks).toBe(2);
    expect(stats.totalBlocks).toBeGreaterThanOrEqual(10);
    expect(stats.utilization).toBeGreaterThan(0);

    (pool as any).release(block1!);
    (pool as any).release(block2!);
  });

  test('should clear all blocks', () => {
    (pool as any).acquire(5);
    (pool as any).acquire(5);
    
    (pool as any).clear();
    
    const stats = (pool as any).getStats();
    expect(stats.usedBlocks).toBe(0);
  });
});

describe('VirtualScrollHandler', () => {
  let handler: VirtualScrollHandler;

  beforeEach(() => {
    handler = new VirtualScrollHandler(200, 20); // 200px viewport, 20px per item
  });

  test('should calculate visible range', () => {
    const range = handler.calculateVisibleRange(100, 1000);
    
    // At 100px scroll with 20px items, start should be 5
    expect(range.start).toBe(5);
    expect(range.end).toBeGreaterThan(range.start);
    expect(range.end).toBeLessThanOrEqual(1000);
  });

  test('should handle scroll at top', () => {
    const range = handler.calculateVisibleRange(0, 100);
    expect(range.start).toBe(0);
  });

  test('should handle scroll at bottom', () => {
    const range = handler.calculateVisibleRange(2000, 100);
    expect(range.end).toBe(100);
  });

  test('should update viewport size', () => {
    handler.setViewportSize(400);
    const range = handler.calculateVisibleRange(0, 100);
    
    // With 400px viewport, should show more items
    expect(range.end - range.start).toBeGreaterThan(10);
  });

  test('should return current visible range', () => {
    handler.calculateVisibleRange(50, 100);
    const range = handler.getVisibleRange();
    
    expect(range.start).toBeGreaterThanOrEqual(0);
    expect(range.end).toBeGreaterThan(range.start);
  });
});

describe('LargeFileHandler', () => {
  let handler: LargeFileHandler;

  beforeEach(() => {
    const options: LargeFileOptions = {
      chunkSize: 50,
      maxChunksInMemory: 5,
    };
    handler = new LargeFileHandler(options);
  });

  test('should initialize all components', () => {
    expect(handler).toBeDefined();
    expect(typeof handler.getStats).toBe('function');
  });

  test('should provide combined statistics', () => {
    const stats = handler.getStats();
    
    expect(stats).toBeDefined();
    expect(stats.chunkLoader).toBeDefined();
    expect(stats.memoryPool).toBeDefined();
    expect(stats.visibleRange).toBeDefined();
    
    expect(typeof stats.chunkLoader.totalAtoms).toBe('number');
    expect(typeof stats.memoryPool.totalBlocks).toBe('number');
  });
});
