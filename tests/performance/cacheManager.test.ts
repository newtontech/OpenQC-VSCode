/**
 * Cache Manager Tests
 */

import { LRUCache, CacheManager, CacheKeyGenerator, CacheStats } from '../../src/performance/cacheManager';

describe('LRUCache', () => {
  let cache: LRUCache<string>;

  beforeEach(() => {
    cache = new LRUCache<string>({ maxSize: 1000, maxEntries: 5 });
  });

  afterEach(() => {
    cache.dispose();
  });

  describe('Basic Operations', () => {
    it('should set and get values', () => {
      cache.set('key1', 'value1');
      expect(cache.get('key1')).toBe('value1');
    });

    it('should return undefined for missing keys', () => {
      expect(cache.get('nonexistent')).toBeUndefined();
    });

    it('should check if key exists', () => {
      cache.set('key1', 'value1');
      expect(cache.has('key1')).toBe(true);
      expect(cache.has('nonexistent')).toBe(false);
    });

    it('should delete entries', () => {
      cache.set('key1', 'value1');
      expect(cache.delete('key1')).toBe(true);
      expect(cache.get('key1')).toBeUndefined();
    });

    it('should clear all entries', () => {
      cache.set('key1', 'value1');
      cache.set('key2', 'value2');
      cache.clear();
      expect(cache.size()).toBe(0);
    });
  });

  describe('LRU Eviction', () => {
    it('should evict least recently used when max entries reached', () => {
      // Add 5 entries (max)
      for (let i = 0; i < 5; i++) {
        cache.set(`key${i}`, `value${i}`);
      }

      // Add one more, should evict first
      cache.set('key5', 'value5');

      expect(cache.has('key0')).toBe(false);
      expect(cache.has('key5')).toBe(true);
    });

    it('should update access order on get', () => {
      // Add 5 entries
      for (let i = 0; i < 5; i++) {
        cache.set(`key${i}`, `value${i}`);
      }

      // Access key0 to make it recently used
      cache.get('key0');

      // Add one more, should evict key1 (not key0)
      cache.set('key5', 'value5');

      expect(cache.has('key0')).toBe(true);
      expect(cache.has('key1')).toBe(false);
    });
  });

  describe('TTL Support', () => {
    it('should expire entries after TTL', (done) => {
      const shortCache = new LRUCache<string>({
        maxSize: 1000,
        defaultTTL: 100, // 100ms
        cleanupInterval: 50,
      });

      shortCache.set('key1', 'value1');

      // Should be present immediately
      expect(shortCache.get('key1')).toBe('value1');

      // Should be expired after TTL
      setTimeout(() => {
        expect(shortCache.get('key1')).toBeUndefined();
        shortCache.dispose();
        done();
      }, 150);
    });

    it('should accept custom TTL per entry', () => {
      cache.set('key1', 'value1', 5000); // 5 seconds
      cache.set('key2', 'value2', 100); // 100ms

      expect(cache.get('key1')).toBe('value1');
      expect(cache.get('key2')).toBe('value2');
    });
  });

  describe('Statistics', () => {
    it('should track cache hits and misses', () => {
      cache.set('key1', 'value1');

      cache.get('key1'); // Hit
      cache.get('key1'); // Hit
      cache.get('nonexistent'); // Miss

      const stats = cache.getStats();
      expect(stats.hits).toBe(2);
      expect(stats.misses).toBe(1);
    });

    it('should calculate hit rate', () => {
      cache.set('key1', 'value1');

      cache.get('key1'); // Hit
      cache.get('nonexistent'); // Miss

      const stats = cache.getStats();
      expect(stats.hitRate).toBe(0.5); // 1 hit, 1 miss = 50%
    });

    it('should track evictions', () => {
      // Fill cache
      for (let i = 0; i < 5; i++) {
        cache.set(`key${i}`, `value${i}`);
      }

      // Trigger eviction
      cache.set('key5', 'value5');

      const stats = cache.getStats();
      expect(stats.evictions).toBe(1);
    });

    it('should track cache size', () => {
      cache.set('key1', 'value1');
      cache.set('key2', 'value2');

      const stats = cache.getStats();
      expect(stats.entries).toBe(2);
      expect(stats.size).toBeGreaterThan(0);
    });
  });
});

describe('CacheKeyGenerator', () => {
  describe('forFile', () => {
    it('should generate unique keys for different files', () => {
      const key1 = CacheKeyGenerator.forFile('/path/to/file1.txt', 'parse');
      const key2 = CacheKeyGenerator.forFile('/path/to/file2.txt', 'parse');
      expect(key1).not.toBe(key2);
    });

    it('should generate same key for same file and operation', () => {
      const key1 = CacheKeyGenerator.forFile('/path/to/file.txt', 'parse');
      const key2 = CacheKeyGenerator.forFile('/path/to/file.txt', 'parse');
      expect(key1).toBe(key2);
    });
  });

  describe('forConversion', () => {
    it('should generate unique keys for different conversions', () => {
      const key1 = CacheKeyGenerator.forConversion('/path/to/file.txt', 'vasp', 'cp2k');
      const key2 = CacheKeyGenerator.forConversion('/path/to/file.txt', 'vasp', 'qe');
      expect(key1).not.toBe(key2);
    });
  });

  describe('forValidation', () => {
    it('should generate unique keys for different validation checks', () => {
      const key1 = CacheKeyGenerator.forValidation('/path/to/file.txt', ['bond_lengths']);
      const key2 = CacheKeyGenerator.forValidation('/path/to/file.txt', ['atom_overlap']);
      expect(key1).not.toBe(key2);
    });
  });

  describe('forProperties', () => {
    it('should generate unique keys for different properties', () => {
      const key1 = CacheKeyGenerator.forProperties('/path/to/file.txt', ['center_of_mass']);
      const key2 = CacheKeyGenerator.forProperties('/path/to/file.txt', ['moment_of_inertia']);
      expect(key1).not.toBe(key2);
    });
  });
});

describe('CacheManager', () => {
  let manager: CacheManager;

  beforeEach(() => {
    manager = CacheManager.getInstance();
  });

  afterEach(() => {
    manager.clearAll();
    manager.dispose();
  });

  describe('Singleton Pattern', () => {
    it('should return same instance', () => {
      const instance1 = CacheManager.getInstance();
      const instance2 = CacheManager.getInstance();
      expect(instance1).toBe(instance2);
    });
  });

  describe('Multiple Caches', () => {
    it('should provide separate caches for different types', () => {
      const structureCache = manager.getStructureCache();
      const conversionCache = manager.getConversionCache();

      structureCache.set('key1', { data: 'structure' });
      conversionCache.set('key1', 'converted data');

      expect(structureCache.get('key1')).toEqual({ data: 'structure' });
      expect(conversionCache.get('key1')).toBe('converted data');
    });
  });

  describe('Statistics', () => {
    it('should provide combined statistics', () => {
      const structureCache = manager.getStructureCache();
      structureCache.set('key1', { data: 'test' });
      structureCache.get('key1');

      const allStats = manager.getAllStats();
      expect(allStats.structure).toBeDefined();
      expect(allStats.conversion).toBeDefined();
      expect(allStats.validation).toBeDefined();
      expect(allStats.property).toBeDefined();
    });
  });

  describe('Clear All', () => {
    it('should clear all caches', () => {
      manager.getStructureCache().set('key1', 'value1');
      manager.getConversionCache().set('key2', 'value2');

      manager.clearAll();

      expect(manager.getStructureCache().size()).toBe(0);
      expect(manager.getConversionCache().size()).toBe(0);
    });
  });
});
