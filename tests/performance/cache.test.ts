/**
 * Unit tests for Cache module
 */

import {
  IncrementalCache,
  CacheEntry,
  CacheConfig,
  DEFAULT_CACHE_CONFIG,
  getGlobalCache,
  resetGlobalCache,
} from '../../src/performance/cache';

describe('IncrementalCache', () => {
  let cache: IncrementalCache;

  beforeEach(() => {
    cache = new IncrementalCache();
  });

  afterEach(() => {
    cache.clear();
  });

  describe('Basic Operations', () => {
    it('should set and get cache entries', () => {
      const data = { test: 'value', number: 42 };
      const hash = 'abc123';

      cache.set('key1', data, hash);
      const entry = cache.get('key1');

      expect(entry).toBeDefined();
      expect(entry?.data).toEqual(data);
      expect(entry?.hash).toBe(hash);
      expect(entry?.timestamp).toBeGreaterThan(0);
      expect(entry?.size).toBeGreaterThan(0);
    });

    it('should return undefined for non-existent keys', () => {
      const entry = cache.get('nonexistent');
      expect(entry).toBeUndefined();
    });

    it('should check if key exists', () => {
      cache.set('key1', 'data', 'hash1');

      expect(cache.has('key1')).toBe(true);
      expect(cache.has('key2')).toBe(false);
    });

    it('should delete entries', () => {
      cache.set('key1', 'data', 'hash1');

      expect(cache.delete('key1')).toBe(true);
      expect(cache.has('key1')).toBe(false);
      expect(cache.delete('key1')).toBe(false);
    });

    it('should clear all entries', () => {
      cache.set('key1', 'data1', 'hash1');
      cache.set('key2', 'data2', 'hash2');

      cache.clear();

      expect(cache.has('key1')).toBe(false);
      expect(cache.has('key2')).toBe(false);
      expect(cache.getEntryCount()).toBe(0);
    });
  });

  describe('LRU Eviction', () => {
    it('should evict least recently used entries when size limit exceeded', () => {
      const smallCache = new IncrementalCache({ maxSize: 100 });

      smallCache.set('key1', 'small', 'hash1');
      smallCache.set('key2', 'small', 'hash2');

      // Access key1 to make it more recently used
      smallCache.get('key1');

      // Add more data to trigger eviction
      smallCache.set('key3', 'larger data here', 'hash3');

      // key2 should be evicted (least recently used)
      expect(smallCache.has('key2')).toBe(false);
      expect(smallCache.has('key1')).toBe(true);
    });

    it('should update access order on get', () => {
      cache.set('key1', 'data1', 'hash1');
      cache.set('key2', 'data2', 'hash2');

      // Access key1 to update its position
      cache.get('key1');

      const keys = cache.keys();
      expect(keys[keys.length - 1]).toBe('key1');
    });
  });

  describe('Hash-based Validation', () => {
    it('should return data when hash matches', () => {
      const data = { test: 'value' };
      const hash = 'abc123';

      cache.set('key1', data, hash);
      const result = cache.getIfValid('key1', hash);

      expect(result).toEqual(data);
    });

    it('should return null when hash does not match', () => {
      cache.set('key1', 'data', 'hash1');
      const result = cache.getIfValid('key1', 'different_hash');

      expect(result).toBeNull();
    });

    it('should return null for non-existent key', () => {
      const result = cache.getIfValid('nonexistent', 'hash');
      expect(result).toBeNull();
    });
  });

  describe('Statistics', () => {
    it('should track hit and miss counts', () => {
      cache.set('key1', 'data', 'hash');

      // Miss
      cache.get('nonexistent');

      // Hit
      cache.get('key1');
      cache.get('key1');

      const stats = cache.getStats();
      expect(stats.hits).toBe(2);
      expect(stats.misses).toBe(1);
      expect(stats.hitRate).toBeCloseTo(66.67, 1);
    });

    it('should reset statistics', () => {
      cache.set('key1', 'data', 'hash');
      cache.get('key1');
      cache.get('nonexistent');

      cache.resetStats();

      const stats = cache.getStats();
      expect(stats.hits).toBe(0);
      expect(stats.misses).toBe(0);
      expect(stats.hitRate).toBe(0);
    });

    it('should report correct entry count and size', () => {
      cache.set('key1', 'data1', 'hash1');
      cache.set('key2', 'data2', 'hash2');

      const stats = cache.getStats();
      expect(stats.entryCount).toBe(2);
      expect(stats.currentSize).toBeGreaterThan(0);
      expect(stats.maxSize).toBe(DEFAULT_CACHE_CONFIG.maxSize);
    });
  });

  describe('Invalidation', () => {
    it('should invalidate entries by age', () => {
      cache.set('key1', 'data1', 'hash1');
      cache.set('key2', 'data2', 'hash2');

      // Wait a bit
      jest.advanceTimersByTime(100);

      const invalidated = cache.invalidateByAge(50);

      expect(invalidated).toBe(2);
      expect(cache.getEntryCount()).toBe(0);
    });

    it('should not invalidate recent entries', () => {
      cache.set('key1', 'data1', 'hash1');

      const invalidated = cache.invalidateByAge(10000);

      expect(invalidated).toBe(0);
      expect(cache.has('key1')).toBe(true);
    });

    it('should invalidate entries using custom validator', () => {
      cache.set('key1', 'data1', 'hash1');
      cache.set('key2', 'data2', 'hash2');

      const invalidated = cache.invalidateStale((key: string, entry: CacheEntry) => key !== 'key1');

      expect(invalidated).toBe(1);
      expect(cache.has('key1')).toBe(true);
      expect(cache.has('key2')).toBe(false);
    });
  });

  describe('Hash Calculation', () => {
    it('should calculate consistent hashes for strings', () => {
      const content = 'test content';
      const hash1 = IncrementalCache.calculateHash(content);
      const hash2 = IncrementalCache.calculateHash(content);

      expect(hash1).toBe(hash2);
      expect(hash1.length).toBe(64); // SHA-256 hex length
    });

    it('should calculate different hashes for different content', () => {
      const hash1 = IncrementalCache.calculateHash('content1');
      const hash2 = IncrementalCache.calculateHash('content2');

      expect(hash1).not.toBe(hash2);
    });

    it('should calculate hash for buffers', () => {
      const buffer = Buffer.from('test content');
      const hash = IncrementalCache.calculateHashBuffer(buffer);

      expect(hash.length).toBe(64);
    });
  });

  describe('Size Estimation', () => {
    it('should estimate size of primitive types', () => {
      cache.set('bool', true, 'hash');
      cache.set('num', 123, 'hash');
      cache.set('str', 'hello', 'hash');

      const stats = cache.getStats();
      expect(stats.currentSize).toBeGreaterThan(0);
    });

    it('should estimate size of objects', () => {
      cache.set('obj', { a: 1, b: 'test', c: true }, 'hash');

      const stats = cache.getStats();
      expect(stats.currentSize).toBeGreaterThan(0);
    });

    it('should estimate size of arrays', () => {
      cache.set('arr', [1, 2, 3, 4, 5], 'hash');

      const stats = cache.getStats();
      expect(stats.currentSize).toBeGreaterThan(0);
    });
  });

  describe('Configuration', () => {
    it('should use default configuration', () => {
      const defaultCache = new IncrementalCache();
      const stats = defaultCache.getStats();

      expect(stats.maxSize).toBe(DEFAULT_CACHE_CONFIG.maxSize);
    });

    it('should accept custom configuration', () => {
      const customCache = new IncrementalCache({
        maxSize: 1024,
        defaultTTL: 5000,
        enableStats: false,
      });

      const stats = customCache.getStats();
      expect(stats.maxSize).toBe(1024);
    });
  });
});

describe('Global Cache', () => {
  beforeEach(() => {
    resetGlobalCache();
  });

  afterEach(() => {
    resetGlobalCache();
  });

  it('should return same instance', () => {
    const cache1 = getGlobalCache();
    const cache2 = getGlobalCache();

    expect(cache1).toBe(cache2);
  });

  it('should create new instance after reset', () => {
    const cache1 = getGlobalCache();
    resetGlobalCache();
    const cache2 = getGlobalCache();

    expect(cache1).not.toBe(cache2);
  });

  it('should persist data across get calls', () => {
    const cache = getGlobalCache();
    cache.set('test', 'value', 'hash');

    const cache2 = getGlobalCache();
    expect(cache2.get('test')?.data).toBe('value');
  });
});
