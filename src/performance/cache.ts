/**
 * Cache System for Performance Optimization
 *
 * Implements LRU (Least Recently Used) cache with file content hashing
 * and incremental parsing support for quantum chemistry files.
 */

import * as crypto from 'crypto';

/**
 * Cache entry metadata and data
 */
export interface CacheEntry {
  /** Content hash for cache validation */
  hash: string;
  /** Cached data */
  data: any;
  /** Timestamp when entry was created/updated */
  timestamp: number;
  /** Size of cached data in bytes (approximate) */
  size: number;
}

/**
 * Cache statistics for monitoring
 */
export interface CacheStats {
  /** Total number of entries */
  entryCount: number;
  /** Current cache size in bytes */
  currentSize: number;
  /** Maximum allowed cache size */
  maxSize: number;
  /** Cache hit count */
  hits: number;
  /** Cache miss count */
  misses: number;
  /** Hit rate percentage */
  hitRate: number;
}

/**
 * Cache configuration options
 */
export interface CacheConfig {
  /** Maximum cache size in bytes (default: 50MB) */
  maxSize: number;
  /** Default TTL in milliseconds (default: 1 hour) */
  defaultTTL: number;
  /** Whether to enable cache statistics */
  enableStats: boolean;
}

/** Default cache configuration */
export const DEFAULT_CACHE_CONFIG: CacheConfig = {
  maxSize: 50 * 1024 * 1024, // 50MB
  defaultTTL: 60 * 60 * 1000, // 1 hour
  enableStats: true,
};

/**
 * LRU Cache with file content hashing
 *
 * Provides efficient caching with:
 * - LRU eviction policy
 * - Content-based cache invalidation
 * - Size-aware memory management
 * - Statistics tracking
 */
export class IncrementalCache {
  private cache: Map<string, CacheEntry>;
  private accessOrder: string[];
  private config: CacheConfig;
  private currentSize: number;
  private hits: number;
  private misses: number;

  constructor(config: Partial<CacheConfig> = {}) {
    this.config = { ...DEFAULT_CACHE_CONFIG, ...config };
    this.cache = new Map();
    this.accessOrder = [];
    this.currentSize = 0;
    this.hits = 0;
    this.misses = 0;
  }

  /**
   * Calculate hash of file content
   * @param content File content string
   * @returns SHA-256 hash
   */
  static calculateHash(content: string): string {
    return crypto.createHash('sha256').update(content).digest('hex');
  }

  /**
   * Calculate hash of buffer content
   * @param buffer Buffer content
   * @returns SHA-256 hash
   */
  static calculateHashBuffer(buffer: Buffer): string {
    return crypto.createHash('sha256').update(buffer).digest('hex');
  }

  /**
   * Get cached entry by key
   * @param key Cache key
   * @returns Cache entry or undefined
   */
  get(key: string): CacheEntry | undefined {
    const entry = this.cache.get(key);

    if (entry) {
      // Update access order (LRU)
      this.updateAccessOrder(key);
      this.hits++;
      return entry;
    }

    this.misses++;
    return undefined;
  }

  /**
   * Get data from cache if hash matches
   * @param key Cache key
   * @param contentHash Expected content hash
   * @returns Cached data or null if stale/missing
   */
  getIfValid(key: string, contentHash: string): any | null {
    const entry = this.cache.get(key);

    if (entry && entry.hash === contentHash) {
      // Update access order (LRU)
      this.updateAccessOrder(key);
      this.hits++;
      return entry.data;
    }

    this.misses++;
    return null;
  }

  /**
   * Set cache entry
   * @param key Cache key
   * @param data Data to cache
   * @param contentHash Content hash for validation
   */
  set(key: string, data: any, contentHash: string): void {
    const size = this.estimateSize(data);

    // Evict entries if necessary
    while (this.currentSize + size > this.config.maxSize && this.accessOrder.length > 0) {
      this.evictLRU();
    }

    const entry: CacheEntry = {
      hash: contentHash,
      data,
      timestamp: Date.now(),
      size,
    };

    // Remove old size from current size if updating
    const oldEntry = this.cache.get(key);
    if (oldEntry) {
      this.currentSize -= oldEntry.size;
      // Remove from access order to re-add at end
      const index = this.accessOrder.indexOf(key);
      if (index >= 0) {
        this.accessOrder.splice(index, 1);
      }
    }

    this.cache.set(key, entry);
    this.accessOrder.push(key);
    this.currentSize += size;
  }

  /**
   * Check if key exists in cache
   * @param key Cache key
   * @returns True if exists
   */
  has(key: string): boolean {
    return this.cache.has(key);
  }

  /**
   * Delete entry from cache
   * @param key Cache key
   * @returns True if deleted
   */
  delete(key: string): boolean {
    const entry = this.cache.get(key);
    if (entry) {
      this.currentSize -= entry.size;
      this.cache.delete(key);
      const index = this.accessOrder.indexOf(key);
      if (index >= 0) {
        this.accessOrder.splice(index, 1);
      }
      return true;
    }
    return false;
  }

  /**
   * Clear all cache entries
   */
  clear(): void {
    this.cache.clear();
    this.accessOrder = [];
    this.currentSize = 0;
  }

  /**
   * Invalidate stale entries based on hash validation function
   * @param validator Function to check if entry is valid
   * @returns Number of entries invalidated
   */
  invalidateStale(validator: (key: string, entry: CacheEntry) => boolean): number {
    let invalidated = 0;

    for (const [key, entry] of this.cache.entries()) {
      if (!validator(key, entry)) {
        this.delete(key);
        invalidated++;
      }
    }

    return invalidated;
  }

  /**
   * Invalidate entries older than specified age
   * @param maxAge Maximum age in milliseconds
   * @returns Number of entries invalidated
   */
  invalidateByAge(maxAge: number): number {
    const now = Date.now();
    let invalidated = 0;

    for (const [key, entry] of this.cache.entries()) {
      if (now - entry.timestamp > maxAge) {
        this.delete(key);
        invalidated++;
      }
    }

    return invalidated;
  }

  /**
   * Get cache statistics
   */
  getStats(): CacheStats {
    const total = this.hits + this.misses;
    return {
      entryCount: this.cache.size,
      currentSize: this.currentSize,
      maxSize: this.config.maxSize,
      hits: this.hits,
      misses: this.misses,
      hitRate: total > 0 ? (this.hits / total) * 100 : 0,
    };
  }

  /**
   * Reset statistics counters
   */
  resetStats(): void {
    this.hits = 0;
    this.misses = 0;
  }

  /**
   * Get all keys in cache
   */
  keys(): string[] {
    return Array.from(this.cache.keys());
  }

  /**
   * Get current cache size
   */
  getCurrentSize(): number {
    return this.currentSize;
  }

  /**
   * Get number of entries
   */
  getEntryCount(): number {
    return this.cache.size;
  }

  /**
   * Update access order for LRU tracking
   */
  private updateAccessOrder(key: string): void {
    const index = this.accessOrder.indexOf(key);
    if (index >= 0) {
      this.accessOrder.splice(index, 1);
    }
    this.accessOrder.push(key);
  }

  /**
   * Evict least recently used entry
   */
  private evictLRU(): void {
    if (this.accessOrder.length === 0) {
      return;
    }

    const lruKey = this.accessOrder.shift()!;
    const entry = this.cache.get(lruKey);

    if (entry) {
      this.currentSize -= entry.size;
      this.cache.delete(lruKey);
    }
  }

  /**
   * Estimate size of data in bytes
   */
  private estimateSize(data: any): number {
    if (data === null || data === undefined) {
      return 0;
    }

    const type = typeof data;

    switch (type) {
      case 'boolean':
        return 4;
      case 'number':
        return 8;
      case 'string':
        return data.length * 2;
      case 'object':
        if (Array.isArray(data)) {
          return data.reduce((sum, item) => sum + this.estimateSize(item), 0);
        }
        if (Buffer.isBuffer(data)) {
          return data.length;
        }
        return Object.keys(data).reduce((sum, key) => {
          return sum + key.length * 2 + this.estimateSize(data[key]);
        }, 0);
      default:
        return 0;
    }
  }
}

/**
 * Singleton cache instance for global use
 */
let globalCache: IncrementalCache | null = null;

/**
 * Get or create global cache instance
 */
export function getGlobalCache(): IncrementalCache {
  if (!globalCache) {
    globalCache = new IncrementalCache();
  }
  return globalCache;
}

/**
 * Reset global cache instance
 */
export function resetGlobalCache(): void {
  globalCache = null;
}
