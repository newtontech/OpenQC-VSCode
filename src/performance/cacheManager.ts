/**
 * Cache Manager - Intelligent Caching System for OpenQC-VSCode
 *
 * Provides LRU caching with TTL support for:
 * - Parsed structures
 * - Converted formats
 * - Validation results
 * - Computed properties
 */

import * as vscode from 'vscode';
import * as crypto from 'crypto';
import * as fs from 'fs';
import * as path from 'path';

/**
 * Cache entry with metadata
 */
export interface CacheEntry<T> {
  data: T;
  timestamp: number;
  ttl: number;
  accessCount: number;
  lastAccess: number;
  size: number;
  hash: string;
}

/**
 * Cache statistics
 */
export interface CacheStats {
  hits: number;
  misses: number;
  evictions: number;
  size: number;
  maxSize: number;
  entries: number;
  hitRate: number;
}

/**
 * Cache configuration
 */
export interface CacheConfig {
  maxSize: number; // Maximum cache size in bytes
  maxEntries: number; // Maximum number of entries
  defaultTTL: number; // Default time-to-live in milliseconds
  cleanupInterval: number; // Cleanup interval in milliseconds
  persistToDisk: boolean; // Whether to persist cache to disk
  cacheDir?: string; // Directory for persistent cache
}

export const DEFAULT_CACHE_CONFIG: CacheConfig = {
  maxSize: 100 * 1024 * 1024, // 100 MB
  maxEntries: 1000,
  defaultTTL: 30 * 60 * 1000, // 30 minutes
  cleanupInterval: 5 * 60 * 1000, // 5 minutes
  persistToDisk: true,
};

/**
 * LRU Cache with TTL support
 */
export class LRUCache<T> {
  private cache: Map<string, CacheEntry<T>> = new Map();
  private config: CacheConfig;
  private stats: CacheStats = {
    hits: 0,
    misses: 0,
    evictions: 0,
    size: 0,
    maxSize: 0,
    entries: 0,
    hitRate: 0,
  };
  private cleanupTimer?: NodeJS.Timeout;
  private context?: vscode.ExtensionContext;

  constructor(config: Partial<CacheConfig> = {}, context?: vscode.ExtensionContext) {
    this.config = { ...DEFAULT_CACHE_CONFIG, ...config };
    this.stats.maxSize = this.config.maxSize;
    this.context = context;

    if (this.config.persistToDisk && context) {
      this.loadFromDisk();
    }

    this.startCleanup();
  }

  /**
   * Get item from cache
   */
  get(key: string): T | undefined {
    const entry = this.cache.get(key);

    if (!entry) {
      this.stats.misses++;
      this.updateHitRate();
      return undefined;
    }

    // Check TTL
    if (Date.now() - entry.timestamp > entry.ttl) {
      this.delete(key);
      this.stats.misses++;
      this.updateHitRate();
      return undefined;
    }

    // Update access metadata
    entry.accessCount++;
    entry.lastAccess = Date.now();

    // Move to end (most recently used)
    this.cache.delete(key);
    this.cache.set(key, entry);

    this.stats.hits++;
    this.updateHitRate();

    return entry.data;
  }

  /**
   * Set item in cache
   */
  set(key: string, value: T, ttl?: number): void {
    const size = this.estimateSize(value);
    const entryTTL = ttl || this.config.defaultTTL;

    // Check if we need to evict entries
    while (
      (this.stats.size + size > this.config.maxSize || this.cache.size >= this.config.maxEntries) &&
      this.cache.size > 0
    ) {
      this.evictLRU();
    }

    // Remove old entry if exists
    if (this.cache.has(key)) {
      this.delete(key);
    }

    const entry: CacheEntry<T> = {
      data: value,
      timestamp: Date.now(),
      ttl: entryTTL,
      accessCount: 0,
      lastAccess: Date.now(),
      size,
      hash: this.hashValue(value),
    };

    this.cache.set(key, entry);
    this.stats.size += size;
    this.stats.entries = this.cache.size;
  }

  /**
   * Delete item from cache
   */
  delete(key: string): boolean {
    const entry = this.cache.get(key);
    if (!entry) {
      return false;
    }

    this.cache.delete(key);
    this.stats.size -= entry.size;
    this.stats.entries = this.cache.size;

    return true;
  }

  /**
   * Check if key exists and is valid
   */
  has(key: string): boolean {
    const entry = this.cache.get(key);
    if (!entry) {
      return false;
    }

    // Check TTL
    if (Date.now() - entry.timestamp > entry.ttl) {
      this.delete(key);
      return false;
    }

    return true;
  }

  /**
   * Clear all entries
   */
  clear(): void {
    this.cache.clear();
    this.stats.size = 0;
    this.stats.entries = 0;
  }

  /**
   * Get cache statistics
   */
  getStats(): CacheStats {
    return { ...this.stats };
  }

  /**
   * Get all keys
   */
  keys(): string[] {
    return Array.from(this.cache.keys());
  }

  /**
   * Get cache size
   */
  size(): number {
    return this.cache.size;
  }

  /**
   * Persist cache to disk
   */
  async saveToDisk(): Promise<void> {
    if (!this.config.persistToDisk || !this.context) {
      return;
    }

    try {
      const cacheDir = this.config.cacheDir || path.join(this.context.globalStorageUri.fsPath, 'cache');
      await fs.promises.mkdir(cacheDir, { recursive: true });

      const cacheFile = path.join(cacheDir, 'openqc-cache.json');
      const data: Record<string, CacheEntry<T>> = {};

      this.cache.forEach((entry, key) => {
        data[key] = entry;
      });

      await fs.promises.writeFile(cacheFile, JSON.stringify(data, null, 2));
    } catch (error) {
      console.error('Failed to save cache to disk:', error);
    }
  }

  /**
   * Load cache from disk
   */
  private async loadFromDisk(): Promise<void> {
    if (!this.config.persistToDisk || !this.context) {
      return;
    }

    try {
      const cacheDir = this.config.cacheDir || path.join(this.context.globalStorageUri.fsPath, 'cache');
      const cacheFile = path.join(cacheDir, 'openqc-cache.json');

      if (!fs.existsSync(cacheFile)) {
        return;
      }

      const content = await fs.promises.readFile(cacheFile, 'utf-8');
      const data: Record<string, CacheEntry<T>> = JSON.parse(content);

      // Load entries that haven't expired
      const now = Date.now();
      for (const [key, entry] of Object.entries(data)) {
        if (now - entry.timestamp <= entry.ttl) {
          this.cache.set(key, entry);
          this.stats.size += entry.size;
        }
      }

      this.stats.entries = this.cache.size;
    } catch (error) {
      console.error('Failed to load cache from disk:', error);
    }
  }

  /**
   * Evict least recently used entry
   */
  private evictLRU(): void {
    // First entry is least recently used (Map maintains insertion order)
    const firstKey = this.cache.keys().next().value;
    if (firstKey) {
      this.delete(firstKey);
      this.stats.evictions++;
    }
  }

  /**
   * Start periodic cleanup
   */
  private startCleanup(): void {
    this.cleanupTimer = setInterval(() => {
      this.cleanup();
    }, this.config.cleanupInterval);
  }

  /**
   * Cleanup expired entries
   */
  private cleanup(): void {
    const now = Date.now();
    const keysToDelete: string[] = [];

    this.cache.forEach((entry, key) => {
      if (now - entry.timestamp > entry.ttl) {
        keysToDelete.push(key);
      }
    });

    keysToDelete.forEach(key => this.delete(key));

    // Save to disk after cleanup
    if (this.config.persistToDisk) {
      this.saveToDisk();
    }
  }

  /**
   * Update hit rate statistic
   */
  private updateHitRate(): void {
    const total = this.stats.hits + this.stats.misses;
    this.stats.hitRate = total > 0 ? this.stats.hits / total : 0;
  }

  /**
   * Estimate size of value in bytes
   */
  private estimateSize(value: T): number {
    try {
      const json = JSON.stringify(value);
      return json.length * 2; // UTF-16 encoding
    } catch {
      return 1024; // Default 1KB
    }
  }

  /**
   * Hash value for integrity check
   */
  private hashValue(value: T): string {
    try {
      const json = JSON.stringify(value);
      return crypto.createHash('md5').update(json).digest('hex');
    } catch {
      return '';
    }
  }

  /**
   * Dispose resources
   */
  dispose(): void {
    if (this.cleanupTimer) {
      clearInterval(this.cleanupTimer);
    }
    this.saveToDisk();
  }
}

/**
 * File-based cache key generator
 */
export class CacheKeyGenerator {
  /**
   * Generate cache key for file content
   */
  static forFile(filepath: string, operation: string): string {
    try {
      const stats = fs.statSync(filepath);
      const content = `${filepath}:${stats.mtimeMs}:${operation}`;
      return crypto.createHash('sha256').update(content).digest('hex');
    } catch (error) {
      // Fallback for tests where file doesn't exist
      const content = `${filepath}:${operation}`;
      return crypto.createHash('sha256').update(content).digest('hex');
    }
  }

  /**
   * Generate cache key for structure conversion
   */
  static forConversion(
    filepath: string,
    sourceFormat: string,
    targetFormat: string
  ): string {
    return this.forFile(filepath, `convert:${sourceFormat}:${targetFormat}`);
  }

  /**
   * Generate cache key for validation
   */
  static forValidation(filepath: string, checks: string[]): string {
    return this.forFile(filepath, `validate:${checks.join(',')}`);
  }

  /**
   * Generate cache key for property calculation
   */
  static forProperties(filepath: string, properties: string[]): string {
    return this.forFile(filepath, `properties:${properties.join(',')}`);
  }
}

/**
 * Global cache instances
 */
export class CacheManager {
  private static instance: CacheManager;
  private structureCache: LRUCache<any>;
  private conversionCache: LRUCache<string>;
  private validationCache: LRUCache<any>;
  private propertyCache: LRUCache<any>;
  private context?: vscode.ExtensionContext;

  private constructor(context?: vscode.ExtensionContext) {
    this.context = context;

    this.structureCache = new LRUCache<any>(
      { maxSize: 50 * 1024 * 1024, maxEntries: 200 },
      context
    );

    this.conversionCache = new LRUCache<string>(
      { maxSize: 20 * 1024 * 1024, maxEntries: 500 },
      context
    );

    this.validationCache = new LRUCache<any>(
      { maxSize: 10 * 1024 * 1024, maxEntries: 300 },
      context
    );

    this.propertyCache = new LRUCache<any>(
      { maxSize: 10 * 1024 * 1024, maxEntries: 300 },
      context
    );
  }

  static getInstance(context?: vscode.ExtensionContext): CacheManager {
    if (!CacheManager.instance) {
      CacheManager.instance = new CacheManager(context);
    }
    return CacheManager.instance;
  }

  getStructureCache(): LRUCache<any> {
    return this.structureCache;
  }

  getConversionCache(): LRUCache<string> {
    return this.conversionCache;
  }

  getValidationCache(): LRUCache<any> {
    return this.validationCache;
  }

  getPropertyCache(): LRUCache<any> {
    return this.propertyCache;
  }

  /**
   * Clear all caches
   */
  clearAll(): void {
    this.structureCache.clear();
    this.conversionCache.clear();
    this.validationCache.clear();
    this.propertyCache.clear();
  }

  /**
   * Get combined statistics
   */
  getAllStats(): Record<string, CacheStats> {
    return {
      structure: this.structureCache.getStats(),
      conversion: this.conversionCache.getStats(),
      validation: this.validationCache.getStats(),
      property: this.propertyCache.getStats(),
    };
  }

  /**
   * Dispose all caches
   */
  dispose(): void {
    this.structureCache.dispose();
    this.conversionCache.dispose();
    this.validationCache.dispose();
    this.propertyCache.dispose();
  }
}
