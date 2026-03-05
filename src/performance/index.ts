/**
 * Performance Module
 *
 * Provides performance optimization features for OpenQC-VSCode:
 * - Lazy loading for large structures
 * - WebWorker for background computations
 * - LRU caching with TTL support
 * - Incremental parsing with change detection
 */

export * from './lazyLoading';
export * from './computeWorker';
export * from './workerManager';
export * from './cacheManager';
export * from './incrementalParser';
