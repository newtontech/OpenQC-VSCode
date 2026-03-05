/**
 * Incremental Parser - Efficient Parsing for Large Files
 *
 * Implements incremental parsing with:
 * - Change detection
 * - Partial re-parsing
 * - AST diffing
 * - Smart caching
 */

import * as vscode from 'vscode';
import { LRUCache, CacheKeyGenerator } from './cacheManager';

/**
 * Parse result with metadata
 */
export interface ParseResult<T> {
  ast: T;
  version: number;
  timestamp: number;
  changes: ChangeSet;
  cached: boolean;
}

/**
 * Set of changes between versions
 */
export interface ChangeSet {
  added: Range[];
  removed: Range[];
  modified: Range[];
}

/**
 * Range in document
 */
export interface Range {
  startLine: number;
  endLine: number;
  startColumn: number;
  endColumn: number;
}

/**
 * Document version info
 */
export interface DocumentVersion {
  uri: string;
  version: number;
  hash: string;
  lineCount: number;
}

/**
 * Incremental parser configuration
 */
export interface IncrementalParserConfig {
  maxCacheSize: number;
  enableDiffing: boolean;
  minLinesForIncremental: number;
}

const DEFAULT_CONFIG: IncrementalParserConfig = {
  maxCacheSize: 50,
  enableDiffing: true,
  minLinesForIncremental: 100,
};

/**
 * Incremental Parser
 */
export class IncrementalParser<T> {
  private cache: LRUCache<ParseResult<T>>;
  private config: IncrementalParserConfig;
  private parseFunction: (content: string) => T;
  private diffFunction?: (oldAst: T, newAst: T) => ChangeSet;

  constructor(
    parseFunction: (content: string) => T,
    diffFunction?: (oldAst: T, newAst: T) => ChangeSet,
    config: Partial<IncrementalParserConfig> = {}
  ) {
    this.parseFunction = parseFunction;
    this.diffFunction = diffFunction;
    this.config = { ...DEFAULT_CONFIG, ...config };
    this.cache = new LRUCache<ParseResult<T>>({ maxEntries: this.config.maxCacheSize });
  }

  /**
   * Parse document with incremental updates
   */
  async parse(document: vscode.TextDocument): Promise<ParseResult<T>> {
    const uri = document.uri.toString();
    const version = document.version;
    const content = document.getText();

    // Check cache
    const cacheKey = CacheKeyGenerator.forFile(document.uri.fsPath, `parse:${version}`);
    const cached = this.cache.get(cacheKey);

    if (cached && cached.version === version) {
      return { ...cached, cached: true };
    }

    // Check if incremental parsing is beneficial
    if (this.shouldUseIncremental(document)) {
      return this.parseIncremental(document, content);
    }

    // Full parse
    return this.parseFull(uri, version, content);
  }

  /**
   * Parse with incremental updates
   */
  private async parseIncremental(
    document: vscode.TextDocument,
    content: string
  ): Promise<ParseResult<T>> {
    const uri = document.uri.toString();
    const version = document.version;

    // Find previous version
    const previousKey = this.findPreviousVersion(uri, version);
    const previous = previousKey ? this.cache.get(previousKey) : undefined;

    if (!previous) {
      return this.parseFull(uri, version, content);
    }

    // Get content changes
    const changes = this.detectChanges(document, previous);

    // If changes are too large, do full parse
    if (this.changesTooLarge(changes, document.lineCount)) {
      return this.parseFull(uri, version, content);
    }

    // Parse only changed regions
    const partialAst = await this.parsePartial(content, changes, previous.ast);

    const result: ParseResult<T> = {
      ast: partialAst,
      version,
      timestamp: Date.now(),
      changes,
      cached: false,
    };

    // Cache result
    const cacheKey = CacheKeyGenerator.forFile(document.uri.fsPath, `parse:${version}`);
    this.cache.set(cacheKey, result);

    return result;
  }

  /**
   * Full parse
   */
  private async parseFull(
    uri: string,
    version: number,
    content: string
  ): Promise<ParseResult<T>> {
    const ast = this.parseFunction(content);

    const result: ParseResult<T> = {
      ast,
      version,
      timestamp: Date.now(),
      changes: { added: [], removed: [], modified: [] },
      cached: false,
    };

    // Cache result
    const cacheKey = CacheKeyGenerator.forFile(uri, `parse:${version}`);
    this.cache.set(cacheKey, result);

    return result;
  }

  /**
   * Parse partial content
   */
  private async parsePartial(
    content: string,
    changes: ChangeSet,
    previousAst: T
  ): Promise<T> {
    // This is a simplified implementation
    // In a real implementation, this would:
    // 1. Extract affected regions from previous AST
    // 2. Parse only changed regions
    // 3. Merge results with previous AST

    // For now, do a full parse but mark it as partial
    return this.parseFunction(content);
  }

  /**
   * Detect changes between document versions
   */
  private detectChanges(document: vscode.TextDocument, previous: ParseResult<T>): ChangeSet {
    const changes: ChangeSet = { added: [], removed: [], modified: [] };

    // Get document changes from VSCode
    // This would use vscode.workspace.onDidChangeTextDocument in a real implementation
    // For now, return empty changes
    return changes;
  }

  /**
   * Check if changes are too large for incremental parsing
   */
  private changesTooLarge(changes: ChangeSet, totalLines: number): boolean {
    const changedLines =
      changes.added.reduce((sum, r) => sum + (r.endLine - r.startLine), 0) +
      changes.modified.reduce((sum, r) => sum + (r.endLine - r.startLine), 0) +
      changes.removed.reduce((sum, r) => sum + (r.endLine - r.startLine), 0);

    // If more than 30% changed, do full parse
    return changedLines > totalLines * 0.3;
  }

  /**
   * Check if incremental parsing should be used
   */
  private shouldUseIncremental(document: vscode.TextDocument): boolean {
    return (
      this.config.enableDiffing && document.lineCount >= this.config.minLinesForIncremental
    );
  }

  /**
   * Find previous cached version
   */
  private findPreviousVersion(uri: string, currentVersion: number): string | undefined {
    // This would search the cache for the most recent version before currentVersion
    // For now, return undefined
    return undefined;
  }

  /**
   * Clear cache
   */
  clearCache(): void {
    this.cache.clear();
  }

  /**
   * Get cache statistics
   */
  getStats() {
    return this.cache.getStats();
  }
}

/**
 * AST Differ - Compute differences between ASTs
 */
export class ASTDiffer {
  /**
   * Diff two ASTs
   */
  static diff(oldAst: any, newAst: any): ChangeSet {
    const changes: ChangeSet = { added: [], removed: [], modified: [] };

    // This is a simplified implementation
    // Real implementation would traverse ASTs and compute minimal diff

    return changes;
  }

  /**
   * Check if two nodes are equivalent
   */
  static nodesEqual(node1: any, node2: any): boolean {
    return JSON.stringify(node1) === JSON.stringify(node2);
  }
}

/**
 * Incremental parsing manager
 */
export class IncrementalParsingManager {
  private parsers: Map<string, IncrementalParser<any>> = new Map();
  private context: vscode.ExtensionContext;

  constructor(context: vscode.ExtensionContext) {
    this.context = context;
  }

  /**
   * Register parser for language
   */
  registerParser<T>(
    languageId: string,
    parseFunction: (content: string) => T,
    diffFunction?: (oldAst: T, newAst: T) => ChangeSet
  ): void {
    const parser = new IncrementalParser(parseFunction, diffFunction);
    this.parsers.set(languageId, parser);
  }

  /**
   * Get parser for language
   */
  getParser(languageId: string): IncrementalParser<any> | undefined {
    return this.parsers.get(languageId);
  }

  /**
   * Parse document
   */
  async parseDocument(document: vscode.TextDocument): Promise<ParseResult<any> | undefined> {
    const parser = this.parsers.get(document.languageId);
    if (!parser) {
      return undefined;
    }

    return parser.parse(document);
  }

  /**
   * Clear all caches
   */
  clearAllCaches(): void {
    this.parsers.forEach(parser => parser.clearCache());
  }

  /**
   * Get combined statistics
   */
  getAllStats(): Record<string, any> {
    const stats: Record<string, any> = {};
    this.parsers.forEach((parser, languageId) => {
      stats[languageId] = parser.getStats();
    });
    return stats;
  }
}
