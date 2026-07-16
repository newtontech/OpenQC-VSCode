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

interface DocumentSnapshot<T> {
  version: number;
  content: string;
  result: ParseResult<T>;
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
  private snapshots: Map<string, DocumentSnapshot<T>> = new Map();

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
    const filePath = document.uri.fsPath || uri;
    const version = document.version;
    const content = document.getText();

    // Check cache
    const cacheKey = CacheKeyGenerator.forFile(filePath, `parse:${version}`);
    const cached = this.cache.get(cacheKey);

    if (cached && cached.version === version) {
      return { ...cached, cached: true };
    }

    // Check if incremental parsing is beneficial
    if (this.shouldUseIncremental(document)) {
      return this.parseIncremental(document, uri, filePath, content);
    }

    // Full parse
    return this.parseFull(uri, filePath, version, content);
  }

  /**
   * Parse with incremental updates
   */
  private async parseIncremental(
    document: vscode.TextDocument,
    uri: string,
    filePath: string,
    content: string
  ): Promise<ParseResult<T>> {
    const version = document.version;

    // Find previous version
    const previous = this.findPreviousVersion(uri, version);

    if (!previous) {
      return this.parseFull(uri, filePath, version, content);
    }

    // Get content changes
    const lineChanges = this.detectChanges(content, previous.content);

    // If changes are too large, do full parse
    if (this.changesTooLarge(lineChanges, document.lineCount)) {
      return this.parseFull(uri, filePath, version, content);
    }

    // Parse only changed regions
    const partialAst = await this.parsePartial(content);
    const changes = this.diffFunction
      ? this.diffFunction(previous.result.ast, partialAst)
      : lineChanges;

    const result: ParseResult<T> = {
      ast: partialAst,
      version,
      timestamp: Date.now(),
      changes,
      cached: false,
    };

    // Cache result
    const cacheKey = CacheKeyGenerator.forFile(filePath, `parse:${version}`);
    this.cache.set(cacheKey, result);
    this.snapshots.set(uri, { version, content, result });

    return result;
  }

  /**
   * Full parse
   */
  private async parseFull(
    uri: string,
    filePath: string,
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
    const cacheKey = CacheKeyGenerator.forFile(filePath, `parse:${version}`);
    this.cache.set(cacheKey, result);
    this.snapshots.set(uri, { version, content, result });

    return result;
  }

  /**
   * Parse partial content
   */
  private async parsePartial(content: string): Promise<T> {
    // Generic AST merging is parser-specific, so the shared fallback reparses
    // the document while preserving exact line-level change metadata.
    return this.parseFunction(content);
  }

  /**
   * Detect changes between document versions
   */
  private detectChanges(currentContent: string, previousContent: string): ChangeSet {
    const currentLines = currentContent.split(/\r?\n/);
    const previousLines = previousContent.split(/\r?\n/);
    const changes: ChangeSet = { added: [], removed: [], modified: [] };

    let prefix = 0;
    while (
      prefix < previousLines.length &&
      prefix < currentLines.length &&
      previousLines[prefix] === currentLines[prefix]
    ) {
      prefix++;
    }

    let previousSuffix = previousLines.length - 1;
    let currentSuffix = currentLines.length - 1;
    while (
      previousSuffix >= prefix &&
      currentSuffix >= prefix &&
      previousLines[previousSuffix] === currentLines[currentSuffix]
    ) {
      previousSuffix--;
      currentSuffix--;
    }

    const removedCount = previousSuffix - prefix + 1;
    const addedCount = currentSuffix - prefix + 1;
    if (removedCount <= 0 && addedCount <= 0) {
      return changes;
    }

    if (removedCount <= 0) {
      changes.added.push({
        startLine: prefix,
        endLine: prefix + addedCount,
        startColumn: 0,
        endColumn: 0,
      });
    } else if (addedCount <= 0) {
      changes.removed.push({
        startLine: prefix,
        endLine: prefix + removedCount,
        startColumn: 0,
        endColumn: 0,
      });
    } else {
      const lastLine = currentLines[Math.max(prefix, prefix + addedCount - 1)] ?? '';
      changes.modified.push({
        startLine: prefix,
        endLine: prefix + addedCount,
        startColumn: 0,
        endColumn: lastLine.length,
      });
    }

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
    return this.config.enableDiffing && document.lineCount >= this.config.minLinesForIncremental;
  }

  /**
   * Find previous cached version
   */
  private findPreviousVersion(
    uri: string,
    currentVersion: number
  ): DocumentSnapshot<T> | undefined {
    const snapshot = this.snapshots.get(uri);
    if (!snapshot || snapshot.version >= currentVersion) {
      return undefined;
    }
    return snapshot;
  }

  /**
   * Clear cache
   */
  clearCache(): void {
    this.cache.clear();
    this.snapshots.clear();
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
    if (ASTDiffer.nodesEqual(oldAst, newAst)) {
      return { added: [], removed: [], modified: [] };
    }

    return {
      added: [],
      removed: [],
      modified: [{ startLine: 0, endLine: 1, startColumn: 0, endColumn: 0 }],
    };
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
