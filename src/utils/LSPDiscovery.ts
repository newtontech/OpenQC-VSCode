/**
 * LSP Discovery Module
 *
 * Provides explicit (opt-in) GitHub-based LSP discovery for refreshing
 * the server list. Normal document-open flows use the bundled registry
 * from `src/lsp/registry.ts` instead.
 *
 * @module LSPDiscovery
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/96
 */

import * as vscode from 'vscode';
import { listBundledLspServers } from '../lsp/registry';
import { createComponentLogger } from './Logger';

export interface LSPServerDefinition {
  /** Repository ID, e.g., "vasp-lsp" */
  id: string;

  /** Human-readable name, e.g., "VASP" */
  name: string;

  /** Full GitHub repository path, e.g., "newtontech/VASP-LSP" */
  repository: string;

  /** Executable name, e.g., "vasp-lsp" */
  executable: string;

  /** VSCode language ID, e.g., "vasp" */
  languageId: string;

  /** File extensions/patterns this LSP handles */
  fileExtensions: string[];

  /** File names (for files without extensions like INCAR, POSCAR) */
  fileNames?: string[];

  /** Whether this LSP is enabled by default */
  enabled: boolean;

  /** Repository URL */
  repositoryUrl: string;

  /** Description from GitHub */
  description?: string;

  /** Last updated timestamp from GitHub */
  lastUpdated?: string;
}

interface GitHubRepo {
  name: string;
  full_name: string;
  description: string | null;
  html_url: string;
  updated_at: string;
  pushed_at: string;
}

interface CacheEntry {
  data: LSPServerDefinition[];
  timestamp: number;
}

export const DEFAULT_LSP_SERVER_DEFINITIONS: readonly LSPServerDefinition[] =
  listBundledLspServers().map(entry => ({
    id: entry.id,
    name: entry.name,
    repository: entry.repository,
    executable: entry.executable,
    languageId: entry.languageId,
    fileExtensions: [...entry.fileExtensions],
    fileNames: [...entry.fileNames],
    enabled: entry.enabled,
    repositoryUrl: entry.repositoryUrl,
  }));

export class LSPDiscovery {
  private static readonly CACHE_KEY = 'openqc.lsp.discovery.cache';
  private static readonly CACHE_TTL_MS = 60 * 60 * 1000; // 1 hour
  private static readonly GITHUB_API_URL = 'https://api.github.com/orgs/newtontech/repos';

  private context: vscode.ExtensionContext | undefined;
  private cache: CacheEntry | null = null;
  private logger = createComponentLogger('LSPDiscovery');

  /**
   * Known LSP to language mappings
   * These are used when auto-detecting from repository names
   */
  private static readonly KNOWN_MAPPINGS = LSPDiscovery.createKnownMappings();

  constructor(context?: vscode.ExtensionContext) {
    this.context = context;
    this.loadCacheFromStorage();
  }

  /**
   * Fetch LSP repositories from GitHub API
   * Uses caching to avoid rate limiting
   */
  async fetchLSPRepositories(forceRefresh = false): Promise<LSPServerDefinition[]> {
    // Check cache first
    if (
      !forceRefresh &&
      this.cache &&
      Date.now() - this.cache.timestamp < LSPDiscovery.CACHE_TTL_MS
    ) {
      this.logger.debug('Using cached LSP list');
      return this.cache.data;
    }

    try {
      console.log('[LSPDiscovery] Fetching LSP repositories from GitHub...');
      const repos = await this.fetchGitHubRepos();
      const lspRepos = repos.filter(repo => repo.name.toLowerCase().includes('lsp'));

      const definitions = lspRepos.map(repo => this.convertToDefinition(repo));

      // Update cache
      this.cache = {
        data: definitions,
        timestamp: Date.now(),
      };
      this.saveCacheToStorage();

      this.logger.info(`Found ${definitions.length} LSP repositories`);
      return definitions;
    } catch (error) {
      this.logger.error('Failed to fetch repositories from GitHub', error as Error);

      // Return cached data if available, even if expired
      if (this.cache) {
        this.logger.warn('Returning stale cache due to error');
        return this.cache.data;
      }

      // Return hardcoded fallback
      this.logger.warn('Returning hardcoded fallback definitions');
      const fallback = this.getFallbackDefinitions();
      this.cache = { data: fallback, timestamp: Date.now() };
      return fallback;
    }
  }

  static getDefaultDefinitions(): LSPServerDefinition[] {
    return DEFAULT_LSP_SERVER_DEFINITIONS.map(definition => ({
      ...definition,
      fileExtensions: [...definition.fileExtensions],
      fileNames: definition.fileNames ? [...definition.fileNames] : undefined,
    }));
  }

  private static createKnownMappings(): Record<
    string,
    {
      name: string;
      languageId: string;
      extensions: string[];
      fileNames?: string[];
      executable: string;
    }
  > {
    const mappings: Record<
      string,
      {
        name: string;
        languageId: string;
        extensions: string[];
        fileNames?: string[];
        executable: string;
      }
    > = {};

    for (const definition of DEFAULT_LSP_SERVER_DEFINITIONS) {
      mappings[definition.id] = {
        name: definition.name,
        languageId: definition.languageId,
        extensions: [...definition.fileExtensions],
        fileNames: definition.fileNames ? [...definition.fileNames] : undefined,
        executable: definition.executable,
      };
    }

    mappings['cp2k-lsp'] = { ...mappings['cp2k-lsp-enhanced'] };
    return mappings;
  }

  /**
   * Fetch repositories from GitHub API
   */
  private async fetchGitHubRepos(): Promise<GitHubRepo[]> {
    const url = `${LSPDiscovery.GITHUB_API_URL}?per_page=100&sort=updated`;

    const response = await fetch(url, {
      headers: {
        Accept: 'application/vnd.github.v3+json',
        'User-Agent': 'OpenQC-VSCode-Extension',
      },
    });

    if (!response.ok) {
      if (response.status === 403) {
        throw new Error('GitHub API rate limit exceeded. Please try again later.');
      }
      throw new Error(`GitHub API error: ${response.status} ${response.statusText}`);
    }

    return response.json();
  }

  /**
   * Convert GitHub repository to LSP server definition
   */
  private convertToDefinition(repo: GitHubRepo): LSPServerDefinition {
    const mapping = LSPDiscovery.KNOWN_MAPPINGS[repo.name] || this.inferMapping(repo.name);

    return {
      id: repo.name,
      name: mapping.name,
      repository: repo.full_name,
      executable: mapping.executable,
      languageId: mapping.languageId,
      fileExtensions: mapping.extensions,
      fileNames: mapping.fileNames,
      enabled: true,
      repositoryUrl: repo.html_url,
      description: repo.description || undefined,
      lastUpdated: repo.pushed_at,
    };
  }

  /**
   * Infer mapping from repository name for unknown LSPs
   */
  private inferMapping(repoName: string): {
    name: string;
    languageId: string;
    extensions: string[];
    fileNames?: string[];
    executable: string;
  } {
    // Remove '-lsp' suffix and convert to title case
    const baseName = repoName.replace(/-lsp(-enhanced)?$/, '');
    const name = baseName
      .split('-')
      .map(word => word.charAt(0).toUpperCase() + word.slice(1))
      .join(' ');

    return {
      name,
      languageId: baseName.toLowerCase(),
      extensions: ['inp'], // Default extension
      fileNames: [],
      executable: repoName,
    };
  }

  /**
   * Get fallback definitions when API fails
   */
  private getFallbackDefinitions(): LSPServerDefinition[] {
    return LSPDiscovery.getDefaultDefinitions();
  }

  /**
   * Load cache from VSCode global storage
   */
  private loadCacheFromStorage(): void {
    if (this.context) {
      const stored = this.context.globalState.get<CacheEntry>(LSPDiscovery.CACHE_KEY);
      if (stored) {
        this.cache = stored;
      }
    }
  }

  /**
   * Save cache to VSCode global storage
   */
  private saveCacheToStorage(): void {
    if (this.context && this.cache) {
      this.context.globalState.update(LSPDiscovery.CACHE_KEY, this.cache);
    }
  }

  /**
   * Clear the cache
   */
  clearCache(): void {
    this.cache = null;
    if (this.context) {
      this.context.globalState.update(LSPDiscovery.CACHE_KEY, undefined);
    }
  }

  /**
   * Get the last cache update time
   */
  getLastCacheTime(): Date | null {
    return this.cache ? new Date(this.cache.timestamp) : null;
  }

  /**
   * Check if new LSPs have been added since last fetch
   */
  async checkForNewLSPs(): Promise<LSPServerDefinition[]> {
    const previousIds = this.cache
      ? this.cache.data.map(l => l.id)
      : this.getFallbackDefinitions().map(l => l.id);
    const currentLSPs = await this.fetchLSPRepositories(true);

    return currentLSPs.filter(lsp => !previousIds.includes(lsp.id));
  }
}

export default LSPDiscovery;
