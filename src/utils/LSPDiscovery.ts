/**
 * LSP Discovery Module
 * 
 * Dynamically fetches and manages Language Server Protocol (LSP) repositories
 * from the OpenQuantumChemistry GitHub organization.
 * 
 * This eliminates the need for hardcoded LSP lists in package.json and LSPManager.ts
 * 
 * @module LSPDiscovery
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/13
 */

import * as vscode from 'vscode';

export interface LSPServerDefinition {
  /** Repository ID, e.g., "vasp-lsp" */
  id: string;
  
  /** Human-readable name, e.g., "VASP" */
  name: string;
  
  /** Full GitHub repository path, e.g., "OpenQuantumChemistry/vasp-lsp" */
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

export class LSPDiscovery {
  private static readonly CACHE_KEY = 'openqc.lsp.discovery.cache';
  private static readonly CACHE_TTL_MS = 60 * 60 * 1000; // 1 hour
  private static readonly GITHUB_API_URL = 'https://api.github.com/orgs/OpenQuantumChemistry/repos';
  
  private context: vscode.ExtensionContext | undefined;
  private cache: CacheEntry | null = null;

  /**
   * Known LSP to language mappings
   * These are used when auto-detecting from repository names
   */
  private static readonly KNOWN_MAPPINGS: Record<string, {
    name: string;
    languageId: string;
    extensions: string[];
    fileNames?: string[];
  }> = {
    'vasp-lsp': {
      name: 'VASP',
      languageId: 'vasp',
      extensions: [],
      fileNames: ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'CONTCAR', 'OSZICAR', 'OUTCAR', 'vasprun.xml']
    },
    'gaussian-lsp': {
      name: 'Gaussian',
      languageId: 'gaussian',
      extensions: ['gjf', 'com'],
      fileNames: []
    },
    'orca-lsp': {
      name: 'ORCA',
      languageId: 'orca',
      extensions: ['inp'],
      fileNames: []
    },
    'cp2k-lsp': {
      name: 'CP2K',
      languageId: 'cp2k',
      extensions: ['inp'],
      fileNames: []
    },
    'cp2k-lsp-enhanced': {
      name: 'CP2K',
      languageId: 'cp2k',
      extensions: ['inp'],
      fileNames: []
    },
    'qe-lsp': {
      name: 'Quantum ESPRESSO',
      languageId: 'qe',
      extensions: ['in', 'pw.in', 'relax.in', 'vc-relax.in', 'scf.in', 'nscf.in', 'bands.in', 'ph.in', 'dos.in'],
      fileNames: []
    },
    'gamess-lsp': {
      name: 'GAMESS',
      languageId: 'gamess',
      extensions: ['inp'],
      fileNames: []
    },
    'nwchem-lsp': {
      name: 'NWChem',
      languageId: 'nwchem',
      extensions: ['nw', 'nwinp'],
      fileNames: []
    }
  };

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
    if (!forceRefresh && this.cache && Date.now() - this.cache.timestamp < LSPDiscovery.CACHE_TTL_MS) {
      console.log('[LSPDiscovery] Using cached LSP list');
      return this.cache.data;
    }

    try {
      console.log('[LSPDiscovery] Fetching LSP repositories from GitHub...');
      const repos = await this.fetchGitHubRepos();
      const lspRepos = repos.filter(repo => 
        repo.name.toLowerCase().includes('lsp')
      );
      
      const definitions = lspRepos.map(repo => this.convertToDefinition(repo));
      
      // Update cache
      this.cache = {
        data: definitions,
        timestamp: Date.now()
      };
      this.saveCacheToStorage();
      
      console.log(`[LSPDiscovery] Found ${definitions.length} LSP repositories`);
      return definitions;
    } catch (error) {
      console.error('[LSPDiscovery] Failed to fetch repositories:', error);
      
      // Return cached data if available, even if expired
      if (this.cache) {
        console.log('[LSPDiscovery] Returning stale cache due to error');
        return this.cache.data;
      }
      
      // Return hardcoded fallback
      console.log('[LSPDiscovery] Returning hardcoded fallback');
      return this.getFallbackDefinitions();
    }
  }

  /**
   * Fetch repositories from GitHub API
   */
  private async fetchGitHubRepos(): Promise<GitHubRepo[]> {
    const url = `${LSPDiscovery.GITHUB_API_URL}?per_page=100&sort=updated`;
    
    const response = await fetch(url, {
      headers: {
        'Accept': 'application/vnd.github.v3+json',
        'User-Agent': 'OpenQC-VSCode-Extension'
      }
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
      executable: this.getExecutableName(repo.name),
      languageId: mapping.languageId,
      fileExtensions: mapping.extensions,
      fileNames: mapping.fileNames,
      enabled: true,
      repositoryUrl: repo.html_url,
      description: repo.description || undefined,
      lastUpdated: repo.pushed_at
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
      fileNames: []
    };
  }

  /**
   * Get executable name from repository name
   */
  private getExecutableName(repoName: string): string {
    // Use repo name as executable (e.g., "vasp-lsp" -> "vasp-lsp")
    // For cp2k-lsp-enhanced, use "cp2k-lsp-enhanced"
    return repoName;
  }

  /**
   * Get fallback definitions when API fails
   */
  private getFallbackDefinitions(): LSPServerDefinition[] {
    return [
      {
        id: 'vasp-lsp',
        name: 'VASP',
        repository: 'OpenQuantumChemistry/vasp-lsp',
        executable: 'vasp-lsp',
        languageId: 'vasp',
        fileExtensions: [],
        fileNames: ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR'],
        enabled: true,
        repositoryUrl: 'https://github.com/OpenQuantumChemistry/vasp-lsp'
      },
      {
        id: 'gaussian-lsp',
        name: 'Gaussian',
        repository: 'OpenQuantumChemistry/gaussian-lsp',
        executable: 'gaussian-lsp',
        languageId: 'gaussian',
        fileExtensions: ['gjf', 'com'],
        enabled: true,
        repositoryUrl: 'https://github.com/OpenQuantumChemistry/gaussian-lsp'
      },
      {
        id: 'orca-lsp',
        name: 'ORCA',
        repository: 'OpenQuantumChemistry/orca-lsp',
        executable: 'orca-lsp',
        languageId: 'orca',
        fileExtensions: ['inp'],
        enabled: true,
        repositoryUrl: 'https://github.com/OpenQuantumChemistry/orca-lsp'
      },
      {
        id: 'cp2k-lsp-enhanced',
        name: 'CP2K',
        repository: 'OpenQuantumChemistry/cp2k-lsp-enhanced',
        executable: 'cp2k-lsp-enhanced',
        languageId: 'cp2k',
        fileExtensions: ['inp'],
        enabled: true,
        repositoryUrl: 'https://github.com/OpenQuantumChemistry/cp2k-lsp-enhanced'
      },
      {
        id: 'qe-lsp',
        name: 'Quantum ESPRESSO',
        repository: 'OpenQuantumChemistry/qe-lsp',
        executable: 'qe-lsp',
        languageId: 'qe',
        fileExtensions: ['in', 'pw.in', 'relax.in'],
        enabled: true,
        repositoryUrl: 'https://github.com/OpenQuantumChemistry/qe-lsp'
      },
      {
        id: 'gamess-lsp',
        name: 'GAMESS',
        repository: 'OpenQuantumChemistry/gamess-lsp',
        executable: 'gamess-lsp',
        languageId: 'gamess',
        fileExtensions: ['inp'],
        enabled: true,
        repositoryUrl: 'https://github.com/OpenQuantumChemistry/gamess-lsp'
      },
      {
        id: 'nwchem-lsp',
        name: 'NWChem',
        repository: 'OpenQuantumChemistry/nwchem-lsp',
        executable: 'nwchem-lsp',
        languageId: 'nwchem',
        fileExtensions: ['nw', 'nwinp'],
        enabled: true,
        repositoryUrl: 'https://github.com/OpenQuantumChemistry/nwchem-lsp'
      }
    ];
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
    const currentLSPs = await this.fetchLSPRepositories();
    const previousIds = this.cache ? this.cache.data.map(l => l.id) 
      : this.getFallbackDefinitions().map(l => l.id);
    
    return currentLSPs.filter(lsp => !previousIds.includes(lsp.id));
  }
}

export default LSPDiscovery;
