import { LSPDiscovery, LSPServerDefinition } from '../../../src/utils/LSPDiscovery';

// Mock vscode module
jest.mock('vscode', () => ({
  workspace: {
    getConfiguration: jest.fn(() => ({ get: jest.fn() })),
  },
}));

// Mock global fetch
const mockFetch = jest.fn();
global.fetch = mockFetch;

describe('LSPDiscovery', () => {
  let discovery: LSPDiscovery;

  beforeEach(() => {
    jest.clearAllMocks();
    discovery = new LSPDiscovery();
    discovery.clearCache();
  });

  describe('fetchLSPRepositories', () => {
    it('should return fallback definitions when API fails', async () => {
      mockFetch.mockRejectedValueOnce(new Error('Network error'));
      const definitions = await discovery.fetchLSPRepositories();
      expect(definitions.length).toBeGreaterThan(0);
      expect(definitions.some((d: LSPServerDefinition) => d.id === 'vasp-lsp')).toBe(true);
    });

    it('should return fallback definitions when API returns non-ok', async () => {
      mockFetch.mockResolvedValueOnce({ ok: false, status: 403 } as any);
      const definitions = await discovery.fetchLSPRepositories();
      expect(definitions.length).toBeGreaterThan(0);
    });

    it('should parse GitHub API response correctly', async () => {
      const mockRepos = [
        {
          name: 'test-lsp',
          full_name: 'OpenQuantumChemistry/test-lsp',
          description: 'Test LSP',
          html_url: 'https://github.com/OpenQuantumChemistry/test-lsp',
          pushed_at: '2024-01-01T00:00:00Z',
        },
      ];
      mockFetch.mockResolvedValueOnce({
        ok: true,
        json: async () => mockRepos,
      } as any);
      const definitions = await discovery.fetchLSPRepositories();
      expect(definitions.length).toBeGreaterThan(0);
    });
  });

  describe('caching', () => {
    it('should cache results and use cache on subsequent calls', async () => {
      mockFetch.mockRejectedValueOnce(new Error('Network error'));
      await discovery.fetchLSPRepositories();
      expect(discovery.getLastCacheTime()).toBeInstanceOf(Date);

      // Second call should use cache
      const secondCall = await discovery.fetchLSPRepositories();
      expect(secondCall).toBeDefined();
      expect(mockFetch).toHaveBeenCalledTimes(1);
    });

    it('should allow force refresh', async () => {
      mockFetch.mockRejectedValueOnce(new Error('Network error'));
      await discovery.fetchLSPRepositories();

      mockFetch.mockRejectedValueOnce(new Error('Network error'));
      const refreshed = await discovery.fetchLSPRepositories(true);
      expect(refreshed).toBeDefined();
      expect(mockFetch).toHaveBeenCalledTimes(2);
    });
  });

  describe('clearCache', () => {
    it('should clear the cache', async () => {
      mockFetch.mockRejectedValueOnce(new Error('Network error'));
      await discovery.fetchLSPRepositories();
      expect(discovery.getLastCacheTime()).not.toBeNull();

      discovery.clearCache();
      expect(discovery.getLastCacheTime()).toBeNull();
    });
  });

  describe('VASP LSP fallback', () => {
    it('should have correct VASP configuration in fallback', async () => {
      mockFetch.mockRejectedValueOnce(new Error('Network error'));
      const definitions = await discovery.fetchLSPRepositories();
      const vasp = definitions.find((d: LSPServerDefinition) => d.id === 'vasp-lsp');

      expect(vasp).toBeDefined();
      expect(vasp!.name).toBe('VASP');
      expect(vasp!.languageId).toBe('vasp');
      expect(vasp!.executable).toBe('vasp-lsp');
    });
  });

  describe('Gaussian LSP fallback', () => {
    it('should have correct Gaussian configuration in fallback', async () => {
      mockFetch.mockRejectedValueOnce(new Error('Network error'));
      const definitions = await discovery.fetchLSPRepositories();
      const gaussian = definitions.find((d: LSPServerDefinition) => d.id === 'gaussian-lsp');

      expect(gaussian).toBeDefined();
      expect(gaussian!.name).toBe('Gaussian');
      expect(gaussian!.languageId).toBe('gaussian');
    });
  });
});
