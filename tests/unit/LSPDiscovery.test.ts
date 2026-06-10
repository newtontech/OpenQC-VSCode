import packageJson from '../../package.json';
import { LSPDiscovery, LSPServerDefinition } from '../../src/utils/LSPDiscovery';
import { Logger } from '../../src/utils/Logger';

type MockContext = {
  globalState: {
    get: jest.Mock;
    update: jest.Mock;
  };
};

const now = 1_800_000_000_000;

const createContext = (cached?: { data: LSPServerDefinition[]; timestamp: number }) =>
  ({
    globalState: {
      get: jest.fn(() => cached),
      update: jest.fn(),
    },
  }) as MockContext;

const mockGitHubResponse = (body: unknown, status = 200, statusText = 'OK') => ({
  ok: status >= 200 && status < 300,
  status,
  statusText,
  json: jest.fn().mockResolvedValue(body),
});

const cp2kRepo = {
  name: 'cp2k-lsp-enhanced',
  full_name: 'newtontech/cp2k-lsp-enhanced',
  description: 'CP2K language server',
  html_url: 'https://github.com/newtontech/cp2k-lsp-enhanced',
  updated_at: '2026-06-01T00:00:00Z',
  pushed_at: '2026-06-02T00:00:00Z',
};

const gaussianRepo = {
  name: 'gaussian-lsp',
  full_name: 'newtontech/gaussian-lsp',
  description: null,
  html_url: 'https://github.com/newtontech/gaussian-lsp',
  updated_at: '2026-06-01T00:00:00Z',
  pushed_at: '2026-06-02T00:00:00Z',
};

const unknownRepo = {
  name: 'my-qc-lsp',
  full_name: 'newtontech/my-qc-lsp',
  description: 'Custom LSP',
  html_url: 'https://github.com/newtontech/my-qc-lsp',
  updated_at: '2026-06-01T00:00:00Z',
  pushed_at: '2026-06-02T00:00:00Z',
};

describe('LSPDiscovery', () => {
  let logSpy: jest.SpyInstance;
  let errorSpy: jest.SpyInstance;
  let dateSpy: jest.SpyInstance;
  let loggerErrorSpy: jest.SpyInstance;
  let loggerWarnSpy: jest.SpyInstance;
  let loggerInfoSpy: jest.SpyInstance;

  beforeEach(() => {
    logSpy = jest.spyOn(console, 'log').mockImplementation(() => undefined);
    errorSpy = jest.spyOn(console, 'error').mockImplementation(() => undefined);
    dateSpy = jest.spyOn(Date, 'now').mockReturnValue(now);
    loggerErrorSpy = jest.spyOn(Logger.getInstance(), 'error').mockImplementation(() => undefined);
    loggerWarnSpy = jest.spyOn(Logger.getInstance(), 'warn').mockImplementation(() => undefined);
    loggerInfoSpy = jest.spyOn(Logger.getInstance(), 'info').mockImplementation(() => undefined);
    (global as any).fetch = jest.fn();
  });

  afterEach(() => {
    logSpy.mockRestore();
    errorSpy.mockRestore();
    dateSpy.mockRestore();
    loggerErrorSpy.mockRestore();
    loggerWarnSpy.mockRestore();
    loggerInfoSpy.mockRestore();
    delete (global as any).fetch;
  });

  it('fetches LSP repositories, maps executables, and stores the cache', async () => {
    const context = createContext();
    (global as any).fetch.mockResolvedValue(
      mockGitHubResponse([
        cp2kRepo,
        gaussianRepo,
        unknownRepo,
        { ...unknownRepo, name: 'not-a-language-server' },
      ])
    );

    const definitions = await new LSPDiscovery(context as any).fetchLSPRepositories();

    expect((global as any).fetch).toHaveBeenCalledWith(
      'https://api.github.com/orgs/newtontech/repos?per_page=100&sort=updated',
      expect.objectContaining({
        headers: expect.objectContaining({
          Accept: 'application/vnd.github.v3+json',
          'User-Agent': 'OpenQC-VSCode-Extension',
        }),
      })
    );
    expect(definitions).toHaveLength(3);
    expect(definitions[0]).toMatchObject({
      id: 'cp2k-lsp-enhanced',
      executable: 'cp2k-language-server',
      languageId: 'cp2k',
      description: 'CP2K language server',
      lastUpdated: '2026-06-02T00:00:00Z',
    });
    expect(definitions[1]).toMatchObject({
      id: 'gaussian-lsp',
      executable: 'gaussian-lsp',
      description: undefined,
    });
    expect(definitions[2]).toMatchObject({
      id: 'my-qc-lsp',
      name: 'My Qc',
      languageId: 'my-qc',
      fileExtensions: ['inp'],
    });
    expect(context.globalState.update).toHaveBeenCalledWith(
      'openqc.lsp.discovery.cache',
      expect.objectContaining({ data: definitions, timestamp: now })
    );
  });

  it('keeps canonical LSP defaults aligned with contributed configuration defaults', () => {
    const properties = packageJson.contributes.configuration.properties as Record<
      string,
      { default: unknown }
    >;

    for (const definition of LSPDiscovery.getDefaultDefinitions()) {
      expect(properties[`openqc.lsp.${definition.languageId}.enabled`]?.default).toBe(
        definition.enabled
      );
      expect(properties[`openqc.lsp.${definition.languageId}.path`]?.default).toBe(
        definition.executable
      );
    }
  });

  it('keeps canonical LSP defaults aligned with contributed language patterns', () => {
    const languages = packageJson.contributes.languages as Array<{
      id: string;
      extensions?: string[];
      filenames?: string[];
    }>;

    for (const definition of LSPDiscovery.getDefaultDefinitions()) {
      const language = languages.find(candidate => candidate.id === definition.languageId);

      expect(language).toBeDefined();
      expect(language?.extensions || []).toEqual(
        definition.fileExtensions.map(extension => `.${extension}`)
      );
      expect(language?.filenames || []).toEqual(definition.fileNames || []);
    }
  });

  it('returns cloned canonical defaults to avoid shared-state drift', () => {
    const first = LSPDiscovery.getDefaultDefinitions();
    const second = LSPDiscovery.getDefaultDefinitions();

    first[0].fileExtensions.push('mutated');

    expect(second[0].fileExtensions).not.toContain('mutated');
  });

  it('uses fresh persisted cache without calling GitHub', async () => {
    const cached = {
      data: [
        {
          id: 'cached-lsp',
          name: 'Cached',
          repository: 'newtontech/cached-lsp',
          executable: 'cached-lsp',
          languageId: 'cached',
          fileExtensions: ['inp'],
          enabled: true,
          repositoryUrl: 'https://github.com/newtontech/cached-lsp',
        },
      ],
      timestamp: now - 1000,
    };

    const discovery = new LSPDiscovery(createContext(cached) as any);
    const definitions = await discovery.fetchLSPRepositories();

    expect(definitions).toBe(cached.data);
    expect((global as any).fetch).not.toHaveBeenCalled();
    expect(discovery.getLastCacheTime()).toEqual(new Date(cached.timestamp));
  });

  it('returns stale cache when a refresh fails', async () => {
    const cached = {
      data: [
        {
          id: 'stale-lsp',
          name: 'Stale',
          repository: 'newtontech/stale-lsp',
          executable: 'stale-lsp',
          languageId: 'stale',
          fileExtensions: ['inp'],
          enabled: true,
          repositoryUrl: 'https://github.com/newtontech/stale-lsp',
        },
      ],
      timestamp: now - 2 * 60 * 60 * 1000,
    };
    (global as any).fetch.mockResolvedValue(mockGitHubResponse([], 500, 'Server Error'));

    const definitions = await new LSPDiscovery(createContext(cached) as any).fetchLSPRepositories();

    expect(definitions).toBe(cached.data);
    expect(loggerErrorSpy).toHaveBeenCalled();
    expect(loggerWarnSpy).toHaveBeenCalled();
  });

  it('uses fallback definitions when GitHub fails and no cache exists', async () => {
    (global as any).fetch.mockResolvedValue(mockGitHubResponse([], 403, 'Forbidden'));

    const definitions = await new LSPDiscovery(createContext() as any).fetchLSPRepositories();
    const defaultDefinitions = LSPDiscovery.getDefaultDefinitions();

    expect(definitions).toHaveLength(defaultDefinitions.length);
    expect(definitions.map(lsp => lsp.id)).toEqual(defaultDefinitions.map(lsp => lsp.id));
    expect(definitions.find(lsp => lsp.id === 'cp2k-lsp-enhanced')).toMatchObject({
      executable: 'cp2k-language-server',
      languageId: 'cp2k',
    });
  });

  it('clears cache state from memory and persisted storage', async () => {
    const cached = {
      data: [
        {
          id: 'cached-lsp',
          name: 'Cached',
          repository: 'newtontech/cached-lsp',
          executable: 'cached-lsp',
          languageId: 'cached',
          fileExtensions: ['inp'],
          enabled: true,
          repositoryUrl: 'https://github.com/newtontech/cached-lsp',
        },
      ],
      timestamp: now,
    };
    const context = createContext(cached);
    const discovery = new LSPDiscovery(context as any);

    discovery.clearCache();

    expect(discovery.getLastCacheTime()).toBeNull();
    expect(context.globalState.update).toHaveBeenCalledWith(
      'openqc.lsp.discovery.cache',
      undefined
    );
  });

  it('detects newly discovered LSP repositories relative to the previous cache', async () => {
    const context = createContext({
      data: [
        {
          id: 'gaussian-lsp',
          name: 'Gaussian',
          repository: 'newtontech/gaussian-lsp',
          executable: 'gaussian-lsp',
          languageId: 'gaussian',
          fileExtensions: ['gjf', 'com'],
          enabled: true,
          repositoryUrl: 'https://github.com/newtontech/gaussian-lsp',
        },
      ],
      timestamp: now - 1000,
    });
    (global as any).fetch.mockResolvedValue(mockGitHubResponse([gaussianRepo, cp2kRepo]));

    const newLsps = await new LSPDiscovery(context as any).checkForNewLSPs();

    expect(newLsps).toHaveLength(1);
    expect(newLsps[0]).toMatchObject({
      id: 'cp2k-lsp-enhanced',
      executable: 'cp2k-language-server',
    });
  });
});
