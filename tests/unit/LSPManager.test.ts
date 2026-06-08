import { LSPManager } from '../../src/managers/LSPManager';
import { LSPDiscovery } from '../../src/utils/LSPDiscovery';
import { Logger } from '../../src/utils/Logger';

// Mock functions object - declared with var for hoisting compatibility
const mockFunctions: {
  showInfo: jest.Mock;
  showWarning: jest.Mock;
  showError: jest.Mock;
  languageClient: jest.Mock;
  statusBarItem: {
    show: jest.Mock;
    hide: jest.Mock;
    dispose: jest.Mock;
    text?: string;
    tooltip?: string;
    command?: string;
    name?: string;
  };
} = {
  showInfo: jest.fn(),
  showWarning: jest.fn(),
  showError: jest.fn(),
  languageClient: jest.fn(),
  statusBarItem: {
    show: jest.fn(),
    hide: jest.fn(),
    dispose: jest.fn(),
  },
};

// Mock vscode module
jest.mock('vscode', () => ({
  window: {
    showInformationMessage: (...args: any[]) => mockFunctions.showInfo(...args),
    showWarningMessage: (...args: any[]) => mockFunctions.showWarning(...args),
    showErrorMessage: (...args: any[]) => mockFunctions.showError(...args),
    createOutputChannel: jest.fn(() => ({
      appendLine: jest.fn(),
      append: jest.fn(),
      clear: jest.fn(),
      show: jest.fn(),
      hide: jest.fn(),
      dispose: jest.fn(),
    })),
    createStatusBarItem: jest.fn(() => mockFunctions.statusBarItem),
  },
  StatusBarAlignment: { Left: 1, Right: 2 },
  ConfigurationTarget: { Workspace: 1 },
  workspace: {
    getConfiguration: jest.fn(() => ({
      get: jest.fn((key: string, defaultValue?: any) => {
        if (key.includes('enabled')) return true;
        if (key.includes('path')) {
          const softwareMap: Record<string, string> = {
            cp2k: 'cp2k-language-server',
            gaussian: 'gaussian-lsp',
            vasp: 'vasp-lsp',
            orca: 'orca-lsp',
            qe: 'qe-lsp',
            gamess: 'gamess-lsp',
            nwchem: 'nwchem-lsp',
          };
          for (const [sw, path] of Object.entries(softwareMap)) {
            if (key.includes(sw)) return path;
          }
        }
        return defaultValue;
      }),
    })),
    getWorkspaceFolder: jest.fn(),
    createFileSystemWatcher: jest.fn(() => ({
      onDidCreate: jest.fn(),
      onDidChange: jest.fn(),
      onDidDelete: jest.fn(),
      dispose: jest.fn(),
    })),
  },
  RelativePattern: jest.fn().mockImplementation((base: any, pattern: string) => ({
    base,
    pattern,
  })),
  env: {
    remoteName: undefined,
  },
}));

// Mock the command resolver's isExecutableAvailable to always succeed in tests
jest.mock('../../src/lsp/commandResolver', () => {
  const actual = jest.requireActual('../../src/lsp/commandResolver');
  return {
    ...actual,
    isExecutableAvailable: jest.fn().mockResolvedValue(true),
  };
});

// Mock vscode-languageclient/node
jest.mock('vscode-languageclient/node', () => ({
  LanguageClient: jest
    .fn()
    .mockImplementation((...args: any[]) => mockFunctions.languageClient(...args)),
  TransportKind: { stdio: 0 },
}));

describe('LSPManager', () => {
  let lspManager: LSPManager;
  let mockDocument: any;
  let mockDiscovery: Pick<LSPDiscovery, 'fetchLSPRepositories'>;

  beforeEach(() => {
    jest.clearAllMocks();
    mockFunctions.statusBarItem.text = undefined;
    mockFunctions.statusBarItem.tooltip = undefined;
    mockFunctions.statusBarItem.command = undefined;
    mockFunctions.statusBarItem.name = undefined;
    Logger.resetInstance();
    // Reset the isExecutableAvailable mock to default (true)
    const { isExecutableAvailable } = require('../../src/lsp/commandResolver');
    (isExecutableAvailable as jest.Mock).mockResolvedValue(true);
    const vscode = require('vscode');
    vscode.env.remoteName = undefined;
    vscode.workspace.getWorkspaceFolder.mockReturnValue(undefined);
    mockDiscovery = {
      fetchLSPRepositories: jest.fn().mockResolvedValue(LSPDiscovery.getDefaultDefinitions()),
    };
    lspManager = new LSPManager(undefined, mockDiscovery as LSPDiscovery);
    mockDocument = {
      uri: { fsPath: '/test/file.com' },
      fileName: '/test/file.com',
      languageId: 'gaussian',
      getText: jest.fn().mockReturnValue('%chk=test.chk\n# B3LYP/6-31G(d)\n\n0 1\nH 0 0 0'),
    } as any;
    // Default mock implementations
    mockFunctions.languageClient.mockImplementation(() => ({
      start: jest.fn().mockResolvedValue(undefined),
      stop: jest.fn().mockResolvedValue(undefined),
      needsStop: jest.fn().mockReturnValue(true),
    }));
  });

  afterEach(async () => {
    await lspManager.dispose();
    Logger.resetInstance();
  });

  describe('startLSPForDocument', () => {
    it('should start LSP for a valid document using bundled registry (no GitHub call)', async () => {
      await lspManager.startLSPForDocument(mockDocument);
      // Bundled registry is used instead of GitHub discovery during normal flows
      expect(mockDiscovery.fetchLSPRepositories).not.toHaveBeenCalled();
      expect(mockFunctions.showInfo).not.toHaveBeenCalled();
      expect(mockFunctions.statusBarItem.text).toBe('$(check) OpenQC LSP: Gaussian running');
    });

    it('should not show automatic success popups after LSP startup', async () => {
      await lspManager.startLSPForDocument(mockDocument);

      expect(mockFunctions.showInfo).not.toHaveBeenCalledWith(
        expect.stringContaining('Language Server started')
      );
      expect(mockFunctions.showError).not.toHaveBeenCalled();
    });

    it('should use cross-platform command resolution for executable checks', async () => {
      const { isExecutableAvailable } = require('../../src/lsp/commandResolver');
      await lspManager.startLSPForDocument(mockDocument);

      // isExecutableAvailable should be called (from the mocked module)
      expect(isExecutableAvailable).toHaveBeenCalled();
      expect(mockFunctions.languageClient).toHaveBeenCalledWith(
        expect.stringContaining('openqc-gaussian-'),
        'OpenQC Gaussian Language Server',
        expect.objectContaining({ command: 'gaussian-lsp' }),
        expect.objectContaining({
          documentSelector: [{ scheme: 'file', language: 'gaussian' }],
        })
      );
    });

    it('should not start LSP if software cannot be detected', async () => {
      const unknownDoc = {
        ...mockDocument,
        fileName: '/test/unknown.xyz',
        getText: jest.fn().mockReturnValue('unknown'),
      };
      await lspManager.startLSPForDocument(unknownDoc);
      expect(mockFunctions.showInfo).not.toHaveBeenCalled();
    });

    it('should not start LSP if already running', async () => {
      await lspManager.startLSPForDocument(mockDocument);
      await lspManager.startLSPForDocument(mockDocument);
      expect(mockFunctions.languageClient).toHaveBeenCalledTimes(1);
    });

    it('should start separate clients for the same language in different workspace folders', async () => {
      const vscode = require('vscode');
      const workspaceA = {
        name: 'workspace-a',
        index: 0,
        uri: { fsPath: '/workspace-a', scheme: 'file', toString: () => 'file:///workspace-a' },
      };
      const workspaceB = {
        name: 'workspace-b',
        index: 1,
        uri: { fsPath: '/workspace-b', scheme: 'file', toString: () => 'file:///workspace-b' },
      };
      const docA = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace-a/file.com',
          scheme: 'file',
          toString: () => 'file:///workspace-a/file.com',
        },
        fileName: '/workspace-a/file.com',
      };
      const docB = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace-b/file.com',
          scheme: 'file',
          toString: () => 'file:///workspace-b/file.com',
        },
        fileName: '/workspace-b/file.com',
      };
      vscode.workspace.getWorkspaceFolder.mockImplementation((uri: any) =>
        uri.fsPath.startsWith('/workspace-a') ? workspaceA : workspaceB
      );

      await lspManager.startLSPForDocument(docA);
      await lspManager.startLSPForDocument(docB);

      expect(mockFunctions.languageClient).toHaveBeenCalledTimes(2);
      expect(mockFunctions.languageClient).toHaveBeenNthCalledWith(
        1,
        expect.stringContaining('workspace-a'),
        'OpenQC Gaussian Language Server (workspace-a)',
        expect.any(Object),
        expect.objectContaining({ workspaceFolder: workspaceA })
      );
      expect(mockFunctions.languageClient).toHaveBeenNthCalledWith(
        2,
        expect.stringContaining('workspace-b'),
        'OpenQC Gaussian Language Server (workspace-b)',
        expect.any(Object),
        expect.objectContaining({ workspaceFolder: workspaceB })
      );
    });

    it('should coalesce concurrent starts for the same language and workspace', async () => {
      const vscode = require('vscode');
      const workspaceFolder = {
        name: 'chemistry',
        index: 0,
        uri: { fsPath: '/workspace', scheme: 'file', toString: () => 'file:///workspace' },
      };
      vscode.workspace.getWorkspaceFolder.mockReturnValue(workspaceFolder);
      const mockStop = jest.fn().mockResolvedValue(undefined);
      let resolveStart: () => void = () => {};
      const startPromise = new Promise<void>(resolve => {
        resolveStart = resolve;
      });
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockReturnValue(startPromise),
        stop: mockStop,
        needsStop: jest.fn().mockReturnValue(true),
      }));
      const docA = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace/a.com',
          scheme: 'file',
          toString: () => 'file:///workspace/a.com',
        },
        fileName: '/workspace/a.com',
      };
      const docB = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace/b.com',
          scheme: 'file',
          toString: () => 'file:///workspace/b.com',
        },
        fileName: '/workspace/b.com',
      };

      const firstStart = lspManager.startLSPForDocument(docA);
      const secondStart = lspManager.startLSPForDocument(docB);

      await new Promise(resolve => setImmediate(resolve));

      expect(mockFunctions.languageClient).toHaveBeenCalledTimes(1);

      resolveStart();
      await Promise.all([firstStart, secondStart]);
      await lspManager.stopLSPForDocument(docA);

      expect(mockStop).not.toHaveBeenCalled();

      await lspManager.stopLSPForDocument(docB);

      expect(mockStop).toHaveBeenCalledTimes(1);
    });

    it('should handle errors when LSP executable is not found', async () => {
      const { isExecutableAvailable } = require('../../src/lsp/commandResolver');
      (isExecutableAvailable as jest.Mock).mockResolvedValueOnce(false);
      await lspManager.startLSPForDocument(mockDocument);
      expect(mockFunctions.showError).toHaveBeenCalledWith(
        expect.stringContaining('Failed to start')
      );
    });

    it('should include remote host context in executable startup errors', async () => {
      const vscode = require('vscode');
      const { isExecutableAvailable } = require('../../src/lsp/commandResolver');
      vscode.env.remoteName = 'ssh-remote';
      (isExecutableAvailable as jest.Mock).mockResolvedValueOnce(false);

      await lspManager.startLSPForDocument(mockDocument);

      expect(mockFunctions.showError).toHaveBeenCalledWith(
        expect.stringContaining('remote extension host (Remote SSH)')
      );
      expect(mockFunctions.showError).toHaveBeenCalledWith(
        expect.stringContaining('openqc.lsp.gaussian.command')
      );
      expect(mockFunctions.showError).toHaveBeenCalledWith(
        expect.stringContaining('remote workspace environment')
      );
    });

    it('should preserve user executable path overrides over discovery defaults', async () => {
      const vscode = require('vscode');
      vscode.workspace.getConfiguration.mockReturnValueOnce({
        get: jest.fn((key: string, defaultValue?: any) => {
          if (key === 'gaussian.enabled') return true;
          if (key === 'gaussian.path') return '/opt/openqc/custom-gaussian-lsp';
          return defaultValue;
        }),
      });
      lspManager = new LSPManager(undefined, mockDiscovery as LSPDiscovery);

      await lspManager.startLSPForDocument(mockDocument);

      // The resolved command should be the custom path, not the default executable.
      // The path contains '/' so it is treated as a filesystem path (fs.access check),
      // but in test environment fs.access will fail since the file doesn't exist.
      // We need to make execFile succeed OR check the LanguageClient was constructed
      // with the overridden command.
      expect(mockFunctions.languageClient).toHaveBeenCalledWith(
        expect.stringContaining('openqc-gaussian-'),
        expect.any(String),
        expect.objectContaining({ command: '/opt/openqc/custom-gaussian-lsp' }),
        expect.any(Object)
      );
    });

    it('should handle paths with spaces in user config', async () => {
      const vscode = require('vscode');
      vscode.workspace.getConfiguration.mockReturnValueOnce({
        get: jest.fn((key: string, defaultValue?: any) => {
          if (key === 'gaussian.enabled') return true;
          if (key === 'gaussian.path') return '/opt/my tools/gaussian-lsp';
          return defaultValue;
        }),
      });
      lspManager = new LSPManager(undefined, mockDiscovery as LSPDiscovery);

      await lspManager.startLSPForDocument(mockDocument);

      expect(mockFunctions.languageClient).toHaveBeenCalledWith(
        expect.any(String),
        expect.any(String),
        expect.objectContaining({ command: '/opt/my tools/gaussian-lsp' }),
        expect.any(Object)
      );
    });

    it('should clean up client map on error', async () => {
      mockFunctions.languageClient
        .mockImplementationOnce(() => ({
          start: jest.fn().mockRejectedValue(new Error('Start failed')),
          stop: jest.fn(),
          needsStop: jest.fn().mockReturnValue(false),
        }))
        .mockImplementation(() => ({
          start: jest.fn().mockResolvedValue(undefined),
          stop: jest.fn(),
          needsStop: jest.fn().mockReturnValue(true),
        }));
      await lspManager.startLSPForDocument(mockDocument);
      await lspManager.startLSPForDocument(mockDocument);
      expect(mockFunctions.languageClient).toHaveBeenCalledTimes(2);
    });
  });

  describe('stopLSPForDocument', () => {
    it('should stop running LSP', async () => {
      const mockStop = jest.fn().mockResolvedValue(undefined);
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: mockStop,
        needsStop: jest.fn().mockReturnValue(true),
      }));
      await lspManager.startLSPForDocument(mockDocument);
      await lspManager.stopLSPForDocument(mockDocument);
      expect(mockStop).toHaveBeenCalled();
    });

    it('should handle errors when stopping LSP', async () => {
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: jest.fn().mockRejectedValue(new Error('Stop failed')),
        needsStop: jest.fn().mockReturnValue(true),
      }));
      await lspManager.startLSPForDocument(mockDocument);
      await lspManager.stopLSPForDocument(mockDocument);
      expect(mockFunctions.showWarning).toHaveBeenCalledWith(
        expect.stringContaining('Error stopping')
      );
    });

    it('should not call stop if client does not need stop', async () => {
      const mockStop = jest.fn();
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: mockStop,
        needsStop: jest.fn().mockReturnValue(false),
      }));
      await lspManager.startLSPForDocument(mockDocument);
      await lspManager.stopLSPForDocument(mockDocument);
      expect(mockStop).not.toHaveBeenCalled();
    });

    it('should do nothing if software is not detected', async () => {
      const unknownDoc = {
        ...mockDocument,
        fileName: '/test/unknown.xyz',
        getText: jest.fn().mockReturnValue('unknown'),
      };
      await expect(lspManager.stopLSPForDocument(unknownDoc)).resolves.not.toThrow();
    });

    it('should keep a workspace client running until the last matching document closes', async () => {
      const vscode = require('vscode');
      const workspaceFolder = {
        name: 'chemistry',
        index: 0,
        uri: { fsPath: '/workspace', scheme: 'file', toString: () => 'file:///workspace' },
      };
      vscode.workspace.getWorkspaceFolder.mockReturnValue(workspaceFolder);
      const mockStop = jest.fn().mockResolvedValue(undefined);
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: mockStop,
        needsStop: jest.fn().mockReturnValue(true),
      }));
      const docA = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace/a.com',
          scheme: 'file',
          toString: () => 'file:///workspace/a.com',
        },
        fileName: '/workspace/a.com',
      };
      const docB = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace/b.com',
          scheme: 'file',
          toString: () => 'file:///workspace/b.com',
        },
        fileName: '/workspace/b.com',
      };

      await lspManager.startLSPForDocument(docA);
      await lspManager.startLSPForDocument(docB);
      await lspManager.stopLSPForDocument(docA);

      expect(mockStop).not.toHaveBeenCalled();

      await lspManager.stopLSPForDocument(docB);

      expect(mockStop).toHaveBeenCalledTimes(1);
    });

    it('should stop only the client for the closed document workspace', async () => {
      const vscode = require('vscode');
      const workspaceA = {
        name: 'workspace-a',
        index: 0,
        uri: { fsPath: '/workspace-a', scheme: 'file', toString: () => 'file:///workspace-a' },
      };
      const workspaceB = {
        name: 'workspace-b',
        index: 1,
        uri: { fsPath: '/workspace-b', scheme: 'file', toString: () => 'file:///workspace-b' },
      };
      const stopA = jest.fn().mockResolvedValue(undefined);
      const stopB = jest.fn().mockResolvedValue(undefined);
      mockFunctions.languageClient
        .mockImplementationOnce(() => ({
          start: jest.fn().mockResolvedValue(undefined),
          stop: stopA,
          needsStop: jest.fn().mockReturnValue(true),
        }))
        .mockImplementationOnce(() => ({
          start: jest.fn().mockResolvedValue(undefined),
          stop: stopB,
          needsStop: jest.fn().mockReturnValue(true),
        }));
      const docA = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace-a/file.com',
          scheme: 'file',
          toString: () => 'file:///workspace-a/file.com',
        },
        fileName: '/workspace-a/file.com',
      };
      const docB = {
        ...mockDocument,
        uri: {
          fsPath: '/workspace-b/file.com',
          scheme: 'file',
          toString: () => 'file:///workspace-b/file.com',
        },
        fileName: '/workspace-b/file.com',
      };
      vscode.workspace.getWorkspaceFolder.mockImplementation((uri: any) =>
        uri.fsPath.startsWith('/workspace-a') ? workspaceA : workspaceB
      );

      await lspManager.startLSPForDocument(docA);
      await lspManager.startLSPForDocument(docB);
      await lspManager.stopLSPForDocument(docA);

      expect(stopA).toHaveBeenCalledTimes(1);
      expect(stopB).not.toHaveBeenCalled();
    });
  });

  describe('restartLSPForDocument', () => {
    it('should successfully restart LSP', async () => {
      const mockStart = jest.fn().mockResolvedValue(undefined);
      const mockStop = jest.fn().mockResolvedValue(undefined);
      mockFunctions.languageClient.mockImplementation(() => ({
        start: mockStart,
        stop: mockStop,
        needsStop: jest.fn().mockReturnValue(true),
      }));
      // First start the LSP, then restart it
      await lspManager.startLSPForDocument(mockDocument);
      await lspManager.restartLSPForDocument(mockDocument);
      expect(mockStop).toHaveBeenCalled();
      expect(mockStart).toHaveBeenCalledTimes(2); // Once for initial start, once for restart
      expect(mockFunctions.showInfo).not.toHaveBeenCalledWith(
        expect.stringContaining('Language Server started')
      );
    });

    it('should show warning if software cannot be detected', async () => {
      const unknownDoc = {
        ...mockDocument,
        fileName: '/test/unknown.xyz',
        getText: jest.fn().mockReturnValue('unknown'),
      };
      await lspManager.restartLSPForDocument(unknownDoc);
      expect(mockFunctions.showWarning).toHaveBeenCalledWith(
        'Could not detect quantum chemistry software for this file'
      );
    });

    it('should show error if restart fails', async () => {
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockRejectedValue(new Error('Restart failed')),
        stop: jest.fn(),
        needsStop: jest.fn().mockReturnValue(true),
      }));
      await lspManager.restartLSPForDocument(mockDocument);
      // When restart fails during start, it shows 'Failed to start'
      expect(mockFunctions.showError).toHaveBeenCalledWith(
        expect.stringContaining('Failed to start')
      );
    });

    it('should not use a fixed delay between stop and start', async () => {
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: jest.fn().mockResolvedValue(undefined),
        needsStop: jest.fn().mockReturnValue(true),
      }));
      const startTime = Date.now();
      await lspManager.restartLSPForDocument(mockDocument);
      const elapsed = Date.now() - startTime;
      // Should complete well under 500ms since there is no longer a fixed sleep.
      // Allow generous overhead for test execution but assert no 500ms delay.
      expect(elapsed).toBeLessThan(500);
    });
  });

  describe('dispose', () => {
    it('should stop all clients on dispose', async () => {
      const mockStop = jest.fn().mockResolvedValue(undefined);
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: mockStop,
        needsStop: jest.fn().mockReturnValue(true),
      }));
      const gaussianDoc = {
        ...mockDocument,
        fileName: '/test/file.com',
        getText: jest.fn().mockReturnValue('%chk=test\n# B3LYP'),
      };
      const vaspDoc = {
        ...mockDocument,
        fileName: '/test/INCAR',
        getText: jest.fn().mockReturnValue('ENCUT=520'),
      };
      await lspManager.startLSPForDocument(gaussianDoc);
      await lspManager.startLSPForDocument(vaspDoc);
      await lspManager.dispose();
      expect(mockStop).toHaveBeenCalledTimes(2);
    });

    it('should handle errors during dispose', async () => {
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: jest.fn().mockRejectedValue(new Error('Dispose failed')),
        needsStop: jest.fn().mockReturnValue(true),
      }));
      await lspManager.startLSPForDocument(mockDocument);
      await lspManager.dispose();
      // Logger error is called via output channel, no console.error needed
      // The key behavior is that dispose completes without throwing
    });
  });

  describe('getLanguageId and getExtensions', () => {
    it('should map all software types correctly', async () => {
      // Use unique file extensions/names for each software type to ensure
      // each one creates a new LSP client
      // Content must match enough patterns for FileTypeDetector (confidence > 0.5)
      const testCases = [
        { fileName: '/test/cp2k.inp', content: '&GLOBAL\n&FORCE_EVAL\nPROJECT_NAME test' },
        { fileName: '/test/INCAR', content: 'ENCUT=520' },
        { fileName: '/test/gaussian.com', content: '%chk=test\n# B3LYP\n\n0 1' },
        { fileName: '/test/orca.inp', content: '! HF\n%pal nprocs 4 end\n%maxcore 2000' },
        {
          fileName: '/test/qe.in',
          content: '&CONTROL\n&SYSTEM\ncalculation = "scf"\npseudo_dir = "./"',
        },
        { fileName: '/test/gamess.inp', content: '$BASIS\n$CONTRL\n$SYSTEM\nruntyp=energy' },
        {
          fileName: '/test/nwchem.nw',
          content: 'geometry\n  H 0 0 0\nend\n\nbasis\n  * library 6-31G*\nend',
        },
      ];
      for (const tc of testCases) {
        const doc: any = {
          uri: { fsPath: tc.fileName },
          fileName: tc.fileName,
          getText: () => tc.content,
        };
        await lspManager.startLSPForDocument(doc);
      }
      expect(mockFunctions.languageClient).toHaveBeenCalledTimes(7);
    });
  });
});
