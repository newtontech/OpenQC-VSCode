import { LSPManager } from '../../src/managers/LSPManager';

// Mock functions object - declared with var for hoisting compatibility
const mockFunctions: {
  showInfo: jest.Mock;
  showWarning: jest.Mock;
  showError: jest.Mock;
  exec: jest.Mock;
  languageClient: jest.Mock;
} = {
  showInfo: jest.fn(),
  showWarning: jest.fn(),
  showError: jest.fn(),
  exec: jest.fn(),
  languageClient: jest.fn(),
};

// Mock vscode module
jest.mock('vscode', () => ({
  window: {
    showInformationMessage: (...args: any[]) => mockFunctions.showInfo(...args),
    showWarningMessage: (...args: any[]) => mockFunctions.showWarning(...args),
    showErrorMessage: (...args: any[]) => mockFunctions.showError(...args),
  },
  workspace: {
    getConfiguration: jest.fn(() => ({
      get: jest.fn((key: string, defaultValue?: any) => {
        if (key.includes('enabled')) return true;
        if (key.includes('path')) {
          const softwareMap: Record<string, string> = {
            cp2k: 'cp2k-lsp-enhanced',
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
    createFileSystemWatcher: jest.fn(() => ({
      onDidCreate: jest.fn(),
      onDidChange: jest.fn(),
      onDidDelete: jest.fn(),
      dispose: jest.fn(),
    })),
  },
}));

// Mock child_process
jest.mock('child_process', () => ({
  exec: (...args: any[]) => mockFunctions.exec(...args),
}));

// Mock util - promisify should return a function that returns a Promise
jest.mock('util', () => ({
  promisify: jest.fn((fn: Function) => {
    return (...args: any[]) => {
      return new Promise((resolve, reject) => {
        fn(...args, (err: any, result: any) => {
          if (err) reject(err);
          else resolve(result);
        });
      });
    };
  }),
}));

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

  beforeEach(() => {
    jest.clearAllMocks();
    lspManager = new LSPManager();
    mockDocument = {
      uri: { fsPath: '/test/file.com' },
      fileName: '/test/file.com',
      languageId: 'gaussian',
      getText: jest.fn().mockReturnValue('%chk=test.chk\n# B3LYP/6-31G(d)\n\n0 1\nH 0 0 0'),
    } as any;
    // Default mock implementations
    mockFunctions.exec.mockImplementation((cmd: string, callback: Function) => {
      callback(null, { stdout: '/usr/bin/test' });
    });
    mockFunctions.languageClient.mockImplementation(() => ({
      start: jest.fn().mockResolvedValue(undefined),
      stop: jest.fn().mockResolvedValue(undefined),
      needsStop: jest.fn().mockReturnValue(true),
    }));
  });

  afterEach(() => {
    lspManager.dispose();
  });

  describe('startLSPForDocument', () => {
    it('should start LSP for a valid document', async () => {
      await lspManager.startLSPForDocument(mockDocument);
      expect(mockFunctions.showInfo).toHaveBeenCalledWith(
        expect.stringContaining('Language Server started')
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

    it('should handle errors when LSP executable is not found', async () => {
      mockFunctions.exec.mockImplementationOnce((cmd: string, callback: Function) => {
        callback(new Error('command not found'), null);
      });
      await lspManager.startLSPForDocument(mockDocument);
      expect(mockFunctions.showError).toHaveBeenCalledWith(
        expect.stringContaining('Failed to start')
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
      expect(mockFunctions.showInfo).toHaveBeenCalledWith(
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

    it('should wait between stop and start', async () => {
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: jest.fn().mockResolvedValue(undefined),
        needsStop: jest.fn().mockReturnValue(true),
      }));
      const startTime = Date.now();
      await lspManager.restartLSPForDocument(mockDocument);
      expect(Date.now() - startTime).toBeGreaterThanOrEqual(500);
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
      lspManager.dispose();
      await new Promise(resolve => setTimeout(resolve, 100));
      expect(mockStop).toHaveBeenCalledTimes(2);
    });

    it('should handle errors during dispose', async () => {
      const consoleSpy = jest.spyOn(console, 'error').mockImplementation();
      mockFunctions.languageClient.mockImplementation(() => ({
        start: jest.fn().mockResolvedValue(undefined),
        stop: jest.fn().mockRejectedValue(new Error('Dispose failed')),
        needsStop: jest.fn().mockReturnValue(true),
      }));
      await lspManager.startLSPForDocument(mockDocument);
      lspManager.dispose();
      await new Promise(resolve => setTimeout(resolve, 100));
      expect(consoleSpy).toHaveBeenCalledWith(
        expect.stringContaining('Error stopping'),
        expect.any(Error)
      );
      consoleSpy.mockRestore();
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
        { fileName: '/test/qe.in', content: '&CONTROL\n&SYSTEM\ncalculation = "scf"\npseudo_dir = "./"' },
        { fileName: '/test/gamess.inp', content: '$BASIS\n$CONTRL\n$SYSTEM\nruntyp=energy' },
        { fileName: '/test/nwchem.nw', content: 'geometry\n  H 0 0 0\nend\n\nbasis\n  * library 6-31G*\nend' },
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
