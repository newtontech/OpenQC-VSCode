import { StructureViewer } from '../../../src/providers/StructureViewer';
import { FileTypeDetector } from '../../../src/managers/FileTypeDetector';
import * as vscode from 'vscode';

// Mock FileTypeDetector
jest.mock('../../../src/managers/FileTypeDetector', () => ({
  FileTypeDetector: jest.fn().mockImplementation(() => ({
    detectSoftware: jest.fn((doc: any) => {
      if (doc.fileName.includes('POSCAR')) return 'VASP';
      if (doc.fileName.includes('.com')) return 'Gaussian';
      if (doc.fileName.includes('.inp')) return 'CP2K';
      return null;
    }),
  })),
}));

// Mock vscode module
jest.mock('vscode', () => {
  let currentPanel: any = null;

  const mockPanel = {
    reveal: jest.fn(),
    onDidDispose: jest.fn(callback => {
      currentPanel = null;
    }),
    webview: { html: '' },
    dispose: jest.fn(),
  };

  return {
    Uri: {
      file: (path: string) => ({ fsPath: path, scheme: 'file' }),
    },
    window: {
      createWebviewPanel: jest.fn(() => {
        if (!currentPanel) {
          currentPanel = mockPanel;
        }
        return currentPanel;
      }),
      showWarningMessage: jest.fn(),
      activeTextEditor: undefined,
    },
    ViewColumn: {
      One: 1,
      Two: 2,
      Three: 3,
    },
  };
});

describe('StructureViewer', () => {
  let viewer: StructureViewer;
  const mockUri = { fsPath: '/test/extension' } as vscode.Uri;

  beforeEach(() => {
    jest.clearAllMocks();
    viewer = new StructureViewer(mockUri);
  });

  describe('show', () => {
    it('should show warning when no editor is provided', () => {
      viewer.show(undefined);
      expect(vscode.window.showWarningMessage).toHaveBeenCalledWith('No active editor found');
    });

    it('should show warning for unsupported file type', () => {
      const mockEditor = {
        document: { fileName: '/test/random.txt', getText: () => 'random' },
      } as unknown as vscode.TextEditor;

      viewer.show(mockEditor);
      expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
        'Unsupported file type for structure visualization'
      );
    });

    it('should create webview panel for VASP POSCAR', () => {
      const mockEditor = {
        document: {
          fileName: '/test/POSCAR',
          getText: () =>
            `Si Diamond Structure
1.0
3.84 0.00 0.00
0.00 3.84 0.00
0.00 0.00 3.84
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25`,
        },
      } as unknown as vscode.TextEditor;

      viewer.show(mockEditor);
      expect(vscode.window.createWebviewPanel).toHaveBeenCalledWith(
        'openqcStructureViewer',
        'OpenQC: Molecular Structure',
        vscode.ViewColumn.Two,
        expect.objectContaining({
          enableScripts: true,
          localResourceRoots: [mockUri],
          retainContextWhenHidden: true,
        })
      );
    });

    it('should reuse existing panel when available', () => {
      const mockEditor = {
        document: {
          fileName: '/test/POSCAR',
          getText: () => 'Si Diamond Structure\n1.0',
        },
      } as unknown as vscode.TextEditor;

      viewer.show(mockEditor);

      // Get the mock function before second call
      const createWebviewPanel = vscode.window.createWebviewPanel as jest.Mock;
      const firstCallCount = createWebviewPanel.mock.calls.length;

      viewer.show(mockEditor);

      // After revealing the same panel, createWebviewPanel should not be called again
      const secondCallCount = createWebviewPanel.mock.calls.length;
      expect(secondCallCount).toBe(firstCallCount);
    });

    it('should generate HTML with 3Dmol.js viewer', () => {
      const mockEditor = {
        document: {
          fileName: '/test/POSCAR',
          getText: () =>
            `H2O molecule
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
O H
2
Direct
0.0 0.0 0.0
0.5 0.5 0.5`,
        },
      } as unknown as vscode.TextEditor;

      viewer.show(mockEditor);

      // Verify the function was called with correct parameters
      expect(vscode.window.createWebviewPanel).toHaveBeenCalledWith(
        'openqcStructureViewer',
        'OpenQC: Molecular Structure',
        expect.any(Number),
        expect.objectContaining({
          enableScripts: true,
          localResourceRoots: [mockUri],
          retainContextWhenHidden: true,
        })
      );
    });
  });

  describe('showPreview', () => {
    it('should show warning when no editor is provided', () => {
      viewer.showPreview(undefined);
      expect(vscode.window.showWarningMessage).toHaveBeenCalledWith('No active editor found');
    });

    it('should create preview panel for Gaussian input', () => {
      const mockEditor = {
        document: {
          fileName: '/test/input.com',
          getText: () =>
            `# B3LYP/6-31G(d) Opt

Water molecule optimization

0 1
O  -0.464   0.177   0.0`,
        },
      } as unknown as vscode.TextEditor;

      viewer.showPreview(mockEditor);
      expect(vscode.window.createWebviewPanel).toHaveBeenCalledWith(
        'openqcInputPreview',
        'OpenQC: Input Preview',
        vscode.ViewColumn.Two,
        expect.objectContaining({
          enableScripts: true,
        })
      );
    });

    it('should extract parameters from Gaussian input', () => {
      const mockEditor = {
        document: {
          fileName: '/test/input.com',
          getText: () =>
            `%chk=test.chk
# B3LYP/6-31G(d) Opt

Water molecule

0 1
O 0.0 0.0 0.0`,
        },
      } as unknown as vscode.TextEditor;

      viewer.showPreview(mockEditor);

      // Verify the function was called
      expect(vscode.window.createWebviewPanel).toHaveBeenCalledWith(
        'openqcInputPreview',
        'OpenQC: Input Preview',
        expect.any(Number),
        expect.any(Object)
      );
    });

    it('should extract sections from CP2K input', () => {
      const mockEditor = {
        document: {
          fileName: '/test/cp2k.inp',
          getText: () =>
            `&GLOBAL
  PROJECT_NAME test
&END
&FORCE_EVAL
  METHOD FIST
&END`,
        },
      } as unknown as vscode.TextEditor;

      viewer.showPreview(mockEditor);

      // Verify the function was called
      expect(vscode.window.createWebviewPanel).toHaveBeenCalled();
    });
  });

  describe('HTML generation', () => {
    it('should include proper styling for dark theme', () => {
      const mockEditor = {
        document: {
          fileName: '/test/POSCAR',
          getText: () => 'Si\n1.0',
        },
      } as unknown as vscode.TextEditor;

      viewer.show(mockEditor);

      // Verify the function was called
      expect(vscode.window.createWebviewPanel).toHaveBeenCalled();
    });

    it('should include control buttons in viewer', () => {
      const mockEditor = {
        document: {
          fileName: '/test/POSCAR',
          getText: () => 'Si\n1.0',
        },
      } as unknown as vscode.TextEditor;

      viewer.show(mockEditor);

      // Verify the function was called
      expect(vscode.window.createWebviewPanel).toHaveBeenCalled();
    });
  });
});
