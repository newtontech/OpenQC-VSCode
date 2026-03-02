import * as vscode from 'vscode';
import {
  MoleculeViewerPanel,
  MoleculeViewerState,
} from '../../../src/visualizers/MoleculeViewerPanel';

// Mock vscode module
jest.mock('vscode', () => ({
  Uri: {
    file: (path: string) => ({ fsPath: path, scheme: 'file' }),
  },
  window: {
    createWebviewPanel: jest.fn(),
    activeTextEditor: undefined,
    showSaveDialog: jest.fn(),
    showInformationMessage: jest.fn(),
    showErrorMessage: jest.fn(),
  },
  workspace: {
    fs: {
      writeFile: jest.fn(),
    },
  },
  ViewColumn: {
    One: 1,
    Two: 2,
  },
  WebviewPanel: class MockWebviewPanel {
    onDidDispose = jest.fn();
    onDidChangeViewState = jest.fn();
    webview = {
      postMessage: jest.fn(),
      onDidReceiveMessage: jest.fn(),
    };
    reveal = jest.fn();
    dispose = jest.fn();
  },
}));

describe('MoleculeViewerPanel', () => {
  let mockExtensionUri: vscode.Uri;
  let mockWebviewPanel: any;

  beforeEach(() => {
    mockExtensionUri = vscode.Uri.file('/path/to/extension');
    mockWebviewPanel = {
      webview: {
        postMessage: jest.fn(),
        onDidReceiveMessage: jest.fn(),
        html: '',
      },
      onDidDispose: jest.fn(),
      onDidChangeViewState: jest.fn(),
      reveal: jest.fn(),
      dispose: jest.fn(),
    };

    const createWebviewPanelMock = vscode.window.createWebviewPanel as jest.Mock;
    createWebviewPanelMock.mockReturnValue(mockWebviewPanel);

    // Clear any existing panel
    (MoleculeViewerPanel as any).currentPanel = undefined;
  });

  afterEach(() => {
    jest.clearAllMocks();
    (MoleculeViewerPanel as any).currentPanel = undefined;
  });

  describe('createOrShow', () => {
    it('creates a new webview panel when none exists', () => {
      const createSpy = vscode.window.createWebviewPanel as jest.Mock;

      MoleculeViewerPanel.createOrShow(mockExtensionUri, 'xyz-content', 'water.xyz');

      expect(createSpy).toHaveBeenCalledWith(
        'openqc.moleculeViewer',
        '3D Molecule Viewer',
        1, // ViewColumn.One when no active editor
        expect.objectContaining({
          enableScripts: true,
        })
      );
    });

    it('sends initialize message to webview on creation', () => {
      // Ensure postMessage is a jest mock function
      mockWebviewPanel.webview.postMessage = jest.fn();

      MoleculeViewerPanel.createOrShow(mockExtensionUri, 'xyz-content', 'water.xyz');

      // Check that postMessage was called
      expect(mockWebviewPanel.webview.postMessage).toHaveBeenCalled();
      const calls = mockWebviewPanel.webview.postMessage.mock.calls;
      const initCall = calls.find((call: any[]) => call[0] && call[0].type === 'initialize');
      expect(initCall).toBeDefined();
    });

    it('reveals existing panel instead of creating new one', () => {
      // First call creates the panel
      MoleculeViewerPanel.createOrShow(mockExtensionUri, 'xyz-content', 'water.xyz');
      const createSpy = vscode.window.createWebviewPanel as jest.Mock;
      createSpy.mockClear();

      // Second call should reveal existing panel (not create new one)
      MoleculeViewerPanel.createOrShow(mockExtensionUri, 'new-xyz', 'new.xyz');

      expect(createSpy).not.toHaveBeenCalled();
      // Note: reveal is called on the internal panel reference
      // We verify by checking that createWebviewPanel wasn't called again
    });
  });

  describe('handleMessage', () => {
    it('handles exportImage message and shows save dialog', async () => {
      const mockUri = vscode.Uri.file('/path/to/image.png');
      (vscode.window.showSaveDialog as jest.Mock).mockResolvedValue(mockUri);
      (vscode.workspace.fs.writeFile as jest.Mock).mockResolvedValue(undefined);

      const blob = new Blob(['fake-image-data'], { type: 'image/png' });
      const message = {
        type: 'exportImage',
        data: blob,
      };

      await MoleculeViewerPanel.handleMessage(message);

      expect(vscode.window.showSaveDialog).toHaveBeenCalledWith({
        filters: { Images: ['png'] },
        defaultUri: expect.any(Object),
      });
    });

    it('handles exportImage when user cancels save dialog', async () => {
      (vscode.window.showSaveDialog as jest.Mock).mockResolvedValue(undefined);

      const blob = new Blob(['fake-image-data'], { type: 'image/png' });
      const message = {
        type: 'exportImage',
        data: blob,
      };

      await MoleculeViewerPanel.handleMessage(message);

      // Should not try to write file when dialog is cancelled
      expect(vscode.workspace.fs.writeFile).not.toHaveBeenCalled();
    });

    it('handles unknown message types gracefully', () => {
      const message = {
        type: 'unknownType',
        data: 'some data',
      };

      // Should not throw error
      expect(() => {
        MoleculeViewerPanel.handleMessage(message);
      }).not.toThrow();
    });

    it('handles error messages', () => {
      const message = {
        type: 'error',
        message: 'Test error',
      };

      MoleculeViewerPanel.handleMessage(message);

      expect(vscode.window.showErrorMessage).toHaveBeenCalledWith('Molecule Viewer: Test error');
    });

    it('handles info messages', () => {
      const message = {
        type: 'info',
        message: 'Test info',
      };

      MoleculeViewerPanel.handleMessage(message);

      expect(vscode.window.showInformationMessage).toHaveBeenCalledWith('Test info');
    });
  });

  describe('MoleculeViewerState type', () => {
    it('has correct structure', () => {
      const state: MoleculeViewerState = {
        filename: 'test.xyz',
        xyz: 'test-xyz-content',
      };

      expect(state).toHaveProperty('filename');
      expect(state).toHaveProperty('xyz');
    });
  });

  describe('viewType', () => {
    it('has correct viewType', () => {
      expect(MoleculeViewerPanel.viewType).toBe('openqc.moleculeViewer');
    });
  });

  describe('dispose', () => {
    it('clears currentPanel and disposes resources', () => {
      MoleculeViewerPanel.createOrShow(mockExtensionUri, 'xyz-content', 'water.xyz');

      // Verify panel was created
      expect((MoleculeViewerPanel as any).currentPanel).toBeDefined();

      // Get reference before disposing
      const panel = (MoleculeViewerPanel as any).currentPanel;
      panel.dispose();

      // Verify panel was cleared
      expect((MoleculeViewerPanel as any).currentPanel).toBeUndefined();
      expect(mockWebviewPanel.dispose).toHaveBeenCalled();
    });

    it('disposes all tracked disposables', () => {
      // Create mock disposables
      const mockDisposable1 = { dispose: jest.fn() };
      const mockDisposable2 = { dispose: jest.fn() };

      MoleculeViewerPanel.createOrShow(mockExtensionUri, 'xyz-content', 'water.xyz');
      const panel = (MoleculeViewerPanel as any).currentPanel;

      // Add mock disposables to the internal array
      panel._disposables.push(mockDisposable1, mockDisposable2);

      panel.dispose();

      expect(mockDisposable1.dispose).toHaveBeenCalled();
      expect(mockDisposable2.dispose).toHaveBeenCalled();
    });
  });

  describe('serialize', () => {
    it('returns state with filename', () => {
      MoleculeViewerPanel.createOrShow(mockExtensionUri, 'xyz-content', 'test.xyz');

      const panel = (MoleculeViewerPanel as any).currentPanel;
      const state = panel.serialize();

      expect(state).toEqual({ filename: 'test.xyz' });
    });
  });

  describe('deserialize', () => {
    it('creates panel from deserialized state', () => {
      const state: MoleculeViewerState = {
        filename: 'deserialized.xyz',
        xyz: 'deserialized-xyz',
      };

      const createSpy = vscode.window.createWebviewPanel as jest.Mock;

      MoleculeViewerPanel.deserialize(state, mockExtensionUri);

      expect(createSpy).toHaveBeenCalled();
    });
  });
});
