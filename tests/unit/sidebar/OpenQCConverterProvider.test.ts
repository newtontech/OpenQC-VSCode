import { OpenQCConverterProvider } from '../../../src/sidebar/OpenQCConverterProvider';
import * as vscode from 'vscode';

// Mock vscode module
jest.mock('vscode', () => ({
  WebviewView: class WebviewView {
    webview = {
      options: {},
      html: '',
      onDidReceiveMessage: jest.fn(),
      postMessage: jest.fn(),
      asWebviewUri: jest.fn((uri: any) => uri),
      cspSource: 'vscode-webview://test',
    };
    visible = true;
    onDidDispose = jest.fn();
    show = jest.fn();
    hide = jest.fn();
  },
  Uri: {
    file: jest.fn((path: string) => ({
      path,
      fsPath: path,
      scheme: 'file',
      toString: () => path,
    })),
    joinPath: jest.fn((uri: any, ...paths: string[]) => ({
      path: [uri.path, ...paths].join('/'),
      fsPath: [uri.path, ...paths].join('/'),
      scheme: 'file',
      toString: () => [uri.path, ...paths].join('/'),
    })),
  },
  commands: {
    executeCommand: jest.fn(),
  },
}));

describe('OpenQCConverterProvider', () => {
  let provider: OpenQCConverterProvider;
  let mockExtensionUri: vscode.Uri;
  let messageHandler: Function | undefined;
  let mockWebviewView: any;
  let mockContext: any;
  let mockToken: any;

  beforeEach(() => {
    jest.clearAllMocks();
    mockExtensionUri = {
      path: '/test/extension',
      fsPath: '/test/extension',
      scheme: 'file',
      toString: () => '/test/extension',
    } as vscode.Uri;
    provider = new OpenQCConverterProvider(mockExtensionUri);

    messageHandler = undefined;
    mockWebviewView = {
      webview: {
        options: {},
        html: '',
        onDidReceiveMessage: jest.fn((handler: Function) => {
          messageHandler = handler;
        }),
        postMessage: jest.fn(),
        asWebviewUri: jest.fn((uri: any) => uri),
        cspSource: 'vscode-webview://test',
      },
      visible: true,
      onDidDispose: jest.fn(),
      show: jest.fn(),
      hide: jest.fn(),
    };
    mockContext = {
      state: {},
    };
    mockToken = {
      isCancellationRequested: false,
      onCancellationRequested: jest.fn(),
    };
  });

  describe('constructor', () => {
    it('should create provider with correct view type', () => {
      expect(OpenQCConverterProvider.viewType).toBe('openqc.converterSidebar');
      expect(provider).toBeInstanceOf(OpenQCConverterProvider);
    });

    it('should store extension URI', () => {
      expect(provider).toBeDefined();
    });
  });

  describe('resolveWebviewView', () => {
    it('should set webview options correctly', () => {
      provider.resolveWebviewView(mockWebviewView, mockContext, mockToken);
      expect(mockWebviewView.webview.options).toEqual({
        enableScripts: true,
        localResourceRoots: [mockExtensionUri],
      });
    });

    it('should set webview HTML content', () => {
      provider.resolveWebviewView(mockWebviewView, mockContext, mockToken);
      expect(mockWebviewView.webview.html).toContain('<!DOCTYPE html>');
      expect(mockWebviewView.webview.html).toContain('OpenQC ASE Integration');
    });

    it('should register message handler', () => {
      provider.resolveWebviewView(mockWebviewView, mockContext, mockToken);
      expect(mockWebviewView.webview.onDidReceiveMessage).toHaveBeenCalled();
      expect(messageHandler).toBeDefined();
    });
  });

  describe('message handling', () => {
    beforeEach(() => {
      provider.resolveWebviewView(mockWebviewView, mockContext, mockToken);
    });

    it('should handle convertToASE message', () => {
      expect(messageHandler).toBeDefined();
      messageHandler!({ type: 'convertToASE' });
      expect(vscode.commands.executeCommand).toHaveBeenCalledWith('openqc.convertToASE');
    });

    it('should handle convertFromASE message', () => {
      messageHandler!({ type: 'convertFromASE' });
      expect(vscode.commands.executeCommand).toHaveBeenCalledWith('openqc.convertFromASE');
    });

    it('should handle migrateFormat message', () => {
      messageHandler!({ type: 'migrateFormat' });
      expect(vscode.commands.executeCommand).toHaveBeenCalledWith('openqc.migrateFormat');
    });

    it('should handle quickConvert message with from and to formats', () => {
      messageHandler!({ type: 'quickConvert', from: 'vasp', to: 'cp2k' });
      expect(vscode.commands.executeCommand).toHaveBeenCalledWith(
        'openqc.quickConvert',
        'vasp',
        'cp2k'
      );
    });

    it('should handle openSettings message', () => {
      messageHandler!({ type: 'openSettings' });
      expect(vscode.commands.executeCommand).toHaveBeenCalledWith(
        'workbench.action.openSettings',
        'openqc'
      );
    });

    it('should ignore unknown message types', () => {
      expect(() => {
        messageHandler!({ type: 'unknownType' });
      }).not.toThrow();
      expect(vscode.commands.executeCommand).not.toHaveBeenCalled();
    });
  });

  describe('HTML content', () => {
    beforeEach(() => {
      provider.resolveWebviewView(mockWebviewView, mockContext, mockToken);
    });

    it('should include CSP nonce in HTML', () => {
      const html = mockWebviewView.webview.html;
      expect(html).toContain('nonce-');
      expect(html).toContain("script-src 'nonce-");
    });

    it('should include supported format tags', () => {
      const html = mockWebviewView.webview.html;
      expect(html).toContain('VASP');
      expect(html).toContain('CP2K');
      expect(html).toContain('Gaussian');
      expect(html).toContain('ORCA');
      expect(html).toContain('XYZ');
    });

    it('should include ASE converter buttons', () => {
      const html = mockWebviewView.webview.html;
      expect(html).toContain('Convert to ASE Atoms');
      expect(html).toContain('Convert from ASE Atoms');
      expect(html).toContain('Migrate Format');
    });

    it('should include quick convert UI elements', () => {
      const html = mockWebviewView.webview.html;
      expect(html).toContain('fromFormat');
      expect(html).toContain('toFormat');
    });

    it('should use VSCode CSS variables', () => {
      const html = mockWebviewView.webview.html;
      expect(html).toContain('--vscode-font-family');
      expect(html).toContain('--vscode-foreground');
      expect(html).toContain('--vscode-button-background');
    });

    it('should include JavaScript for message passing', () => {
      const html = mockWebviewView.webview.html;
      expect(html).toContain('acquireVsCodeApi');
      expect(html).toContain('vscode.postMessage');
    });
  });
});
