import * as vscode from 'vscode';
import { MoleculeViewerWebview } from '../../../src/visualizers/MoleculeViewerWebview';

describe('MoleculeViewerWebview', () => {
  const mockExtensionUri = vscode.Uri.file('/test/openqc');
  const mockWebview = {
    asWebviewUri: jest.fn((uri: vscode.Uri) => uri),
    cspSource: 'vscode-resource:',
  } as unknown as vscode.Webview;

  beforeEach(() => {
    jest.clearAllMocks();
  });

  describe('generateWebviewHTML', () => {
    it('generates HTML with local NGL Viewer script', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML(mockWebview, mockExtensionUri);

      expect(html).toContain('ngl');
      expect(html).toContain('<script');
      expect(html).toContain('</script>');
      expect(html).toContain('media/vendor/ngl/ngl.js');
      expect(html).not.toContain('https://unpkg.com');
      expect(html).not.toContain('https://cdn.jsdelivr.net');
    });

    it('includes VSCode webview API script', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML(mockWebview, mockExtensionUri);

      expect(html).toContain('acquireVsCodeApi');
    });

    it('includes canvas element for 3D rendering', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML(mockWebview, mockExtensionUri);

      expect(html).toContain('<canvas');
      expect(html).toContain('id="ngl-canvas"');
    });

    it('includes message listener for extension communication', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML(mockWebview, mockExtensionUri);

      expect(html).toContain('window.addEventListener');
      expect(html).toContain('message');
    });

    it('handles initialize message with structure data', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML(mockWebview, mockExtensionUri);

      expect(html).toContain('initialize');
      expect(html).toContain('loadStructure');
    });

    it('adds nonces to inline scripts and styles', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML(mockWebview, mockExtensionUri);

      expect(html).toMatch(/<style nonce="[A-Za-z0-9]{32}">/);
      expect(html).toMatch(/<script nonce="[A-Za-z0-9]{32}" src="[^"]*ngl\.js"><\/script>/);
      expect(html).toMatch(/<script nonce="[A-Za-z0-9]{32}">/);
    });
  });

  describe('getWebviewOptions', () => {
    it('returns correct webview options', () => {
      const options = MoleculeViewerWebview.getWebviewOptions(mockExtensionUri);

      expect(options).toHaveProperty('enableScripts', true);
      expect(options.localResourceRoots).toHaveLength(1);
      expect(options.localResourceRoots?.[0].toString()).toContain('media/vendor/ngl');
    });
  });

  describe('getCSP', () => {
    it('generates valid Content Security Policy', () => {
      const csp = MoleculeViewerWebview.getCSP('vscode-resource:', 'testnonce');

      expect(csp).toContain('default-src');
      expect(csp).toContain('script-src');
      expect(csp).toContain('style-src');
      expect(csp).toContain('img-src');
      expect(csp).toContain("'nonce-testnonce'");
      expect(csp).toContain('vscode-resource:');
    });

    it('does not allow remote CDN resources', () => {
      const csp = MoleculeViewerWebview.getCSP('vscode-resource:', 'testnonce');

      expect(csp).not.toContain('https://unpkg.com');
      expect(csp).not.toContain('https://cdn.jsdelivr.net');
      expect(csp).not.toContain('https:');
      expect(csp).not.toContain("'unsafe-inline'");
    });
  });
});
