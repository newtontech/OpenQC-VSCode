/**
 * Unit tests for OpenQCViewerWebview (issue #107).
 * @module tests/unit/webviews/openqcViewerWebview.test
 */

import { OpenQCViewerWebview } from '../../../src/webviews/openqcViewerWebview';

// Mock vscode
jest.mock('vscode', () => ({
  Uri: {
    joinPath: jest.fn((base: any, ...parts: string[]) => ({
      fsPath: `${base.fsPath}/${parts.join('/')}`,
      path: `${base.path}/${parts.join('/')}`,
    })),
  },
}));

describe('OpenQCViewerWebview', () => {
  const mockWebview = {
    asWebviewUri: jest.fn((uri: any) => `vscode-webview://${uri.path}`),
    cspSource: 'https://test-host',
  };

  const mockExtensionUri = { fsPath: '/ext', path: '/ext' };

  beforeEach(() => {
    jest.clearAllMocks();
  });

  describe('generateHTML', () => {
    it('generates valid HTML with bundled assets', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('<!DOCTYPE html>');
      expect(html).toContain('<html');
      expect(html).toContain('</html>');
    });

    it('includes CSP with nonce', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('Content-Security-Policy');
      expect(html).toContain("script-src 'nonce-");
      expect(html).toContain("style-src 'nonce-");
    });

    it('loads bundled viewer JS', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('openqc-viewer.js');
    });

    it('loads bundled CSS', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('openqc-viewer.css');
    });

    it('loads 3Dmol from bundled assets, not CDN', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('3Dmol-min.js');
      expect(html).not.toContain('cdnjs.cloudflare.com');
      expect(html).not.toContain('unpkg.com');
    });

    it('includes viewer-canvas div', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="viewer-canvas"');
    });

    it('includes toolbar with style selector', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="style-select"');
      expect(html).toContain('ball-stick');
      expect(html).toContain('spacefill');
    });

    it('includes periodic controls', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="periodic-controls"');
      expect(html).toContain('id="sc-na"');
    });

    it('includes trajectory controls', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="trajectory-controls"');
      expect(html).toContain('id="frame-slider"');
    });

    it('includes WebGL fallback', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="webgl-fallback"');
    });

    it('includes status bar', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="status-bar"');
      expect(html).toContain('id="atom-count"');
    });
  });

  describe('getWebviewOptions', () => {
    it('enables scripts', () => {
      const options = OpenQCViewerWebview.getWebviewOptions(mockExtensionUri as any);
      expect(options.enableScripts).toBe(true);
    });

    it('sets local resource roots for media and 3dmol', () => {
      const options = OpenQCViewerWebview.getWebviewOptions(mockExtensionUri as any);
      expect(options.localResourceRoots).toHaveLength(2);
    });
  });
});
