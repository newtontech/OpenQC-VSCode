import { MoleculeViewerWebview } from '../../../src/visualizers/MoleculeViewerWebview';

describe('MoleculeViewerWebview', () => {
  describe('generateWebviewHTML', () => {
    it('generates HTML with NGL Viewer script', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML();

      expect(html).toContain('ngl');
      expect(html).toContain('<script');
      expect(html).toContain('</script>');
    });

    it('includes VSCode webview API script', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML();

      expect(html).toContain('acquireVsCodeApi');
    });

    it('includes canvas element for 3D rendering', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML();

      expect(html).toContain('<canvas');
      expect(html).toContain('id="ngl-canvas"');
    });

    it('includes message listener for extension communication', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML();

      expect(html).toContain('window.addEventListener');
      expect(html).toContain('message');
    });

    it('handles initialize message with structure data', () => {
      const html = MoleculeViewerWebview.generateWebviewHTML();

      expect(html).toContain('initialize');
      expect(html).toContain('loadStructure');
    });
  });

  describe('getWebviewOptions', () => {
    it('returns correct webview options', () => {
      const options = MoleculeViewerWebview.getWebviewOptions();

      expect(options).toHaveProperty('enableScripts', true);
    });
  });

  describe('getCSP', () => {
    it('generates valid Content Security Policy', () => {
      const csp = MoleculeViewerWebview.getCSP();

      expect(csp).toContain('default-src');
      expect(csp).toContain('script-src');
      expect(csp).toContain('style-src');
      expect(csp).toContain('img-src');
    });

    it('allows NGL Viewer CDN', () => {
      const csp = MoleculeViewerWebview.getCSP();

      expect(csp).toContain('https://unpkg.com');
      expect(csp).toContain('https://cdn.jsdelivr.net');
    });
  });
});
