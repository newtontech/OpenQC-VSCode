import { generateNonce } from '../../../src/utils/nonce';
import * as fs from 'fs';
import * as path from 'path';

const REPO_ROOT = path.resolve(__dirname, '../../..');
const WEBVIEW_NONCE_SOURCES = [
  'src/providers/StructureViewer.ts',
  'src/providers/DataPlotter.ts',
  'src/webviews/openqcViewerWebview.ts',
  'src/webviews/resultsWebview.ts',
  'src/sidebar/OpenQCConverterProvider.ts',
  'src/visualizers/ThreeJsWebview.ts',
  'src/visualizers/MoleculeViewerWebview.ts',
];

describe('generateNonce', () => {
  it('generates a 32-character CSP nonce by default', () => {
    expect(generateNonce()).toMatch(/^[A-Za-z0-9]{32}$/);
  });

  it('does not depend on Math.random for CSP nonce generation', () => {
    const randomSpy = jest.spyOn(Math, 'random').mockImplementation(() => {
      throw new Error('Math.random must not be used for CSP nonces');
    });

    expect(() => generateNonce()).not.toThrow();

    randomSpy.mockRestore();
  });

  it('supports custom nonce lengths for webview callers', () => {
    expect(generateNonce(48)).toMatch(/^[A-Za-z0-9]{48}$/);
  });

  it('keeps webview CSP nonce sources off Math.random', () => {
    for (const source of WEBVIEW_NONCE_SOURCES) {
      const text = fs.readFileSync(path.join(REPO_ROOT, source), 'utf8');

      expect(text).not.toContain('Math.random');
    }
  });
});
