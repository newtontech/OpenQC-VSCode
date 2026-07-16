import * as path from 'path';
import {
  buildVisualSmokePaths,
  isVisualSmokeReadyMarker,
  parseVisualSmokeArgs,
} from '../../../src/smoke/visualSmokePaths';

describe('visual smoke path options', () => {
  const workspaceRoot = path.resolve('/repo/openqc');

  it('uses source workspace paths by default', () => {
    const options = parseVisualSmokeArgs([], workspaceRoot);

    expect(options.extensionRoot).toBe(workspaceRoot);
    expect(options.outputPrefix).toBe('openqc-viewer-smoke');
    expect(options.outputDir).toBe(path.join(workspaceRoot, 'output', 'playwright'));
  });

  it('accepts an extracted VSIX extension root and output prefix', () => {
    const options = parseVisualSmokeArgs(
      [
        '--extension-root',
        '/tmp/openqc-vsix/extension',
        '--output-prefix',
        'openqc-viewer-vsix-smoke',
      ],
      workspaceRoot
    );

    expect(options.extensionRoot).toBe(path.resolve('/tmp/openqc-vsix/extension'));
    expect(options.outputPrefix).toBe('openqc-viewer-vsix-smoke');
    expect(buildVisualSmokePaths(options).htmlPath).toBe(
      path.join(workspaceRoot, 'output', 'playwright', 'openqc-viewer-vsix-smoke.html')
    );
  });

  it('rejects flags without values', () => {
    expect(() => parseVisualSmokeArgs(['--extension-root'], workspaceRoot)).toThrow(
      'Missing value for --extension-root'
    );
  });

  it('rejects unknown flags', () => {
    expect(() => parseVisualSmokeArgs(['--unknown'], workspaceRoot)).toThrow(
      'Unknown visual smoke option: --unknown'
    );
  });

  it('accepts only a verified ready marker after a browser polling timeout', () => {
    expect(isVisualSmokeReadyMarker({ state: 'ready', text: 'ready' })).toBe(true);
    expect(isVisualSmokeReadyMarker({ state: 'error', text: 'viewer failed' })).toBe(false);
    expect(isVisualSmokeReadyMarker({ state: 'booting', text: 'booting' })).toBe(false);
    expect(isVisualSmokeReadyMarker({ state: undefined, text: null })).toBe(false);
  });
});
