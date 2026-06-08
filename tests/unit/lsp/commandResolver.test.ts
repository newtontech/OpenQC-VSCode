import {
  resolveLspCommand,
  readCommandOverrides,
  isExecutableAvailable,
} from '../../../src/lsp/commandResolver';
import { LSPServerRegistryEntry } from '../../../src/lsp/types';

// ---------------------------------------------------------------------------
// Mock vscode for readCommandOverrides
// ---------------------------------------------------------------------------

const mockConfigGet = jest.fn();

jest.mock('vscode', () => ({
  workspace: {
    getConfiguration: jest.fn(() => ({
      get: mockConfigGet,
    })),
  },
}));

// ---------------------------------------------------------------------------
// Test data
// ---------------------------------------------------------------------------

function makeEntry(overrides: Partial<LSPServerRegistryEntry> = {}): LSPServerRegistryEntry {
  return {
    id: 'gaussian-lsp',
    name: 'Gaussian',
    repository: 'newtontech/gaussian-lsp',
    executable: 'gaussian-lsp',
    languageId: 'gaussian',
    fileExtensions: ['gjf', 'com'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/gaussian-lsp',
    stability: 'stable',
    ...overrides,
  };
}

// ---------------------------------------------------------------------------
// resolveLspCommand
// ---------------------------------------------------------------------------

describe('resolveLspCommand', () => {
  const entry = makeEntry();

  it('uses registry executable when no overrides are provided', () => {
    const result = resolveLspCommand(entry, {});
    expect(result).toEqual({
      kind: 'pathOrCommand',
      command: 'gaussian-lsp',
      args: ['--stdio'],
      env: undefined,
    });
  });

  it('uses deprecated path override when command is not set', () => {
    const result = resolveLspCommand(entry, { path: '/opt/custom/gaussian-lsp' });
    expect(result).toEqual({
      kind: 'pathOrCommand',
      command: '/opt/custom/gaussian-lsp',
      args: ['--stdio'],
      env: undefined,
    });
  });

  it('prefers command over path when both are set', () => {
    const result = resolveLspCommand(entry, {
      path: '/opt/old/gaussian-lsp',
      command: '/opt/new/gaussian-lsp',
    });
    expect(result).toEqual(
      expect.objectContaining({
        kind: 'pathOrCommand',
        command: '/opt/new/gaussian-lsp',
      })
    );
  });

  it('uses custom args when provided', () => {
    const result = resolveLspCommand(entry, { args: ['--stdio', '--verbose'] });
    expect(result.args).toEqual(['--stdio', '--verbose']);
  });

  it('uses registry args before falling back to --stdio', () => {
    const result = resolveLspCommand(makeEntry({ args: ['serve', '--stdio'] }), {});
    expect(result.args).toEqual(['serve', '--stdio']);
  });

  it('passes through env variables', () => {
    const env = { PYTHONPATH: '/opt/lib', DEBUG: '1' };
    const result = resolveLspCommand(entry, { env });
    expect(result.env).toEqual(env);
  });

  it('omits env when empty or undefined', () => {
    expect(resolveLspCommand(entry, { env: undefined }).env).toBeUndefined();
    expect(resolveLspCommand(entry, { env: {} }).env).toBeUndefined();
  });

  it('handles paths with spaces', () => {
    const result = resolveLspCommand(entry, { path: '/opt/my tools/gaussian-lsp' });
    expect(result).toEqual(
      expect.objectContaining({
        kind: 'pathOrCommand',
        command: '/opt/my tools/gaussian-lsp',
      })
    );
  });

  it('defaults to --stdio args when no overrides are given', () => {
    const result = resolveLspCommand(entry, {});
    expect(result.args).toEqual(['--stdio']);
  });
});

// ---------------------------------------------------------------------------
// readCommandOverrides
// ---------------------------------------------------------------------------

describe('readCommandOverrides', () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  it('returns undefined for all fields when no config is set', () => {
    mockConfigGet.mockReturnValue(undefined);
    const vscode = require('vscode');
    const config = vscode.workspace.getConfiguration('openqc.lsp');
    const overrides = readCommandOverrides(config, 'gaussian');

    expect(overrides.path).toBeUndefined();
    expect(overrides.command).toBeUndefined();
    expect(overrides.args).toBeUndefined();
    expect(overrides.env).toBeUndefined();
  });

  it('reads deprecated path setting', () => {
    mockConfigGet.mockImplementation((key: string) => {
      if (key === 'gaussian.path') return '/opt/gaussian-lsp';
      return undefined;
    });
    const vscode = require('vscode');
    const config = vscode.workspace.getConfiguration('openqc.lsp');
    const overrides = readCommandOverrides(config, 'gaussian');
    expect(overrides.path).toBe('/opt/gaussian-lsp');
  });

  it('reads command, args, and env settings', () => {
    mockConfigGet.mockImplementation((key: string) => {
      if (key === 'cp2k.command') return '/usr/local/bin/cp2k-language-server';
      if (key === 'cp2k.args') return ['--stdio', '--log-level', 'debug'];
      if (key === 'cp2k.env') return { CP2K_ROOT: '/opt/cp2k' };
      return undefined;
    });
    const vscode = require('vscode');
    const config = vscode.workspace.getConfiguration('openqc.lsp');
    const overrides = readCommandOverrides(config, 'cp2k');

    expect(overrides.command).toBe('/usr/local/bin/cp2k-language-server');
    expect(overrides.args).toEqual(['--stdio', '--log-level', 'debug']);
    expect(overrides.env).toEqual({ CP2K_ROOT: '/opt/cp2k' });
  });
});

// ---------------------------------------------------------------------------
// isExecutableAvailable
// ---------------------------------------------------------------------------

describe('isExecutableAvailable', () => {
  it('returns false for a clearly nonexistent path', async () => {
    const result = await isExecutableAvailable('/no/such/path/binary');
    expect(result).toBe(false);
  });

  it('returns false for a nonexistent relative path', async () => {
    const result = await isExecutableAvailable('./nonexistent-binary-xyz');
    expect(result).toBe(false);
  });

  it('treats Windows-style paths as paths without shell lookup', async () => {
    const result = await isExecutableAvailable('C:\\no\\such\\binary.exe');
    expect(result).toBe(false);
  });

  it('returns false for a nonexistent command on PATH', async () => {
    const result = await isExecutableAvailable('no-such-executable-xyz-12345');
    expect(result).toBe(false);
  });

  it('returns true for node (commonly available on PATH)', async () => {
    const result = await isExecutableAvailable('node');
    expect(result).toBe(true);
  });
});
