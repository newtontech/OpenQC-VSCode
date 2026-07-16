/**
 * Unit tests for PythonBridge (issue #74).
 * @module tests/unit/python/PythonBridge.test
 */

import { execPythonJson, checkBackend } from '../../../src/python/PythonBridge';
import type { BackendCheckResult } from '../../../src/python/PythonBackendStatus';

// Mock child_process
jest.mock('child_process', () => ({
  execFile: jest.fn(),
}));

// Mock vscode
jest.mock('vscode', () => ({
  workspace: {
    getConfiguration: jest.fn(() => ({
      get: jest.fn((key: string, def: any) => {
        if (key === 'pythonPath') return 'python3';
        if (key === 'timeoutMs') return 30000;
        return def;
      }),
    })),
  },
  window: {
    createOutputChannel: jest.fn(() => ({
      appendLine: jest.fn(),
      show: jest.fn(),
      dispose: jest.fn(),
    })),
  },
  CancellationTokenSource: class {
    token = { isCancellationRequested: false, onCancellationRequested: jest.fn() };
    cancel() {
      this.token.isCancellationRequested = true;
    }
  },
}));

// Mock Logger
jest.mock('../../../src/utils/Logger', () => ({
  Logger: {
    getInstance: jest.fn(() => ({
      debug: jest.fn(),
      info: jest.fn(),
      warn: jest.fn(),
      error: jest.fn(),
      setConfig: jest.fn(),
    })),
  },
}));

import { execFile } from 'child_process';
const mockExecFile = execFile as any as jest.Mock;

describe('PythonBridge', () => {
  beforeEach(() => {
    mockExecFile.mockReset();
  });

  describe('execPythonJson', () => {
    it('parses valid JSON response from Python', async () => {
      const mockData = { success: true, python: { version: '3.11.0' } };
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(null, JSON.stringify(mockData), '');
        return {} as any;
      });

      const result = await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'python3',
        timeoutMs: 5000,
      });

      expect(result.success).toBe(true);
      expect(result.data).toEqual(mockData);
    });

    it('adds the bundled Python package directory to PYTHONPATH', async () => {
      const mockData = { success: true };
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(null, JSON.stringify(mockData), '');
        return {} as any;
      });

      await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'python3',
        timeoutMs: 5000,
      });

      const options = mockExecFile.mock.calls[0][2];
      expect(options.env.PYTHONPATH).toContain('/python');
    });

    it('handles missing Python executable (ENOENT)', async () => {
      const error: any = new Error('spawn python3 ENOENT');
      error.code = 'ENOENT';
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(error, '', '');
        return {} as any;
      });

      const result = await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'bad-python',
        timeoutMs: 5000,
      });

      expect(result.success).toBe(false);
      expect(result.error?.message).toContain('not found');
    });

    it('handles subprocess timeout', async () => {
      const error: any = new Error('Timed out');
      error.killed = true;
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(error, '', '');
        return {} as any;
      });

      const result = await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'python3',
        timeoutMs: 1000,
      });

      expect(result.success).toBe(false);
      expect(result.timedOut).toBe(true);
    });

    it('handles invalid JSON response', async () => {
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(null, 'not valid json', '');
        return {} as any;
      });

      const result = await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'python3',
        timeoutMs: 5000,
      });

      expect(result.success).toBe(false);
      expect(result.error?.message).toContain('Invalid JSON');
    });

    it('handles non-zero exit code', async () => {
      const error: any = new Error('Command failed');
      error.code = 1;
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(error, '', 'ModuleNotFoundError: No module named xyz');
        return {} as any;
      });

      const result = await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'python3',
        timeoutMs: 5000,
      });

      expect(result.success).toBe(false);
      expect(result.exitCode).toBe(1);
    });

    it('handles empty output', async () => {
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(null, '', '');
        return {} as any;
      });

      const result = await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'python3',
        timeoutMs: 5000,
      });

      expect(result.success).toBe(false);
      expect(result.error?.message).toContain('empty output');
    });

    it('handles structured stderr error', async () => {
      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(null, 'not json', JSON.stringify({ error: { message: 'custom error' } }));
        return {} as any;
      });

      const result = await execPythonJson('openqc.bridge.check_backend', undefined, {
        pythonPath: 'python3',
        timeoutMs: 5000,
      });

      expect(result.success).toBe(false);
      expect(result.error?.message).toBe('custom error');
    });
  });

  describe('checkBackend', () => {
    it('returns BackendCheckResult on success', async () => {
      const mockData: BackendCheckResult = {
        success: true,
        python: {
          executable: '/usr/bin/python3',
          version: '3.11.0',
        },
        packages: {
          pymatgen: { available: true, version: '2024.1.0' },
          ase: { available: false, installHint: 'pip install ase' },
        },
        externalTools: {
          multiwfn: { available: false, path: null },
        },
      };

      mockExecFile.mockImplementation((_cmd: any, _args: any, _opts: any, cb: any) => {
        cb(null, JSON.stringify(mockData), '');
        return {} as any;
      });

      const result = await checkBackend({ pythonPath: 'python3' });

      expect(result.success).toBe(true);
      expect(result.data?.python.version).toBe('3.11.0');
      expect(result.data?.packages.pymatgen.available).toBe(true);
      expect(result.data?.packages.ase.available).toBe(false);
      expect(result.data?.packages.ase.installHint).toBe('pip install ase');
    });
  });
});
