/**
 * Complex Property Mapper Unit Tests
 *
 * Tests for complex property mapping between quantum chemistry codes
 */

import * as vscode from 'vscode';
import {
  ComplexPropertyMapper,
  initializeComplexPropertyMapper,
} from '../../../src/ase/ComplexPropertyMapper';
import { spawn } from 'child_process';
import * as path from 'path';

// Mock child_process
jest.mock('child_process');
const mockSpawn = spawn as jest.MockedFunction<typeof spawn>;
const mockExec = jest.fn();
jest.mock('child_process', () => ({
  spawn: jest.fn(),
  exec: jest
    .fn()
    .mockImplementation(
      (cmd: string, callback: (error: any, stdout: { stdout: string; stderr: string }) => void) => {
        return mockExec(cmd, callback);
      }
    ),
}));

// Mock vscode
jest.mock('vscode', () => ({
  workspace: {
    getConfiguration: jest.fn(),
  },
  Uri: {
    file: jest.fn((path: string) => ({ fsPath: path, scheme: 'file' })),
  },
}));

describe('ComplexPropertyMapper', () => {
  let mapper: ComplexPropertyMapper;
  let mockContext: jest.Mocked<vscode.ExtensionContext>;
  let mockProcess: any;

  beforeEach(() => {
    // Reset mocks
    jest.clearAllMocks();

    // Create mock context
    mockContext = {
      extensionPath: '/test/extension/path',
      subscriptions: [],
      globalState: {
        get: jest.fn(),
        update: jest.fn(),
        keys: jest.fn(),
        setKeysForSync: jest.fn(),
      },
      workspaceState: {
        get: jest.fn(),
        update: jest.fn(),
        keys: jest.fn(),
        setKeysForSync: jest.fn(),
      },
      asAbsolutePath: jest.fn(p => path.join('/test/extension/path', p)),
      asAbsolutePathUri: jest.fn(),
      extensionUri: { fsPath: '/test/extension/path', scheme: 'file' } as any,
      storageUri: { fsPath: '/test/storage', scheme: 'file' } as any,
      globalStorageUri: { fsPath: '/test/global-storage', scheme: 'file' } as any,
      logPath: '/test/logs',
    } as any;

    // Mock workspace configuration
    (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
      get: jest.fn((key: string, defaultValue: any) => defaultValue || 'python3'),
    });

    // Mock spawn process
    mockProcess = {
      stdout: {
        on: jest.fn(),
      },
      stderr: {
        on: jest.fn(),
      },
      on: jest.fn(),
    };
    mockSpawn.mockReturnValue(mockProcess);

    // Create mapper instance
    mapper = new ComplexPropertyMapper(mockContext);
  });

  describe('Constructor', () => {
    it('should initialize with default python path', () => {
      expect(mapper).toBeDefined();
      expect(vscode.workspace.getConfiguration).toHaveBeenCalledWith('openqc');
    });

    it('should use custom python path from configuration', () => {
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: jest.fn((key: string) => {
          if (key === 'pythonPath') return '/custom/python';
          return undefined;
        }),
      });

      mapper = new ComplexPropertyMapper(mockContext);
      expect(mapper).toBeDefined();
    });
  });

  describe('convertHubbard', () => {
    it('should convert Hubbard U parameters successfully', async () => {
      const hubbard = {
        enabled: true,
        method: 'Dudarev' as const,
        uValues: { Fe: 4.0, Co: 3.5 },
        jValues: { Fe: 0.9, Co: 0.8 },
      };

      const mockResult = {
        enabled: true,
        method: 'Dudarev',
        u_values: { Fe: 4.0, Co: 3.5 },
        j_values: { Fe: 0.9, Co: 0.8 },
      };

      // Setup spawn mock to call success callback
      mockSpawn.mockImplementation(() => {
        const proc = {
          stdout: {
            on: jest.fn((event: string, callback: (data: Buffer) => void) => {
              if (event === 'data') {
                callback(Buffer.from(JSON.stringify(mockResult)));
              }
            }),
          },
          stderr: {
            on: jest.fn(),
          },
          on: jest.fn((event: string, callback: (code: number) => void) => {
            if (event === 'close') {
              callback(0);
            }
          }),
        };
        return proc as any;
      });

      const result = await mapper.convertHubbard(hubbard, 'vasp', 'qe');

      expect(result.enabled).toBe(true);
      expect(result.method).toBe('Dudarev');
      expect(result.uValues).toEqual({ Fe: 4.0, Co: 3.5 });
      expect(result.jValues).toEqual({ Fe: 0.9, Co: 0.8 });
    });

    it('should handle conversion errors', async () => {
      const hubbard = {
        enabled: true,
        method: 'Dudarev' as const,
        uValues: { Fe: 4.0 },
      };

      mockSpawn.mockImplementation(() => {
        const proc = {
          stdout: { on: jest.fn() },
          stderr: {
            on: jest.fn((event: string, callback: (data: Buffer) => void) => {
              if (event === 'data') {
                callback(Buffer.from('Conversion error'));
              }
            }),
          },
          on: jest.fn((event: string, callback: (code: number) => void) => {
            if (event === 'close') {
              callback(1);
            }
          }),
        };
        return proc as any;
      });

      await expect(mapper.convertHubbard(hubbard, 'vasp', 'qe')).rejects.toThrow(
        'Mapper failed with code 1: Conversion error'
      );
    });

    it('should handle process errors', async () => {
      const hubbard = {
        enabled: true,
        method: 'Dudarev' as const,
        uValues: { Fe: 4.0 },
      };

      mockSpawn.mockImplementation(() => {
        const proc = {
          stdout: { on: jest.fn() },
          stderr: { on: jest.fn() },
          on: jest.fn((event: string, callback: (error: Error) => void) => {
            if (event === 'error') {
              callback(new Error('Process spawn failed'));
            }
          }),
        };
        return proc as any;
      });

      await expect(mapper.convertHubbard(hubbard, 'vasp', 'qe')).rejects.toThrow(
        'Failed to execute mapper: Process spawn failed'
      );
    });
  });

  describe('isAvailable', () => {
    it('should return true when ASE is available', async () => {
      const mockExec = require('child_process').exec;
      mockExec.mockImplementation(
        (
          cmd: string,
          callback: (error: any, result: { stdout: string; stderr: string }) => void
        ) => {
          callback(null, { stdout: '3.22.0\n', stderr: '' });
        }
      );

      const result = await mapper.isAvailable();

      expect(result).toBe(true);
    });

    it('should return false when ASE is not available', async () => {
      const mockExec = require('child_process').exec;
      mockExec.mockImplementation(
        (
          cmd: string,
          callback: (error: any, result: { stdout: string; stderr: string }) => void
        ) => {
          callback(new Error('ASE not found'), { stdout: '', stderr: 'ASE not found' });
        }
      );

      const result = await mapper.isAvailable();

      expect(result).toBe(false);
    });
  });

  describe('initializeComplexPropertyMapper', () => {
    it('should create and return a ComplexPropertyMapper instance', () => {
      const instance = initializeComplexPropertyMapper(mockContext);

      expect(instance).toBeInstanceOf(ComplexPropertyMapper);
    });

    it('should return valid ComplexPropertyMapper instance', () => {
      const instance1 = initializeComplexPropertyMapper(mockContext);
      expect(instance1).toBeInstanceOf(ComplexPropertyMapper);
    });
  });
});
