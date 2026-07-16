import {
  checkAnalyzer,
  executeAnalyzer,
  generateC2xCommand,
  generateMultiwfnScript,
  generateOpenBabelCommand,
  getPathLookupCommand,
  OPENBABEL_CONFIG,
  previewCommand,
} from '../../../src/analyzers/ExternalAnalyzer';
import * as vscode from 'vscode';
import { execFile } from 'child_process';
import { access } from 'fs/promises';

jest.mock('child_process', () => ({
  execFile: jest.fn(),
}));
jest.mock('fs/promises', () => ({
  access: jest.fn(),
}));

describe('ExternalAnalyzer', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    (access as unknown as jest.Mock).mockResolvedValue(undefined);
  });

  it('generates non-interactive Multiwfn commands with operation stdin', () => {
    const command = generateMultiwfnScript(
      '/tmp/water.fchk',
      'Electron density cube',
      '/tmp/water.cube'
    );

    expect(command.args).toEqual(['/tmp/water.fchk']);
    expect(command.stdin).toContain('/tmp/water.cube');
    expect(command.stdin).toContain('q');
    expect(command.expectedOutputPath).toBe('/tmp/water.cube');
  });

  it('previews command stdin before external execution', async () => {
    (vscode.window.showWarningMessage as jest.Mock).mockResolvedValue('Execute');
    const command = generateMultiwfnScript('/tmp/water.fchk', 'ELF cube', '/tmp/elf.cube');

    const accepted = await previewCommand(command);

    expect(accepted).toBe(true);
    expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
      'Confirm external analyzer execution?',
      expect.objectContaining({
        detail: expect.stringContaining('Standard input:'),
      }),
      'Execute',
      'Cancel'
    );
  });

  it('writes analyzer stdin to the spawned process', async () => {
    const stdin = { write: jest.fn(), end: jest.fn() };
    (execFile as unknown as jest.Mock).mockImplementation((_exe, _args, _opts, cb) => {
      cb(null, 'done', '');
      return { stdin };
    });
    const command = generateMultiwfnScript(
      '/tmp/water.fchk',
      'Population analysis',
      '/tmp/pop.out'
    );

    const result = await executeAnalyzer(command, 1000);

    expect(result.success).toBe(true);
    expect(stdin.write).toHaveBeenCalledWith(command.stdin);
    expect(stdin.end).toHaveBeenCalled();
  });

  it('fails successful analyzer exits when the expected output file is missing', async () => {
    (access as unknown as jest.Mock).mockRejectedValue(new Error('missing'));
    (execFile as unknown as jest.Mock).mockImplementation((_exe, _args, _opts, cb) => {
      cb(null, 'done', '');
      return {};
    });
    const command = generateOpenBabelCommand('/tmp/water.xyz', '/tmp/water.pdb');

    const result = await executeAnalyzer(command, 1000);

    expect(result.success).toBe(false);
    expect(result.exitCode).toBe(0);
    expect(result.stderr).toContain('Expected analyzer output was not created');
  });

  it('keeps c2x conversion as argv-only execution', () => {
    const command = generateC2xCommand('/tmp/density.check', 'convert', '/tmp/density.cube');

    expect(command.args).toEqual(['/tmp/density.check', '/tmp/density.cube']);
    expect(command.stdin).toBeUndefined();
    expect(command.expectedOutputPath).toBe('/tmp/density.cube');
  });

  it('generates Open Babel conversion commands with output inference', () => {
    const command = generateOpenBabelCommand('/tmp/water.xyz', '/tmp/water.pdb');

    expect(command).toMatchObject({
      executable: 'obabel',
      args: ['/tmp/water.xyz', '-O', '/tmp/water.pdb'],
      expectedOutputPath: '/tmp/water.pdb',
      description: expect.stringContaining('Open Babel'),
    });
    expect(command.stdin).toBeUndefined();
  });

  it('uses the Windows where command for PATH analyzer lookup', () => {
    expect(getPathLookupCommand('obabel', 'win32')).toEqual({
      executable: 'where',
      args: ['obabel'],
    });
  });

  it('uses the POSIX which command for PATH analyzer lookup', () => {
    expect(getPathLookupCommand('obabel', 'darwin')).toEqual({
      executable: 'which',
      args: ['obabel'],
    });
  });

  it('reports PATH analyzer lookup using the first resolved executable path', async () => {
    (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
      get: jest.fn((key: string, defaultValue?: unknown) => {
        if (key === 'openqc.external.openBabelPath') {
          return '';
        }
        if (key === 'openqc.external.allowExternalAnalyzers') {
          return true;
        }
        return defaultValue;
      }),
    });
    (execFile as unknown as jest.Mock).mockImplementation((_exe, _args, cb) => {
      cb(null, '/usr/local/bin/obabel\n/opt/homebrew/bin/obabel\n', '');
    });

    const status = await checkAnalyzer(OPENBABEL_CONFIG);

    expect(status).toEqual({
      id: 'openbabel',
      available: true,
      path: '/usr/local/bin/obabel',
      enabled: true,
    });
  });
});
