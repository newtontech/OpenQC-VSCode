import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import * as vscode from 'vscode';
import {
  checkAnalyzer,
  executeAnalyzer,
  generateMultiwfnScript,
  OPENBABEL_CONFIG,
} from '../../../src/analyzers/ExternalAnalyzer';

describe('ExternalAnalyzer real process integration', () => {
  let tempDir: string;

  beforeEach(() => {
    tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-external-analyzer-'));
    jest.clearAllMocks();
  });

  afterEach(() => {
    fs.rmSync(tempDir, { recursive: true, force: true });
  });

  it('runs a fake analyzer executable and validates a relative output path under cwd', async () => {
    const executable = writeExecutable(
      'fake-c2x.sh',
      [
        '#!/bin/sh',
        'printf "converted %s -> %s\\n" "$1" "$2"',
        'printf "cube data\\n" > "$2"',
      ].join('\n')
    );

    const result = await executeAnalyzer(
      {
        executable,
        args: ['density.check', 'density.cube'],
        cwd: tempDir,
        expectedOutputPath: 'density.cube',
        description: 'fake c2x relative output smoke',
      },
      1000
    );

    expect(result).toMatchObject({
      success: true,
      exitCode: 0,
      stderr: '',
    });
    expect(result.stdout).toContain('converted density.check -> density.cube');
    expect(fs.readFileSync(path.join(tempDir, 'density.cube'), 'utf8')).toContain('cube data');
  });

  it('passes generated Multiwfn stdin to a real child process', async () => {
    const stdinCapture = path.join(tempDir, 'multiwfn.stdin');
    const outputFile = path.join(tempDir, 'density.cube');
    const executable = writeExecutable(
      'fake-multiwfn.sh',
      ['#!/bin/sh', `cat > "${stdinCapture}"`, `printf "cube\\n" > "${outputFile}"`].join('\n')
    );
    const command = generateMultiwfnScript('water.fchk', 'ELF cube', outputFile);
    command.executable = executable;
    command.cwd = tempDir;

    const result = await executeAnalyzer(command, 1000);

    expect(result.success).toBe(true);
    expect(fs.readFileSync(stdinCapture, 'utf8')).toBe(command.stdin);
    expect(fs.readFileSync(outputFile, 'utf8')).toContain('cube');
  });

  it('treats a configured non-executable analyzer path as unavailable', async () => {
    const notExecutable = path.join(tempDir, 'obabel-not-executable');
    fs.writeFileSync(notExecutable, '#!/bin/sh\nexit 0\n', { mode: 0o644 });
    mockWorkspaceSettings({
      'openqc.external.openBabelPath': notExecutable,
      'openqc.external.allowExternalAnalyzers': true,
    });

    const status = await checkAnalyzer(OPENBABEL_CONFIG);

    expect(status).toEqual({
      id: 'openbabel',
      available: false,
      path: notExecutable,
      enabled: true,
    });
  });

  function writeExecutable(filename: string, content: string): string {
    const executable = path.join(tempDir, filename);
    fs.writeFileSync(executable, `${content}\n`, { mode: 0o755 });
    return executable;
  }
});

function mockWorkspaceSettings(values: Record<string, unknown>): void {
  (vscode.workspace.getConfiguration as jest.Mock).mockImplementation((section?: string) => ({
    get: jest.fn((key: string, defaultValue?: unknown) => {
      const fullKey = section ? `${section}.${key}` : key;
      if (fullKey in values) {
        return values[fullKey];
      }
      if (key in values) {
        return values[key];
      }
      return defaultValue;
    }),
  }));
}
