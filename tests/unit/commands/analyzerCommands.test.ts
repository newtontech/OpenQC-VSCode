const mockCheckAnalyzer = jest.fn();
const mockGenerateMultiwfnScript = jest.fn();
const mockGenerateC2xCommand = jest.fn();
const mockGenerateOpenBabelCommand = jest.fn();
const mockPreviewCommand = jest.fn();
const mockExecuteAnalyzer = jest.fn();

jest.mock('../../../src/analyzers/ExternalAnalyzer', () => ({
  MULTIWFN_CONFIG: {
    id: 'multiwfn',
    displayName: 'Multiwfn',
    settingKey: 'openqc.external.multiwfnPath',
    defaultCommand: 'Multiwfn',
    description: 'Wavefunction analysis',
  },
  C2X_CONFIG: {
    id: 'c2x',
    displayName: 'c2x',
    settingKey: 'openqc.external.c2xPath',
    defaultCommand: 'c2x',
    description: 'Density conversion',
  },
  OPENBABEL_CONFIG: {
    id: 'openbabel',
    displayName: 'Open Babel',
    settingKey: 'openqc.external.openBabelPath',
    defaultCommand: 'obabel',
    description: 'Structure and molecular format conversion',
  },
  checkAnalyzer: mockCheckAnalyzer,
  generateMultiwfnScript: mockGenerateMultiwfnScript,
  generateC2xCommand: mockGenerateC2xCommand,
  generateOpenBabelCommand: mockGenerateOpenBabelCommand,
  previewCommand: mockPreviewCommand,
  executeAnalyzer: mockExecuteAnalyzer,
}));

import * as vscode from 'vscode';
import { registerAnalyzerCommands } from '../../../src/commands/analyzerCommands';

describe('analyzer commands', () => {
  const channel = {
    clear: jest.fn(),
    appendLine: jest.fn(),
    show: jest.fn(),
    dispose: jest.fn(),
  };
  const externalSettings: Record<string, unknown> = {
    allowExternalAnalyzers: false,
    timeoutMs: 60000,
  };

  beforeEach(() => {
    jest.clearAllMocks();
    externalSettings.allowExternalAnalyzers = false;
    externalSettings.timeoutMs = 60000;

    (vscode.commands.registerCommand as jest.Mock).mockImplementation(
      (command: string, callback: (...args: unknown[]) => unknown) => ({
        command,
        callback,
        dispose: jest.fn(),
      })
    );
    (vscode.workspace.getConfiguration as jest.Mock).mockImplementation((section?: string) => ({
      get: jest.fn((key: string, defaultValue?: unknown) => {
        if (section === 'openqc.external' && key in externalSettings) {
          return externalSettings[key];
        }
        return defaultValue;
      }),
    }));
    (vscode.window.createOutputChannel as jest.Mock).mockReturnValue(channel);
    (vscode.window.withProgress as jest.Mock).mockImplementation((_options, task) => task());
    (vscode.workspace.openTextDocument as jest.Mock).mockResolvedValue({
      uri: { fsPath: 'untitled' },
    });
    mockPreviewCommand.mockResolvedValue(true);
  });

  it('registers external analyzer commands', () => {
    const context = { subscriptions: [] } as unknown as vscode.ExtensionContext;

    registerAnalyzerCommands(context);

    expect(vscode.commands.registerCommand).toHaveBeenCalledWith(
      'openqc.checkExternalAnalyzers',
      expect.any(Function)
    );
    expect(vscode.commands.registerCommand).toHaveBeenCalledWith(
      'openqc.generateMultiwfnScript',
      expect.any(Function)
    );
    expect(vscode.commands.registerCommand).toHaveBeenCalledWith(
      'openqc.runMultiwfnAnalysis',
      expect.any(Function)
    );
    expect(vscode.commands.registerCommand).toHaveBeenCalledWith(
      'openqc.convertDensityC2x',
      expect.any(Function)
    );
    expect(vscode.commands.registerCommand).toHaveBeenCalledWith(
      'openqc.convertStructureOpenBabel',
      expect.any(Function)
    );
    expect(context.subscriptions).toHaveLength(5);
  });

  it('reports analyzer availability and global enablement in the output channel', async () => {
    externalSettings.allowExternalAnalyzers = true;
    mockCheckAnalyzer
      .mockResolvedValueOnce({
        id: 'multiwfn',
        available: true,
        path: '/opt/Multiwfn',
        enabled: true,
      })
      .mockResolvedValueOnce({
        id: 'c2x',
        available: false,
        path: null,
        enabled: true,
      })
      .mockResolvedValueOnce({
        id: 'openbabel',
        available: true,
        path: '/opt/obabel',
        enabled: true,
      });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.checkExternalAnalyzers')();

    expect(channel.appendLine).toHaveBeenCalledWith(expect.stringContaining('multiwfn'));
    expect(channel.appendLine).toHaveBeenCalledWith(expect.stringContaining('/opt/Multiwfn'));
    expect(channel.appendLine).toHaveBeenCalledWith(expect.stringContaining('c2x'));
    expect(channel.appendLine).toHaveBeenCalledWith(expect.stringContaining('openbabel'));
    expect(channel.appendLine).toHaveBeenCalledWith(expect.stringContaining('/opt/obabel'));
    expect(channel.appendLine).toHaveBeenCalledWith(
      expect.stringContaining('Analyzers are ENABLED')
    );
    expect(channel.show).toHaveBeenCalledWith(true);
  });

  it('generates a reviewable Multiwfn script without executing the analyzer', async () => {
    (vscode.window.showInputBox as jest.Mock).mockResolvedValue('/tmp/water.fchk');
    (vscode.window.showQuickPick as jest.Mock).mockResolvedValue('Population analysis');
    mockGenerateMultiwfnScript.mockReturnValue({
      executable: 'Multiwfn',
      args: ['/tmp/water.fchk'],
      cwd: '',
      stdin: '7\n1\n0\nq\n',
      expectedOutputPath: '/tmp/water.fchk.population_analysis.out',
      description: 'Multiwfn population analysis',
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.generateMultiwfnScript')();

    expect(mockGenerateMultiwfnScript).toHaveBeenCalledWith(
      '/tmp/water.fchk',
      'Population analysis',
      '/tmp/water.fchk.population_analysis.out'
    );
    expect(vscode.workspace.openTextDocument).toHaveBeenCalledWith(
      expect.objectContaining({
        language: 'shellscript',
        content: expect.stringContaining("<<'EOF'"),
      })
    );
    expect(mockExecuteAnalyzer).not.toHaveBeenCalled();
  });

  it('generates cube-oriented Multiwfn scripts with .cube output paths', async () => {
    (vscode.window.showInputBox as jest.Mock).mockResolvedValue('/tmp/water.fchk');
    (vscode.window.showQuickPick as jest.Mock).mockResolvedValue('ELF cube');
    mockGenerateMultiwfnScript.mockReturnValue({
      executable: 'Multiwfn',
      args: ['/tmp/water.fchk'],
      cwd: '',
      stdin: '5\n9\n2\n/tmp/water.fchk.elf_cube.cube\n0\nq\n',
      expectedOutputPath: '/tmp/water.fchk.elf_cube.cube',
      description: 'Multiwfn ELF cube',
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.generateMultiwfnScript')();

    expect(mockGenerateMultiwfnScript).toHaveBeenCalledWith(
      '/tmp/water.fchk',
      'ELF cube',
      '/tmp/water.fchk.elf_cube.cube'
    );
    expect(vscode.workspace.openTextDocument).toHaveBeenCalledWith(
      expect.objectContaining({
        content: expect.stringContaining('/tmp/water.fchk.elf_cube.cube'),
      })
    );
  });

  it('stops Multiwfn execution when external analyzers are disabled', async () => {
    mockCheckAnalyzer.mockResolvedValue({
      id: 'multiwfn',
      available: true,
      path: '/opt/Multiwfn',
      enabled: false,
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.runMultiwfnAnalysis')();

    expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
      expect.stringContaining('External analyzers are disabled')
    );
    expect(vscode.window.showInputBox).not.toHaveBeenCalled();
    expect(mockExecuteAnalyzer).not.toHaveBeenCalled();
  });

  it('stops Multiwfn execution when the executable is unavailable', async () => {
    mockCheckAnalyzer.mockResolvedValue({
      id: 'multiwfn',
      available: false,
      path: null,
      enabled: true,
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.runMultiwfnAnalysis')();

    expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
      expect.stringContaining('Multiwfn executable was not found')
    );
    expect(mockExecuteAnalyzer).not.toHaveBeenCalled();
  });

  it('runs Multiwfn with configured path, timeout, and failure output', async () => {
    externalSettings.timeoutMs = 1234;
    mockCheckAnalyzer.mockResolvedValue({
      id: 'multiwfn',
      available: true,
      path: '/opt/Multiwfn',
      enabled: true,
    });
    (vscode.window.showInputBox as jest.Mock).mockResolvedValue('/tmp/water.fchk');
    (vscode.window.showQuickPick as jest.Mock).mockResolvedValue('ELF cube');
    mockGenerateMultiwfnScript.mockReturnValue({
      executable: 'Multiwfn',
      args: ['/tmp/water.fchk'],
      cwd: '',
      stdin: '5\n9\n2\n/tmp/water.fchk.elf_cube.cube\n0\nq\n',
      expectedOutputPath: '/tmp/water.fchk.elf_cube.cube',
      description: 'Multiwfn ELF cube',
    });
    mockExecuteAnalyzer.mockResolvedValue({
      success: false,
      stdout: 'partial stdout',
      stderr: 'failed stderr',
      exitCode: 2,
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.runMultiwfnAnalysis')();

    expect(mockGenerateMultiwfnScript).toHaveBeenCalledWith(
      '/tmp/water.fchk',
      'ELF cube',
      '/tmp/water.fchk.elf_cube.cube'
    );
    expect(mockExecuteAnalyzer).toHaveBeenCalledWith(
      expect.objectContaining({ executable: '/opt/Multiwfn' }),
      1234
    );
    expect(channel.appendLine).toHaveBeenCalledWith('Exit code: 2');
    expect(channel.appendLine).toHaveBeenCalledWith('partial stdout');
    expect(channel.appendLine).toHaveBeenCalledWith('failed stderr');
    expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
      'External analyzer failed; see OpenQC External Analyzer Run output'
    );
  });

  it('runs c2x conversion with configured path and success message', async () => {
    mockCheckAnalyzer.mockResolvedValue({
      id: 'c2x',
      available: true,
      path: '/opt/c2x',
      enabled: true,
    });
    (vscode.window.showInputBox as jest.Mock)
      .mockResolvedValueOnce('/tmp/density.check')
      .mockResolvedValueOnce('/tmp/density.cube');
    mockGenerateC2xCommand.mockReturnValue({
      executable: 'c2x',
      args: ['/tmp/density.check', '/tmp/density.cube'],
      cwd: '',
      description: 'c2x conversion',
    });
    mockExecuteAnalyzer.mockResolvedValue({
      success: true,
      stdout: 'converted',
      stderr: '',
      exitCode: 0,
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.convertDensityC2x')();

    expect(mockExecuteAnalyzer).toHaveBeenCalledWith(
      expect.objectContaining({ executable: '/opt/c2x' }),
      60000
    );
    expect(vscode.window.showInformationMessage).toHaveBeenCalledWith(
      'External analyzer completed successfully'
    );
  });

  it('stops Open Babel conversion when external analyzers are disabled', async () => {
    mockCheckAnalyzer.mockResolvedValue({
      id: 'openbabel',
      available: true,
      path: '/opt/obabel',
      enabled: false,
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.convertStructureOpenBabel')();

    expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
      expect.stringContaining('External analyzers are disabled')
    );
    expect(vscode.window.showInputBox).not.toHaveBeenCalled();
    expect(mockExecuteAnalyzer).not.toHaveBeenCalled();
  });

  it('runs Open Babel conversion with configured path and preview confirmation', async () => {
    externalSettings.timeoutMs = 9876;
    mockCheckAnalyzer.mockResolvedValue({
      id: 'openbabel',
      available: true,
      path: '/opt/obabel',
      enabled: true,
    });
    (vscode.window.showInputBox as jest.Mock)
      .mockResolvedValueOnce('/tmp/water.xyz')
      .mockResolvedValueOnce('/tmp/water.pdb');
    mockGenerateOpenBabelCommand.mockReturnValue({
      executable: 'obabel',
      args: ['/tmp/water.xyz', '-O', '/tmp/water.pdb'],
      cwd: '',
      expectedOutputPath: '/tmp/water.pdb',
      description: 'Open Babel: convert /tmp/water.xyz → /tmp/water.pdb',
    });
    mockExecuteAnalyzer.mockResolvedValue({
      success: true,
      stdout: '1 molecule converted',
      stderr: '',
      exitCode: 0,
    });
    registerAnalyzerCommands({ subscriptions: [] } as unknown as vscode.ExtensionContext);

    await commandHandler('openqc.convertStructureOpenBabel')();

    expect(mockGenerateOpenBabelCommand).toHaveBeenCalledWith('/tmp/water.xyz', '/tmp/water.pdb');
    expect(mockPreviewCommand).toHaveBeenCalledWith(
      expect.objectContaining({ executable: '/opt/obabel' })
    );
    expect(mockExecuteAnalyzer).toHaveBeenCalledWith(
      expect.objectContaining({ executable: '/opt/obabel' }),
      9876
    );
    expect(channel.appendLine).toHaveBeenCalledWith('1 molecule converted');
  });
});

function commandHandler(command: string): () => Promise<void> {
  const handler = (vscode.commands.registerCommand as jest.Mock).mock.calls.find(
    ([registered]) => registered === command
  )?.[1];
  if (!handler) {
    throw new Error(`Command not registered: ${command}`);
  }
  return handler;
}
