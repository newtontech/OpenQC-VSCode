const mockCheckBackend = jest.fn();
const mockLogger = {
  debug: jest.fn(),
  info: jest.fn(),
  warn: jest.fn(),
  error: jest.fn(),
};

jest.mock('../../../src/python/PythonBridge', () => ({
  checkBackend: mockCheckBackend,
}));

jest.mock('../../../src/utils/Logger', () => ({
  Logger: {
    getInstance: jest.fn(() => mockLogger),
  },
}));

import * as vscode from 'vscode';
import { registerPythonBackendCommands } from '../../../src/commands/pythonBackendCommands';

describe('python backend commands', () => {
  const channel = {
    clear: jest.fn(),
    appendLine: jest.fn(),
    show: jest.fn(),
    dispose: jest.fn(),
  };

  beforeEach(() => {
    jest.clearAllMocks();
    (vscode.commands.registerCommand as jest.Mock).mockImplementation(
      (command: string, callback: (...args: unknown[]) => unknown) => ({
        command,
        callback,
        dispose: jest.fn(),
      })
    );
    (vscode.window.withProgress as jest.Mock).mockImplementation((_options, task) =>
      task(undefined, { isCancellationRequested: false, onCancellationRequested: jest.fn() })
    );
    (vscode.window.createOutputChannel as jest.Mock).mockReturnValue(channel);
  });

  it('shows degraded backend detail and labeled capability readiness', async () => {
    mockCheckBackend.mockResolvedValue({
      success: true,
      data: {
        success: true,
        status: 'degraded',
        statusDetail: 'Python is available, but these core packages are missing: ase, cclib.',
        python: {
          executable: '/usr/bin/python3',
          version: '3.11.0',
          platform: 'darwin',
        },
        packages: {
          numpy: { available: true, version: '2.0.0' },
          ase: { available: false, installHint: 'pip install ase' },
          cclib: { available: false, installHint: 'pip install cclib' },
        },
        externalTools: {
          multiwfn: { available: false, path: null },
        },
        capabilities: {
          formatConversion: {
            label: 'Format conversion',
            status: 'degraded',
            detail: 'dpdata is available but ASE is missing.',
            requires: ['ase'],
          },
          outputParsing: {
            label: 'Calculation output parsing',
            status: 'missing',
            detail: 'cclib parses output files.',
            requires: ['cclib'],
          },
        },
        missingPackages: ['ase', 'cclib'],
      },
    });

    const context = { subscriptions: [] } as unknown as vscode.ExtensionContext;
    registerPythonBackendCommands(context);

    await commandHandler('openqc.checkPythonBackend')();

    expect(channel.appendLine).toHaveBeenCalledWith('Status: DEGRADED');
    expect(channel.appendLine).toHaveBeenCalledWith(
      'Detail: Python is available, but these core packages are missing: ase, cclib.'
    );
    expect(channel.appendLine).toHaveBeenCalledWith(
      'Capabilities: 0 available, 1 degraded, 1 missing'
    );
    expect(channel.appendLine).toHaveBeenCalledWith('   ⚠️ Format conversion: DEGRADED');
    expect(channel.appendLine).toHaveBeenCalledWith('      Install for full support: ase');
    expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
      'Python backend degraded: missing ase, cclib. 1/3 scientific packages available'
    );
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
