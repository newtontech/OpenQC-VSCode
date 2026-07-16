import * as vscode from 'vscode';
import { registerExportCommands } from '../../../src/commands/exportCommands';

describe('exportCommands', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    (vscode.commands.registerCommand as jest.Mock).mockImplementation(
      (_command: string, _handler: (...args: any[]) => unknown) => ({ dispose: jest.fn() })
    );
    (vscode.window.showQuickPick as jest.Mock).mockResolvedValue({
      label: 'XYZ',
      description: 'xyz',
    });
    (vscode.window.showSaveDialog as jest.Mock).mockResolvedValue({
      fsPath: '/tmp/bonded.xyz',
    });
  });

  it('shows exporter warnings after a successful command-palette export', async () => {
    const context = { subscriptions: [] } as unknown as vscode.ExtensionContext;
    registerExportCommands(context, () => ({
      atoms: [
        { element: 'C', x: 0, y: 0, z: 0 },
        { element: 'O', x: 1.2, y: 0, z: 0 },
      ],
      bonds: [{ from: 0, to: 1, order: 2 }],
    }));

    const handler = (vscode.commands.registerCommand as jest.Mock).mock.calls.find(
      ([command]) => command === 'openqc.exportStructure'
    )?.[1] as () => Promise<void>;

    await handler();

    expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
      expect.stringContaining('does not preserve edited bond topology or bond order')
    );
  });
});
