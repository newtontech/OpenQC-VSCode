const mockConverter = {
  readToAtoms: jest.fn(),
  writeFromAtoms: jest.fn(),
  convertFormat: jest.fn(),
  getSupportedFormats: jest.fn(),
};

jest.mock('../../../src/ase/ASEConverter', () => ({
  ASEFormat: {
    VASP: 'vasp',
    CP2K: 'cp2k',
    QE: 'qe',
    Gaussian: 'gaussian',
    ORCA: 'orca',
    NWChem: 'nwchem',
    GAMESS: 'gamess',
    LAMMPS: 'lammps',
    XYZ: 'xyz',
    PDB: 'pdb',
    CIF: 'cif',
  },
  ASEConverter: jest.fn(() => mockConverter),
}));

import * as vscode from 'vscode';
import { registerASECommands } from '../../../src/ase/commands';

describe('ASE command handlers', () => {
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
      task({ report: jest.fn() })
    );
    (vscode.workspace.openTextDocument as jest.Mock).mockResolvedValue({
      uri: { fsPath: '/tmp/output.xyz' },
    });
    (vscode.window.showTextDocument as jest.Mock).mockResolvedValue(undefined);
    (vscode.window.showInputBox as jest.Mock).mockResolvedValue('/tmp/water.xyz');
    (vscode.window.showQuickPick as jest.Mock).mockImplementation(
      async (items: any[]) => items.find(item => item.format === 'xyz') ?? items[0]
    );
    mockConverter.getSupportedFormats.mockReturnValue({
      vasp: { name: 'VASP', extensions: ['POSCAR'], description: 'VASP structure format' },
      cp2k: { name: 'CP2K', extensions: ['.inp'], description: 'CP2K input format' },
      qe: {
        name: 'Quantum ESPRESSO',
        extensions: ['.in'],
        description: 'Quantum ESPRESSO input format',
      },
      xyz: { name: 'XYZ', extensions: ['.xyz'], description: 'Generic XYZ format' },
    });
    mockConverter.writeFromAtoms.mockResolvedValue({
      success: true,
      warnings: [],
      metadata: {},
    });
    mockConverter.convertFormat.mockResolvedValue({
      success: true,
      warnings: [],
      metadata: {},
    });
  });

  it('normalizes ASE Atoms JSON before writing from Atoms', async () => {
    const context = {
      subscriptions: [],
      extensionPath: '/ext',
    } as unknown as vscode.ExtensionContext;
    registerASECommands(context);
    (vscode.window as any).activeTextEditor = {
      document: {
        getText: jest.fn(() =>
          JSON.stringify({
            data: {
              atoms: {
                numbers: [8, 1, 1],
                positions: [
                  [0, 0, 0],
                  [0.75, 0, 0.5],
                  [-0.75, 0, 0.5],
                ],
                cell: [12, 13, 14],
                pbc: true,
                info: { source: 'unit-test' },
              },
            },
          })
        ),
      },
    };

    await commandHandler('openqc.convertFromASE')();

    expect(mockConverter.writeFromAtoms).toHaveBeenCalledWith(
      {
        chemical_symbols: ['O', 'H', 'H'],
        positions: [
          [0, 0, 0],
          [0.75, 0, 0.5],
          [-0.75, 0, 0.5],
        ],
        cell: [
          [12, 0, 0],
          [0, 13, 0],
          [0, 0, 14],
        ],
        pbc: [true, true, true],
        info: { source: 'unit-test' },
      },
      '/tmp/water.xyz',
      'xyz'
    );
    expect(vscode.window.showInformationMessage).toHaveBeenCalledWith(
      'Successfully wrote /tmp/water.xyz'
    );
  });

  it('rejects malformed ASE Atoms JSON before invoking the converter', async () => {
    const context = {
      subscriptions: [],
      extensionPath: '/ext',
    } as unknown as vscode.ExtensionContext;
    registerASECommands(context);
    (vscode.window as any).activeTextEditor = {
      document: {
        getText: jest.fn(() =>
          JSON.stringify({
            aseAtoms: {
              chemical_symbols: ['H', 'H'],
              positions: [[0, 0, 0]],
            },
          })
        ),
      },
    };

    await commandHandler('openqc.convertFromASE')();

    expect(mockConverter.writeFromAtoms).not.toHaveBeenCalled();
    expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
      expect.stringContaining('chemical_symbols has 2 entries but positions has 1')
    );
  });

  it('validates quick-convert format arguments from the sidebar', async () => {
    const context = {
      subscriptions: [],
      extensionPath: '/ext',
    } as unknown as vscode.ExtensionContext;
    registerASECommands(context);
    (vscode.window as any).activeTextEditor = {
      document: {
        uri: { fsPath: '/tmp/POSCAR' },
      },
    };

    await commandHandler('openqc.quickConvert')('POSCAR', 'quantum espresso');

    expect(mockConverter.convertFormat).toHaveBeenCalledWith(
      '/tmp/POSCAR',
      '/tmp/POSCAR_qe.in',
      'vasp',
      'qe'
    );

    await commandHandler('openqc.quickConvert')('bad-format', 'cp2k');

    expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
      expect.stringContaining('Unsupported ASE quick conversion')
    );
  });
});

function commandHandler(command: string): (...args: unknown[]) => Promise<void> {
  const handler = (vscode.commands.registerCommand as jest.Mock).mock.calls.find(
    ([registered]) => registered === command
  )?.[1];
  if (!handler) {
    throw new Error(`Command not registered: ${command}`);
  }
  return handler;
}
