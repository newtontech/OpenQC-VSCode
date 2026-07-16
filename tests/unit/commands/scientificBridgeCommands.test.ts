jest.mock(
  'vscode',
  () => ({
    commands: {
      registerCommand: jest.fn((_command: string, callback: (...args: unknown[]) => unknown) => ({
        dispose: jest.fn(),
        callback,
      })),
    },
    window: {
      activeTextEditor: undefined,
      showErrorMessage: jest.fn(),
      showWarningMessage: jest.fn(),
      showInformationMessage: jest.fn(),
      showTextDocument: jest.fn(),
      showInputBox: jest.fn(),
      createWebviewPanel: jest.fn(() => ({
        webview: { html: '' },
      })),
      withProgress: jest.fn((_options, task) => task()),
    },
    workspace: {
      openTextDocument: jest.fn(),
    },
    ProgressLocation: {
      Notification: 15,
    },
    ViewColumn: {
      Two: 2,
    },
  }),
  { virtual: true }
);

jest.mock('../../../src/python/StructureBridge', () => ({
  parseStructure: jest.fn(),
  generateSupercell: jest.fn(),
}));

jest.mock('../../../src/results/OutputParserBridge', () => ({
  parseOutput: jest.fn(),
  extractTrajectory: jest.fn(),
}));

jest.mock('../../../src/visualizers/OpenQCViewerPanel', () => ({
  OpenQCViewerPanel: {
    createOrShow: jest.fn(),
  },
}));

import * as vscode from 'vscode';
import { registerScientificBridgeCommands } from '../../../src/commands/scientificBridgeCommands';
import { parseStructure, generateSupercell } from '../../../src/python/StructureBridge';
import { extractTrajectory, parseOutput } from '../../../src/results/OutputParserBridge';
import { OpenQCViewerPanel } from '../../../src/visualizers/OpenQCViewerPanel';

const mockParseStructure = parseStructure as jest.Mock;
const mockGenerateSupercell = generateSupercell as jest.Mock;
const mockParseOutput = parseOutput as jest.Mock;
const mockExtractTrajectory = extractTrajectory as jest.Mock;

describe('scientific bridge commands', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    (vscode.window as any).activeTextEditor = {
      document: {
        uri: { fsPath: '/tmp/POSCAR' },
        fileName: '/tmp/POSCAR',
        getText: jest.fn(() => 'Si'),
      },
    };
    (vscode.workspace.openTextDocument as jest.Mock).mockResolvedValue({
      uri: { fsPath: 'untitled' },
    });
    (vscode.window.showInputBox as jest.Mock).mockResolvedValue('2 2 1');
  });

  it('registers and invokes real bridge handlers', async () => {
    const context = { subscriptions: [], extensionUri: { fsPath: '/ext' } } as any;
    const structure = {
      name: 'Si',
      kind: 'periodic',
      atoms: [{ element: 'Si', x: 0, y: 0, z: 0 }],
    };
    mockParseStructure.mockResolvedValue({ success: true, data: structure });
    mockGenerateSupercell.mockResolvedValue({
      success: true,
      data: { ...structure, atoms: [...structure.atoms, ...structure.atoms] },
    });
    mockParseOutput.mockResolvedValue({
      success: true,
      data: { software: 'vasp', finalEnergy: { value: -1, unit: 'eV' }, scfEnergies: [-0.5, -1] },
    });
    mockExtractTrajectory.mockResolvedValue({
      success: true,
      data: {
        supported: true,
        frameCount: 2,
        frames: [
          {
            schemaVersion: 'openqc.structure.v1',
            kind: 'molecule',
            name: 'Frame 1',
            atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
          },
          {
            schemaVersion: 'openqc.structure.v1',
            kind: 'molecule',
            name: 'Frame 2',
            atoms: [{ element: 'H', x: 0, y: 0, z: 0.1 }],
          },
        ],
        energies: [-0.5, -1],
      },
    });

    registerScientificBridgeCommands(context);

    await commandHandler('openqc.parseStructurePython')();
    expect(mockParseStructure).toHaveBeenCalledWith('/tmp/POSCAR', undefined, 'Si');
    expect(vscode.workspace.openTextDocument).toHaveBeenCalledWith(
      expect.objectContaining({ language: 'json' })
    );

    await commandHandler('openqc.generateSupercell')();
    expect(mockGenerateSupercell).toHaveBeenCalledWith(structure, 2, 2, 1);
    expect(OpenQCViewerPanel.createOrShow).toHaveBeenCalled();

    await commandHandler('openqc.parseCalculationOutput')();
    expect(mockParseOutput).toHaveBeenCalledWith('/tmp/POSCAR');

    await commandHandler('openqc.showSCFConvergence')();
    expect(vscode.window.createWebviewPanel).toHaveBeenCalledWith(
      'openqc.scfConvergence',
      'OpenQC: SCF Convergence',
      vscode.ViewColumn.Two,
      expect.objectContaining({ enableScripts: false })
    );

    await commandHandler('openqc.showOptimizationTrajectory')();
    expect(mockExtractTrajectory).toHaveBeenCalledWith('/tmp/POSCAR');
    expect(OpenQCViewerPanel.createOrShow).toHaveBeenCalledWith(
      context.extensionUri,
      expect.objectContaining({
        schemaVersion: 'openqc.structure.v1',
        kind: 'trajectory',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
        frames: expect.arrayContaining([
          expect.objectContaining({ name: 'Frame 1' }),
          expect.objectContaining({ name: 'Frame 2' }),
        ]),
      }),
      '/tmp/POSCAR trajectory'
    );
  });

  it('shows a structured warning when an output has no extractable trajectory', async () => {
    const context = { subscriptions: [], extensionUri: { fsPath: '/ext' } } as any;
    mockExtractTrajectory.mockResolvedValue({
      success: true,
      data: {
        supported: false,
        frameCount: 0,
        frames: [],
        warnings: ['No atom coordinate trajectory was found in this output file'],
      },
    });
    registerScientificBridgeCommands(context);

    await commandHandler('openqc.showOptimizationTrajectory')();

    expect(OpenQCViewerPanel.createOrShow).not.toHaveBeenCalled();
    expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
      'No atom coordinate trajectory was found in this output file'
    );
  });

  it('does not announce success for unrecognized calculation output', async () => {
    const context = { subscriptions: [], extensionUri: { fsPath: '/ext' } } as any;
    mockParseOutput.mockResolvedValue({
      success: true,
      data: {
        schemaVersion: 'openqc.results.v1',
        software: 'auto',
        success: false,
        cclibAvailable: false,
        warnings: ['No calculation output data could be extracted from this file.'],
      },
    });
    registerScientificBridgeCommands(context);

    await commandHandler('openqc.parseCalculationOutput')();

    expect(vscode.workspace.openTextDocument).not.toHaveBeenCalled();
    expect(vscode.window.showInformationMessage).not.toHaveBeenCalledWith(
      expect.stringContaining('Parsed calculation output')
    );
    expect(vscode.window.showWarningMessage).toHaveBeenCalledWith(
      'No calculation output data could be extracted from this file.'
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
