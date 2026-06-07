const mockConverter = {
  checkBackend: jest.fn(),
  getExtensionForFormat: jest.fn(),
  convert: jest.fn(),
  batchConvert: jest.fn(),
};

jest.mock(
  'vscode',
  () => ({
    commands: {
      registerCommand: jest.fn((_command: string, callback: (...args: unknown[]) => unknown) => ({
        command: _command,
        callback,
        dispose: jest.fn(),
      })),
    },
    window: {
      activeTextEditor: undefined,
      showErrorMessage: jest.fn(),
      showWarningMessage: jest.fn(),
      showInformationMessage: jest.fn(),
      showQuickPick: jest.fn(),
      showOpenDialog: jest.fn(),
      withProgress: jest.fn((_options, task) => task()),
    },
    workspace: {
      openTextDocument: jest.fn(),
      fs: {
        readFile: jest.fn(),
      },
    },
    env: {
      clipboard: {
        writeText: jest.fn(),
      },
      openExternal: jest.fn(),
    },
    Uri: {
      file: jest.fn((filePath: string) => ({ fsPath: filePath })),
    },
    ProgressLocation: {
      Notification: 15,
    },
  }),
  { virtual: true }
);

jest.mock('../../../src/converters', () => ({
  FormatConverter: jest.fn(() => mockConverter),
  SupportedFormat: {
    XYZ: 'xyz',
    PDB: 'pdb',
    CIF: 'cif',
    POSCAR: 'poscar',
    Gaussian: 'gaussian',
    ORCA: 'orca',
  },
}));

import * as vscode from 'vscode';
import { registerFormatConversionCommands } from '../../../src/commands/formatConversionCommands';

describe('format conversion commands', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    mockConverter.checkBackend.mockResolvedValue(true);
    mockConverter.getExtensionForFormat.mockImplementation((format: string) =>
      format === 'gaussian' ? 'gjf' : format
    );
    mockConverter.convert.mockResolvedValue({ success: true });
  });

  it('registers conversion commands with the extension context', () => {
    const context = { subscriptions: [] } as unknown as vscode.ExtensionContext;

    registerFormatConversionCommands(context);

    expect(vscode.commands.registerCommand).toHaveBeenCalledWith(
      'openqc.convertFormat',
      expect.any(Function)
    );
    expect(vscode.commands.registerCommand).toHaveBeenCalledWith(
      'openqc.convertToGaussian',
      expect.any(Function)
    );
    expect(context.subscriptions).toHaveLength(7);
  });

  it('passes the selected quick-convert target format to the backend', async () => {
    const context = { subscriptions: [] } as unknown as vscode.ExtensionContext;
    registerFormatConversionCommands(context);

    const convertToGaussian = (vscode.commands.registerCommand as jest.Mock).mock.calls.find(
      ([command]) => command === 'openqc.convertToGaussian'
    )?.[1] as (uri?: vscode.Uri) => Promise<void>;

    await convertToGaussian({ fsPath: '/tmp/water.xyz' } as vscode.Uri);

    expect(mockConverter.convert).toHaveBeenCalledWith(
      '/tmp/water.xyz',
      '/tmp/water.gjf',
      undefined,
      'gaussian'
    );
  });
});
