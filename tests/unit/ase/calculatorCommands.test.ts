const mockRegisterCommand = jest.fn(
  (command: string, callback: (...args: unknown[]) => unknown) => ({
    command,
    callback,
    dispose: jest.fn(),
  })
);
const mockShowInputBox = jest.fn();
const mockShowQuickPick = jest.fn();
const mockShowWarningMessage = jest.fn();
const mockShowErrorMessage = jest.fn();
const mockShowInformationMessage = jest.fn();
const mockWithProgress = jest.fn((_options, task) => task({ report: jest.fn() }));
const mockOpenTextDocument = jest.fn();
const mockShowTextDocument = jest.fn();
const mockParseOutput = jest.fn();

jest.mock(
  'vscode',
  () => ({
    commands: {
      registerCommand: mockRegisterCommand,
    },
    window: {
      activeTextEditor: undefined,
      showInputBox: mockShowInputBox,
      showQuickPick: mockShowQuickPick,
      showWarningMessage: mockShowWarningMessage,
      showErrorMessage: mockShowErrorMessage,
      showInformationMessage: mockShowInformationMessage,
      showTextDocument: mockShowTextDocument,
      withProgress: mockWithProgress,
    },
    workspace: {
      getConfiguration: jest.fn(() => ({
        get: jest.fn((_key: string, defaultValue?: unknown) => defaultValue),
      })),
      openTextDocument: mockOpenTextDocument,
    },
    ProgressLocation: {
      Notification: 15,
    },
  }),
  { virtual: true }
);

const mockCalculator = {
  isAvailable: jest.fn(),
  isCalculatorAvailable: jest.fn(),
  generateInput: jest.fn(),
  runCalculation: jest.fn(),
  readResults: jest.fn(),
  calculate: jest.fn(),
};
const mockCreateCalculator = jest.fn(() => mockCalculator);

jest.mock('../../../src/ase/ASECalculator', () => ({
  CalculatorFactory: jest.fn().mockImplementation(() => ({
    createCalculator: mockCreateCalculator,
  })),
  CalculatorType: {
    VASP: 'vasp',
    CP2K: 'cp2k',
    QE: 'qe',
  },
  getCalculatorConfiguration: jest.fn(() => ({
    vasp: { command: 'vasp', defaultParams: {} },
    cp2k: { command: 'cp2k.popt', defaultParams: {} },
    qe: { command: 'pw.x', defaultParams: {} },
  })),
}));

const mockReadToAtoms = jest.fn();

jest.mock('../../../src/ase/ASEConverter', () => ({
  ASEConverter: jest.fn().mockImplementation(() => ({
    readToAtoms: mockReadToAtoms,
  })),
}));

jest.mock('../../../src/results/OutputParserBridge', () => ({
  parseOutput: mockParseOutput,
}));

import { registerCalculatorCommands } from '../../../src/ase/calculatorCommands';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

describe('calculator commands', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    mockRegisterCommand.mockImplementation(
      (command: string, callback: (...args: unknown[]) => unknown) => ({
        command,
        callback,
        dispose: jest.fn(),
      })
    );
    mockCalculator.isAvailable.mockResolvedValue(true);
    mockCalculator.isCalculatorAvailable.mockResolvedValue(true);
    mockCalculator.runCalculation.mockResolvedValue({
      success: true,
      energy: -1.23,
      warnings: [],
      metadata: {},
      outputFiles: [],
    });
    mockCalculator.readResults.mockResolvedValue({
      success: true,
      energy: -2.34,
      warnings: [],
      metadata: {},
      outputFiles: ['/tmp/calc/OUTCAR'],
    });
    mockParseOutput.mockResolvedValue({
      success: true,
      data: {
        schemaVersion: 'openqc.results.v1',
        success: true,
        software: 'orca',
        finalEnergy: { value: -75.98, unit: 'hartree' },
        scfEnergies: [-75.98],
        warnings: [],
      },
    });
    mockOpenTextDocument.mockResolvedValue({ uri: { fsPath: 'results.json' } });
  });

  it('registers the calculator command surface', () => {
    const context = createContext();

    registerCalculatorCommands(context);

    expect(registeredCommandIds()).toEqual([
      'openqc.generateCalculatorInput',
      'openqc.runCalculation',
      'openqc.readResults',
      'openqc.quickVASP',
      'openqc.quickCP2K',
      'openqc.quickQE',
      'openqc.checkCalculators',
    ]);
    expect(context.subscriptions).toHaveLength(7);
  });

  it('stops execution with a degraded backend message when ASE is unavailable', async () => {
    registerCalculatorCommands(createContext());
    mockShowInputBox.mockResolvedValue('/tmp/calc');
    mockShowQuickPick
      .mockResolvedValueOnce({ label: 'VASP', type: 'vasp' })
      .mockResolvedValueOnce({ label: 'Use defaults', value: false });
    mockCalculator.isAvailable.mockResolvedValue(false);

    await commandHandler('openqc.runCalculation')();

    expect(mockShowWarningMessage).toHaveBeenCalledWith(
      expect.stringContaining('ASE Python backend is unavailable')
    );
    expect(mockCalculator.isCalculatorAvailable).not.toHaveBeenCalled();
    expect(mockCalculator.runCalculation).not.toHaveBeenCalled();
    expect(mockWithProgress).not.toHaveBeenCalled();
  });

  it('passes customized calculator parameters from JSON input into execution config', async () => {
    registerCalculatorCommands(createContext());
    mockShowInputBox
      .mockResolvedValueOnce('/tmp/calc')
      .mockResolvedValueOnce('{"encut":520,"kpts":[2,2,2]}');
    mockShowQuickPick
      .mockResolvedValueOnce({ label: 'VASP', type: 'vasp' })
      .mockResolvedValueOnce({ label: 'Customize parameters', value: true });

    await commandHandler('openqc.runCalculation')();

    expect(mockCreateCalculator).toHaveBeenCalledWith(
      expect.objectContaining({
        type: 'vasp',
        command: 'vasp',
        parameters: {
          encut: 520,
          kpts: [2, 2, 2],
        },
      })
    );
    expect(mockCalculator.runCalculation).toHaveBeenCalledWith('/tmp/calc');
  });

  it('stops calculator execution on invalid custom parameter JSON', async () => {
    registerCalculatorCommands(createContext());
    mockShowInputBox.mockResolvedValueOnce('/tmp/calc').mockResolvedValueOnce('{bad json');
    mockShowQuickPick
      .mockResolvedValueOnce({ label: 'VASP', type: 'vasp' })
      .mockResolvedValueOnce({ label: 'Customize parameters', value: true });

    await commandHandler('openqc.runCalculation')();

    expect(mockShowErrorMessage).toHaveBeenCalledWith(
      expect.stringContaining('Invalid calculator parameter JSON')
    );
    expect(mockCreateCalculator).not.toHaveBeenCalled();
    expect(mockCalculator.runCalculation).not.toHaveBeenCalled();
  });

  it('routes detected ORCA output directories through the output bridge', async () => {
    const tmpdir = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-orca-results-'));
    const outputPath = path.join(tmpdir, 'ORCA.out');
    fs.writeFileSync(
      outputPath,
      '* O   R   C   A *\nFINAL SINGLE POINT ENERGY     -75.980000\n',
      'utf8'
    );
    try {
      registerCalculatorCommands(createContext());
      mockShowInputBox.mockResolvedValue(tmpdir);

      await commandHandler('openqc.readResults')();

      expect(mockParseOutput).toHaveBeenCalledWith(outputPath, 'orca');
      expect(mockShowQuickPick).not.toHaveBeenCalled();
      expect(mockCreateCalculator).not.toHaveBeenCalled();
      expect(mockOpenTextDocument).toHaveBeenCalledWith(
        expect.objectContaining({
          language: 'json',
          content: expect.stringContaining('"software": "orca"'),
        })
      );
    } finally {
      fs.rmSync(tmpdir, { recursive: true, force: true });
    }
  });

  it('keeps the ASE calculator result reader for directories without bridge outputs', async () => {
    const tmpdir = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-vasp-results-'));
    try {
      registerCalculatorCommands(createContext());
      mockShowInputBox.mockResolvedValue(tmpdir);
      mockShowQuickPick.mockResolvedValueOnce({ label: 'VASP', type: 'vasp' });

      await commandHandler('openqc.readResults')();

      expect(mockParseOutput).not.toHaveBeenCalled();
      expect(mockCreateCalculator).toHaveBeenCalledWith({ type: 'vasp', parameters: {} });
      expect(mockCalculator.readResults).toHaveBeenCalledWith(tmpdir);
      expect(mockOpenTextDocument).toHaveBeenCalledWith(
        expect.objectContaining({
          language: 'json',
          content: expect.stringContaining('"energy": -2.34'),
        })
      );
    } finally {
      fs.rmSync(tmpdir, { recursive: true, force: true });
    }
  });
});

function createContext(): any {
  return {
    extensionPath: '/tmp/openqc-extension',
    subscriptions: [],
  };
}

function registeredCommandIds(): string[] {
  return mockRegisterCommand.mock.calls.map(([command]) => command);
}

function commandHandler(command: string): () => Promise<void> {
  const handler = mockRegisterCommand.mock.calls.find(
    ([registered]) => registered === command
  )?.[1];
  if (!handler) {
    throw new Error(`Command not registered: ${command}`);
  }
  return handler as () => Promise<void>;
}
