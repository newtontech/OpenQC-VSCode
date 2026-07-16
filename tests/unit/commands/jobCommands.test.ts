import * as vscode from 'vscode';
import { registerJobCommands, createJobExportPayload } from '../../../src/commands/jobCommands';
import { JobItem, JobTreeProvider } from '../../../src/sidebar/JobTreeProvider';

const runCalculation = jest.fn();
const createCalculator = jest.fn(() => ({ runCalculation }));

jest.mock('../../../src/ase/ASECalculator', () => ({
  CalculatorType: {
    VASP: 'vasp',
    CP2K: 'cp2k',
    QE: 'qe',
  },
  CalculatorFactory: jest.fn().mockImplementation(() => ({ createCalculator })),
  getCalculatorConfiguration: jest.fn(() => ({
    vasp: { command: '/usr/local/bin/vasp_std', defaultParams: { encut: 520 } },
    cp2k: { command: 'cp2k.popt', defaultParams: { cutoff: 400 } },
    qe: { command: 'pw.x', defaultParams: { ecutwfc: 40 } },
  })),
}));

describe('jobCommands', () => {
  let provider: jest.Mocked<JobTreeProvider>;

  beforeEach(() => {
    jest.clearAllMocks();
    jest.spyOn(Date, 'now').mockReturnValue(1234567890);
    (vscode.commands.registerCommand as jest.Mock).mockImplementation(
      (_command: string, _handler: (...args: any[]) => unknown) => ({ dispose: jest.fn() })
    );
    (vscode.window.showQuickPick as jest.Mock).mockResolvedValue({
      label: 'VASP',
      description: 'vasp',
    });
    (vscode.window.showOpenDialog as jest.Mock).mockResolvedValue([{ fsPath: '/tmp/vasp-run' }]);
    (vscode.window.showInputBox as jest.Mock).mockResolvedValue('Static VASP');
    (vscode.window.showSaveDialog as jest.Mock).mockResolvedValue({
      fsPath: '/tmp/static_vasp_results.json',
    });
    provider = {
      refresh: jest.fn(),
      addJob: jest.fn(),
      updateJobStatus: jest.fn(),
      updateJobResult: jest.fn(),
      cancelJob: jest.fn(() => true),
      restartJob: jest.fn(),
      getJob: jest.fn(),
    } as unknown as jest.Mocked<JobTreeProvider>;
  });

  afterEach(() => {
    jest.restoreAllMocks();
  });

  it('runs an existing calculator directory and stores the real successful result', async () => {
    runCalculation.mockResolvedValue({
      success: true,
      energy: -10.5,
      forces: [[0.1, 0.2, 0.3]],
      warnings: ['converged with warning'],
      metadata: { parser: 'ase' },
      outputFiles: ['OUTCAR'],
    });
    registerJobCommands({ subscriptions: [] } as any, provider);

    await getCommandHandler('openqc.sidebar.runCalculation')();

    expect(createCalculator).toHaveBeenCalledWith({
      type: 'vasp',
      command: '/usr/local/bin/vasp_std',
      parameters: { encut: 520 },
    });
    expect(runCalculation).toHaveBeenCalledWith('/tmp/vasp-run', undefined, expect.any(Object));
    expect(provider.addJob).toHaveBeenCalledWith(
      expect.objectContaining({
        id: 'job-1234567890',
        label: 'Static VASP',
        status: 'running',
        progress: 5,
        software: 'VASP',
        workingDirectory: '/tmp/vasp-run',
      })
    );
    expect(provider.updateJobResult).toHaveBeenCalledWith(
      'job-1234567890',
      'completed',
      expect.objectContaining({
        success: true,
        energy: -10.5,
        outputFiles: ['OUTCAR'],
      }),
      100
    );
  });

  it('stores failed calculator execution results instead of fabricating output', async () => {
    runCalculation.mockResolvedValue({
      success: false,
      error: 'vasp_std exited with code 1',
      warnings: [],
      metadata: {},
      outputFiles: ['stdout.log'],
    });
    registerJobCommands({ subscriptions: [] } as any, provider);

    await getCommandHandler('openqc.sidebar.runCalculation')();

    expect(provider.updateJobResult).toHaveBeenCalledWith(
      'job-1234567890',
      'failed',
      expect.objectContaining({
        success: false,
        error: 'vasp_std exited with code 1',
        outputFiles: ['stdout.log'],
      }),
      100
    );
    expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
      'Calculation failed: vasp_std exited with code 1'
    );
  });

  it('cancels the running calculator and ignores late results', async () => {
    let resolveRun: (value: unknown) => void;
    runCalculation.mockImplementation(
      () =>
        new Promise(resolve => {
          resolveRun = resolve;
        })
    );
    registerJobCommands({ subscriptions: [] } as any, provider);

    const runPromise = getCommandHandler('openqc.sidebar.runCalculation')();
    await flushPromises();

    const signal = runCalculation.mock.calls[0][2] as AbortSignal;
    expect(signal.aborted).toBe(false);
    (provider.getJob as jest.Mock).mockReturnValue(
      new JobItem('job-1234567890', 'Static VASP', 'cancelled', 25, 'VASP')
    );

    await getCommandHandler('openqc.sidebar.cancelJob')(
      new JobItem('job-1234567890', 'Static VASP', 'running', 25, 'VASP')
    );

    expect(signal.aborted).toBe(true);
    expect(provider.cancelJob).toHaveBeenCalledWith('job-1234567890');
    resolveRun!({
      success: true,
      energy: -99,
      warnings: [],
      metadata: {},
      outputFiles: ['OUTCAR'],
    });
    await runPromise;

    expect(provider.updateJobResult).not.toHaveBeenCalled();
  });

  it('restarts a job by creating and executing a concrete restart job', async () => {
    runCalculation.mockResolvedValue({
      success: true,
      energy: -11,
      warnings: [],
      metadata: {},
      outputFiles: ['OUTCAR'],
    });
    const restartJob = new JobItem(
      'job-old-restart-1234567890',
      'Old Job (restart)',
      'running',
      5,
      'VASP',
      undefined,
      undefined,
      '/tmp/vasp-run'
    );
    (provider.restartJob as jest.Mock).mockReturnValue(restartJob);
    (provider.getJob as jest.Mock).mockReturnValue(restartJob);
    registerJobCommands({ subscriptions: [] } as any, provider);

    await getCommandHandler('openqc.sidebar.restartJob')(
      new JobItem(
        'job-old',
        'Old Job',
        'completed',
        100,
        'VASP',
        undefined,
        undefined,
        '/tmp/vasp-run'
      )
    );

    expect(provider.restartJob).toHaveBeenCalledWith('job-old');
    expect(runCalculation).toHaveBeenCalledWith('/tmp/vasp-run', undefined, expect.any(Object));
    expect(provider.updateJobResult).toHaveBeenCalledWith(
      'job-old-restart-1234567890',
      'completed',
      expect.objectContaining({ energy: -11 }),
      100
    );
  });

  it('does not start a job when the calculation name prompt is cancelled', async () => {
    (vscode.window.showInputBox as jest.Mock).mockResolvedValue(undefined);
    registerJobCommands({ subscriptions: [] } as any, provider);

    await getCommandHandler('openqc.sidebar.runCalculation')();

    expect(provider.addJob).not.toHaveBeenCalled();
    expect(runCalculation).not.toHaveBeenCalled();
  });

  it('renders captured job result data in the results webview', async () => {
    registerJobCommands({ subscriptions: [] } as any, provider);
    const item = new JobItem(
      'job-real',
      'Real Job',
      'completed',
      100,
      'VASP',
      new Date('2026-07-03T00:00:00Z'),
      new Date('2026-07-03T00:00:01Z'),
      '/tmp/real-job',
      {
        success: true,
        energy: -10.5,
        warnings: ['low precision'],
        metadata: { source: 'OUTCAR' },
        outputFiles: ['OUTCAR'],
      }
    );

    await getCommandHandler('openqc.sidebar.viewResults')(item);

    const panel = (vscode.window.createWebviewPanel as jest.Mock).mock.results[0].value;
    expect(panel.webview.html).toContain('Final energy: -10.5');
    expect(panel.webview.html).toContain('Output files: OUTCAR');
    expect(panel.webview.html).not.toContain('Sample output data');
  });

  it('exports captured result data as JSON', async () => {
    registerJobCommands({ subscriptions: [] } as any, provider);
    const item = new JobItem(
      'job-export',
      'Export Job',
      'completed',
      100,
      'CP2K',
      undefined,
      undefined,
      '/tmp/cp2k',
      {
        success: true,
        energy: -123.4,
        outputFiles: ['openqc_cp2k_smoke.out'],
      }
    );

    await getCommandHandler('openqc.sidebar.exportData')(item);

    const buffer = (vscode.workspace.fs.writeFile as jest.Mock).mock.calls[0][1] as Buffer;
    const payload = JSON.parse(buffer.toString());
    expect(payload.result.energy).toBe(-123.4);
    expect(payload.result.outputFiles).toEqual(['openqc_cp2k_smoke.out']);
    expect(payload).not.toHaveProperty('data.energies');
  });

  it('builds export payloads from JobItem result data', () => {
    const item = new JobItem(
      'job-payload',
      'Payload',
      'completed',
      100,
      'QE',
      undefined,
      undefined,
      '/tmp/qe',
      {
        success: true,
        energy: -3.14,
      }
    );

    expect(createJobExportPayload(item)).toMatchObject({
      jobId: 'job-payload',
      workingDirectory: '/tmp/qe',
      result: { success: true, energy: -3.14 },
    });
  });
});

function getCommandHandler(command: string): (...args: any[]) => Promise<void> {
  const handler = (vscode.commands.registerCommand as jest.Mock).mock.calls.find(
    ([registeredCommand]) => registeredCommand === command
  )?.[1];
  if (!handler) {
    throw new Error(`Command not registered: ${command}`);
  }
  return handler;
}

function flushPromises(): Promise<void> {
  return new Promise(resolve => setImmediate(resolve));
}
