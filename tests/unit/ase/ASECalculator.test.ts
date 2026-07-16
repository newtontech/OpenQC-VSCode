import { EventEmitter } from 'events';
import { spawn } from 'child_process';
import { ASECalculator, CalculatorType } from '../../../src/ase/ASECalculator';

jest.mock('child_process', () => ({
  spawn: jest.fn(),
  exec: jest.fn(),
}));

jest.mock(
  'vscode',
  () => ({
    workspace: {
      getConfiguration: jest.fn(() => ({
        get: jest.fn((key: string, defaultValue?: unknown) => {
          if (key === 'pythonPath') {
            return 'python3';
          }
          return defaultValue;
        }),
      })),
    },
  }),
  { virtual: true }
);

describe('ASECalculator Python CLI contract', () => {
  const context = {
    extensionPath: '/ext',
  } as any;
  const atoms = {
    chemical_symbols: ['H', 'H'],
    positions: [
      [0, 0, 0],
      [0, 0, 0.74],
    ],
    cell: [
      [8, 0, 0],
      [0, 8, 0],
      [0, 0, 8],
    ],
    pbc: [false, false, false],
  };

  beforeEach(() => {
    jest.clearAllMocks();
  });

  it('sends generate as a single CalculatorInput JSON payload and maps returned files', async () => {
    mockSpawnJson({
      success: true,
      status: 'completed',
      work_dir: '/tmp/vasp',
      input_files: ['POSCAR', 'INCAR', 'KPOINTS'],
      output_files: [],
      warnings: [],
      metadata: { calculator: 'vasp' },
    });
    const calculator = new ASECalculator(context, {
      type: CalculatorType.VASP,
      command: 'vasp_std',
      parameters: { encut: 520 },
    });

    const result = await calculator.generateInput(atoms, '/tmp/vasp');

    const args = (spawn as jest.Mock).mock.calls[0][1];
    expect(args[1]).toBe('generate');
    expect(args).toHaveLength(3);
    expect(JSON.parse(args[2])).toMatchObject({
      atoms,
      calculator: 'vasp',
      parameters: { encut: 520 },
      work_dir: '/tmp/vasp',
    });
    expect(result.files).toEqual({
      POSCAR: '/tmp/vasp/POSCAR',
      INCAR: '/tmp/vasp/INCAR',
      KPOINTS: '/tmp/vasp/KPOINTS',
    });
  });

  it('runs existing input directories through the run-existing CLI contract', async () => {
    mockSpawnJson({
      success: true,
      status: 'completed',
      work_dir: '/tmp/vasp',
      input_files: ['INCAR'],
      output_files: ['/tmp/vasp/OUTCAR'],
      energy: -10.5,
      warnings: [],
      metadata: { return_code: 0 },
    });
    const calculator = new ASECalculator(context, {
      type: CalculatorType.VASP,
      command: 'vasp_std',
      environment: { OMP_NUM_THREADS: '1' },
      parameters: {},
    });

    const result = await calculator.runCalculation('/tmp/vasp');

    const args = (spawn as jest.Mock).mock.calls[0][1];
    expect(args[1]).toBe('run-existing');
    expect(args[2]).toBe('vasp');
    expect(args[3]).toBe('/tmp/vasp');
    expect(JSON.parse(args[4])).toMatchObject({
      executable: 'vasp_std',
      command: 'vasp_std',
      environment: { OMP_NUM_THREADS: '1' },
    });
    expect(result.energy).toBe(-10.5);
  });

  it('sends read as read <work_dir> <calculator>', async () => {
    mockSpawnJson({
      success: true,
      status: 'completed',
      work_dir: '/tmp/qe',
      input_files: ['qe.in'],
      output_files: ['qe.out'],
      energy: -20,
      warnings: [],
      metadata: { calculator: 'qe' },
    });
    const calculator = new ASECalculator(context, {
      type: CalculatorType.QE,
      command: 'pw.x',
      parameters: {},
    });

    await calculator.readResults('/tmp/qe');

    const args = (spawn as jest.Mock).mock.calls[0][1];
    expect(args.slice(1)).toEqual(['read', '/tmp/qe', 'qe']);
  });

  it('kills the Python calculator process when an abort signal is raised', async () => {
    const child = mockSpawnHanging();
    const calculator = new ASECalculator(context, {
      type: CalculatorType.VASP,
      command: 'vasp_std',
      parameters: {},
    });
    const controller = new AbortController();

    const resultPromise = calculator.runCalculation('/tmp/vasp', undefined, controller.signal);
    controller.abort();
    const result = await resultPromise;

    expect(child.kill).toHaveBeenCalledWith('SIGTERM');
    expect(result.success).toBe(false);
    expect(result.error).toContain('cancelled');
    expect(result.metadata.cancelled).toBe(true);
  });
});

function mockSpawnJson(payload: unknown): void {
  (spawn as jest.Mock).mockImplementation(() => {
    const child = new EventEmitter() as any;
    child.stdout = new EventEmitter();
    child.stderr = new EventEmitter();
    child.kill = jest.fn();
    setImmediate(() => {
      child.stdout.emit('data', JSON.stringify(payload));
      child.emit('close', 0);
    });
    return child;
  });
}

function mockSpawnHanging(): any {
  const child = new EventEmitter() as any;
  child.stdout = new EventEmitter();
  child.stderr = new EventEmitter();
  child.kill = jest.fn(() => {
    setImmediate(() => child.emit('close', null));
    return true;
  });
  (spawn as jest.Mock).mockReturnValue(child);
  return child;
}
