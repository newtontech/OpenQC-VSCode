const mockExecPythonJson = jest.fn();

jest.mock('../../../src/python/PythonBridge', () => ({
  execPythonJson: mockExecPythonJson,
}));

jest.mock('../../../src/utils/Logger', () => ({
  Logger: {
    getInstance: () => ({
      info: jest.fn(),
      warn: jest.fn(),
      debug: jest.fn(),
      error: jest.fn(),
    }),
  },
}));

import {
  extractTrajectory,
  parseOutput,
  summarizeOutput,
} from '../../../src/results/OutputParserBridge';

describe('OutputParserBridge', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    mockExecPythonJson.mockResolvedValue({ success: true, data: {} });
  });

  it('parses calculation output through the Python output bridge', async () => {
    await parseOutput('/tmp/gaussian.log', 'gaussian');

    expect(mockExecPythonJson).toHaveBeenCalledWith('openqc.bridge.output_bridge', {
      command: 'parse',
      args: { path: '/tmp/gaussian.log', software: 'gaussian' },
    });
  });

  it('summarizes calculation output through the Python output bridge', async () => {
    await summarizeOutput('/tmp/orca.out');

    expect(mockExecPythonJson).toHaveBeenCalledWith('openqc.bridge.output_bridge', {
      command: 'summarize',
      args: { path: '/tmp/orca.out', software: 'auto' },
    });
  });

  it('extracts optimization trajectories through the Python output bridge', async () => {
    await extractTrajectory('/tmp/gaussian.log', 'gaussian');

    expect(mockExecPythonJson).toHaveBeenCalledWith('openqc.bridge.output_bridge', {
      command: 'trajectory',
      args: { path: '/tmp/gaussian.log', software: 'gaussian' },
    });
  });
});
