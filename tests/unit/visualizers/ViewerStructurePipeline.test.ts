import { ViewerStructurePipeline } from '../../../src/visualizers/ViewerStructurePipeline';
import type { WorkerManager } from '../../../src/performance/workerManager';

jest.mock('vscode', () => require('../../mocks/vscode'));

function createContext(): any {
  return {
    globalStorageUri: { fsPath: '/tmp/openqc-viewer-pipeline-tests' },
    subscriptions: [],
  };
}

function createDocument(fileName: string, content: string, version = 1): any {
  return {
    fileName,
    version,
    lineCount: content.split(/\r?\n/).length,
    uri: {
      fsPath: fileName,
      toString: () => `file://${fileName}`,
    },
    getText: () => content,
  };
}

describe('ViewerStructurePipeline', () => {
  it('routes POSCAR parsing through the real worker and then serves the viewer cache', async () => {
    const submitTask = jest.fn().mockResolvedValue({ id: 'viewer-worker-task' });
    const waitForTask = jest.fn().mockResolvedValue({
      status: 'completed',
      result: {
        chemical_symbols: ['Si'],
        positions: [[0, 0, 0]],
        pbc: [true, true, true],
        cell: [
          [5, 0, 0],
          [0, 5, 0],
          [0, 0, 5],
        ],
      },
    });
    const worker = { submitTask, waitForTask } as unknown as WorkerManager;
    const pipeline = new ViewerStructurePipeline(createContext(), worker);
    const document = createDocument('/tmp/openqc-viewer-worker/POSCAR', 'fixture');

    const first = await pipeline.parse(document, 'VASP');
    const second = await pipeline.parse(document, 'VASP');

    expect(first.source).toBe('worker');
    expect(first.structure.atoms).toEqual([{ element: 'Si', x: 0, y: 0, z: 0 }]);
    expect(first.structure.metadata?.source?.parser).toBe('compute-worker-native');
    expect(second.source).toBe('cache');
    expect(submitTask).toHaveBeenCalledTimes(1);
    expect(waitForTask).toHaveBeenCalledWith('viewer-worker-task', 10_000);
  });

  it('uses the incremental parser for supported non-worker documents', async () => {
    const worker = {
      submitTask: jest.fn(),
      waitForTask: jest.fn(),
    } as unknown as WorkerManager;
    const pipeline = new ViewerStructurePipeline(createContext(), worker);
    const prefix = Array.from({ length: 120 }, () => '# padding').join('\n');
    const firstContent = `${prefix}\n&COORD\nH 0 0 0\n&END COORD`;
    const secondContent = `${prefix}\n&COORD\nH 1 0 0\n&END COORD`;

    const first = await pipeline.parse(createDocument('/tmp/viewer.cp2k', firstContent, 1), 'CP2K');
    const second = await pipeline.parse(
      createDocument('/tmp/viewer.cp2k', secondContent, 2),
      'CP2K'
    );

    expect(first.source).toBe('incremental');
    expect(first.incremental).toBe(false);
    expect(second.source).toBe('incremental');
    expect(second.incremental).toBe(true);
    expect(second.structure.atoms[0].x).toBe(1);
    expect(worker.submitTask).not.toHaveBeenCalled();
  });
});
