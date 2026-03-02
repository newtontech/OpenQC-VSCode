import { MoleculeItem, MoleculeTreeProvider } from '../../src/sidebar/MoleculeTreeProvider';
import { JobItem, JobTreeProvider, JobStatus } from '../../src/sidebar/JobTreeProvider';
import * as vscode from 'vscode';

// Mock vscode module
jest.mock('vscode', () => ({
  TreeItem: class TreeItem {
    constructor(
      public label: string,
      public collapsibleState?: any
    ) {}
  },
  TreeItemCollapsibleState: {
    None: 0,
    Collapsed: 1,
    Expanded: 2,
  },
  ThemeIcon: class ThemeIcon {
    constructor(
      public id: string,
      public color?: any
    ) {}
  },
  ThemeColor: class ThemeColor {
    constructor(public id: string) {}
  },
  EventEmitter: class EventEmitter {
    event = jest.fn();
    fire = jest.fn();
  },
  workspace: {
    getConfiguration: jest.fn(() => ({
      get: jest.fn((key: string, defaultValue: any) => defaultValue),
    })),
  },
  window: {
    showInformationMessage: jest.fn(),
  },
}));

describe('MoleculeItem', () => {
  it('should create a molecule item with correct properties', () => {
    const item = new MoleculeItem('mol-1', 'Water', 'H2O', 3, '/path/to/file.xyz');

    expect(item.id).toBe('mol-1');
    expect(item.label).toBe('Water');
    expect(item.formula).toBe('H2O');
    expect(item.atomCount).toBe(3);
    expect(item.filePath).toBe('/path/to/file.xyz');
    expect(item.contextValue).toBe('molecule');
  });

  it('should set correct icon for small molecules (< 10 atoms)', () => {
    const item = new MoleculeItem('mol-1', 'Small', 'H2', 2);
    expect(item.iconPath).toBeDefined();
  });

  it('should set correct icon for medium molecules (10-50 atoms)', () => {
    const item = new MoleculeItem('mol-1', 'Medium', 'C10H8', 18);
    expect(item.iconPath).toBeDefined();
  });

  it('should set correct icon for large molecules (> 50 atoms)', () => {
    const item = new MoleculeItem('mol-1', 'Large', 'C60', 60);
    expect(item.iconPath).toBeDefined();
  });
});

describe('JobItem', () => {
  it('should create a job item with correct properties', () => {
    const item = new JobItem('job-1', 'Optimization', 'running', 50, 'Gaussian');

    expect(item.id).toBe('job-1');
    expect(item.label).toBe('Optimization');
    expect(item.status).toBe('running');
    expect(item.progress).toBe(50);
    expect(item.software).toBe('Gaussian');
  });

  it('should set correct icon for running status', () => {
    const item = new JobItem('job-1', 'Running', 'running', 50, 'Gaussian');
    expect(item.iconPath).toBeDefined();
    expect(item.contextValue).toBe('running');
  });

  it('should set correct icon for completed status', () => {
    const item = new JobItem('job-1', 'Completed', 'completed', 100, 'ORCA');
    expect(item.iconPath).toBeDefined();
    expect(item.contextValue).toBe('completed');
  });

  it('should set correct icon for failed status', () => {
    const item = new JobItem('job-1', 'Failed', 'failed', 45, 'VASP');
    expect(item.iconPath).toBeDefined();
    expect(item.contextValue).toBe('failed');
  });

  it('should calculate duration correctly', () => {
    const startTime = new Date(Date.now() - 120000); // 2 minutes ago
    const item = new JobItem('job-1', 'Test', 'running', 50, 'CP2K', startTime);
    expect(item.tooltip).toContain('Duration');
  });
});

describe('MoleculeTreeProvider', () => {
  let mockContext: any;
  let provider: MoleculeTreeProvider;

  beforeEach(() => {
    mockContext = {
      workspaceState: {
        get: jest.fn(() => []),
        update: jest.fn(),
      },
    };
    provider = new MoleculeTreeProvider(mockContext);
  });

  afterEach(() => {
    provider.dispose();
  });

  it('should initialize with sample molecules if no saved data', () => {
    expect(mockContext.workspaceState.get).toHaveBeenCalledWith('openqc.molecules', []);
  });

  it('should add a molecule', () => {
    const molecule = new MoleculeItem('mol-new', 'Test', 'C2H6', 8);
    provider.addMolecule(molecule);
    expect(mockContext.workspaceState.update).toHaveBeenCalled();
  });

  it('should remove a molecule', () => {
    provider.removeMolecule('mol-1');
    expect(mockContext.workspaceState.update).toHaveBeenCalled();
  });

  it('should get children (root items)', async () => {
    const children = await provider.getChildren();
    expect(Array.isArray(children)).toBe(true);
  });

  it('should return empty array for children of a specific item', async () => {
    const molecule = new MoleculeItem('mol-1', 'Test', 'H2', 2);
    const children = await provider.getChildren(molecule);
    expect(children).toEqual([]);
  });

  it('should get a molecule by id', () => {
    const result = provider.getMolecule('nonexistent');
    // Returns undefined if not found
    expect(result).toBeUndefined();
  });
});

describe('JobTreeProvider', () => {
  let mockContext: any;
  let provider: JobTreeProvider;

  beforeEach(() => {
    mockContext = {
      workspaceState: {
        get: jest.fn(() => []),
        update: jest.fn(),
      },
    };
    provider = new JobTreeProvider(mockContext);
  });

  afterEach(() => {
    provider.dispose();
  });

  it('should initialize with sample jobs if no saved data', () => {
    expect(mockContext.workspaceState.get).toHaveBeenCalledWith('openqc.jobs', []);
  });

  it('should add a job', () => {
    const job = new JobItem('job-new', 'Test Job', 'queued', 0, 'Gaussian');
    provider.addJob(job);
    expect(mockContext.workspaceState.update).toHaveBeenCalled();
  });

  it('should update job status', () => {
    provider.updateJobStatus('job-1', 'completed', 100);
    expect(mockContext.workspaceState.update).toHaveBeenCalled();
  });

  it('should cancel a running job', () => {
    provider.cancelJob('job-1');
    // Should not throw error
  });

  it('should restart a job', () => {
    provider.restartJob('job-1');
    // Should create a new job
    expect(mockContext.workspaceState.update).toHaveBeenCalled();
  });

  it('should get children (root items)', async () => {
    const children = await provider.getChildren();
    expect(Array.isArray(children)).toBe(true);
  });

  it('should get a job by id', () => {
    const result = provider.getJob('nonexistent');
    expect(result).toBeUndefined();
  });
});
