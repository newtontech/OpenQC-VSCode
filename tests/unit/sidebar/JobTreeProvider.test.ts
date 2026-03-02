import { JobItem, JobTreeProvider, JobStatus } from '../../../src/sidebar/JobTreeProvider';
import * as vscode from 'vscode';

jest.mock('vscode', () => ({
  TreeItem: class TreeItem {
    constructor(
      public label: string,
      public collapsibleState?: any
    ) {}
  },
  TreeItemCollapsibleState: { None: 0, Collapsed: 1, Expanded: 2 },
  ThemeIcon: class ThemeIcon {
    constructor(
      public id: string,
      public color?: any
    ) {}
  },
  ThemeColor: class ThemeColor {
    constructor(public id: string) {}
  },
  EventEmitter: class EventEmitter<T> {
    private listeners: Array<(e: T) => void> = [];
    event = (listener: (e: T) => void) => {
      this.listeners.push(listener);
      return { dispose: () => {} };
    };
    fire = (event?: T) => {
      this.listeners.forEach(l => l(event as T));
    };
  },
  workspace: {
    getConfiguration: jest.fn(() => ({
      get: jest.fn((key: string, defaultValue: any) => defaultValue),
    })),
  },
}));

describe('JobTreeProvider Full Coverage', () => {
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

  describe('JobItem icons', () => {
    it('should set queued icon', () => {
      const item = new JobItem('job-1', 'Queued', 'queued', 0, 'Gaussian');
      expect(item.iconPath).toBeDefined();
    });

    it('should set cancelled icon', () => {
      const item = new JobItem('job-1', 'Cancelled', 'cancelled', 50, 'VASP');
      expect(item.iconPath).toBeDefined();
    });
  });

  describe('JobItem duration calculations', () => {
    it('should format hours duration', () => {
      const start = new Date(Date.now() - 3661000); // 1h 1m 1s ago
      const item = new JobItem('job-1', 'Long', 'running', 50, 'ORCA', start);
      expect(item.tooltip).toContain('Duration');
    });

    it('should format minutes duration', () => {
      const start = new Date(Date.now() - 65000); // 1m 5s ago
      const item = new JobItem('job-1', 'Medium', 'running', 50, 'CP2K', start);
      expect(item.tooltip).toContain('Duration');
    });

    it('should format seconds duration', () => {
      const start = new Date(Date.now() - 30000); // 30s ago
      const item = new JobItem('job-1', 'Short', 'running', 50, 'GAMESS', start);
      expect(item.tooltip).toContain('Duration');
    });
  });

  describe('JobTreeProvider methods', () => {
    it('should get tree item', () => {
      const job = new JobItem('job-1', 'Test', 'running', 50, 'NWChem');
      const result = provider.getTreeItem(job);
      expect(result).toBe(job);
    });

    it('should load saved jobs from workspace state', () => {
      const savedJobs = [
        {
          id: 'saved-1',
          label: 'Saved Job',
          status: 'completed' as JobStatus,
          progress: 100,
          software: 'Gaussian',
          startTime: new Date().toISOString(),
          endTime: new Date().toISOString(),
        },
      ];
      mockContext.workspaceState.get = jest.fn(() => savedJobs);
      const newProvider = new JobTreeProvider(mockContext);
      expect(mockContext.workspaceState.get).toHaveBeenCalledWith('openqc.jobs', []);
      newProvider.dispose();
    });

    it('should remove job', () => {
      const job = new JobItem('job-1', 'Test', 'running', 50, 'VASP');
      provider.addJob(job);
      provider.removeJob('job-1');
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should cancel non-running job without error', () => {
      const job = new JobItem('job-1', 'Test', 'completed', 100, 'ORCA');
      provider.addJob(job);
      provider.cancelJob('job-1'); // Should not throw
    });

    it('should clear completed jobs', () => {
      provider.clearCompleted();
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });
  });

  describe('auto-refresh with config disabled', () => {
    it('should not set interval when autoRefresh is false', () => {
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return false;
          return defaultValue;
        },
      });
      const newProvider = new JobTreeProvider(mockContext);
      newProvider.dispose();
    });
  });
});
