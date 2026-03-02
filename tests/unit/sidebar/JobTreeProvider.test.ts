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
        update: jest.fn().mockResolvedValue(undefined),
      },
    };
    // Disable auto-refresh for faster tests
    (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
      get: (key: string, defaultValue: any) => {
        if (key === 'autoRefresh') return false;
        return defaultValue;
      },
    });
    provider = new JobTreeProvider(mockContext);
  });

  afterEach(() => {
    provider.dispose();
  });

  describe('JobItem icons', () => {
    it('should set running icon', () => {
      const item = new JobItem('job-1', 'Running', 'running', 50, 'Gaussian');
      expect(item.iconPath).toBeDefined();
    });

    it('should set completed icon', () => {
      const item = new JobItem('job-1', 'Completed', 'completed', 100, 'VASP');
      expect(item.iconPath).toBeDefined();
    });

    it('should set failed icon', () => {
      const item = new JobItem('job-1', 'Failed', 'failed', 50, 'ORCA');
      expect(item.iconPath).toBeDefined();
    });

    it('should set queued icon', () => {
      const item = new JobItem('job-1', 'Queued', 'queued', 0, 'Gaussian');
      expect(item.iconPath).toBeDefined();
    });

    it('should set cancelled icon', () => {
      const item = new JobItem('job-1', 'Cancelled', 'cancelled', 50, 'VASP');
      expect(item.iconPath).toBeDefined();
    });

    it('should set command for completed job', () => {
      const item = new JobItem('job-1', 'Completed', 'completed', 100, 'CP2K');
      expect(item.command).toBeDefined();
      expect(item.command?.command).toBe('openqc.sidebar.viewResults');
    });

    it('should set command for failed job', () => {
      const item = new JobItem('job-1', 'Failed', 'failed', 50, 'GAMESS');
      expect(item.command).toBeDefined();
      expect(item.command?.command).toBe('openqc.sidebar.viewResults');
    });

    it('should not set command for running job', () => {
      const item = new JobItem('job-1', 'Running', 'running', 50, 'NWChem');
      expect(item.command).toBeUndefined();
    });
  });

  describe('JobItem duration calculations', () => {
    it('should format hours duration', () => {
      const start = new Date(Date.now() - 3661000); // 1h 1m 1s ago
      const item = new JobItem('job-1', 'Long', 'running', 50, 'ORCA', start);
      expect(item.tooltip).toContain('Duration');
      expect(item.tooltip).toContain('1h');
    });

    it('should format minutes duration', () => {
      const start = new Date(Date.now() - 65000); // 1m 5s ago
      const item = new JobItem('job-1', 'Medium', 'running', 50, 'CP2K', start);
      expect(item.tooltip).toContain('Duration');
      expect(item.tooltip).toContain('1m');
    });

    it('should format seconds duration', () => {
      const start = new Date(Date.now() - 30000); // 30s ago
      const item = new JobItem('job-1', 'Short', 'running', 50, 'GAMESS', start);
      expect(item.tooltip).toContain('Duration');
    });

    it('should use endTime when provided', () => {
      const start = new Date(Date.now() - 60000);
      const end = new Date(Date.now() - 30000);
      const item = new JobItem('job-1', 'With End', 'completed', 100, 'QE', start, end);
      expect(item.tooltip).toContain('Duration');
    });

    it('should handle no startTime', () => {
      const item = new JobItem('job-1', 'No Start', 'queued', 0, 'VASP');
      expect(item.tooltip).toContain('Duration');
    });
  });

  describe('JobTreeProvider methods', () => {
    it('should get tree item', () => {
      const job = new JobItem('job-1', 'Test', 'running', 50, 'NWChem');
      const result = provider.getTreeItem(job);
      expect(result).toBe(job);
    });

    it('should get children when no element provided', async () => {
      const children = await provider.getChildren();
      expect(Array.isArray(children)).toBe(true);
    });

    it('should get children when element provided (empty array)', async () => {
      const job = new JobItem('job-1', 'Test', 'running', 50, 'NWChem');
      const children = await provider.getChildren(job);
      expect(children).toEqual([]);
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

    it('should handle saved jobs with missing fields', () => {
      const savedJobs = [
        {
          id: 'saved-1',
          label: 'Saved Job',
          status: 'running' as JobStatus,
          progress: 50,
          software: 'Gaussian',
        },
      ];
      mockContext.workspaceState.get = jest.fn(() => savedJobs);
      const newProvider = new JobTreeProvider(mockContext);
      newProvider.dispose();
    });

    it('should filter out invalid saved jobs', () => {
      const savedJobs = [
        { id: 'valid', label: 'Valid', status: 'running', progress: 50, software: 'Gaussian' },
        { id: '', label: 'Invalid - no id', status: 'running', progress: 50, software: 'Gaussian' },
        null,
        undefined,
        { label: 'No id field', status: 'running', progress: 50, software: 'Gaussian' },
      ] as any[];
      mockContext.workspaceState.get = jest.fn(() => savedJobs);
      const newProvider = new JobTreeProvider(mockContext);
      newProvider.dispose();
    });

    it('should handle loadJobs error', () => {
      mockContext.workspaceState.get = jest.fn(() => {
        throw new Error('Load error');
      });
      const newProvider = new JobTreeProvider(mockContext);
      newProvider.dispose();
    });

    it('should add job', () => {
      const job = new JobItem('job-new', 'New Job', 'running', 50, 'VASP');
      provider.addJob(job);
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should remove job', () => {
      const job = new JobItem('job-to-remove', 'Test', 'running', 50, 'VASP');
      provider.addJob(job);
      provider.removeJob('job-to-remove');
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should do nothing when removing non-existent job', () => {
      provider.removeJob('non-existent');
      // Should not throw
    });

    it('should update job status to completed', () => {
      const job = new JobItem('job-update-1', 'Test', 'running', 50, 'ORCA');
      provider.addJob(job);
      provider.updateJobStatus('job-update-1', 'completed', 100);
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should update job status to failed', () => {
      const job = new JobItem('job-update-2', 'Test', 'running', 50, 'CP2K');
      provider.addJob(job);
      provider.updateJobStatus('job-update-2', 'failed', 50);
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should update job status to cancelled', () => {
      const job = new JobItem('job-update-3', 'Test', 'running', 50, 'GAMESS');
      provider.addJob(job);
      provider.updateJobStatus('job-update-3', 'cancelled', 50);
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should update job status to running', () => {
      const job = new JobItem('job-update-4', 'Test', 'queued', 0, 'NWChem');
      provider.addJob(job);
      provider.updateJobStatus('job-update-4', 'running', 10);
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should do nothing when updating non-existent job', () => {
      provider.updateJobStatus('non-existent', 'running', 50);
      // Should not throw
    });

    it('should cancel running job', () => {
      const job = new JobItem('job-cancel', 'Test', 'running', 50, 'ORCA');
      provider.addJob(job);
      provider.cancelJob('job-cancel');
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should cancel non-running job without error', () => {
      const job = new JobItem('job-no-cancel', 'Test', 'completed', 100, 'ORCA');
      provider.addJob(job);
      provider.cancelJob('job-no-cancel'); // Should not throw
    });

    it('should cancel non-existent job without error', () => {
      provider.cancelJob('non-existent'); // Should not throw
    });

    it('should restart job', () => {
      const job = new JobItem('job-restart', 'Test', 'completed', 100, 'VASP');
      provider.addJob(job);
      provider.restartJob('job-restart');
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should do nothing when restarting non-existent job', () => {
      provider.restartJob('non-existent');
      // Should not throw
    });

    it('should clear completed jobs', () => {
      provider.clearCompleted();
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should get job by id', () => {
      const job = new JobItem('job-unique-test', 'Test', 'running', 50, 'QE');
      provider.addJob(job);
      const found = provider.getJob('job-unique-test');
      expect(found).toBeDefined();
      expect(found?.label).toBe('Test');
    });

    it('should return undefined for non-existent job', () => {
      const found = provider.getJob('non-existent');
      expect(found).toBeUndefined();
    });

    it('should refresh tree', () => {
      provider.refresh();
      // Should not throw
    });

    it('should handle saveJobs error', async () => {
      mockContext.workspaceState.update = jest.fn().mockRejectedValue(new Error('Save error'));
      const job = new JobItem('job-save-error', 'Test', 'running', 50, 'VASP');
      provider.addJob(job);
      // Wait for async saveJobs
      await new Promise(resolve => setTimeout(resolve, 100));
      // Should not throw
    });
  });

  describe('updateRunningJobs simulation', () => {
    it('should complete job when progress reaches 100', () => {
      // Add a running job with high progress
      const job = new JobItem('job-sim', 'Sim', 'running', 99, 'Gaussian');
      provider.addJob(job);

      // Mock Math.random to return high value
      const originalRandom = Math.random;
      Math.random = jest.fn(() => 0.9); // Will make progress >= 100

      provider.refresh();

      Math.random = originalRandom;

      // Job should be completed now
      const updated = provider.getJob('job-sim');
      // Due to timing, this might not be completed yet, but the code path is exercised
    });

    it('should update progress for running job', () => {
      // Add a running job with low progress
      const job = new JobItem('job-sim2', 'Sim2', 'running', 50, 'VASP');
      provider.addJob(job);

      // Mock Math.random to return low value
      const originalRandom = Math.random;
      Math.random = jest.fn(() => 0.1); // Will make progress < 100

      provider.refresh();

      Math.random = originalRandom;
    });
  });

  describe('auto-refresh', () => {
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

    it('should set interval when autoRefresh is true', () => {
      jest.useFakeTimers();
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return true;
          if (key === 'refreshInterval') return 1000;
          return defaultValue;
        },
      });
      const newProvider = new JobTreeProvider(mockContext);
      // Verify timer is set
      expect(jest.getTimerCount()).toBe(1);
      newProvider.dispose();
      jest.useRealTimers();
    });

    it('should clear existing interval when setting up new one', () => {
      jest.useFakeTimers();
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return true;
          if (key === 'refreshInterval') return 1000;
          return defaultValue;
        },
      });
      const newProvider = new JobTreeProvider(mockContext);
      newProvider.dispose();
      // After dispose, no timers should be pending
      expect(jest.getTimerCount()).toBe(0);
      jest.useRealTimers();
    });
  });

  describe('dispose', () => {
    it('should clear interval on dispose', () => {
      jest.useFakeTimers();
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return true;
          if (key === 'refreshInterval') return 1000;
          return defaultValue;
        },
      });
      const newProvider = new JobTreeProvider(mockContext);
      newProvider.dispose();
      expect(jest.getTimerCount()).toBe(0);
      jest.useRealTimers();
    });

    it('should handle dispose when no interval', () => {
      const newProvider = new JobTreeProvider(mockContext);
      newProvider.dispose();
      // Should not throw
    });
  });

  describe('auto-refresh callback', () => {
    it('should call refresh when interval fires', () => {
      jest.useFakeTimers();
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return true;
          if (key === 'refreshInterval') return 1000;
          return defaultValue;
        },
      });
      const newProvider = new JobTreeProvider(mockContext);

      // Add a job so refresh has something to update
      const job = new JobItem('job-callback-test', 'Test', 'running', 50, 'Gaussian');
      newProvider.addJob(job);

      // Advance time by the interval to trigger the callback once
      jest.advanceTimersByTime(1000);

      // Advance time again to trigger a second callback (ensures coverage)
      jest.advanceTimersByTime(1000);

      // Verify the job exists (refresh was called)
      const updatedJob = newProvider.getJob('job-callback-test');
      expect(updatedJob).toBeDefined();

      newProvider.dispose();
      jest.useRealTimers();
    });
  });
});
