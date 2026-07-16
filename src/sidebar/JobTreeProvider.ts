import * as vscode from 'vscode';
import { createComponentLogger } from '../utils/Logger';

/**
 * Job status types
 */
export type JobStatus = 'running' | 'completed' | 'failed' | 'queued' | 'cancelled';

export interface JobResultData {
  success: boolean;
  energy?: number;
  forces?: number[][];
  stress?: number[][];
  error?: string;
  warnings?: string[];
  metadata?: Record<string, any>;
  outputFiles?: string[];
}

interface StoredJob {
  id: string;
  label: string;
  status: JobStatus;
  progress: number;
  software: string;
  startTime?: string;
  endTime?: string;
  workingDirectory?: string;
  result?: JobResultData;
}

/**
 * Represents a calculation job item in the tree view
 */
export class JobItem extends vscode.TreeItem {
  constructor(
    public readonly id: string,
    public readonly label: string,
    public readonly status: JobStatus,
    public readonly progress: number,
    public readonly software: string,
    public readonly startTime?: Date,
    public readonly endTime?: Date,
    public readonly workingDirectory?: string,
    public readonly result?: JobResultData,
    public readonly collapsibleState: vscode.TreeItemCollapsibleState = vscode
      .TreeItemCollapsibleState.None
  ) {
    super(label, collapsibleState);

    const duration = this.calculateDuration();
    this.tooltip = `${label}\nStatus: ${status}\nSoftware: ${software}\nProgress: ${progress}%\nDuration: ${duration}${workingDirectory ? `\nDirectory: ${workingDirectory}` : ''}`;
    this.description = `${software} - ${progress}%`;
    this.contextValue = status;

    // Set icon based on status
    switch (status) {
      case 'running':
        this.iconPath = new vscode.ThemeIcon('sync~spin', new vscode.ThemeColor('charts.yellow'));
        break;
      case 'completed':
        this.iconPath = new vscode.ThemeIcon('check', new vscode.ThemeColor('charts.green'));
        break;
      case 'failed':
        this.iconPath = new vscode.ThemeIcon('error', new vscode.ThemeColor('charts.red'));
        break;
      case 'queued':
        this.iconPath = new vscode.ThemeIcon('clock', new vscode.ThemeColor('charts.blue'));
        break;
      case 'cancelled':
        this.iconPath = new vscode.ThemeIcon('circle-slash', new vscode.ThemeColor('charts.gray'));
        break;
    }

    // Add command for completed/failed jobs to view results
    if (status === 'completed' || status === 'failed') {
      this.command = {
        command: 'openqc.sidebar.viewResults',
        title: 'View Results',
        arguments: [this],
      };
    }
  }

  /**
   * Calculate job duration
   */
  private calculateDuration(): string {
    const end = this.endTime || new Date();
    const start = this.startTime || end;
    const diff = end.getTime() - start.getTime();
    const seconds = Math.floor(diff / 1000);
    const minutes = Math.floor(seconds / 60);
    const hours = Math.floor(minutes / 60);

    if (hours > 0) {
      return `${hours}h ${minutes % 60}m`;
    } else if (minutes > 0) {
      return `${minutes}m ${seconds % 60}s`;
    } else {
      return `${seconds}s`;
    }
  }
}

/**
 * Tree data provider for the Calculation Jobs view
 */
export class JobTreeProvider implements vscode.TreeDataProvider<JobItem> {
  private _onDidChangeTreeData: vscode.EventEmitter<JobItem | undefined | null | void> =
    new vscode.EventEmitter<JobItem | undefined | null | void>();
  readonly onDidChangeTreeData: vscode.Event<JobItem | undefined | null | void> =
    this._onDidChangeTreeData.event;

  private jobs: JobItem[] = [];
  private autoRefreshInterval: ReturnType<typeof setInterval> | undefined;
  private logger = createComponentLogger('JobTreeProvider');

  constructor(private context: vscode.ExtensionContext) {
    this.loadJobs();
    this.setupAutoRefresh();
  }

  /**
   * Get a tree item for an element
   */
  getTreeItem(element: JobItem): vscode.TreeItem {
    return element;
  }

  /**
   * Get children of an element (root items if no element provided)
   */
  getChildren(element?: JobItem): Thenable<JobItem[]> {
    if (element) {
      // No children for job items (flat list)
      return Promise.resolve([]);
    }
    return Promise.resolve(this.jobs);
  }

  /**
   * Refresh the tree view
   */
  refresh(): void {
    this.updateRunningJobs();
    this._onDidChangeTreeData.fire();
  }

  /**
   * Add a new job to the view
   */
  addJob(job: JobItem): void {
    this.jobs.push(job);
    this.saveJobs();
    this._onDidChangeTreeData.fire();
  }

  /**
   * Update job status
   */
  updateJobStatus(id: string, status: JobStatus, progress: number): boolean {
    const index = this.jobs.findIndex(j => j.id === id);
    if (index >= 0) {
      const job = this.jobs[index];
      // Create updated job with new status
      const updatedJob = new JobItem(
        job.id,
        job.label,
        status,
        progress,
        job.software,
        job.startTime,
        status === 'completed' || status === 'failed' || status === 'cancelled'
          ? new Date()
          : undefined,
        job.workingDirectory,
        job.result
      );
      this.jobs[index] = updatedJob;
      this.saveJobs();
      this._onDidChangeTreeData.fire();
      return true;
    }
    return false;
  }

  updateJobResult(
    id: string,
    status: 'completed' | 'failed',
    result: JobResultData,
    progress = status === 'completed' ? 100 : 0
  ): boolean {
    const index = this.jobs.findIndex(j => j.id === id);
    if (index >= 0) {
      const job = this.jobs[index];
      if (job.status === 'cancelled') {
        return false;
      }
      this.jobs[index] = new JobItem(
        job.id,
        job.label,
        status,
        progress,
        job.software,
        job.startTime,
        new Date(),
        job.workingDirectory,
        result
      );
      this.saveJobs();
      this._onDidChangeTreeData.fire();
      return true;
    }
    return false;
  }

  /**
   * Remove a job from the view
   */
  removeJob(id: string): void {
    const index = this.jobs.findIndex(j => j.id === id);
    if (index >= 0) {
      this.jobs.splice(index, 1);
      this.saveJobs();
      this._onDidChangeTreeData.fire();
    }
  }

  /**
   * Cancel a running job
   */
  cancelJob(id: string): boolean {
    const job = this.jobs.find(j => j.id === id);
    if (job && (job.status === 'running' || job.status === 'queued')) {
      return this.updateJobStatus(id, 'cancelled', job.progress);
    }
    return false;
  }

  /**
   * Restart a job
   */
  restartJob(id: string): JobItem | undefined {
    const job = this.jobs.find(j => j.id === id);
    if (job) {
      const newJob = new JobItem(
        `${job.id}-restart-${Date.now()}`,
        `${job.label} (restart)`,
        'running',
        5,
        job.software,
        new Date(),
        undefined,
        job.workingDirectory
      );
      this.addJob(newJob);
      return newJob;
    }
    return undefined;
  }

  /**
   * Clear completed/failed jobs
   */
  clearCompleted(): void {
    this.jobs = this.jobs.filter(j => j.status === 'running' || j.status === 'queued');
    this.saveJobs();
    this._onDidChangeTreeData.fire();
  }

  /**
   * Get job by ID
   */
  getJob(id: string): JobItem | undefined {
    return this.jobs.find(j => j.id === id);
  }

  /**
   * Load jobs from workspace state
   */
  private loadJobs(): void {
    try {
      const saved = this.context.workspaceState.get<StoredJob[]>('openqc.jobs', []);
      this.jobs = saved
        .filter(j => j && j.id)
        .map(
          j =>
            new JobItem(
              j.id,
              j.label,
              j.status,
              j.progress,
              j.software,
              j.startTime ? new Date(j.startTime) : undefined,
              j.endTime ? new Date(j.endTime) : undefined,
              j.workingDirectory,
              j.result
            )
        );
    } catch (error) {
      this.logger.error('Failed to load jobs from workspace state', error as Error);
      this.jobs = [];
    }
  }

  /**
   * Save jobs to workspace state
   */
  private async saveJobs(): Promise<void> {
    try {
      const stored: StoredJob[] = this.jobs.map(job => ({
        id: job.id,
        label: job.label,
        status: job.status,
        progress: job.progress,
        software: job.software,
        startTime: job.startTime?.toISOString(),
        endTime: job.endTime?.toISOString(),
        workingDirectory: job.workingDirectory,
        result: job.result,
      }));

      await this.context.workspaceState.update('openqc.jobs', stored);
    } catch (error) {
      this.logger.error('Failed to save jobs to workspace state', error as Error);
    }
  }

  private updateRunningJobs(): void {
    // Jobs are advanced by concrete command handlers when calculator execution returns.
    // Refresh only redraws persisted state.
  }

  /**
   * Setup automatic refresh
   */
  private setupAutoRefresh(): void {
    const config = vscode.workspace.getConfiguration('openqc.sidebar');
    const autoRefresh = config.get<boolean>('autoRefresh', true);
    const interval = config.get<number>('refreshInterval', 5000);

    if (this.autoRefreshInterval) {
      clearInterval(this.autoRefreshInterval);
      this.autoRefreshInterval = undefined;
    }

    if (autoRefresh) {
      this.autoRefreshInterval = setInterval(() => {
        this.refresh();
      }, interval);
    }
  }

  /**
   * Dispose of resources
   */
  dispose(): void {
    if (this.autoRefreshInterval) {
      clearInterval(this.autoRefreshInterval);
      this.autoRefreshInterval = undefined;
    }
  }
}
