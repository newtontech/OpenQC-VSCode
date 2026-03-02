import * as vscode from 'vscode';

/**
 * Job status types
 */
export type JobStatus = 'running' | 'completed' | 'failed' | 'queued' | 'cancelled';

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
    public readonly collapsibleState: vscode.TreeItemCollapsibleState = vscode
      .TreeItemCollapsibleState.None
  ) {
    super(label, collapsibleState);

    const duration = this.calculateDuration();
    this.tooltip = `${label}\nStatus: ${status}\nSoftware: ${software}\nProgress: ${progress}%\nDuration: ${duration}`;
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
  updateJobStatus(id: string, status: JobStatus, progress: number): void {
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
          : undefined
      );
      this.jobs[index] = updatedJob;
      this.saveJobs();
      this._onDidChangeTreeData.fire();
    }
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
  cancelJob(id: string): void {
    const job = this.jobs.find(j => j.id === id);
    if (job && job.status === 'running') {
      this.updateJobStatus(id, 'cancelled', job.progress);
    }
  }

  /**
   * Restart a job
   */
  restartJob(id: string): void {
    const job = this.jobs.find(j => j.id === id);
    if (job) {
      const newJob = new JobItem(
        `${job.id}-restart-${Date.now()}`,
        `${job.label} (restart)`,
        'queued',
        0,
        job.software,
        undefined
      );
      this.addJob(newJob);
    }
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
      const saved = this.context.workspaceState.get<any[]>('openqc.jobs', []);
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
              j.endTime ? new Date(j.endTime) : undefined
            )
        );

      // If no saved jobs, add some sample data for demonstration
      if (this.jobs.length === 0) {
        this.addSampleJobs();
      }
    } catch (error) {
      console.error('Failed to load jobs:', error);
      this.jobs = [];
      // Add sample jobs as fallback
      this.addSampleJobs();
    }
  }

  /**
   * Save jobs to workspace state
   */
  private async saveJobs(): Promise<void> {
    try {
      await this.context.workspaceState.update('openqc.jobs', this.jobs);
    } catch (error) {
      console.error('Failed to save jobs:', error);
    }
  }

  /**
   * Add sample jobs for demonstration
   */
  private addSampleJobs(): void {
    const samples = [
      new JobItem(
        'job-1',
        'Water Optimization',
        'running',
        65,
        'Gaussian',
        new Date(Date.now() - 120000)
      ),
      new JobItem(
        'job-2',
        'Benzene Frequency',
        'completed',
        100,
        'ORCA',
        new Date(Date.now() - 3600000),
        new Date(Date.now() - 3000000)
      ),
      new JobItem('job-3', 'Surface Calculation', 'queued', 0, 'VASP', undefined),
      new JobItem(
        'job-4',
        'MD Simulation',
        'failed',
        45,
        'CP2K',
        new Date(Date.now() - 7200000),
        new Date(Date.now() - 7000000)
      ),
    ];
    this.jobs = samples;
    this.saveJobs();
  }

  /**
   * Update progress for running jobs (simulation)
   */
  private updateRunningJobs(): void {
    this.jobs.forEach(job => {
      if (job.status === 'running' && job.progress < 100) {
        const increment = Math.random() * 5;
        const newProgress = Math.min(job.progress + increment, 100);
        if (newProgress >= 100) {
          this.updateJobStatus(job.id, 'completed', 100);
        } else {
          this.updateJobStatus(job.id, 'running', Math.floor(newProgress));
        }
      }
    });
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
