/**
 * Worker Manager - Manages WebWorker lifecycle and communication
 *
 * Provides a high-level API for offloading heavy computations to background threads.
 */

import * as vscode from 'vscode';
import { ComputeWorker, WorkerMessage, WorkerMessageType } from './computeWorker';

export interface WorkerTask {
  id: string;
  type: WorkerMessageType;
  payload: any;
  priority: 'low' | 'normal' | 'high';
  status: 'pending' | 'running' | 'completed' | 'failed';
  progress?: number;
  result?: any;
  error?: string;
  startTime?: number;
  endTime?: number;
}

export interface WorkerStats {
  activeWorkers: number;
  pendingTasks: number;
  completedTasks: number;
  failedTasks: number;
  averageDuration: number;
}

/**
 * Worker Manager
 *
 * Manages a pool of WebWorkers for parallel computation
 */
export class WorkerManager {
  private workers: Worker[] = [];
  private taskQueue: WorkerTask[] = [];
  private taskMap: Map<string, WorkerTask> = new Map();
  private maxWorkers: number;
  private taskCounter: number = 0;
  private activeTasks: number = 0;
  private completedDurations: number[] = [];
  private readonly computeWorker = new ComputeWorker();
  private stats: WorkerStats = {
    activeWorkers: 0,
    pendingTasks: 0,
    completedTasks: 0,
    failedTasks: 0,
    averageDuration: 0,
  };

  constructor(maxWorkers: number = 4) {
    this.maxWorkers = maxWorkers;
  }

  /**
   * Initialize worker pool
   */
  async initialize(): Promise<void> {
    // In VSCode extension context, we might use a different approach
    // For now, we'll create workers on-demand
    console.log('WorkerManager initialized');
  }

  /**
   * Submit a task to the worker pool
   */
  async submitTask(
    type: WorkerMessageType,
    payload: any,
    priority: 'low' | 'normal' | 'high' = 'normal'
  ): Promise<WorkerTask> {
    const task: WorkerTask = {
      id: this.generateTaskId(),
      type,
      payload,
      priority,
      status: 'pending',
    };

    this.taskMap.set(task.id, task);

    // Add to queue with priority
    this.addToQueue(task);

    // Try to process queue
    this.processQueue();

    return task;
  }

  /**
   * Get task status
   */
  getTask(taskId: string): WorkerTask | undefined {
    return this.taskMap.get(taskId);
  }

  /**
   * Wait for task completion
   */
  async waitForTask(taskId: string, timeout: number = 30000): Promise<WorkerTask> {
    return new Promise((resolve, reject) => {
      const startTime = Date.now();

      const check = () => {
        const task = this.taskMap.get(taskId);

        if (!task) {
          reject(new Error(`Task ${taskId} not found`));
          return;
        }

        if (task.status === 'completed') {
          resolve(task);
          return;
        }

        if (task.status === 'failed') {
          reject(new Error(task.error || 'Task failed'));
          return;
        }

        if (Date.now() - startTime > timeout) {
          reject(new Error('Task timeout'));
          return;
        }

        // Check again in 100ms
        setTimeout(check, 100);
      };

      check();
    });
  }

  /**
   * Cancel a task
   */
  cancelTask(taskId: string): boolean {
    const task = this.taskMap.get(taskId);
    if (!task) {
      return false;
    }

    if (task.status === 'pending') {
      // Remove from queue
      const index = this.taskQueue.indexOf(task);
      if (index >= 0) {
        this.taskQueue.splice(index, 1);
      }
      task.status = 'failed';
      task.error = 'Cancelled by user';
      return true;
    }

    // Cannot cancel running tasks
    return false;
  }

  /**
   * Get worker statistics
   */
  getStats(): WorkerStats {
    return { ...this.stats };
  }

  /**
   * Shutdown all workers
   */
  async shutdown(): Promise<void> {
    // Terminate all workers
    this.workers.forEach(worker => worker.terminate());
    this.workers = [];

    // Clear queues
    this.taskQueue = [];
    this.taskMap.clear();

    console.log('WorkerManager shutdown complete');
  }

  /**
   * Generate unique task ID
   */
  private generateTaskId(): string {
    return `task-${++this.taskCounter}-${Date.now()}`;
  }

  /**
   * Add task to queue with priority
   */
  private addToQueue(task: WorkerTask): void {
    // Priority order: high > normal > low
    const priorityOrder = { high: 0, normal: 1, low: 2 };

    let inserted = false;
    for (let i = 0; i < this.taskQueue.length; i++) {
      if (priorityOrder[task.priority] < priorityOrder[this.taskQueue[i].priority]) {
        this.taskQueue.splice(i, 0, task);
        inserted = true;
        break;
      }
    }

    if (!inserted) {
      this.taskQueue.push(task);
    }

    this.stats.pendingTasks = this.taskQueue.length;
  }

  /**
   * Process task queue
   */
  private processQueue(): void {
    while (this.taskQueue.length > 0 && this.activeTasks < this.maxWorkers) {
      void this.executeTask(this.taskQueue.shift()!);
    }
    this.updateQueueStats();
  }

  /**
   * Execute a task
   */
  private async executeTask(task: WorkerTask): Promise<void> {
    task.status = 'running';
    task.startTime = Date.now();
    this.activeTasks++;
    this.updateQueueStats();

    try {
      const response = await this.computeWorker.processMessage({
        id: task.id,
        type: task.type,
        payload: task.payload,
      });
      if (!response.success) {
        throw new Error(response.error || 'Worker task failed');
      }

      task.result = response.result;
      task.status = 'completed';
      task.endTime = Date.now();

      this.stats.completedTasks++;
    } catch (error) {
      task.error = error instanceof Error ? error.message : String(error);
      task.status = 'failed';
      task.endTime = Date.now();

      this.stats.failedTasks++;
    } finally {
      this.activeTasks = Math.max(0, this.activeTasks - 1);
      if (task.startTime && task.endTime) {
        this.completedDurations.push(task.endTime - task.startTime);
        this.stats.averageDuration =
          this.completedDurations.reduce((sum, duration) => sum + duration, 0) /
          this.completedDurations.length;
      }
      this.updateQueueStats();
    }

    // Process next task
    this.processQueue();
  }

  private updateQueueStats(): void {
    this.stats.activeWorkers = this.activeTasks;
    this.stats.pendingTasks = this.taskQueue.length;
  }
}

// Singleton instance
let workerManager: WorkerManager | null = null;

/**
 * Get or create WorkerManager instance
 */
export function getWorkerManager(): WorkerManager {
  if (!workerManager) {
    workerManager = new WorkerManager();
  }
  return workerManager;
}
