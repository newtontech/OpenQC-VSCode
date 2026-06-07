/**
 * Worker Manager - Manages WebWorker lifecycle and communication
 *
 * Provides a high-level API for offloading heavy computations to background threads.
 */

import * as vscode from 'vscode';
import { WorkerMessage, WorkerResponse, WorkerMessageType } from './computeWorker';

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
    // Check if we can start more tasks
    if (this.taskQueue.length === 0) {
      return;
    }

    // For now, execute tasks synchronously
    // In production, this would use actual WebWorkers
    this.executeTask(this.taskQueue.shift()!);
  }

  /**
   * Execute a task
   */
  private async executeTask(task: WorkerTask): Promise<void> {
    task.status = 'running';
    task.startTime = Date.now();

    try {
      // Simulate worker execution
      // In production, this would send message to WebWorker
      const result = await this.simulateWorkerExecution(task);

      task.result = result;
      task.status = 'completed';
      task.endTime = Date.now();

      this.stats.completedTasks++;
    } catch (error) {
      task.error = error instanceof Error ? error.message : String(error);
      task.status = 'failed';
      task.endTime = Date.now();

      this.stats.failedTasks++;
    }

    this.stats.pendingTasks = this.taskQueue.length;

    // Process next task
    this.processQueue();
  }

  /**
   * Simulate worker execution (for testing)
   */
  private async simulateWorkerExecution(task: WorkerTask): Promise<any> {
    // Simulate computation time
    await new Promise(resolve => setTimeout(resolve, 100));

    // Return mock result
    switch (task.type) {
      case WorkerMessageType.PARSE_STRUCTURE:
        return {
          chemical_symbols: ['C', 'C', 'C'],
          positions: [
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
          ],
          pbc: [false, false, false],
        };

      case WorkerMessageType.CONVERT_FORMAT:
        return '# Converted structure\nC 0 0 0\nC 1 0 0\n';

      case WorkerMessageType.VALIDATE_STRUCTURE:
        return { valid: true, warnings: [] };

      case WorkerMessageType.CALCULATE_PROPERTIES:
        return { center_of_mass: [1, 0, 0] };

      case WorkerMessageType.MIGRATE_PARAMETERS:
        return { migrated: {}, warnings: [] };

      default:
        throw new Error(`Unknown task type: ${task.type}`);
    }
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
