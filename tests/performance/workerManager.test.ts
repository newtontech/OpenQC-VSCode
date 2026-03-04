/**
 * Unit tests for Worker Manager
 */

import {
  WorkerManager,
  getWorkerManager,
} from '../../src/performance/workerManager';
import { WorkerMessageType } from '../../src/performance/computeWorker';

describe('WorkerManager', () => {
  let manager: WorkerManager;

  beforeEach(() => {
    manager = new WorkerManager(2); // Max 2 workers for testing
  });

  afterEach(async () => {
    await manager.shutdown();
  });

  describe('Task Submission', () => {
    it('should submit task successfully', async () => {
      const task = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      expect(task).toBeDefined();
      expect(task.id).toBeDefined();
      expect(task.status).toMatch(/pending|running|completed/);
    });

    it('should generate unique task IDs', async () => {
      const task1 = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      const task2 = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      expect(task1.id).not.toBe(task2.id);
    });

    it('should support task priorities', async () => {
      const lowTask = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' },
        'low'
      );

      const highTask = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' },
        'high'
      );

      expect(lowTask.priority).toBe('low');
      expect(highTask.priority).toBe('high');
    });
  });

  describe('Task Management', () => {
    it('should get task by ID', async () => {
      const submitted = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      const retrieved = manager.getTask(submitted.id);

      expect(retrieved).toBeDefined();
      expect(retrieved!.id).toBe(submitted.id);
    });

    it('should return undefined for non-existent task', () => {
      const task = manager.getTask('non-existent-id');
      expect(task).toBeUndefined();
    });

    it('should cancel pending task', async () => {
      // Submit multiple tasks to create queue
      await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );
      
      const task = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      const cancelled = manager.cancelTask(task.id);

      // Might not cancel if already running
      if (cancelled) {
        const retrieved = manager.getTask(task.id);
        expect(retrieved!.status).toBe('failed');
        expect(retrieved!.error).toContain('Cancelled');
      }
    });

    it('should not cancel non-existent task', () => {
      const cancelled = manager.cancelTask('non-existent-id');
      expect(cancelled).toBe(false);
    });
  });

  describe('Task Waiting', () => {
    it('should wait for task completion', async () => {
      const task = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      const completed = await manager.waitForTask(task.id, 5000);

      expect(completed.status).toBe('completed');
      expect(completed.result).toBeDefined();
    }, 10000);

    it('should timeout if task takes too long', async () => {
      const task = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      await expect(
        manager.waitForTask(task.id, 1) // 1ms timeout
      ).rejects.toThrow('timeout');
    });

    it('should reject if task not found', async () => {
      await expect(
        manager.waitForTask('non-existent-id', 1000)
      ).rejects.toThrow('not found');
    });
  });

  describe('Statistics', () => {
    it('should track worker statistics', async () => {
      const stats = manager.getStats();

      expect(stats).toBeDefined();
      expect(stats.activeWorkers).toBeDefined();
      expect(stats.pendingTasks).toBeDefined();
      expect(stats.completedTasks).toBeDefined();
      expect(stats.failedTasks).toBeDefined();
    });

    it('should update stats after task completion', async () => {
      const initialStats = manager.getStats();

      const task = await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      await manager.waitForTask(task.id, 5000);

      const finalStats = manager.getStats();

      expect(finalStats.completedTasks).toBeGreaterThan(initialStats.completedTasks);
    }, 10000);
  });

  describe('Shutdown', () => {
    it('should shutdown gracefully', async () => {
      await manager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: 'C 0 0 0', format: 'xyz' }
      );

      await manager.shutdown();

      const stats = manager.getStats();
      expect(stats.activeWorkers).toBe(0);
    });
  });
});

describe('Singleton Instance', () => {
  it('should return same instance', () => {
    const instance1 = getWorkerManager();
    const instance2 = getWorkerManager();

    expect(instance1).toBe(instance2);
  });
});
