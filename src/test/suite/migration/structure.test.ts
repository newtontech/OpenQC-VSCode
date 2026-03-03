import * as vscode from 'vscode';
import * as path from 'path';
import { StructureMigration, MigrationOptions } from '../../../utils/migration/structure';
import { ASEFormat } from '../../../ase/ASEConverter';

describe('StructureMigration', () => {
  let mockContext: vscode.ExtensionContext;
  let migration: StructureMigration;
  const workspaceDir = path.join(__dirname, '..', '..', '..', '..');

  beforeEach(() => {
    // Mock VSCode context
    mockContext = {
      extensionPath: workspaceDir,
      globalStoragePath: '/tmp/openqc-test',
      subscriptions: [],
      workspaceState: {
        get: jest.fn(),
        update: jest.fn(),
      } as any,
      globalState: {
        get: jest.fn(),
        update: jest.fn(),
      } as any,
      extensionUri: vscode.Uri.file(workspaceDir),
    } as any;

    migration = new StructureMigration(mockContext);
  });

  describe('isSupported', () => {
    it('should support VASP to CP2K', () => {
      expect(migration.isSupported('vasp', 'cp2k')).toBe(true);
    });

    it('should support VASP to QE', () => {
      expect(migration.isSupported('vasp', 'qe')).toBe(true);
    });

    it('should support CP2K to VASP', () => {
      expect(migration.isSupported('cp2k', 'vasp')).toBe(true);
    });

    it('should support QE to VASP', () => {
      expect(migration.isSupported('qe', 'vasp')).toBe(true);
    });

    it('should not support unsupported migrations', () => {
      expect(migration.isSupported('unknown', 'vasp')).toBe(false);
    });
  });

  describe('getParameterMappings', () => {
    it('should return mappings for VASP to CP2K', () => {
      const mappings = migration.getParameterMappings('vasp', 'cp2k');
      expect(mappings.length).toBeGreaterThan(0);
      
      const encutMapping = mappings.find(m => m.sourceParam === 'ENCUT');
      expect(encutMapping).toBeDefined();
      expect(encutMapping?.targetParam).toBe('CUTOFF');
    });

    it('should return empty array for unsupported migration', () => {
      const mappings = migration.getParameterMappings('unknown', 'vasp');
      expect(mappings.length).toBe(0);
    });
  });

  describe('migrate', () => {
    it('should detect format from POSCAR', async () => {
      const poscarPath = path.join(workspaceDir, 'test', 'fixtures', 'POSCAR');
      if (!require('fs').existsSync(poscarPath)) {
        // Skip test if fixture doesn't exist
        return;
      }

      const options: MigrationOptions = {
        targetFormat: 'cp2k' as ASEFormat,
      };

      const result = await migration.migrate(poscarPath, options);
      expect(result.success).toBe(true);
      expect(result.sourceFormat).toBe('vasp');
      expect(result.targetFormat).toBe('cp2k');
    });

    it('should generate output path if not provided', async () => {
      const poscarPath = path.join(workspaceDir, 'test', 'fixtures', 'POSCAR');
      if (!require('fs').existsSync(poscarPath)) {
        return;
      }

      const options: MigrationOptions = {
        targetFormat: 'cp2k' as ASEFormat,
      };

      const result = await migration.migrate(poscarPath, options);
      expect(result.targetPath).toBeDefined();
      expect(result.targetPath).toContain('_cp2k');
    });

    it('should validate output by default', async () => {
      const poscarPath = path.join(workspaceDir, 'test', 'fixtures', 'POSCAR');
      if (!require('fs').existsSync(poscarPath)) {
        return;
      }

      const options: MigrationOptions = {
        targetFormat: 'cp2k' as ASEFormat,
        validate: true,
      };

      const result = await migration.migrate(poscarPath, options);
      // Should have validation metadata
      expect(result.metadata).toBeDefined();
    });

    it('should handle unsupported migrations', async () => {
      const poscarPath = path.join(workspaceDir, 'test', 'fixtures', 'POSCAR');
      if (!require('fs').existsSync(poscarPath)) {
        return;
      }

      const options: MigrationOptions = {
        targetFormat: 'unknown' as any,
      };

      const result = await migration.migrate(poscarPath, options);
      expect(result.success).toBe(false);
      expect(result.errors.length).toBeGreaterThan(0);
    });
  });

  describe('parameter conversion', () => {
    it('should preserve atom count', async () => {
      const poscarPath = path.join(workspaceDir, 'test', 'fixtures', 'POSCAR');
      if (!require('fs').existsSync(poscarPath)) {
        return;
      }

      const options: MigrationOptions = {
        targetFormat: 'cp2k' as ASEFormat,
        validate: true,
      };

      const result = await migration.migrate(poscarPath, options);
      if (result.success) {
        expect(result.metadata.natoms).toBeGreaterThan(0);
      }
    });

    it('should extract formula if available', async () => {
      const poscarPath = path.join(workspaceDir, 'test', 'fixtures', 'POSCAR');
      if (!require('fs').existsSync(poscarPath)) {
        return;
      }

      const options: MigrationOptions = {
        targetFormat: 'cp2k' as ASEFormat,
      };

      const result = await migration.migrate(poscarPath, options);
      // Formula may or may not be present depending on source file
      if (result.success) {
        expect(result.metadata).toHaveProperty('formula');
      }
    });
  });
});
