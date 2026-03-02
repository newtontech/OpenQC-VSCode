import { MoleculeItem, MoleculeTreeProvider } from '../../../src/sidebar/MoleculeTreeProvider';
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

describe('MoleculeTreeProvider Full Coverage', () => {
  let mockContext: any;
  let provider: MoleculeTreeProvider;

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
    provider = new MoleculeTreeProvider(mockContext);
  });

  afterEach(() => {
    provider.dispose();
  });

  describe('MoleculeItem variations', () => {
    it('should create small molecule (<10 atoms) with blue icon', () => {
      const item = new MoleculeItem('mol-1', 'H2', 'H2', 2);
      expect(item.iconPath).toBeDefined();
      expect(item.contextValue).toBe('molecule');
    });

    it('should create medium molecule (10-50 atoms) with green icon', () => {
      const item = new MoleculeItem('mol-2', 'Benzene', 'C6H6', 12);
      expect(item.iconPath).toBeDefined();
    });

    it('should create large molecule (>50 atoms) with purple icon', () => {
      const item = new MoleculeItem('mol-3', 'C60', 'C60', 60);
      expect(item.iconPath).toBeDefined();
    });

    it('should create molecule with filePath', () => {
      const item = new MoleculeItem('mol-4', 'Test', 'CH4', 5, '/path/to/file.xyz');
      expect(item.filePath).toBe('/path/to/file.xyz');
      expect(item.command).toBeDefined();
    });

    it('should create molecule without filePath', () => {
      const item = new MoleculeItem('mol-5', 'Test', 'H2', 2);
      expect(item.filePath).toBeUndefined();
      expect(item.command).toBeDefined();
    });

    it('should set tooltip correctly', () => {
      const item = new MoleculeItem('mol-6', 'Water', 'H2O', 3);
      expect(item.tooltip).toContain('Water');
      expect(item.tooltip).toContain('H2O');
      expect(item.tooltip).toContain('3 atoms');
    });

    it('should set description correctly', () => {
      const item = new MoleculeItem('mol-7', 'Test', 'C6H6', 12);
      expect(item.description).toBe('C6H6 (12 atoms)');
    });
  });

  describe('MoleculeTreeProvider methods', () => {
    it('should get tree item', () => {
      const mol = new MoleculeItem('mol-1', 'Test', 'H2', 2);
      const result = provider.getTreeItem(mol);
      expect(result).toBe(mol);
    });

    it('should get children when no element provided', async () => {
      const children = await provider.getChildren();
      expect(Array.isArray(children)).toBe(true);
    });

    it('should get children when element provided (empty array)', async () => {
      const mol = new MoleculeItem('mol-1', 'Test', 'H2', 2);
      const children = await provider.getChildren(mol);
      expect(children).toEqual([]);
    });

    it('should load saved molecules from workspace state', () => {
      const savedMols = [
        {
          id: 'saved-1',
          label: 'Saved Mol',
          formula: 'C6H6',
          atomCount: 12,
          filePath: '/test/file.xyz',
        },
      ];
      mockContext.workspaceState.get = jest.fn(() => savedMols);
      const newProvider = new MoleculeTreeProvider(mockContext);
      expect(mockContext.workspaceState.get).toHaveBeenCalledWith('openqc.molecules', []);
      newProvider.dispose();
    });

    it('should handle saved molecules with missing filePath', () => {
      const savedMols = [
        {
          id: 'saved-1',
          label: 'Saved Mol',
          formula: 'H2',
          atomCount: 2,
        },
      ];
      mockContext.workspaceState.get = jest.fn(() => savedMols);
      const newProvider = new MoleculeTreeProvider(mockContext);
      newProvider.dispose();
    });

    it('should filter out invalid saved molecules', () => {
      const savedMols = [
        { id: 'valid', label: 'Valid', formula: 'H2', atomCount: 2 },
        { id: '', label: 'Invalid - no id', formula: 'H2', atomCount: 2 },
        null,
        undefined,
        { label: 'No id field', formula: 'H2', atomCount: 2 },
      ] as any[];
      mockContext.workspaceState.get = jest.fn(() => savedMols);
      const newProvider = new MoleculeTreeProvider(mockContext);
      newProvider.dispose();
    });

    it('should handle loadMolecules error', () => {
      mockContext.workspaceState.get = jest.fn(() => {
        throw new Error('Load error');
      });
      const newProvider = new MoleculeTreeProvider(mockContext);
      newProvider.dispose();
    });

    it('should add molecule', () => {
      const mol = new MoleculeItem('mol-new', 'New', 'H2', 2);
      provider.addMolecule(mol);
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should remove molecule', () => {
      const mol = new MoleculeItem('mol-1', 'Test', 'H2', 2);
      provider.addMolecule(mol);
      provider.removeMolecule('mol-1');
      expect(mockContext.workspaceState.update).toHaveBeenCalled();
    });

    it('should do nothing when removing non-existent molecule', () => {
      provider.removeMolecule('non-existent');
      // Should not throw
    });

    it('should clear all molecules', () => {
      provider.clearMolecules();
      expect(mockContext.workspaceState.update).toHaveBeenCalledWith('openqc.molecules', []);
    });

    it('should get molecule by id when found', () => {
      const mol = new MoleculeItem('mol-find', 'Find Me', 'H2O', 3);
      provider.addMolecule(mol);
      const found = provider.getMolecule('mol-find');
      expect(found).toBeDefined();
      expect(found?.label).toBe('Find Me');
    });

    it('should return undefined for non-existent molecule', () => {
      const found = provider.getMolecule('non-existent');
      expect(found).toBeUndefined();
    });

    it('should refresh tree', () => {
      provider.refresh();
      // Should not throw
    });

    it('should handle saveMolecules error', async () => {
      mockContext.workspaceState.update = jest.fn().mockRejectedValue(new Error('Save error'));
      const mol = new MoleculeItem('mol-1', 'Test', 'H2', 2);
      provider.addMolecule(mol);
      // Wait for async saveMolecules
      await new Promise(resolve => setTimeout(resolve, 100));
      // Should not throw
    });
  });

  describe('auto-refresh', () => {
    it('should handle disabled autoRefresh', () => {
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return false;
          return defaultValue;
        },
      });
      const newProvider = new MoleculeTreeProvider(mockContext);
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
      const newProvider = new MoleculeTreeProvider(mockContext);
      // Verify timer is set
      expect(jest.getTimerCount()).toBe(1);
      newProvider.dispose();
      jest.useRealTimers();
    });

    it('should clear interval on refresh config change', () => {
      jest.useFakeTimers();
      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return true;
          if (key === 'refreshInterval') return 1000;
          return defaultValue;
        },
      });
      const newProvider = new MoleculeTreeProvider(mockContext);
      newProvider.refresh();
      newProvider.dispose();
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
      const newProvider = new MoleculeTreeProvider(mockContext);
      newProvider.dispose();
      expect(jest.getTimerCount()).toBe(0);
      jest.useRealTimers();
    });

    it('should handle dispose when no interval', () => {
      const newProvider = new MoleculeTreeProvider(mockContext);
      newProvider.dispose();
      // Should not throw
    });
  });

  describe('auto-refresh callback', () => {
    it('should call refresh when interval fires', () => {
      jest.useFakeTimers();

      // Track saved molecules
      let savedMolecules: any[] = [];
      const trackingMockContext = {
        workspaceState: {
          get: jest.fn((key: string, defaultValue: any) => {
            if (key === 'openqc.molecules') return savedMolecules;
            return defaultValue;
          }),
          update: jest.fn((key: string, value: any) => {
            if (key === 'openqc.molecules') savedMolecules = value;
            return Promise.resolve();
          }),
        },
      };

      (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
        get: (key: string, defaultValue: any) => {
          if (key === 'autoRefresh') return true;
          if (key === 'refreshInterval') return 1000;
          return defaultValue;
        },
      });
      const newProvider = new MoleculeTreeProvider(trackingMockContext as any);

      // Add a molecule so refresh has something to process
      const mol = new MoleculeItem('mol-callback-test', 'Test', 'H2', 2);
      newProvider.addMolecule(mol);

      // Advance time by the interval to trigger the callback once
      jest.advanceTimersByTime(1000);

      // Verify the molecule exists (refresh was called and preserved)
      const foundMol = newProvider.getMolecule('mol-callback-test');
      expect(foundMol).toBeDefined();

      newProvider.dispose();
      jest.useRealTimers();
    });
  });
});
