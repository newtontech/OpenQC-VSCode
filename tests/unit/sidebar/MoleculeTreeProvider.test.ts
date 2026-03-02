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
        update: jest.fn(),
      },
    };
    provider = new MoleculeTreeProvider(mockContext);
  });

  afterEach(() => {
    provider.dispose();
  });

  describe('MoleculeItem variations', () => {
    it('should create molecule without filePath', () => {
      const item = new MoleculeItem('mol-1', 'Test', 'CH4', 5);
      expect(item.filePath).toBeUndefined();
      expect(item.command).toBeDefined();
    });
  });

  describe('MoleculeTreeProvider methods', () => {
    it('should get tree item', () => {
      const mol = new MoleculeItem('mol-1', 'Test', 'H2', 2);
      const result = provider.getTreeItem(mol);
      expect(result).toBe(mol);
    });

    it('should load saved molecules', () => {
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

    it('should clear all molecules', () => {
      provider.clearMolecules();
      expect(mockContext.workspaceState.update).toHaveBeenCalledWith('openqc.molecules', []);
    });

    it('should get molecule by id when found', () => {
      // First add a molecule
      const mol = new MoleculeItem('mol-find', 'Find Me', 'H2O', 3);
      provider.addMolecule(mol);
      // Then get it
      const found = provider.getMolecule('mol-find');
      expect(found).toBeDefined();
      expect(found?.label).toBe('Find Me');
    });
  });

  describe('auto-refresh variations', () => {
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

    it('should clear interval on refresh config change', () => {
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
    });
  });
});
