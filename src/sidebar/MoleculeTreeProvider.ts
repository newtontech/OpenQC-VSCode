import * as vscode from 'vscode';

/**
 * Represents a molecule item in the tree view
 */
export class MoleculeItem extends vscode.TreeItem {
  constructor(
    public readonly id: string,
    public readonly label: string,
    public readonly formula: string,
    public readonly atomCount: number,
    public readonly filePath?: string,
    public readonly collapsibleState: vscode.TreeItemCollapsibleState = vscode
      .TreeItemCollapsibleState.None
  ) {
    super(label, collapsibleState);

    this.tooltip = `${label} (${formula}) - ${atomCount} atoms`;
    this.description = `${formula} (${atomCount} atoms)`;
    this.contextValue = 'molecule';

    // Set icon based on molecule size
    if (atomCount < 10) {
      this.iconPath = new vscode.ThemeIcon('symbol-color', new vscode.ThemeColor('charts.blue'));
    } else if (atomCount < 50) {
      this.iconPath = new vscode.ThemeIcon('symbol-color', new vscode.ThemeColor('charts.green'));
    } else {
      this.iconPath = new vscode.ThemeIcon('symbol-color', new vscode.ThemeColor('charts.purple'));
    }

    // Command to open the molecule when clicked
    this.command = {
      command: 'openqc.sidebar.openMolecule',
      title: 'Open Molecule',
      arguments: [this],
    };
  }
}

/**
 * Tree data provider for the Molecules view
 */
export class MoleculeTreeProvider implements vscode.TreeDataProvider<MoleculeItem> {
  private _onDidChangeTreeData: vscode.EventEmitter<MoleculeItem | undefined | null | void> =
    new vscode.EventEmitter<MoleculeItem | undefined | null | void>();
  readonly onDidChangeTreeData: vscode.Event<MoleculeItem | undefined | null | void> =
    this._onDidChangeTreeData.event;

  private molecules: MoleculeItem[] = [];
  private autoRefreshInterval: ReturnType<typeof setInterval> | undefined;

  constructor(private context: vscode.ExtensionContext) {
    this.loadMolecules();
    this.setupAutoRefresh();
  }

  /**
   * Get a tree item for an element
   */
  getTreeItem(element: MoleculeItem): vscode.TreeItem {
    return element;
  }

  /**
   * Get children of an element (root items if no element provided)
   */
  getChildren(element?: MoleculeItem): Thenable<MoleculeItem[]> {
    if (element) {
      // No children for molecule items (flat list)
      return Promise.resolve([]);
    }
    return Promise.resolve(this.molecules);
  }

  /**
   * Refresh the tree view
   */
  refresh(): void {
    this.loadMolecules();
    this._onDidChangeTreeData.fire();
  }

  /**
   * Add a new molecule to the view
   */
  addMolecule(molecule: MoleculeItem): void {
    this.molecules.push(molecule);
    this.saveMolecules();
    this._onDidChangeTreeData.fire();
  }

  /**
   * Remove a molecule from the view
   */
  removeMolecule(id: string): void {
    const index = this.molecules.findIndex(m => m.id === id);
    if (index >= 0) {
      this.molecules.splice(index, 1);
      this.saveMolecules();
      this._onDidChangeTreeData.fire();
    }
  }

  /**
   * Clear all molecules
   */
  clearMolecules(): void {
    this.molecules = [];
    this.saveMolecules();
    this._onDidChangeTreeData.fire();
  }

  /**
   * Get molecule by ID
   */
  getMolecule(id: string): MoleculeItem | undefined {
    return this.molecules.find(m => m.id === id);
  }

  /**
   * Load molecules from workspace state
   */
  private loadMolecules(): void {
    try {
      const saved = this.context.workspaceState.get<MoleculeItem[]>('openqc.molecules', []);
      this.molecules = saved
        .filter(m => m && m.id)
        .map(m => new MoleculeItem(m.id, m.label, m.formula, m.atomCount, m.filePath));

      // If no saved molecules, add some sample data for demonstration
      if (this.molecules.length === 0) {
        this.addSampleMolecules();
      }
    } catch (error) {
      console.error('Failed to load molecules:', error);
      this.molecules = [];
      // Add sample molecules as fallback
      this.addSampleMolecules();
    }
  }

  /**
   * Save molecules to workspace state
   */
  private async saveMolecules(): Promise<void> {
    try {
      await this.context.workspaceState.update('openqc.molecules', this.molecules);
    } catch (error) {
      console.error('Failed to save molecules:', error);
    }
  }

  /**
   * Add sample molecules for demonstration
   */
  private addSampleMolecules(): void {
    const samples = [
      new MoleculeItem('mol-1', 'Water', 'H2O', 3, undefined),
      new MoleculeItem('mol-2', 'Benzene', 'C6H6', 12, undefined),
      new MoleculeItem('mol-3', 'Caffeine', 'C8H10N4O2', 24, undefined),
      new MoleculeItem('mol-4', 'Aspirin', 'C9H8O4', 21, undefined),
      new MoleculeItem('mol-5', 'Buckminsterfullerene', 'C60', 60, undefined),
    ];
    this.molecules = samples;
    this.saveMolecules();
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
