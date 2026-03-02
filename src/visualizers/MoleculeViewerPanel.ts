import * as vscode from 'vscode';
import { MoleculeViewerWebview } from './MoleculeViewerWebview';

export interface MoleculeViewerState {
  xyz?: string;
  json?: string;
  filename: string;
}

export class MoleculeViewerPanel {
  public static currentPanel: MoleculeViewerPanel | undefined;
  public static readonly viewType = 'openqc.moleculeViewer';

  private readonly _panel: vscode.WebviewPanel;
  private _disposables: vscode.Disposable[] = [];

  /**
   * Create or show the molecule viewer panel
   */
  public static createOrShow(extensionUri: vscode.Uri, xyzContent: string, filename: string): void {
    const column = vscode.window.activeTextEditor
      ? vscode.window.activeTextEditor.viewColumn
      : undefined;

    // If we already have a panel, show it
    if (MoleculeViewerPanel.currentPanel) {
      MoleculeViewerPanel.currentPanel._panel.reveal(column);
      MoleculeViewerPanel.currentPanel._sendInitializeMessage({ xyz: xyzContent });
      return;
    }

    // Otherwise, create a new panel
    const panel = vscode.window.createWebviewPanel(
      MoleculeViewerPanel.viewType,
      '3D Molecule Viewer',
      column || vscode.ViewColumn.One,
      MoleculeViewerWebview.getWebviewOptions()
    );

    MoleculeViewerPanel.currentPanel = new MoleculeViewerPanel(
      panel,
      extensionUri,
      xyzContent,
      filename
    );
  }

  /**
   * Handle incoming messages from the webview
   */
  public static handleMessage(message: any): void {
    switch (message.type) {
      case 'exportImage':
        MoleculeViewerPanel._handleExportImage(message.data);
        break;
      case 'error':
        vscode.window.showErrorMessage(`Molecule Viewer: ${message.message}`);
        break;
      case 'info':
        vscode.window.showInformationMessage(message.message);
        break;
    }
  }

  private constructor(
    panel: vscode.WebviewPanel,
    private readonly _extensionUri: vscode.Uri,
    xyzContent: string,
    private readonly _filename: string
  ) {
    this._panel = panel;

    // Set the webview's initial HTML content
    this._panel.webview.html = MoleculeViewerWebview.generateWebviewHTML();

    // Listen for when the panel is disposed
    this._panel.onDidDispose(() => this.dispose(), null, this._disposables);

    // Update the content when the panel becomes visible
    this._panel.onDidChangeViewState(
      () => {
        if (this._panel.visible) {
          this._sendInitializeMessage({ xyz: xyzContent });
        }
      },
      null,
      this._disposables
    );

    // Handle messages from the webview
    this._panel.webview.onDidReceiveMessage(
      message => MoleculeViewerPanel.handleMessage(message),
      null,
      this._disposables
    );

    // Send initial structure data
    this._sendInitializeMessage({ xyz: xyzContent });
  }

  private _sendInitializeMessage(structure: { xyz?: string; json?: string }): void {
    this._panel.webview.postMessage({
      type: 'initialize',
      structure,
      filename: this._filename,
    });
  }

  private static async _handleExportImage(blob: Blob): Promise<void> {
    const uri = await vscode.window.showSaveDialog({
      filters: {
        Images: ['png'],
      },
      defaultUri: vscode.Uri.file('molecule.png'),
    });

    if (!uri) {
      return;
    }

    try {
      const arrayBuffer = await blob.arrayBuffer();
      const buffer = Buffer.from(arrayBuffer);
      await vscode.workspace.fs.writeFile(uri, buffer);
      vscode.window.showInformationMessage(`Image saved to ${uri.fsPath}`);
    } catch (error) {
      vscode.window.showErrorMessage(`Failed to save image: ${error}`);
    }
  }

  public dispose(): void {
    MoleculeViewerPanel.currentPanel = undefined;

    this._panel.dispose();

    while (this._disposables.length) {
      const disposable = this._disposables.pop();
      if (disposable) {
        disposable.dispose();
      }
    }
  }

  /**
   * Serialize the panel state for restoration
   */
  public serialize(): MoleculeViewerState {
    return {
      filename: this._filename,
    };
  }

  /**
   * Deserialize and restore panel state
   */
  public static deserialize(state: MoleculeViewerState, extensionUri: vscode.Uri): void {
    // Note: In a real implementation, we'd need to store the XYZ content
    // This is a simplified version that shows the pattern
    MoleculeViewerPanel.createOrShow(extensionUri, state.xyz || '', state.filename);
  }
}
