/**
 * OpenQCViewerPanel - DTO-driven viewer panel for molecular/crystal structures.
 *
 * Accepts `OpenQCStructure` as its primary input contract and forwards the
 * serialized DTO to the webview for rendering.  Falls back to the legacy
 * XYZ path when needed.
 *
 * Replaces `MoleculeViewerPanel` as the canonical entry-point for the
 * `openqc.visualizeStructure` command.
 *
 * @module visualizers/OpenQCViewerPanel
 */

import * as vscode from 'vscode';
import type { OpenQCStructure } from '../structures/OpenQCStructure';
import { validateOpenQCStructure, type ValidationResult } from '../structures/validation';
import { openQCStructureToXYZ } from '../structures/converters';
import { MoleculeViewerWebview } from './MoleculeViewerWebview';
import { Logger } from '../utils/Logger';

const logger = Logger.getInstance();

// ---------------------------------------------------------------------------
// Output channel (lazy-created)
// ---------------------------------------------------------------------------

let _outputChannel: vscode.OutputChannel | undefined;

function getOutputChannel(): vscode.OutputChannel {
  if (!_outputChannel) {
    _outputChannel = vscode.window.createOutputChannel('OpenQC Viewer');
  }
  return _outputChannel;
}

// ---------------------------------------------------------------------------
// Viewer state
// ---------------------------------------------------------------------------

export interface OpenQCViewerState {
  structure?: OpenQCStructure;
  filename: string;
}

// ---------------------------------------------------------------------------
// Panel
// ---------------------------------------------------------------------------

export class OpenQCViewerPanel {
  public static currentPanel: OpenQCViewerPanel | undefined;
  public static readonly viewType = 'openqc.openqcViewer';

  private readonly _panel: vscode.WebviewPanel;
  private _disposables: vscode.Disposable[] = [];
  private _currentStructure?: OpenQCStructure;

  // -----------------------------------------------------------------------
  // Factory
  // -----------------------------------------------------------------------

  /**
   * Create or reveal the DTO-driven viewer.
   *
   * @param extensionUri - Extension root URI.
   * @param structure    - The canonical structure DTO.
   * @param filename     - Human-readable filename for the panel title.
   */
  public static createOrShow(
    extensionUri: vscode.Uri,
    structure: OpenQCStructure,
    filename: string
  ): void {
    const column = vscode.window.activeTextEditor
      ? vscode.window.activeTextEditor.viewColumn
      : undefined;

    if (OpenQCViewerPanel.currentPanel) {
      OpenQCViewerPanel.currentPanel._panel.reveal(column);
      OpenQCViewerPanel.currentPanel._sendStructure(structure);
      return;
    }

    const panel = vscode.window.createWebviewPanel(
      OpenQCViewerPanel.viewType,
      `OpenQC: ${filename}`,
      column || vscode.ViewColumn.One,
      MoleculeViewerWebview.getWebviewOptions(extensionUri)
    );

    OpenQCViewerPanel.currentPanel = new OpenQCViewerPanel(
      panel,
      extensionUri,
      structure,
      filename
    );
  }

  // -----------------------------------------------------------------------
  // Instance
  // -----------------------------------------------------------------------

  private constructor(
    panel: vscode.WebviewPanel,
    private readonly _extensionUri: vscode.Uri,
    structure: OpenQCStructure,
    private readonly _filename: string
  ) {
    this._panel = panel;

    // Validate before sending
    const validation = validateOpenQCStructure(structure);
    if (!validation.valid) {
      const msg = `Invalid OpenQCStructure: ${validation.errors.join('; ')}`;
      logger.warn(msg);
      getOutputChannel().appendLine(msg);
    }

    // Render webview HTML
    this._panel.webview.html = MoleculeViewerWebview.generateWebviewHTML(
      this._panel.webview,
      this._extensionUri
    );

    this._panel.onDidDispose(() => this.dispose(), null, this._disposables);

    this._panel.onDidChangeViewState(
      () => {
        if (this._panel.visible && this._currentStructure) {
          this._sendStructure(this._currentStructure);
        }
      },
      null,
      this._disposables
    );

    this._panel.webview.onDidReceiveMessage(
      message => OpenQCViewerPanel._handleMessage(message),
      null,
      this._disposables
    );

    // Send initial structure
    this._sendStructure(structure);
  }

  // -----------------------------------------------------------------------
  // Public helpers
  // -----------------------------------------------------------------------

  /** Validate and return diagnostics for a proposed structure. */
  public static validate(structure: unknown): ValidationResult {
    return validateOpenQCStructure(structure);
  }

  /** Return the current structure (used for testing / state serialization). */
  public getCurrentStructure(): OpenQCStructure | undefined {
    return this._currentStructure;
  }

  public dispose(): void {
    OpenQCViewerPanel.currentPanel = undefined;
    this._panel.dispose();
    while (this._disposables.length) {
      const d = this._disposables.pop();
      if (d) {
        d.dispose();
      }
    }
  }

  // -----------------------------------------------------------------------
  // Internal messaging
  // -----------------------------------------------------------------------

  private _sendStructure(structure: OpenQCStructure): void {
    this._currentStructure = structure;

    // Validate
    const validation = validateOpenQCStructure(structure);

    // Send the full DTO as JSON so the webview can render it
    this._panel.webview.postMessage({
      type: 'loadStructure',
      structure: JSON.stringify(structure),
      filename: this._filename,
      validation,
    });

    // Also send as XYZ for backward compatibility with NGL-based webview
    try {
      const xyz = openQCStructureToXYZ(structure);
      this._panel.webview.postMessage({
        type: 'initialize',
        structure: { xyz },
        filename: this._filename,
      });
    } catch {
      // If XYZ conversion fails, the DTO path is still available
    }

    logger.info(`Viewer: loaded ${structure.atoms.length} atoms from ${this._filename}`);
    if (validation.warnings.length > 0) {
      getOutputChannel().appendLine(`Warnings: ${validation.warnings.join('; ')}`);
    }
  }

  // -----------------------------------------------------------------------
  // Message handlers
  // -----------------------------------------------------------------------

  private static _handleMessage(message: any): void {
    switch (message.type) {
      case 'exportImage':
        OpenQCViewerPanel._handleExportImage(message.data);
        break;
      case 'error':
        vscode.window.showErrorMessage(`OpenQC Viewer: ${message.message}`);
        break;
      case 'info':
        vscode.window.showInformationMessage(message.message);
        break;
      case 'viewerReady':
        logger.debug('Viewer webview reports ready');
        break;
    }
  }

  private static async _handleExportImage(data: any): Promise<void> {
    // Accept both Blob-like objects and base64 data URLs
    let buffer: Buffer;

    if (typeof data === 'string' && data.startsWith('data:')) {
      // base64 data URL
      const base64 = data.split(',')[1];
      buffer = Buffer.from(base64, 'base64');
    } else if (typeof data === 'string') {
      // raw base64
      buffer = Buffer.from(data, 'base64');
    } else {
      vscode.window.showErrorMessage('Image export received unsupported data format');
      return;
    }

    const uri = await vscode.window.showSaveDialog({
      filters: { Images: ['png'] },
      defaultUri: vscode.Uri.file('openqc-snapshot.png'),
    });

    if (!uri) {
      return;
    }

    try {
      await vscode.workspace.fs.writeFile(uri, buffer);
      vscode.window.showInformationMessage(`Image saved to ${uri.fsPath}`);
    } catch (error) {
      vscode.window.showErrorMessage(`Failed to save image: ${error}`);
    }
  }
}
