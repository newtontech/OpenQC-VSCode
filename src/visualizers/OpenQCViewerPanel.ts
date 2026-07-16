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
import { OpenQCViewerWebview } from '../webviews/openqcViewerWebview';
import { Logger } from '../utils/Logger';
import { StructureExporter, type StructureData } from '../utils/structureExporter';

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
  private _dirty = false;

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
      OpenQCViewerPanel.currentPanel._replaceStructure(structure, filename);
      return;
    }

    const panel = vscode.window.createWebviewPanel(
      OpenQCViewerPanel.viewType,
      `OpenQC: ${filename}`,
      column || vscode.ViewColumn.One,
      OpenQCViewerWebview.getWebviewOptions(extensionUri)
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
    private _filename: string
  ) {
    this._panel = panel;

    // Validate before sending
    const validation = validateOpenQCStructure(structure);
    if (!validation.valid) {
      const msg = `Invalid OpenQCStructure: ${validation.errors.join('; ')}`;
      logger.warn(msg);
      getOutputChannel().appendLine(msg);
    }

    // Render webview HTML using bundled assets (no CDN dependency)
    this._panel.webview.html = OpenQCViewerWebview.generateHTML(
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
      message => this._handleMessage(message),
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

  /** Return whether the viewer has unsaved webview-side edits. */
  public isDirty(): boolean {
    return this._dirty;
  }

  public static getCurrentStructureData(): StructureData | undefined {
    const structure = OpenQCViewerPanel.currentPanel?._currentStructure;
    if (!structure) {
      return undefined;
    }

    return {
      atoms: structure.atoms.map(atom => ({
        element: atom.element,
        x: atom.x,
        y: atom.y,
        z: atom.z,
      })),
      bonds: structure.bonds,
      cell: structure.cell,
      pbc: structure.cell?.pbc,
      metadata: structure.metadata,
    };
  }

  public static async saveCurrentStructureToSource(): Promise<boolean> {
    if (!OpenQCViewerPanel.currentPanel) {
      vscode.window.showErrorMessage('No OpenQC structure viewer is open');
      return false;
    }

    return OpenQCViewerPanel.currentPanel._saveCurrentStructureToSource();
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
    this._dirty = false;
    this._updateTitle();

    // Validate
    const validation = validateOpenQCStructure(structure);

    // Send the full DTO as JSON so the webview can render it
    this._panel.webview.postMessage({
      type: 'loadStructure',
      structure: JSON.stringify(structure),
      filename: this._filename,
      validation,
    });

    logger.info(`Viewer: loaded ${structure.atoms.length} atoms from ${this._filename}`);
    if (validation.warnings.length > 0) {
      getOutputChannel().appendLine(`Warnings: ${validation.warnings.join('; ')}`);
    }
  }

  private _replaceStructure(structure: OpenQCStructure, filename: string): void {
    this._filename = filename;
    this._sendStructure(structure);
  }

  // -----------------------------------------------------------------------
  // Message handlers
  // -----------------------------------------------------------------------

  private _handleMessage(message: any): void {
    switch (message.type) {
      case 'exportImage':
        OpenQCViewerPanel._handleExportImage(message.data);
        break;
      case 'structureEdited':
      case 'structureUpdated':
        this._handleStructureUpdate(message.structure, message.dirty ?? true);
        break;
      case 'exportEditedStructure':
      case 'exportStructure':
        this._handleStructureExportRequest(message.structure, message.dirty);
        break;
      case 'saveEditedStructureToSource':
        if (
          message.structure &&
          !this._handleStructureUpdate(message.structure, message.dirty ?? true)
        ) {
          break;
        }
        void this._saveCurrentStructureToSource();
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

  private _handleStructureUpdate(structurePayload: unknown, dirty: boolean): boolean {
    const structure = parseStructurePayload(structurePayload);
    if (!structure) {
      vscode.window.showErrorMessage('OpenQC Viewer: edited structure message was invalid');
      return false;
    }

    const validation = validateOpenQCStructure(structure);
    if (!validation.valid) {
      vscode.window.showErrorMessage(
        `OpenQC Viewer: edited structure is invalid: ${validation.errors.join('; ')}`
      );
      return false;
    }

    this._currentStructure = structure;
    this._dirty = dirty;
    this._updateTitle();
    if (validation.warnings.length > 0) {
      getOutputChannel().appendLine(`Edited structure warnings: ${validation.warnings.join('; ')}`);
    }
    return true;
  }

  private _handleStructureExportRequest(structurePayload: unknown, dirty?: unknown): void {
    const nextDirty = typeof dirty === 'boolean' ? dirty : this._dirty;
    if (structurePayload && !this._handleStructureUpdate(structurePayload, nextDirty)) {
      return;
    }
    void vscode.commands.executeCommand('openqc.exportStructureWithPicker');
  }

  private _updateTitle(): void {
    this._panel.title = `OpenQC: ${this._dirty ? '• ' : ''}${this._filename}`;
  }

  private async _saveCurrentStructureToSource(): Promise<boolean> {
    const structureData = OpenQCViewerPanel.getCurrentStructureData();
    if (!this._currentStructure || !structureData || structureData.atoms.length === 0) {
      vscode.window.showErrorMessage('No edited structure is available to save');
      return false;
    }

    const format = StructureExporter.inferNativeFormatFromPath(this._filename);
    if (!format) {
      vscode.window.showErrorMessage(
        'OpenQC can only write edited structures back to XYZ, Extended XYZ, PDB, CIF, POSCAR, or CONTCAR source files. Use Export Structure for other formats.'
      );
      return false;
    }

    if (!this._dirty) {
      vscode.window.showInformationMessage('No OpenQC viewer edits to save');
      return false;
    }

    const confirmation = await vscode.window.showWarningMessage(
      `Overwrite source file with the edited structure?\n${this._filename}`,
      { modal: true },
      'Save'
    );
    if (confirmation !== 'Save') {
      return false;
    }

    const exporter = new StructureExporter();
    const result = await exporter.overwriteStructureFile(structureData, this._filename);
    if (!result.success) {
      vscode.window.showErrorMessage(`Save failed: ${result.error}`);
      return false;
    }

    this._dirty = false;
    this._updateTitle();
    await this._panel.webview.postMessage({
      type: 'markStructureSaved',
      structure: JSON.stringify(this._currentStructure),
    });
    vscode.window.showInformationMessage(`Edited structure saved to ${this._filename}`);
    return true;
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

function parseStructurePayload(payload: unknown): OpenQCStructure | undefined {
  if (!payload) {
    return undefined;
  }

  if (typeof payload === 'string') {
    try {
      return JSON.parse(payload) as OpenQCStructure;
    } catch {
      return undefined;
    }
  }

  if (typeof payload === 'object') {
    return payload as OpenQCStructure;
  }

  return undefined;
}
