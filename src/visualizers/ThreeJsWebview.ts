/**
 * VSCode Webview integration for Three.js molecular visualization
 *
 * Phase 2: 3D Visualization - Webview Integration
 *
 * This module provides the webview HTML/JS for Three.js rendering
 * within the VSCode extension context.
 */

import * as vscode from 'vscode';
import { ThreeJsRenderer } from './ThreeJsRenderer';
import { Atom, RepresentationMode } from './types';

export interface ThreeJsWebviewOptions {
  title: string;
  onViewChange?: (state: any) => void;
}

export class ThreeJsWebview {
  private static currentPanel: ThreeJsWebview | undefined;
  private readonly panel: vscode.WebviewPanel;
  private disposables: vscode.Disposable[] = [];
  private renderer: ThreeJsRenderer | null = null;

  /**
   * Create or show the webview panel
   */
  public static createOrShow(
    extensionUri: vscode.Uri,
    options: ThreeJsWebviewOptions
  ): ThreeJsWebview {
    const column = vscode.window.activeTextEditor
      ? vscode.window.activeTextEditor.viewColumn
      : undefined;

    if (ThreeJsWebview.currentPanel) {
      ThreeJsWebview.currentPanel.panel.reveal(column);
      return ThreeJsWebview.currentPanel;
    }

    const panel = vscode.window.createWebviewPanel(
      'openqc.threeJsViewer',
      options.title,
      column || vscode.ViewColumn.One,
      {
        enableScripts: true,
        retainContextWhenHidden: true,
        localResourceRoots: [vscode.Uri.joinPath(extensionUri, 'media')],
      }
    );

    ThreeJsWebview.currentPanel = new ThreeJsWebview(panel, extensionUri, options);
    return ThreeJsWebview.currentPanel;
  }

  private constructor(
    panel: vscode.WebviewPanel,
    private readonly extensionUri: vscode.Uri,
    private readonly options: ThreeJsWebviewOptions
  ) {
    this.panel = panel;

    // Set the webview's initial HTML content
    this.panel.webview.html = this.getWebviewContent();

    // Handle messages from the webview
    this.panel.webview.onDidReceiveMessage(
      (message) => this.handleMessage(message),
      null,
      this.disposables
    );

    // Clean up when panel is disposed
    this.panel.onDidDispose(() => this.dispose(), null, this.disposables);
  }

  /**
   * Get webview HTML content
   */
  private getWebviewContent(): string {
    const nonce = getNonce();

    return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'nonce-${nonce}'; style-src 'unsafe-inline';">
    <title>3D Molecular Visualization</title>
    <style>
        body {
            margin: 0;
            padding: 0;
            overflow: hidden;
            font-family: var(--vscode-font-family);
            background-color: var(--vscode-editor-background);
            color: var(--vscode-editor-foreground);
        }
        #container {
            width: 100%;
            height: 100vh;
            position: relative;
        }
        #three-canvas {
            width: 100%;
            height: 100%;
            display: block;
        }
        .controls {
            position: absolute;
            top: 10px;
            right: 10px;
            background: var(--vscode-editor-background);
            border: 1px solid var(--vscode-panel-border);
            border-radius: 4px;
            padding: 10px;
            display: flex;
            flex-direction: column;
            gap: 8px;
            z-index: 1000;
        }
        .control-group {
            display: flex;
            flex-direction: column;
            gap: 4px;
        }
        .control-label {
            font-size: 12px;
            color: var(--vscode-descriptionForeground);
            font-weight: 600;
        }
        .control-row {
            display: flex;
            gap: 8px;
        }
        button {
            background: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            padding: 6px 12px;
            border-radius: 2px;
            cursor: pointer;
            font-size: 13px;
            min-width: 80px;
        }
        button:hover {
            background: var(--vscode-button-hoverBackground);
        }
        select {
            background: var(--vscode-dropdown-background);
            color: var(--vscode-dropdown-foreground);
            border: 1px solid var(--vscode-dropdown-border);
            padding: 4px 8px;
            border-radius: 2px;
            font-size: 13px;
            min-width: 120px;
        }
        .info-panel {
            position: absolute;
            bottom: 10px;
            left: 10px;
            background: var(--vscode-editor-background);
            border: 1px solid var(--vscode-panel-border);
            border-radius: 4px;
            padding: 8px 12px;
            font-size: 13px;
            max-width: 300px;
        }
        .loading {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            text-align: center;
        }
        .spinner {
            border: 3px solid var(--vscode-panel-border);
            border-top: 3px solid var(--vscode-button-background);
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
            margin: 0 auto 10px;
        }
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
    </style>
</head>
<body>
    <div id="container">
        <div id="loading" class="loading" style="display: none;">
            <div class="spinner"></div>
            <div>Loading structure...</div>
        </div>
        <canvas id="three-canvas"></canvas>

        <div class="controls">
            <div class="control-group">
                <div class="control-label">Representation</div>
                <select id="representation-mode">
                    <option value="ball-and-stick">Ball and Stick</option>
                    <option value="space-filling">Space Filling</option>
                    <option value="wireframe">Wireframe</option>
                    <option value="stick">Stick</option>
                </select>
            </div>

            <div class="control-group">
                <div class="control-label">View</div>
                <div class="control-row">
                    <button id="reset-view">Reset</button>
                    <button id="toggle-axes">Axes</button>
                </div>
            </div>

            <div class="control-group">
                <div class="control-label">Export</div>
                <button id="export-image">Save Image</button>
            </div>
        </div>

        <div class="info-panel" id="info-panel" style="display: none;">
            <div id="atom-count">Atoms: 0</div>
            <div id="bond-count">Bonds: 0</div>
        </div>
    </div>

    <script nonce="${nonce}">
        // Three.js will be loaded here via postMessage communication
        const vscode = acquireVsCodeApi();

        // DOM elements
        const canvas = document.getElementById('three-canvas');
        const loading = document.getElementById('loading');
        const infoPanel = document.getElementById('info-panel');
        const atomCountEl = document.getElementById('atom-count');
        const bondCountEl = document.getElementById('bond-count');

        // Controls
        const representationMode = document.getElementById('representation-mode');
        const resetViewBtn = document.getElementById('reset-view');
        const toggleAxesBtn = document.getElementById('toggle-axes');
        const exportImageBtn = document.getElementById('export-image');

        // State
        let currentStructure = null;
        let currentConfig = {
            representationMode: 'ball-and-stick',
            showAxes: false
        };

        // Event handlers
        representationMode.addEventListener('change', () => {
            currentConfig.representationMode = representationMode.value;
            vscode.postMessage({
                type: 'changeRepresentation',
                mode: representationMode.value
            });
        });

        resetViewBtn.addEventListener('click', () => {
            vscode.postMessage({ type: 'resetView' });
        });

        toggleAxesBtn.addEventListener('click', () => {
            currentConfig.showAxes = !currentConfig.showAxes;
            vscode.postMessage({
                type: 'toggleAxes',
                show: currentConfig.showAxes
            });
        });

        exportImageBtn.addEventListener('click', () => {
            vscode.postMessage({ type: 'exportImage', format: 'png' });
        });

        // Handle messages from extension
        window.addEventListener('message', event => {
            const message = event.data;

            switch (message.type) {
                case 'initialize':
                    loading.style.display = 'block';
                    currentStructure = message.structure;
                    updateInfoPanel(message.structure);
                    vscode.postMessage({
                        type: 'render',
                        structure: message.structure,
                        config: currentConfig
                    });
                    break;

                case 'renderComplete':
                    loading.style.display = 'none';
                    infoPanel.style.display = 'block';
                    break;

                case 'exportComplete':
                    vscode.postMessage({
                        type: 'saveImage',
                        data: message.data
                    });
                    break;

                case 'error':
                    loading.style.display = 'none';
                    showError(message.message);
                    break;
            }
        });

        function updateInfoPanel(structure) {
            if (structure && structure.atoms) {
                atomCountEl.textContent = 'Atoms: ' + structure.atoms.length;
                bondCountEl.textContent = 'Bonds: ' + (structure.bonds ? structure.bonds.length : 0);
            }
        }

        function showError(message) {
            const errorDiv = document.createElement('div');
            errorDiv.className = 'error-message';
            errorDiv.textContent = message;
            errorDiv.style.cssText = 'position: absolute; top: 10px; left: 50%; transform: translateX(-50%); background: var(--vscode-errorBackground); color: var(--vscode-errorForeground); padding: 8px 16px; border-radius: 4px;';
            document.body.appendChild(errorDiv);
            setTimeout(() => errorDiv.remove(), 5000);
        }

        // Notify extension that webview is ready
        vscode.postMessage({ type: 'ready' });
    </script>
</body>
</html>`;
  }

  /**
   * Handle messages from webview
   */
  private handleMessage(message: any): void {
    switch (message.type) {
      case 'ready':
        // Webview is ready to receive data
        break;

      case 'render':
        this.renderStructure(message.structure, message.config);
        break;

      case 'changeRepresentation':
        this.changeRepresentation(message.mode);
        break;

      case 'resetView':
        this.resetView();
        break;

      case 'toggleAxes':
        this.toggleAxes(message.show);
        break;

      case 'exportImage':
        this.exportImage(message.format);
        break;

      case 'saveImage':
        this.saveImage(message.data);
        break;
    }
  }

  /**
   * Render molecular structure
   */
  private renderStructure(structure: any, config: any): void {
    // This will be handled by the actual Three.js implementation
    // in the renderer class
    if (this.options.onViewChange) {
      this.options.onViewChange({ structure, config });
    }
  }

  /**
   * Change representation mode
   */
  private changeRepresentation(mode: RepresentationMode): void {
    if (this.options.onViewChange) {
      this.options.onViewChange({ representationMode: mode });
    }
  }

  /**
   * Reset camera view
   */
  private resetView(): void {
    if (this.options.onViewChange) {
      this.options.onViewChange({ action: 'resetView' });
    }
  }

  /**
   * Toggle axes display
   */
  private toggleAxes(show: boolean): void {
    if (this.options.onViewChange) {
      this.options.onViewChange({ showAxes: show });
    }
  }

  /**
   * Export image
   */
  private exportImage(format: string): void {
    if (this.options.onViewChange) {
      this.options.onViewChange({ action: 'export', format });
    }
  }

  /**
   * Save image to disk
   */
  private saveImage(data: string): void {
    vscode.window
      .showSaveDialog({
        filters: { Images: ['png'] },
        defaultUri: vscode.Uri.file('molecule.png'),
      })
      .then((uri) => {
        if (uri) {
          vscode.workspace.fs
            .writeFile(uri, Buffer.from(data, 'base64'))
            .then(
              () => {
                vscode.window.showInformationMessage(`Image saved to ${uri.fsPath}`);
              },
              (error) => {
                vscode.window.showErrorMessage(`Failed to save image: ${error}`);
              }
            );
        }
      });
  }

  /**
   * Update webview with new structure
   */
  public update(structure: { atoms: Atom[] }): void {
    this.panel.webview.postMessage({
      type: 'initialize',
      structure,
    });
  }

  /**
   * Notify webview of render completion
   */
  public notifyRenderComplete(): void {
    this.panel.webview.postMessage({
      type: 'renderComplete',
    });
  }

  /**
   * Notify webview of error
   */
  public notifyError(message: string): void {
    this.panel.webview.postMessage({
      type: 'error',
      message,
    });
  }

  /**
   * Dispose of resources
   */
  public dispose(): void {
    ThreeJsWebview.currentPanel = undefined;

    this.panel.dispose();

    while (this.disposables.length) {
      const disposable = this.disposables.pop();
      if (disposable) {
        disposable.dispose();
      }
    }

    if (this.renderer) {
      this.renderer.dispose();
      this.renderer = null;
    }
  }
}

/**
 * Generate a random nonce for CSP
 */
function getNonce(): string {
  let text = '';
  const possible =
    'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
  for (let i = 0; i < 32; i++) {
    text += possible.charAt(Math.floor(Math.random() * possible.length));
  }
  return text;
}
