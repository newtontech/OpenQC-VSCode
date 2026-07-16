/**
 * OpenQC bundled webview HTML generator.
 *
 * Generates a self-contained webview that loads 3Dmol.js and viewer scripts
 * from bundled extension assets via webview.asWebviewUri().
 * No CDN dependency — safe for Remote SSH and offline use.
 *
 * @module webviews/openqcViewerWebview
 */

import * as vscode from 'vscode';
import { generateNonce } from '../utils/nonce';

export class OpenQCViewerWebview {
  /**
   * Generate the complete HTML for the bundled OpenQC viewer.
   */
  static generateHTML(webview: vscode.Webview, extensionUri: vscode.Uri): string {
    const nonce = getNonce();

    // Bundled asset URIs
    const viewerJsUri = webview.asWebviewUri(
      vscode.Uri.joinPath(extensionUri, 'media', 'openqc-viewer.js')
    );
    const viewerCssUri = webview.asWebviewUri(
      vscode.Uri.joinPath(extensionUri, 'media', 'openqc-viewer.css')
    );
    const threeDmolUri = webview.asWebviewUri(
      vscode.Uri.joinPath(extensionUri, 'media', 'vendor', '3dmol', '3Dmol-min.js')
    );

    const csp = [
      `default-src 'none'`,
      `script-src 'nonce-${nonce}'`,
      `style-src 'nonce-${nonce}'`,
      `img-src data: blob:`,
      `font-src data:`,
    ].join('; ');

    return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta http-equiv="Content-Security-Policy" content="${csp}">
  <title>OpenQC Structure Viewer</title>
  <link rel="stylesheet" nonce="${nonce}" href="${viewerCssUri}">
</head>
<body>
  <div id="openqc-container">
    <!-- Toolbar -->
    <div id="toolbar">
      <div class="toolbar-group">
        <select id="style-select" class="select">
          <option value="ball-stick">Ball &amp; Stick</option>
          <option value="spacefill">Space Filling</option>
          <option value="line">Line</option>
          <option value="cartoon">Cartoon</option>
        </select>
      </div>
      <div class="toolbar-separator"></div>
      <div class="toolbar-group">
        <button id="btn-reset" class="btn" title="Reset Camera">Reset</button>
        <button id="btn-rotate" class="btn secondary" title="Auto Rotate">Rotate</button>
        <button id="btn-labels" class="btn secondary" title="Toggle Atom Labels">Labels</button>
        <button id="btn-measure" class="btn secondary" title="Measure Distance, Angle, or Dihedral">Measure</button>
        <button id="btn-export" class="btn secondary" title="Export PNG">Export</button>
      </div>
      <div class="toolbar-separator"></div>
      <div class="toolbar-group">
        <button id="btn-toggle-cell" class="btn secondary" title="Toggle Unit Cell">Cell</button>
      </div>
    </div>

    <!-- Editing controls -->
    <div id="editing-controls">
      <label for="edit-element">Atom</label>
      <input type="text" id="edit-element" value="C" maxlength="2" title="Element symbol">
      <input type="number" id="edit-x" value="0" step="0.1" title="Add X coordinate">
      <input type="number" id="edit-y" value="0" step="0.1" title="Add Y coordinate">
      <input type="number" id="edit-z" value="0" step="0.1" title="Add Z coordinate">
      <button id="btn-add-atom" class="btn secondary" title="Add Atom">Add</button>
      <button id="btn-delete-atom" class="btn secondary" title="Delete Selected Atom">Delete</button>
      <span class="editing-separator"></span>
      <label for="move-dx">Move</label>
      <input type="number" id="move-dx" value="0" step="0.1" title="Delta X">
      <input type="number" id="move-dy" value="0" step="0.1" title="Delta Y">
      <input type="number" id="move-dz" value="0" step="0.1" title="Delta Z">
      <button id="btn-move-atom" class="btn secondary" title="Move Selected Atom">Move</button>
      <span class="editing-separator"></span>
      <label for="bond-from">Bond</label>
      <input type="number" id="bond-from" value="1" min="1" step="1" title="First atom index">
      <input type="number" id="bond-to" value="2" min="1" step="1" title="Second atom index">
      <select id="bond-order" class="select" title="Bond order">
        <option value="1">1</option>
        <option value="2">2</option>
        <option value="3">3</option>
      </select>
      <button id="btn-add-bond" class="btn secondary" title="Add or Update Bond">Bond</button>
      <button id="btn-delete-bond" class="btn secondary" title="Delete Bond">Unbond</button>
      <span id="bond-count">0 bonds</span>
      <span class="editing-separator"></span>
      <label>Free</label>
      <label class="constraint-toggle" title="Free selected atom along X">
        <input type="checkbox" id="constraint-x" checked>
        X
      </label>
      <label class="constraint-toggle" title="Free selected atom along Y">
        <input type="checkbox" id="constraint-y" checked>
        Y
      </label>
      <label class="constraint-toggle" title="Free selected atom along Z">
        <input type="checkbox" id="constraint-z" checked>
        Z
      </label>
      <button id="btn-set-constraints" class="btn secondary" title="Apply Selected Atom Constraints">Apply</button>
      <span class="editing-separator"></span>
      <button id="btn-undo-edit" class="btn secondary" title="Undo Edit">Undo</button>
      <button id="btn-redo-edit" class="btn secondary" title="Redo Edit">Redo</button>
      <button id="btn-save-source" class="btn secondary" title="Save Edited Structure to Source File">Save Source</button>
      <button id="btn-export-structure" class="btn secondary" title="Export Edited Structure">Export Structure</button>
      <span id="dirty-indicator">Clean</span>
    </div>

    <!-- Periodic controls (hidden unless periodic structure) -->
    <div id="periodic-controls">
      <label>Supercell:</label>
      <input type="number" id="sc-na" value="2" min="1" max="5">
      <span>×</span>
      <input type="number" id="sc-nb" value="2" min="1" max="5">
      <span>×</span>
      <input type="number" id="sc-nc" value="2" min="1" max="5">
      <button id="btn-supercell" class="btn secondary">Generate</button>
    </div>

    <!-- Trajectory controls (hidden unless trajectory) -->
    <div id="trajectory-controls">
      <button id="btn-play" class="btn secondary">▶ Play</button>
      <input type="range" id="frame-slider" min="0" max="0" value="0">
      <span id="frame-info">Frame 0/0</span>
    </div>

    <!-- Loading overlay -->
    <div id="loading-overlay">
      <div class="spinner"></div>
      <span>Loading OpenQC Viewer...</span>
    </div>

    <!-- 3D canvas -->
    <div id="viewer-canvas"></div>

    <!-- WebGL fallback -->
    <div id="webgl-fallback">
      <h2>WebGL Not Available</h2>
      <p>WebGL is required for 3D molecular visualization.</p>
      <p>This may occur in some Remote SSH configurations or older browsers.</p>
    </div>

    <!-- Error panel -->
    <div id="error-panel"></div>

    <!-- Status bar -->
    <div id="status-bar">
      <div class="status-left">
        <span id="status-text">Ready</span>
      </div>
      <span id="atom-count"></span>
    </div>
  </div>

  <!-- Load 3Dmol.js from bundled extension assets -->
  <script nonce="${nonce}" src="${threeDmolUri}"></script>
  <!-- Load OpenQC viewer script -->
  <script nonce="${nonce}" src="${viewerJsUri}"></script>
</body>
</html>`;
  }

  /**
   * Get webview options for the bundled viewer.
   */
  static getWebviewOptions(extensionUri: vscode.Uri): vscode.WebviewOptions {
    return {
      enableScripts: true,
      localResourceRoots: [vscode.Uri.joinPath(extensionUri, 'media')],
    };
  }
}

const getNonce = generateNonce;
