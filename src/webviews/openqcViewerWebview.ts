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
      vscode.Uri.joinPath(extensionUri, 'node_modules', '3dmol', 'build', '3Dmol-min.js')
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
        <button id="btn-export" class="btn secondary" title="Export PNG">Export</button>
      </div>
      <div class="toolbar-separator"></div>
      <div class="toolbar-group">
        <button id="btn-toggle-cell" class="btn secondary" title="Toggle Unit Cell">Cell</button>
      </div>
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
      localResourceRoots: [
        vscode.Uri.joinPath(extensionUri, 'media'),
        vscode.Uri.joinPath(extensionUri, 'node_modules', '3dmol'),
      ],
    };
  }
}

function getNonce(): string {
  let text = '';
  const chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
  for (let i = 0; i < 32; i++) {
    text += chars.charAt(Math.floor(Math.random() * chars.length));
  }
  return text;
}
