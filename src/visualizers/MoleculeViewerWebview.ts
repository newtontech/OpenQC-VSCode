import * as vscode from 'vscode';

export class MoleculeViewerWebview {
  /**
   * Generate the HTML content for the molecule viewer webview
   * @returns HTML string with embedded JavaScript
   */
  static generateWebviewHTML(): string {
    const csp = this.getCSP();
    const nonce = this.getNonce();

    return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta http-equiv="Content-Security-Policy" content="${csp}">
  <title>3D Molecular Structure Viewer</title>
  <style>
    body {
      margin: 0;
      padding: 0;
      overflow: hidden;
      font-family: var(--vscode-font-family);
      background-color: var(--vscode-editor-background);
      color: var(--vscode-editor-foreground);
    }
    #ngl-container {
      width: 100vw;
      height: 100vh;
      display: flex;
      flex-direction: column;
    }
    #ngl-canvas {
      flex: 1;
      width: 100%;
      height: 100%;
    }
    #controls {
      padding: 10px;
      background-color: var(--vscode-editorHoverWidget-background);
      border-bottom: 1px solid var(--vscode-panel-border);
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
    }
    .btn {
      background-color: var(--vscode-button-background);
      color: var(--vscode-button-foreground);
      border: none;
      padding: 6px 12px;
      cursor: pointer;
      font-size: 13px;
      border-radius: 2px;
    }
    .btn:hover {
      background-color: var(--vscode-button-hoverBackground);
    }
    .select {
      background-color: var(--vscode-dropdown-background);
      color: var(--vscode-dropdown-foreground);
      border: 1px solid var(--vscode-dropdown-border);
      padding: 4px 8px;
      font-size: 13px;
      border-radius: 2px;
    }
    #status {
      padding: 5px 10px;
      font-size: 12px;
      color: var(--vscode-descriptionForeground);
    }
    #error {
      display: none;
      padding: 10px;
      background-color: var(--vscode-errorBackground);
      color: var(--vscode-errorForeground);
      border-radius: 3px;
      margin: 10px;
    }
    #loading {
      display: flex;
      align-items: center;
      justify-content: center;
      height: 100%;
      font-size: 14px;
      color: var(--vscode-descriptionForeground);
    }
    .spinner {
      border: 3px solid var(--vscode-editorWidget-border);
      border-top: 3px solid var(--vscode-button-background);
      border-radius: 50%;
      width: 24px;
      height: 24px;
      animation: spin 1s linear infinite;
      margin-right: 10px;
    }
    @keyframes spin {
      0% { transform: rotate(0deg); }
      100% { transform: rotate(360deg); }
    }
  </style>
</head>
<body>
  <div id="ngl-container">
    <div id="controls">
      <select id="representation" class="select">
        <option value="ball-stick">Ball & Stick</option>
        <option value="space-filling">Space Filling</option>
        <option value="licorice">Licorice</option>
        <option value="cartoon">Cartoon</option>
        <option value="wireframe">Wireframe</option>
      </select>
      <button id="reset" class="btn">Reset View</button>
      <button id="auto-rotate" class="btn">Auto Rotate</button>
      <button id="export-image" class="btn">Export Image</button>
      <span id="status">Ready</span>
    </div>
    <div id="loading"><div class="spinner"></div><span>Loading 3D viewer...</span></div>
    <canvas id="ngl-canvas"></canvas>
    <div id="error"></div>
  </div>

  <script src="https://unpkg.com/ngl@2.2.1/dist/ngl.js"></script>
  <script>
    (function() {
      const vscode = acquireVsCodeApi();
      let stage = null;
      let structureComp = null;
      let autoRotate = false;

      // Initialize NGL Stage
      function initStage() {
        document.getElementById('loading').style.display = 'none';

        stage = new NGL.Stage('ngl-canvas', {
          backgroundColor: 'var(--vscode-editor-background)',
          backgroundColor2: 'var(--vscode-editor-background)',
          impostor: true,
          quality: 'medium',
          lightIntensity: 1.5
        });

        stage.handleResize();
        window.addEventListener('resize', () => stage.handleResize());

        // Set up mouse interactions
        stage.mouseControls.add('cameraPick', (stage, pickingProxy) => {
          if (pickingProxy && pickingProxy.atom) {
            const atom = pickingProxy.atom;
            const atomName = atom.qualifiedName();
            const pos = atom.position;
            showStatus(\`\${atomName} at (\${pos.x.toFixed(3)}, \${pos.y.toFixed(3)}, \${pos.z.toFixed(3)})\`);
          }
        });

        setupControls();
        showStatus('3D viewer loaded. Waiting for structure...');
      }

      // Setup UI controls
      function setupControls() {
        document.getElementById('representation').addEventListener('change', (e) => {
          if (structureComp) {
            changeRepresentation(e.target.value);
          }
        });

        document.getElementById('reset').addEventListener('click', () => {
          if (stage) {
            stage.autoView();
          }
        });

        document.getElementById('auto-rotate').addEventListener('click', () => {
          autoRotate = !autoRotate;
          if (stage) {
            if (autoRotate) {
              stage.setSpin(true);
              showStatus('Auto rotation: ON');
            } else {
              stage.setSpin(false);
              showStatus('Auto rotation: OFF');
            }
          }
        });

        document.getElementById('export-image').addEventListener('click', () => {
          if (stage) {
            stage.makeImage({
              factor: 2,
              antialias: true,
              trim: false,
              transparent: false
            }).then((blob) => {
              vscode.postMessage({
                type: 'exportImage',
                data: blob
              });
              showStatus('Image exported');
            });
          }
        });
      }

      // Load structure from data
      function loadStructure(data) {
        if (!stage) return;

        // Clear existing structure
        if (structureComp) {
          stage.removeComponent(structureComp);
          structureComp = null;
        }

        // Parse structure data
        let structure;
        if (data.xyz) {
          structure = NGL.Structure.parseString(data.xyz, { ext: 'xyz' });
        } else if (data.json) {
          const parsed = JSON.parse(data.json);
          const xyz = jsonToXYZ(parsed);
          structure = NGL.Structure.parseString(xyz, { ext: 'xyz' });
        } else {
          showError('No structure data provided');
          return;
        }

        // Add to stage
        structureComp = stage.addComponentFromObject(structure);

        // Apply default representation
        structureComp.addRepresentation('ball-stick', {
          aspectRatio: 1.2,
          radiusSize: 0.15,
          multipleBond: true,
          bondScale: 0.5
        });

        // Auto-center and zoom
        stage.autoView();
        showStatus(\`Loaded \${structure.atomCount} atoms\`);
      }

      // Convert JSON data to XYZ format
      function jsonToXYZ(data) {
        const atoms = data.atoms || [];
        const comment = data.comment || 'molecule';
        let xyz = atoms.length + '\\n';
        xyz += comment + '\\n';
        atoms.forEach(atom => {
          xyz += \`\${atom.elem} \${atom.x.toFixed(6)} \${atom.y.toFixed(6)} \${atom.z.toFixed(6)}\\n\`;
        });
        return xyz;
      }

      // Change representation style
      function changeRepresentation(type) {
        if (!structureComp) return;

        structureComp.clearRepresentations();

        const props = {
          aspectRatio: type === 'space-filling' ? 1 : 1.2,
          radiusSize: type === 'space-filling' ? 1 : 0.15,
          multipleBond: true,
          bondScale: type === 'licorice' ? 0.3 : 0.5
        };

        structureComp.addRepresentation(type, props);
        showStatus(\`Representation: \${type}\`);
      }

      // Show status message
      function showStatus(msg) {
        document.getElementById('status').textContent = msg;
      }

      // Show error message
      function showError(msg) {
        const errorEl = document.getElementById('error');
        errorEl.textContent = msg;
        errorEl.style.display = 'block';
        setTimeout(() => {
          errorEl.style.display = 'none';
        }, 5000);
      }

      // Handle messages from extension
      window.addEventListener('message', (event) => {
        const message = event.data;

        switch (message.type) {
          case 'initialize':
            if (!stage) {
              initStage();
            }
            if (message.structure) {
              loadStructure(message.structure);
            }
            break;

          case 'loadStructure':
            loadStructure(message.structure);
            break;

          case 'changeRepresentation':
            if (message.representation) {
              changeRepresentation(message.representation);
            }
            break;

          case 'resetView':
            if (stage) {
              stage.autoView();
            }
            break;

          case 'setSpin':
            autoRotate = message.spin;
            if (stage) {
              stage.setSpin(autoRotate);
            }
            break;
        }
      });

      // Initialize when DOM is ready
      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', initStage);
      } else {
        initStage();
      }
    })();
  </script>
</body>
</html>`;
  }

  /**
   * Get webview options
   */
  static getWebviewOptions(): vscode.WebviewOptions {
    return {
      enableScripts: true,
    };
  }

  /**
   * Generate Content Security Policy for the webview
   */
  static getCSP(): string {
    return [
      `default-src 'none';`,
      `script-src https://unpkg.com https://cdn.jsdelivr.net 'nonce-${this.getNonce()}';`,
      `style-src 'unsafe-inline';`,
      `img-src https: data:;`,
      `font-src https: data:;`,
    ].join(' ');
  }

  /**
   * Generate a random nonce for CSP
   */
  private static getNonce(): string {
    let text = '';
    const possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    for (let i = 0; i < 32; i++) {
      text += possible.charAt(Math.floor(Math.random() * possible.length));
    }
    return text;
  }
}
