/**
 * Visualizer View - 3D molecular visualization using 3Dmol.js
 */

import * as vscode from 'vscode';

export class VisualizerView {
    public static currentPanel: VisualizerView | undefined;
    public static readonly viewType = 'openqc.visualizer';
    
    private readonly _panel: vscode.WebviewPanel;
    private readonly _extensionUri: vscode.Uri;
    private _disposables: vscode.Disposable[] = [];
    
    public static createOrShow(extensionUri: vscode.Uri, document: vscode.TextDocument) {
        const column = vscode.window.activeTextEditor
            ? vscode.window.activeTextEditor.viewColumn
            : undefined;
        
        // If we already have a panel, show it
        if (VisualizerView.currentPanel) {
            VisualizerView.currentPanel._panel.reveal(column);
            VisualizerView.currentPanel._update(document);
            return;
        }
        
        // Create a new panel
        const panel = vscode.window.createWebviewPanel(
            VisualizerView.viewType,
            'OpenQC Visualizer',
            column || vscode.ViewColumn.One,
            {
                enableScripts: true,
                retainContextWhenHidden: true
            }
        );
        
        VisualizerView.currentPanel = new VisualizerView(panel, extensionUri, document);
    }
    
    private constructor(
        panel: vscode.WebviewPanel,
        extensionUri: vscode.Uri,
        document: vscode.TextDocument
    ) {
        this._panel = panel;
        this._extensionUri = extensionUri;
        
        // Set the webview's initial html content
        this._update(document);
        
        // Listen for when the panel is disposed
        this._panel.onDidDispose(() => this.dispose(), null, this._disposables);
        
        // Handle messages from the webview
        this._panel.webview.onDidReceiveMessage(
            message => {
                switch (message.command) {
                    case 'alert':
                        vscode.window.showInformationMessage(message.text);
                        break;
                    case 'export':
                        this._handleExport(message.format, message.data);
                        break;
                }
            },
            null,
            this._disposables
        );
    }
    
    private async _update(document: vscode.TextDocument) {
        const content = document.getText();
        const filename = document.fileName.split('/').pop() || 'structure';
        
        // Parse structure based on file type
        const structure = this._parseStructure(content, document.languageId);
        
        this._panel.webview.html = this._getHtmlForWebview(
            this._panel.webview,
            structure,
            filename
        );
    }
    
    private _parseStructure(content: string, languageId: string): any {
        // Basic structure parsing
        // In production, this would call the Python parser
        
        const atoms: { element: string; x: number; y: number; z: number }[] = [];
        const lines = content.split('\n');
        
        switch (languageId) {
            case 'gaussian':
                // Find the geometry section
                let inGeometry = false;
                for (const line of lines) {
                    if (line.trim().match(/^[A-Za-z]+\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+$/)) {
                        const parts = line.trim().split(/\s+/);
                        if (parts.length >= 4) {
                            atoms.push({
                                element: parts[0],
                                x: parseFloat(parts[1]),
                                y: parseFloat(parts[2]),
                                z: parseFloat(parts[3])
                            });
                        }
                    }
                }
                break;
                
            case 'vasp':
                // Parse POSCAR format
                if (lines.length > 8) {
                    const scale = parseFloat(lines[1]);
                    const cell = [
                        lines[2].trim().split(/\s+/).map(Number),
                        lines[3].trim().split(/\s+/).map(Number),
                        lines[4].trim().split(/\s+/).map(Number)
                    ];
                    // ... more parsing needed
                }
                break;
                
            default:
                // Try to find XYZ-like coordinates
                for (const line of lines) {
                    const match = line.match(/^([A-Za-z]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)$/);
                    if (match) {
                        atoms.push({
                            element: match[1],
                            x: parseFloat(match[2]),
                            y: parseFloat(match[3]),
                            z: parseFloat(match[4])
                        });
                    }
                }
        }
        
        return { atoms };
    }
    
    private _handleExport(format: string, data: string) {
        // Handle export requests from webview
        vscode.window.showInformationMessage(`Export as ${format} requested`);
    }
    
    private _getHtmlForWebview(webview: vscode.Webview, structure: any, filename: string): string {
        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>OpenQC Visualizer - ${filename}</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.1/3Dmol.min.js"></script>
    <style>
        body {
            margin: 0;
            padding: 0;
            font-family: var(--vscode-font-family);
            background-color: var(--vscode-editor-background);
            color: var(--vscode-editor-foreground);
        }
        #viewer-container {
            width: 100%;
            height: 100vh;
            position: relative;
        }
        #viewer {
            width: 100%;
            height: 100%;
        }
        .toolbar {
            position: absolute;
            top: 10px;
            left: 10px;
            z-index: 100;
            display: flex;
            flex-direction: column;
            gap: 5px;
        }
        .toolbar button {
            padding: 8px 12px;
            background: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            border-radius: 3px;
            cursor: pointer;
            font-size: 12px;
        }
        .toolbar button:hover {
            background: var(--vscode-button-hoverBackground);
        }
        .info-panel {
            position: absolute;
            bottom: 10px;
            left: 10px;
            background: rgba(0,0,0,0.7);
            padding: 10px;
            border-radius: 5px;
            font-size: 12px;
        }
        .controls {
            position: absolute;
            top: 10px;
            right: 10px;
            background: rgba(0,0,0,0.7);
            padding: 10px;
            border-radius: 5px;
        }
        .controls label {
            display: block;
            margin: 5px 0;
        }
    </style>
</head>
<body>
    <div id="viewer-container">
        <div id="viewer"></div>
        <div class="toolbar">
            <button onclick="setStyle('stick')">Stick</button>
            <button onclick="setStyle('sphere')">Sphere</button>
            <button onclick="setStyle('cartoon')">Cartoon</button>
            <button onclick="resetView()">Reset View</button>
            <button onclick="toggleSpin()">Spin</button>
            <button onclick="exportImage()">Export PNG</button>
        </div>
        <div class="info-panel" id="info">
            <strong>${filename}</strong><br>
            Atoms: <span id="atom-count">${structure.atoms.length}</span>
        </div>
        <div class="controls">
            <label>
                Background:
                <select onchange="setBackground(this.value)">
                    <option value="black">Black</option>
                    <option value="white">White</option>
                    <option value="gray">Gray</option>
                </select>
            </label>
        </div>
    </div>
    
    <script>
        let viewer;
        let spinning = false;
        
        function init() {
            viewer = $3Dmol.createViewer('viewer', {
                defaultcolors: $3Dmol.rasmolElementColors
            });
            
            // Add atoms
            const structureData = ${JSON.stringify(structure)};
            const atoms = structureData.atoms.map((a, i) => ({
                elem: a.element,
                x: a.x,
                y: a.y,
                z: a.z,
                serial: i + 1
            }));
            
            viewer.addModel(atoms, 'xyz');
            viewer.setStyle({}, {stick: {}, sphere: {scale: 0.3}});
            viewer.zoomTo();
            viewer.render();
        }
        
        function setStyle(style) {
            const styles = {
                stick: {stick: {}},
                sphere: {sphere: {}},
                cartoon: {cartoon: {}}
            };
            viewer.setStyle({}, styles[style] || styles.stick);
            viewer.render();
        }
        
        function resetView() {
            viewer.zoomTo();
            viewer.render();
        }
        
        function toggleSpin() {
            spinning = !spinning;
            viewer.spin(spinning);
        }
        
        function setBackground(color) {
            viewer.setBackgroundColor(color);
            viewer.render();
        }
        
        function exportImage() {
            const img = viewer.pngURI();
            const link = document.createElement('a');
            link.download = 'structure.png';
            link.href = img;
            link.click();
        }
        
        // Initialize on load
        window.onload = init;
    </script>
</body>
</html>`;
    }
    
    public dispose() {
        VisualizerView.currentPanel = undefined;
        
        this._panel.dispose();
        
        while (this._disposables.length) {
            const x = this._disposables.pop();
            if (x) {
                x.dispose();
            }
        }
    }
}
