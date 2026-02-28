import * as vscode from 'vscode';
import * as path from 'path';
import { Molecule3D } from '../visualizers/Molecule3D';
import { FileTypeDetector, QuantumChemistrySoftware } from '../managers/FileTypeDetector';

export class StructureViewer {
    private panel: vscode.WebviewPanel | undefined;
    private extensionUri: vscode.Uri;
    private molecule3D: Molecule3D;
    private fileTypeDetector: FileTypeDetector;

    constructor(extensionUri: vscode.Uri) {
        this.extensionUri = extensionUri;
        this.molecule3D = new Molecule3D();
        this.fileTypeDetector = new FileTypeDetector();
    }

    show(editor: vscode.TextEditor | undefined): void {
        if (!editor) {
            vscode.window.showWarningMessage('No active editor found');
            return;
        }

        const software = this.fileTypeDetector.detectSoftware(editor.document);
        if (!software) {
            vscode.window.showWarningMessage('Unsupported file type for structure visualization');
            return;
        }

        const column = vscode.ViewColumn.Two;
        
        if (this.panel) {
            this.panel.reveal(column);
        } else {
            this.panel = vscode.window.createWebviewPanel(
                'openqcStructureViewer',
                'OpenQC: Molecular Structure',
                column,
                {
                    enableScripts: true,
                    localResourceRoots: [this.extensionUri],
                    retainContextWhenHidden: true
                }
            );

            this.panel.onDidDispose(() => {
                this.panel = undefined;
            });
        }

        this.updateContent(editor.document, software);
    }

    showPreview(editor: vscode.TextEditor | undefined): void {
        if (!editor) {
            vscode.window.showWarningMessage('No active editor found');
            return;
        }

        const software = this.fileTypeDetector.detectSoftware(editor.document);
        if (!software) {
            vscode.window.showWarningMessage('Unsupported file type for input preview');
            return;
        }

        const column = vscode.ViewColumn.Two;
        
        if (this.panel) {
            this.panel.reveal(column);
        } else {
            this.panel = vscode.window.createWebviewPanel(
                'openqcInputPreview',
                'OpenQC: Input Preview',
                column,
                {
                    enableScripts: true,
                    localResourceRoots: [this.extensionUri],
                    retainContextWhenHidden: true
                }
            );

            this.panel.onDidDispose(() => {
                this.panel = undefined;
            });
        }

        this.updatePreviewContent(editor.document, software);
    }

    private updateContent(document: vscode.TextDocument, software: QuantumChemistrySoftware): void {
        if (!this.panel) return;

        const content = document.getText();
        const atoms = this.molecule3D.parseAtoms(content, software);
        
        this.panel.webview.html = this.get3DViewerHtml(atoms, software);
    }

    private updatePreviewContent(document: vscode.TextDocument, software: QuantumChemistrySoftware): void {
        if (!this.panel) return;

        const content = document.getText();
        const previewData = this.parseInputPreview(content, software);
        
        this.panel.webview.html = this.getPreviewHtml(previewData, software);
    }

    private parseInputPreview(content: string, software: QuantumChemistrySoftware): any {
        // Parse input file and extract structured data
        const lines = content.split('\n');
        const data: any = {
            software,
            sections: [],
            parameters: [],
            atoms: []
        };

        // Extract parameters based on software type
        switch (software) {
            case 'CP2K':
                data.sections = this.extractCp2kSections(content);
                break;
            case 'VASP':
                data.parameters = this.extractVaspParameters(content);
                break;
            case 'Gaussian':
                data.parameters = this.extractGaussianParameters(content);
                break;
            case 'ORCA':
                data.parameters = this.extractOrcaParameters(content);
                break;
            case 'Quantum ESPRESSO':
                data.sections = this.extractQeSections(content);
                break;
            case 'GAMESS':
                data.parameters = this.extractGamessParameters(content);
                break;
            case 'NWChem':
                data.sections = this.extractNwchemSections(content);
                break;
        }

        return data;
    }

    private extractCp2kSections(content: string): any[] {
        const sections: any[] = [];
        const sectionRegex = /&(\w+)[\s\S]*?&END\s*\1/gi;
        let match;
        while ((match = sectionRegex.exec(content)) !== null) {
            sections.push({
                name: match[1],
                content: match[0].substring(0, 200) + '...'
            });
        }
        return sections;
    }

    private extractVaspParameters(content: string): any[] {
        const params: any[] = [];
        const lines = content.split('\n');
        for (const line of lines) {
            const match = line.match(/^(\w+)\s*=\s*(.+)/);
            if (match) {
                params.push({ name: match[1], value: match[2] });
            }
        }
        return params;
    }

    private extractGaussianParameters(content: string): any[] {
        const params: any[] = [];
        const routeMatch = content.match(/^#(.+)$/m);
        if (routeMatch) {
            params.push({ name: 'Route Section', value: routeMatch[1].trim() });
        }
        const chkMatch = content.match(/^%chk=(.+)$/m);
        if (chkMatch) {
            params.push({ name: 'Checkpoint File', value: chkMatch[1].trim() });
        }
        return params;
    }

    private extractOrcaParameters(content: string): any[] {
        const params: any[] = [];
        const simpleInputMatch = content.match(/^!(.+)$/m);
        if (simpleInputMatch) {
            params.push({ name: 'Simple Input', value: simpleInputMatch[1].trim() });
        }
        return params;
    }

    private extractQeSections(content: string): any[] {
        const sections: any[] = [];
        const sectionNames = ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'];
        for (const name of sectionNames) {
            const regex = new RegExp(`&${name}[\\s\\S]*?&END`, 'i');
            const match = content.match(regex);
            if (match) {
                sections.push({ name, content: 'Present' });
            }
        }
        return sections;
    }

    private extractGamessParameters(content: string): any[] {
        const params: any[] = [];
        const groupRegex = /^\s*\$(\w+)[\s\S]*?\$END/gi;
        let match;
        while ((match = groupRegex.exec(content)) !== null) {
            params.push({ name: `Group: ${match[1]}`, value: 'Present' });
        }
        return params;
    }

    private extractNwchemSections(content: string): any[] {
        const sections: any[] = [];
        const blockRegex = /^(geometry|basis|scf|dft|mp2|ccsd|title)\s+\w*/gim;
        let match;
        while ((match = blockRegex.exec(content)) !== null) {
            sections.push({ name: match[1], content: 'Present' });
        }
        return sections;
    }

    private get3DViewerHtml(atoms: any[], software: QuantumChemistrySoftware): string {
        const atomsJson = JSON.stringify(atoms);
        
        return `<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>OpenQC: Molecular Structure - ${software}</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.6/3Dmol-min.js"></script>
    <style>
        body { 
            margin: 0; 
            padding: 0; 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1e1e1e;
            color: #cccccc;
        }
        #header {
            padding: 10px 20px;
            background: #252526;
            border-bottom: 1px solid #3c3c3c;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        #title {
            font-size: 14px;
            font-weight: 600;
        }
        #software-badge {
            background: #0e639c;
            color: white;
            padding: 4px 8px;
            border-radius: 3px;
            font-size: 11px;
            text-transform: uppercase;
        }
        #container {
            position: relative;
            width: 100%;
            height: calc(100vh - 50px);
        }
        #viewer {
            width: 100%;
            height: 100%;
        }
        #controls {
            position: absolute;
            bottom: 20px;
            left: 20px;
            background: rgba(30, 30, 30, 0.9);
            padding: 10px;
            border-radius: 5px;
            border: 1px solid #3c3c3c;
        }
        .control-btn {
            background: #0e639c;
            color: white;
            border: none;
            padding: 6px 12px;
            margin: 2px;
            border-radius: 3px;
            cursor: pointer;
            font-size: 11px;
        }
        .control-btn:hover {
            background: #1177bb;
        }
        #info {
            position: absolute;
            top: 20px;
            right: 20px;
            background: rgba(30, 30, 30, 0.9);
            padding: 10px;
            border-radius: 5px;
            border: 1px solid #3c3c3c;
            font-size: 12px;
        }
    </style>
</head>
<body>
    <div id="header">
        <span id="title">Molecular Structure Viewer</span>
        <span id="software-badge">${software}</span>
    </div>
    <div id="container">
        <div id="viewer"></div>
        <div id="controls">
            <button class="control-btn" onclick="setStyle('stick')">Stick</button>
            <button class="control-btn" onclick="setStyle('sphere')">Sphere</button>
            <button class="control-btn" onclick="setStyle('line')">Line</button>
            <button class="control-btn" onclick="setStyle('cartoon')">Cartoon</button>
            <br>
            <button class="control-btn" onclick="viewer.spin(true)">Spin</button>
            <button class="control-btn" onclick="viewer.spin(false)">Stop</button>
            <button class="control-btn" onclick="viewer.zoomTo()">Reset View</button>
        </div>
        <div id="info">
            <div>Atoms: <span id="atom-count">0</span></div>
        </div>
    </div>
    <script>
        let viewer = null;
        const atoms = ${atomsJson};
        
        document.addEventListener('DOMContentLoaded', function() {
            const element = document.getElementById('viewer');
            const config = { backgroundColor: '#1e1e1e' };
            viewer = $3Dmol.createViewer(element, config);
            
            if (atoms && atoms.length > 0) {
                const m = viewer.addModel();
                m.addAtoms(atoms);
                viewer.zoomTo();
                viewer.setStyle({}, {stick: {}});
                viewer.render();
                
                document.getElementById('atom-count').textContent = atoms.length;
            } else {
                element.innerHTML = '<div style="padding: 50px; text-align: center; color: #888;">' +
                    'No atomic coordinates found in file.<br>' +
                    'Supported formats: XYZ, POSCAR, Gaussian, etc.</div>';
            }
        });
        
        function setStyle(style) {
            if (!viewer) return;
            viewer.removeAllModels();
            const m = viewer.addModel();
            m.addAtoms(atoms);
            
            if (style === 'stick') {
                viewer.setStyle({}, {stick: {}});
            } else if (style === 'sphere') {
                viewer.setStyle({}, {sphere: {}});
            } else if (style === 'line') {
                viewer.setStyle({}, {line: {}});
            } else if (style === 'cartoon') {
                viewer.setStyle({}, {cartoon: {}});
            }
            viewer.render();
        }
    </script>
</body>
</html>`;
    }

    private getPreviewHtml(data: any, software: QuantumChemistrySoftware): string {
        return `<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>OpenQC: Input Preview - ${software}</title>
    <style>
        body { 
            margin: 0; 
            padding: 0; 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1e1e1e;
            color: #cccccc;
        }
        #header {
            padding: 15px 20px;
            background: #252526;
            border-bottom: 1px solid #3c3c3c;
        }
        #title {
            font-size: 16px;
            font-weight: 600;
            margin-bottom: 5px;
        }
        #software {
            color: #0e639c;
            font-size: 12px;
        }
        #content {
            padding: 20px;
        }
        .section {
            background: #252526;
            border: 1px solid #3c3c3c;
            border-radius: 5px;
            margin-bottom: 15px;
            overflow: hidden;
        }
        .section-header {
            background: #2d2d30;
            padding: 10px 15px;
            font-weight: 600;
            font-size: 13px;
            border-bottom: 1px solid #3c3c3c;
        }
        .section-content {
            padding: 15px;
        }
        .param-row {
            display: flex;
            padding: 5px 0;
            border-bottom: 1px solid #3c3c3c;
        }
        .param-row:last-child {
            border-bottom: none;
        }
        .param-name {
            width: 200px;
            color: #9cdcfe;
            font-family: 'Consolas', monospace;
            font-size: 12px;
        }
        .param-value {
            flex: 1;
            color: #ce9178;
            font-family: 'Consolas', monospace;
            font-size: 12px;
        }
        .badge {
            display: inline-block;
            background: #0e639c;
            color: white;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 10px;
            margin-right: 5px;
        }
    </style>
</head>
<body>
    <div id="header">
        <div id="title">Input File Preview</div>
        <div id="software">${software}</div>
    </div>
    <div id="content">
        ${data.sections.length > 0 ? `
        <div class="section">
            <div class="section-header">Input Sections</div>
            <div class="section-content">
                ${data.sections.map((s: any) => `
                    <div class="param-row">
                        <span class="badge">SECTION</span>
                        <span class="param-name">${s.name}</span>
                    </div>
                `).join('')}
            </div>
        </div>
        ` : ''}
        
        ${data.parameters.length > 0 ? `
        <div class="section">
            <div class="section-header">Parameters</div>
            <div class="section-content">
                ${data.parameters.map((p: any) => `
                    <div class="param-row">
                        <span class="param-name">${p.name}</span>
                        <span class="param-value">${p.value}</span>
                    </div>
                `).join('')}
            </div>
        </div>
        ` : ''}
    </div>
</body>
</html>`;
    }
}