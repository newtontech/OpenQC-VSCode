import * as vscode from 'vscode';
import { FileTypeDetector, QuantumChemistrySoftware } from '../managers/FileTypeDetector';

export class DataPlotter {
    private panel: vscode.WebviewPanel | undefined;
    private extensionUri: vscode.Uri;
    private fileTypeDetector: FileTypeDetector;

    constructor(extensionUri: vscode.Uri) {
        this.extensionUri = extensionUri;
        this.fileTypeDetector = new FileTypeDetector();
    }

    show(editor: vscode.TextEditor | undefined): void {
        if (!editor) {
            vscode.window.showWarningMessage('No active editor found');
            return;
        }

        const software = this.fileTypeDetector.detectSoftware(editor.document);
        if (!software) {
            vscode.window.showWarningMessage('Unsupported file type for data plotting');
            return;
        }

        const column = vscode.ViewColumn.Two;
        
        if (this.panel) {
            this.panel.reveal(column);
        } else {
            this.panel = vscode.window.createWebviewPanel(
                'openqcDataPlotter',
                'OpenQC: Data Plotter',
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

    private updateContent(document: vscode.TextDocument, software: QuantumChemistrySoftware): void {
        if (!this.panel) return;

        const content = document.getText();
        const plotData = this.extractPlotData(content, software);
        
        this.panel.webview.html = this.getPlotterHtml(plotData, software);
    }

    private extractPlotData(content: string, software: QuantumChemistrySoftware): any {
        const data: any = {
            software,
            plots: [],
            metadata: {}
        };

        // Extract data based on software type
        switch (software) {
            case 'CP2K':
                data.plots = this.extractCp2kData(content);
                break;
            case 'VASP':
                data.plots = this.extractVaspData(content);
                break;
            case 'Gaussian':
                data.plots = this.extractGaussianData(content);
                break;
            case 'ORCA':
                data.plots = this.extractOrcaData(content);
                break;
            case 'Quantum ESPRESSO':
                data.plots = this.extractQeData(content);
                break;
            case 'GAMESS':
                data.plots = this.extractGamessData(content);
                break;
            case 'NWChem':
                data.plots = this.extractNwchemData(content);
                break;
        }

        return data;
    }

    private extractCp2kData(content: string): any[] {
        // Extract energy convergence data
        const plots: any[] = [];
        const energyPattern = /Total energy\s*:\s*([-\d.]+)/gi;
        const energies: number[] = [];
        let match;
        while ((match = energyPattern.exec(content)) !== null) {
            energies.push(parseFloat(match[1]));
        }
        if (energies.length > 0) {
            plots.push({
                title: 'Energy Convergence',
                type: 'line',
                x: energies.map((_, i) => i + 1),
                y: energies,
                xLabel: 'Iteration',
                yLabel: 'Energy (a.u.)'
            });
        }
        return plots;
    }

    private extractVaspData(content: string): any[] {
        const plots: any[] = [];
        // Extract KPOINTS data
        const kpointsMatch = content.match(/KPOINTS[\s\S]*?^\s*$/m);
        if (kpointsMatch) {
            const lines = kpointsMatch[0].split('\n');
            if (lines.length > 3) {
                const grid = lines[3].trim().split(/\s+/).map(Number);
                plots.push({
                    title: 'K-point Grid',
                    type: 'bar',
                    x: ['Kx', 'Ky', 'Kz'],
                    y: grid,
                    xLabel: 'Direction',
                    yLabel: 'K-points'
                });
            }
        }
        return plots;
    }

    private extractGaussianData(content: string): any[] {
        const plots: any[] = [];
        // Extract SCF energies
        const scfPattern = /SCF Done:\s*E\([^)]+\)\s*=\s*([-\d.]+)/gi;
        const energies: number[] = [];
        let match;
        while ((match = scfPattern.exec(content)) !== null) {
            energies.push(parseFloat(match[1]));
        }
        if (energies.length > 0) {
            plots.push({
                title: 'SCF Energy Convergence',
                type: 'line',
                x: energies.map((_, i) => i + 1),
                y: energies,
                xLabel: 'Cycle',
                yLabel: 'Energy (a.u.)'
            });
        }
        return plots;
    }

    private extractOrcaData(content: string): any[] {
        const plots: any[] = [];
        // Extract energy data from output-like content
        const energyPattern = /FINAL SINGLE POINT ENERGY\s+([-\d.]+)/gi;
        const energies: number[] = [];
        let match;
        while ((match = energyPattern.exec(content)) !== null) {
            energies.push(parseFloat(match[1]));
        }
        if (energies.length > 0) {
            plots.push({
                title: 'Energy',
                type: 'bar',
                x: ['Final Energy'],
                y: [energies[energies.length - 1]],
                xLabel: '',
                yLabel: 'Energy (a.u.)'
            });
        }
        return plots;
    }

    private extractQeData(content: string): any[] {
        const plots: any[] = [];
        // Extract SCF convergence
        const energyPattern = /total energy\s*=\s*([-\d.]+)/gi;
        const energies: number[] = [];
        let match;
        while ((match = energyPattern.exec(content)) !== null) {
            energies.push(parseFloat(match[1]));
        }
        if (energies.length > 0) {
            plots.push({
                title: 'Total Energy Convergence',
                type: 'line',
                x: energies.map((_, i) => i + 1),
                y: energies,
                xLabel: 'Iteration',
                yLabel: 'Energy (Ry)'
            });
        }
        return plots;
    }

    private extractGamessData(content: string): any[] {
        const plots: any[] = [];
        // Extract SCF data
        const energyPattern = /TOTAL ENERGY\s*=\s*([-\d.]+)/gi;
        const energies: number[] = [];
        let match;
        while ((match = energyPattern.exec(content)) !== null) {
            energies.push(parseFloat(match[1]));
        }
        if (energies.length > 0) {
            plots.push({
                title: 'Total Energy',
                type: 'line',
                x: energies.map((_, i) => i + 1),
                y: energies,
                xLabel: 'Iteration',
                yLabel: 'Energy (a.u.)'
            });
        }
        return plots;
    }

    private extractNwchemData(content: string): any[] {
        const plots: any[] = [];
        // Extract energy data
        const energyPattern = /Total SCF energy\s*=\s*([-\d.]+)/gi;
        const energies: number[] = [];
        let match;
        while ((match = energyPattern.exec(content)) !== null) {
            energies.push(parseFloat(match[1]));
        }
        if (energies.length > 0) {
            plots.push({
                title: 'SCF Energy',
                type: 'line',
                x: energies.map((_, i) => i + 1),
                y: energies,
                xLabel: 'Iteration',
                yLabel: 'Energy (a.u.)'
            });
        }
        return plots;
    }

    private getPlotterHtml(data: any, software: QuantumChemistrySoftware): string {
        const plotsJson = JSON.stringify(data.plots);
        
        return `<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>OpenQC: Data Plotter - ${software}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
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
            padding: 20px;
        }
        .plot-container {
            background: #252526;
            border: 1px solid #3c3c3c;
            border-radius: 5px;
            margin-bottom: 20px;
            overflow: hidden;
        }
        .plot-title {
            padding: 10px 15px;
            background: #2d2d30;
            font-weight: 600;
            font-size: 13px;
            border-bottom: 1px solid #3c3c3c;
        }
        .plot-area {
            height: 400px;
            padding: 10px;
        }
        .no-data {
            padding: 50px;
            text-align: center;
            color: #888;
        }
    </style>
</head>
<body>
    <div id="header">
        <span id="title">Data Plotter</span>
        <span id="software-badge">${software}</span>
    </div>
    <div id="container">
        <div id="plots"></div>
    </div>
    <script>
        const plots = ${plotsJson};
        const container = document.getElementById('plots');
        
        if (plots && plots.length > 0) {
            plots.forEach((plot, index) => {
                const plotDiv = document.createElement('div');
                plotDiv.className = 'plot-container';
                plotDiv.innerHTML = '<div class="plot-title">' + plot.title + '</div>' +
                    '<div class="plot-area" id="plot-' + index + '"></div>';
                container.appendChild(plotDiv);
                
                const layout = {
                    paper_bgcolor: '#252526',
                    plot_bgcolor: '#1e1e1e',
                    font: { color: '#cccccc' },
                    xaxis: { 
                        title: plot.xLabel,
                        gridcolor: '#3c3c3c',
                        linecolor: '#3c3c3c'
                    },
                    yaxis: { 
                        title: plot.yLabel,
                        gridcolor: '#3c3c3c',
                        linecolor: '#3c3c3c'
                    },
                    margin: { t: 20, r: 20, b: 60, l: 80 }
                };
                
                const trace = {
                    x: plot.x,
                    y: plot.y,
                    type: plot.type,
                    mode: 'lines+markers',
                    line: { color: '#0e639c', width: 2 },
                    marker: { color: '#0e639c', size: 8 }
                };
                
                if (plot.type === 'bar') {
                    trace.type = 'bar';
                    trace.marker = { color: '#0e639c' };
                }
                
                Plotly.newPlot('plot-' + index, [trace], layout, { responsive: true });
            });
        } else {
            container.innerHTML = '<div class="no-data">' +
                '<h3>No Plot Data Found</h3>' +
                '<p>This file does not contain extractable plot data.</p>' +
                '<p>Supported data types: SCF energies, convergence data, k-points, etc.</p>' +
                '</div>';
        }
    </script>
</body>
</html>`;
    }
}