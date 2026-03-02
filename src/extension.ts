import * as vscode from 'vscode';
import { LSPManager } from './managers/LSPManager';
import { StructureViewer } from './providers/StructureViewer';
import { DataPlotter } from './providers/DataPlotter';
import {
  CompletionProvider,
  DiagnosticsProvider,
  HoverProvider,
  DefinitionProvider,
} from './providers/lsp';
import { FileTypeDetector } from './managers/FileTypeDetector';
import { MoleculeTreeProvider, JobTreeProvider, MoleculeItem, JobItem } from './sidebar';
import { MoleculeViewerPanel } from './visualizers/MoleculeViewerPanel';
import { StructureConverter } from './visualizers/StructureConverter';
import { Molecule3D } from './visualizers/Molecule3D';
import { createParser } from './parsers';

let lspManager: LSPManager;
let structureViewer: StructureViewer;
let dataPlotter: DataPlotter;
let diagnosticsProvider: DiagnosticsProvider;
let fileTypeDetector: FileTypeDetector;
let moleculeProvider: MoleculeTreeProvider;
let jobProvider: JobTreeProvider;

export function activate(context: vscode.ExtensionContext) {
  console.log('OpenQC-VSCode extension is now active!');

  // Set sidebar enabled context
  vscode.commands.executeCommand('setContext', 'openqc.sidebar.enabled', true);

  // Initialize FileTypeDetector
  fileTypeDetector = new FileTypeDetector();

  // Initialize LSP Manager
  lspManager = new LSPManager();

  // Initialize visualization providers
  structureViewer = new StructureViewer(context.extensionUri);
  dataPlotter = new DataPlotter(context.extensionUri);

  // Initialize sidebar providers
  moleculeProvider = new MoleculeTreeProvider(context);
  jobProvider = new JobTreeProvider(context);

  // Initialize LSP providers
  diagnosticsProvider = new DiagnosticsProvider();
  const completionProvider = new CompletionProvider();
  const hoverProvider = new HoverProvider();
  const definitionProvider = new DefinitionProvider();

  // Language IDs for quantum chemistry software
  const languageIds = ['cp2k', 'vasp', 'gaussian', 'orca', 'qe', 'gamess', 'nwchem'];

  // Register language providers
  const disposables = [
    // Completion provider
    vscode.languages.registerCompletionItemProvider(
      languageIds,
      completionProvider,
      '=',
      ' ',
      '\n'
    ),

    // Hover provider
    vscode.languages.registerHoverProvider(languageIds, hoverProvider),

    // Definition provider
    vscode.languages.registerDefinitionProvider(languageIds, definitionProvider),

    // Validation on document change
    vscode.workspace.onDidChangeTextDocument(event => {
      diagnosticsProvider.validateDocument(event.document);
    }),

    // Validation on document save
    vscode.workspace.onDidSaveTextDocument(document => {
      diagnosticsProvider.validateDocument(document);
    }),

    // Clear diagnostics on document close
    vscode.workspace.onDidCloseTextDocument(document => {
      diagnosticsProvider.clearDiagnostics(document);
    }),

    // Visualization commands
    vscode.commands.registerCommand('openqc.visualizeStructure', async () => {
      const editor = vscode.window.activeTextEditor;
      if (!editor) {
        vscode.window.showErrorMessage('No active text editor');
        return;
      }

      const document = editor.document;
      const software = fileTypeDetector.detectSoftware(document);
      if (!software) {
        vscode.window.showErrorMessage('Unsupported file type for visualization');
        return;
      }

      try {
        // Parse the file content
        const content = document.getText();
        const parser = createParser(software, content, document.fileName);

        // Extract atoms using Molecule3D
        const molecule3D = new Molecule3D();
        const atoms = molecule3D.parseAtoms(content, software);

        if (atoms.length === 0) {
          vscode.window.showErrorMessage('No molecular structure found in file');
          return;
        }

        // Convert to XYZ format
        const xyzContent = StructureConverter.atomsToXYZ(atoms, document.fileName);

        // Show the 3D viewer
        MoleculeViewerPanel.createOrShow(context.extensionUri, xyzContent, document.fileName);
      } catch (error) {
        vscode.window.showErrorMessage(`Failed to visualize structure: ${error}`);
      }
    }),

    vscode.commands.registerCommand('openqc.plotData', () => {
      dataPlotter.show(vscode.window.activeTextEditor);
    }),

    vscode.commands.registerCommand('openqc.previewInput', () => {
      structureViewer.showPreview(vscode.window.activeTextEditor);
    }),

    // LSP management commands
    vscode.commands.registerCommand('openqc.startLSP', async () => {
      const editor = vscode.window.activeTextEditor;
      if (editor) {
        await lspManager.startLSPForDocument(editor.document);
      }
    }),

    vscode.commands.registerCommand('openqc.stopLSP', async () => {
      const editor = vscode.window.activeTextEditor;
      if (editor) {
        await lspManager.stopLSPForDocument(editor.document);
      }
    }),

    vscode.commands.registerCommand('openqc.restartLSP', async () => {
      const editor = vscode.window.activeTextEditor;
      if (editor) {
        await lspManager.restartLSPForDocument(editor.document);
      }
    }),

    // Validate current document
    vscode.commands.registerCommand('openqc.validate', async () => {
      const editor = vscode.window.activeTextEditor;
      if (editor) {
        await diagnosticsProvider.validateDocument(editor.document);
        vscode.window.showInformationMessage('Input file validated');
      }
    }),

    // Sidebar: Refresh molecules view
    vscode.commands.registerCommand('openqc.sidebar.refreshMolecules', () => {
      moleculeProvider.refresh();
      vscode.window.showInformationMessage('Molecules refreshed');
    }),

    // Sidebar: Refresh jobs view
    vscode.commands.registerCommand('openqc.sidebar.refreshJobs', () => {
      jobProvider.refresh();
      vscode.window.showInformationMessage('Jobs refreshed');
    }),

    // Sidebar: Open molecule
    vscode.commands.registerCommand('openqc.sidebar.openMolecule', async (item: MoleculeItem) => {
      if (item.filePath) {
        try {
          const doc = await vscode.workspace.openTextDocument(item.filePath);
          await vscode.window.showTextDocument(doc);
        } catch (error) {
          vscode.window.showErrorMessage(`Failed to open file: ${error}`);
        }
      } else {
        vscode.window.showInformationMessage(`Selected molecule: ${item.label} (${item.formula})`);
      }
    }),

    // Sidebar: Delete molecule
    vscode.commands.registerCommand('openqc.sidebar.deleteMolecule', (item: MoleculeItem) => {
      moleculeProvider.removeMolecule(item.id);
      vscode.window.showInformationMessage(`Removed molecule: ${item.label}`);
    }),

    // Sidebar: Run calculation
    vscode.commands.registerCommand('openqc.sidebar.runCalculation', () => {
      vscode.window
        .showInputBox({
          prompt: 'Enter calculation name',
          placeHolder: 'e.g., Geometry Optimization',
        })
        .then(name => {
          if (name) {
            jobProvider.addJob(new JobItem(`job-${Date.now()}`, name, 'queued', 0, 'Gaussian'));
            vscode.window.showInformationMessage(`Started calculation: ${name}`);
          }
        });
    }),

    // Sidebar: View results
    vscode.commands.registerCommand('openqc.sidebar.viewResults', async (item: JobItem) => {
      if (item.status !== 'completed' && item.status !== 'failed') {
        vscode.window.showInformationMessage(`Results not available for ${item.status} jobs`);
        return;
      }

      // Create a results panel
      const panel = vscode.window.createWebviewPanel(
        'openqc.results',
        `Results: ${item.label}`,
        vscode.ViewColumn.Two,
        {
          enableScripts: true,
          retainContextWhenHidden: true,
        }
      );

      // Generate results content
      const resultsData = {
        jobId: item.id,
        jobName: item.label,
        software: item.software,
        status: item.status,
        startTime: item.startTime?.toISOString(),
        endTime: item.endTime?.toISOString(),
        duration:
          typeof item.tooltip === 'string' ? item.tooltip.split('Duration: ')[1] : undefined,
        // In a real implementation, this would come from actual job output
        output: `Results for ${item.label}\n\nSoftware: ${item.software}\nStatus: ${item.status}\n\nSample output data would appear here.\n\nFor completed jobs, this would include:\n- Final energies\n- Optimized geometries\n- Convergence data\n- Properties calculated\n\nFor failed jobs, this would include:\n- Error messages\n- Stack traces\n- Diagnostic information`,
      };

      panel.webview.html = getResultsHtml(resultsData);
    }),

    // Sidebar: Export data
    vscode.commands.registerCommand('openqc.sidebar.exportData', async (item: JobItem) => {
      if (item.status !== 'completed') {
        vscode.window.showWarningMessage('Can only export data from completed jobs');
        return;
      }

      const uri = await vscode.window.showSaveDialog({
        filters: {
          JSON: ['json'],
          CSV: ['csv'],
          'All Files': ['*'],
        },
        defaultUri: vscode.Uri.file(`${item.label.replace(/\s+/g, '_')}_results.json`),
        saveLabel: 'Export Results',
      });

      if (uri) {
        try {
          const exportData = {
            jobId: item.id,
            jobName: item.label,
            software: item.software,
            status: item.status,
            startTime: item.startTime?.toISOString(),
            endTime: item.endTime?.toISOString(),
            timestamp: new Date().toISOString(),
            // In a real implementation, this would include actual job results
            data: {
              energies: [-76.0, -76.1, -76.2],
              forces: [
                [0.01, 0.02, 0.03],
                [-0.01, -0.02, -0.03],
              ],
              converged: true,
            },
          };

          const content = JSON.stringify(exportData, null, 2);
          await vscode.workspace.fs.writeFile(uri, Buffer.from(content));
          vscode.window.showInformationMessage(`Exported ${item.label} to ${uri.fsPath}`);
        } catch (error) {
          vscode.window.showErrorMessage(`Failed to export data: ${error}`);
        }
      }
    }),

    // Sidebar: Cancel job
    vscode.commands.registerCommand('openqc.sidebar.cancelJob', (item: JobItem) => {
      jobProvider.cancelJob(item.id);
      vscode.window.showInformationMessage(`Cancelled job: ${item.label}`);
    }),

    // Sidebar: Restart job
    vscode.commands.registerCommand('openqc.sidebar.restartJob', (item: JobItem) => {
      jobProvider.restartJob(item.id);
      vscode.window.showInformationMessage(`Restarted job: ${item.label}`);
    }),

    // Auto-start LSP on document open
    vscode.workspace.onDidOpenTextDocument(async document => {
      await lspManager.startLSPForDocument(document);
      // Also validate the document
      if (fileTypeDetector.detectSoftware(document)) {
        await diagnosticsProvider.validateDocument(document);
      }
    }),

    // Clean up LSP on document close
    vscode.workspace.onDidCloseTextDocument(async document => {
      await lspManager.stopLSPForDocument(document);
      diagnosticsProvider.clearDiagnostics(document);
    }),
  ];

  context.subscriptions.push(...disposables);
  context.subscriptions.push(diagnosticsProvider);
  context.subscriptions.push(moleculeProvider);
  context.subscriptions.push(jobProvider);

  // Register tree views
  const moleculeTreeView = vscode.window.createTreeView('openqc.molecules', {
    treeDataProvider: moleculeProvider,
    showCollapseAll: true,
  });

  const jobTreeView = vscode.window.createTreeView('openqc.jobs', {
    treeDataProvider: jobProvider,
    showCollapseAll: true,
  });

  context.subscriptions.push(moleculeTreeView);
  context.subscriptions.push(jobTreeView);

  // Start LSP and validate for already open documents
  vscode.window.visibleTextEditors.forEach(async editor => {
    await lspManager.startLSPForDocument(editor.document);
    if (fileTypeDetector.detectSoftware(editor.document)) {
      await diagnosticsProvider.validateDocument(editor.document);
    }
  });
}

function getResultsHtml(data: any): string {
  return `<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Results: ${data.jobName}</title>
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
        #status {
            display: inline-block;
            padding: 2px 8px;
            border-radius: 3px;
            font-size: 11px;
            text-transform: uppercase;
            background: ${data.status === 'completed' ? '#3c9' : '#f66'};
            color: white;
        }
        #content {
            padding: 20px;
        }
        .info-section {
            background: #252526;
            border: 1px solid #3c3c3c;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 15px;
        }
        .info-section h3 {
            margin-top: 0;
            margin-bottom: 10px;
            font-size: 14px;
            color: #0e639c;
        }
        .info-row {
            display: flex;
            padding: 5px 0;
            border-bottom: 1px solid #3c3c3c;
        }
        .info-row:last-child {
            border-bottom: none;
        }
        .info-label {
            width: 150px;
            color: #9cdcfe;
            font-weight: 500;
        }
        .info-value {
            flex: 1;
            color: #ce9178;
        }
        #output {
            background: #1e1e1e;
            border: 1px solid #3c3c3c;
            border-radius: 5px;
            padding: 15px;
            font-family: 'Consolas', monospace;
            font-size: 12px;
            white-space: pre-wrap;
            color: #cccccc;
        }
    </style>
</head>
<body>
    <div id="header">
        <div id="title">${data.jobName}</div>
        <div><span id="status">${data.status}</span> ${data.software}</div>
    </div>
    <div id="content">
        <div class="info-section">
            <h3>Job Information</h3>
            <div class="info-row">
                <span class="info-label">Job ID:</span>
                <span class="info-value">${data.jobId}</span>
            </div>
            <div class="info-row">
                <span class="info-label">Software:</span>
                <span class="info-value">${data.software}</span>
            </div>
            <div class="info-row">
                <span class="info-label">Status:</span>
                <span class="info-value">${data.status}</span>
            </div>
            <div class="info-row">
                <span class="info-label">Duration:</span>
                <span class="info-value">${data.duration || 'N/A'}</span>
            </div>
            ${
              data.startTime
                ? `
            <div class="info-row">
                <span class="info-label">Started:</span>
                <span class="info-value">${new Date(data.startTime).toLocaleString()}</span>
            </div>
            `
                : ''
            }
            ${
              data.endTime
                ? `
            <div class="info-row">
                <span class="info-label">Completed:</span>
                <span class="info-value">${new Date(data.endTime).toLocaleString()}</span>
            </div>
            `
                : ''
            }
        </div>
        <div class="info-section">
            <h3>Output</h3>
            <div id="output">${data.output}</div>
        </div>
    </div>
</body>
</html>`;
}

export function deactivate() {
  if (lspManager) {
    lspManager.dispose();
  }
  if (diagnosticsProvider) {
    diagnosticsProvider.dispose();
  }
}
