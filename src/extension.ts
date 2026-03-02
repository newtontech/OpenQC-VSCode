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
    vscode.commands.registerCommand('openqc.visualizeStructure', () => {
      structureViewer.show(vscode.window.activeTextEditor);
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
    vscode.commands.registerCommand('openqc.sidebar.viewResults', (item: JobItem) => {
      vscode.window.showInformationMessage(`Viewing results for: ${item.label} (${item.status})`);
    }),

    // Sidebar: Export data
    vscode.commands.registerCommand('openqc.sidebar.exportData', (item: JobItem) => {
      vscode.window
        .showSaveDialog({
          filters: {
            JSON: ['json'],
            CSV: ['csv'],
            'All Files': ['*'],
          },
        })
        .then(uri => {
          if (uri) {
            vscode.window.showInformationMessage(`Exported ${item.label} to ${uri.fsPath}`);
          }
        });
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

export function deactivate() {
  if (lspManager) {
    lspManager.dispose();
  }
  if (diagnosticsProvider) {
    diagnosticsProvider.dispose();
  }
}
