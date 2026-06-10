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
import { registerASECommands } from './ase/commands';
import { registerMigrationCommands } from './commands/migrationCommands';
import { registerPythonBackendCommands } from './commands/pythonBackendCommands';
import { registerScientificBridgeCommands } from './commands/scientificBridgeCommands';
import { registerAnalyzerCommands } from './commands/analyzerCommands';
import { registerAICommands } from './ai/aiCommands';
import { registerExportCommands } from './commands/exportCommands';
import { FileTypeDetector } from './managers/FileTypeDetector';
import { MoleculeTreeProvider, JobTreeProvider, MoleculeItem, JobItem } from './sidebar';
import { OpenQCConverterProvider } from './sidebar/OpenQCConverterProvider';
import { MoleculeViewerPanel } from './visualizers/MoleculeViewerPanel';
import { OpenQCViewerPanel } from './visualizers/OpenQCViewerPanel';
import { StructureConverter } from './visualizers/StructureConverter';
import { Molecule3D } from './visualizers/Molecule3D';
import {
  createOpenQCStructure,
  molecularStructureToOpenQCStructure,
  poscarToOpenQCStructure,
} from './structures/converters';
import type { OpenQCStructure } from './structures/OpenQCStructure';
import { createParser } from './parsers';
import { registerFormatConversionCommands } from './commands/formatConversionCommands';
import { renderResultsWebviewHtml } from './webviews/resultsWebview';
import { Logger, LogLevel } from './utils/Logger';
import { getLanguageFeaturePolicy, readLanguageFeatureMode } from './languageFeatures';
import { listBundledLspServers } from './lsp/registry';
import { registerDslAuthoringContextCommand } from './lsp/dslAuthoringContext';

let lspManager: LSPManager;
let structureViewer: StructureViewer;
let dataPlotter: DataPlotter;
let completionProvider: CompletionProvider;
let hoverProvider: HoverProvider;
let definitionProvider: DefinitionProvider;
let diagnosticsProvider: DiagnosticsProvider;
let fileTypeDetector: FileTypeDetector;
let moleculeProvider: MoleculeTreeProvider;
let jobProvider: JobTreeProvider;
let converterProvider: OpenQCConverterProvider;
const logger = Logger.getInstance();

export function activate(context: vscode.ExtensionContext) {
  // Configure logger from VS Code settings
  const config = vscode.workspace.getConfiguration('openqc.logging');
  logger.setConfig({
    level: Logger.parseLogLevel(config.get<string>('level', 'info')),
    showUserMessages: config.get<boolean>('showUserMessages', true),
  });

  logger.info('OpenQC-VSCode extension activated');
  const languageFeaturePolicy = getLanguageFeaturePolicy(readLanguageFeatureMode());

  // Set sidebar enabled context
  vscode.commands.executeCommand('setContext', 'openqc.sidebar.enabled', true);

  // Initialize FileTypeDetector
  fileTypeDetector = new FileTypeDetector();

  // Initialize LSP Manager
  lspManager = new LSPManager(context);

  // Initialize visualization providers
  structureViewer = new StructureViewer(context.extensionUri);
  dataPlotter = new DataPlotter(context.extensionUri);

  // Initialize sidebar providers
  moleculeProvider = new MoleculeTreeProvider(context);
  jobProvider = new JobTreeProvider(context);

  // Initialize OpenQC Converter Sidebar Provider
  converterProvider = new OpenQCConverterProvider(context.extensionUri);
  context.subscriptions.push(
    vscode.window.registerWebviewViewProvider(OpenQCConverterProvider.viewType, converterProvider)
  );

  // Set converter enabled context
  vscode.commands.executeCommand('setContext', 'openqc.converterEnabled', true);

  diagnosticsProvider = new DiagnosticsProvider();

  registerFormatConversionCommands(context);

  // Language IDs for bundled OpenQC-managed LSP servers.
  const languageIds = listBundledLspServers().map(server => server.languageId);

  // Register language providers
  const disposables: vscode.Disposable[] = [];

  if (languageFeaturePolicy.registerLocalProviders) {
    completionProvider = new CompletionProvider();
    hoverProvider = new HoverProvider();
    definitionProvider = new DefinitionProvider();

    disposables.push(
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
      vscode.languages.registerDefinitionProvider(languageIds, definitionProvider)
    );
  }

  if (languageFeaturePolicy.registerLocalDiagnostics) {
    disposables.push(
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
      })
    );
  }

  disposables.push(
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
        const content = document.getText();
        const fileName = document.fileName;
        let structure: OpenQCStructure | undefined;

        // 1. Try format-specific DTO path for periodic formats
        if (software === 'VASP') {
          const basename = fileName.split(/[/\\]/).pop()?.toUpperCase() ?? '';
          if (basename === 'POSCAR' || basename === 'CONTCAR') {
            try {
              structure = poscarToOpenQCStructure(content, fileName);
            } catch {
              logger.warn('POSCAR DTO parse failed, falling back to legacy path');
            }
          }
        }

        // 2. Try StructureConverter DTO path for other formats
        if (!structure) {
          try {
            const converter = new StructureConverter();
            const molecular = converter.autoConvert(content, fileName);
            if (molecular.atoms.length > 0) {
              structure = molecularStructureToOpenQCStructure(molecular, {
                sourceSoftware: software,
                sourceParser: 'native',
              });
            }
          } catch {
            logger.warn('StructureConverter DTO parse failed, falling back to legacy path');
          }
        }

        // 3. Final legacy fallback: Molecule3D → atoms → createOpenQCStructure
        if (!structure) {
          const molecule3D = new Molecule3D();
          const atoms = molecule3D.parseAtoms(content, software);
          if (atoms.length > 0) {
            structure = createOpenQCStructure(atoms, {
              name: fileName,
              sourceSoftware: software,
            });
          }
        }

        if (!structure || structure.atoms.length === 0) {
          vscode.window.showErrorMessage('No molecular structure found in file');
          return;
        }

        // Show the DTO-driven viewer
        OpenQCViewerPanel.createOrShow(context.extensionUri, structure, fileName);
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

    vscode.commands.registerCommand('openqc.lsp.showStatus', () => {
      lspManager.showStatus();
    }),

    vscode.commands.registerCommand('openqc.lsp.showLogs', () => {
      lspManager.showLogs();
    }),

    vscode.commands.registerCommand('openqc.lsp.restartCurrent', async () => {
      await lspManager.restartCurrent();
    }),

    vscode.commands.registerCommand('openqc.lsp.selectExecutable', async () => {
      await lspManager.selectExecutable();
    }),

    vscode.commands.registerCommand('openqc.lsp.generateCompatibilityReport', async () => {
      await lspManager.generateCompatibilityReport();
    }),

    // Validate current document
    vscode.commands.registerCommand('openqc.validate', async () => {
      const editor = vscode.window.activeTextEditor;
      if (editor) {
        if (languageFeaturePolicy.validateWithLocalDiagnostics) {
          await diagnosticsProvider.validateDocument(editor.document);
        } else {
          await lspManager.startLSPForDocument(editor.document);
        }
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
          retainContextWhenHidden: true,
          localResourceRoots: [],
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

      panel.webview.html = renderResultsWebviewHtml(resultsData, panel.webview.cspSource);
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
      if (languageFeaturePolicy.autoStartLsp) {
        await lspManager.startLSPForDocument(document);
      }
      if (
        languageFeaturePolicy.validateWithLocalDiagnostics &&
        fileTypeDetector.detectSoftware(document)
      ) {
        await diagnosticsProvider.validateDocument(document);
      }
    }),

    // Clean up LSP on document close
    vscode.workspace.onDidCloseTextDocument(async document => {
      await lspManager.stopLSPForDocument(document);
      if (languageFeaturePolicy.registerLocalDiagnostics) {
        diagnosticsProvider.clearDiagnostics(document);
      }
    })
  );

  // Register ASE commands
  registerPythonBackendCommands(context);
  registerScientificBridgeCommands(context);
  registerAnalyzerCommands(context);
  registerMigrationCommands(context);
  registerAICommands(context);
  registerExportCommands(context);

  // Register DSL authoring context command for agent workflows
  registerDslAuthoringContextCommand(context, () => lspManager);

  // Show Converter Panel command
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.showConverterPanel', () => {
      vscode.commands.executeCommand('openqc-converter-panel.focus');
    })
  );

  context.subscriptions.push(...disposables);

  console.log('OpenQC-VSCode: All providers registered successfully!');
}

export async function deactivate() {
  if (lspManager) {
    await lspManager.dispose();
  }
  if (diagnosticsProvider) {
    diagnosticsProvider.dispose();
  }
}
