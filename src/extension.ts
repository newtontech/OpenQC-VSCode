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
import { registerCalculatorCommands } from './ase/calculatorCommands';
import { registerMigrationCommands } from './commands/migrationCommands';
import { registerPythonBackendCommands } from './commands/pythonBackendCommands';
import { registerScientificBridgeCommands } from './commands/scientificBridgeCommands';
import { registerAnalyzerCommands } from './commands/analyzerCommands';
import { registerJobCommands } from './commands/jobCommands';
import { registerAICommands } from './ai/aiCommands';
import { registerExportCommands } from './commands/exportCommands';
import { FileTypeDetector } from './managers/FileTypeDetector';
import { MoleculeTreeProvider, JobTreeProvider, MoleculeItem } from './sidebar';
import { OpenQCConverterProvider } from './sidebar/OpenQCConverterProvider';
import { MoleculeViewerPanel } from './visualizers/MoleculeViewerPanel';
import { OpenQCViewerPanel } from './visualizers/OpenQCViewerPanel';
import { ViewerStructurePipeline } from './visualizers/ViewerStructurePipeline';
import { registerFormatConversionCommands } from './commands/formatConversionCommands';
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
let viewerStructurePipeline: ViewerStructurePipeline;
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
  viewerStructurePipeline = new ViewerStructurePipeline(context);

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
        const fileName = document.fileName;
        const structure = (await viewerStructurePipeline.parse(document, software)).structure;

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

  for (const document of vscode.workspace.textDocuments) {
    if (fileTypeDetector.detectSoftware(document)) {
      if (languageFeaturePolicy.autoStartLsp) {
        void lspManager.startLSPForDocument(document);
      }
      if (languageFeaturePolicy.validateWithLocalDiagnostics) {
        void diagnosticsProvider.validateDocument(document);
      }
    }
  }

  // Register ASE commands
  registerASECommands(context);
  registerCalculatorCommands(context);
  registerPythonBackendCommands(context);
  registerScientificBridgeCommands(context);
  registerAnalyzerCommands(context);
  registerJobCommands(context, jobProvider);
  registerMigrationCommands(context);
  registerAICommands(context);
  registerExportCommands(
    context,
    () => OpenQCViewerPanel.getCurrentStructureData(),
    () => OpenQCViewerPanel.saveCurrentStructureToSource()
  );

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
  viewerStructurePipeline?.clear();
  if (lspManager) {
    await lspManager.dispose();
  }
  if (diagnosticsProvider) {
    diagnosticsProvider.dispose();
  }
}
