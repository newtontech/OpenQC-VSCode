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
import { registerAICommands } from './ai/aiCommands';
import { FileTypeDetector } from './managers/FileTypeDetector';
import { MoleculeTreeProvider, JobTreeProvider, MoleculeItem, JobItem } from './sidebar';
import { OpenQCConverterProvider } from './sidebar/OpenQCConverterProvider';
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
let converterProvider: OpenQCConverterProvider;

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

  // Initialize OpenQC Converter Sidebar Provider
  converterProvider = new OpenQCConverterProvider(context.extensionUri);
  context.subscriptions.push(
    vscode.window.registerWebviewViewProvider(
      OpenQCConverterProvider.viewType,
      converterProvider
    )
  );
  
  // Set converter enabled context
  vscode.commands.executeCommand('setContext', 'openqc.converterEnabled', true);

  // Initialize LSP providers
  diagnosticsProvider = new DiagnosticsProvider();

  // Register all commands
  registerCommands(context);
  registerASECommands(context);
  registerMigrationCommands(context);
  registerAICommands(context);
  
  console.log('OpenQC-VSCode: All providers registered successfully!');
}

function registerCommands(context: vscode.ExtensionContext): void {
  // ... 其他命令注册代码 ...
  
  // Show Converter Panel command
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.showConverterPanel', () => {
      vscode.commands.executeCommand('openqc-converter-panel.focus');
    })
  );
}

export function deactivate() {
  console.log('OpenQC-VSCode extension is now deactivated');
}
