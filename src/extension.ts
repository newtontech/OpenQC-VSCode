import * as vscode from 'vscode';
import { LSPManager } from './managers/LSPManager';
import { StructureViewer } from './providers/StructureViewer';
import { DataPlotter } from './providers/DataPlotter';

let lspManager: LSPManager;
let structureViewer: StructureViewer;
let dataPlotter: DataPlotter;

export function activate(context: vscode.ExtensionContext) {
    console.log('OpenQC-VSCode extension is now active!');

    // Initialize LSP Manager
    lspManager = new LSPManager();
    
    // Initialize visualization providers
    structureViewer = new StructureViewer(context.extensionUri);
    dataPlotter = new DataPlotter(context.extensionUri);

    // Register commands
    const disposables = [
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

        // Auto-start LSP on document open
        vscode.workspace.onDidOpenTextDocument(async (document) => {
            await lspManager.startLSPForDocument(document);
        }),

        // Clean up LSP on document close
        vscode.workspace.onDidCloseTextDocument(async (document) => {
            await lspManager.stopLSPForDocument(document);
        })
    ];

    context.subscriptions.push(...disposables);

    // Start LSP for already open documents
    if (vscode.window.activeTextEditor) {
        lspManager.startLSPForDocument(vscode.window.activeTextEditor.document);
    }
}

export function deactivate() {
    if (lspManager) {
        lspManager.dispose();
    }
}