/**
 * Visualize Command - Display molecular structure in 3D viewer
 */

import * as vscode from 'vscode';
import { VisualizerView } from '../views/visualizer';

export class VisualizeCommand {
    constructor(
        private context: vscode.ExtensionContext,
        private outputChannel: vscode.OutputChannel
    ) {}
    
    register() {
        const disposable = vscode.commands.registerCommand(
            'openqc.visualize',
            async () => {
                const editor = vscode.window.activeTextEditor;
                if (!editor) {
                    vscode.window.showErrorMessage('No active editor');
                    return;
                }
                
                const document = editor.document;
                const languageId = document.languageId;
                
                // Validate language
                const supportedLanguages = ['gaussian', 'vasp', 'quantumespresso', 'orca'];
                if (!supportedLanguages.includes(languageId)) {
                    vscode.window.showErrorMessage(
                        `Unsupported file type: ${languageId}. ` +
                        `Supported types: ${supportedLanguages.join(', ')}`
                    );
                    return;
                }
                
                this.outputChannel.appendLine(`Visualizing ${document.fileName}`);
                
                // Create or show the visualizer panel
                VisualizerView.createOrShow(this.context.extensionUri, document);
            }
        );
        
        this.context.subscriptions.push(disposable);
    }
}
