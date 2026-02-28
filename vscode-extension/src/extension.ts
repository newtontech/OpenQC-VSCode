/**
 * OpenQC - VSCode Extension for Quantum Chemistry
 * Main extension entry point
 */

import * as vscode from 'vscode';
import { VisualizeCommand } from './commands/visualize';
import { ConvertCommand } from './commands/convert';
import { ValidateCommand } from './commands/validate';
import { AIAssistCommand } from './commands/aiAssist';
import { ConnectServerCommand } from './commands/connectServer';
import { SubmitJobCommand } from './commands/submitJob';
import { GaussianCompletionProvider } from './providers/completion';
import { QCDocumentSymbolProvider } from './providers/symbols';
import { QCHoverProvider } from './providers/hover';

let outputChannel: vscode.OutputChannel;

export function activate(context: vscode.ExtensionContext) {
    outputChannel = vscode.window.createOutputChannel('OpenQC');
    outputChannel.appendLine('OpenQC extension is activating...');
    
    // Register commands
    const commands = [
        new VisualizeCommand(context, outputChannel),
        new ConvertCommand(context, outputChannel),
        new ValidateCommand(context, outputChannel),
        new AIAssistCommand(context, outputChannel),
        new ConnectServerCommand(context, outputChannel),
        new SubmitJobCommand(context, outputChannel),
    ];
    
    commands.forEach(cmd => cmd.register());
    
    // Register language providers
    const languages = ['gaussian', 'vasp', 'quantumespresso', 'orca'];
    
    languages.forEach(lang => {
        // Completion provider
        context.subscriptions.push(
            vscode.languages.registerCompletionItemProvider(
                { language: lang },
                new GaussianCompletionProvider(),
                ' ', '\t'
            )
        );
        
        // Hover provider
        context.subscriptions.push(
            vscode.languages.registerHoverProvider(
                { language: lang },
                new QCHoverProvider(lang)
            )
        );
        
        // Document symbol provider
        context.subscriptions.push(
            vscode.languages.registerDocumentSymbolProvider(
                { language: lang },
                new QCDocumentSymbolProvider(lang)
            )
        );
    });
    
    outputChannel.appendLine('OpenQC extension activated successfully!');
    
    // Show welcome message on first activation
    const hasShownWelcome = context.globalState.get<boolean>('openqc.welcomeShown');
    if (!hasShownWelcome) {
        vscode.window.showInformationMessage(
            'Welcome to OpenQC! Your quantum chemistry toolkit is ready.',
            'Open Documentation',
            'View Examples'
        ).then(selection => {
            if (selection === 'Open Documentation') {
                vscode.env.openExternal(
                    vscode.Uri.parse('https://github.com/newtontech/OpenQC#readme')
                );
            } else if (selection === 'View Examples') {
                vscode.commands.executeCommand('vscode.openFolder');
            }
        });
        context.globalState.update('openqc.welcomeShown', true);
    }
}

export function deactivate() {
    if (outputChannel) {
        outputChannel.dispose();
    }
}
