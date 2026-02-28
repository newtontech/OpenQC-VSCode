/**
 * AI Assist Command - AI-powered parameter modification and suggestions
 */

import * as vscode from 'vscode';

export class AIAssistCommand {
    private suggestions = [
        'Change functional to B3LYP',
        'Change functional to PBE0',
        'Change functional to M06-2X',
        'Change basis set to 6-31G(d)',
        'Change basis set to def2-TZVP',
        'Add solvation (PCM)',
        'Add dispersion correction (D3)',
        'Optimize geometry',
        'Calculate frequency',
        'Calculate NMR',
        'Calculate UV-Vis',
        'Add counterpoise correction',
        'Increase integration grid',
        'Tighten convergence criteria'
    ];
    
    constructor(
        private context: vscode.ExtensionContext,
        private outputChannel: vscode.OutputChannel
    ) {}
    
    register() {
        const disposable = vscode.commands.registerCommand(
            'openqc.aiAssist',
            async () => {
                const editor = vscode.window.activeTextEditor;
                if (!editor) {
                    vscode.window.showErrorMessage('No active editor');
                    return;
                }
                
                const config = vscode.workspace.getConfiguration('openqc');
                if (!config.get<boolean>('enableAI')) {
                    const enable = await vscode.window.showInformationMessage(
                        'AI features are disabled. Enable them?',
                        'Yes', 'No'
                    );
                    if (enable !== 'Yes') {
                        return;
                    }
                    config.update('enableAI', true, vscode.ConfigurationTarget.Global);
                }
                
                // Show quick pick for common modifications
                const action = await vscode.window.showQuickPick(
                    [
                        { label: '$(sparkle) Quick Modification', description: 'Apply common parameter changes' },
                        { label: '$(comment) Natural Language', description: 'Describe what you want to change' },
                        { label: '$(check) Analyze Input', description: 'Get AI analysis of current input' },
                        { label: '$(lightbulb) Suggest Improvements', description: 'Get optimization suggestions' }
                    ],
                    { placeHolder: 'Select AI action' }
                );
                
                if (!action) {
                    return;
                }
                
                switch (action.label) {
                    case '$(sparkle) Quick Modification':
                        await this.quickModification(editor);
                        break;
                    case '$(comment) Natural Language':
                        await this.naturalLanguageModification(editor);
                        break;
                    case '$(check) Analyze Input':
                        await this.analyzeInput(editor);
                        break;
                    case '$(lightbulb) Suggest Improvements':
                        await this.suggestImprovements(editor);
                        break;
                }
            }
        );
        
        this.context.subscriptions.push(disposable);
    }
    
    private async quickModification(editor: vscode.TextEditor): Promise<void> {
        const selection = await vscode.window.showQuickPick(
            this.suggestions,
            { placeHolder: 'Select modification' }
        );
        
        if (!selection) {
            return;
        }
        
        this.outputChannel.appendLine(`Applying modification: ${selection}`);
        
        // TODO: Call MCP server to apply modification
        vscode.window.showInformationMessage(
            `Modification "${selection}" will be applied via MCP server. ` +
            `This feature requires the OpenQC MCP server to be running.`
        );
    }
    
    private async naturalLanguageModification(editor: vscode.TextEditor): Promise<void> {
        const input = await vscode.window.showInputBox({
            prompt: 'Describe the modification you want',
            placeHolder: 'e.g., "Change the functional to B3LYP and add solvation"'
        });
        
        if (!input) {
            return;
        }
        
        this.outputChannel.appendLine(`Natural language request: ${input}`);
        
        // TODO: Call MCP server with natural language
        vscode.window.showInformationMessage(
            `Processing: "${input}". ` +
            `This feature requires the OpenQC MCP server.`
        );
    }
    
    private async analyzeInput(editor: vscode.TextEditor): Promise<void> {
        const content = editor.document.getText();
        this.outputChannel.appendLine('Analyzing input file...');
        
        // Basic analysis
        const analysis: string[] = [];
        
        if (content.includes('B3LYP')) {
            analysis.push('• Using B3LYP functional (hybrid GGA)');
        } else if (content.includes('PBE0')) {
            analysis.push('• Using PBE0 functional (hybrid GGA)');
        } else if (content.includes('M06')) {
            analysis.push('• Using M06 functional (meta-GGA)');
        }
        
        if (content.includes('opt')) {
            analysis.push('• Geometry optimization requested');
        }
        if (content.includes('freq')) {
            analysis.push('• Frequency calculation requested');
        }
        if (content.includes('scrf')) {
            analysis.push('• Solvation model included');
        }
        
        if (analysis.length === 0) {
            analysis.push('• Basic calculation setup detected');
        }
        
        vscode.window.showInformationMessage(
            'Analysis Complete',
            { modal: true, detail: analysis.join('\n') }
        );
    }
    
    private async suggestImprovements(editor: vscode.TextEditor): Promise<void> {
        const content = editor.document.getText();
        
        const suggestions: string[] = [];
        
        // Basic suggestions
        if (content.includes('B3LYP') && !content.includes('empirical')) {
            suggestions.push('Consider adding D3 dispersion correction for better non-covalent interactions');
        }
        
        if (!content.includes('int=ultrafine')) {
            suggestions.push('Consider using int=ultrafine for more accurate integration');
        }
        
        if (content.includes('opt') && !content.includes('freq')) {
            suggestions.push('Consider adding freq to verify the optimized structure is a minimum');
        }
        
        if (suggestions.length === 0) {
            suggestions.push('Your input looks well-optimized!');
        }
        
        vscode.window.showInformationMessage(
            'Suggestions',
            { modal: true, detail: suggestions.join('\n') }
        );
    }
}
