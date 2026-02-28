/**
 * Validate Command - Check quantum chemistry input files for errors
 */

import * as vscode from 'vscode';
import { spawn } from 'child_process';

interface ValidationResult {
    valid: boolean;
    errors: string[];
    warnings: string[];
    suggestions: string[];
}

export class ValidateCommand {
    constructor(
        private context: vscode.ExtensionContext,
        private outputChannel: vscode.OutputChannel
    ) {}
    
    register() {
        const disposable = vscode.commands.registerCommand(
            'openqc.validate',
            async () => {
                const editor = vscode.window.activeTextEditor;
                if (!editor) {
                    vscode.window.showErrorMessage('No active editor');
                    return;
                }
                
                const document = editor.document;
                const content = document.getText();
                
                this.outputChannel.appendLine(`Validating ${document.fileName}`);
                
                try {
                    const result = await this.runValidator(
                        content,
                        document.languageId
                    );
                    
                    this.showValidationResult(result);
                } catch (error) {
                    vscode.window.showErrorMessage(
                        `Validation failed: ${error}`
                    );
                }
            }
        );
        
        this.context.subscriptions.push(disposable);
    }
    
    private async runValidator(
        content: string, 
        languageId: string
    ): Promise<ValidationResult> {
        // For now, return basic validation
        // TODO: Call Python validator for full validation
        const result: ValidationResult = {
            valid: true,
            errors: [],
            warnings: [],
            suggestions: []
        };
        
        // Basic checks
        if (!content.trim()) {
            result.valid = false;
            result.errors.push('File is empty');
        }
        
        // Language-specific checks
        switch (languageId) {
            case 'gaussian':
                if (!content.includes('#')) {
                    result.warnings.push('No route section found');
                }
                break;
            case 'vasp':
                const lines = content.split('\n');
                if (lines.length < 8) {
                    result.warnings.push('POSCAR may be incomplete');
                }
                break;
        }
        
        return result;
    }
    
    private showValidationResult(result: ValidationResult): void {
        if (result.valid && result.errors.length === 0 && result.warnings.length === 0) {
            vscode.window.showInformationMessage('✅ Input file is valid!');
        } else {
            const messages: string[] = [];
            
            if (result.errors.length > 0) {
                messages.push('❌ Errors:\n' + result.errors.map(e => `  • ${e}`).join('\n'));
            }
            
            if (result.warnings.length > 0) {
                messages.push('⚠️ Warnings:\n' + result.warnings.map(w => `  • ${w}`).join('\n'));
            }
            
            vscode.window.showWarningMessage(
                result.valid ? 'Validation completed with warnings' : 'Validation failed',
                { modal: true, detail: messages.join('\n\n') },
                'View Details'
            ).then(selection => {
                if (selection === 'View Details') {
                    this.outputChannel.show();
                }
            });
        }
    }
}
