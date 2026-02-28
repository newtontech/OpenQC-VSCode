/**
 * Convert Command - Convert between quantum chemistry file formats
 */

import * as vscode from 'vscode';
import { spawn } from 'child_process';

export class ConvertCommand {
    private formats = [
        'Gaussian (.gjf)',
        'VASP POSCAR',
        'Quantum ESPRESSO (.pw)',
        'ORCA (.inp)',
        'XYZ',
        'CIF',
        'PDB',
        'mol2',
        'ASE JSON'
    ];
    
    constructor(
        private context: vscode.ExtensionContext,
        private outputChannel: vscode.OutputChannel
    ) {}
    
    register() {
        const disposable = vscode.commands.registerCommand(
            'openqc.convert',
            async () => {
                const editor = vscode.window.activeTextEditor;
                if (!editor) {
                    vscode.window.showErrorMessage('No active editor');
                    return;
                }
                
                // Ask for target format
                const targetFormat = await vscode.window.showQuickPick(
                    this.formats,
                    { placeHolder: 'Select target format' }
                );
                
                if (!targetFormat) {
                    return;
                }
                
                const sourceFile = editor.document.uri.fsPath;
                const targetFile = await vscode.window.showInputBox({
                    prompt: 'Enter output file path',
                    value: sourceFile.replace(/\.[^.]+$/, this.getExtension(targetFormat))
                });
                
                if (!targetFile) {
                    return;
                }
                
                this.outputChannel.appendLine(
                    `Converting ${sourceFile} to ${targetFormat} -> ${targetFile}`
                );
                
                // Call Python converter
                await this.runConverter(sourceFile, targetFile, targetFormat);
            }
        );
        
        this.context.subscriptions.push(disposable);
    }
    
    private getExtension(format: string): string {
        const extMap: Record<string, string> = {
            'Gaussian (.gjf)': '.gjf',
            'VASP POSCAR': '.vasp',
            'Quantum ESPRESSO (.pw)': '.pw',
            'ORCA (.inp)': '.inp',
            'XYZ': '.xyz',
            'CIF': '.cif',
            'PDB': '.pdb',
            'mol2': '.mol2',
            'ASE JSON': '.json'
        };
        return extMap[format] || '.xyz';
    }
    
    private async runConverter(
        source: string, 
        target: string, 
        format: string
    ): Promise<void> {
        const config = vscode.workspace.getConfiguration('openqc');
        const pythonPath = config.get<string>('pythonPath') || 'python';
        
        return new Promise((resolve, reject) => {
            const proc = spawn(pythonPath, [
                '-m', 'openqc.converters',
                '--input', source,
                '--output', target,
                '--format', format
            ]);
            
            proc.on('close', (code) => {
                if (code === 0) {
                    vscode.window.showInformationMessage(
                        `Successfully converted to ${format}`
                    );
                    vscode.workspace.openTextDocument(target).then(doc => {
                        vscode.window.showTextDocument(doc);
                    });
                    resolve();
                } else {
                    vscode.window.showErrorMessage(
                        `Conversion failed with code ${code}`
                    );
                    reject(new Error(`Conversion failed: ${code}`));
                }
            });
            
            proc.on('error', (err) => {
                vscode.window.showErrorMessage(
                    `Failed to run converter: ${err.message}`
                );
                reject(err);
            });
        });
    }
}
