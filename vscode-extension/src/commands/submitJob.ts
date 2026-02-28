/**
 * Submit Job Command - Submit jobs to remote compute servers
 */

import * as vscode from 'vscode';

interface JobConfig {
    inputFile: string;
    jobName: string;
    queue: string;
    nodes: number;
    cpusPerNode: number;
    wallTime: string;
    memory?: string;
}

export class SubmitJobCommand {
    constructor(
        private context: vscode.ExtensionContext,
        private outputChannel: vscode.OutputChannel
    ) {}
    
    register() {
        const disposable = vscode.commands.registerCommand(
            'openqc.submitJob',
            async () => {
                // Check if connected to a server
                const connectedServer = this.context.globalState.get('openqc.connectedServer');
                
                if (!connectedServer) {
                    const connect = await vscode.window.showInformationMessage(
                        'Not connected to any server. Would you like to connect?',
                        'Yes', 'No'
                    );
                    
                    if (connect === 'Yes') {
                        await vscode.commands.executeCommand('openqc.connectServer');
                        return;
                    }
                    return;
                }
                
                const editor = vscode.window.activeTextEditor;
                if (!editor) {
                    vscode.window.showErrorMessage('No active editor');
                    return;
                }
                
                const inputFile = editor.document.uri.fsPath;
                
                // Gather job configuration
                const jobConfig = await this.getJobConfig(inputFile);
                
                if (!jobConfig) {
                    return;
                }
                
                // Submit job
                await this.submitJob(jobConfig);
            }
        );
        
        this.context.subscriptions.push(disposable);
    }
    
    private async getJobConfig(inputFile: string): Promise<JobConfig | undefined> {
        const jobName = await vscode.window.showInputBox({
            prompt: 'Enter job name',
            value: inputFile.split('/').pop()?.replace(/\.[^.]+$/, '') || 'qc_job',
            placeHolder: 'my_qc_job'
        });
        
        if (!jobName) {
            return undefined;
        }
        
        const queue = await vscode.window.showQuickPick(
            ['debug', 'short', 'normal', 'long', 'gpu'],
            { placeHolder: 'Select queue/partition' }
        );
        
        if (!queue) {
            return undefined;
        }
        
        const nodesStr = await vscode.window.showInputBox({
            prompt: 'Number of nodes',
            value: '1',
            placeHolder: '1'
        });
        
        const nodes = parseInt(nodesStr || '1', 10);
        
        const cpusStr = await vscode.window.showInputBox({
            prompt: 'CPUs per node',
            value: '32',
            placeHolder: '32'
        });
        
        const cpusPerNode = parseInt(cpusStr || '32', 10);
        
        const wallTime = await vscode.window.showInputBox({
            prompt: 'Wall time (HH:MM:SS)',
            value: '24:00:00',
            placeHolder: '24:00:00'
        });
        
        const memory = await vscode.window.showInputBox({
            prompt: 'Memory per node (optional)',
            placeHolder: 'e.g., 64GB'
        });
        
        return {
            inputFile,
            jobName,
            queue,
            nodes,
            cpusPerNode,
            wallTime: wallTime || '24:00:00',
            memory
        };
    }
    
    private async submitJob(config: JobConfig): Promise<void> {
        this.outputChannel.appendLine(`Submitting job: ${config.jobName}`);
        this.outputChannel.appendLine(`  Queue: ${config.queue}`);
        this.outputChannel.appendLine(`  Nodes: ${config.nodes}`);
        this.outputChannel.appendLine(`  CPUs/node: ${config.cpusPerNode}`);
        this.outputChannel.appendLine(`  Wall time: ${config.wallTime}`);
        
        // TODO: Generate and submit Slurm/PBS script
        const slurmScript = this.generateSlurmScript(config);
        
        this.outputChannel.appendLine('\nGenerated Slurm script:');
        this.outputChannel.appendLine(slurmScript);
        
        const action = await vscode.window.showInformationMessage(
            `Job "${config.jobName}" ready to submit`,
            'Submit', 'View Script', 'Cancel'
        );
        
        if (action === 'View Script') {
            const doc = await vscode.workspace.openTextDocument({
                content: slurmScript,
                language: 'shellscript'
            });
            vscode.window.showTextDocument(doc);
        } else if (action === 'Submit') {
            vscode.window.showInformationMessage(
                `Submitting job to queue... ` +
                `(Job submission feature coming soon)`
            );
        }
    }
    
    private generateSlurmScript(config: JobConfig): string {
        return `#!/bin/bash
#SBATCH --job-name=${config.jobName}
#SBATCH --partition=${config.queue}
#SBATCH --nodes=${config.nodes}
#SBATCH --ntasks-per-node=${config.cpusPerNode}
#SBATCH --time=${config.wallTime}
${config.memory ? `#SBATCH --mem=${config.memory}` : '#SBATCH --mem=64GB'}
#SBATCH --output=${config.jobName}_%j.out
#SBATCH --error=${config.jobName}_%j.err

# Load necessary modules
# module load gaussian/16
# module load vasp/6.3

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Run the calculation
# mpirun -np $SLURM_NTASKS vasp_std
# g16 < ${config.inputFile}

echo "Job completed at $(date)"
`;
    }
}
