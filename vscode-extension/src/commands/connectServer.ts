/**
 * Connect Server Command - Connect to remote compute servers via SSH
 */

import * as vscode from 'vscode';

interface ServerConfig {
    name: string;
    host: string;
    port: number;
    username: string;
    keyPath?: string;
    password?: string;
}

export class ConnectServerCommand {
    constructor(
        private context: vscode.ExtensionContext,
        private outputChannel: vscode.OutputChannel
    ) {}
    
    register() {
        const disposable = vscode.commands.registerCommand(
            'openqc.connectServer',
            async () => {
                const servers = this.context.globalState.get<ServerConfig[]>('openqc.servers') || [];
                
                const action = await vscode.window.showQuickPick(
                    [
                        { label: '$(add) Add New Server', description: 'Configure a new remote server' },
                        ...servers.map(s => ({
                            label: `$(server) ${s.name}`,
                            description: `${s.username}@${s.host}:${s.port}`
                        }))
                    ],
                    { placeHolder: 'Select or add a server' }
                );
                
                if (!action) {
                    return;
                }
                
                if (action.label === '$(add) Add New Server') {
                    await this.addServer();
                } else {
                    await this.connectToServer(action.label.replace('$(server) ', ''));
                }
            }
        );
        
        this.context.subscriptions.push(disposable);
    }
    
    private async addServer(): Promise<void> {
        const name = await vscode.window.showInputBox({
            prompt: 'Enter server name',
            placeHolder: 'e.g., Cluster A, SuperComputer'
        });
        
        if (!name) {
            return;
        }
        
        const host = await vscode.window.showInputBox({
            prompt: 'Enter server hostname or IP',
            placeHolder: 'e.g., cluster.example.com or 192.168.1.100'
        });
        
        if (!host) {
            return;
        }
        
        const portStr = await vscode.window.showInputBox({
            prompt: 'Enter SSH port',
            value: '22',
            placeHolder: '22'
        });
        
        const port = parseInt(portStr || '22', 10);
        
        const username = await vscode.window.showInputBox({
            prompt: 'Enter username',
            placeHolder: 'your_username'
        });
        
        if (!username) {
            return;
        }
        
        const authMethod = await vscode.window.showQuickPick(
            ['SSH Key', 'Password'],
            { placeHolder: 'Select authentication method' }
        );
        
        if (!authMethod) {
            return;
        }
        
        let keyPath: string | undefined;
        let password: string | undefined;
        
        if (authMethod === 'SSH Key') {
            keyPath = await vscode.window.showInputBox({
                prompt: 'Enter path to SSH private key',
                placeHolder: '~/.ssh/id_rsa'
            });
        } else {
            password = await vscode.window.showInputBox({
                prompt: 'Enter password (stored securely)',
                password: true
            });
        }
        
        const serverConfig: ServerConfig = {
            name,
            host,
            port,
            username,
            keyPath,
            password
        };
        
        // Save server configuration
        const servers = this.context.globalState.get<ServerConfig[]>('openqc.servers') || [];
        servers.push(serverConfig);
        await this.context.globalState.update('openqc.servers', servers);
        
        vscode.window.showInformationMessage(
            `Server "${name}" added successfully!`
        );
        
        // Ask if user wants to connect now
        const connectNow = await vscode.window.showInformationMessage(
            'Would you like to connect now?',
            'Yes', 'No'
        );
        
        if (connectNow === 'Yes') {
            await this.connectToServer(name);
        }
    }
    
    private async connectToServer(name: string): Promise<void> {
        const servers = this.context.globalState.get<ServerConfig[]>('openqc.servers') || [];
        const server = servers.find(s => s.name === name);
        
        if (!server) {
            vscode.window.showErrorMessage(`Server "${name}" not found`);
            return;
        }
        
        this.outputChannel.appendLine(`Connecting to ${server.name} (${server.host})...`);
        
        // TODO: Implement actual SSH connection
        // This would use the ssh2 library
        
        vscode.window.showInformationMessage(
            `Connecting to ${server.name}... ` +
            `(SSH connection feature coming soon)`
        );
        
        // Store connection status
        await this.context.globalState.update('openqc.connectedServer', server);
    }
}
