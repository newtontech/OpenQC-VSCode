import * as vscode from 'vscode';
import { LanguageClient, LanguageClientOptions, ServerOptions, TransportKind } from 'vscode-languageclient/node';
import { FileTypeDetector, QuantumChemistrySoftware } from './FileTypeDetector';

interface LSPServerConfig {
    name: string;
    enabled: boolean;
    path: string;
    args: string[];
}

export class LSPManager {
    private clients: Map<string, LanguageClient> = new Map();
    private fileTypeDetector: FileTypeDetector;
    private config: vscode.WorkspaceConfiguration;

    constructor() {
        this.fileTypeDetector = new FileTypeDetector();
        this.config = vscode.workspace.getConfiguration('openqc.lsp');
    }

    async startLSPForDocument(document: vscode.TextDocument): Promise<void> {
        const software = this.fileTypeDetector.detectSoftware(document);
        if (!software) {
            return;
        }

        const languageId = this.getLanguageId(software);
        if (this.clients.has(languageId)) {
            return; // LSP already running
        }

        const serverConfig = this.getServerConfig(software);
        if (!serverConfig.enabled) {
            return;
        }

        const client = await this.createLanguageClient(software, serverConfig, document);
        if (client) {
            this.clients.set(languageId, client);
            await client.start();
            vscode.window.showInformationMessage(`${software} Language Server started`);
        }
    }

    async stopLSPForDocument(document: vscode.TextDocument): Promise<void> {
        const software = this.fileTypeDetector.detectSoftware(document);
        if (!software) {
            return;
        }

        const languageId = this.getLanguageId(software);
        const client = this.clients.get(languageId);
        if (client) {
            await client.stop();
            this.clients.delete(languageId);
        }
    }

    async restartLSPForDocument(document: vscode.TextDocument): Promise<void> {
        await this.stopLSPForDocument(document);
        await this.startLSPForDocument(document);
    }

    private async createLanguageClient(
        software: QuantumChemistrySoftware,
        config: LSPServerConfig,
        document: vscode.TextDocument
    ): Promise<LanguageClient | undefined> {
        const serverOptions: ServerOptions = {
            command: config.path,
            args: config.args,
            transport: TransportKind.stdio
        };

        const clientOptions: LanguageClientOptions = {
            documentSelector: [{ scheme: 'file', language: this.getLanguageId(software) }],
            synchronize: {
                fileEvents: vscode.workspace.createFileSystemWatcher(`**/*.{${this.getExtensions(software).join(',')}}`)
            }
        };

        return new LanguageClient(
            `openqc-${software.toLowerCase()}`,
            `OpenQC ${software} Language Server`,
            serverOptions,
            clientOptions
        );
    }

    private getServerConfig(software: QuantumChemistrySoftware): LSPServerConfig {
        const softwareKey = software.toLowerCase();
        return {
            name: `${software} LSP`,
            enabled: this.config.get<boolean>(`${softwareKey}.enabled`, true),
            path: this.config.get<string>(`${softwareKey}.path`, `${softwareKey}-lsp`),
            args: ['--stdio']
        };
    }

    private getLanguageId(software: QuantumChemistrySoftware): string {
        const mapping: Record<QuantumChemistrySoftware, string> = {
            'CP2K': 'cp2k',
            'VASP': 'vasp',
            'Gaussian': 'gaussian',
            'ORCA': 'orca',
            'Quantum ESPRESSO': 'qe',
            'GAMESS': 'gamess',
            'NWChem': 'nwchem'
        };
        return mapping[software];
    }

    private getExtensions(software: QuantumChemistrySoftware): string[] {
        const mapping: Record<QuantumChemistrySoftware, string[]> = {
            'CP2K': ['inp'],
            'VASP': ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR'],
            'Gaussian': ['gjf', 'com'],
            'ORCA': ['inp'],
            'Quantum ESPRESSO': ['in', 'pw.in', 'relax.in'],
            'GAMESS': ['inp'],
            'NWChem': ['nw']
        };
        return mapping[software];
    }

    dispose(): void {
        this.clients.forEach(async (client) => {
            await client.stop();
        });
        this.clients.clear();
    }
}