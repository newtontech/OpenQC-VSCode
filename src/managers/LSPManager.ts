import * as vscode from 'vscode';
import {
  LanguageClient,
  LanguageClientOptions,
  ServerOptions,
  TransportKind,
} from 'vscode-languageclient/node';
import { FileTypeDetector, QuantumChemistrySoftware } from './FileTypeDetector';
import { LSPDiscovery, LSPServerDefinition } from '../utils/LSPDiscovery';

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
  private lspDiscovery: LSPDiscovery;

  constructor(context?: vscode.ExtensionContext) {
    this.fileTypeDetector = new FileTypeDetector();
    this.config = vscode.workspace.getConfiguration('openqc.lsp');
    this.lspDiscovery = new LSPDiscovery(context);
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

    try {
      const client = await this.createLanguageClient(software, serverConfig, document);
      if (client) {
        this.clients.set(languageId, client);
        await client.start();
        vscode.window.showInformationMessage(`${software} Language Server started`);
      }
    } catch (error) {
      console.error(`Error starting ${software} Language Server:`, error);
      vscode.window.showErrorMessage(`Failed to start ${software} Language Server: ${error}`);
      // Clean up the client if it was added
      this.clients.delete(languageId);
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
      try {
        // Check if client is running before stopping
        if (client.needsStop()) {
          await client.stop();
        }
      } catch (error) {
        console.error(`Error stopping ${software} Language Server:`, error);
        vscode.window.showWarningMessage(`Error stopping ${software} Language Server: ${error}`);
      } finally {
        // Always remove from clients map to allow restart
        this.clients.delete(languageId);
      }
    }
  }

  async restartLSPForDocument(document: vscode.TextDocument): Promise<void> {
    const software = this.fileTypeDetector.detectSoftware(document);
    if (!software) {
      vscode.window.showWarningMessage('Could not detect quantum chemistry software for this file');
      return;
    }

    try {
      await this.stopLSPForDocument(document);
      // Add a small delay to ensure the process is fully terminated
      await new Promise(resolve => setTimeout(resolve, 500));
      await this.startLSPForDocument(document);
    } catch (error) {
      console.error(`Error restarting ${software} Language Server:`, error);
      vscode.window.showErrorMessage(`Failed to restart ${software} Language Server: ${error}`);
    }
  }

  private async createLanguageClient(
    software: QuantumChemistrySoftware,
    config: LSPServerConfig,
    document: vscode.TextDocument
  ): Promise<LanguageClient | undefined> {
    try {
      // Verify the LSP executable exists
      const { exec } = require('child_process');
      const { promisify } = require('util');
      const execAsync = promisify(exec);

      try {
        await execAsync(`which ${config.path}`);
      } catch {
        throw new Error(
          `LSP executable '${config.path}' not found in PATH. Please install ${software} LSP server.`
        );
      }

      const serverOptions: ServerOptions = {
        command: config.path,
        args: config.args,
        transport: TransportKind.stdio,
      };

      const clientOptions: LanguageClientOptions = {
        documentSelector: [{ scheme: 'file', language: this.getLanguageId(software) }],
        synchronize: {
          fileEvents: vscode.workspace.createFileSystemWatcher(
            `**/*.{${this.getExtensions(software).join(',')}}`
          ),
        },
      };

      return new LanguageClient(
        `openqc-${software.toLowerCase()}`,
        `OpenQC ${software} Language Server`,
        serverOptions,
        clientOptions
      );
    } catch (error) {
      console.error(`Failed to create LanguageClient for ${software}:`, error);
      throw error;
    }
  }

  private getServerConfig(software: QuantumChemistrySoftware): LSPServerConfig {
    const softwareKey = software.toLowerCase();
    return {
      name: `${software} LSP`,
      enabled: this.config.get<boolean>(`${softwareKey}.enabled`, true),
      path: this.config.get<string>(`${softwareKey}.path`, `${softwareKey}-lsp`),
      args: ['--stdio'],
    };
  }

  private getLanguageId(software: QuantumChemistrySoftware): string {
    const mapping: Record<QuantumChemistrySoftware, string> = {
      CP2K: 'cp2k',
      VASP: 'vasp',
      Gaussian: 'gaussian',
      ORCA: 'orca',
      'Quantum ESPRESSO': 'qe',
      GAMESS: 'gamess',
      NWChem: 'nwchem',
    };
    return mapping[software];
  }

  private getExtensions(software: QuantumChemistrySoftware): string[] {
    const mapping: Record<QuantumChemistrySoftware, string[]> = {
      CP2K: ['inp'],
      VASP: ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR'],
      Gaussian: ['gjf', 'com'],
      ORCA: ['inp'],
      'Quantum ESPRESSO': ['in', 'pw.in', 'relax.in'],
      GAMESS: ['inp'],
      NWChem: ['nw'],
    };
    return mapping[software];
  }

  /**
   * Dynamically discover available LSPs from OpenQuantumChemistry GitHub organization
   */
  async discoverAvailableLSPs(): Promise<LSPServerDefinition[]> {
    try {
      return await this.lspDiscovery.fetchLSPRepositories();
    } catch (error) {
      console.error('[LSPManager] Failed to discover LSPs:', error);
      return [];
    }
  }

  /**
   * Force refresh LSP list from GitHub
   */
  async refreshLSPList(): Promise<void> {
    try {
      await this.lspDiscovery.fetchLSPRepositories(true);
      vscode.window.showInformationMessage('LSP list refreshed successfully');
    } catch (error) {
      vscode.window.showErrorMessage(`Failed to refresh LSP list: ${error}`);
    }
  }

  dispose(): void {
    const stopPromises = Array.from(this.clients.entries()).map(async ([languageId, client]) => {
      try {
        if (client.needsStop()) {
          await client.stop();
        }
      } catch (error) {
        console.error(`Error stopping ${languageId} Language Server during dispose:`, error);
      }
    });

    Promise.all(stopPromises).then(() => {
      this.clients.clear();
    });
  }
}
