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
  definition: LSPServerDefinition;
}

export class LSPManager {
  private clients: Map<string, LanguageClient> = new Map();
  private fileTypeDetector: FileTypeDetector;
  private config: vscode.WorkspaceConfiguration;
  private discovery: LSPDiscovery;
  private serverDefinitions: LSPServerDefinition[] = LSPDiscovery.getDefaultDefinitions();

  constructor(context?: vscode.ExtensionContext, discovery?: LSPDiscovery) {
    this.fileTypeDetector = new FileTypeDetector();
    this.config = vscode.workspace.getConfiguration('openqc.lsp');
    this.discovery = discovery || new LSPDiscovery(context);
  }

  async startLSPForDocument(document: vscode.TextDocument): Promise<void> {
    const software = this.fileTypeDetector.detectSoftware(document);
    if (!software) {
      return;
    }

    const serverConfig = await this.getServerConfig(software);
    if (!serverConfig) {
      return;
    }

    const languageId = serverConfig.definition.languageId;
    if (this.clients.has(languageId)) {
      return; // LSP already running
    }

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

    const definition = this.findServerDefinition(software);
    if (!definition) {
      return;
    }

    const languageId = definition.languageId;
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
        documentSelector: [{ scheme: 'file', language: config.definition.languageId }],
        synchronize: {
          fileEvents: vscode.workspace.createFileSystemWatcher(
            this.createWatcherPattern(config.definition)
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

  private async getServerConfig(
    software: QuantumChemistrySoftware
  ): Promise<LSPServerConfig | undefined> {
    await this.refreshDiscoveredDefinitions();

    const definition = this.findServerDefinition(software);
    if (!definition) {
      return undefined;
    }

    const softwareKey = definition.languageId;
    return {
      name: `${definition.name} LSP`,
      enabled: this.config.get<boolean>(`${softwareKey}.enabled`, definition.enabled),
      path: this.config.get<string>(`${softwareKey}.path`, definition.executable),
      args: ['--stdio'],
      definition,
    };
  }

  private async refreshDiscoveredDefinitions(): Promise<void> {
    this.serverDefinitions = await this.discovery.fetchLSPRepositories();
  }

  private findServerDefinition(
    software: QuantumChemistrySoftware
  ): LSPServerDefinition | undefined {
    return this.serverDefinitions.find(definition => definition.name === software);
  }

  private createWatcherPattern(definition: LSPServerDefinition): string {
    const filePatterns = [
      ...definition.fileExtensions.map(extension => `*.${extension}`),
      ...(definition.fileNames || []),
    ];

    if (filePatterns.length === 0) {
      return '**/*';
    }

    if (filePatterns.length === 1) {
      return `**/${filePatterns[0]}`;
    }

    return `**/{${filePatterns.join(',')}}`;
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
