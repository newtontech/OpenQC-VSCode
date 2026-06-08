import * as vscode from 'vscode';
import {
  LanguageClient,
  LanguageClientOptions,
  ServerOptions,
  TransportKind,
} from 'vscode-languageclient/node';
import { FileTypeDetector, QuantumChemistrySoftware } from './FileTypeDetector';
import { LSPDiscovery, LSPServerDefinition } from '../utils/LSPDiscovery';
import { createComponentLogger } from '../utils/Logger';

interface LSPServerConfig {
  name: string;
  enabled: boolean;
  path: string;
  args: string[];
  definition: LSPServerDefinition;
}

interface LSPClientIdentity {
  key: string;
  label: string;
  workspaceFolder?: vscode.WorkspaceFolder;
}

interface LSPClientRecord {
  client: LanguageClient;
  documents: Set<string>;
  languageId: string;
}

export class LSPManager {
  private clients: Map<string, LSPClientRecord> = new Map();
  private startingClients: Map<string, Promise<void>> = new Map();
  private fileTypeDetector: FileTypeDetector;
  private config: vscode.WorkspaceConfiguration;
  private discovery: LSPDiscovery;
  private serverDefinitions: LSPServerDefinition[] = LSPDiscovery.getDefaultDefinitions();
  private logger = createComponentLogger('LSPManager');

  constructor(context?: vscode.ExtensionContext, discovery?: LSPDiscovery) {
    this.fileTypeDetector = new FileTypeDetector();
    this.config = vscode.workspace.getConfiguration('openqc.lsp');
    this.discovery = discovery || new LSPDiscovery(context);
  }

  /**
   * Start the matching LSP client for a quantum chemistry document.
   *
   * Detects the document's software format, resolves the configured language server,
   * and reuses an existing workspace-scoped client when one is already running.
   *
   * @param document - VS Code document that should be attached to an LSP client.
   */
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
    const identity = this.getClientIdentity(languageId, document);
    const documentKey = this.getDocumentKey(document);
    const existingRecord = this.clients.get(identity.key);
    if (existingRecord) {
      existingRecord.documents.add(documentKey);
      return; // LSP already running
    }

    if (!serverConfig.enabled) {
      return;
    }

    const pendingStart = this.startingClients.get(identity.key);
    if (pendingStart) {
      await pendingStart;
      const startedRecord = this.clients.get(identity.key);
      if (startedRecord) {
        startedRecord.documents.add(documentKey);
      }
      return;
    }

    const startPromise = this.startClient(software, serverConfig, document, identity, documentKey);
    this.startingClients.set(identity.key, startPromise);
    try {
      await startPromise;
    } finally {
      this.startingClients.delete(identity.key);
    }
  }

  private async startClient(
    software: QuantumChemistrySoftware,
    serverConfig: LSPServerConfig,
    document: vscode.TextDocument,
    identity: LSPClientIdentity,
    documentKey: string
  ): Promise<void> {
    try {
      const client = await this.createLanguageClient(software, serverConfig, document, identity);
      if (client) {
        await client.start();
        this.clients.set(identity.key, {
          client,
          documents: new Set([documentKey]),
          languageId: serverConfig.definition.languageId,
        });
        this.logger.info(`${software} Language Server started`);
        vscode.window.showInformationMessage(`${software} Language Server started`);
      }
    } catch (error) {
      this.logger.error(`Failed to start ${software} Language Server`, error as Error);
      vscode.window.showErrorMessage(`Failed to start ${software} Language Server: ${error}`);
      // Clean up the client if it was added
      this.clients.delete(identity.key);
    }
  }

  /**
   * Stop tracking a document and shut down its LSP client when no documents remain.
   *
   * @param document - VS Code document that should be detached from its LSP client.
   */
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
    const identity = this.getClientIdentity(languageId, document);
    const record = this.clients.get(identity.key);
    if (!record) {
      return;
    }

    record.documents.delete(this.getDocumentKey(document));
    if (record.documents.size > 0) {
      return;
    }

    try {
      await this.stopClientRecord(identity.key, record);
    } catch (error) {
      this.logger.error(`Error stopping ${software} Language Server`, error as Error);
      vscode.window.showWarningMessage(`Error stopping ${software} Language Server: ${error}`);
    }
  }

  /**
   * Restart the LSP client associated with a quantum chemistry document.
   *
   * @param document - VS Code document whose detected language server should be restarted.
   */
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
      this.logger.error(`Error restarting ${software} Language Server`, error as Error);
      vscode.window.showErrorMessage(`Failed to restart ${software} Language Server: ${error}`);
    }
  }

  private async createLanguageClient(
    software: QuantumChemistrySoftware,
    config: LSPServerConfig,
    document: vscode.TextDocument,
    identity: LSPClientIdentity
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

      const watcherPattern = this.createWatcherPattern(config.definition);
      const fileEvents = identity.workspaceFolder
        ? vscode.workspace.createFileSystemWatcher(
            new vscode.RelativePattern(identity.workspaceFolder, watcherPattern)
          )
        : vscode.workspace.createFileSystemWatcher(watcherPattern);

      const documentSelector = identity.workspaceFolder
        ? [
            {
              scheme: document.uri.scheme || 'file',
              language: config.definition.languageId,
              pattern: this.createWorkspaceDocumentPattern(
                identity.workspaceFolder,
                watcherPattern
              ),
            },
          ]
        : [{ scheme: document.uri.scheme || 'file', language: config.definition.languageId }];

      const clientOptions: LanguageClientOptions = {
        documentSelector,
        workspaceFolder: identity.workspaceFolder,
        synchronize: {
          fileEvents,
        },
      };

      return new LanguageClient(
        `openqc-${software.toLowerCase()}-${this.sanitizeClientId(identity.key)}`,
        `OpenQC ${software} Language Server${identity.label ? ` (${identity.label})` : ''}`,
        serverOptions,
        clientOptions
      );
    } catch (error) {
      this.logger.error(`Failed to create LanguageClient for ${software}`, error as Error);
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

  private createWorkspaceDocumentPattern(
    workspaceFolder: vscode.WorkspaceFolder,
    watcherPattern: string
  ): string {
    const workspacePath = workspaceFolder.uri.fsPath.replace(/\\/g, '/').replace(/\/+$/, '');
    return `${workspacePath}/${watcherPattern}`;
  }

  private getClientIdentity(languageId: string, document: vscode.TextDocument): LSPClientIdentity {
    const workspaceFolder = vscode.workspace.getWorkspaceFolder(document.uri);
    if (workspaceFolder) {
      return {
        key: `${languageId}:workspace:${this.uriToString(workspaceFolder.uri)}`,
        label: workspaceFolder.name,
        workspaceFolder,
      };
    }

    return {
      key: `${languageId}:document:${this.getDocumentKey(document)}`,
      label: '',
    };
  }

  private getDocumentKey(document: vscode.TextDocument): string {
    return this.uriToString(document.uri) || document.fileName;
  }

  private uriToString(uri: vscode.Uri): string {
    return typeof uri.toString === 'function' ? uri.toString() : uri.fsPath;
  }

  private sanitizeClientId(value: string): string {
    return value.replace(/[^a-zA-Z0-9_-]+/g, '-').replace(/^-+|-+$/g, '');
  }

  private async stopClientRecord(clientKey: string, record: LSPClientRecord): Promise<void> {
    try {
      if (record.client.needsStop()) {
        await record.client.stop();
      }
    } finally {
      this.clients.delete(clientKey);
    }
  }

  /**
   * Dynamically discover available LSPs from OpenQuantumChemistry GitHub organization
   */
  async discoverAvailableLSPs(): Promise<LSPServerDefinition[]> {
    try {
      return await this.discovery.fetchLSPRepositories();
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
      await this.discovery.fetchLSPRepositories(true);
      vscode.window.showInformationMessage('LSP list refreshed successfully');
    } catch (error) {
      vscode.window.showErrorMessage(`Failed to refresh LSP list: ${error}`);
    }
  }

  /**
   * Stop all managed LSP clients and clear their records.
   */
  async dispose(): Promise<void> {
    const stopPromises = Array.from(this.clients.entries()).map(async ([clientKey, record]) => {
      try {
        await this.stopClientRecord(clientKey, record);
      } catch (error) {
        this.logger.error(
          `Error stopping ${record.languageId} Language Server during dispose`,
          error as Error
        );
      }
    });

    await Promise.all(stopPromises);
  }
}
