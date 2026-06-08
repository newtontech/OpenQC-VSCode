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
import { listBundledLspServers } from '../lsp/registry';
import {
  resolveLspCommand,
  readCommandOverrides,
  isExecutableAvailable,
  ResolvedLspCommand,
} from '../lsp/commandResolver';
import { LSPServerRegistryEntry } from '../lsp/types';
import {
  describeExtensionHostContext,
  getExtensionHostContext,
} from '../utils/extensionHostContext';

interface LSPServerConfig {
  name: string;
  enabled: boolean;
  resolvedCommand: ResolvedLspCommand;
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
  private serverDefinitions: LSPServerDefinition[] = LSPManager.bundledDefinitions();
  private logger = createComponentLogger('LSPManager');

  constructor(context?: vscode.ExtensionContext, discovery?: LSPDiscovery) {
    this.fileTypeDetector = new FileTypeDetector();
    this.config = vscode.workspace.getConfiguration('openqc.lsp');
    this.discovery = discovery || new LSPDiscovery(context);
  }

  /**
   * Build `LSPServerDefinition[]` from the bundled registry.
   *
   * This is the default source of server definitions used during normal
   * document-open flows. No network calls are required.
   */
  private static bundledDefinitions(): LSPServerDefinition[] {
    return listBundledLspServers().map(entry => ({
      id: entry.id,
      name: entry.name,
      repository: entry.repository,
      executable: entry.executable,
      languageId: entry.languageId,
      fileExtensions: [...entry.fileExtensions],
      fileNames: [...entry.fileNames],
      enabled: entry.enabled,
      repositoryUrl: entry.repositoryUrl,
    }));
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
    const hostContextDescription = describeExtensionHostContext();
    try {
      const client = await this.createLanguageClient(software, serverConfig, document, identity);
      if (client) {
        await client.start();
        this.clients.set(identity.key, {
          client,
          documents: new Set([documentKey]),
          languageId: serverConfig.definition.languageId,
        });
        this.logger.info(`${software} Language Server started in ${hostContextDescription}`);
        vscode.window.showInformationMessage(`${software} Language Server started`);
      }
    } catch (error) {
      const message = `Failed to start ${software} Language Server in ${hostContextDescription}: ${error}`;
      this.logger.error(message, error as Error);
      vscode.window.showErrorMessage(message);
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
      const resolved = config.resolvedCommand;
      const hostContext = getExtensionHostContext();
      const hostContextDescription = describeExtensionHostContext(hostContext);
      const commandSummary =
        resolved.kind === 'pythonModule'
          ? `${resolved.python} -m ${resolved.module}`
          : resolved.command;

      // Verify the executable is available using a cross-platform check.
      // For pythonModule kind, check the python binary.
      const executableToCheck =
        resolved.kind === 'pythonModule' ? resolved.python : resolved.command;

      this.logger.info(
        `Starting ${software} Language Server with '${commandSummary}' in ${hostContextDescription}`
      );

      const available = await isExecutableAvailable(executableToCheck);
      if (!available) {
        const settingKey = `openqc.lsp.${config.definition.languageId}.command`;
        throw new Error(
          `LSP executable '${executableToCheck}' not found in ${hostContext.label}. ` +
            `Install ${software} LSP server or update ${settingKey} for that environment. ` +
            hostContext.commandResolution
        );
      }

      let serverCommand: string;
      let serverArgs: string[];
      const serverEnv: Record<string, string> | undefined = resolved.env;

      if (resolved.kind === 'pythonModule') {
        serverCommand = resolved.python;
        serverArgs = ['-m', resolved.module, ...resolved.args];
      } else {
        serverCommand = resolved.command;
        serverArgs = resolved.args;
      }

      const serverOptions: ServerOptions = {
        command: serverCommand,
        args: serverArgs,
        transport: TransportKind.stdio,
        ...(serverEnv ? { options: { env: { ...process.env, ...serverEnv } } } : {}),
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

    const languageId = definition.languageId;
    const overrides = readCommandOverrides(this.config, languageId);

    // Find the matching registry entry for command resolution
    const registryEntry = this.findRegistryEntry(languageId);
    const resolvedCommand = registryEntry
      ? resolveLspCommand(registryEntry, overrides)
      : {
          kind: 'pathOrCommand' as const,
          command: overrides.command || overrides.path || definition.executable,
          args: overrides.args ?? ['--stdio'],
          env: overrides.env,
        };

    return {
      name: `${definition.name} LSP`,
      enabled: this.config.get<boolean>(`${languageId}.enabled`, definition.enabled),
      resolvedCommand,
      definition,
    };
  }

  private findRegistryEntry(languageId: string): LSPServerRegistryEntry | undefined {
    const allEntries = listBundledLspServers();
    return allEntries.find(entry => entry.languageId === languageId);
  }

  private async refreshDiscoveredDefinitions(): Promise<void> {
    this.serverDefinitions = LSPManager.bundledDefinitions();
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
   * Explicitly discover available LSPs from GitHub.
   * This is an opt-in operation and does NOT run during normal document-open flows.
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
