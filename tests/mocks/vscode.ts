// Mock vscode module for tests
export const window = {
  showInformationMessage: jest.fn(),
  showWarningMessage: jest.fn(),
  showErrorMessage: jest.fn(),
  activeTextEditor: undefined,
  createWebviewPanel: jest.fn(() => ({
    webview: {
      html: '',
      onDidReceiveMessage: jest.fn(),
      postMessage: jest.fn(),
      asWebviewUri: jest.fn((uri: any) => uri),
    },
    onDidDispose: jest.fn((cb: () => void) => ({ dispose: jest.fn() })),
    reveal: jest.fn(),
    dispose: jest.fn(),
    visible: true,
  })),
  createOutputChannel: jest.fn(() => ({
    appendLine: jest.fn(),
    append: jest.fn(),
    clear: jest.fn(),
    show: jest.fn(),
    hide: jest.fn(),
    dispose: jest.fn(),
  })),
};

export const workspace = {
  getConfiguration: jest.fn(() => ({
    get: jest.fn((key: string, defaultValue?: any) => defaultValue),
    update: jest.fn(),
    has: jest.fn(() => true),
  })),
  createFileSystemWatcher: jest.fn(() => ({
    onDidCreate: jest.fn(),
    onDidChange: jest.fn(),
    onDidDelete: jest.fn(),
    dispose: jest.fn(),
  })),
  onDidOpenTextDocument: jest.fn(() => ({ dispose: jest.fn() })),
  onDidCloseTextDocument: jest.fn(() => ({ dispose: jest.fn() })),
  workspaceFolders: [],
  getWorkspaceFolder: jest.fn(),
  asRelativePath: jest.fn((path: string) => path),
  fs: {
    readFile: jest.fn(),
    writeFile: jest.fn(),
    stat: jest.fn(),
  },
};

export const commands = {
  registerCommand: jest.fn(() => ({ dispose: jest.fn() })),
};

export const ViewColumn = {
  One: 1,
  Two: 2,
  Three: 3,
};

export const Uri = {
  file: jest.fn((path: string) => ({
    path,
    fsPath: path,
    scheme: 'file',
    toString: () => path,
  })),
  parse: jest.fn((uri: string) => ({
    path: uri,
    fsPath: uri,
    scheme: 'file',
    toString: () => uri,
  })),
};

export const EventEmitter = jest.fn().mockImplementation(() => ({
  event: jest.fn(),
  fire: jest.fn(),
  dispose: jest.fn(),
}));

export const ExtensionContext = jest.fn();
export const TextDocument = jest.fn();
export const TextEditor = jest.fn();

export default {
  window,
  workspace,
  commands,
  ViewColumn,
  Uri,
  EventEmitter,
  ExtensionContext,
  TextDocument,
  TextEditor,
};
