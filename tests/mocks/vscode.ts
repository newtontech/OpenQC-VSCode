// Mock vscode module for tests
export const version = '1.120.0';

export const window = {
  showInformationMessage: jest.fn(() => Promise.resolve(undefined)),
  showWarningMessage: jest.fn(() => Promise.resolve(undefined)),
  showErrorMessage: jest.fn(() => Promise.resolve(undefined)),
  showTextDocument: jest.fn(() => Promise.resolve(undefined)),
  activeTextEditor: undefined,
  visibleTextEditors: [],
  tabGroups: {
    onDidChangeTabs: jest.fn(() => ({ dispose: jest.fn() })),
  },
  onDidChangeActiveTextEditor: jest.fn(() => ({ dispose: jest.fn() })),
  onDidChangeVisibleTextEditors: jest.fn(() => ({ dispose: jest.fn() })),
  onDidChangeTextEditorSelection: jest.fn(() => ({ dispose: jest.fn() })),
  onDidChangeTextEditorVisibleRanges: jest.fn(() => ({ dispose: jest.fn() })),
  createWebviewPanel: jest.fn(() => ({
    webview: {
      html: '',
      onDidReceiveMessage: jest.fn(),
      postMessage: jest.fn(),
      asWebviewUri: jest.fn((uri: any) => uri),
      cspSource: 'vscode-resource:',
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
    replace: jest.fn(),
    name: 'mock',
    error: jest.fn(),
    warn: jest.fn(),
    info: jest.fn(),
    log: jest.fn(),
    trace: jest.fn(),
    logLevel: 3, // LogLevel.Info
    onDidChangeLogLevel: jest.fn(() => ({ dispose: jest.fn() })),
  })),
};

export const workspace = {
  getConfiguration: jest.fn((section?: string) => {
    if (section === 'openqc.logging') {
      return {
        get: jest.fn((key: string, defaultValue?: any) => {
          if (key === 'level') return 'info';
          if (key === 'showUserMessages') return true;
          return defaultValue;
        }),
        update: jest.fn(),
        has: jest.fn(() => true),
      };
    }
    return {
      get: jest.fn((key: string, defaultValue?: any) => defaultValue),
      update: jest.fn(),
      has: jest.fn(() => true),
    };
  }),
  createFileSystemWatcher: jest.fn(() => ({
    onDidCreate: jest.fn(),
    onDidChange: jest.fn(),
    onDidDelete: jest.fn(),
    dispose: jest.fn(),
  })),
  onDidOpenTextDocument: jest.fn(() => ({ dispose: jest.fn() })),
  onDidChangeTextDocument: jest.fn(() => ({ dispose: jest.fn() })),
  onDidCloseTextDocument: jest.fn(() => ({ dispose: jest.fn() })),
  onDidSaveTextDocument: jest.fn(() => ({ dispose: jest.fn() })),
  onWillSaveTextDocument: jest.fn(() => ({ dispose: jest.fn() })),
  onDidChangeConfiguration: jest.fn(() => ({ dispose: jest.fn() })),
  onDidChangeWorkspaceFolders: jest.fn(() => ({ dispose: jest.fn() })),
  textDocuments: [],
  workspaceFolders: [],
  getWorkspaceFolder: jest.fn(),
  asRelativePath: jest.fn((path: string) => path),
  fs: {
    readFile: jest.fn(),
    writeFile: jest.fn(),
    stat: jest.fn(),
  },
};

export const languages = {
  registerCompletionItemProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerHoverProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerDefinitionProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerReferenceProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerDocumentHighlightProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerDocumentSymbolProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerCodeActionsProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerCodeLensProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerDocumentFormattingEditProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerDocumentRangeFormattingEditProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerOnTypeFormattingEditProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerRenameProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerDocumentLinkProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerColorProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerFoldingRangeProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerImplementationProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerSelectionRangeProvider: jest.fn(() => ({ dispose: jest.fn() })),
  registerTypeDefinitionProvider: jest.fn(() => ({ dispose: jest.fn() })),
  createDiagnosticCollection: jest.fn(() => ({
    set: jest.fn(),
    delete: jest.fn(),
    clear: jest.fn(),
    dispose: jest.fn(),
  })),
  match: jest.fn(() => 10),
};

export class Position {
  constructor(
    public line: number,
    public character: number
  ) {}
}

export class Range {
  constructor(
    public start: Position,
    public end: Position
  ) {}
}

export class CompletionItem {
  constructor(public label: string) {}
}

export class CodeLens {
  constructor(public range: Range) {}
}

export class CodeAction {
  constructor(
    public title: string,
    public kind?: unknown
  ) {}
}

export class Diagnostic {
  constructor(
    public range: Range,
    public message: string,
    public severity?: unknown
  ) {}
}

export class DocumentLink {
  constructor(
    public range: Range,
    public target?: unknown
  ) {}
}

export class InlayHint {
  constructor(
    public position: Position,
    public label: unknown,
    public kind?: unknown
  ) {}
}

export class SymbolInformation {
  constructor(
    public name: string,
    public kind: unknown,
    public rangeOrContainer: unknown,
    public locationOrUri?: unknown,
    public containerName?: string
  ) {}
}

export class CallHierarchyItem {
  constructor(
    public kind: unknown,
    public name: string,
    public detail: string,
    public uri: unknown,
    public range: Range,
    public selectionRange: Range
  ) {}
}

export class TypeHierarchyItem {
  constructor(
    public kind: unknown,
    public name: string,
    public detail: string,
    public uri: unknown,
    public range: Range,
    public selectionRange: Range
  ) {}
}

export class Location {
  constructor(
    public uri: unknown,
    public range: Range
  ) {}
}

export class DocumentSymbol {
  children: DocumentSymbol[] = [];
  tags?: unknown[];

  constructor(
    public name: string,
    public detail: string,
    public kind: unknown,
    public range: Range,
    public selectionRange: Range
  ) {}
}

export class CompletionList {
  constructor(
    public items: unknown[],
    public isIncomplete?: boolean
  ) {}
}

export class MarkdownString {
  value = '';
  supportHtml = false;

  constructor(value = '') {
    this.value = value;
  }

  appendMarkdown(value: string): MarkdownString {
    this.value += value;
    return this;
  }

  appendCodeblock(value: string): MarkdownString {
    this.value += value;
    return this;
  }
}

export class SnippetString {
  constructor(public value = '') {}
}

export class TextEdit {
  constructor(
    public range: Range,
    public newText: string
  ) {}
}

export class CancellationError extends Error {
  constructor() {
    super('Canceled');
    this.name = 'CancellationError';
  }
}

export const DiagnosticSeverity = {
  Error: 0,
  Warning: 1,
  Information: 2,
  Hint: 3,
};

export const LogLevel = {
  Off: 0,
  Trace: 1,
  Debug: 2,
  Info: 3,
  Warning: 4,
  Error: 5,
};

export class Disposable {
  static from(...disposables: { dispose: () => unknown }[]): Disposable {
    return new Disposable(() => disposables.forEach(d => d.dispose()));
  }
  constructor(private callOnDispose: () => unknown) {}
  dispose(): any {
    this.callOnDispose();
  }
}

export class CancellationTokenSource {
  token = {
    isCancellationRequested: false,
    onCancellationRequested: jest.fn(() => ({ dispose: jest.fn() })),
  };
  cancel(): void {
    this.token.isCancellationRequested = true;
  }
  dispose(): void {}
}

export class TabInputText {
  constructor(public uri: unknown) {}
}
export class TabInputTextDiff {
  constructor(
    public original: unknown,
    public modified: unknown
  ) {}
}
export class TabInputNotebook {
  constructor(public uri: unknown) {}
}
export class TabInputCustom {
  constructor(public uri: unknown) {} // eslint-disable-line @typescript-eslint/no-useless-constructor
}

export const CodeActionKind = {
  Empty: { value: '' },
  QuickFix: { value: 'quickfix' },
  Refactor: { value: 'refactor' },
  RefactorExtract: { value: 'refactor.extract' },
  RefactorInline: { value: 'refactor.inline' },
  RefactorRewrite: { value: 'refactor.rewrite' },
  Source: { value: 'source' },
  SourceOrganizeImports: { value: 'source.organizeImports' },
};

export const CompletionItemKind = {
  Text: 1,
};

export const CompletionItemTag = {
  Deprecated: 1,
};

export const DiagnosticTag = {
  Unnecessary: 1,
  Deprecated: 2,
};

export const DocumentHighlightKind = {
  Text: 1,
  Read: 2,
  Write: 3,
};

export const SymbolKind = {
  File: 1,
  Module: 2,
  Namespace: 3,
  Package: 4,
  Class: 5,
  Method: 6,
  Property: 7,
  Field: 8,
  Constructor: 9,
  Enum: 10,
  Interface: 11,
  Function: 12,
  Variable: 13,
  Constant: 14,
  String: 15,
  Number: 16,
  Boolean: 17,
  Array: 18,
  Object: 19,
  Key: 20,
  Null: 21,
  EnumMember: 22,
  Struct: 23,
  Event: 24,
  Operator: 25,
  TypeParameter: 26,
};

export const SymbolTag = {
  Deprecated: 1,
};

export const FoldingRangeKind = {
  Comment: 1,
  Imports: 2,
  Region: 3,
};

export const commands = {
  registerCommand: jest.fn(() => ({ dispose: jest.fn() })),
};

export const env = {
  remoteName: undefined,
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
  joinPath: jest.fn((base: any, ...paths: string[]) => {
    const basePath = base.path || base.fsPath || base.toString();
    const joinedPath = [basePath, ...paths].join('/').replace(/\/+/g, '/');
    return {
      path: joinedPath,
      fsPath: joinedPath,
      scheme: base.scheme || 'file',
      toString: () => joinedPath,
    };
  }),
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
  version,
  window,
  workspace,
  commands,
  env,
  languages,
  Position,
  Range,
  CompletionItem,
  CodeLens,
  CodeAction,
  Diagnostic,
  DocumentLink,
  InlayHint,
  SymbolInformation,
  CallHierarchyItem,
  TypeHierarchyItem,
  Location,
  DocumentSymbol,
  CompletionList,
  MarkdownString,
  SnippetString,
  TextEdit,
  CancellationError,
  CancellationTokenSource,
  DiagnosticSeverity,
  LogLevel,
  Disposable,
  CodeActionKind,
  CompletionItemKind,
  CompletionItemTag,
  DiagnosticTag,
  DocumentHighlightKind,
  SymbolKind,
  SymbolTag,
  FoldingRangeKind,
  ViewColumn,
  Uri,
  EventEmitter,
  ExtensionContext,
  TextDocument,
  TextEditor,
  TabInputText,
  TabInputTextDiff,
  TabInputNotebook,
  TabInputCustom,
};
