import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {
  CompletionRequest,
  DidChangeTextDocumentNotification,
  DidOpenTextDocumentNotification,
  LanguageClient,
  PublishDiagnosticsNotification,
  TransportKind,
} from 'vscode-languageclient/node';

interface FakeLspEvent {
  method: string;
  id?: number | string;
  uri?: string;
  text?: string;
  version?: number;
  diagnosticCount?: number;
  sequence: number;
}

const fixturePath = path.join(__dirname, '../fixtures/fake-lsp-server.js');

function readEvents(eventsPath: string): FakeLspEvent[] {
  if (!fs.existsSync(eventsPath)) {
    return [];
  }

  return fs
    .readFileSync(eventsPath, 'utf8')
    .split('\n')
    .filter(Boolean)
    .map(line => JSON.parse(line) as FakeLspEvent);
}

async function waitForEvent(
  eventsPath: string,
  predicate: (event: FakeLspEvent) => boolean
): Promise<FakeLspEvent> {
  const timeoutAt = Date.now() + 5000;

  while (Date.now() < timeoutAt) {
    const match = readEvents(eventsPath).find(predicate);
    if (match) {
      return match;
    }
    await new Promise(resolve => setTimeout(resolve, 25));
  }

  throw new Error(`Timed out waiting for fake LSP event in ${eventsPath}`);
}

describe('fake stdio LSP client integration', () => {
  let client: LanguageClient | undefined;
  let tempDir: string;
  let eventsPath: string;

  beforeEach(() => {
    tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-fake-lsp-'));
    eventsPath = path.join(tempDir, 'events.ndjson');
  });

  afterEach(async () => {
    if (client?.needsStop()) {
      await client.stop();
    }
    client = undefined;
    fs.rmSync(tempDir, { recursive: true, force: true });
  });

  it('starts, exchanges document lifecycle messages, completes, and shuts down', async () => {
    const diagnostics: unknown[] = [];
    client = new LanguageClient(
      'openqc-fake-lsp-integration',
      'OpenQC Fake LSP Integration',
      {
        command: process.execPath,
        args: [fixturePath],
        transport: TransportKind.stdio,
        options: {
          env: {
            ...process.env,
            FAKE_LSP_EVENTS_PATH: eventsPath,
          },
        },
      },
      {
        documentSelector: [{ scheme: 'file', language: 'gaussian' }],
        initializationOptions: {
          testName: 'fakeLspClient',
        },
      }
    );

    client.onNotification(PublishDiagnosticsNotification.type, params => {
      diagnostics.push(params);
    });

    await client.start();
    await waitForEvent(eventsPath, event => event.method === 'initialize');
    await waitForEvent(eventsPath, event => event.method === 'initialized');

    const uri = 'file:///tmp/openqc-fake-lsp-test.gjf';
    client.sendNotification(DidOpenTextDocumentNotification.type, {
      textDocument: {
        uri,
        languageId: 'gaussian',
        version: 1,
        text: '# hf/sto-3g\n\n0 1\nH 0 0 0\n',
      },
    });

    await waitForEvent(
      eventsPath,
      event =>
        event.method === 'fake/documentState' &&
        event.uri === uri &&
        event.version === 1 &&
        event.text?.includes('hf/sto-3g') === true
    );

    client.sendNotification(DidChangeTextDocumentNotification.type, {
      textDocument: { uri, version: 2 },
      contentChanges: [{ text: '# problem\n\n0 1\nH 0 0 0\n' }],
    });

    await waitForEvent(
      eventsPath,
      event =>
        event.method === 'fake/documentState' &&
        event.uri === uri &&
        event.version === 2 &&
        event.text?.includes('problem') === true
    );
    await waitForEvent(
      eventsPath,
      event =>
        event.method === 'textDocument/publishDiagnostics' &&
        event.uri === uri &&
        event.diagnosticCount === 1
    );

    const completion = await client.sendRequest(CompletionRequest.type, {
      textDocument: { uri },
      position: { line: 0, character: 1 },
      context: { triggerKind: 1 },
    });

    expect(completion).toEqual(
      expect.objectContaining({
        isIncomplete: false,
        items: [
          expect.objectContaining({
            label: 'fake-lsp-completion',
            detail: 'Fake LSP completion item',
          }),
        ],
      })
    );
    expect(diagnostics).toEqual([
      expect.objectContaining({
        uri,
        diagnostics: [],
      }),
      expect.objectContaining({
        uri,
        diagnostics: [expect.objectContaining({ source: 'fake-lsp' })],
      }),
    ]);

    await client.stop();
    await waitForEvent(eventsPath, event => event.method === 'shutdown');
    await waitForEvent(eventsPath, event => event.method === 'exit');
  });
});
