#!/usr/bin/env node

const fs = require('fs');

let buffer = Buffer.alloc(0);
let nextVersion = 0;
const documents = new Map();
const eventsPath = process.env.FAKE_LSP_EVENTS_PATH;

function logEvent(event) {
  if (!eventsPath) {
    return;
  }

  fs.appendFileSync(
    eventsPath,
    `${JSON.stringify({ ...event, sequence: ++nextVersion })}\n`,
    'utf8'
  );
}

function writeMessage(message) {
  const body = JSON.stringify(message);
  process.stdout.write(`Content-Length: ${Buffer.byteLength(body, 'utf8')}\r\n\r\n${body}`);
}

function sendResponse(id, result) {
  writeMessage({ jsonrpc: '2.0', id, result });
}

function sendNotification(method, params) {
  writeMessage({ jsonrpc: '2.0', method, params });
}

function textHasProblem(text) {
  return text.toLowerCase().includes('problem');
}

function publishDiagnostics(uri) {
  const text = documents.get(uri) || '';
  const diagnostics = textHasProblem(text)
    ? [
        {
          range: {
            start: { line: 0, character: 0 },
            end: { line: 0, character: 7 },
          },
          severity: 2,
          source: 'fake-lsp',
          message: 'Fake diagnostic for problem marker',
        },
      ]
    : [];

  sendNotification('textDocument/publishDiagnostics', { uri, diagnostics });
  logEvent({ method: 'textDocument/publishDiagnostics', uri, diagnosticCount: diagnostics.length });
}

function handleMessage(message) {
  const { id, method, params } = message;
  logEvent({ method, id });

  if (method === 'initialize') {
    sendResponse(id, {
      capabilities: {
        textDocumentSync: {
          openClose: true,
          change: 1,
        },
        completionProvider: {
          triggerCharacters: ['#'],
        },
      },
    });
    return;
  }

  if (method === 'initialized') {
    return;
  }

  if (method === 'textDocument/didOpen') {
    const document = params.textDocument;
    documents.set(document.uri, document.text);
    logEvent({
      method: 'fake/documentState',
      uri: document.uri,
      text: document.text,
      version: document.version,
    });
    publishDiagnostics(document.uri);
    return;
  }

  if (method === 'textDocument/didChange') {
    const uri = params.textDocument.uri;
    const text = params.contentChanges[params.contentChanges.length - 1].text;
    documents.set(uri, text);
    logEvent({
      method: 'fake/documentState',
      uri,
      text,
      version: params.textDocument.version,
    });
    publishDiagnostics(uri);
    return;
  }

  if (method === 'textDocument/completion') {
    sendResponse(id, {
      isIncomplete: false,
      items: [
        {
          label: 'fake-lsp-completion',
          kind: 14,
          detail: 'Fake LSP completion item',
        },
      ],
    });
    return;
  }

  if (method === 'shutdown') {
    sendResponse(id, null);
    return;
  }

  if (method === 'exit') {
    process.exit(0);
  }

  if (id !== undefined) {
    writeMessage({
      jsonrpc: '2.0',
      id,
      error: {
        code: -32601,
        message: `Method not found: ${method}`,
      },
    });
  }
}

function readMessages() {
  while (true) {
    const headerEnd = buffer.indexOf('\r\n\r\n');
    if (headerEnd === -1) {
      return;
    }

    const header = buffer.subarray(0, headerEnd).toString('utf8');
    const lengthMatch = /Content-Length: (\d+)/i.exec(header);
    if (!lengthMatch) {
      throw new Error(`Missing Content-Length header: ${header}`);
    }

    const length = Number(lengthMatch[1]);
    const bodyStart = headerEnd + 4;
    const bodyEnd = bodyStart + length;
    if (buffer.length < bodyEnd) {
      return;
    }

    const body = buffer.subarray(bodyStart, bodyEnd).toString('utf8');
    buffer = buffer.subarray(bodyEnd);
    handleMessage(JSON.parse(body));
  }
}

process.stdin.on('data', chunk => {
  buffer = Buffer.concat([buffer, chunk]);
  readMessages();
});

process.stdin.on('end', () => process.exit(0));
