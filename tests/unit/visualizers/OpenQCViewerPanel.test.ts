/**
 * Unit tests for OpenQCViewerPanel (issue #106).
 * @module tests/unit/visualizers/OpenQCViewerPanel.test
 */

import { OpenQCViewerPanel } from '../../../src/visualizers/OpenQCViewerPanel';
import { validateOpenQCStructure } from '../../../src/structures/validation';
import { WATER_STRUCTURE, SILICON_STRUCTURE } from '../../../src/structures/fixtures';
import { createOpenQCStructure } from '../../../src/structures/converters';

// Mock vscode
jest.mock('vscode', () => {
  const messages: any[] = [];
  const panel = {
    webview: {
      html: '',
      postMessage: jest.fn((msg: any) => {
        messages.push(msg);
        return true;
      }),
      onDidReceiveMessage: jest.fn(),
      asWebviewUri: jest.fn((uri: any) => uri),
      cspSource: 'https://test-host',
    },
    reveal: jest.fn(),
    onDidDispose: jest.fn((cb: any) => {
      cb();
      return { dispose: jest.fn() };
    }),
    onDidChangeViewState: jest.fn(() => ({ dispose: jest.fn() })),
    dispose: jest.fn(),
  };

  return {
    window: {
      activeTextEditor: undefined,
      createWebviewPanel: jest.fn(() => panel),
      createOutputChannel: jest.fn(() => ({
        appendLine: jest.fn(),
        dispose: jest.fn(),
      })),
      showErrorMessage: jest.fn(),
      showInformationMessage: jest.fn(),
      showSaveDialog: jest.fn(),
    },
    ViewColumn: { One: 1, Two: 2 },
    Uri: {
      joinPath: jest.fn((...parts: string[]) => ({ fsPath: parts.join('/') })),
      file: jest.fn((p: string) => ({ fsPath: p })),
    },
    Disposable: { from: jest.fn() },
    _messages: messages,
    _panel: panel,
  };
});

const vscode = require('vscode');

describe('OpenQCViewerPanel', () => {
  beforeEach(() => {
    // Reset singleton between tests
    (OpenQCViewerPanel as any).currentPanel = undefined;
    vscode._messages.length = 0;
    jest.clearAllMocks();
  });

  describe('validate (static)', () => {
    it('accepts a valid molecule structure', () => {
      const result = OpenQCViewerPanel.validate(WATER_STRUCTURE);
      expect(result.valid).toBe(true);
      expect(result.errors).toHaveLength(0);
    });

    it('accepts a valid periodic structure', () => {
      const result = OpenQCViewerPanel.validate(SILICON_STRUCTURE);
      expect(result.valid).toBe(true);
    });

    it('rejects an empty atoms array', () => {
      const result = OpenQCViewerPanel.validate({
        ...WATER_STRUCTURE,
        atoms: [],
      });
      expect(result.valid).toBe(false);
      expect(result.errors).toContain('atoms array must not be empty');
    });

    it('rejects invalid structure', () => {
      const result = OpenQCViewerPanel.validate(null);
      expect(result.valid).toBe(false);
    });
  });

  describe('createOrShow', () => {
    it('creates a new panel when none exists', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      expect(vscode.window.createWebviewPanel).toHaveBeenCalledTimes(1);
      expect(vscode._panel.webview.postMessage).toHaveBeenCalled();
    });

    it('reveals existing panel when one already exists', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      // Second call should reveal, not create
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, SILICON_STRUCTURE, 'POSCAR');

      expect(vscode.window.createWebviewPanel).toHaveBeenCalledTimes(1);
      expect(vscode._panel.reveal).toHaveBeenCalled();
    });

    it('sends both loadStructure and initialize messages', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      const messages = vscode._messages;
      const loadMsg = messages.find((m: any) => m.type === 'loadStructure');
      const initMsg = messages.find((m: any) => m.type === 'initialize');

      expect(loadMsg).toBeDefined();
      expect(loadMsg.structure).toBeDefined();
      expect(loadMsg.filename).toBe('water.xyz');
      expect(loadMsg.validation).toBeDefined();

      expect(initMsg).toBeDefined();
      expect(initMsg.structure).toBeDefined();
      expect(initMsg.structure.xyz).toBeDefined();
    });

    it('sends periodic structure with cell data', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, SILICON_STRUCTURE, 'POSCAR');

      const loadMsg = vscode._messages.find((m: any) => m.type === 'loadStructure');
      const parsedStructure = JSON.parse(loadMsg.structure);
      expect(parsedStructure.cell).toBeDefined();
      expect(parsedStructure.cell.a).toEqual([5.43, 0.0, 0.0]);
    });
  });

  describe('createOpenQCStructure integration', () => {
    it('creates a valid molecule structure from atoms', () => {
      const atoms = [
        { element: 'C', x: 0, y: 0, z: 0 },
        { element: 'H', x: 0.6276, y: 0.6276, z: 0.6276 },
      ];

      const structure = createOpenQCStructure(atoms, {
        name: 'methane-partial',
        sourceSoftware: 'test',
      });

      expect(structure.atoms).toHaveLength(2);
      expect(structure.kind).toBe('molecule');
      expect(structure.name).toBe('methane-partial');
      expect(structure.metadata?.source?.software).toBe('test');

      const validation = validateOpenQCStructure(structure);
      expect(validation.valid).toBe(true);
    });

    it('creates a valid periodic structure with cell', () => {
      const atoms = [
        { element: 'Si', x: 0, y: 0, z: 0 },
        { element: 'Si', x: 0.25, y: 0.25, z: 0.25 },
      ];

      const structure = createOpenQCStructure(atoms, {
        name: 'Si',
        cell: {
          a: [5.43, 0, 0],
          b: [0, 5.43, 0],
          c: [0, 0, 5.43],
          pbc: [true, true, true],
        },
      });

      expect(structure.kind).toBe('periodic');
      expect(structure.cell).toBeDefined();
    });
  });
});
