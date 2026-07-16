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
  let messageHandler: ((message: any) => void) | undefined;
  const panel = {
    title: '',
    webview: {
      html: '',
      postMessage: jest.fn((msg: any) => {
        messages.push(msg);
        return true;
      }),
      onDidReceiveMessage: jest.fn((handler: (message: any) => void) => {
        messageHandler = handler;
        return { dispose: jest.fn() };
      }),
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
      showWarningMessage: jest.fn(),
      showErrorMessage: jest.fn(),
      showInformationMessage: jest.fn(),
      showSaveDialog: jest.fn(),
    },
    workspace: {
      fs: {
        writeFile: jest.fn(),
      },
    },
    commands: {
      executeCommand: jest.fn(),
    },
    ViewColumn: { One: 1, Two: 2 },
    Uri: {
      joinPath: jest.fn((...parts: string[]) => ({ fsPath: parts.join('/') })),
      file: jest.fn((p: string) => ({ fsPath: p })),
    },
    Disposable: { from: jest.fn() },
    _messages: messages,
    _panel: panel,
    _getMessageHandler: () => messageHandler,
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

    it('updates the source filename when reusing an existing viewer panel', async () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, '/tmp/a.xyz');
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, '/tmp/b.xyz');
      vscode.window.showWarningMessage.mockResolvedValue('Save');

      vscode._getMessageHandler()({
        type: 'structureUpdated',
        structure: JSON.stringify({
          ...WATER_STRUCTURE,
          atoms: [{ ...WATER_STRUCTURE.atoms[0], x: 7 }],
        }),
      });

      const saved = await OpenQCViewerPanel.saveCurrentStructureToSource();

      expect(saved).toBe(true);
      expect(vscode.workspace.fs.writeFile).toHaveBeenCalledWith(
        expect.objectContaining({ fsPath: '/tmp/b.xyz' }),
        expect.any(Buffer)
      );
      expect(vscode.workspace.fs.writeFile).not.toHaveBeenCalledWith(
        expect.objectContaining({ fsPath: '/tmp/a.xyz' }),
        expect.any(Buffer)
      );
    });

    it('sends only the canonical DTO load message to the production viewer', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      const messages = vscode._messages;
      const loadMsg = messages.find((m: any) => m.type === 'loadStructure');
      const initMsg = messages.find((m: any) => m.type === 'initialize');

      expect(loadMsg).toBeDefined();
      expect(loadMsg.structure).toBeDefined();
      expect(loadMsg.filename).toBe('water.xyz');
      expect(loadMsg.validation).toBeDefined();

      expect(initMsg).toBeUndefined();
    });

    it('sends periodic structure with cell data', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, SILICON_STRUCTURE, 'POSCAR');

      const loadMsg = vscode._messages.find((m: any) => m.type === 'loadStructure');
      const parsedStructure = JSON.parse(loadMsg.structure);
      expect(parsedStructure.cell).toBeDefined();
      expect(parsedStructure.cell.a).toEqual([5.43, 0.0, 0.0]);
      expect(vscode._messages.find((m: any) => m.type === 'initialize')).toBeUndefined();
    });

    it('updates current export data from webview-side structure edits', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');
      const edited = {
        ...WATER_STRUCTURE,
        atoms: [
          WATER_STRUCTURE.atoms[0],
          { ...WATER_STRUCTURE.atoms[1], x: 2.5 },
          WATER_STRUCTURE.atoms[2],
        ],
        bonds: [{ from: 0, to: 1, order: 2 }],
      };

      vscode._getMessageHandler()({
        type: 'structureUpdated',
        structure: JSON.stringify(edited),
      });

      expect(OpenQCViewerPanel.currentPanel?.isDirty()).toBe(true);
      expect(vscode._panel.title).toBe('OpenQC: • water.xyz');
      expect(OpenQCViewerPanel.getCurrentStructureData()?.atoms[1].x).toBe(2.5);
      expect(OpenQCViewerPanel.getCurrentStructureData()?.bonds).toEqual([
        { from: 0, to: 1, order: 2 },
      ]);
    });

    it('updates structure and invokes export picker on webview export request', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      vscode._getMessageHandler()({
        type: 'exportEditedStructure',
        structure: WATER_STRUCTURE,
      });

      expect(vscode.commands.executeCommand).toHaveBeenCalledWith(
        'openqc.exportStructureWithPicker'
      );
    });

    it('does not mark the panel dirty when exporting a clean webview payload', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      vscode._getMessageHandler()({
        type: 'exportEditedStructure',
        structure: JSON.stringify(WATER_STRUCTURE),
        dirty: false,
      });

      expect(OpenQCViewerPanel.currentPanel?.isDirty()).toBe(false);
      expect(vscode._panel.title).toBe('OpenQC: water.xyz');
      expect(vscode.commands.executeCommand).toHaveBeenCalledWith(
        'openqc.exportStructureWithPicker'
      );
    });

    it('rejects invalid edited structures without replacing current structure', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      vscode._getMessageHandler()({
        type: 'structureEdited',
        structure: { ...WATER_STRUCTURE, atoms: [] },
      });

      expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
        expect.stringContaining('edited structure is invalid')
      );
      expect(OpenQCViewerPanel.getCurrentStructureData()?.atoms).toHaveLength(3);
    });

    it('rejects edited structures with invalid bonds without replacing current structure', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      vscode._getMessageHandler()({
        type: 'structureUpdated',
        structure: JSON.stringify({
          ...WATER_STRUCTURE,
          bonds: [{ from: 0, to: 99, order: 1 }],
        }),
      });

      expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
        expect.stringContaining('edited structure is invalid')
      );
      expect(OpenQCViewerPanel.getCurrentStructureData()?.bonds).toBeUndefined();
    });

    it('does not invoke export picker when edited export payload is invalid', () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, 'water.xyz');

      vscode._getMessageHandler()({
        type: 'exportEditedStructure',
        structure: { ...WATER_STRUCTURE, atoms: [] },
      });

      expect(vscode.commands.executeCommand).not.toHaveBeenCalled();
    });

    it('saves edited structures back to native writable source files after confirmation', async () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, '/tmp/water.xyz');
      const edited = {
        ...WATER_STRUCTURE,
        atoms: [
          WATER_STRUCTURE.atoms[0],
          { ...WATER_STRUCTURE.atoms[1], x: 2.5 },
          WATER_STRUCTURE.atoms[2],
        ],
      };
      vscode.window.showWarningMessage.mockResolvedValue('Save');

      vscode._getMessageHandler()({
        type: 'structureUpdated',
        structure: JSON.stringify(edited),
      });

      const saved = await OpenQCViewerPanel.saveCurrentStructureToSource();

      expect(saved).toBe(true);
      expect(vscode.workspace.fs.writeFile).toHaveBeenCalledWith(
        expect.objectContaining({ fsPath: '/tmp/water.xyz' }),
        expect.any(Buffer)
      );
      const savedBuffer = vscode.workspace.fs.writeFile.mock.calls[0][1] as Buffer;
      expect(savedBuffer.toString()).toContain('H 2.500000 0.757200 -0.469200');
      expect(vscode._panel.webview.postMessage).toHaveBeenCalledWith(
        expect.objectContaining({ type: 'markStructureSaved' })
      );
      expect(OpenQCViewerPanel.currentPanel?.isDirty()).toBe(false);
      expect(vscode._panel.title).toBe('OpenQC: /tmp/water.xyz');
    });

    it('writes edited bonds back to PDB CONECT records for native PDB sources', async () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, '/tmp/water.pdb');
      vscode.window.showWarningMessage.mockResolvedValue('Save');

      vscode._getMessageHandler()({
        type: 'structureUpdated',
        structure: JSON.stringify({
          ...WATER_STRUCTURE,
          bonds: [{ from: 0, to: 1, order: 1 }],
        }),
      });

      const saved = await OpenQCViewerPanel.saveCurrentStructureToSource();

      expect(saved).toBe(true);
      const savedBuffer = vscode.workspace.fs.writeFile.mock.calls[0][1] as Buffer;
      expect(savedBuffer.toString()).toContain('CONECT    1    2');
    });

    it('rejects source writeback for non-native source formats', async () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, '/tmp/cp2k.inp');

      vscode._getMessageHandler()({
        type: 'structureUpdated',
        structure: JSON.stringify({
          ...WATER_STRUCTURE,
          atoms: [{ ...WATER_STRUCTURE.atoms[0], x: 1 }],
        }),
      });

      const saved = await OpenQCViewerPanel.saveCurrentStructureToSource();

      expect(saved).toBe(false);
      expect(vscode.window.showErrorMessage).toHaveBeenCalledWith(
        expect.stringContaining('only write edited structures back')
      );
      expect(vscode.workspace.fs.writeFile).not.toHaveBeenCalled();
    });

    it('handles webview save-source requests using the latest edited payload', async () => {
      OpenQCViewerPanel.createOrShow({ fsPath: '/ext' } as any, WATER_STRUCTURE, '/tmp/water.xyz');
      vscode.window.showWarningMessage.mockResolvedValue('Save');

      vscode._getMessageHandler()({
        type: 'saveEditedStructureToSource',
        structure: JSON.stringify({
          ...WATER_STRUCTURE,
          atoms: [{ ...WATER_STRUCTURE.atoms[0], x: 4 }],
        }),
      });

      await new Promise(resolve => setImmediate(resolve));

      const savedBuffer = vscode.workspace.fs.writeFile.mock.calls[0][1] as Buffer;
      expect(savedBuffer.toString()).toContain('O 4.000000 0.000000 0.117300');
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
