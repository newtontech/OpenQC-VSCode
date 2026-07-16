/**
 * Unit tests for OpenQCViewerWebview (issue #107).
 * @module tests/unit/webviews/openqcViewerWebview.test
 */

import { OpenQCViewerWebview } from '../../../src/webviews/openqcViewerWebview';
import * as fs from 'fs';
import * as path from 'path';
import * as vm from 'vm';

// Mock vscode
jest.mock('vscode', () => ({
  Uri: {
    joinPath: jest.fn((base: any, ...parts: string[]) => ({
      fsPath: `${base.fsPath}/${parts.join('/')}`,
      path: `${base.path}/${parts.join('/')}`,
    })),
  },
}));

function loadViewerTestApi(): any {
  const script = fs.readFileSync(
    path.resolve(__dirname, '../../../media/openqc-viewer.js'),
    'utf8'
  );
  const context: any = {
    window: { addEventListener: jest.fn() },
    document: {
      readyState: 'loading',
      addEventListener: jest.fn(),
      getElementById: jest.fn(),
    },
    acquireVsCodeApi: () => ({ postMessage: jest.fn() }),
    setInterval,
    clearInterval,
    setTimeout,
    JSON,
    parseInt,
  };

  vm.runInNewContext(script, context);
  return context.window.__openqcViewerTest;
}

describe('OpenQCViewerWebview', () => {
  const mockWebview = {
    asWebviewUri: jest.fn((uri: any) => `vscode-webview://${uri.path}`),
    cspSource: 'https://test-host',
  };

  const mockExtensionUri = { fsPath: '/ext', path: '/ext' };

  beforeEach(() => {
    jest.clearAllMocks();
  });

  describe('generateHTML', () => {
    it('generates valid HTML with bundled assets', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('<!DOCTYPE html>');
      expect(html).toContain('<html');
      expect(html).toContain('</html>');
    });

    it('includes CSP with nonce', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('Content-Security-Policy');
      expect(html).toContain("script-src 'nonce-");
      expect(html).toContain("style-src 'nonce-");
    });

    it('loads bundled viewer JS', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('openqc-viewer.js');
    });

    it('loads bundled CSS', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('openqc-viewer.css');
    });

    it('loads 3Dmol from bundled assets, not CDN', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('3Dmol-min.js');
      expect(html).not.toContain('cdnjs.cloudflare.com');
      expect(html).not.toContain('unpkg.com');
    });

    it('includes viewer-canvas div', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="viewer-canvas"');
    });

    it('includes toolbar with style selector', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="style-select"');
      expect(html).toContain('ball-stick');
      expect(html).toContain('spacefill');
    });

    it('includes selection, measurement, and export controls', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="btn-labels"');
      expect(html).toContain('id="btn-measure"');
      expect(html).toContain('id="btn-export"');
    });

    it('includes bounded atom editing controls', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="editing-controls"');
      expect(html).toContain('id="edit-element"');
      expect(html).toContain('id="btn-add-atom"');
      expect(html).toContain('id="btn-delete-atom"');
      expect(html).toContain('id="btn-move-atom"');
      expect(html).toContain('id="bond-from"');
      expect(html).toContain('id="bond-to"');
      expect(html).toContain('id="bond-order"');
      expect(html).toContain('id="btn-add-bond"');
      expect(html).toContain('id="btn-delete-bond"');
      expect(html).toContain('id="bond-count"');
      expect(html).toContain('id="constraint-x"');
      expect(html).toContain('id="constraint-y"');
      expect(html).toContain('id="constraint-z"');
      expect(html).toContain('id="btn-set-constraints"');
      expect(html).toContain('id="btn-undo-edit"');
      expect(html).toContain('id="btn-redo-edit"');
      expect(html).toContain('id="btn-save-source"');
      expect(html).toContain('id="btn-export-structure"');
      expect(html).toContain('id="dirty-indicator"');
      expect(html).toContain('>Clean</span>');
    });

    it('includes periodic controls', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="periodic-controls"');
      expect(html).toContain('id="sc-na"');
    });

    it('includes trajectory controls', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="trajectory-controls"');
      expect(html).toContain('id="frame-slider"');
    });

    it('includes WebGL fallback', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="webgl-fallback"');
    });

    it('includes status bar', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('id="status-bar"');
      expect(html).toContain('id="atom-count"');
    });

    it('loads 3Dmol from vendored media assets', () => {
      const html = OpenQCViewerWebview.generateHTML(mockWebview as any, mockExtensionUri as any);

      expect(html).toContain('media/vendor/3dmol/3Dmol-min.js');
    });
  });

  describe('getWebviewOptions', () => {
    it('enables scripts', () => {
      const options = OpenQCViewerWebview.getWebviewOptions(mockExtensionUri as any);
      expect(options.enableScripts).toBe(true);
    });

    it('sets local resource roots to media assets only', () => {
      const options = OpenQCViewerWebview.getWebviewOptions(mockExtensionUri as any);
      expect(options.localResourceRoots).toHaveLength(1);
      expect((options.localResourceRoots?.[0] as any).fsPath).toContain('media');
    });
  });

  describe('bundled viewer script', () => {
    it('keeps the 3Dmol canvas constrained below toolbar controls', () => {
      const css = fs.readFileSync(
        path.resolve(__dirname, '../../../media/openqc-viewer.css'),
        'utf8'
      );

      expect(css).toMatch(/#viewer-canvas\s*{[^}]*position:\s*relative;/s);
      expect(css).toMatch(/#viewer-canvas\s*{[^}]*overflow:\s*hidden;/s);
    });

    it('keeps deterministic edit helpers uniquely declared', () => {
      const script = fs.readFileSync(
        path.resolve(__dirname, '../../../media/openqc-viewer.js'),
        'utf8'
      );

      expect(script.match(/function normalizeSelectionIndex\(/g)).toHaveLength(1);
      expect(script.match(/function applyEditingAction\(/g)).toHaveLength(1);
      expect(script.match(/function calculateDihedral\(/g)).toHaveLength(1);
      expect(script.match(/function createEditedStructureExportMessage\(/g)).toHaveLength(1);
      expect(script.match(/function createEditedStructureSaveMessage\(/g)).toHaveLength(1);
      expect(script.match(/function createStructureUpdateMessage\(/g)).toHaveLength(1);
      expect(script.match(/function normalizeStructureBond\(/g)).toHaveLength(1);
      expect(script.match(/function normalizeBondPair\(/g)).toHaveLength(1);
      expect(script.match(/function normalizeSelectiveDynamics\(/g)).toHaveLength(1);
      expect(script.match(/function renderExplicitBonds\(/g)).toHaveLength(1);
      expect(script.match(/function replicateBondsForSupercell\(/g)).toHaveLength(1);
      expect(script.match(/function buildSupercellStructure\(/g)).toHaveLength(1);
      expect(script).toContain("type: 'structureUpdated'");
      expect(script).toContain("type: 'saveEditedStructureToSource'");
      expect(script).toContain("case 'markStructureSaved'");
      expect(script).toContain("getElementById('btn-save-source')");
      expect(script).toContain("getElementById('btn-add-bond')");
      expect(script).toContain("getElementById('btn-set-constraints')");
    });

    it('converts fractional periodic coordinates to Cartesian XYZ before 3Dmol rendering', () => {
      const api = loadViewerTestApi();

      const xyz = api.structureToXYZ({
        name: 'Si',
        kind: 'periodic',
        atoms: [
          { element: 'Si', x: 0, y: 0, z: 0 },
          { element: 'Si', x: 0.25, y: 0.25, z: 0.25 },
        ],
        cell: {
          a: [2.715, 2.715, 0],
          b: [0, 2.715, 2.715],
          c: [2.715, 0, 2.715],
          pbc: [true, true, true],
          coordinateMode: 'fractional',
        },
      });

      expect(xyz).toContain('Si 1.357500 1.357500 1.357500');
    });

    it('exposes deterministic distance, angle, and dihedral helpers for measurement mode', () => {
      const api = loadViewerTestApi();

      expect(api.calculateDistance({ x: 0, y: 0, z: 0 }, { x: 0, y: 3, z: 4 })).toBe(5);
      expect(
        api.calculateAngle({ x: 1, y: 0, z: 0 }, { x: 0, y: 0, z: 0 }, { x: 0, y: 1, z: 0 })
      ).toBeCloseTo(90);
      expect(
        api.calculateDihedral(
          { x: 1, y: 0, z: 0 },
          { x: 0, y: 0, z: 0 },
          { x: 0, y: 1, z: 0 },
          { x: 0, y: 1, z: 1 }
        )
      ).toBeCloseTo(90);
      expect(typeof api.handleAtomClick).toBe('function');
    });

    it('applies add, move, delete, undo, and redo as deterministic edit transitions', () => {
      const api = loadViewerTestApi();
      let state = api.createEditorState({
        name: 'H',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      });

      expect(state.dirty).toBe(false);

      state = api.applyEditingAction(state, {
        type: 'addAtom',
        atom: { element: 'cl', x: 1, y: 2, z: 3 },
      });
      expect(state.structure.atoms).toHaveLength(2);
      expect(state.structure.atoms[1]).toMatchObject({ element: 'Cl', x: 1, y: 2, z: 3 });
      expect(state.selectedIndex).toBe(1);
      expect(state.dirty).toBe(true);
      expect(state.undoStack).toHaveLength(1);

      state = api.applyEditingAction(state, { type: 'moveSelectedAtom', dx: 0.5, dy: -1, dz: 2 });
      expect(state.structure.atoms[1]).toMatchObject({ element: 'Cl', x: 1.5, y: 1, z: 5 });

      state = api.applyEditingAction(state, { type: 'deleteSelectedAtom' });
      expect(state.structure.atoms).toHaveLength(1);
      expect(state.selectedIndex).toBe(0);

      state = api.applyEditingAction(state, { type: 'undo' });
      expect(state.structure.atoms).toHaveLength(2);
      expect(state.structure.atoms[1]).toMatchObject({ element: 'Cl', x: 1.5, y: 1, z: 5 });

      state = api.applyEditingAction(state, { type: 'undo' });
      expect(state.structure.atoms[1]).toMatchObject({ element: 'Cl', x: 1, y: 2, z: 3 });

      state = api.applyEditingAction(state, { type: 'redo' });
      expect(state.structure.atoms[1]).toMatchObject({ element: 'Cl', x: 1.5, y: 1, z: 5 });
    });

    it('applies bond add, order update, delete, and atom-delete remapping deterministically', () => {
      const api = loadViewerTestApi();
      let state = api.createEditorState({
        name: 'linear',
        atoms: [
          { element: 'C', x: 0, y: 0, z: 0 },
          { element: 'O', x: 1, y: 0, z: 0 },
          { element: 'H', x: 2, y: 0, z: 0 },
        ],
      });

      state = api.applyEditingAction(state, {
        type: 'addOrUpdateBond',
        from: 0,
        to: 1,
        order: 2,
      });
      expect(state.structure.bonds).toEqual([{ from: 0, to: 1, order: 2 }]);
      expect(state.dirty).toBe(true);
      expect(state.undoStack).toHaveLength(1);

      state = api.applyEditingAction(state, {
        type: 'addOrUpdateBond',
        from: 1,
        to: 0,
        order: 3,
      });
      expect(state.structure.bonds).toEqual([{ from: 0, to: 1, order: 3 }]);

      state = api.applyEditingAction(state, { type: 'deleteBond', from: 1, to: 0 });
      expect(state.structure.bonds).toEqual([]);

      state = api.applyEditingAction(state, { type: 'undo' });
      expect(state.structure.bonds).toEqual([{ from: 0, to: 1, order: 3 }]);

      state = api.createEditorState(
        {
          name: 'remap',
          atoms: [
            { element: 'C', x: 0, y: 0, z: 0 },
            { element: 'O', x: 1, y: 0, z: 0 },
            { element: 'H', x: 2, y: 0, z: 0 },
          ],
          bonds: [
            { from: 0, to: 2, order: 2 },
            { from: 1, to: 2, order: 1 },
          ],
        },
        1
      );
      state = api.applyEditingAction(state, { type: 'deleteSelectedAtom' });
      expect(state.structure.atoms).toHaveLength(2);
      expect(state.structure.bonds).toEqual([{ from: 0, to: 1, order: 2 }]);
    });

    it('applies selected-atom selective dynamics constraints with undo/redo', () => {
      const api = loadViewerTestApi();
      let state = api.createEditorState(
        {
          name: 'slab',
          atoms: [
            { element: 'Fe', x: 0, y: 0, z: 0, selectiveDynamics: [true, true, true] },
            { element: 'Fe', x: 0, y: 0, z: 1 },
          ],
        },
        1
      );

      state = api.applyEditingAction(state, {
        type: 'setSelectedAtomConstraints',
        selectiveDynamics: [false, false, true],
      });

      expect(state.structure.atoms[0].selectiveDynamics).toEqual([true, true, true]);
      expect(state.structure.atoms[1].selectiveDynamics).toEqual([false, false, true]);
      expect(state.dirty).toBe(true);

      state = api.applyEditingAction(state, { type: 'undo' });
      expect(state.structure.atoms[1].selectiveDynamics).toBeUndefined();
      expect(state.dirty).toBe(false);

      state = api.applyEditingAction(state, { type: 'redo' });
      expect(state.structure.atoms[1].selectiveDynamics).toEqual([false, false, true]);
    });

    it('normalizes VASP-style selective dynamics flags without coercing invalid triples', () => {
      const api = loadViewerTestApi();

      expect(api.normalizeSelectiveDynamics([true, 'F', 't'])).toEqual([true, false, true]);
      expect(api.normalizeSelectiveDynamics(['T', 'false', 1])).toBeNull();
      expect(api.normalizeSelectiveDynamics([true, false])).toBeNull();
    });

    it('rejects invalid bond orders and replicates valid bonds for supercells', () => {
      const api = loadViewerTestApi();
      const structure = {
        atoms: [
          { element: 'C', x: 0, y: 0, z: 0 },
          { element: 'O', x: 1, y: 0, z: 0 },
        ],
      };

      expect(api.normalizeStructureBond({ from: 0, to: 1, order: 0 }, structure)).toBeNull();
      expect(api.normalizeStructureBond({ from: 0, to: 1, order: 4 }, structure)).toBeNull();
      expect(api.normalizeStructureBond({ from: 0, to: 1, order: 1.5 }, structure)).toBeNull();
      expect(api.normalizeStructureBond({ from: 0, to: 1, order: '2' }, structure)).toBeNull();

      expect(api.replicateBondsForSupercell([{ from: 0, to: 1, order: 2 }], 2, 3)).toEqual([
        { from: 0, to: 1, order: 2 },
        { from: 2, to: 3, order: 2 },
        { from: 4, to: 5, order: 2 },
      ]);
    });

    it('builds supercells without dropping atom metadata or edited bonds', () => {
      const api = loadViewerTestApi();
      const supercell = api.buildSupercellStructure(
        {
          name: 'constrained slab',
          kind: 'periodic',
          atoms: [
            {
              element: 'Fe',
              x: 0,
              y: 0,
              z: 0,
              selectiveDynamics: [false, false, true],
              tag: 'bottom',
            },
            {
              element: 'O',
              x: 0.5,
              y: 0.5,
              z: 0.5,
              selectiveDynamics: [true, true, true],
              charge: -1,
            },
          ],
          bonds: [{ from: 0, to: 1, order: 2 }],
          cell: {
            a: [2, 0, 0],
            b: [0, 2, 0],
            c: [0, 0, 4],
            pbc: [true, true, true],
            coordinateMode: 'fractional',
          },
        },
        2,
        1,
        1
      );

      expect(supercell.atoms).toHaveLength(4);
      expect(supercell.atoms[0]).toMatchObject({
        element: 'Fe',
        x: 0,
        y: 0,
        z: 0,
        selectiveDynamics: [false, false, true],
        tag: 'bottom',
      });
      expect(supercell.atoms[1]).toMatchObject({
        element: 'O',
        x: 1,
        y: 1,
        z: 2,
        selectiveDynamics: [true, true, true],
        charge: -1,
      });
      expect(supercell.atoms[2]).toMatchObject({
        element: 'Fe',
        x: 2,
        y: 0,
        z: 0,
        selectiveDynamics: [false, false, true],
        tag: 'bottom',
      });
      expect(supercell.cell).toMatchObject({
        a: [4, 0, 0],
        b: [0, 2, 0],
        c: [0, 0, 4],
        pbc: [true, true, true],
        coordinateMode: 'cartesian',
      });
      expect(supercell.bonds).toEqual([
        { from: 0, to: 1, order: 2 },
        { from: 2, to: 3, order: 2 },
      ]);
    });

    it('clears dirty state when undo returns to the loaded baseline', () => {
      const api = loadViewerTestApi();
      let state = api.createEditorState({
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      });

      state = api.applyEditingAction(state, {
        type: 'addAtom',
        atom: { element: 'O', x: 0, y: 0, z: 1 },
      });
      expect(state.dirty).toBe(true);

      state = api.applyEditingAction(state, { type: 'undo' });
      expect(state.structure.atoms).toHaveLength(1);
      expect(state.dirty).toBe(false);
      expect(state.redoStack).toHaveLength(1);
    });

    it('creates an edited-structure export message without extension-side dependencies', () => {
      const api = loadViewerTestApi();
      const message = api.createEditedStructureExportMessage(
        {
          name: 'Edited',
          atoms: [
            { element: 'h', x: 0, y: 0, z: 0 },
            { element: 'o', x: 1, y: 0, z: 0 },
          ],
          bonds: [{ from: 0, to: 1, order: 2 }],
        },
        true
      );

      expect(message).toMatchObject({
        type: 'exportEditedStructure',
        dirty: true,
        atomCount: 2,
        bondCount: 1,
      });
      expect(JSON.parse(message.structure)).toMatchObject({
        name: 'Edited',
        atoms: [
          { element: 'H', x: 0, y: 0, z: 0 },
          { element: 'O', x: 1, y: 0, z: 0 },
        ],
        bonds: [{ from: 0, to: 1, order: 2 }],
      });
    });

    it('creates an edited-structure save message without extension-side dependencies', () => {
      const api = loadViewerTestApi();
      const message = api.createEditedStructureSaveMessage(
        {
          name: 'Edited',
          atoms: [
            { element: 'cl', x: 1, y: 2, z: 3 },
            { element: 'na', x: 2, y: 2, z: 3 },
          ],
          bonds: [{ from: 1, to: 0, order: 3 }],
        },
        true
      );

      expect(message).toMatchObject({
        type: 'saveEditedStructureToSource',
        dirty: true,
        atomCount: 2,
        bondCount: 1,
      });
      expect(JSON.parse(message.structure)).toMatchObject({
        name: 'Edited',
        atoms: [
          { element: 'Cl', x: 1, y: 2, z: 3 },
          { element: 'Na', x: 2, y: 2, z: 3 },
        ],
        bonds: [{ from: 0, to: 1, order: 3 }],
      });
    });

    it('includes selected atom constraints in edited-structure save messages', () => {
      const api = loadViewerTestApi();
      let state = api.createEditorState(
        {
          name: 'surface',
          atoms: [{ element: 'Pt', x: 0, y: 0, z: 0 }],
        },
        0
      );

      state = api.applyEditingAction(state, {
        type: 'setSelectedAtomConstraints',
        selectiveDynamics: [false, false, true],
      });
      const message = api.createEditedStructureSaveMessage(state.structure, state.dirty);

      expect(message).toMatchObject({
        type: 'saveEditedStructureToSource',
        dirty: true,
        atomCount: 1,
      });
      expect(JSON.parse(message.structure).atoms[0].selectiveDynamics).toEqual([
        false,
        false,
        true,
      ]);
    });

    it('keeps the complete editing surface enabled through exactly 10,000 atoms', () => {
      const api = loadViewerTestApi();
      const structure = createDeterministicStructure(10_000);
      const plan = api.createViewerLoadPlan(structure);

      expect(api.FULL_EDIT_ATOM_LIMIT).toBe(10_000);
      expect(plan.mode).toBe('editable');
      expect(plan.renderedAtomCount).toBe(10_000);
      expect(plan.renderStructure).toBe(structure);
      expect(plan.status).toContain('full editing enabled');
    });

    it('creates a deterministic bounded read-only LOD for a 100,000-atom fixture', () => {
      const api = loadViewerTestApi();
      const structure = createDeterministicStructure(100_000);
      const startedAt = Date.now();
      const first = api.createViewerLoadPlan(structure);
      const second = api.createViewerLoadPlan(structure);
      const elapsedMs = Date.now() - startedAt;

      expect(first.mode).toBe('readonly-lod');
      expect(first.atomCount).toBe(100_000);
      expect(first.renderedAtomCount).toBe(10_000);
      expect(first.renderStructure.atoms).toHaveLength(10_000);
      expect(first.renderStructure.bonds).toEqual([]);
      expect(first.renderStructure.frames).toBeUndefined();
      expect(first.renderStructure.atoms[0]).toBe(structure.atoms[0]);
      expect(first.renderStructure.atoms.at(-1)).toBe(structure.atoms.at(-1));
      expect(first.renderStructure.atoms.map((atom: any) => atom.x)).toEqual(
        second.renderStructure.atoms.map((atom: any) => atom.x)
      );
      expect(first.status).toContain('Read-only LOD');
      expect(elapsedMs).toBeLessThan(2_000);
    }, 5_000);
  });
});

function createDeterministicStructure(atomCount: number): any {
  return {
    schemaVersion: 'openqc.structure.v1',
    kind: 'molecule',
    name: `${atomCount}-atom-fixture`,
    atoms: Array.from({ length: atomCount }, (_, index) => ({
      element: index % 2 === 0 ? 'C' : 'H',
      x: index,
      y: index % 97,
      z: index % 193,
    })),
    bonds: [{ from: 0, to: atomCount - 1, order: 1 }],
    frames: [{ atoms: [] }],
  };
}
