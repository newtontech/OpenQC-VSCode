/**
 * OpenQC Bundled Viewer - Remote-safe Webview viewer shell.
 *
 * This viewer is fully self-contained and does not depend on any CDN.
 * It loads 3Dmol.js from bundled extension assets via webview.asWebviewUri.
 * It receives structure data through postMessage and renders molecular/crystal structures.
 *
 * Message protocol:
 *   loadStructure { structure: JSON string of OpenQCStructure }
 *   setStyle      { style: 'ball-stick' | 'spacefill' | 'line' | 'cartoon' }
 *   showCell      { enabled: boolean }
 *   resetCamera   {}
 *   exportSnapshot {}
 *   exportEditedStructure {}
 *   saveEditedStructureToSource {}
 */

(function () {
  'use strict';

  // -----------------------------------------------------------------------
  // Globals
  // -----------------------------------------------------------------------
  const vscode = acquireVsCodeApi();
  let viewer = null; // 3Dmol viewer
  let currentModel = null;
  let currentStructure = null;
  let cellShape = null;
  let autoRotate = false;
  let labelsVisible = false;
  let measurementMode = false;
  let selectedAtoms = [];
  let playbackInterval = null;
  let editorState = null;
  let renderedStructure = null;
  let viewerMode = 'editable';

  const EDIT_HISTORY_LIMIT = 50;
  const FULL_EDIT_ATOM_LIMIT = 10000;
  const LOD_RENDER_ATOM_LIMIT = 10000;

  // -----------------------------------------------------------------------
  // Initialization
  // -----------------------------------------------------------------------
  function init() {
    const canvas = document.getElementById('viewer-canvas');
    if (!canvas) return;

    // Check WebGL availability
    const testCanvas = document.createElement('canvas');
    const gl = testCanvas.getContext('webgl') || testCanvas.getContext('experimental-webgl');
    if (!gl) {
      showWebGLFallback();
      return;
    }

    try {
      // Use 3Dmol from bundled assets (window.$3Dmol loaded via script tag)
      if (typeof $3Dmol === 'undefined') {
        showError('3Dmol.js library not loaded. Extension assets may be missing.');
        return;
      }

      viewer = $3Dmol.createViewer(canvas, {
        backgroundColor: '0x1e1e1e',
        antialias: true,
      });

      setupControls();
      hideLoading();
      showStatus('Viewer ready. Waiting for structure...');
      vscode.postMessage({ type: 'viewerReady' });
    } catch (e) {
      showError('Failed to initialize viewer: ' + e.message);
    }
  }

  // -----------------------------------------------------------------------
  // Structure loading
  // -----------------------------------------------------------------------
  async function loadStructure(structureJson, options = {}) {
    if (!viewer) return;

    try {
      if (!options.keepPlayback) {
        stopPlayback();
      }

      const structure = normalizeEditableStructure(JSON.parse(structureJson));
      const loadPlan = createViewerLoadPlan(structure);
      viewerMode = loadPlan.mode;
      renderedStructure = loadPlan.renderStructure;
      const selectedIndex = options.keepSelection
        ? normalizeSelectionIndex(options.selectedIndex ?? getPrimarySelectedIndex(), structure)
        : -1;

      currentStructure = structure;
      if (viewerMode === 'readonly-lod') {
        editorState = null;
      } else if (!options.preserveEditing) {
        editorState = createEditorState(structure, selectedIndex);
      } else if (editorState) {
        editorState = finalizeEditingState({
          ...editorState,
          structure: cloneValue(structure),
          selectedIndex,
        });
      }
      selectedAtoms = viewerMode === 'editable' ? selectedAtomArray(structure, selectedIndex) : [];

      // Convert OpenQCStructure to XYZ for 3Dmol
      if (viewerMode === 'readonly-lod') {
        await yieldToMainThread();
      }
      const xyz = structureToXYZ(renderedStructure);

      viewer.removeAllModels();
      viewer.removeAllSurfaces();
      viewer.removeAllLabels();
      viewer.removeAllShapes();

      // Add model from XYZ
      currentModel = viewer.addModel(xyz, 'xyz');
      setupAtomPicking(viewerMode === 'editable');

      // Apply default style
      applyStyle(viewerMode === 'editable' ? 'ball-stick' : 'lod');

      // Show cell if periodic
      if (structure.cell) {
        showPeriodicControls(true);
        renderCell(structure.cell);
      } else {
        showPeriodicControls(false);
        renderSceneShapes(false);
      }

      // Show trajectory controls if frames present
      if (viewerMode === 'editable' && structure.frames && structure.frames.length > 0) {
        showTrajectoryControls(structure);
      } else if (!options.keepPlayback) {
        hideTrajectoryControls();
      }

      renderLabels();
      viewer.zoomTo();
      viewer.render();

      updateAtomCount(structure.atoms.length, loadPlan.renderedAtomCount);
      updateEditingControls();
      showStatus(loadPlan.status);
    } catch (e) {
      showError('Failed to load structure: ' + e.message);
    }
  }

  function structureToXYZ(structure) {
    const atoms = cartesianAtomsForRendering(structure);
    const comment = structure.name || 'OpenQC';
    let xyz = atoms.length + '\n' + comment + '\n';
    for (const atom of atoms) {
      xyz += `${atom.element.padEnd(2)} ${atom.x.toFixed(6)} ${atom.y.toFixed(6)} ${atom.z.toFixed(6)}\n`;
    }
    return xyz;
  }

  function createViewerLoadPlan(structure) {
    const atomCount = structure?.atoms?.length || 0;
    if (atomCount <= FULL_EDIT_ATOM_LIMIT) {
      return {
        mode: 'editable',
        atomCount,
        renderedAtomCount: atomCount,
        renderStructure: structure,
        status: `Loaded ${atomCount} atoms${structure.kind === 'periodic' ? ' (periodic)' : ''}; full editing enabled`,
      };
    }

    const sampledAtoms = deterministicAtomSample(structure.atoms, LOD_RENDER_ATOM_LIMIT);
    return {
      mode: 'readonly-lod',
      atomCount,
      renderedAtomCount: sampledAtoms.length,
      renderStructure: {
        ...structure,
        atoms: sampledAtoms,
        bonds: [],
        frames: undefined,
      },
      status: `Read-only LOD: rendered ${sampledAtoms.length} of ${atomCount} atoms; editing, labels, measurements, trajectory, and supercell are disabled`,
    };
  }

  function deterministicAtomSample(atoms, limit) {
    if (!Array.isArray(atoms) || atoms.length <= limit) return Array.isArray(atoms) ? atoms : [];
    const sampleSize = Math.max(2, Math.floor(limit));
    return Array.from({ length: sampleSize }, (_, index) => {
      const sourceIndex = Math.floor((index * (atoms.length - 1)) / (sampleSize - 1));
      return atoms[sourceIndex];
    });
  }

  function yieldToMainThread() {
    return new Promise(resolve => setTimeout(resolve, 0));
  }

  function cartesianAtomsForRendering(structure) {
    const atoms = structure.atoms || [];
    if (!structure.cell || structure.cell.coordinateMode !== 'fractional') {
      return atoms.map(atom => ({ ...atom }));
    }

    return atoms.map(atom => {
      const [x, y, z] = fractionalToCartesian([atom.x, atom.y, atom.z], structure.cell);
      return { ...atom, x, y, z };
    });
  }

  function fractionalToCartesian(fractional, cell) {
    const [u, v, w] = fractional;
    return [
      u * cell.a[0] + v * cell.b[0] + w * cell.c[0],
      u * cell.a[1] + v * cell.b[1] + w * cell.c[1],
      u * cell.a[2] + v * cell.b[2] + w * cell.c[2],
    ];
  }

  // -----------------------------------------------------------------------
  // Rendering
  // -----------------------------------------------------------------------
  function applyStyle(style) {
    if (!currentModel) return;

    const explicitBonds = hasExplicitBonds(renderedStructure);
    viewer.setStyle({}, {}); // Clear

    switch (style) {
      case 'ball-stick':
        viewer.setStyle(
          {},
          explicitBonds
            ? { sphere: { scale: 0.3, colorscheme: 'Jmol' } }
            : {
                stick: { radius: 0.15, colorscheme: 'Jmol' },
                sphere: { scale: 0.3, colorscheme: 'Jmol' },
              }
        );
        break;
      case 'spacefill':
        viewer.setStyle(
          {},
          {
            sphere: { scale: 1.0, colorscheme: 'Jmol' },
          }
        );
        break;
      case 'line':
        viewer.setStyle(
          {},
          explicitBonds
            ? { sphere: { scale: 0.14, colorscheme: 'Jmol' } }
            : { line: { colorscheme: 'Jmol' } }
        );
        break;
      case 'cartoon':
        viewer.setStyle({}, { cartoon: { colorscheme: 'Jmol' } });
        break;
      case 'lod':
        viewer.setStyle({}, { sphere: { scale: 0.16, colorscheme: 'Jmol' } });
        break;
      default:
        viewer.setStyle(
          {},
          {
            stick: { radius: 0.15, colorscheme: 'Jmol' },
            sphere: { scale: 0.3, colorscheme: 'Jmol' },
          }
        );
    }

    viewer.render();
  }

  function renderSceneShapes(showCell = Boolean(cellShape)) {
    if (!viewer) return;

    viewer.removeAllShapes();
    cellShape = false;

    if (showCell && currentStructure?.cell) {
      drawCellShape(currentStructure.cell);
      cellShape = true;
    }

    renderExplicitBonds(renderedStructure);
  }

  function renderCell(cell) {
    if (!viewer) return;

    renderSceneShapes(true);
  }

  function drawCellShape(cell) {
    if (!viewer || !cell) return;

    const a = cell.a;
    const b = cell.b;
    const c = cell.c;

    const corners = [
      [0, 0, 0],
      a,
      b,
      c,
      addVec(a, b),
      addVec(a, c),
      addVec(b, c),
      addVec(addVec(a, b), c),
    ];

    const edges = [
      [0, 1],
      [0, 2],
      [0, 3],
      [1, 4],
      [1, 5],
      [2, 4],
      [2, 6],
      [3, 5],
      [3, 6],
      [4, 7],
      [5, 7],
      [6, 7],
    ];

    for (const [i, j] of edges) {
      const start = corners[i];
      const end = corners[j];
      viewer.addCylinder({
        start: { x: start[0], y: start[1], z: start[2] },
        end: { x: end[0], y: end[1], z: end[2] },
        radius: 0.05,
        color: 'yellow',
        opacity: 0.6,
      });
    }
  }

  function renderExplicitBonds(structure) {
    if (!viewer || !hasExplicitBonds(structure)) return;

    const atoms = cartesianAtomsForRendering(structure);
    for (const bond of structure.bonds) {
      const normalizedBond = normalizeStructureBond(bond, structure);
      if (!normalizedBond) continue;
      const start = atoms[normalizedBond.from];
      const end = atoms[normalizedBond.to];
      if (!start || !end) continue;
      viewer.addCylinder({
        start: { x: start.x, y: start.y, z: start.z },
        end: { x: end.x, y: end.y, z: end.z },
        radius: 0.06 + (normalizedBond.order - 1) * 0.025,
        color: 'white',
        opacity: 0.9,
      });
    }
  }

  function hasExplicitBonds(structure) {
    return Array.isArray(structure?.bonds) && structure.bonds.length > 0;
  }

  function addVec(a, b) {
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  }

  // -----------------------------------------------------------------------
  // Selection, labels, and measurements
  // -----------------------------------------------------------------------
  function setupAtomPicking(enabled = true) {
    if (!viewer || !enabled) return;
    try {
      viewer.setClickable({}, true, atom => {
        handleAtomClick(atom);
      });
    } catch (e) {
      showStatus('Atom selection is unavailable in this 3Dmol runtime');
    }
  }

  function handleAtomClick(atom) {
    const selected = normalizePickedAtom(atom);
    if (!selected) return;

    if (!measurementMode) {
      const selectedIndex = normalizeSelectionIndex(selected.index, currentStructure);
      if (editorState && selectedIndex >= 0) {
        editorState = applyEditingAction(editorState, { type: 'selectAtom', index: selectedIndex });
      }
      selectedAtoms =
        selectedIndex >= 0 ? selectedAtomArray(currentStructure, selectedIndex) : [selected];
      renderLabels();
      updateEditingControls();
      showStatus(`Selected ${formatAtomLabel(selected)}`);
      return;
    }

    selectedAtoms.push(selected);
    if (selectedAtoms.length > 4) {
      selectedAtoms = [selected];
    }

    renderLabels();

    if (selectedAtoms.length === 1) {
      showStatus(`Measurement start: ${formatAtomLabel(selected)}`);
    } else if (selectedAtoms.length === 2) {
      showStatus(`Distance: ${calculateDistance(selectedAtoms[0], selectedAtoms[1]).toFixed(3)} Å`);
    } else if (selectedAtoms.length === 3) {
      showStatus(
        `Angle: ${calculateAngle(selectedAtoms[0], selectedAtoms[1], selectedAtoms[2]).toFixed(2)}°`
      );
    } else {
      showStatus(
        `Dihedral: ${calculateDihedral(
          selectedAtoms[0],
          selectedAtoms[1],
          selectedAtoms[2],
          selectedAtoms[3]
        ).toFixed(2)}°`
      );
    }
  }

  function normalizePickedAtom(atom) {
    if (!atom) return null;
    const index = Number.isInteger(atom.index)
      ? atom.index
      : Number.isInteger(atom.serial)
        ? atom.serial - 1
        : selectedAtoms.length;
    return {
      element: atom.elem || atom.element || atom.atom || 'X',
      x: Number(atom.x),
      y: Number(atom.y),
      z: Number(atom.z),
      index,
    };
  }

  function renderLabels() {
    if (!viewer) return;
    viewer.removeAllLabels();
    if (viewerMode === 'readonly-lod') return;

    if (labelsVisible && currentStructure) {
      const atoms = cartesianAtomsForRendering(currentStructure);
      atoms.forEach((atom, index) => {
        viewer.addLabel(`${atom.element}${index + 1}`, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          fontColor: 'white',
          backgroundColor: 'black',
          backgroundOpacity: 0.45,
          fontSize: 11,
          inFront: true,
        });
      });
    }

    selectedAtoms.forEach(atom => {
      viewer.addLabel(formatAtomLabel(atom), {
        position: { x: atom.x, y: atom.y, z: atom.z },
        fontColor: 'black',
        backgroundColor: '#f2cc60',
        backgroundOpacity: 0.9,
        fontSize: 12,
        inFront: true,
      });
    });

    if (selectedAtoms.length >= 2) {
      const [a, b] = selectedAtoms;
      const midpoint = midpointBetween(a, b);
      viewer.addLabel(`${calculateDistance(a, b).toFixed(3)} Å`, {
        position: midpoint,
        fontColor: 'black',
        backgroundColor: '#89d185',
        backgroundOpacity: 0.9,
        fontSize: 12,
        inFront: true,
      });
    }

    if (selectedAtoms.length === 3) {
      const [, center] = selectedAtoms;
      viewer.addLabel(`${calculateAngle(selectedAtoms[0], center, selectedAtoms[2]).toFixed(2)}°`, {
        position: { x: center.x, y: center.y, z: center.z },
        fontColor: 'black',
        backgroundColor: '#75beff',
        backgroundOpacity: 0.9,
        fontSize: 12,
        inFront: true,
      });
    }

    if (selectedAtoms.length === 4) {
      const [, b, c] = selectedAtoms;
      viewer.addLabel(
        `${calculateDihedral(
          selectedAtoms[0],
          selectedAtoms[1],
          selectedAtoms[2],
          selectedAtoms[3]
        ).toFixed(2)}°`,
        {
          position: {
            x: (b.x + c.x) / 2,
            y: (b.y + c.y) / 2,
            z: (b.z + c.z) / 2,
          },
          fontColor: '#ffffff',
          backgroundColor: '#1f6feb',
          backgroundOpacity: 0.9,
          fontSize: 12,
          inFront: true,
        }
      );
    }
  }

  function toggleLabels() {
    if (viewerMode === 'readonly-lod') {
      showStatus('Labels are disabled in read-only LOD mode');
      return;
    }
    labelsVisible = !labelsVisible;
    document.getElementById('btn-labels')?.classList.toggle('active', labelsVisible);
    renderLabels();
    viewer?.render();
    showStatus(labelsVisible ? 'Atom labels shown' : 'Atom labels hidden');
  }

  function toggleMeasurementMode() {
    if (viewerMode === 'readonly-lod') {
      showStatus('Measurements are disabled in read-only LOD mode');
      return;
    }
    measurementMode = !measurementMode;
    selectedAtoms = [];
    document.getElementById('btn-measure')?.classList.toggle('active', measurementMode);
    renderLabels();
    viewer?.render();
    showStatus(measurementMode ? 'Measurement mode active' : 'Measurement mode off');
  }

  function formatAtomLabel(atom) {
    return `${atom.element}${Number.isInteger(atom.index) ? atom.index + 1 : ''}`;
  }

  function midpointBetween(a, b) {
    return {
      x: (a.x + b.x) / 2,
      y: (a.y + b.y) / 2,
      z: (a.z + b.z) / 2,
    };
  }

  function calculateDistance(a, b) {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    const dz = a.z - b.z;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  function calculateAngle(a, b, c) {
    const ab = { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
    const cb = { x: c.x - b.x, y: c.y - b.y, z: c.z - b.z };
    const dot = ab.x * cb.x + ab.y * cb.y + ab.z * cb.z;
    const denom =
      Math.sqrt(ab.x * ab.x + ab.y * ab.y + ab.z * ab.z) *
      Math.sqrt(cb.x * cb.x + cb.y * cb.y + cb.z * cb.z);
    if (denom === 0) return 0;
    const cosine = Math.min(1, Math.max(-1, dot / denom));
    return (Math.acos(cosine) * 180) / Math.PI;
  }

  function calculateDihedral(a, b, c, d) {
    const b0 = vectorBetween(b, a);
    const b1 = normalizePointVector(vectorBetween(b, c));
    const b2 = vectorBetween(c, d);
    const v = subtractPointVector(b0, scalePointVector(b1, dotPointVector(b0, b1)));
    const w = subtractPointVector(b2, scalePointVector(b1, dotPointVector(b2, b1)));
    const x = dotPointVector(v, w);
    const y = dotPointVector(crossPointVector(v, b1), w);
    if (!Number.isFinite(x) || !Number.isFinite(y)) return 0;
    return (Math.atan2(y, x) * 180) / Math.PI;
  }

  function vectorBetween(from, to) {
    return {
      x: finiteNumber(to.x) - finiteNumber(from.x),
      y: finiteNumber(to.y) - finiteNumber(from.y),
      z: finiteNumber(to.z) - finiteNumber(from.z),
    };
  }

  function normalizePointVector(vector) {
    const length = Math.sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    if (!length) return { x: 0, y: 0, z: 0 };
    return {
      x: vector.x / length,
      y: vector.y / length,
      z: vector.z / length,
    };
  }

  function scalePointVector(vector, scale) {
    return {
      x: vector.x * scale,
      y: vector.y * scale,
      z: vector.z * scale,
    };
  }

  function subtractPointVector(a, b) {
    return {
      x: a.x - b.x,
      y: a.y - b.y,
      z: a.z - b.z,
    };
  }

  function dotPointVector(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  function crossPointVector(a, b) {
    return {
      x: a.y * b.z - a.z * b.y,
      y: a.z * b.x - a.x * b.z,
      z: a.x * b.y - a.y * b.x,
    };
  }

  // -----------------------------------------------------------------------
  // Editing state and actions
  // -----------------------------------------------------------------------
  function createEditorState(structure, selectedIndex = -1) {
    const normalizedStructure = normalizeEditableStructure(structure || { atoms: [] });
    const baselineJson = stableStructureJson(normalizedStructure);
    return finalizeEditingState({
      structure: normalizedStructure,
      selectedIndex,
      undoStack: [],
      redoStack: [],
      baselineJson,
      dirty: false,
    });
  }

  function applyEditingAction(state, action) {
    const current = normalizeEditorState(state);
    const editAction = action || {};

    switch (editAction.type) {
      case 'selectAtom':
        return finalizeEditingState({
          ...current,
          selectedIndex: normalizeSelectionIndex(editAction.index, current.structure),
        });

      case 'addAtom': {
        const nextStructure = structureForEditMutation(current.structure);
        const atom = normalizeStructureAtom(editAction.atom || editAction);
        nextStructure.atoms.push(atom);
        return commitEditingState(current, nextStructure, nextStructure.atoms.length - 1);
      }

      case 'deleteSelectedAtom': {
        const selectedIndex = selectedIndexForAction(editAction, current);
        if (selectedIndex < 0) return current;

        const nextStructure = structureForEditMutation(current.structure);
        nextStructure.atoms.splice(selectedIndex, 1);
        nextStructure.bonds = remapBondsAfterAtomDelete(nextStructure.bonds, selectedIndex);
        const nextSelectedIndex = Math.min(selectedIndex, nextStructure.atoms.length - 1);
        return commitEditingState(current, nextStructure, nextSelectedIndex);
      }

      case 'moveSelectedAtom': {
        const selectedIndex = selectedIndexForAction(editAction, current);
        if (selectedIndex < 0) return current;

        const nextStructure = structureForEditMutation(current.structure);
        const atom = nextStructure.atoms[selectedIndex];
        const delta = vectorFromAction(editAction);
        nextStructure.atoms[selectedIndex] = {
          ...atom,
          x: atom.x + delta.x,
          y: atom.y + delta.y,
          z: atom.z + delta.z,
        };
        return commitEditingState(current, nextStructure, selectedIndex);
      }

      case 'setSelectedAtomConstraints': {
        const selectedIndex = selectedIndexForAction(editAction, current);
        const selectiveDynamics = normalizeSelectiveDynamics(
          editAction.selectiveDynamics || editAction.constraints || editAction
        );
        if (selectedIndex < 0 || !selectiveDynamics) return current;

        const nextStructure = structureForEditMutation(current.structure);
        const atom = nextStructure.atoms[selectedIndex];
        nextStructure.atoms[selectedIndex] = {
          ...atom,
          selectiveDynamics,
        };
        return commitEditingState(current, nextStructure, selectedIndex);
      }

      case 'addOrUpdateBond': {
        const nextStructure = structureForEditMutation(current.structure);
        const bond = normalizeStructureBond(editAction, nextStructure);
        if (!bond) return current;

        const bondIndex = findBondIndex(nextStructure.bonds, bond.from, bond.to);
        if (bondIndex >= 0) {
          nextStructure.bonds[bondIndex] = bond;
        } else {
          nextStructure.bonds.push(bond);
        }
        return commitEditingState(current, nextStructure, current.selectedIndex);
      }

      case 'deleteBond': {
        const nextStructure = structureForEditMutation(current.structure);
        const pair = normalizeBondPair(editAction, nextStructure);
        if (!pair) return current;

        const bondIndex = findBondIndex(nextStructure.bonds, pair.from, pair.to);
        if (bondIndex < 0) return current;

        nextStructure.bonds.splice(bondIndex, 1);
        return commitEditingState(current, nextStructure, current.selectedIndex);
      }

      case 'undo':
        return popEditingHistory(current, 'undo');

      case 'redo':
        return popEditingHistory(current, 'redo');

      default:
        return current;
    }
  }

  function normalizeEditorState(state) {
    if (!state || !state.structure) {
      return createEditorState({ atoms: [] });
    }

    const structure = normalizeEditableStructure(state.structure);
    const baselineJson =
      typeof state.baselineJson === 'string' ? state.baselineJson : stableStructureJson(structure);

    return finalizeEditingState({
      structure,
      selectedIndex: state.selectedIndex,
      undoStack: normalizeHistoryStack(state.undoStack),
      redoStack: normalizeHistoryStack(state.redoStack),
      baselineJson,
      dirty: false,
    });
  }

  function commitEditingState(state, structure, selectedIndex) {
    return finalizeEditingState({
      ...state,
      structure: normalizeEditableStructure(structure),
      selectedIndex,
      undoStack: pushHistorySnapshot(state.undoStack, snapshotEditingState(state)),
      redoStack: [],
    });
  }

  function popEditingHistory(state, direction) {
    const fromStack = direction === 'undo' ? state.undoStack : state.redoStack;
    if (!fromStack.length) return state;

    const toStack = direction === 'undo' ? state.redoStack : state.undoStack;
    const snapshot = fromStack[fromStack.length - 1];
    const nextState = {
      ...state,
      structure: normalizeEditableStructure(snapshot.structure),
      selectedIndex: snapshot.selectedIndex,
      undoStack:
        direction === 'undo'
          ? fromStack.slice(0, -1)
          : pushHistorySnapshot(toStack, snapshotEditingState(state)),
      redoStack:
        direction === 'undo'
          ? pushHistorySnapshot(toStack, snapshotEditingState(state))
          : fromStack.slice(0, -1),
    };

    return finalizeEditingState(nextState);
  }

  function finalizeEditingState(state) {
    const structure = normalizeEditableStructure(state.structure);
    const baselineJson =
      typeof state.baselineJson === 'string' ? state.baselineJson : stableStructureJson(structure);

    return {
      structure,
      selectedIndex: normalizeSelectionIndex(state.selectedIndex, structure),
      undoStack: normalizeHistoryStack(state.undoStack),
      redoStack: normalizeHistoryStack(state.redoStack),
      baselineJson,
      dirty: stableStructureJson(structure) !== baselineJson,
    };
  }

  function snapshotEditingState(state) {
    return {
      structure: normalizeEditableStructure(state.structure),
      selectedIndex: normalizeSelectionIndex(state.selectedIndex, state.structure),
    };
  }

  function normalizeHistoryStack(stack) {
    if (!Array.isArray(stack)) return [];
    return stack.slice(-EDIT_HISTORY_LIMIT).map(snapshot => ({
      structure: normalizeEditableStructure(snapshot.structure),
      selectedIndex: snapshot.selectedIndex,
    }));
  }

  function pushHistorySnapshot(stack, snapshot) {
    return normalizeHistoryStack([...(Array.isArray(stack) ? stack : []), snapshot]);
  }

  function normalizeEditableStructure(structure) {
    const normalized = cloneValue(structure || {});
    normalized.atoms = Array.isArray(normalized.atoms)
      ? normalized.atoms.map(atom => normalizeStructureAtom(atom))
      : [];
    normalized.bonds = Array.isArray(normalized.bonds)
      ? normalized.bonds
          .map(bond => normalizeStructureBond(bond, normalized))
          .filter(bond => Boolean(bond))
      : [];
    return normalized;
  }

  function normalizeStructureAtom(atom) {
    const source = atom || {};
    const normalized = {
      ...source,
      element: normalizeElement(source.element || source.elem || source.atom || 'X'),
      x: finiteNumber(source.x),
      y: finiteNumber(source.y),
      z: finiteNumber(source.z),
    };
    const selectiveDynamics = normalizeSelectiveDynamics(source.selectiveDynamics);
    if (selectiveDynamics) {
      normalized.selectiveDynamics = selectiveDynamics;
    } else {
      delete normalized.selectiveDynamics;
    }
    return normalized;
  }

  function normalizeElement(value) {
    const letters = String(value || 'X')
      .replace(/[^a-z]/gi, '')
      .slice(0, 2);
    if (!letters) return 'X';
    return letters.charAt(0).toUpperCase() + letters.slice(1).toLowerCase();
  }

  function normalizeStructureBond(bond, structure) {
    const pair = normalizeBondPair(bond, structure);
    if (!pair) return null;
    const order = normalizeBondOrder(bond?.order);
    if (!order) return null;
    return {
      from: pair.from,
      to: pair.to,
      order,
    };
  }

  function normalizeBondPair(action, structure) {
    const from = normalizeSelectionIndex(action?.from, structure);
    const to = normalizeSelectionIndex(action?.to, structure);
    if (from < 0 || to < 0 || from === to) return null;
    return from < to ? { from, to } : { from: to, to: from };
  }

  function normalizeBondOrder(value) {
    if (value === undefined || value === null || value === '') {
      return 1;
    }
    if (typeof value !== 'number' || !Number.isInteger(value)) {
      return null;
    }
    const order = value;
    if (order < 1 || order > 3) return null;
    return order;
  }

  function normalizeSelectiveDynamics(value) {
    if (!Array.isArray(value) || value.length !== 3) {
      return null;
    }

    const flags = value.map(normalizeSelectiveFlag);
    return flags.every(flag => typeof flag === 'boolean') ? flags : null;
  }

  function normalizeSelectiveFlag(value) {
    if (typeof value === 'boolean') {
      return value;
    }
    if (typeof value === 'string') {
      const normalized = value.trim().toLowerCase();
      if (normalized === 't' || normalized === 'true') {
        return true;
      }
      if (normalized === 'f' || normalized === 'false') {
        return false;
      }
    }
    return null;
  }

  function findBondIndex(bonds, from, to) {
    const pair = from < to ? { from, to } : { from: to, to: from };
    return (Array.isArray(bonds) ? bonds : []).findIndex(
      bond => bond.from === pair.from && bond.to === pair.to
    );
  }

  function remapBondsAfterAtomDelete(bonds, deletedIndex) {
    if (!Array.isArray(bonds)) return [];
    return bonds
      .filter(bond => bond.from !== deletedIndex && bond.to !== deletedIndex)
      .map(bond => ({
        ...bond,
        from: bond.from > deletedIndex ? bond.from - 1 : bond.from,
        to: bond.to > deletedIndex ? bond.to - 1 : bond.to,
      }));
  }

  function structureForEditMutation(structure) {
    const nextStructure = normalizeEditableStructure(structure);
    delete nextStructure.frames;
    return nextStructure;
  }

  function selectedIndexForAction(action, state) {
    const explicitIndex = action && action.index;
    const index = explicitIndex === 0 || explicitIndex ? explicitIndex : state.selectedIndex;
    return normalizeSelectionIndex(index, state.structure);
  }

  function normalizeSelectionIndex(index, structure) {
    const numericIndex = Number(index);
    const atoms = structure?.atoms || [];
    if (!Number.isInteger(numericIndex) || numericIndex < 0 || numericIndex >= atoms.length) {
      return -1;
    }
    return numericIndex;
  }

  function selectedAtomArray(structure, selectedIndex) {
    const atom = atomSelectionAtIndex(structure, selectedIndex);
    return atom ? [atom] : [];
  }

  function atomSelectionAtIndex(structure, selectedIndex) {
    const index = normalizeSelectionIndex(selectedIndex, structure);
    if (index < 0) return null;

    const atoms = cartesianAtomsForRendering(structure);
    const atom = atoms[index];
    if (!atom) return null;

    return {
      element: atom.element,
      x: finiteNumber(atom.x),
      y: finiteNumber(atom.y),
      z: finiteNumber(atom.z),
      index,
    };
  }

  function getPrimarySelectedIndex() {
    return selectedAtoms.length ? selectedAtoms[0].index : -1;
  }

  function vectorFromAction(action) {
    const delta = action.delta || action;
    return {
      x: finiteNumber(delta.x ?? delta.dx),
      y: finiteNumber(delta.y ?? delta.dy),
      z: finiteNumber(delta.z ?? delta.dz),
    };
  }

  function cartesianToStructureCoordinates(point, structure) {
    if (!usesFractionalCoordinates(structure)) {
      return {
        x: finiteNumber(point.x),
        y: finiteNumber(point.y),
        z: finiteNumber(point.z),
      };
    }

    const fractional = cartesianVectorToFractional(
      [finiteNumber(point.x), finiteNumber(point.y), finiteNumber(point.z)],
      structure.cell
    );
    return { x: fractional[0], y: fractional[1], z: fractional[2] };
  }

  function cartesianDeltaToStructureDelta(delta, structure) {
    return cartesianToStructureCoordinates(delta, structure);
  }

  function usesFractionalCoordinates(structure) {
    return Boolean(structure?.cell && structure.cell.coordinateMode === 'fractional');
  }

  function cartesianVectorToFractional(vector, cell) {
    const a = cell.a;
    const b = cell.b;
    const c = cell.c;
    const determinant = dotVec(a, crossVec(b, c));
    if (!determinant) return vector;

    return [
      dotVec(vector, crossVec(b, c)) / determinant,
      dotVec(vector, crossVec(c, a)) / determinant,
      dotVec(vector, crossVec(a, b)) / determinant,
    ];
  }

  function crossVec(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }

  function dotVec(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  function finiteNumber(value, fallback = 0) {
    const numericValue = Number(value);
    return Number.isFinite(numericValue) ? numericValue : fallback;
  }

  function stableStructureJson(structure) {
    return JSON.stringify(normalizeEditableStructure(structure));
  }

  function cloneValue(value) {
    return JSON.parse(JSON.stringify(value || {}));
  }

  function applyViewerEdit(action, statusMessage) {
    if (viewerMode !== 'editable') {
      showStatus('Editing is disabled above 10,000 atoms (read-only LOD mode)');
      return;
    }
    if (!currentStructure || !editorState) {
      showStatus('Load a structure before editing');
      return;
    }

    const previousJson = stableStructureJson(editorState.structure);
    const previousSelectedIndex = editorState.selectedIndex;
    const previousUndoCount = editorState.undoStack.length;
    const previousRedoCount = editorState.redoStack.length;
    const nextState = applyEditingAction(editorState, action);

    editorState = nextState;
    currentStructure = normalizeEditableStructure(nextState.structure);
    selectedAtoms = selectedAtomArray(currentStructure, nextState.selectedIndex);

    const changed =
      previousJson !== stableStructureJson(nextState.structure) ||
      previousSelectedIndex !== nextState.selectedIndex ||
      previousUndoCount !== nextState.undoStack.length ||
      previousRedoCount !== nextState.redoStack.length;

    if (changed) {
      loadStructure(JSON.stringify(currentStructure), {
        preserveEditing: true,
        keepSelection: true,
        selectedIndex: nextState.selectedIndex,
      });
      postStructureUpdate(nextState);
    } else {
      renderLabels();
      viewer?.render();
      updateEditingControls();
    }

    if (statusMessage) {
      showStatus(statusMessage);
    }
  }

  function addAtomFromInputs() {
    if (!currentStructure) {
      showStatus('Load a structure before adding atoms');
      return;
    }

    const element = normalizeElement(document.getElementById('edit-element')?.value || 'C');
    const point = cartesianToStructureCoordinates(
      {
        x: readNumberInput('edit-x'),
        y: readNumberInput('edit-y'),
        z: readNumberInput('edit-z'),
      },
      currentStructure
    );
    const label = `${element}${currentStructure.atoms.length + 1}`;
    applyViewerEdit({ type: 'addAtom', atom: { element, ...point } }, `Added ${label}`);
  }

  function deleteSelectedAtom() {
    const selected = atomSelectionAtIndex(currentStructure, editorState?.selectedIndex);
    if (!selected) {
      showStatus('Select an atom before deleting');
      return;
    }

    applyViewerEdit({ type: 'deleteSelectedAtom' }, `Deleted ${formatAtomLabel(selected)}`);
  }

  function moveSelectedAtomFromInputs() {
    const selected = atomSelectionAtIndex(currentStructure, editorState?.selectedIndex);
    if (!selected) {
      showStatus('Select an atom before moving');
      return;
    }

    const cartesianDelta = {
      x: readNumberInput('move-dx'),
      y: readNumberInput('move-dy'),
      z: readNumberInput('move-dz'),
    };
    const structureDelta = cartesianDeltaToStructureDelta(cartesianDelta, currentStructure);

    applyViewerEdit(
      {
        type: 'moveSelectedAtom',
        dx: structureDelta.x,
        dy: structureDelta.y,
        dz: structureDelta.z,
      },
      `Moved ${formatAtomLabel(selected)} by ${cartesianDelta.x}, ${cartesianDelta.y}, ${cartesianDelta.z}`
    );
  }

  function addOrUpdateBondFromInputs() {
    if (!currentStructure) {
      showStatus('Load a structure before editing bonds');
      return;
    }

    const pair = readBondInputPair();
    if (!pair) {
      showStatus('Enter two different valid atom indices before editing bonds');
      return;
    }

    const order = readBondOrderInput();
    applyViewerEdit(
      { type: 'addOrUpdateBond', from: pair.from, to: pair.to, order },
      `Bond ${pair.from + 1}-${pair.to + 1} set to order ${order}`
    );
  }

  function deleteBondFromInputs() {
    if (!currentStructure) {
      showStatus('Load a structure before editing bonds');
      return;
    }

    const pair = readBondInputPair();
    if (!pair) {
      showStatus('Enter two different valid atom indices before deleting a bond');
      return;
    }

    applyViewerEdit(
      { type: 'deleteBond', from: pair.from, to: pair.to },
      `Bond ${pair.from + 1}-${pair.to + 1} deleted`
    );
  }

  function setSelectedAtomConstraintsFromInputs() {
    const selected = atomSelectionAtIndex(currentStructure, editorState?.selectedIndex);
    if (!selected) {
      showStatus('Select an atom before setting constraints');
      return;
    }

    const selectiveDynamics = readConstraintInputs();
    applyViewerEdit(
      { type: 'setSelectedAtomConstraints', selectiveDynamics },
      `Constraints updated for ${formatAtomLabel(selected)}`
    );
  }

  function undoEdit() {
    if (!editorState?.undoStack.length) {
      showStatus('Nothing to undo');
      return;
    }
    applyViewerEdit({ type: 'undo' }, 'Undo');
  }

  function redoEdit() {
    if (!editorState?.redoStack.length) {
      showStatus('Nothing to redo');
      return;
    }
    applyViewerEdit({ type: 'redo' }, 'Redo');
  }

  function createEditedStructureExportMessage(structure, dirty) {
    const normalizedStructure = normalizeEditableStructure(structure);
    return {
      type: 'exportEditedStructure',
      structure: JSON.stringify(normalizedStructure),
      dirty: Boolean(dirty),
      atomCount: normalizedStructure.atoms.length,
      bondCount: normalizedStructure.bonds.length,
    };
  }

  function createEditedStructureSaveMessage(structure, dirty) {
    const normalizedStructure = normalizeEditableStructure(structure);
    return {
      type: 'saveEditedStructureToSource',
      structure: JSON.stringify(normalizedStructure),
      dirty: Boolean(dirty),
      atomCount: normalizedStructure.atoms.length,
      bondCount: normalizedStructure.bonds.length,
    };
  }

  function createStructureUpdateMessage(state) {
    const normalizedStructure = normalizeEditableStructure(state?.structure || currentStructure);
    return {
      type: 'structureUpdated',
      structure: JSON.stringify(normalizedStructure),
      dirty: Boolean(state?.dirty),
      atomCount: normalizedStructure.atoms.length,
      bondCount: normalizedStructure.bonds.length,
    };
  }

  function postStructureUpdate(state) {
    vscode.postMessage(createStructureUpdateMessage(state));
  }

  function exportEditedStructure() {
    if (!currentStructure) {
      showStatus('Load a structure before exporting edits');
      return;
    }

    vscode.postMessage(createEditedStructureExportMessage(currentStructure, editorState?.dirty));
    showStatus(editorState?.dirty ? 'Edited structure exported' : 'Structure exported');
  }

  function saveEditedStructureToSource() {
    if (!currentStructure) {
      showStatus('Load a structure before saving edits');
      return;
    }

    if (!editorState?.dirty) {
      showStatus('No edits to save');
      return;
    }

    vscode.postMessage(createEditedStructureSaveMessage(currentStructure, editorState?.dirty));
    showStatus('Saving edited structure...');
  }

  function markStructureSaved(structureJson) {
    const structure = structureJson
      ? normalizeEditableStructure(JSON.parse(structureJson))
      : currentStructure;
    currentStructure = structure;
    renderedStructure = structure;
    viewerMode = 'editable';
    editorState = createEditorState(structure, editorState?.selectedIndex);
    selectedAtoms = selectedAtomArray(currentStructure, editorState.selectedIndex);
    updateEditingControls();
    renderLabels();
    viewer?.render();
    showStatus('Edited structure saved to source');
  }

  function readNumberInput(id, fallback = 0) {
    return finiteNumber(document.getElementById(id)?.value, fallback);
  }

  function readBondInputPair() {
    const from = Math.round(readNumberInput('bond-from', 1)) - 1;
    const to = Math.round(readNumberInput('bond-to', 2)) - 1;
    const pair = normalizeBondPair({ from, to }, currentStructure);
    return pair;
  }

  function readBondOrderInput() {
    const order = normalizeBondOrder(readNumberInput('bond-order', 1));
    return order || 1;
  }

  function readConstraintInputs() {
    return ['constraint-x', 'constraint-y', 'constraint-z'].map(id =>
      Boolean(document.getElementById(id)?.checked)
    );
  }

  function updateEditingControls() {
    const hasStructure = Boolean(currentStructure?.atoms);
    const isEditable = viewerMode === 'editable';
    const hasSelection =
      isEditable && normalizeSelectionIndex(editorState?.selectedIndex, currentStructure) >= 0;
    const hasBondCapableStructure = isEditable && Boolean(currentStructure?.atoms?.length >= 2);

    setDisabled('btn-add-atom', !hasStructure || !isEditable);
    setDisabled('btn-delete-atom', !hasStructure || !hasSelection);
    setDisabled('btn-move-atom', !hasStructure || !hasSelection);
    setDisabled('btn-add-bond', !hasBondCapableStructure);
    setDisabled('btn-delete-bond', !hasBondCapableStructure);
    setDisabled('constraint-x', !hasStructure || !hasSelection);
    setDisabled('constraint-y', !hasStructure || !hasSelection);
    setDisabled('constraint-z', !hasStructure || !hasSelection);
    setDisabled('btn-set-constraints', !hasStructure || !hasSelection);
    setDisabled('btn-undo-edit', !isEditable || !editorState?.undoStack.length);
    setDisabled('btn-redo-edit', !isEditable || !editorState?.redoStack.length);
    setDisabled('btn-save-source', !hasStructure || !editorState?.dirty);
    setDisabled('btn-export-structure', !hasStructure);
    setDisabled('btn-labels', !isEditable);
    setDisabled('btn-measure', !isEditable);
    setDisabled('btn-supercell', !isEditable);
    syncConstraintInputs();

    const bondCount = document.getElementById('bond-count');
    if (bondCount) {
      const count = currentStructure?.bonds?.length || 0;
      bondCount.textContent = `${count} ${count === 1 ? 'bond' : 'bonds'}`;
    }

    const indicator = document.getElementById('dirty-indicator');
    if (indicator) {
      indicator.textContent = isEditable ? (editorState?.dirty ? 'Edited' : 'Clean') : 'Read-only LOD';
      indicator.classList.toggle('dirty', Boolean(editorState?.dirty));
      indicator.title = !isEditable
        ? 'Structures above 10,000 atoms use a deterministic read-only level of detail'
        : editorState?.dirty
          ? 'Structure has webview edits that are not written back to the source file'
          : 'No webview edits relative to the loaded structure';
    }
  }

  function setDisabled(id, disabled) {
    const element = document.getElementById(id);
    if (element) element.disabled = disabled;
  }

  function syncConstraintInputs() {
    const selected = normalizeSelectionIndex(editorState?.selectedIndex, currentStructure);
    const atom = selected >= 0 ? currentStructure?.atoms?.[selected] : null;
    const flags = normalizeSelectiveDynamics(atom?.selectiveDynamics) || [true, true, true];
    ['constraint-x', 'constraint-y', 'constraint-z'].forEach((id, index) => {
      const input = document.getElementById(id);
      if (input) {
        input.checked = flags[index];
      }
    });
  }

  // -----------------------------------------------------------------------
  // Export
  // -----------------------------------------------------------------------
  function exportSnapshot() {
    if (!viewer) return;

    try {
      const imgData = viewer.pngURI();
      if (imgData) {
        // Send base64 data URL instead of Blob
        vscode.postMessage({
          type: 'exportImage',
          data: imgData,
        });
        showStatus('Image exported');
      }
    } catch (e) {
      showError('Export failed: ' + e.message);
    }
  }

  // -----------------------------------------------------------------------
  // UI Controls
  // -----------------------------------------------------------------------
  function setupControls() {
    // Style selector
    document.getElementById('style-select')?.addEventListener('change', e => {
      applyStyle(e.target.value);
    });

    // Reset camera
    document.getElementById('btn-reset')?.addEventListener('click', () => {
      if (viewer) {
        viewer.zoomTo();
        viewer.render();
      }
    });

    // Auto rotate
    document.getElementById('btn-rotate')?.addEventListener('click', () => {
      autoRotate = !autoRotate;
      // 3Dmol doesn't have built-in auto-rotate, simulate with spin
      if (viewer) {
        viewer.spin(autoRotate);
      }
      const btn = document.getElementById('btn-rotate');
      if (btn) {
        btn.classList.toggle('active', autoRotate);
      }
    });

    // Atom labels
    document.getElementById('btn-labels')?.addEventListener('click', () => {
      toggleLabels();
    });

    // Distance/angle measurement
    document.getElementById('btn-measure')?.addEventListener('click', () => {
      toggleMeasurementMode();
    });

    // Export image
    document.getElementById('btn-export')?.addEventListener('click', () => {
      exportSnapshot();
    });

    // Editing controls
    document.getElementById('btn-add-atom')?.addEventListener('click', () => {
      addAtomFromInputs();
    });

    document.getElementById('btn-delete-atom')?.addEventListener('click', () => {
      deleteSelectedAtom();
    });

    document.getElementById('btn-move-atom')?.addEventListener('click', () => {
      moveSelectedAtomFromInputs();
    });

    document.getElementById('btn-add-bond')?.addEventListener('click', () => {
      addOrUpdateBondFromInputs();
    });

    document.getElementById('btn-delete-bond')?.addEventListener('click', () => {
      deleteBondFromInputs();
    });

    document.getElementById('btn-set-constraints')?.addEventListener('click', () => {
      setSelectedAtomConstraintsFromInputs();
    });

    document.getElementById('btn-undo-edit')?.addEventListener('click', () => {
      undoEdit();
    });

    document.getElementById('btn-redo-edit')?.addEventListener('click', () => {
      redoEdit();
    });

    document.getElementById('btn-save-source')?.addEventListener('click', () => {
      saveEditedStructureToSource();
    });

    document.getElementById('btn-export-structure')?.addEventListener('click', () => {
      exportEditedStructure();
    });

    // Toggle unit cell
    document.getElementById('btn-toggle-cell')?.addEventListener('click', () => {
      if (currentStructure?.cell) {
        if (cellShape) {
          renderSceneShapes(false);
          cellShape = null;
          showStatus('Unit cell hidden');
        } else {
          renderSceneShapes(true);
          showStatus('Unit cell shown');
        }
        viewer.render();
      }
    });

    // Supercell controls
    document.getElementById('btn-supercell')?.addEventListener('click', () => {
      generateSupercell();
    });

    updateEditingControls();
  }

  function generateSupercell() {
    if (!currentStructure?.cell || viewerMode !== 'editable') return;

    const na = parseInt(document.getElementById('sc-na')?.value || '2');
    const nb = parseInt(document.getElementById('sc-nb')?.value || '2');
    const nc = parseInt(document.getElementById('sc-nc')?.value || '2');

    const maxAtoms = 10000;
    const currentCount = currentStructure.atoms.length;
    const totalAtoms = currentCount * na * nb * nc;

    if (totalAtoms > maxAtoms) {
      showStatus(`⚠️ Supercell would have ${totalAtoms} atoms (max ${maxAtoms}). Reduce size.`);
      return;
    }

    const superStructure = buildSupercellStructure(currentStructure, na, nb, nc);

    loadStructure(JSON.stringify(superStructure));
    showStatus(`Supercell ${na}x${nb}x${nc}: ${totalAtoms} atoms`);
  }

  function buildSupercellStructure(structure, na, nb, nc) {
    const normalizedStructure = normalizeEditableStructure(structure);
    const cell = normalizedStructure.cell;
    if (!cell) {
      return normalizedStructure;
    }

    const repeatA = positiveInteger(na, 1);
    const repeatB = positiveInteger(nb, 1);
    const repeatC = positiveInteger(nc, 1);
    const baseAtoms = cartesianAtomsForRendering(normalizedStructure);
    const superAtoms = [];
    const imageCount = repeatA * repeatB * repeatC;
    const superBonds = replicateBondsForSupercell(
      normalizedStructure.bonds,
      normalizedStructure.atoms.length,
      imageCount
    );

    for (let ia = 0; ia < repeatA; ia++) {
      for (let ib = 0; ib < repeatB; ib++) {
        for (let ic = 0; ic < repeatC; ic++) {
          for (const atom of baseAtoms) {
            const dx = ia * cell.a[0] + ib * cell.b[0] + ic * cell.c[0];
            const dy = ia * cell.a[1] + ib * cell.b[1] + ic * cell.c[1];
            const dz = ia * cell.a[2] + ib * cell.b[2] + ic * cell.c[2];
            superAtoms.push({
              ...atom,
              x: atom.x + dx,
              y: atom.y + dy,
              z: atom.z + dz,
            });
          }
        }
      }
    }

    return {
      ...normalizedStructure,
      atoms: superAtoms,
      bonds: superBonds,
      cell: {
        ...cell,
        a: [cell.a[0] * repeatA, cell.a[1] * repeatA, cell.a[2] * repeatA],
        b: [cell.b[0] * repeatB, cell.b[1] * repeatB, cell.b[2] * repeatB],
        c: [cell.c[0] * repeatC, cell.c[1] * repeatC, cell.c[2] * repeatC],
        pbc: cell.pbc,
        coordinateMode: 'cartesian',
      },
    };
  }

  function positiveInteger(value, fallback) {
    const numericValue = Number(value);
    if (!Number.isInteger(numericValue) || numericValue < 1) {
      return fallback;
    }
    return numericValue;
  }

  function replicateBondsForSupercell(bonds, atomCount, imageCount) {
    if (!Array.isArray(bonds) || atomCount <= 0 || imageCount <= 0) {
      return [];
    }

    const normalizedBonds = bonds
      .map(bond => normalizeStructureBond(bond, { atoms: new Array(atomCount) }))
      .filter(bond => Boolean(bond));
    const replicated = [];
    for (let imageIndex = 0; imageIndex < imageCount; imageIndex++) {
      const offset = imageIndex * atomCount;
      for (const bond of normalizedBonds) {
        replicated.push({
          from: bond.from + offset,
          to: bond.to + offset,
          order: bond.order,
        });
      }
    }
    return replicated;
  }

  // -----------------------------------------------------------------------
  // Trajectory controls
  // -----------------------------------------------------------------------
  function showTrajectoryControls(structure) {
    const controls = document.getElementById('trajectory-controls');
    if (!controls || !structure.frames) return;

    controls.classList.add('visible');
    const slider = document.getElementById('frame-slider');
    if (slider) {
      slider.max = Math.max(0, structure.frames.length - 1);
      slider.value = 0;
      slider.oninput = () => {
        const idx = parseInt(slider.value);
        loadFrame(idx, structure);
      };
    }

    // Play/pause button
    const playButton = document.getElementById('btn-play');
    if (playButton) {
      playButton.onclick = () => {
        togglePlayback(structure);
      };
      playButton.textContent = '▶ Play';
    }

    const info = document.getElementById('frame-info');
    if (info) {
      info.textContent = `Frame 1/${structure.frames.length}`;
    }
  }

  function hideTrajectoryControls() {
    const controls = document.getElementById('trajectory-controls');
    if (controls) {
      controls.classList.remove('visible');
    }
  }

  function loadFrame(index, structure) {
    if (!structure.frames || index >= structure.frames.length) return;
    const frame = structure.frames[index];
    const combined = { ...structure, ...frame, frames: undefined };
    loadStructure(JSON.stringify(combined), { keepPlayback: true });

    const info = document.getElementById('frame-info');
    if (info) {
      info.textContent = `Frame ${index + 1}/${structure.frames.length}`;
    }
  }

  function togglePlayback(structure) {
    if (playbackInterval) {
      stopPlayback();
    } else {
      const slider = document.getElementById('frame-slider');
      const playButton = document.getElementById('btn-play');
      let idx = parseInt(slider?.value || '0');
      if (playButton) {
        playButton.textContent = 'Pause';
      }
      playbackInterval = setInterval(() => {
        if (!structure.frames) return;
        idx = (idx + 1) % structure.frames.length;
        if (slider) slider.value = idx;
        loadFrame(idx, structure);
      }, 200);
      showStatus('Playing trajectory...');
    }
  }

  function stopPlayback() {
    if (!playbackInterval) return;
    clearInterval(playbackInterval);
    playbackInterval = null;
    const playButton = document.getElementById('btn-play');
    if (playButton) {
      playButton.textContent = '▶ Play';
    }
    showStatus('Playback paused');
  }

  // -----------------------------------------------------------------------
  // Periodic controls visibility
  // -----------------------------------------------------------------------
  function showPeriodicControls(visible) {
    const controls = document.getElementById('periodic-controls');
    if (controls) {
      controls.classList.toggle('visible', visible);
    }
  }

  // -----------------------------------------------------------------------
  // UI helpers
  // -----------------------------------------------------------------------
  function showStatus(msg) {
    const el = document.getElementById('status-text');
    if (el) el.textContent = msg;
  }

  function updateAtomCount(count, renderedCount = count) {
    const el = document.getElementById('atom-count');
    if (el) {
      el.textContent = renderedCount === count ? `${count} atoms` : `${count} atoms (${renderedCount} LOD)`;
    }
  }

  function showError(msg) {
    const el = document.getElementById('error-panel');
    if (el) {
      el.textContent = msg;
      el.style.display = 'block';
      setTimeout(() => {
        el.style.display = 'none';
      }, 5000);
    }
  }

  function hideLoading() {
    const el = document.getElementById('loading-overlay');
    if (el) el.classList.add('hidden');
  }

  function showWebGLFallback() {
    hideLoading();
    const fallback = document.getElementById('webgl-fallback');
    if (fallback) fallback.classList.add('visible');
  }

  // -----------------------------------------------------------------------
  // Message handler
  // -----------------------------------------------------------------------
  window.addEventListener('message', event => {
    const message = event.data;

    switch (message.type) {
      case 'loadStructure':
        if (message.structure) {
          loadStructure(message.structure);
        }
        break;

      case 'initialize':
        // Backward compat: XYZ format
        if (message.structure?.xyz) {
          if (!viewer) init();
          viewer.removeAllModels();
          currentModel = viewer.addModel(message.structure.xyz, 'xyz');
          applyStyle('ball-stick');
          viewer.zoomTo();
          viewer.render();
          hideLoading();
        }
        break;

      case 'setStyle':
        if (message.style) {
          applyStyle(message.style);
        }
        break;

      case 'showCell':
        if (message.enabled && currentStructure?.cell) {
          renderSceneShapes(true);
        } else {
          renderSceneShapes(false);
          cellShape = null;
        }
        viewer?.render();
        break;

      case 'resetCamera':
        viewer?.zoomTo();
        viewer?.render();
        break;

      case 'exportSnapshot':
        exportSnapshot();
        break;

      case 'markStructureSaved':
        if (message.structure) {
          markStructureSaved(message.structure);
        }
        break;
    }
  });

  if (typeof window !== 'undefined') {
    window.__openqcViewerTest = {
      cartesianAtomsForRendering,
      fractionalToCartesian,
      structureToXYZ,
      calculateDistance,
      calculateAngle,
      calculateDihedral,
      createEditorState,
      applyEditingAction,
      createEditedStructureExportMessage,
      createEditedStructureSaveMessage,
      normalizeElement,
      normalizeStructureBond,
      normalizeBondPair,
      normalizeBondOrder,
      normalizeSelectiveDynamics,
      buildSupercellStructure,
      replicateBondsForSupercell,
      cartesianToStructureCoordinates,
      cartesianDeltaToStructureDelta,
      markStructureSaved,
      handleAtomClick,
      createViewerLoadPlan,
      deterministicAtomSample,
      FULL_EDIT_ATOM_LIMIT,
      LOD_RENDER_ATOM_LIMIT,
    };
  }

  // Initialize when DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
