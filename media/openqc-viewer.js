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
  function loadStructure(structureJson) {
    if (!viewer) return;

    try {
      const structure = JSON.parse(structureJson);
      currentStructure = structure;

      // Convert OpenQCStructure to XYZ for 3Dmol
      const xyz = structureToXYZ(structure);

      viewer.removeAllModels();
      viewer.removeAllSurfaces();
      viewer.removeAllLabels();
      viewer.removeAllShapes();

      // Add model from XYZ
      currentModel = viewer.addModel(xyz, 'xyz');

      // Apply default style
      applyStyle('ball-stick');

      // Show cell if periodic
      if (structure.cell) {
        showPeriodicControls(true);
        renderCell(structure.cell);
      } else {
        showPeriodicControls(false);
        if (cellShape) {
          viewer.removeAllShapes();
          cellShape = null;
        }
      }

      // Show trajectory controls if frames present
      if (structure.frames && structure.frames.length > 0) {
        showTrajectoryControls(structure);
      }

      viewer.zoomTo();
      viewer.render();

      updateAtomCount(structure.atoms.length);
      showStatus(`Loaded ${structure.atoms.length} atoms${structure.kind === 'periodic' ? ' (periodic)' : ''}`);
    } catch (e) {
      showError('Failed to load structure: ' + e.message);
    }
  }

  function structureToXYZ(structure) {
    const atoms = structure.atoms || [];
    const comment = structure.name || 'OpenQC';
    let xyz = atoms.length + '\n' + comment + '\n';
    for (const atom of atoms) {
      xyz += `${atom.element.padEnd(2)} ${atom.x.toFixed(6)} ${atom.y.toFixed(6)} ${atom.z.toFixed(6)}\n`;
    }
    return xyz;
  }

  // -----------------------------------------------------------------------
  // Rendering
  // -----------------------------------------------------------------------
  function applyStyle(style) {
    if (!currentModel) return;

    viewer.setStyle({}, {}); // Clear

    switch (style) {
      case 'ball-stick':
        viewer.setStyle({}, {
          stick: { radius: 0.15, colorscheme: 'Jmol' },
          sphere: { scale: 0.3, colorscheme: 'Jmol' },
        });
        break;
      case 'spacefill':
        viewer.setStyle({}, {
          sphere: { scale: 1.0, colorscheme: 'Jmol' },
        });
        break;
      case 'line':
        viewer.setStyle({}, { line: { colorscheme: 'Jmol' } });
        break;
      case 'cartoon':
        viewer.setStyle({}, { cartoon: { colorscheme: 'Jmol' } });
        break;
      default:
        viewer.setStyle({}, {
          stick: { radius: 0.15, colorscheme: 'Jmol' },
          sphere: { scale: 0.3, colorscheme: 'Jmol' },
        });
    }

    viewer.render();
  }

  function renderCell(cell) {
    if (!viewer) return;

    viewer.removeAllShapes();

    const a = cell.a;
    const b = cell.b;
    const c = cell.c;

    const corners = [
      [0, 0, 0],
      a, b, c,
      addVec(a, b),
      addVec(a, c),
      addVec(b, c),
      addVec(addVec(a, b), c),
    ];

    const edges = [
      [0, 1], [0, 2], [0, 3],
      [1, 4], [1, 5],
      [2, 4], [2, 6],
      [3, 5], [3, 6],
      [4, 7], [5, 7], [6, 7],
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

    cellShape = true;
  }

  function addVec(a, b) {
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
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
    document.getElementById('style-select')?.addEventListener('change', (e) => {
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

    // Export image
    document.getElementById('btn-export')?.addEventListener('click', () => {
      exportSnapshot();
    });

    // Toggle unit cell
    document.getElementById('btn-toggle-cell')?.addEventListener('click', () => {
      if (currentStructure?.cell) {
        if (cellShape) {
          viewer.removeAllShapes();
          cellShape = null;
          showStatus('Unit cell hidden');
        } else {
          renderCell(currentStructure.cell);
          showStatus('Unit cell shown');
        }
        viewer.render();
      }
    });

    // Supercell controls
    document.getElementById('btn-supercell')?.addEventListener('click', () => {
      generateSupercell();
    });
  }

  function generateSupercell() {
    if (!currentStructure?.cell) return;

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

    // Generate supercell by replicating atoms
    const cell = currentStructure.cell;
    const baseAtoms = currentStructure.atoms;
    const superAtoms = [];

    for (let ia = 0; ia < na; ia++) {
      for (let ib = 0; ib < nb; ib++) {
        for (let ic = 0; ic < nc; ic++) {
          for (const atom of baseAtoms) {
            const dx = ia * cell.a[0] + ib * cell.b[0] + ic * cell.c[0];
            const dy = ia * cell.a[1] + ib * cell.b[1] + ic * cell.c[1];
            const dz = ia * cell.a[2] + ib * cell.b[2] + ic * cell.c[2];
            superAtoms.push({
              element: atom.element,
              x: atom.x + dx,
              y: atom.y + dy,
              z: atom.z + dz,
            });
          }
        }
      }
    }

    // Create supercell structure and display it
    const superStructure = {
      ...currentStructure,
      atoms: superAtoms,
      cell: {
        a: [cell.a[0] * na, cell.a[1] * na, cell.a[2] * na],
        b: [cell.b[0] * nb, cell.b[1] * nb, cell.b[2] * nb],
        c: [cell.c[0] * nc, cell.c[1] * nc, cell.c[2] * nc],
        pbc: cell.pbc,
      },
    };

    loadStructure(JSON.stringify(superStructure));
    showStatus(`Supercell ${na}x${nb}x${nc}: ${totalAtoms} atoms`);
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
      slider.max = structure.frames.length;
      slider.value = 0;
      slider.oninput = () => {
        const idx = parseInt(slider.value);
        loadFrame(idx, structure);
      };
    }

    // Play/pause button
    document.getElementById('btn-play')?.addEventListener('click', () => {
      togglePlayback(structure);
    });
  }

  function loadFrame(index, structure) {
    if (!structure.frames || index >= structure.frames.length) return;
    const frame = structure.frames[index];
    const combined = { ...structure, ...frame, frames: undefined };
    loadStructure(JSON.stringify(combined));

    const info = document.getElementById('frame-info');
    if (info) {
      info.textContent = `Frame ${index + 1}/${structure.frames.length}`;
    }
  }

  let playbackInterval = null;
  function togglePlayback(structure) {
    if (playbackInterval) {
      clearInterval(playbackInterval);
      playbackInterval = null;
      showStatus('Playback paused');
    } else {
      const slider = document.getElementById('frame-slider');
      let idx = 0;
      playbackInterval = setInterval(() => {
        if (!structure.frames) return;
        idx = (idx + 1) % structure.frames.length;
        if (slider) slider.value = idx;
        loadFrame(idx, structure);
      }, 200);
      showStatus('Playing trajectory...');
    }
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

  function updateAtomCount(count) {
    const el = document.getElementById('atom-count');
    if (el) el.textContent = `${count} atoms`;
  }

  function showError(msg) {
    const el = document.getElementById('error-panel');
    if (el) {
      el.textContent = msg;
      el.style.display = 'block';
      setTimeout(() => { el.style.display = 'none'; }, 5000);
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
  window.addEventListener('message', (event) => {
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
          renderCell(currentStructure.cell);
        } else {
          viewer?.removeAllShapes();
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
    }
  });

  // Initialize when DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
