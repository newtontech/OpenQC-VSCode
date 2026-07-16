#!/usr/bin/env node

import fs from 'fs';
import { createRequire } from 'module';
import path from 'path';
import { pathToFileURL } from 'url';
import { chromium } from 'playwright';
import sharp from 'sharp';

const require = createRequire(import.meta.url);
const root = path.resolve(new URL('..', import.meta.url).pathname);
const { buildVisualSmokePaths, parseVisualSmokeArgs } = require(
  path.join(root, 'out', 'smoke', 'visualSmokePaths.js')
);
const options = parseVisualSmokeArgs(process.argv.slice(2), root);
const { extensionRoot, outputDir, outputPrefix } = options;
const { htmlPath, screenshotPath, narrowScreenshotPath, reportPath } =
  buildVisualSmokePaths(options);

fs.mkdirSync(outputDir, { recursive: true });

function requireProductionWebview(targetExtensionRoot) {
  const webviewModulePath = path.join(
    targetExtensionRoot,
    'out',
    'webviews',
    'openqcViewerWebview.js'
  );
  if (!fs.existsSync(webviewModulePath)) {
    throw new Error(
      `Compiled webview generator not found in ${targetExtensionRoot}. Run \`npm run compile\` before \`npm run test:visual\`, or pass an extracted VSIX extension root.`
    );
  }

  const moduleApi = require('module');
  const originalLoad = moduleApi._load;
  const vscodeMock = {
    Uri: {
      joinPath(base, ...parts) {
        const basePath = base.fsPath || base.path;
        const fsPath = path.join(basePath, ...parts);
        return { fsPath, path: fsPath };
      },
    },
  };

  try {
    moduleApi._load = function loadWithVscodeMock(request, parent, isMain) {
      if (request === 'vscode') {
        return vscodeMock;
      }
      return originalLoad.call(this, request, parent, isMain);
    };
    return require(webviewModulePath).OpenQCViewerWebview;
  } finally {
    moduleApi._load = originalLoad;
  }
}

const OpenQCViewerWebview = requireProductionWebview(extensionRoot);
const productionHtml = OpenQCViewerWebview.generateHTML(
  {
    asWebviewUri(uri) {
      return pathToFileURL(uri.fsPath || uri.path).href;
    },
  },
  { fsPath: extensionRoot, path: extensionRoot }
);

const nonce = productionHtml.match(/script-src 'nonce-([^']+)'/)?.[1];
if (!nonce) {
  throw new Error('Unable to extract CSP nonce from production OpenQC viewer HTML.');
}

const smokeStyle = `<style nonce="${nonce}">
  #visual-smoke-ready {
    position: fixed;
    right: 8px;
    top: 8px;
    z-index: 20;
    padding: 4px 8px;
    border-radius: 4px;
    background: #1f6feb;
    color: #fff;
    font: 12px sans-serif;
  }
</style>`;

const bootstrapScript = `<script nonce="${nonce}">
    window.__openqcMessages = [];
    window.acquireVsCodeApi = () => ({
      postMessage(message) {
        window.__openqcMessages.push(message);
      },
    });
    window.addEventListener('error', event => {
      const marker = document.getElementById('visual-smoke-ready');
      if (marker) {
        marker.dataset.state = 'error';
        marker.textContent = event.message || 'visual smoke error';
      }
    });
</script>`;

const driverScript = `<script nonce="${nonce}">
    const smokeStructure = {
      schemaVersion: 'openqc.structure.v1',
      name: 'OpenQC visual smoke NaCl',
      kind: 'periodic',
      atoms: [
        { element: 'Na', x: 0, y: 0, z: 0 },
        { element: 'Cl', x: 0.5, y: 0.5, z: 0.5 },
        { element: 'Na', x: 0, y: 0.5, z: 0.5 },
        { element: 'Cl', x: 0.5, y: 0, z: 0 },
      ],
      cell: {
        a: [5.64, 0, 0],
        b: [0, 5.64, 0],
        c: [0, 0, 5.64],
        pbc: [true, true, true],
        coordinateMode: 'fractional',
      },
      frames: [
        { atoms: [{ element: 'Na', x: 0, y: 0, z: 0 }, { element: 'Cl', x: 0.5, y: 0.5, z: 0.5 }] },
        { atoms: [{ element: 'Na', x: 0.01, y: 0, z: 0 }, { element: 'Cl', x: 0.5, y: 0.49, z: 0.5 }] },
      ],
    };

    window.addEventListener('load', () => {
      setTimeout(() => {
        window.postMessage({ type: 'loadStructure', structure: JSON.stringify(smokeStructure) }, '*');
      }, 100);

      setTimeout(() => {
        document.getElementById('btn-labels')?.click();
        const element = document.getElementById('edit-element');
        const x = document.getElementById('edit-x');
        const y = document.getElementById('edit-y');
        const z = document.getElementById('edit-z');
        if (element) element.value = 'O';
        if (x) x.value = '2.8';
        if (y) y.value = '1.4';
        if (z) z.value = '1.4';
        document.getElementById('btn-add-atom')?.click();

        window.__openqcViewerTest?.handleAtomClick({ elem: 'O', index: 4, x: 2.8, y: 1.4, z: 1.4 });
        const dx = document.getElementById('move-dx');
        const dy = document.getElementById('move-dy');
        const dz = document.getElementById('move-dz');
        if (dx) dx.value = '0.2';
        if (dy) dy.value = '0.1';
        if (dz) dz.value = '0';
        document.getElementById('btn-undo-edit')?.click();
        document.getElementById('btn-redo-edit')?.click();

        const moveEnabledBeforeClick = !document.getElementById('btn-move-atom')?.disabled;
        document.getElementById('btn-move-atom')?.click();
        document.getElementById('btn-delete-atom')?.click();
        document.getElementById('btn-undo-edit')?.click();
        document.getElementById('btn-undo-edit')?.click();
        document.getElementById('btn-redo-edit')?.click();

        const constraintX = document.getElementById('constraint-x');
        const constraintY = document.getElementById('constraint-y');
        const constraintZ = document.getElementById('constraint-z');
        if (constraintX) constraintX.checked = false;
        if (constraintY) constraintY.checked = false;
        if (constraintZ) constraintZ.checked = true;
        const constraintsEnabledBeforeClick = !document.getElementById('btn-set-constraints')?.disabled;
        document.getElementById('btn-set-constraints')?.click();

        const beforeBondMessageCount = window.__openqcMessages.length;
        const bondFrom = document.getElementById('bond-from');
        const bondTo = document.getElementById('bond-to');
        const bondOrder = document.getElementById('bond-order');
        if (bondFrom) bondFrom.value = '1';
        if (bondTo) bondTo.value = '5';
        if (bondOrder) bondOrder.value = '2';
        document.getElementById('btn-add-bond')?.click();
        if (bondOrder) bondOrder.value = '3';
        document.getElementById('btn-add-bond')?.click();
        document.getElementById('btn-delete-bond')?.click();
        document.getElementById('btn-add-bond')?.click();

        document.getElementById('btn-measure')?.click();
        window.__openqcViewerTest?.handleAtomClick({ elem: 'Na', index: 0, x: 0, y: 0, z: 0 });
        window.__openqcViewerTest?.handleAtomClick({ elem: 'Cl', index: 1, x: 2.82, y: 2.82, z: 2.82 });
        window.__openqcViewerTest?.handleAtomClick({ elem: 'Na', index: 2, x: 0, y: 2.82, z: 2.82 });
        window.__openqcViewerTest?.handleAtomClick({ elem: 'Cl', index: 3, x: 2.82, y: 0, z: 0 });
        const measurementStatus = document.getElementById('status-text')?.textContent || '';
        document.getElementById('btn-measure')?.click();

        document.getElementById('btn-export-structure')?.click();
        const dirtyIndicatorBeforeSave = document.getElementById('dirty-indicator')?.textContent || '';
        document.getElementById('btn-save-source')?.click();
        const saveMessage = window.__openqcMessages.findLast?.(
          message => message.type === 'saveEditedStructureToSource'
        ) || [...window.__openqcMessages].reverse().find(
          message => message.type === 'saveEditedStructureToSource'
        );
        if (saveMessage?.structure) {
          window.__openqcViewerTest?.markStructureSaved(saveMessage.structure);
        }
        document.getElementById('btn-export')?.click();
        const atomCount = document.getElementById('atom-count')?.textContent || '';
        const bondCount = document.getElementById('bond-count')?.textContent || '';
        const dirtyIndicatorAfterSave = document.getElementById('dirty-indicator')?.textContent || '';
        const scNa = document.getElementById('sc-na');
        const scNb = document.getElementById('sc-nb');
        const scNc = document.getElementById('sc-nc');
        if (scNa) scNa.value = '2';
        if (scNb) scNb.value = '1';
        if (scNc) scNc.value = '1';
        document.getElementById('btn-supercell')?.click();
        const supercellAtomCount = document.getElementById('atom-count')?.textContent || '';
        const supercellBondCount = document.getElementById('bond-count')?.textContent || '';
        document.getElementById('btn-export-structure')?.click();
        const structureUpdates = window.__openqcMessages
          .filter(message => message.type === 'structureUpdated')
          .map(message => ({
            ...message,
            parsed: JSON.parse(message.structure),
          }));
        const bondEditUpdates = window.__openqcMessages
          .slice(beforeBondMessageCount)
          .filter(message => message.type === 'structureUpdated')
          .map(message => ({
            ...message,
            parsed: JSON.parse(message.structure),
          }));
        const exported = window.__openqcMessages.some(message => {
          if (message.type !== 'exportEditedStructure') return false;
          const parsed = JSON.parse(message.structure);
          return (
            parsed.atoms.length === 5 &&
            parsed.bonds?.length === 1 &&
            parsed.bonds[0]?.from === 0 &&
            parsed.bonds[0]?.to === 4 &&
            parsed.bonds[0]?.order === 3 &&
            JSON.stringify(parsed.atoms[4]?.selectiveDynamics) === JSON.stringify([false, false, true]) &&
            message.bondCount === 1 &&
            message.dirty === true
          );
        });
        const saved = Boolean(
          saveMessage &&
          (() => {
            const parsed = JSON.parse(saveMessage.structure);
            return (
              parsed.atoms.length === 5 &&
              parsed.bonds?.length === 1 &&
              parsed.bonds[0]?.from === 0 &&
              parsed.bonds[0]?.to === 4 &&
              parsed.bonds[0]?.order === 3 &&
              JSON.stringify(parsed.atoms[4]?.selectiveDynamics) === JSON.stringify([false, false, true]) &&
              saveMessage.bondCount === 1 &&
              saveMessage.dirty === true
            );
          })()
        );
        const supercellExported = window.__openqcMessages.some(message => {
          if (message.type !== 'exportEditedStructure') return false;
          const parsed = JSON.parse(message.structure);
          return (
            parsed.atoms.length === 10 &&
            parsed.bonds?.length === 2 &&
            parsed.bonds[0]?.from === 0 &&
            parsed.bonds[0]?.to === 4 &&
            parsed.bonds[0]?.order === 3 &&
            parsed.bonds[1]?.from === 5 &&
            parsed.bonds[1]?.to === 9 &&
            parsed.bonds[1]?.order === 3 &&
            JSON.stringify(parsed.atoms[4]?.selectiveDynamics) === JSON.stringify([false, false, true]) &&
            JSON.stringify(parsed.atoms[9]?.selectiveDynamics) === JSON.stringify([false, false, true]) &&
            parsed.cell?.coordinateMode === 'cartesian'
          );
        });
        const exportedImage = window.__openqcMessages.some(message =>
          message.type === 'exportImage' &&
          typeof message.data === 'string' &&
          message.data.startsWith('data:image/png')
        );
        const sawEditedUpdate = structureUpdates.some(message =>
          message.parsed.atoms.length === 5 && message.dirty === true
        );
        const sawUndoUpdate = structureUpdates.some(message =>
          message.parsed.atoms.length === 4 && message.dirty === false
        );
        const sawMoveUpdate = structureUpdates.some(message =>
          message.parsed.atoms.length === 5 &&
          message.dirty === true &&
          Number(message.parsed.atoms[4]?.x) > 0.49
        );
        const sawConstraintUpdate = structureUpdates.some(message =>
          message.parsed.atoms.length === 5 &&
          message.dirty === true &&
          JSON.stringify(message.parsed.atoms[4]?.selectiveDynamics) === JSON.stringify([false, false, true])
        );
        const sawDeleteUpdate = structureUpdates.some(message =>
          message.parsed.atoms.length === 4 && message.dirty === true
        );
        const sawBondOrderTwo = bondEditUpdates.some(message =>
          message.bondCount === 1 &&
          message.parsed.bonds?.length === 1 &&
          message.parsed.bonds[0]?.from === 0 &&
          message.parsed.bonds[0]?.to === 4 &&
          message.parsed.bonds[0]?.order === 2
        );
        const sawBondOrderThree = bondEditUpdates.some(message =>
          message.bondCount === 1 &&
          message.parsed.bonds?.length === 1 &&
          message.parsed.bonds[0]?.from === 0 &&
          message.parsed.bonds[0]?.to === 4 &&
          message.parsed.bonds[0]?.order === 3
        );
        const sawBondDelete = bondEditUpdates.some(message =>
          message.bondCount === 0 && Array.isArray(message.parsed.bonds) && message.parsed.bonds.length === 0
        );
        const ready =
          atomCount.includes('5 atoms') &&
          bondCount.includes('1 bond') &&
          supercellAtomCount.includes('10 atoms') &&
          supercellBondCount.includes('2 bonds') &&
          dirtyIndicatorBeforeSave === 'Edited' &&
          dirtyIndicatorAfterSave === 'Clean' &&
          moveEnabledBeforeClick &&
          constraintsEnabledBeforeClick &&
          measurementStatus.startsWith('Dihedral:') &&
          sawEditedUpdate &&
          sawUndoUpdate &&
          sawMoveUpdate &&
          sawConstraintUpdate &&
          sawDeleteUpdate &&
          sawBondOrderTwo &&
          sawBondOrderThree &&
          sawBondDelete &&
          exported &&
          saved &&
          supercellExported &&
          exportedImage;
        const marker = document.getElementById('visual-smoke-ready');
        const checks = {
          atomCount: atomCount.includes('5 atoms'),
          bondCount: bondCount.includes('1 bond'),
          supercellAtomCount: supercellAtomCount.includes('10 atoms'),
          supercellBondCount: supercellBondCount.includes('2 bonds'),
          dirtyBefore: dirtyIndicatorBeforeSave === 'Edited',
          cleanAfter: dirtyIndicatorAfterSave === 'Clean',
          moveEnabledBeforeClick,
          constraintsEnabledBeforeClick,
          measurement: measurementStatus.startsWith('Dihedral:'),
          sawEditedUpdate,
          sawUndoUpdate,
          sawMoveUpdate,
          sawConstraintUpdate,
          sawDeleteUpdate,
          sawBondOrderTwo,
          sawBondOrderThree,
          sawBondDelete,
          exported,
          saved,
          supercellExported,
          exportedImage,
        };
        const failedChecks = Object.entries(checks)
          .filter(([, passed]) => !passed)
          .map(([name]) => name);
        marker.dataset.state = ready ? 'ready' : 'error';
        marker.textContent = ready ? 'ready' : 'viewer smoke failed: ' + failedChecks.join(', ');
      }, 1800);
    });
  </script>`;

const html = productionHtml
  .replace('<title>OpenQC Structure Viewer</title>', '<title>OpenQC Viewer Visual Smoke</title>')
  .replace('</head>', `${smokeStyle}\n${bootstrapScript}\n</head>`)
  .replace('<body>', '<body>\n  <div id="visual-smoke-ready" data-state="booting">booting</div>')
  .replace('</body>', `${driverScript}\n</body>`);

fs.writeFileSync(htmlPath, html, 'utf8');

const screenshots = [
  {
    name: 'desktop',
    path: screenshotPath,
    viewport: '1280,800',
    minWidth: 1000,
    minHeight: 700,
    minBytes: 20000,
  },
  {
    name: 'narrow',
    path: narrowScreenshotPath,
    viewport: '430,780',
    minWidth: 390,
    minHeight: 700,
    minBytes: 12000,
  },
];

async function captureAndInspect(browser, screenshot) {
  const [width, height] = screenshot.viewport.split(',').map(Number);
  const page = await browser.newPage({ viewport: { width, height } });
  const diagnostics = [];
  page.on('console', message => diagnostics.push(`console.${message.type()}: ${message.text()}`));
  page.on('pageerror', error => diagnostics.push(`pageerror: ${error.message}`));
  page.on('requestfailed', request => {
    diagnostics.push(
      `requestfailed: ${request.url()} (${request.failure()?.errorText || 'unknown error'})`
    );
  });
  page.setDefaultTimeout(20_000);
  try {
    await page.goto(pathToFileURL(htmlPath).href, { waitUntil: 'load', timeout: 20_000 });
    try {
      await page.waitForFunction(() => {
        const state = document.getElementById('visual-smoke-ready')?.dataset.state;
        return state === 'ready' || state === 'error';
      });
    } catch (error) {
      const state = await page.locator('#visual-smoke-ready').evaluate(element => ({
        state: element.dataset.state,
        text: element.textContent,
      }));
      throw new Error(
        `Viewer smoke timed out in state ${state.state || 'missing'} (${state.text || 'no status'}). ${diagnostics.join(' | ') || 'No browser diagnostics were emitted.'}`,
        { cause: error }
      );
    }
    const marker = await page.locator('#visual-smoke-ready').evaluate(element => ({
      state: element.dataset.state,
      text: element.textContent,
    }));
    if (marker.state !== 'ready') {
      throw new Error(`Viewer smoke reported an error: ${marker.text || 'unknown error'}`);
    }
    await page.waitForTimeout(500);
    await page.screenshot({ path: screenshot.path });
  } finally {
    await page.close();
  }

  const image = sharp(screenshot.path).removeAlpha();
  const { data, info } = await image.raw().toBuffer({ resolveWithObject: true });
  const uniqueSamples = new Set();
  for (
    let offset = 0;
    offset < data.length;
    offset += Math.max(3, Math.floor(data.length / 6000))
  ) {
    const aligned = offset - (offset % 3);
    uniqueSamples.add(`${data[aligned]},${data[aligned + 1]},${data[aligned + 2]}`);
  }

  const entry = {
    name: screenshot.name,
    screenshotPath: screenshot.path,
    viewport: screenshot.viewport,
    width: info.width,
    height: info.height,
    channels: info.channels,
    uniqueSampleCount: uniqueSamples.size,
    bytes: fs.statSync(screenshot.path).size,
  };

  if (info.width < screenshot.minWidth || info.height < screenshot.minHeight) {
    throw new Error(
      `${screenshot.name} screenshot dimensions are too small: ${info.width}x${info.height}`
    );
  }

  if (entry.bytes < screenshot.minBytes || uniqueSamples.size < 16) {
    throw new Error(
      `${screenshot.name} screenshot appears blank or under-rendered: ${entry.bytes} bytes, ${uniqueSamples.size} sampled colors`
    );
  }

  return entry;
}

const screenshotReports = [];
const browser = await chromium.launch({ headless: true, timeout: 20_000 });
try {
  for (const screenshot of screenshots) {
    screenshotReports.push(await captureAndInspect(browser, screenshot));
  }
} finally {
  await browser.close();
}

const report = {
  extensionRoot,
  htmlPath,
  outputPrefix,
  screenshotPath,
  screenshots: screenshotReports,
};

fs.writeFileSync(reportPath, `${JSON.stringify(report, null, 2)}\n`, 'utf8');

console.log(`OpenQC viewer visual smoke passed (${outputPrefix}): ${screenshotPath}`);
console.log(`OpenQC viewer narrow smoke passed (${outputPrefix}): ${narrowScreenshotPath}`);
console.log(`Visual smoke report: ${reportPath}`);
