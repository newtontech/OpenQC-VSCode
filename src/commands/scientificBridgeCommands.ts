/**
 * VS Code commands for scientific Python bridges.
 * @module commands/scientificBridgeCommands
 */

import * as vscode from 'vscode';
import {
  parseStructure,
  generateSupercell as generateSupercellBridge,
} from '../python/StructureBridge';
import { extractTrajectory, parseOutput } from '../results/OutputParserBridge';
import { OpenQCViewerPanel } from '../visualizers/OpenQCViewerPanel';
import type { OpenQCResults } from '../results/OpenQCResults';
import type { OpenQCStructure } from '../structures/OpenQCStructure';

export function registerScientificBridgeCommands(context: vscode.ExtensionContext): void {
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.parseStructurePython', async () => {
      await parseStructureCommand();
    }),
    vscode.commands.registerCommand('openqc.generateSupercell', async () => {
      await generateSupercellCommand(context);
    }),
    vscode.commands.registerCommand('openqc.parseCalculationOutput', async () => {
      await parseCalculationOutputCommand();
    }),
    vscode.commands.registerCommand('openqc.showSCFConvergence', async () => {
      await showScfConvergenceCommand(context);
    }),
    vscode.commands.registerCommand('openqc.showOptimizationTrajectory', async () => {
      await showOptimizationTrajectoryCommand(context);
    })
  );
}

async function parseStructureCommand(): Promise<void> {
  const document = getActiveDocument();
  if (!document) {
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Parsing structure with OpenQC Python bridge...',
      cancellable: false,
    },
    async () => {
      const result = await parseStructure(document.uri.fsPath, undefined, document.getText());
      if (!result.success || !result.data) {
        vscode.window.showErrorMessage(`Structure parse failed: ${formatBridgeError(result)}`);
        return;
      }

      await showJsonDocument(result.data, 'json');
      vscode.window.showInformationMessage(
        `Parsed structure: ${result.data.atoms.length} atoms (${result.data.kind})`
      );
    }
  );
}

async function generateSupercellCommand(context: vscode.ExtensionContext): Promise<void> {
  const document = getActiveDocument();
  if (!document) {
    return;
  }

  const dims = await vscode.window.showInputBox({
    prompt: 'Supercell dimensions',
    value: '2 2 2',
    placeHolder: 'na nb nc',
    validateInput: value => {
      const parsed = parseSupercellDims(value);
      return parsed ? undefined : 'Enter three positive integers, for example: 2 2 2';
    },
  });
  if (!dims) {
    return;
  }

  const [na, nb, nc] = parseSupercellDims(dims)!;

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Generating ${na}x${nb}x${nc} supercell...`,
      cancellable: false,
    },
    async () => {
      const parsed = await parseStructure(document.uri.fsPath, undefined, document.getText());
      if (!parsed.success || !parsed.data) {
        vscode.window.showErrorMessage(`Structure parse failed: ${formatBridgeError(parsed)}`);
        return;
      }

      const result = await generateSupercellBridge(parsed.data, na, nb, nc);
      if (!result.success || !result.data) {
        vscode.window.showErrorMessage(`Supercell generation failed: ${formatBridgeError(result)}`);
        return;
      }

      OpenQCViewerPanel.createOrShow(
        context.extensionUri,
        result.data,
        `${document.fileName} ${na}x${nb}x${nc}`
      );
      vscode.window.showInformationMessage(
        `Generated supercell with ${result.data.atoms.length} atoms`
      );
    }
  );
}

async function parseCalculationOutputCommand(): Promise<void> {
  const document = getActiveDocument();
  if (!document) {
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Parsing calculation output...',
      cancellable: false,
    },
    async () => {
      const result = await parseOutput(document.uri.fsPath);
      if (!result.success || !result.data) {
        vscode.window.showErrorMessage(`Output parse failed: ${formatBridgeError(result)}`);
        return;
      }
      if (result.data.success === false) {
        vscode.window.showWarningMessage(
          result.data.warnings?.[0] ??
            result.data.errors?.[0] ??
            'No calculation output data found.'
        );
        return;
      }

      await showJsonDocument(result.data, 'json');
      const energy = result.data.finalEnergy
        ? ` Final energy: ${result.data.finalEnergy.value} ${result.data.finalEnergy.unit}.`
        : '';
      vscode.window.showInformationMessage(`Parsed calculation output.${energy}`);
    }
  );
}

async function showScfConvergenceCommand(context: vscode.ExtensionContext): Promise<void> {
  const document = getActiveDocument();
  if (!document) {
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Parsing SCF convergence...',
      cancellable: false,
    },
    async () => {
      const result = await parseOutput(document.uri.fsPath);
      if (!result.success || !result.data) {
        vscode.window.showErrorMessage(`Output parse failed: ${formatBridgeError(result)}`);
        return;
      }

      const energies = result.data.scfEnergies ?? result.data.optimization?.energies ?? [];
      if (energies.length === 0) {
        vscode.window.showWarningMessage(
          'No SCF or optimization energy series found in this output.'
        );
        return;
      }

      const panel = vscode.window.createWebviewPanel(
        'openqc.scfConvergence',
        'OpenQC: SCF Convergence',
        vscode.ViewColumn.Two,
        { enableScripts: false, localResourceRoots: [context.extensionUri] }
      );
      panel.webview.html = renderScfConvergenceHtml(result.data, energies);
    }
  );
}

async function showOptimizationTrajectoryCommand(context: vscode.ExtensionContext): Promise<void> {
  const document = getActiveDocument();
  if (!document) {
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Extracting optimization trajectory...',
      cancellable: false,
    },
    async () => {
      const result = await extractTrajectory(document.uri.fsPath);
      if (!result.success || !result.data) {
        vscode.window.showErrorMessage(
          `Trajectory extraction failed: ${formatBridgeError(result)}`
        );
        return;
      }

      if (!result.data.supported || result.data.frames.length === 0) {
        vscode.window.showWarningMessage(
          result.data.warnings?.[0] ?? 'No optimization trajectory frames found in this output.'
        );
        return;
      }

      const trajectory = toTrajectoryStructure(result.data.frames, document.fileName);
      OpenQCViewerPanel.createOrShow(
        context.extensionUri,
        trajectory,
        `${document.fileName} trajectory`
      );
      vscode.window.showInformationMessage(
        `Loaded optimization trajectory with ${trajectory.frames?.length ?? 0} frames`
      );
    }
  );
}

function getActiveDocument(): vscode.TextDocument | undefined {
  const document = vscode.window.activeTextEditor?.document;
  if (!document) {
    vscode.window.showErrorMessage('No active text editor');
  }
  return document;
}

async function showJsonDocument(data: unknown, language: string): Promise<void> {
  const doc = await vscode.workspace.openTextDocument({
    content: JSON.stringify(data, null, 2),
    language,
  });
  await vscode.window.showTextDocument(doc);
}

function formatBridgeError(result: { error?: { message?: string }; stderr?: string }): string {
  return result.error?.message ?? result.stderr ?? 'Unknown error';
}

function parseSupercellDims(value: string): [number, number, number] | undefined {
  const values = value
    .trim()
    .split(/[,\s]+/)
    .map(part => Number.parseInt(part, 10));
  if (values.length !== 3 || values.some(value => !Number.isInteger(value) || value < 1)) {
    return undefined;
  }
  return values as [number, number, number];
}

function toTrajectoryStructure(frames: unknown[], filename: string): OpenQCStructure {
  const normalizedFrames = frames.map((frame, index) => normalizeTrajectoryFrame(frame, index));
  const firstFrame = normalizedFrames[0];

  return {
    ...firstFrame,
    kind: 'trajectory',
    name: firstFrame.name || `${filename} trajectory`,
    frames: normalizedFrames,
  };
}

function normalizeTrajectoryFrame(frame: unknown, index: number): OpenQCStructure {
  const source = (frame && typeof frame === 'object' ? frame : {}) as Partial<OpenQCStructure>;
  return {
    schemaVersion: 'openqc.structure.v1',
    kind: source.kind === 'periodic' || source.kind === 'surface' ? source.kind : 'molecule',
    name: source.name || `Frame ${index + 1}`,
    atoms: Array.isArray(source.atoms) ? source.atoms : [],
    ...(source.cell ? { cell: source.cell } : {}),
    ...(Array.isArray(source.bonds) ? { bonds: source.bonds } : {}),
    ...(source.metadata ? { metadata: source.metadata } : {}),
  };
}

function renderScfConvergenceHtml(results: OpenQCResults, energies: number[]): string {
  const points = toPolylinePoints(energies, 760, 320);
  const min = Math.min(...energies);
  const max = Math.max(...energies);
  const title = escapeHtml(results.software || 'OpenQC');
  const finalEnergy = results.finalEnergy
    ? `${results.finalEnergy.value} ${results.finalEnergy.unit}`
    : `${energies[energies.length - 1]}`;

  return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline';">
  <title>OpenQC SCF Convergence</title>
  <style>
    body { margin: 0; padding: 16px; font-family: var(--vscode-font-family); color: var(--vscode-foreground); background: var(--vscode-editor-background); }
    h1 { font-size: 16px; margin: 0 0 12px; }
    .meta { color: var(--vscode-descriptionForeground); margin-bottom: 12px; }
    svg { width: 100%; max-width: 820px; height: auto; background: var(--vscode-editorWidget-background); border: 1px solid var(--vscode-panel-border); }
    .axis { stroke: var(--vscode-descriptionForeground); stroke-width: 1; }
    .line { fill: none; stroke: var(--vscode-charts-blue); stroke-width: 2; }
    .point { fill: var(--vscode-charts-blue); }
  </style>
</head>
<body>
  <h1>${title} Energy Convergence</h1>
  <div class="meta">Steps: ${energies.length} · Min: ${formatNumber(min)} · Max: ${formatNumber(max)} · Final: ${escapeHtml(finalEnergy)}</div>
  <svg viewBox="0 0 820 380" role="img" aria-label="SCF convergence plot">
    <line class="axis" x1="40" y1="340" x2="800" y2="340"></line>
    <line class="axis" x1="40" y1="20" x2="40" y2="340"></line>
    <polyline class="line" points="${points}"></polyline>
    ${points
      .split(' ')
      .map(point => {
        const [x, y] = point.split(',');
        return `<circle class="point" cx="${x}" cy="${y}" r="2"></circle>`;
      })
      .join('')}
  </svg>
</body>
</html>`;
}

function toPolylinePoints(values: number[], width: number, height: number): string {
  const min = Math.min(...values);
  const max = Math.max(...values);
  const range = max - min || 1;
  return values
    .map((value, index) => {
      const x = 40 + (values.length === 1 ? 0 : (index / (values.length - 1)) * width);
      const y = 20 + (1 - (value - min) / range) * height;
      return `${formatNumber(x)},${formatNumber(y)}`;
    })
    .join(' ');
}

function formatNumber(value: number): string {
  return Number.isFinite(value) ? value.toFixed(6).replace(/\.?0+$/, '') : String(value);
}

function escapeHtml(value: string): string {
  return value.replace(/[&<>"']/g, char => {
    const escapes: Record<string, string> = {
      '&': '&amp;',
      '<': '&lt;',
      '>': '&gt;',
      '"': '&quot;',
      "'": '&#39;',
    };
    return escapes[char];
  });
}
