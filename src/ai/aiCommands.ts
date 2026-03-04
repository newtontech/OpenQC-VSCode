/**
 * AI Commands
 *
 * VSCode commands for AI-powered features including
 * input optimization, generation, explanation, and debugging.
 */

import * as vscode from 'vscode';
import * as fs from 'fs';
import { AICore, AICoreFactory, AIResponse } from './AICore';

/**
 * Register AI commands
 */
export function registerAICommands(context: vscode.ExtensionContext): void {
  const aiCore = AICoreFactory.getInstance(context);

  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.aiOptimizeInput', async () => {
      await aiOptimizeInput(aiCore);
    })
  );

  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.aiGenerateInput', async () => {
      await aiGenerateInput(aiCore);
    })
  );

  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.aiExplainParameters', async () => {
      await aiExplainParameters(aiCore);
    })
  );

  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.aiDebugCalculation', async () => {
      await aiDebugCalculation(aiCore);
    })
  );

  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.aiSettings', async () => {
      await aiSettings(aiCore);
    })
  );

  context.subscriptions.push(
    vscode.workspace.onDidChangeConfiguration(event => {
      if (event.affectsConfiguration('openqc.ai')) {
        aiCore.refreshConfig();
        vscode.window.showInformationMessage('AI configuration updated');
      }
    })
  );
}

async function aiOptimizeInput(aiCore: AICore): Promise<void> {
  if (!checkAIEnabled(aiCore)) {
    return;
  }

  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const document = editor.document;
  const inputContent = document.getText();
  const software = detectSoftware(document.fileName);

  if (!software) {
    vscode.window.showErrorMessage('Unsupported file type for AI optimization');
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'AI optimizing input file...',
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Analyzing input file...' });

      const response = await aiCore.optimizeInput(inputContent, software);

      if (!response.success) {
        vscode.window.showErrorMessage(`AI optimization failed: ${response.error}`);
        return;
      }

      progress.report({ message: 'Processing suggestions...' });
      displayOptimizationResults(response, document);
    }
  );
}

async function aiGenerateInput(aiCore: AICore): Promise<void> {
  if (!checkAIEnabled(aiCore)) {
    return;
  }

  const description = await vscode.window.showInputBox({
    prompt: 'Describe the calculation you want to perform',
    placeHolder: 'e.g., Geometry optimization of water molecule using DFT with PBE functional',
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a description';
      }
      return null;
    },
  });

  if (!description) {
    return;
  }

  const software = await selectSoftware();
  if (!software) {
    return;
  }

  const structureUri = await vscode.window.showOpenDialog({
    canSelectFiles: true,
    canSelectFolders: false,
    canSelectMany: false,
    filters: {
      'Structure Files': ['xyz', 'pdb', 'cif', 'poscar', 'contcar'],
      'All Files': ['*'],
    },
    openLabel: 'Select Structure (Optional)',
  });

  let ctx: Record<string, any> | undefined;
  if (structureUri && structureUri[0]) {
    try {
      const structureContent = fs.readFileSync(structureUri[0].fsPath, 'utf-8');
      ctx = { structure: structureContent, structureFormat: detectFormat(structureUri[0].fsPath) };
    } catch (error) {
      vscode.window.showWarningMessage(`Failed to read structure file: ${error}`);
    }
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'AI generating input file...',
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Generating input...' });

      const response = await aiCore.generateInput(description, software, ctx);

      if (!response.success) {
        vscode.window.showErrorMessage(`AI generation failed: ${response.error}`);
        return;
      }

      progress.report({ message: 'Opening generated file...' });

      if (response.generatedInput) {
        const doc = await vscode.workspace.openTextDocument({
          content: response.generatedInput,
          language: getLanguageId(software),
        });
        await vscode.window.showTextDocument(doc);
        vscode.window.showInformationMessage('AI-generated input file created');
      }
    }
  );
}

async function aiExplainParameters(aiCore: AICore): Promise<void> {
  if (!checkAIEnabled(aiCore)) {
    return;
  }

  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const document = editor.document;
  const inputContent = document.getText();
  const software = detectSoftware(document.fileName);

  if (!software) {
    vscode.window.showErrorMessage('Unsupported file type for parameter explanation');
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'AI analyzing parameters...',
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Analyzing parameters...' });

      const response = await aiCore.explainParameters(inputContent, software);

      if (!response.success) {
        vscode.window.showErrorMessage(`AI explanation failed: ${response.error}`);
        return;
      }

      progress.report({ message: 'Generating explanation...' });
      displayExplanation(response, software);
    }
  );
}

async function aiDebugCalculation(aiCore: AICore): Promise<void> {
  if (!checkAIEnabled(aiCore)) {
    return;
  }

  const inputUri = await vscode.window.showOpenDialog({
    canSelectFiles: true,
    canSelectFolders: false,
    canSelectMany: false,
    filters: {
      'Input Files': ['inp', 'com', 'gjf', 'in', 'nw'],
      'All Files': ['*'],
    },
    openLabel: 'Select Input File',
  });

  if (!inputUri || !inputUri[0]) {
    return;
  }

  const outputUri = await vscode.window.showOpenDialog({
    canSelectFiles: true,
    canSelectFolders: false,
    canSelectMany: false,
    filters: {
      'Output Files': ['out', 'log', 'stdout'],
      'All Files': ['*'],
    },
    openLabel: 'Select Output/Log File',
  });

  if (!outputUri || !outputUri[0]) {
    return;
  }

  let inputContent: string;
  let outputContent: string;

  try {
    inputContent = fs.readFileSync(inputUri[0].fsPath, 'utf-8');
    outputContent = fs.readFileSync(outputUri[0].fsPath, 'utf-8');
  } catch (error) {
    vscode.window.showErrorMessage(`Failed to read files: ${error}`);
    return;
  }

  const software = detectSoftware(inputUri[0].fsPath);
  if (!software) {
    vscode.window.showErrorMessage('Could not detect software type from input file');
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'AI debugging calculation...',
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Analyzing input and output...' });

      const response = await aiCore.debugCalculation(inputContent, outputContent, software);

      if (!response.success) {
        vscode.window.showErrorMessage(`AI debugging failed: ${response.error}`);
        return;
      }

      progress.report({ message: 'Generating debug report...' });
      displayDebugResults(response, inputUri[0].fsPath, outputUri[0].fsPath);
    }
  );
}

async function aiSettings(aiCore: AICore): Promise<void> {
  const config = aiCore.getConfig();
  const validation = aiCore.validateConfig();

  const items: vscode.QuickPickItem[] = [
    {
      label: '$(gear) Configure AI Settings',
      description: 'Open VSCode settings for AI configuration',
    },
    {
      label: '$(info) AI Status',
      description: `Enabled: ${config.enabled ? 'Yes' : 'No'} | Provider: ${config.provider}`,
    },
    {
      label: '$(symbol-event) Test AI Connection',
      description: 'Test connection to AI backend',
    },
  ];

  const selected = await vscode.window.showQuickPick(items, {
    placeHolder: 'Select AI option',
  });

  if (!selected) {
    return;
  }

  if (selected.label.includes('Configure')) {
    vscode.commands.executeCommand('workbench.action.openSettings', 'openqc.ai');
  } else if (selected.label.includes('Test')) {
    const available = await aiCore.isAvailable();
    if (available) {
      vscode.window.showInformationMessage('AI backend is available and ready');
    } else {
      vscode.window.showErrorMessage(`AI backend is not available. ${validation.errors.join(', ')}`);
    }
  }
}

function checkAIEnabled(aiCore: AICore): boolean {
  if (!aiCore.isEnabled()) {
    vscode.window
      .showWarningMessage(
        'AI features are disabled. Enable them in settings?',
        'Open Settings',
        'Cancel'
      )
      .then(selection => {
        if (selection === 'Open Settings') {
          vscode.commands.executeCommand('workbench.action.openSettings', 'openqc.ai');
        }
      });
    return false;
  }
  return true;
}

function detectSoftware(filename: string): string | null {
  const basename = filename.toLowerCase();
  
  if (basename.includes('incar') || basename.includes('poscar') || basename.includes('vasp')) {
    return 'vasp';
  }
  if (basename.endsWith('.inp')) {
    return 'cp2k';
  }
  if (basename.endsWith('.com') || basename.endsWith('.gjf')) {
    return 'gaussian';
  }
  if (basename.endsWith('.in')) {
    return 'qe';
  }
  if (basename.endsWith('.nw') || basename.endsWith('.nwinp')) {
    return 'nwchem';
  }
  
  return null;
}

function detectFormat(filename: string): string {
  const ext = filename.split('.').pop()?.toLowerCase() || '';
  
  if (['xyz', 'pdb', 'cif'].includes(ext)) {
    return ext;
  }
  if (filename.toLowerCase().includes('poscar') || filename.toLowerCase().includes('contcar')) {
    return 'vasp';
  }
  
  return 'xyz';
}

async function selectSoftware(): Promise<string | undefined> {
  const options = [
    { label: 'VASP', description: 'Vienna Ab initio Simulation Package' },
    { label: 'CP2K', description: 'Car-Parrinello 2K' },
    { label: 'Quantum ESPRESSO', description: 'Plane-wave DFT' },
    { label: 'Gaussian', description: 'Gaussian quantum chemistry' },
    { label: 'ORCA', description: 'ORCA quantum chemistry' },
    { label: 'NWChem', description: 'NWChem computational chemistry' },
    { label: 'GAMESS', description: 'GAMESS quantum chemistry' },
  ];

  const selected = await vscode.window.showQuickPick(options, {
    placeHolder: 'Select quantum chemistry software',
  });

  if (selected) {
    return selected.label.toLowerCase().replace(/\s+/g, '_');
  }

  return undefined;
}

function getLanguageId(software: string): string {
  const mapping: Record<string, string> = {
    vasp: 'vasp',
    cp2k: 'cp2k',
    quantum_espresso: 'qe',
    gaussian: 'gaussian',
    orca: 'orca',
    nwchem: 'nwchem',
    gamess: 'gamess',
  };
  
  return mapping[software] || 'plaintext';
}

function displayOptimizationResults(response: AIResponse, document: vscode.TextDocument): void {
  if (!response.suggestions || response.suggestions.length === 0) {
    vscode.window.showInformationMessage('No optimization suggestions found');
    return;
  }

  const panel = vscode.window.createWebviewPanel(
    'openqc.aiOptimization',
    'AI Optimization Results',
    vscode.ViewColumn.Two,
    { enableScripts: true }
  );

  const suggestionsHtml = response.suggestions
    .map(
      (s, i) => `
    <div class="suggestion ${s.type}">
      <div class="suggestion-header">
        <span class="badge ${s.type}">${s.type.toUpperCase()}</span>
        <span class="parameter">${s.parameter || 'General'}</span>
      </div>
      <div class="message">${s.message}</div>
      <div class="values">
        <div class="current"><strong>Current:</strong> ${s.currentValue || 'N/A'}</div>
        <div class="suggested"><strong>Suggested:</strong> ${s.suggestedValue || 'N/A'}</div>
      </div>
      <div class="explanation">${s.explanation || ''}</div>
    </div>
  `
    )
    .join('');

  panel.webview.html = `
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="UTF-8">
      <title>AI Optimization Results</title>
      <style>
        body { font-family: -apple-system, sans-serif; padding: 20px; background: #1e1e1e; color: #ccc; }
        h2 { margin-top: 0; color: #fff; }
        .suggestion { background: #252526; border: 1px solid #3c3c3c; border-radius: 5px; padding: 15px; margin-bottom: 15px; }
        .suggestion-header { display: flex; align-items: center; margin-bottom: 10px; }
        .badge { padding: 2px 8px; border-radius: 3px; font-size: 11px; margin-right: 10px; }
        .badge.optimization { background: #0e639c; color: white; }
        .badge.warning { background: #f66; color: white; }
        .badge.info { background: #3c9; color: white; }
        .parameter { font-weight: bold; color: #9cdcfe; }
        .values { display: flex; gap: 20px; margin: 10px 0; }
        .current { color: #ce9178; }
        .suggested { color: #b5cea8; }
        .explanation { color: #9cdcfe; font-size: 12px; margin-top: 10px; }
      </style>
    </head>
    <body>
      <h2>AI Optimization Results: ${document.fileName}</h2>
      <div class="suggestions">${suggestionsHtml}</div>
    </body>
    </html>
  `;
}

function displayExplanation(response: AIResponse, software: string): void {
  const panel = vscode.window.createWebviewPanel(
    'openqc.aiExplanation',
    'AI Parameter Explanation',
    vscode.ViewColumn.Two,
    { enableScripts: true }
  );

  panel.webview.html = `
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="UTF-8">
      <title>AI Explanation</title>
      <style>
        body { font-family: -apple-system, sans-serif; padding: 20px; background: #1e1e1e; color: #ccc; }
        h2 { margin-top: 0; color: #fff; }
        .content { line-height: 1.6; }
        .content p { margin-bottom: 15px; }
        .content code { background: #3c3c3c; padding: 2px 5px; border-radius: 3px; }
        .content pre { background: #252526; padding: 15px; border-radius: 5px; overflow-x: auto; }
      </style>
    </head>
    <body>
      <h2>AI Parameter Explanation: ${software.toUpperCase()}</h2>
      <div class="content"><pre>${response.content || 'No explanation available'}</pre></div>
    </body>
    </html>
  `;
}

function displayDebugResults(response: AIResponse, inputFile: string, outputFile: string): void {
  const panel = vscode.window.createWebviewPanel(
    'openqc.aiDebug',
    'AI Debug Results',
    vscode.ViewColumn.Two,
    { enableScripts: true }
  );

  panel.webview.html = `
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="UTF-8">
      <title>AI Debug Results</title>
      <style>
        body { font-family: -apple-system, sans-serif; padding: 20px; background: #1e1e1e; color: #ccc; }
        h2 { margin-top: 0; color: #fff; }
        .files { margin-bottom: 20px; padding: 10px; background: #252526; border-radius: 5px; }
        .content { line-height: 1.6; }
        .content pre { background: #252526; padding: 15px; border-radius: 5px; overflow-x: auto; white-space: pre-wrap; }
        .error { color: #f66; }
        .success { color: #3c9; }
      </style>
    </head>
    <body>
      <h2>AI Debug Analysis</h2>
      <div class="files">
        <div><strong>Input:</strong> ${inputFile}</div>
        <div><strong>Output:</strong> ${outputFile}</div>
      </div>
      <div class="content"><pre>${response.content || 'No analysis available'}</pre></div>
    </body>
    </html>
  `;
}
