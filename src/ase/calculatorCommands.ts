/**
 * ASE Calculator Commands
 *
 * VSCode commands for ASE calculator operations including
 * input generation, calculation execution, and result reading.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import * as fs from 'fs';
import {
  ASECalculator,
  CalculatorFactory,
  CalculatorType,
  CalculatorConfig,
  getCalculatorConfiguration,
} from './ASECalculator';
import { ASEConverter, ASEAtoms } from './ASEConverter';

/**
 * Register calculator commands
 */
export function registerCalculatorCommands(context: vscode.ExtensionContext): void {
  const factory = new CalculatorFactory(context);

  // Command: Generate calculator input
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.generateCalculatorInput',
      async () => await generateCalculatorInput(factory)
    )
  );

  // Command: Run calculation
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.runCalculation',
      async () => await runCalculation(factory)
    )
  );

  // Command: Read calculation results
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.readResults', async () => await readResults())
  );

  // Command: Quick VASP calculation
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.quickVASP',
      async () => await quickCalculation(factory, CalculatorType.VASP)
    )
  );

  // Command: Quick CP2K calculation
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.quickCP2K',
      async () => await quickCalculation(factory, CalculatorType.CP2K)
    )
  );

  // Command: Quick QE calculation
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.quickQE',
      async () => await quickCalculation(factory, CalculatorType.QE)
    )
  );

  // Command: Check calculator availability
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.checkCalculators',
      async () => await checkCalculators(factory)
    )
  );
}

/**
 * Generate input files for a calculator
 */
async function generateCalculatorInput(factory: CalculatorFactory): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  // Select calculator type
  const calculatorType = await selectCalculatorType();
  if (!calculatorType) {
    return;
  }

  // Get output directory
  const outputDir = await vscode.window.showInputBox({
    prompt: 'Enter output directory for input files',
    value: path.dirname(editor.document.uri.fsPath),
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a directory path';
      }
      return null;
    },
  });

  if (!outputDir) {
    return;
  }

  // Create output directory if it doesn't exist
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }

  // Get calculator configuration
  const config = await configureCalculator(calculatorType);
  if (!config) {
    return;
  }

  const calculator = factory.createCalculator(config);

  // Show progress
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Generating ${calculatorType} input`,
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Reading structure...' });

      // Read structure from current file
      const converter = new ASEConverter(
        vscode.extensions.getExtension('newtontech.openqc')!.exports.context
      );
      const filepath = editor.document.uri.fsPath;
      const readResult = await converter.readToAtoms(filepath);

      if (!readResult.success || !readResult.atoms) {
        vscode.window.showErrorMessage(`Failed to read structure: ${readResult.error}`);
        return;
      }

      progress.report({ message: 'Generating input files...' });

      // Generate input
      const result = await calculator.generateInput(readResult.atoms, outputDir);

      if (!result.success) {
        vscode.window.showErrorMessage(`Input generation failed: ${result.error}`);
        return;
      }

      progress.report({ message: 'Opening generated files...' });

      // Open the first generated file
      const files = Object.values(result.files);
      if (files.length > 0) {
        const doc = await vscode.workspace.openTextDocument(files[0]);
        await vscode.window.showTextDocument(doc);
      }

      vscode.window.showInformationMessage(
        `Generated ${Object.keys(result.files).length} input file(s) in ${outputDir}`
      );
    }
  );
}

/**
 * Run a calculation
 */
async function runCalculation(factory: CalculatorFactory): Promise<void> {
  // Select input directory
  const inputDir = await vscode.window.showInputBox({
    prompt: 'Enter input directory containing calculation files',
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a directory path';
      }
      if (!fs.existsSync(value)) {
        return 'Directory does not exist';
      }
      return null;
    },
  });

  if (!inputDir) {
    return;
  }

  // Select calculator type
  const calculatorType = await selectCalculatorType();
  if (!calculatorType) {
    return;
  }

  // Get calculator configuration
  const config = await configureCalculator(calculatorType);
  if (!config) {
    return;
  }

  const calculator = factory.createCalculator(config);

  // Run calculation
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Running ${calculatorType} calculation`,
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Starting calculation...' });

      const result = await calculator.runCalculation(inputDir);

      if (!result.success) {
        vscode.window.showErrorMessage(`Calculation failed: ${result.error}`);
        return;
      }

      progress.report({ message: 'Processing results...' });

      // Display results
      displayCalculationResults(result);
    }
  );
}

/**
 * Read calculation results
 */
async function readResults(): Promise<void> {
  // Select output directory
  const outputDir = await vscode.window.showInputBox({
    prompt: 'Enter output directory containing calculation results',
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a directory path';
      }
      if (!fs.existsSync(value)) {
        return 'Directory does not exist';
      }
      return null;
    },
  });

  if (!outputDir) {
    return;
  }

  // Select calculator type
  const calculatorType = await selectCalculatorType();
  if (!calculatorType) {
    return;
  }

  // Create temporary calculator to read results
  const extension = vscode.extensions.getExtension('newtontech.openqc');
  if (!extension) {
    vscode.window.showErrorMessage('OpenQC extension not found');
    return;
  }

  const { CalculatorFactory } = await import('./ASECalculator');
  const factory = new CalculatorFactory(extension.exports.context);
  const config: CalculatorConfig = {
    type: calculatorType,
    parameters: {},
  };
  const calculator = factory.createCalculator(config);

  // Read results
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Reading ${calculatorType} results`,
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Reading output files...' });

      const result = await calculator.readResults(outputDir);

      if (!result.success) {
        vscode.window.showErrorMessage(`Failed to read results: ${result.error}`);
        return;
      }

      progress.report({ message: 'Displaying results...' });

      displayCalculationResults(result);
    }
  );
}

/**
 * Quick calculation from current file
 */
async function quickCalculation(
  factory: CalculatorFactory,
  calculatorType: CalculatorType
): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  // Get working directory
  const filepath = editor.document.uri.fsPath;
  const defaultDir = path.join(path.dirname(filepath), `${calculatorType}_calc`);

  const workingDir = await vscode.window.showInputBox({
    prompt: 'Enter working directory for calculation',
    value: defaultDir,
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a directory path';
      }
      return null;
    },
  });

  if (!workingDir) {
    return;
  }

  // Create working directory
  if (!fs.existsSync(workingDir)) {
    fs.mkdirSync(workingDir, { recursive: true });
  }

  // Get calculator configuration
  const config = await configureCalculator(calculatorType);
  if (!config) {
    return;
  }

  const calculator = factory.createCalculator(config);

  // Run calculation
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Running ${calculatorType} calculation`,
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Reading structure...' });

      // Read structure
      const converter = new ASEConverter(
        vscode.extensions.getExtension('newtontech.openqc')!.exports.context
      );
      const readResult = await converter.readToAtoms(filepath);

      if (!readResult.success || !readResult.atoms) {
        vscode.window.showErrorMessage(`Failed to read structure: ${readResult.error}`);
        return;
      }

      progress.report({ message: 'Running calculation...' });

      // Run calculation
      const result = await calculator.calculate(readResult.atoms, workingDir);

      if (!result.success) {
        vscode.window.showErrorMessage(`Calculation failed: ${result.error}`);
        return;
      }

      progress.report({ message: 'Processing results...' });

      displayCalculationResults(result);
    }
  );
}

/**
 * Check calculator availability
 */
async function checkCalculators(factory: CalculatorFactory): Promise<void> {
  const config = getCalculatorConfiguration();

  const results: string[] = [];

  for (const [type, calcConfig] of Object.entries(config)) {
    const calculator = factory.createCalculator({
      type: type as CalculatorType,
      command: calcConfig.command,
      parameters: {},
    });

    const backendAvailable = await calculator.isAvailable();
    const calcAvailable = await calculator.isCalculatorAvailable();

    results.push(
      `${type.toUpperCase()}: Backend ${backendAvailable ? '✓' : '✗'}, Calculator ${
        calcAvailable ? '✓' : '✗'
      }`
    );
  }

  vscode.window.showInformationMessage(`Calculator Status: ${results.join(', ')}`, { modal: true });
}

/**
 * Select calculator type
 */
async function selectCalculatorType(): Promise<CalculatorType | undefined> {
  const types = [
    {
      label: 'VASP',
      type: CalculatorType.VASP,
      description: 'Vienna Ab initio Simulation Package',
    },
    { label: 'CP2K', type: CalculatorType.CP2K, description: 'Car-Parrinello 2K' },
    { label: 'Quantum ESPRESSO', type: CalculatorType.QE, description: 'Plane-wave DFT' },
  ];

  const selected = await vscode.window.showQuickPick(types, {
    placeHolder: 'Select calculator type',
  });

  return selected?.type;
}

/**
 * Configure calculator parameters
 */
async function configureCalculator(
  calculatorType: CalculatorType
): Promise<CalculatorConfig | undefined> {
  const config = getCalculatorConfiguration();

  const defaults: CalculatorConfig = {
    type: calculatorType,
    command: config[calculatorType]?.command,
    parameters: config[calculatorType]?.defaultParams || {},
  };

  // Ask for custom parameters (optional)
  const customize = await vscode.window.showQuickPick(
    [
      { label: 'Use defaults', value: false },
      { label: 'Customize parameters', value: true },
    ],
    { placeHolder: 'Calculator configuration' }
  );

  if (!customize) {
    return undefined;
  }

  if (!customize.value) {
    return defaults;
  }

  // TODO: Implement parameter customization UI
  // For now, return defaults
  vscode.window.showInformationMessage('Using default parameters (customization coming soon)');

  return defaults;
}

/**
 * Display calculation results
 */
function displayCalculationResults(result: {
  success: boolean;
  energy?: number;
  forces?: number[][];
  stress?: number[][];
  atoms?: ASEAtoms;
  error?: string;
  warnings: string[];
  metadata: Record<string, any>;
  outputFiles: string[];
}): void {
  if (!result.success) {
    vscode.window.showErrorMessage(`Calculation failed: ${result.error}`);
    return;
  }

  const messages: string[] = [];

  if (result.energy !== undefined) {
    messages.push(`Energy: ${result.energy.toFixed(6)} eV`);
  }

  if (result.forces) {
    messages.push(`Forces: ${result.forces.length} atoms`);
  }

  if (result.stress) {
    messages.push('Stress tensor available');
  }

  if (result.warnings.length > 0) {
    messages.push(`Warnings: ${result.warnings.length}`);
  }

  vscode.window.showInformationMessage(messages.join(', '));

  // Show detailed results in a new document
  const resultContent = formatResults(result);
  vscode.workspace
    .openTextDocument({
      content: resultContent,
      language: 'json',
    })
    .then(doc => vscode.window.showTextDocument(doc));
}

/**
 * Format calculation results for display
 */
function formatResults(result: {
  success: boolean;
  energy?: number;
  forces?: number[][];
  stress?: number[][];
  atoms?: ASEAtoms;
  error?: string;
  warnings: string[];
  metadata: Record<string, any>;
  outputFiles: string[];
}): string {
  const output: Record<string, any> = {
    success: result.success,
    energy: result.energy,
    natoms: result.atoms?.chemical_symbols.length,
    formula: result.atoms?.chemical_symbols
      ? getChemicalFormula(result.atoms.chemical_symbols)
      : undefined,
    has_forces: !!result.forces,
    has_stress: !!result.stress,
    output_files: result.outputFiles,
    warnings: result.warnings,
    metadata: result.metadata,
  };

  if (result.forces) {
    output.max_force = Math.max(
      ...result.forces.map(f => Math.sqrt(f[0] ** 2 + f[1] ** 2 + f[2] ** 2))
    );
  }

  return JSON.stringify(output, null, 2);
}

/**
 * Get chemical formula from element symbols
 */
function getChemicalFormula(symbols: string[]): string {
  const counts: Record<string, number> = {};
  for (const symbol of symbols) {
    counts[symbol] = (counts[symbol] || 0) + 1;
  }

  return Object.entries(counts)
    .map(([element, count]) => `${element}${count > 1 ? count : ''}`)
    .join('');
}
