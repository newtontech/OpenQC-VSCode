/**
 * Format Converter - Quantum chemistry file format conversion adapter
 *
 * This module provides a TypeScript interface to the Python dpdata-based
 * format conversion backend. It supports conversion between VASP, Gaussian,
 * ORCA, XYZ, PDB, CIF, and other quantum chemistry file formats.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { spawn } from 'child_process';

/**
 * Supported file formats for conversion
 */
export enum SupportedFormat {
  VASP = 'vasp',
  POSCAR = 'poscar',
  CONTCAR = 'contcar',
  Gaussian = 'gaussian',
  GJF = 'gjf',
  COM = 'com',
  ORCA = 'orca',
  INP = 'inp',
  XYZ = 'xyz',
  PDB = 'pdb',
  CIF = 'cif',
  CP2K = 'cp2k',
  QE = 'qe',
  QuantumEspresso = 'quantum_espresso',
}

/**
 * Conversion metadata
 */
export interface ConversionMetadata {
  natoms?: number;
  nbonds?: number;
  elements?: string[];
  atom_counts?: number[];
  frames_count?: number;
  original_format?: string;
}

/**
 * Result of a format conversion operation
 */
export interface ConversionResult {
  success: boolean;
  input_format?: string;
  output_format?: string;
  atoms_count?: number;
  frames_count?: number;
  metadata?: ConversionMetadata;
  output_file?: string;
  error?: string;
  error_type?: string;
}

/**
 * Batch conversion result
 */
export interface BatchConversionResult {
  success: boolean;
  total: number;
  successful: number;
  failed: number;
  results: Array<{
    input_file: string;
    output_file: string;
    result: ConversionResult;
  }>;
}

/**
 * Format detection result
 */
export interface FormatDetection {
  format: SupportedFormat;
  confidence: number;
  extension: string;
}

/**
 * Configuration for the FormatConverter
 */
export interface FormatConverterConfig {
  pythonPath?: string;
  scriptPath?: string;
  preserveMetadata?: boolean;
}

/**
 * Main format converter class
 */
export class FormatConverter {
  private config: FormatConverterConfig;
  private outputChannel: vscode.OutputChannel;

  constructor(config: FormatConverterConfig = {}) {
    this.config = {
      pythonPath: config.pythonPath || 'python3',
      scriptPath:
        config.scriptPath || path.resolve(__dirname, '../../../python/format_converter.py'),
      preserveMetadata: config.preserveMetadata ?? true,
    };
    this.outputChannel = vscode.window.createOutputChannel('OpenQC Format Converter');
  }

  /**
   * Detect the format of a file from its extension and content
   */
  static detectFormat(filePath: string): FormatDetection {
    const ext = path.extname(filePath).toLowerCase();
    const basename = path.basename(filePath).toLowerCase();

    // Special VASP files
    if (basename === 'poscar') {
      return { format: SupportedFormat.POSCAR, confidence: 1.0, extension: '' };
    }
    if (basename === 'contcar') {
      return { format: SupportedFormat.CONTCAR, confidence: 1.0, extension: '' };
    }
    if (basename === 'incar') {
      return { format: SupportedFormat.VASP, confidence: 1.0, extension: '' };
    }
    if (basename === 'kpoints') {
      return { format: SupportedFormat.VASP, confidence: 1.0, extension: '' };
    }
    if (basename === 'potcar') {
      return { format: SupportedFormat.VASP, confidence: 1.0, extension: '' };
    }

    // Extension-based detection
    const formatMap: Record<string, SupportedFormat> = {
      '.xyz': SupportedFormat.XYZ,
      '.pdb': SupportedFormat.PDB,
      '.cif': SupportedFormat.CIF,
      '.gjf': SupportedFormat.GJF,
      '.com': SupportedFormat.COM,
      '.inp': SupportedFormat.INP,
    };

    const format = formatMap[ext];
    if (format) {
      return { format, confidence: 0.9, extension: ext };
    }

    // Default to VASP
    return { format: SupportedFormat.POSCAR, confidence: 0.5, extension: ext };
  }

  /**
   * Convert a single file from one format to another
   *
   * @param inputPath Path to input file
   * @param outputPath Path to output file
   * @param fromFormat Source format (auto-detected if not specified)
   * @param toFormat Target format (inferred from output path if not specified)
   */
  async convert(
    inputPath: string,
    outputPath: string,
    fromFormat?: SupportedFormat,
    toFormat?: SupportedFormat
  ): Promise<ConversionResult> {
    try {
      this.log(`Converting ${inputPath} to ${outputPath}`);

      const detection = fromFormat ? undefined : FormatConverter.detectFormat(inputPath);
      const sourceFormat = fromFormat || detection?.format;
      const targetFormat = toFormat || FormatConverter.detectFormat(outputPath).format;

      const args = [this.config.scriptPath!, inputPath, outputPath];

      if (fromFormat) {
        args.push('--from', fromFormat);
      }
      if (toFormat) {
        args.push('--to', toFormat);
      }
      if (!this.config.preserveMetadata) {
        args.push('--no-metadata');
      }
      args.push('--json');

      const stdout = await this.runPython(args);
      const result = JSON.parse(stdout) as ConversionResult;

      if (result.success) {
        this.log(`Conversion successful: ${result.output_file}`);
      } else {
        this.logError(`Conversion failed: ${result.error}`);
      }

      return result;
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : String(error);
      this.logError(`Conversion error: ${errorMessage}`);
      return {
        success: false,
        error: errorMessage,
        error_type: 'ExecutionError',
      };
    }
  }

  /**
   * Batch convert multiple files
   *
   * @param inputPaths Array of input file paths
   * @param outputDir Output directory
   * @param toFormat Target format
   * @param fromFormat Source format (auto-detected if not specified)
   */
  async batchConvert(
    inputPaths: string[],
    outputDir: string,
    toFormat: SupportedFormat,
    fromFormat?: SupportedFormat
  ): Promise<BatchConversionResult> {
    try {
      this.log(`Batch converting ${inputPaths.length} files to ${toFormat}`);

      const args = [
        this.config.scriptPath!,
        '--batch',
        ...inputPaths,
        '--to',
        toFormat,
        '--output-dir',
        outputDir,
        '--json',
      ];

      if (fromFormat) {
        args.push('--from', fromFormat);
      }
      if (!this.config.preserveMetadata) {
        args.push('--no-metadata');
      }

      const stdout = await this.runPython(args);
      const result = JSON.parse(stdout) as BatchConversionResult;

      this.log(`Batch conversion complete: ${result.successful}/${result.total} successful`);

      return result;
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : String(error);
      this.logError(`Batch conversion error: ${errorMessage}`);
      return {
        success: false,
        total: inputPaths.length,
        successful: 0,
        failed: inputPaths.length,
        results: [],
      };
    }
  }

  /**
   * Convert the current editor document to a different format
   */
  async convertCurrentDocument(
    targetFormat: SupportedFormat
  ): Promise<ConversionResult | undefined> {
    const editor = vscode.window.activeTextEditor;
    if (!editor) {
      vscode.window.showErrorMessage('No active document');
      return undefined;
    }

    const inputPath = editor.document.uri.fsPath;
    const inputDir = path.dirname(inputPath);
    const inputBasename = path.basename(inputPath, path.extname(inputPath));
    const extension = this.getExtensionForFormat(targetFormat);
    const outputPath = path.join(inputDir, `${inputBasename}.${extension}`);

    return this.convert(inputPath, outputPath, undefined, targetFormat);
  }

  /**
   * Get file extension for a given format
   */
  getExtensionForFormat(format: SupportedFormat): string {
    const extMap: Record<SupportedFormat, string> = {
      [SupportedFormat.VASP]: 'POSCAR',
      [SupportedFormat.POSCAR]: 'POSCAR',
      [SupportedFormat.CONTCAR]: 'CONTCAR',
      [SupportedFormat.Gaussian]: 'gjf',
      [SupportedFormat.GJF]: 'gjf',
      [SupportedFormat.COM]: 'com',
      [SupportedFormat.ORCA]: 'inp',
      [SupportedFormat.INP]: 'inp',
      [SupportedFormat.XYZ]: 'xyz',
      [SupportedFormat.PDB]: 'pdb',
      [SupportedFormat.CIF]: 'cif',
      [SupportedFormat.CP2K]: 'inp',
      [SupportedFormat.QE]: 'in',
      [SupportedFormat.QuantumEspresso]: 'in',
    };
    return extMap[format];
  }

  /**
   * Check if the Python backend is available
   */
  async checkBackend(): Promise<boolean> {
    try {
      await this.runPython(['--version']);
      await this.runPython(['-c', 'import dpdata']);
      return true;
    } catch {
      return false;
    }
  }

  /**
   * Show setup instructions if backend is not available
   */
  async showSetupInstructions(): Promise<void> {
    const action = await vscode.window.showWarningMessage(
      'OpenQC Format Converter requires Python and dpdata. Install now?',
      'Install',
      'Cancel'
    );

    if (action === 'Install') {
      const terminal = vscode.window.createTerminal('OpenQC Setup');
      terminal.sendText('pip install dpdata');
      terminal.show();
    }
  }

  /**
   * Quick conversion using the existing StructureConverter
   * (for XYZ output without Python backend)
   */
  static convertToXYZ(
    atoms: Array<{ elem: string; x: number; y: number; z: number }>,
    comment: string = 'molecule'
  ): string {
    const lines: string[] = [];

    lines.push(String(atoms.length));
    lines.push(comment);

    for (const atom of atoms) {
      const x = atom.x.toFixed(6);
      const y = atom.y.toFixed(6);
      const z = atom.z.toFixed(6);
      lines.push(`${atom.elem} ${x} ${y} ${z}`);
    }

    return lines.join('\n');
  }

  private log(message: string): void {
    this.outputChannel.appendLine(message);
  }

  private logError(message: string): void {
    this.outputChannel.appendLine(`[ERROR] ${message}`);
  }

  private runPython(args: string[]): Promise<string> {
    return new Promise((resolve, reject) => {
      const child = spawn(this.config.pythonPath!, args, { shell: false });
      let stdout = '';
      let stderr = '';

      child.stdout.on('data', chunk => {
        stdout += chunk.toString();
      });

      child.stderr.on('data', chunk => {
        stderr += chunk.toString();
      });

      child.on('error', reject);
      child.on('close', code => {
        if (code === 0) {
          resolve(stdout);
          return;
        }

        reject(new Error(stderr.trim() || `Python process exited with code ${code}`));
      });
    });
  }

  dispose(): void {
    this.outputChannel.dispose();
  }
}

/**
 * Helper function to create and show a quick format conversion
 */
export async function quickConvert(targetFormat: SupportedFormat): Promise<void> {
  const converter = new FormatConverter();

  if (!(await converter.checkBackend())) {
    await converter.showSetupInstructions();
    return;
  }

  const result = await converter.convertCurrentDocument(targetFormat);

  if (result?.success) {
    const action = await vscode.window.showInformationMessage(
      `Successfully converted to ${targetFormat.toUpperCase()}`,
      'Open File',
      'Copy to Clipboard'
    );

    if (action === 'Open File' && result.output_file) {
      const doc = await vscode.workspace.openTextDocument(result.output_file);
      await vscode.window.showTextDocument(doc);
    } else if (action === 'Copy to Clipboard') {
      // Read and copy file content
      const content = await vscode.workspace.fs.readFile(vscode.Uri.file(result.output_file!));
      await vscode.env.clipboard.writeText(Buffer.from(content).toString('utf8'));
    }
  } else {
    vscode.window.showErrorMessage(`Conversion failed: ${result?.error}`);
  }
}
