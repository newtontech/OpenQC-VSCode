/**
 * ASE Integration Commands
 *
 * VSCode commands for ASE format conversion and migration.
 */

import * as vscode from 'vscode';
import { ASEConverter, ASEFormat, type ASEAtoms, type ConversionResult } from './ASEConverter';

interface AtomsExtraction {
  atoms?: ASEAtoms;
  source?: string;
  error?: string;
}

type ValidationResult<T> = { value: T } | { error: string };

const ELEMENT_SYMBOLS = [
  '',
  ...'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split(
    ' '
  ),
];

const ELEMENT_SYMBOL_SET: ReadonlySet<string> = new Set(ELEMENT_SYMBOLS.slice(1));

/**
 * Register ASE commands
 */
export function registerASECommands(context: vscode.ExtensionContext): void {
  const converter = new ASEConverter(context);

  // Command: Convert to ASE Atoms
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.convertToASE',
      async () => await convertToASE(converter)
    )
  );

  // Command: Convert from ASE Atoms
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.convertFromASE',
      async () => await convertFromASE(converter)
    )
  );

  // Command: Migrate to target format
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.migrateFormat',
      async () => await migrateFormat(converter)
    )
  );

  // Command: Quick convert (VASP -> CP2K, etc.)
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.quickConvert',
      async (from?: string, to?: string) => await quickConvert(converter, from, to)
    )
  );
}

/**
 * Convert current file to ASE Atoms
 */
async function convertToASE(converter: ASEConverter): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const document = editor.document;
  const filepath = document.uri.fsPath;

  // Show progress
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Converting to ASE Atoms',
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Reading file...' });

      const result = await converter.readToAtoms(filepath);

      if (!result.success) {
        vscode.window.showErrorMessage(`Conversion failed: ${result.error}`);
        return;
      }

      progress.report({ message: 'Displaying results...' });

      // Display results in a new document
      const atomsInfo = formatAtomsInfo(result);
      const doc = await vscode.workspace.openTextDocument({
        content: atomsInfo,
        language: 'json',
      });
      await vscode.window.showTextDocument(doc);

      // Show success message
      const formula = result.atoms?.chemical_symbols
        ? getChemicalFormula(result.atoms.chemical_symbols)
        : 'Unknown';
      vscode.window.showInformationMessage(
        `Successfully converted to ASE Atoms: ${formula} (${result.metadata.natoms} atoms)`
      );
    }
  );
}

/**
 * Convert ASE Atoms to target format
 */
async function convertFromASE(converter: ASEConverter): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor with ASE Atoms JSON');
    return;
  }

  let atoms: any;
  try {
    const parsed = JSON.parse(editor.document.getText());
    const extracted = extractASEAtoms(parsed);
    if (!extracted.atoms) {
      vscode.window.showErrorMessage(
        `Current JSON does not contain usable ASE Atoms data: ${extracted.error}`
      );
      return;
    }
    atoms = extracted.atoms;
  } catch (error) {
    vscode.window.showErrorMessage(`Current document is not valid ASE Atoms JSON: ${error}`);
    return;
  }

  const formats = converter.getSupportedFormats();
  const items = Object.entries(formats).map(([key, value]) => ({
    label: value.name,
    description: value.description,
    detail: `Extensions: ${value.extensions.join(', ')}`,
    format: key as ASEFormat,
  }));

  const selected = await vscode.window.showQuickPick(items, {
    placeHolder: 'Select target format',
  });

  if (!selected) {
    return;
  }

  // Ask for output file path
  const outputPath = await vscode.window.showInputBox({
    prompt: 'Enter output file path',
    value: `output${formats[selected.format].extensions[0]}`,
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a file path';
      }
      return null;
    },
  });

  if (!outputPath) {
    return;
  }

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Writing ${selected.format} from ASE Atoms`,
      cancellable: false,
    },
    async () => {
      const result = await converter.writeFromAtoms(atoms, outputPath, selected.format);
      if (!result.success) {
        vscode.window.showErrorMessage(`Conversion failed: ${result.error}`);
        return;
      }

      const doc = await vscode.workspace.openTextDocument(outputPath);
      await vscode.window.showTextDocument(doc);
      vscode.window.showInformationMessage(`Successfully wrote ${outputPath}`);
    }
  );
}

/**
 * Migrate from one format to another
 */
async function migrateFormat(converter: ASEConverter): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const inputPath = editor.document.uri.fsPath;

  // Select target format
  const formats = converter.getSupportedFormats();
  const items = Object.entries(formats).map(([key, value]) => ({
    label: value.name,
    description: value.description,
    detail: `Extensions: ${value.extensions.join(', ')}`,
    format: key as ASEFormat,
  }));

  const selected = await vscode.window.showQuickPick(items, {
    placeHolder: 'Select target format for migration',
  });

  if (!selected) {
    return;
  }

  // Suggest output filename
  const inputBasename = inputPath.replace(/\.[^/.]+$/, '');
  const suggestedOutput = `${inputBasename}_${selected.format}${
    formats[selected.format].extensions[0]
  }`;

  const outputPath = await vscode.window.showInputBox({
    prompt: 'Enter output file path',
    value: suggestedOutput,
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a file path';
      }
      return null;
    },
  });

  if (!outputPath) {
    return;
  }

  // Perform conversion
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Migrating to ${selected.format}`,
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Converting format...' });

      const result = await converter.convertFormat(
        inputPath,
        outputPath,
        undefined,
        selected.format
      );

      if (!result.success) {
        vscode.window.showErrorMessage(`Migration failed: ${result.error}`);
        return;
      }

      progress.report({ message: 'Opening result...' });

      // Open the converted file
      const doc = await vscode.workspace.openTextDocument(outputPath);
      await vscode.window.showTextDocument(doc);

      vscode.window.showInformationMessage(
        `Successfully migrated to ${selected.format}: ${outputPath}`
      );
    }
  );
}

/**
 * Quick convert with common presets
 */
async function quickConvert(
  converter: ASEConverter,
  fromFormat?: string,
  toFormat?: string
): Promise<void> {
  const presets = [
    { label: 'VASP → CP2K', input: ASEFormat.VASP, output: ASEFormat.CP2K },
    { label: 'VASP → Quantum ESPRESSO', input: ASEFormat.VASP, output: ASEFormat.QE },
    { label: 'CP2K → VASP', input: ASEFormat.CP2K, output: ASEFormat.VASP },
    { label: 'Gaussian → ORCA', input: ASEFormat.Gaussian, output: ASEFormat.ORCA },
    { label: 'ORCA → Gaussian', input: ASEFormat.ORCA, output: ASEFormat.Gaussian },
    { label: 'XYZ → VASP', input: ASEFormat.XYZ, output: ASEFormat.VASP },
    { label: 'CIF → VASP', input: ASEFormat.CIF, output: ASEFormat.VASP },
  ];

  const selected =
    fromFormat && toFormat
      ? {
          label: `${fromFormat.toUpperCase()} → ${toFormat.toUpperCase()}`,
          input: normalizeASEFormat(fromFormat, converter.getSupportedFormats()),
          output: normalizeASEFormat(toFormat, converter.getSupportedFormats()),
        }
      : await vscode.window.showQuickPick(
          presets.map(p => p),
          { placeHolder: 'Select conversion preset' }
        );

  if (!selected) {
    return;
  }

  if (!selected.input || !selected.output) {
    const supported = Object.keys(converter.getSupportedFormats()).join(', ');
    vscode.window.showErrorMessage(
      `Unsupported ASE quick conversion: ${fromFormat ?? 'unknown'} → ${
        toFormat ?? 'unknown'
      }. Supported formats: ${supported}`
    );
    return;
  }

  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const inputPath = editor.document.uri.fsPath;
  const inputBasename = inputPath.replace(/\.[^/.]+$/, '');
  const outputPath = `${inputBasename}_${selected.output}${getExtension(selected.output)}`;

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Converting: ${selected.label}`,
      cancellable: false,
    },
    async progress => {
      const result = await converter.convertFormat(
        inputPath,
        outputPath,
        selected.input,
        selected.output
      );

      if (!result.success) {
        vscode.window.showErrorMessage(`Conversion failed: ${result.error}`);
        return;
      }

      const doc = await vscode.workspace.openTextDocument(outputPath);
      await vscode.window.showTextDocument(doc);

      vscode.window.showInformationMessage(`Successfully converted: ${outputPath}`);
    }
  );
}

/**
 * Format ASE Atoms information for display
 */
function formatAtomsInfo(result: ConversionResult): string {
  if (!result.atoms) {
    return JSON.stringify(result, null, 2);
  }

  const extracted = extractASEAtoms(result.atoms);
  const atoms = extracted.atoms ?? result.atoms;
  const formula = getChemicalFormula(atoms.chemical_symbols);

  const info = {
    summary: {
      formula: formula,
      natoms: atoms.chemical_symbols.length,
      periodic: atoms.pbc.some(p => p),
      elements: [...new Set(atoms.chemical_symbols)],
    },
    aseAtoms: atoms,
    structure: {
      chemical_symbols: atoms.chemical_symbols,
      positions: atoms.positions,
      cell: atoms.cell,
      pbc: atoms.pbc,
    },
    metadata: result.metadata,
    warnings: result.warnings,
  };

  return JSON.stringify(info, null, 2);
}

function extractASEAtoms(value: unknown): AtomsExtraction {
  const candidates = getAtomsCandidates(value);
  if (candidates.length === 0) {
    return {
      error:
        'expected a JSON object containing aseAtoms, atoms, data.atoms, structure, or a top-level Atoms object',
    };
  }

  const errors: string[] = [];
  for (const candidate of candidates) {
    const normalized = normalizeAtoms(candidate.value);
    if ('value' in normalized) {
      return { atoms: normalized.value, source: candidate.source };
    }
    errors.push(`${candidate.source}: ${normalized.error}`);
  }

  return {
    error: errors.slice(0, 3).join('; '),
  };
}

function getAtomsCandidates(value: unknown): Array<{ value: unknown; source: string }> {
  if (!isRecord(value)) {
    return [];
  }

  const candidates: Array<{ value: unknown; source: string }> = [];

  addCandidate(candidates, value.aseAtoms, 'aseAtoms');
  addCandidate(candidates, value.atoms, 'atoms');
  if (isRecord(value.data)) {
    addCandidate(candidates, value.data.aseAtoms, 'data.aseAtoms');
    addCandidate(candidates, value.data.atoms, 'data.atoms');
    addCandidate(candidates, value.data.structure, 'data.structure');
  }
  addCandidate(candidates, value.structure, 'structure');
  addCandidate(candidates, value, 'top-level');

  return candidates;
}

function addCandidate(
  candidates: Array<{ value: unknown; source: string }>,
  value: unknown,
  source: string
): void {
  if (value !== undefined && value !== null) {
    candidates.push({ value, source });
  }
}

function normalizeAtoms(value: unknown): ValidationResult<ASEAtoms> {
  if (!isRecord(value)) {
    return { error: 'candidate is not an object' };
  }

  const symbolValue = value.chemical_symbols ?? value.chemicalSymbols ?? value.symbols;
  const symbols = normalizeSymbols(symbolValue, value.numbers);
  if ('error' in symbols) {
    return { error: symbols.error };
  }

  const positions = normalizeVector3List(value.positions, 'positions');
  if ('error' in positions) {
    return { error: positions.error };
  }

  if (symbols.value.length === 0) {
    return { error: 'chemical_symbols must contain at least one atom' };
  }

  if (symbols.value.length !== positions.value.length) {
    return {
      error: `chemical_symbols has ${symbols.value.length} entries but positions has ${positions.value.length}`,
    };
  }

  const pbc = normalizePbc(value.pbc);
  if ('error' in pbc) {
    return { error: pbc.error };
  }

  const cell = normalizeCell(value.cell);
  if ('error' in cell) {
    return { error: cell.error };
  }

  const atoms: ASEAtoms = {
    chemical_symbols: symbols.value,
    positions: positions.value,
    cell: cell.value,
    pbc: pbc.value,
  };

  if (value.info !== undefined) {
    if (!isRecord(value.info)) {
      return { error: 'info must be an object when provided' };
    }
    atoms.info = value.info;
  }

  if (value.masses !== undefined) {
    const masses = normalizeNumberList(value.masses, 'masses');
    if ('error' in masses) {
      return { error: masses.error };
    }
    if (masses.value.length !== symbols.value.length) {
      return {
        error: `masses has ${masses.value.length} entries but chemical_symbols has ${symbols.value.length}`,
      };
    }
    atoms.masses = masses.value;
  }

  if (typeof value.constraints === 'string') {
    atoms.constraints = value.constraints;
  }

  return { value: atoms };
}

function normalizeSymbols(symbolValue: unknown, numbersValue: unknown): ValidationResult<string[]> {
  if (Array.isArray(symbolValue)) {
    const symbols: string[] = [];
    for (const symbol of symbolValue) {
      if (typeof symbol !== 'string' || !ELEMENT_SYMBOL_SET.has(symbol)) {
        return { error: 'chemical_symbols must be valid element symbols' };
      }
      symbols.push(symbol);
    }
    return { value: symbols };
  }

  if (Array.isArray(numbersValue)) {
    const symbols: string[] = [];
    for (const number of numbersValue) {
      if (!Number.isInteger(number) || number < 1 || number >= ELEMENT_SYMBOLS.length) {
        return { error: 'numbers must contain valid atomic numbers' };
      }
      symbols.push(ELEMENT_SYMBOLS[number]);
    }
    return { value: symbols };
  }

  return { error: 'missing chemical_symbols, symbols, or numbers array' };
}

function normalizeVector3List(value: unknown, label: string): ValidationResult<number[][]> {
  if (!Array.isArray(value)) {
    return { error: `${label} must be an array` };
  }

  const rows: number[][] = [];
  for (const row of value) {
    if (!Array.isArray(row) || row.length !== 3 || !row.every(isFiniteNumber)) {
      return { error: `${label} must contain [x, y, z] numeric rows` };
    }
    rows.push([...row]);
  }
  return { value: rows };
}

function normalizeCell(value: unknown): ValidationResult<number[][] | null> {
  if (value === undefined || value === null) {
    return { value: null };
  }

  if (!Array.isArray(value) || value.length !== 3) {
    return { error: 'cell must be a 3-vector or 3x3 numeric matrix' };
  }

  if (value.every(isFiniteNumber)) {
    return {
      value: [
        [value[0], 0, 0],
        [0, value[1], 0],
        [0, 0, value[2]],
      ],
    };
  }

  return normalizeVector3List(value, 'cell');
}

function normalizePbc(value: unknown): ValidationResult<boolean[]> {
  if (value === undefined || value === null) {
    return { value: [false, false, false] };
  }

  if (typeof value === 'boolean') {
    return { value: [value, value, value] };
  }

  if (Array.isArray(value)) {
    if (value.length === 1 && typeof value[0] === 'boolean') {
      return { value: [value[0], value[0], value[0]] };
    }
    if (value.length === 3 && value.every(item => typeof item === 'boolean')) {
      return { value: [...value] };
    }
  }

  return { error: 'pbc must be a boolean or a 3-element boolean array' };
}

function normalizeNumberList(value: unknown, label: string): ValidationResult<number[]> {
  if (!Array.isArray(value) || !value.every(isFiniteNumber)) {
    return { error: `${label} must be a numeric array` };
  }
  return { value: [...value] };
}

function isRecord(value: unknown): value is Record<string, any> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function isFiniteNumber(value: unknown): value is number {
  return typeof value === 'number' && Number.isFinite(value);
}

function normalizeASEFormat(
  value: string,
  formats: Record<string, { name: string }>
): ASEFormat | undefined {
  const normalized = value
    .trim()
    .toLowerCase()
    .replace(/[\s_-]+/g, '');
  const aliases: Record<string, ASEFormat> = {
    poscar: ASEFormat.VASP,
    contcar: ASEFormat.VASP,
    vasp: ASEFormat.VASP,
    espresso: ASEFormat.QE,
    quantumespresso: ASEFormat.QE,
    qe: ASEFormat.QE,
    gaussian: ASEFormat.Gaussian,
    gau: ASEFormat.Gaussian,
    orca: ASEFormat.ORCA,
    nwchem: ASEFormat.NWChem,
    gamess: ASEFormat.GAMESS,
    lammps: ASEFormat.LAMMPS,
    xyz: ASEFormat.XYZ,
    pdb: ASEFormat.PDB,
    cif: ASEFormat.CIF,
    cp2k: ASEFormat.CP2K,
  };

  if (aliases[normalized]) {
    return aliases[normalized];
  }

  const formatKey = Object.keys(formats).find(key => key.toLowerCase() === value.toLowerCase());
  return formatKey as ASEFormat | undefined;
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

/**
 * Get default extension for format
 */
function getExtension(format: ASEFormat): string {
  const extensions: Record<ASEFormat, string> = {
    [ASEFormat.VASP]: '.POSCAR',
    [ASEFormat.CP2K]: '.inp',
    [ASEFormat.QE]: '.in',
    [ASEFormat.Gaussian]: '.com',
    [ASEFormat.ORCA]: '.inp',
    [ASEFormat.NWChem]: '.nw',
    [ASEFormat.GAMESS]: '.inp',
    [ASEFormat.LAMMPS]: '.lmp',
    [ASEFormat.XYZ]: '.xyz',
    [ASEFormat.PDB]: '.pdb',
    [ASEFormat.CIF]: '.cif',
  };
  return extensions[format] || '.txt';
}
