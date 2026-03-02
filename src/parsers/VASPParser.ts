import {
  BaseParser,
  ParseResult,
  ParsedSection,
  ParsedParameter,
  ValidationResult,
  ParseError,
  ParseWarning,
} from './base';

export interface VASPParameters {
  system?: string;
  start?: string;
  charge?: number;
  spin?: number;
  encut?: number;
  prec?: string;
  ibrion?: number;
  nsw?: number;
  isif?: number;
  ismear?: number;
  sigma?: number;
}

export interface POSCARData {
  comment: string;
  scale: number;
  lattice: number[][];
  atomTypes: string[];
  atomCounts: number[];
  coordinateType: 'Direct' | 'Cartesian';
  coordinates: number[][];
}

export class VASPParser extends BaseParser {
  private parsedResult: ParseResult | null = null;
  private filename: string;

  constructor(content: string, filename: string = 'INCAR') {
    super(content);
    this.filename = filename.toUpperCase();
  }

  parseInput(): ParseResult {
    if (this.parsedResult) {
      return this.parsedResult;
    }

    switch (this.filename) {
      case 'POSCAR':
        return this.parsePOSCAR();
      case 'KPOINTS':
        return this.parseKPOINTS();
      case 'POTCAR':
        return this.parsePOTCAR();
      case 'INCAR':
      default:
        return this.parseINCAR();
    }
  }

  private parseINCAR(): ParseResult {
    const sections: ParsedSection[] = [];
    const parameters: ParsedParameter[] = [];
    const errors: ParseError[] = [];
    const warnings: ParseWarning[] = [];

    for (let i = 0; i < this.lines.length; i++) {
      const line = this.lines[i].trim();

      // Skip empty lines and comments
      if (!line || line.startsWith('#')) continue;

      // Remove inline comments
      const lineWithoutComment = line.split('#')[0].trim();
      if (!lineWithoutComment) continue;

      const kv = this.parseKeyValue(lineWithoutComment, '=');
      if (kv) {
        parameters.push({
          name: kv.key.toUpperCase(),
          value: this.convertValue(kv.value),
          line: i,
        });
      } else if (lineWithoutComment && !lineWithoutComment.includes('=')) {
        // Line without '=' is malformed
        warnings.push({
          message: `Possible malformed line: "${line}"`,
          line: i,
        });
      }
    }

    sections.push({
      name: 'INCAR',
      startLine: 0,
      endLine: this.lines.length - 1,
      parameters,
    });

    this.parsedResult = { sections, parameters, errors, warnings };
    return this.parsedResult;
  }

  private parsePOSCAR(): ParseResult {
    const sections: ParsedSection[] = [];
    const parameters: ParsedParameter[] = [];
    const errors: ParseError[] = [];
    const warnings: ParseWarning[] = [];

    if (this.lines.length < 8) {
      errors.push({
        message: 'POSCAR file too short',
        line: 0,
        severity: 'error',
      });
      return { sections, parameters, errors, warnings };
    }

    let lineIdx = 0;

    const comment = this.lines[lineIdx++].trim();
    parameters.push({ name: 'Comment', value: comment, line: 0 });

    const scale = parseFloat(this.lines[lineIdx++].trim());
    if (isNaN(scale)) {
      errors.push({ message: 'Invalid scaling factor', line: 1, severity: 'error' });
    }
    parameters.push({ name: 'SCALE', value: scale, line: 1 });

    const lattice: number[][] = [];
    for (let i = 0; i < 3; i++) {
      const parts = this.lines[lineIdx++].trim().split(/\s+/).map(Number);
      if (parts.length >= 3) {
        lattice.push(parts.slice(0, 3));
      }
    }

    const atomTypes = this.lines[lineIdx++].trim().split(/\s+/);
    const atomCounts = this.lines[lineIdx++].trim().split(/\s+/).map(Number);

    atomTypes.forEach((type, idx) => {
      parameters.push({
        name: `Atom_${type}`,
        value: atomCounts[idx] || 0,
        line: lineIdx - 2 + idx,
      });
    });

    const coordType = this.lines[lineIdx++].trim();

    sections.push({
      name: 'POSCAR',
      startLine: 0,
      endLine: this.lines.length - 1,
      parameters,
    });

    this.parsedResult = { sections, parameters, errors, warnings };
    return this.parsedResult;
  }

  private parseKPOINTS(): ParseResult {
    const sections: ParsedSection[] = [];
    const parameters: ParsedParameter[] = [];
    const errors: ParseError[] = [];
    const warnings: ParseWarning[] = [];

    if (this.lines.length < 4) {
      errors.push({
        message: 'KPOINTS file too short',
        line: 0,
        severity: 'error',
      });
      return { sections, parameters, errors, warnings };
    }

    const comment = this.lines[0].trim();
    const kptStyle = parseInt(this.lines[1].trim());
    const kptScheme = this.lines[2].trim();

    parameters.push({ name: 'Comment', value: comment, line: 0 });
    parameters.push({ name: 'KpointStyle', value: kptStyle, line: 1 });
    parameters.push({ name: 'KpointScheme', value: kptScheme, line: 2 });

    const kptLine = this.lines[3].trim().split(/\s+/).map(Number);
    if (kptLine.length >= 3) {
      parameters.push({ name: 'KX', value: kptLine[0], line: 3 });
      parameters.push({ name: 'KY', value: kptLine[1], line: 3 });
      parameters.push({ name: 'KZ', value: kptLine[2], line: 3 });
    }

    sections.push({
      name: 'KPOINTS',
      startLine: 0,
      endLine: this.lines.length - 1,
      parameters,
    });

    this.parsedResult = { sections, parameters, errors, warnings };
    return this.parsedResult;
  }

  private parsePOTCAR(): ParseResult {
    const sections: ParsedSection[] = [];
    const parameters: ParsedParameter[] = [];
    const errors: ParseError[] = [];
    const warnings: ParseWarning[] = [];

    if (this.lines.length < 1 || !this.lines[0].trim()) {
      errors.push({
        message: 'POTCAR file is empty',
        line: 0,
        severity: 'error',
      });
      return { sections, parameters, errors, warnings };
    }

    // Parse header line: "PAW_PBE H 01Jan2001"
    this.parsePOTCARHeader(parameters);

    // Parse energy cutoffs and other parameters
    this.parsePOTCAREnergyParameters(parameters);

    sections.push({
      name: 'POTCAR',
      startLine: 0,
      endLine: this.lines.length - 1,
      parameters,
    });

    this.parsedResult = { sections, parameters, errors, warnings };
    return this.parsedResult;
  }

  /**
   * Parse POTCAR header line containing type, element, and date
   * Format: "PAW_PBE H 01Jan2001"
   */
  private parsePOTCARHeader(parameters: ParsedParameter[]): void {
    const headerLine = this.lines[0].trim();
    const headerParts = headerLine.split(/\s+/);

    if (headerParts.length >= 3) {
      const [potcarType, element, date] = headerParts;
      parameters.push({ name: 'POTCARType', value: potcarType, line: 0 });
      parameters.push({ name: 'Element', value: element, line: 0 });
      parameters.push({ name: 'Date', value: date, line: 0 });
    }
  }

  /**
   * Parse energy cutoff parameters from POTCAR
   * Extracts ENMAX and ENMIN values
   */
  private parsePOTCAREnergyParameters(parameters: ParsedParameter[]): void {
    const energyParamRegex = /(ENMAX|ENMIN)\s*=\s*([\d.]+)/;

    for (let i = 1; i < this.lines.length; i++) {
      const match = this.lines[i].match(energyParamRegex);
      if (match) {
        const [, paramName, paramValue] = match;
        parameters.push({
          name: paramName,
          value: parseFloat(paramValue),
          line: i,
        });
      }
    }
  }

  private convertValue(value: string): string | number | boolean {
    const lower = value.toLowerCase();
    if (lower === 'true' || lower === '.true.') {
      return true;
    }
    if (lower === 'false' || lower === '.false.') {
      return false;
    }

    const num = Number(value);
    if (!isNaN(num) && value !== '') {
      return num;
    }

    return value;
  }

  validate(): ValidationResult {
    const result = this.parseInput();
    const errors = [...result.errors];
    const warnings = [...result.warnings];

    if (this.filename === 'INCAR') {
      const required = ['ENCUT', 'PREC'];
      const present = result.parameters.map(p => p.name.toUpperCase());

      for (const req of required) {
        if (!present.includes(req)) {
          warnings.push({
            message: `Missing recommended parameter: ${req}`,
            line: 0,
          });
        }
      }
    }

    return {
      valid: errors.length === 0,
      errors,
      warnings,
    };
  }

  getSections(): ParsedSection[] {
    return this.parseInput().sections;
  }

  getParameters(): ParsedParameter[] {
    return this.parseInput().parameters;
  }

  getParameter(name: string): ParsedParameter | undefined {
    return this.parseInput().parameters.find(p => p.name.toUpperCase() === name.toUpperCase());
  }

  getSection(name: string): ParsedSection | undefined {
    return this.parseInput().sections.find(s => s.name.toUpperCase() === name.toUpperCase());
  }

  /**
   * Extract atomic coordinates from POSCAR file
   * Returns array of [x, y, z] coordinates
   */
  getCoordinates(): number[][] {
    if (this.filename !== 'POSCAR') {
      return [];
    }

    const lines = this.lines.filter(line => line.trim() !== '');
    if (lines.length < 8) {
      return [];
    }

    let lineIdx = 0;
    lineIdx++; // Skip comment
    lineIdx++; // Skip scale

    // Skip lattice vectors (3 lines)
    lineIdx += 3;

    // Skip atom types line
    lineIdx++;

    // Get atom counts
    const atomCounts = lines[lineIdx].trim().split(/\s+/).map(Number);
    const totalAtoms = atomCounts.reduce((sum, count) => sum + count, 0);
    lineIdx++;

    // Check for selective dynamics
    if (lines[lineIdx].trim().toLowerCase().startsWith('selective')) {
      lineIdx++;
    }

    // Skip coordinate type line (Direct/Cartesian)
    lineIdx++;

    // Extract coordinates
    const coordinates: number[][] = [];
    for (let i = 0; i < totalAtoms && lineIdx < lines.length; i++) {
      const parts = lines[lineIdx].trim().split(/\s+/);
      if (parts.length >= 3) {
        const x = parseFloat(parts[0]);
        const y = parseFloat(parts[1]);
        const z = parseFloat(parts[2]);
        if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
          coordinates.push([x, y, z]);
        }
      }
      lineIdx++;
    }

    return coordinates;
  }

  /**
   * Extract lattice vectors from POSCAR file
   * Returns array of 3 vectors, each with [x, y, z]
   */
  getLatticeVectors(): number[][] {
    if (this.filename !== 'POSCAR') {
      return [];
    }

    const lines = this.lines.filter(line => line.trim() !== '');
    if (lines.length < 5) {
      return [];
    }

    const lattice: number[][] = [];
    for (let i = 2; i < 5 && i < lines.length; i++) {
      const parts = lines[i].trim().split(/\s+/).map(Number);
      if (parts.length >= 3) {
        lattice.push(parts.slice(0, 3));
      }
    }

    return lattice;
  }

  /**
   * Extract atom types from POSCAR file
   */
  getAtomTypes(): string[] {
    if (this.filename !== 'POSCAR') {
      return [];
    }

    const lines = this.lines.filter(line => line.trim() !== '');
    if (lines.length < 6) {
      return [];
    }

    // Atom types are on line 5 (0-indexed: after comment, scale, 3 lattice lines)
    const atomTypesLine = lines[5].trim();
    if (!atomTypesLine) {
      return [];
    }

    return atomTypesLine.split(/\s+/);
  }

  /**
   * Extract atom counts from POSCAR file
   */
  getAtomCounts(): number[] {
    if (this.filename !== 'POSCAR') {
      return [];
    }

    const lines = this.lines.filter(line => line.trim() !== '');
    if (lines.length < 7) {
      return [];
    }

    // Atom counts are on line 6
    const atomCountsLine = lines[6].trim();
    if (!atomCountsLine) {
      return [];
    }

    return atomCountsLine
      .split(/\s+/)
      .map(Number)
      .filter(n => !isNaN(n));
  }
}
