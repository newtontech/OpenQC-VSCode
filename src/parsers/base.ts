/**
 * Base parser interface for quantum chemistry input files
 *
 * All quantum chemistry parsers (VASP, Gaussian, ORCA, CP2K, QE, GAMESS, NWChem)
 * extend this class to share common parsing utilities.
 *
 * @module parsers/base
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/16
 */

/**
 * Represents an atom extracted from coordinate section of an input file.
 * Parsers that support coordinate extraction should implement `CoordinateCapable`.
 */
export interface ExtractedAtom {
  element: string;
  coords: [number, number, number];
}

/**
 * Interface for parsers that can extract atomic coordinates.
 *
 * Each parser implements this according to its file format conventions.
 * This enables `StructureConverter` and `Molecule3D` to consume coordinates
 * from any parser through a uniform contract.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/16
 */
export interface CoordinateCapable {
  /**
   * Extract atomic coordinates from the parsed input file.
   *
   * @returns Array of atoms with element symbols and xyz coordinates
   */
  getExtractedAtoms(): ExtractedAtom[];
}

export interface ParsedParameter {
  name: string;
  value: string | number | boolean | number[];
  line: number;
  description?: string;
}

export interface ParsedSection {
  name: string;
  startLine: number;
  endLine: number;
  parameters: ParsedParameter[];
  subsections?: ParsedSection[];
}

export interface ParseResult {
  sections: ParsedSection[];
  parameters: ParsedParameter[];
  errors: ParseError[];
  warnings: ParseWarning[];
}

export interface ParseError {
  message: string;
  line: number;
  severity: 'error' | 'warning';
}

export interface ParseWarning {
  message: string;
  line: number;
}

export interface ValidationResult {
  valid: boolean;
  errors: ParseError[];
  warnings: ParseWarning[];
}

export abstract class BaseParser {
  protected content: string;
  protected lines: string[];

  constructor(content: string) {
    this.content = content;
    this.lines = content.split('\n');
  }

  abstract parseInput(): ParseResult;
  abstract validate(): ValidationResult;
  abstract getSections(): ParsedSection[];
  abstract getParameters(): ParsedParameter[];
  abstract getParameter(name: string): ParsedParameter | undefined;
  abstract getSection(name: string): ParsedSection | undefined;

  protected getLineNumber(position: number): number {
    let line = 0;
    let currentPos = 0;
    for (const lineContent of this.lines) {
      if (currentPos + lineContent.length >= position) {
        return line;
      }
      currentPos += lineContent.length + 1;
      line++;
    }
    return line;
  }

  protected stripComments(line: string, commentChars: string[] = ['#', '!']): string {
    let result = line;
    for (const char of commentChars) {
      const index = result.indexOf(char);
      if (index !== -1) {
        result = result.substring(0, index);
      }
    }
    return result.trim();
  }

  protected parseKeyValue(
    line: string,
    delimiter: string | RegExp = /\s*=\s*|\s+/
  ): { key: string; value: string } | null {
    const parts = line.split(delimiter).filter(p => p.trim());
    if (parts.length >= 2) {
      return { key: parts[0].trim(), value: parts.slice(1).join(' ').trim() };
    }
    return null;
  }

  /**
   * Parse a single coordinate line in "Element X Y Z" format.
   *
   * Shared utility for parsers that extract coordinates from XYZ-style
   * coordinate blocks (Gaussian, ORCA, CP2K, GAMESS, NWChem, etc.).
   *
   * @param line - Raw coordinate line text
   * @param startIndex - Column index where coordinates begin (default: 1 for "Element X Y Z")
   * @returns Extracted atom or null if the line is not a valid coordinate
   */
  protected parseCoordinateLine(line: string, startIndex = 1): ExtractedAtom | null {
    const parts = line.trim().split(/\s+/);
    if (parts.length < startIndex + 3) {
      return null;
    }

    const element = parts[0];
    const x = parseFloat(parts[startIndex]);
    const y = parseFloat(parts[startIndex + 1]);
    const z = parseFloat(parts[startIndex + 2]);

    if (isNaN(x) || isNaN(y) || isNaN(z)) {
      return null;
    }

    return { element, coords: [x, y, z] };
  }
}
