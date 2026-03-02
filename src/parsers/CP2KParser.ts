import {
  BaseParser,
  ParseResult,
  ParsedSection,
  ParsedParameter,
  ValidationResult,
  ParseError,
  ParseWarning,
} from './base';

interface CP2KSectionInfo {
  name: string;
  elementName?: string; // For KIND sections, stores the element name
  startLine: number;
  endLine: number;
  parameters: ParsedParameter[];
  subsections?: ParsedSection[];
}

export class CP2KParser extends BaseParser {
  private parsedResult: ParseResult | null = null;
  private sectionInfos: Map<string, CP2KSectionInfo> = new Map();

  parseInput(): ParseResult {
    if (this.parsedResult) {
      return this.parsedResult;
    }

    const sections: ParsedSection[] = [];
    const parameters: ParsedParameter[] = [];
    const errors: ParseError[] = [];
    const warnings: ParseWarning[] = [];
    this.sectionInfos.clear();

    const sectionStack: CP2KSectionInfo[] = [];

    for (let i = 0; i < this.lines.length; i++) {
      const line = this.lines[i];
      const trimmed = line.trim();

      if (!trimmed || trimmed.startsWith('#')) {
        continue;
      }

      if (trimmed.startsWith('&')) {
        const parts = trimmed.substring(1).split(/\s+/);
        const sectionName = parts[0].toUpperCase();

        if (
          trimmed.toUpperCase().startsWith('&END') ||
          trimmed.toUpperCase() === '&' + sectionName + ' END'
        ) {
          if (sectionStack.length > 0) {
            const closedSection = sectionStack.pop()!;
            closedSection.endLine = i;

            if (sectionStack.length > 0) {
              if (!sectionStack[sectionStack.length - 1].subsections) {
                sectionStack[sectionStack.length - 1].subsections = [];
              }
              sectionStack[sectionStack.length - 1].subsections!.push(closedSection);
            } else {
              sections.push(closedSection);
            }
          }
        } else {
          // Extract element name for KIND and XC_FUNCTIONAL sections
          let elementName: string | undefined;
          if ((sectionName === 'KIND' || sectionName === 'XC_FUNCTIONAL') && parts.length > 1) {
            elementName = parts[1];
          }

          const newSection: CP2KSectionInfo = {
            name: sectionName,
            elementName,
            startLine: i,
            endLine: i,
            parameters: [],
          };
          sectionStack.push(newSection);

          // Store section info for later reference
          const key = `${sectionStack.length}-${sectionName}`;
          this.sectionInfos.set(key, newSection);
        }
      } else if (sectionStack.length > 0) {
        const kv = this.parseKeyValueCP2K(trimmed);
        if (kv) {
          const param: ParsedParameter = {
            name: kv.key,
            value: this.convertValue(kv.value),
            line: i,
          };
          sectionStack[sectionStack.length - 1].parameters.push(param);
          parameters.push(param);
        }
      } else {
        const kv = this.parseKeyValueCP2K(trimmed);
        if (kv) {
          parameters.push({
            name: kv.key,
            value: this.convertValue(kv.value),
            line: i,
          });
        }
      }
    }

    while (sectionStack.length > 0) {
      const section = sectionStack.pop()!;
      section.endLine = this.lines.length - 1;

      if (sectionStack.length > 0) {
        if (!sectionStack[sectionStack.length - 1].subsections) {
          sectionStack[sectionStack.length - 1].subsections = [];
        }
        sectionStack[sectionStack.length - 1].subsections!.push(section);
      } else {
        sections.push(section);
      }

      errors.push({
        message: `Section ${section.name} not properly closed`,
        line: section.startLine,
        severity: 'error',
      });
    }

    this.parsedResult = { sections, parameters, errors, warnings };
    return this.parsedResult;
  }

  private parseKeyValueCP2K(line: string): { key: string; value: string } | null {
    const parts = line.split(/\s+/).filter(p => p.trim());
    if (parts.length >= 2) {
      return { key: parts[0], value: parts.slice(1).join(' ') };
    }
    return null;
  }

  private convertValue(value: string): string | number | boolean {
    const lower = value.toLowerCase();
    if (lower === 'true' || lower === 'yes' || lower === 'on') {
      return true;
    }
    if (lower === 'false' || lower === 'no' || lower === 'off') {
      return false;
    }

    const num = Number(value);
    if (!isNaN(num)) {
      return num;
    }

    if (value.startsWith('"') && value.endsWith('"')) {
      return value.slice(1, -1);
    }

    return value;
  }

  validate(): ValidationResult {
    const result = this.parseInput();
    const errors = [...result.errors];
    const warnings = [...result.warnings];

    const hasGlobal = result.sections.some(s => s.name === 'GLOBAL');
    const hasForceEval = result.sections.some(s => s.name === 'FORCE_EVAL');

    if (!hasGlobal) {
      warnings.push({
        message: 'Missing &GLOBAL section',
        line: 0,
      });
    }

    if (!hasForceEval) {
      errors.push({
        message: 'Missing required &FORCE_EVAL section',
        line: 0,
        severity: 'error',
      });
    }

    const globalSection = result.sections.find(s => s.name === 'GLOBAL');
    if (globalSection) {
      const hasProjectName = globalSection.parameters.some(
        p => p.name.toUpperCase() === 'PROJECT_NAME'
      );
      if (!hasProjectName) {
        warnings.push({
          message: 'PROJECT_NAME not set in &GLOBAL',
          line: globalSection.startLine,
        });
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
    return this.parseInput().parameters.find(p => p.name.toLowerCase() === name.toLowerCase());
  }

  getSection(name: string): ParsedSection | undefined {
    return this.parseInput().sections.find(s => s.name.toLowerCase() === name.toLowerCase());
  }

  /**
   * Find a subsection recursively within a section
   */
  private findSubsection(
    section: ParsedSection,
    name: string
  ): ParsedSection | undefined {
    if (section.name.toLowerCase() === name.toLowerCase()) {
      return section;
    }
    if (section.subsections) {
      for (const subsection of section.subsections) {
        const found = this.findSubsection(subsection, name);
        if (found) {return found;}
      }
    }
    return undefined;
  }

  /**
   * Extract atomic coordinates from CP2K input file
   * Returns array of {element: string, coords: number[]} objects
   */
  getCoordinates(): Array<{ element: string; coords: number[] }> {
    const result = this.parseInput();
    const coordinates: Array<{ element: string; coords: number[] }> = [];

    // Find the SUBSYS section
    const subsysSection = result.sections.find(s => s.name === 'FORCE_EVAL');
    if (!subsysSection || !subsysSection.subsections) {
      return coordinates;
    }

    const subsys = this.findSubsection(subsysSection, 'SUBSYS');
    if (!subsys) {
      return coordinates;
    }

    // Find COORD section within SUBSYS
    const coordSection = this.findSubsection(subsys, 'COORD');
    if (!coordSection) {
      return coordinates;
    }

    // Parse coordinate lines
    for (let i = coordSection.startLine + 1; i < coordSection.endLine; i++) {
      const line = this.lines[i].trim();
      if (!line || line.startsWith('&') || line.startsWith('#')) {
        continue;
      }

      // CP2K coordinate format: ELEMENT x y z
      const parts = line.split(/\s+/).filter(p => p.trim());
      if (parts.length >= 4) {
        const element = parts[0];
        const x = parseFloat(parts[1]);
        const y = parseFloat(parts[2]);
        const z = parseFloat(parts[3]);

        if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
          coordinates.push({ element, coords: [x, y, z] });
        }
      }
    }

    return coordinates;
  }

  /**
   * Extract cell parameters from CP2K input file
   * Returns cell definition object
   */
  getCellParameters():
    | { type: 'ABC'; a: number; b: number; c: number }
    | { type: 'vectors'; a: number[]; b: number[]; c: number[] }
    | undefined {
    const result = this.parseInput();

    // Find the SUBSYS section
    const subsysSection = result.sections.find(s => s.name === 'FORCE_EVAL');
    if (!subsysSection || !subsysSection.subsections) {
      return undefined;
    }

    const subsys = this.findSubsection(subsysSection, 'SUBSYS');
    if (!subsys) {
      return undefined;
    }

    // Find CELL section within SUBSYS
    const cellSection = this.findSubsection(subsys, 'CELL');
    if (!cellSection) {
      return undefined;
    }

    // Check for ABC format
    const abcParam = cellSection.parameters.find(p => p.name.toUpperCase() === 'ABC');
    if (abcParam && typeof abcParam.value === 'string') {
      const parts = abcParam.value.split(/\s+/).map(Number);
      if (parts.length === 3 && parts.every(n => !isNaN(n))) {
        return { type: 'ABC', a: parts[0], b: parts[1], c: parts[2] };
      }
    }

    // Check for vector format
    const aParam = cellSection.parameters.find(p => p.name.toUpperCase() === 'A');
    const bParam = cellSection.parameters.find(p => p.name.toUpperCase() === 'B');
    const cParam = cellSection.parameters.find(p => p.name.toUpperCase() === 'C');

    if (aParam && bParam && cParam) {
      const parseVector = (val: any): number[] | undefined => {
        if (typeof val === 'string') {
          const parts = val.split(/\s+/).map(Number);
          if (parts.length >= 3 && parts.slice(0, 3).every(n => !isNaN(n))) {
            return parts.slice(0, 3);
          }
        }
        return undefined;
      };

      const a = parseVector(aParam.value);
      const b = parseVector(bParam.value);
      const c = parseVector(cParam.value);

      if (a && b && c) {
        return { type: 'vectors', a, b, c };
      }
    }

    return undefined;
  }

  /**
   * Extract atom types from KIND sections
   * Returns array of element names
   */
  getAtomTypes(): string[] {
    const result = this.parseInput();
    const atomTypes: string[] = [];

    // Find the FORCE_EVAL section
    const forceEvalSection = result.sections.find(s => s.name === 'FORCE_EVAL');
    if (!forceEvalSection || !forceEvalSection.subsections) {
      return atomTypes;
    }

    // Recursively find all KIND sections
    const findKindSections = (section: ParsedSection | any): ParsedSection[] => {
      const kinds: ParsedSection[] = [];
      if (section.name === 'KIND') {
        kinds.push(section);
      }
      if (section.subsections) {
        for (const subsection of section.subsections) {
          kinds.push(...findKindSections(subsection));
        }
      }
      return kinds;
    };

    const kindSections = findKindSections(forceEvalSection);

    for (const kind of kindSections) {
      // Check if this section has elementName (our custom property)
      const kindInfo = kind as any;
      if (kindInfo.elementName) {
        atomTypes.push(kindInfo.elementName);
      } else {
        // Fallback: try to find element from parameters
        const elementParam = kind.parameters?.find((p: any) =>
          p.name.toUpperCase() === '_SECTION_PARAMETERS'
        );
        if (elementParam && typeof elementParam.value === 'string') {
          atomTypes.push(elementParam.value.split(/\s+/)[0]);
        } else {
          // Last resort: use section name
          atomTypes.push(kind.name);
        }
      }
    }

    return atomTypes;
  }

  /**
   * Extract basis set information from KIND sections
   * Returns map of element -> basis set name
   */
  getBasisSets(): Map<string, string> {
    const result = this.parseInput();
    const basisSets = new Map<string, string>();

    const forceEvalSection = result.sections.find(s => s.name === 'FORCE_EVAL');
    if (!forceEvalSection || !forceEvalSection.subsections) {
      return basisSets;
    }

    const findKindSections = (section: ParsedSection | any): ParsedSection[] => {
      const kinds: ParsedSection[] = [];
      if (section.name === 'KIND') {
        kinds.push(section);
      }
      if (section.subsections) {
        for (const subsection of section.subsections) {
          kinds.push(...findKindSections(subsection));
        }
      }
      return kinds;
    };

    const kindSections = findKindSections(forceEvalSection);

    for (const kind of kindSections) {
      const kindInfo = kind as any;
      let element: string;
      if (kindInfo.elementName) {
        element = kindInfo.elementName;
      } else {
        const elementParam = kind.parameters?.find((p: any) =>
          p.name.toUpperCase() === '_SECTION_PARAMETERS'
        );
        if (elementParam && typeof elementParam.value === 'string') {
          element = elementParam.value.split(/\s+/)[0];
        } else {
          element = kind.name;
        }
      }

      const basisSetParam = kind.parameters?.find(
        (p: any) => p.name.toUpperCase() === 'BASIS_SET'
      );
      if (basisSetParam && typeof basisSetParam.value === 'string') {
        basisSets.set(element, basisSetParam.value);
      }
    }

    return basisSets;
  }

  /**
   * Extract DFT functional information
   * Returns the XC functional name if specified
   */
  getFunctional(): string | undefined {
    const result = this.parseInput();

    const forceEvalSection = result.sections.find(s => s.name === 'FORCE_EVAL');
    if (!forceEvalSection || !forceEvalSection.subsections) {
      return undefined;
    }

    const dftSection = this.findSubsection(forceEvalSection, 'DFT');
    if (!dftSection || !dftSection.subsections) {
      return undefined;
    }

    const xcSection = this.findSubsection(dftSection, 'XC');
    if (!xcSection || !xcSection.subsections) {
      return undefined;
    }

    const xcFunctionalSection = this.findSubsection(xcSection, 'XC_FUNCTIONAL');
    if (!xcFunctionalSection) {
      return undefined;
    }

    // Check for elementName (for sections like &XC_FUNCTIONAL PBE)
    const xcFunctionalInfo = xcFunctionalSection as any;
    if (xcFunctionalInfo.elementName) {
      return xcFunctionalInfo.elementName;
    }

    // The functional name might be in parameters
    const functionalParam = xcFunctionalSection.parameters?.find(
      (p: any) => p.name.toUpperCase() === '_SECTION_PARAMETERS'
    );

    if (functionalParam && typeof functionalParam.value === 'string') {
      // Extract just the functional name from the value
      const parts = functionalParam.value.split(/\s+/);
      return parts[0];
    }

    // Check for subsections like PBE, BLYP, etc.
    if (xcFunctionalSection.subsections && xcFunctionalSection.subsections.length > 0) {
      return xcFunctionalSection.subsections[0].name;
    }

    return undefined;
  }
}
