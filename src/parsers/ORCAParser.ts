import {
  BaseParser,
  ParseResult,
  ParsedSection,
  ParsedParameter,
  ValidationResult,
  ParseError,
  ParseWarning,
} from './base';

export class ORCAParser extends BaseParser {
  private parsedResult: ParseResult | null = null;

  parseInput(): ParseResult {
    if (this.parsedResult) {
      return this.parsedResult;
    }

    const sections: ParsedSection[] = [];
    const parameters: ParsedParameter[] = [];
    const errors: ParseError[] = [];
    const warnings: ParseWarning[] = [];

    let currentSection: ParsedSection | null = null;
    let inBlock = false;
    let blockName = '';

    for (let i = 0; i < this.lines.length; i++) {
      const line = this.lines[i];
      const trimmed = line.trim();

      if (!trimmed || trimmed.startsWith('#')) {
        continue;
      }

      if (trimmed.startsWith('!')) {
        const simpleLine = trimmed.substring(1).trim();
        const simpleParams = this.parseSimpleInputLine(simpleLine, i);
        parameters.push(...simpleParams);
      } else if (trimmed.startsWith('%')) {
        if (inBlock && blockName) {
          if (currentSection) {
            currentSection.endLine = i - 1;
            sections.push(currentSection);
          }
        }

        blockName = trimmed.slice(1).split(/\s+/)[0];
        inBlock = true;
        currentSection = {
          name: blockName,
          startLine: i,
          endLine: i,
          parameters: [],
        };
      } else if (inBlock && trimmed.toLowerCase() === 'end') {
        if (currentSection) {
          currentSection.endLine = i;
          sections.push(currentSection);
        }
        inBlock = false;
        blockName = '';
        currentSection = null;
      } else if (inBlock && currentSection) {
        const kv = this.parseKeyValue(trimmed);
        if (kv) {
          currentSection.parameters.push({
            name: kv.key,
            value: this.convertValue(kv.value),
            line: i,
          });
        } else if (trimmed) {
          currentSection.parameters.push({
            name: 'Value',
            value: trimmed,
            line: i,
          });
        }
      } else if (trimmed && !trimmed.startsWith('*')) {
        const coordMatch = trimmed.match(/^([A-Za-z]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)/);
        if (coordMatch) {
          parameters.push({
            name: `Atom_${coordMatch[1]}`,
            value: `${parseFloat(coordMatch[2])} ${parseFloat(coordMatch[3])} ${parseFloat(
              coordMatch[4]
            )}`,
            line: i,
          });
        }
      }
    }

    if (inBlock && currentSection) {
      currentSection.endLine = this.lines.length - 1;
      sections.push(currentSection);
    }

    if (parameters.length === 0 && sections.length === 0) {
      errors.push({
        message: 'No valid ORCA input found',
        line: 0,
        severity: 'error',
      });
    }

    this.parsedResult = { sections, parameters, errors, warnings };
    return this.parsedResult;
  }

  private parseSimpleInputLine(line: string, lineNum: number): ParsedParameter[] {
    const params: ParsedParameter[] = [];
    const words = line.split(/\s+/).filter(w => w);

    const methodKeywords = [
      'HF',
      'DFT',
      'MP2',
      'CCSD',
      'CCSD(T)',
      'CASSCF',
      'RHF',
      'UHF',
      'ROHF',
      'B3LYP',
      'PBE',
      'BLYP',
      'BP86',
      'M06',
      'WB97X',
      'WB97X-D3',
      'TPSS',
    ];
    const basisKeywords = [
      'def2-SVP',
      'def2-TZVP',
      'def2-TZVPP',
      'def2-QZVP',
      '6-31G',
      '6-311G',
      'cc-pVDZ',
      'cc-pVTZ',
      'cc-pVQZ',
      'STO-3G',
    ];
    const calcTypeKeywords = ['OPT', 'FREQ', 'SP', 'TDDFT', 'COPT', 'NUMFREQ', 'GRAD', 'ENERGY'];

    for (const word of words) {
      const upper = word.toUpperCase();

      if (methodKeywords.some(m => upper === m || upper.includes(m))) {
        params.push({ name: 'Method', value: word, line: lineNum });
      } else if (
        basisKeywords.some(
          b => upper === b.toUpperCase() || upper.includes(b.toUpperCase().replace(/-/g, ''))
        )
      ) {
        params.push({ name: 'Basis', value: word, line: lineNum });
      } else if (calcTypeKeywords.some(c => upper === c || upper.includes(c))) {
        params.push({ name: 'CalculationType', value: word, line: lineNum });
      } else {
        params.push({ name: 'Option', value: word, line: lineNum });
      }
    }

    return params;
  }

  private convertValue(value: string): string | number | boolean {
    const lower = value.toLowerCase();
    if (lower === 'true') {
      return true;
    }
    if (lower === 'false') {
      return false;
    }

    const num = Number(value);
    if (!isNaN(num)) {
      return num;
    }

    return value;
  }

  validate(): ValidationResult {
    const result = this.parseInput();
    const errors = [...result.errors];
    const warnings = [...result.warnings];

    const hasMethod = result.parameters.some(p => p.name === 'Method');
    if (!hasMethod) {
      warnings.push({
        message: 'No calculation method specified in simple input',
        line: 0,
      });
    }

    for (const section of result.sections) {
      if (!section.parameters.length && section.name !== 'PAL') {
        warnings.push({
          message: `Section ${section.name} has no parameters`,
          line: section.startLine,
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
}
