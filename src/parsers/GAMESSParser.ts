import {
  BaseParser,
  ParseResult,
  ParsedSection,
  ParsedParameter,
  ValidationResult,
  ParseError,
  ParseWarning,
} from './base';

export class GAMESSParser extends BaseParser {
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
    let inDataGroup = false;

    for (let i = 0; i < this.lines.length; i++) {
      const line = this.lines[i];
      const trimmed = line.trim();

      if (!trimmed) {
        continue;
      }

      if (trimmed.startsWith('$')) {
        // Extract group name: the first word after $ (up to first whitespace or $END)
        const afterDollar = trimmed.substring(1).trim();
        const firstWord = afterDollar.split(/\s+/)[0].toUpperCase();
        const groupName = firstWord;

        if (groupName === 'END') {
          if (currentSection) {
            currentSection.endLine = i;
            sections.push(currentSection);
            currentSection = null;
          }
          inDataGroup = false;
        } else if (!inDataGroup) {
          if (currentSection) {
            currentSection.endLine = i - 1;
            sections.push(currentSection);
          }

          const groupParams = this.parseDataGroupLine(trimmed, i);
          currentSection = {
            name: groupName,
            startLine: i,
            endLine: i,
            parameters: groupParams,
          };
          parameters.push(...groupParams);
          inDataGroup = true;

          // Check if $END appears on the same line after the group params
          if (trimmed.toUpperCase().includes('$END') && trimmed.toUpperCase().indexOf('$END') > 0) {
            // Strip $END and everything after it from the line content
            currentSection.endLine = i;
            sections.push(currentSection);
            currentSection = null;
            inDataGroup = false;
          }
        }
      } else if (inDataGroup && currentSection) {
        const groupParams = this.parseDataGroupLine(trimmed, i);
        currentSection.parameters.push(...groupParams);
        parameters.push(...groupParams);
        currentSection.endLine = i;
      }
    }

    if (currentSection) {
      currentSection.endLine = this.lines.length - 1;
      sections.push(currentSection);
    }

    this.parsedResult = { sections, parameters, errors, warnings };
    return this.parsedResult;
  }

  private parseDataGroupLine(line: string, lineNum: number): ParsedParameter[] {
    const params: ParsedParameter[] = [];

    if (line.startsWith('$')) {
      const content = line.substring(line.indexOf('$') + 1).trim();
      if (content.toUpperCase() === 'END') {
        return params;
      }

      const groupNameEnd = content.search(/\s/);
      if (groupNameEnd === -1) {
        return params;
      }

      line = content.substring(groupNameEnd).trim();
    }

    // Strip inline $END and everything after it
    const endIdx = line.toUpperCase().indexOf('$END');
    if (endIdx >= 0) {
      line = line.substring(0, endIdx).trim();
    }

    const assignments = line.split(/\s+/).filter(s => s.includes('='));

    for (const assignment of assignments) {
      const kv = this.parseKeyValue(assignment, '=');
      if (kv) {
        params.push({
          name: kv.key.trim().toUpperCase(),
          value: this.convertValue(kv.value.trim()),
          line: lineNum,
        });
      }
    }

    return params;
  }

  private convertValue(value: string): string | number | boolean {
    const lower = value.toLowerCase();
    if (lower === '.true.' || lower === 'true' || lower === 't') {
      return true;
    }
    if (lower === '.false.' || lower === 'false' || lower === 'f') {
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

    const hasContrl = result.sections.some(s => s.name === 'CONTRL');
    const hasSystem = result.sections.some(s => s.name === 'SYSTEM');
    const hasBasis = result.sections.some(s => s.name === 'BASIS');

    if (!hasContrl) {
      errors.push({
        message: 'Missing required $CONTRL data group',
        line: 0,
        severity: 'error',
      });
    }

    if (!hasSystem) {
      errors.push({
        message: 'Missing required $SYSTEM data group',
        line: 0,
        severity: 'error',
      });
    }

    if (!hasBasis) {
      errors.push({
        message: 'Missing required $BASIS data group',
        line: 0,
        severity: 'error',
      });
    }

    const contrlSection = result.sections.find(s => s.name === 'CONTRL');
    if (contrlSection) {
      const hasRuntyp = contrlSection.parameters.some(p => p.name === 'RUNTYP');
      if (!hasRuntyp) {
        errors.push({
          message: 'RUNTYP must be specified in $CONTRL',
          line: contrlSection.startLine,
          severity: 'error',
        });
      }

      const hasScftyp = contrlSection.parameters.some(p => p.name === 'SCFTYP');
      if (!hasScftyp) {
        errors.push({
          message: 'SCFTYP must be specified in $CONTRL',
          line: contrlSection.startLine,
          severity: 'error',
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
