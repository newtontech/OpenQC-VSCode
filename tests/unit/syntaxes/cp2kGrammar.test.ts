import * as fs from 'fs';
import * as path from 'path';

const grammarPath = path.join(__dirname, '../../../syntaxes/cp2k.tmLanguage.json');

interface GrammarPattern {
  name?: string;
  match?: string;
  patterns?: GrammarPattern[];
  captures?: Record<string, { name?: string }>;
}

function collectPatterns(node: unknown, patterns: Array<{ name: string; regex: RegExp }> = []) {
  if (!node || typeof node !== 'object') return patterns;

  const value = node as GrammarPattern & Record<string, unknown>;
  if (typeof value.name === 'string' && typeof value.match === 'string') {
    patterns.push({ name: value.name, regex: new RegExp(value.match) });
  }

  if (value.captures && typeof value.match === 'string') {
    for (const capture of Object.values(value.captures)) {
      if (typeof capture.name === 'string') {
        patterns.push({ name: capture.name, regex: new RegExp(value.match) });
      }
    }
  }

  for (const child of Object.values(value)) {
    if (Array.isArray(child)) {
      child.forEach(item => collectPatterns(item, patterns));
    } else if (child && typeof child === 'object') {
      collectPatterns(child, patterns);
    }
  }

  return patterns;
}

function scopesFor(line: string): string[] {
  const grammar = JSON.parse(fs.readFileSync(grammarPath, 'utf8'));
  return collectPatterns(grammar)
    .filter(pattern => pattern.regex.test(line))
    .map(pattern => pattern.name);
}

describe('CP2K TextMate grammar', () => {
  it('declares the CP2K source scope', () => {
    const grammar = JSON.parse(fs.readFileSync(grammarPath, 'utf8'));
    expect(grammar.scopeName).toBe('source.cp2k');
    expect(grammar.repository).toHaveProperty('preprocessor');
    expect(grammar.repository).toHaveProperty('sections');
    expect(grammar.repository).toHaveProperty('values');
  });

  it('covers full CP2K input syntax used before LSP startup', () => {
    const examples: Record<string, string[]> = {
      '! comment': ['comment.line.cp2k'],
      '&GLOBAL': ['keyword.control.section.begin.cp2k'],
      '&END GLOBAL': ['keyword.control.section.end.cp2k'],
      '  RUN_TYPE ENERGY': ['variable.parameter.cp2k'],
      '  UKS TRUE': ['variable.parameter.cp2k', 'constant.language.boolean.cp2k'],
      '  CUTOFF 400': ['variable.parameter.cp2k', 'constant.numeric.integer.cp2k'],
      '  EPS_SCF 1.0E-7': ['variable.parameter.cp2k', 'constant.numeric.float.cp2k'],
      '  A [angstrom] 10.0 0.0 0.0': ['variable.parameter.cp2k', 'constant.other.unit.cp2k'],
      "@INCLUDE './common.inc'": [
        'keyword.control.preprocessor.cp2k',
        'string.unquoted.include.cp2k',
        'string.quoted.single.cp2k',
      ],
      '@SET PROJECT_NAME water': ['keyword.control.preprocessor.cp2k'],
      '  PROJECT ${PROJECT_NAME}': ['variable.parameter.cp2k', 'variable.other.cp2k'],
      '  COORD_FILE_NAME ./coords.xyz': ['variable.parameter.cp2k', 'string.unquoted.path.cp2k'],
      '  BASIS_SET DZVP-MOLOPT-SR-GTH': ['variable.parameter.cp2k'],
      '  POTENTIAL GTH-PBE': ['variable.parameter.cp2k', 'string.unquoted.datafile.cp2k'],
    };

    for (const [line, expectedScopes] of Object.entries(examples)) {
      const scopes = scopesFor(line);
      for (const expectedScope of expectedScopes) {
        expect(scopes).toContain(expectedScope);
      }
    }
  });
});
