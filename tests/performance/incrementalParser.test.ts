import { IncrementalParser } from '../../src/performance/incrementalParser';

describe('IncrementalParser', () => {
  function createDocument(content: string, version: number, fsPath = '/tmp/input.inp'): any {
    return {
      uri: {
        fsPath,
        toString: () => `file://${fsPath}`,
      },
      version,
      lineCount: content.split('\n').length,
      getText: () => content,
      languageId: 'openqc-test',
    };
  }

  it('returns cached results for the same document version without reparsing', async () => {
    const parse = jest.fn((content: string) => ({ lines: content.split('\n') }));
    const parser = new IncrementalParser(parse, undefined, { minLinesForIncremental: 2 });
    const document = createDocument(['A', 'B', 'C'].join('\n'), 1);

    const first = await parser.parse(document);
    const second = await parser.parse(document);

    expect(first.cached).toBe(false);
    expect(second.cached).toBe(true);
    expect(parse).toHaveBeenCalledTimes(1);
  });

  it('detects added lines between cached document versions', async () => {
    const parser = new IncrementalParser(
      (content: string) => ({ lines: content.split('\n') }),
      undefined,
      { minLinesForIncremental: 2 }
    );

    await parser.parse(createDocument(['A', 'B', 'C', 'D', 'E', 'F'].join('\n'), 1));
    const result = await parser.parse(
      createDocument(['A', 'B', 'B2', 'C', 'D', 'E', 'F'].join('\n'), 2)
    );

    expect(result.changes.added).toEqual([
      { startLine: 2, endLine: 3, startColumn: 0, endColumn: 0 },
    ]);
    expect(result.changes.modified).toEqual([]);
    expect(result.ast.lines).toEqual(['A', 'B', 'B2', 'C', 'D', 'E', 'F']);
  });

  it('detects removed lines between cached document versions', async () => {
    const parser = new IncrementalParser(
      (content: string) => ({ lines: content.split('\n') }),
      undefined,
      { minLinesForIncremental: 2 }
    );

    await parser.parse(createDocument(['A', 'B', 'C', 'D', 'E', 'F'].join('\n'), 1));
    const result = await parser.parse(createDocument(['A', 'C', 'D', 'E', 'F'].join('\n'), 2));

    expect(result.changes.removed).toEqual([
      { startLine: 1, endLine: 2, startColumn: 0, endColumn: 0 },
    ]);
    expect(result.changes.added).toEqual([]);
  });

  it('detects modified lines between cached document versions', async () => {
    const parser = new IncrementalParser(
      (content: string) => ({ lines: content.split('\n') }),
      undefined,
      { minLinesForIncremental: 2 }
    );

    await parser.parse(createDocument(['A', 'B', 'C', 'D', 'E', 'F'].join('\n'), 1));
    const result = await parser.parse(
      createDocument(['A', 'B changed', 'C', 'D', 'E', 'F'].join('\n'), 2)
    );

    expect(result.changes.modified).toEqual([
      { startLine: 1, endLine: 2, startColumn: 0, endColumn: 9 },
    ]);
  });

  it('uses the AST diff function when provided', async () => {
    const diff = jest.fn(() => ({
      added: [],
      removed: [],
      modified: [{ startLine: 10, endLine: 11, startColumn: 0, endColumn: 1 }],
    }));
    const original = ['A', 'B', 'C', 'D', 'E', 'F'].join('\n');
    const updated = ['A', 'B2', 'C', 'D', 'E', 'F'].join('\n');
    const parser = new IncrementalParser((content: string) => ({ text: content }), diff, {
      minLinesForIncremental: 2,
    });

    await parser.parse(createDocument(original, 1));
    const result = await parser.parse(createDocument(updated, 2));

    expect(diff).toHaveBeenCalledWith({ text: original }, { text: updated });
    expect(result.changes.modified).toEqual([
      { startLine: 10, endLine: 11, startColumn: 0, endColumn: 1 },
    ]);
  });
});
