/**
 * Tests for DSL Authoring Context Aggregation
 *
 * Covers full support, partial support, and unsupported language scenarios,
 * plus capability detection and context bundle format validation.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/146
 */

import {
  assembleDSLAuthoringContext,
  formatDSLAuthoringContextMarkdown,
} from '../../../src/lsp/dslAuthoringContext';
import {
  DSLAuthoringContext,
  CapabilityStatus,
  DSLContextOutputFormat,
} from '../../../src/lsp/types';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function makeContext(
  languageId: string,
  serverRunning = false,
  format: DSLContextOutputFormat = 'json'
): DSLAuthoringContext {
  return assembleDSLAuthoringContext(
    `file:///test/example.${languageId}`,
    languageId,
    { line: 0, character: 5 },
    serverRunning,
    format
  );
}

// ---------------------------------------------------------------------------
// Full support scenario
// ---------------------------------------------------------------------------

describe('DSL Authoring Context - Full support (registered language with built-in data)', () => {
  describe('VASP', () => {
    let ctx: DSLAuthoringContext;

    beforeAll(() => {
      ctx = makeContext('vasp');
    });

    it('identifies the correct server', () => {
      expect(ctx.serverId).toBe('vasp-lsp');
      expect(ctx.serverName).toBe('VASP');
      expect(ctx.languageId).toBe('vasp');
      expect(ctx.stability).toBe('stable');
    });

    it('includes the document URI', () => {
      expect(ctx.documentUri).toContain('vasp');
    });

    it('provides a language description', () => {
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data).toBeDefined();
      expect(ctx.languageDescription.data!.name).toBe('VASP INCAR');
      expect(ctx.languageDescription.data!.summary).toBeTruthy();
      expect(ctx.languageDescription.data!.documentationUri).toContain('vasp.at');
    });

    it('provides section/keyword schema', () => {
      expect(ctx.sectionKeywordSchema.status).toBe('available');
      expect(ctx.sectionKeywordSchema.data).toBeDefined();
      expect(ctx.sectionKeywordSchema.data!.name).toBeTruthy();
      expect(ctx.sectionKeywordSchema.data!.description).toBeTruthy();
    });

    it('provides examples', () => {
      expect(ctx.examples.status).toBe('available');
      expect(ctx.examples.items.length).toBeGreaterThan(0);
      for (const ex of ctx.examples.items) {
        expect(ex.title).toBeTruthy();
        expect(ex.code).toBeTruthy();
      }
    });

    it('provides next-token guidance', () => {
      expect(ctx.nextTokenGuidance.status).toBe('available');
      expect(ctx.nextTokenGuidance.data).toBeDefined();
      expect(ctx.nextTokenGuidance.data!.candidates.length).toBeGreaterThan(0);
      expect(ctx.nextTokenGuidance.data!.hint).toBeTruthy();
    });

    it('reports diagnostics as unknown when server is not running', () => {
      expect(ctx.diagnostics.status).toBe('unknown');
      expect(ctx.diagnostics.items).toEqual([]);
    });

    it('reports code actions as unknown when server is not running', () => {
      expect(ctx.codeActions.status).toBe('unknown');
      expect(ctx.codeActions.items).toEqual([]);
    });

    it('includes an ISO-8601 assembledAt timestamp', () => {
      expect(ctx.assembledAt).toMatch(/^\d{4}-\d{2}-\d{2}T/);
    });
  });

  describe('Gaussian', () => {
    let ctx: DSLAuthoringContext;

    beforeAll(() => {
      ctx = makeContext('gaussian');
    });

    it('identifies the correct server', () => {
      expect(ctx.serverId).toBe('gaussian-lsp');
      expect(ctx.serverName).toBe('Gaussian');
      expect(ctx.stability).toBe('stable');
    });

    it('provides a language description', () => {
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('Gaussian Input');
    });

    it('provides examples', () => {
      expect(ctx.examples.status).toBe('available');
      expect(ctx.examples.items.length).toBeGreaterThanOrEqual(2);
    });

    it('provides next-token guidance with keyword candidates', () => {
      expect(ctx.nextTokenGuidance.status).toBe('available');
      const candidates = ctx.nextTokenGuidance.data!.candidates;
      expect(candidates).toContain('opt');
      expect(candidates).toContain('B3LYP');
    });
  });

  describe('Quantum ESPRESSO', () => {
    let ctx: DSLAuthoringContext;

    beforeAll(() => {
      ctx = makeContext('qe');
    });

    it('provides a language description with QE-specific documentation', () => {
      expect(ctx.languageDescription.data!.name).toBe('Quantum ESPRESSO Input');
      expect(ctx.languageDescription.data!.documentationUri).toContain('quantum-espresso');
    });

    it('provides keyword schema with calculation type', () => {
      expect(ctx.sectionKeywordSchema.status).toBe('available');
    });
  });

  describe('when server is running', () => {
    let ctx: DSLAuthoringContext;

    beforeAll(() => {
      ctx = makeContext('vasp', true);
    });

    it('marks serverRunning as true', () => {
      expect(ctx.serverRunning).toBe(true);
    });

    it('reports diagnostics as available (even if empty)', () => {
      expect(ctx.diagnostics.status).toBe('available');
    });

    it('reports code actions as available (even if empty)', () => {
      expect(ctx.codeActions.status).toBe('available');
    });
  });
});

// ---------------------------------------------------------------------------
// Partial support scenario
// ---------------------------------------------------------------------------

describe('DSL Authoring Context - Partial support (registered language without all built-in data)', () => {
  describe('LAMMPS (registered but no keyword schemas or examples)', () => {
    let ctx: DSLAuthoringContext;

    beforeAll(() => {
      ctx = makeContext('lammps');
    });

    it('identifies the correct server', () => {
      expect(ctx.serverId).toBe('lammps-lsp');
      expect(ctx.serverName).toBe('LAMMPS');
    });

    it('provides a language description from built-in or generated data', () => {
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('LAMMPS Input');
    });

    it('reports section/keyword schema as unavailable', () => {
      expect(ctx.sectionKeywordSchema.status).toBe('unavailable');
      expect(ctx.sectionKeywordSchema.data).toBeUndefined();
    });

    it('reports examples as unavailable', () => {
      expect(ctx.examples.status).toBe('unavailable');
      expect(ctx.examples.items).toEqual([]);
    });

    it('reports next-token guidance as unavailable', () => {
      expect(ctx.nextTokenGuidance.status).toBe('unavailable');
    });
  });

  describe('GPUMD (registered, experimental, partial data)', () => {
    let ctx: DSLAuthoringContext;

    beforeAll(() => {
      ctx = makeContext('gpumd');
    });

    it('marks stability as experimental', () => {
      expect(ctx.stability).toBe('experimental');
    });

    it('provides a language description', () => {
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('GPUMD Input');
    });

    it('reports examples as unavailable', () => {
      expect(ctx.examples.status).toBe('unavailable');
    });
  });

  describe('ABACUS (registered, no keyword schemas or examples)', () => {
    let ctx: DSLAuthoringContext;

    beforeAll(() => {
      ctx = makeContext('abacus');
    });

    it('provides a language description', () => {
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('ABACUS Input');
    });

    it('reports section/keyword schema as unavailable', () => {
      expect(ctx.sectionKeywordSchema.status).toBe('unavailable');
    });
  });

  describe('Languages with generated descriptions (no explicit built-in)', () => {
    it('generates description from registry for PySCF', () => {
      const ctx = makeContext('pyscf');
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('PySCF Input');
      expect(ctx.languageDescription.data!.summary).toContain('PySCF');
    });

    it('generates description from registry for MLIP', () => {
      const ctx = makeContext('mlip');
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('MLIP Configuration');
      expect(ctx.languageDescription.data!.summary).toContain('interatomic potential');
    });

    it('generates description from registry for PyATB', () => {
      const ctx = makeContext('pyatb');
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('PyATB Input');
    });

    it('generates description from registry for GROMACS', () => {
      const ctx = makeContext('gromacs');
      expect(ctx.languageDescription.status).toBe('available');
      expect(ctx.languageDescription.data!.name).toBe('GROMACS Input');
    });
  });
});

// ---------------------------------------------------------------------------
// Unsupported language scenario
// ---------------------------------------------------------------------------

describe('DSL Authoring Context - Unsupported language', () => {
  let ctx: DSLAuthoringContext;

  beforeAll(() => {
    ctx = makeContext('plaintext');
  });

  it('uses "unknown" server id', () => {
    expect(ctx.serverId).toBe('unknown');
    expect(ctx.serverName).toBe('plaintext');
  });

  it('reports language description as unavailable', () => {
    expect(ctx.languageDescription.status).toBe('unavailable');
  });

  it('reports section/keyword schema as unavailable', () => {
    expect(ctx.sectionKeywordSchema.status).toBe('unavailable');
  });

  it('reports examples as unavailable', () => {
    expect(ctx.examples.status).toBe('unavailable');
    expect(ctx.examples.items).toEqual([]);
  });

  it('reports next-token guidance as unavailable', () => {
    expect(ctx.nextTokenGuidance.status).toBe('unavailable');
  });

  it('defaults stability to experimental', () => {
    expect(ctx.stability).toBe('experimental');
  });
});

// ---------------------------------------------------------------------------
// Capability detection
// ---------------------------------------------------------------------------

describe('Capability detection', () => {
  it('reports available for all capabilities when VASP server is running', () => {
    const ctx = makeContext('vasp', true);

    expect(ctx.languageDescription.status).toBe('available');
    expect(ctx.sectionKeywordSchema.status).toBe('available');
    expect(ctx.examples.status).toBe('available');
    expect(ctx.nextTokenGuidance.status).toBe('available');
    expect(ctx.diagnostics.status).toBe('available');
    expect(ctx.codeActions.status).toBe('available');
  });

  it('reports unknown for diagnostics and code actions when server is not running', () => {
    const ctx = makeContext('gaussian', false);

    expect(ctx.diagnostics.status).toBe('unknown');
    expect(ctx.codeActions.status).toBe('unknown');
  });

  it('degrades gracefully for an unknown language with server not running', () => {
    const ctx = makeContext('unknown-lang', false);

    expect(ctx.languageDescription.status).toBe('unavailable');
    expect(ctx.sectionKeywordSchema.status).toBe('unavailable');
    expect(ctx.examples.status).toBe('unavailable');
    expect(ctx.nextTokenGuidance.status).toBe('unavailable');
    expect(ctx.diagnostics.status).toBe('unknown');
    expect(ctx.codeActions.status).toBe('unknown');
  });
});

// ---------------------------------------------------------------------------
// Context bundle format
// ---------------------------------------------------------------------------

describe('Context bundle format', () => {
  it('always includes all top-level fields', () => {
    const ctx = makeContext('vasp');

    const requiredKeys: Array<keyof DSLAuthoringContext> = [
      'documentUri',
      'languageId',
      'serverId',
      'serverName',
      'stability',
      'serverRunning',
      'languageDescription',
      'sectionKeywordSchema',
      'examples',
      'nextTokenGuidance',
      'diagnostics',
      'codeActions',
      'assembledAt',
    ];

    for (const key of requiredKeys) {
      expect(ctx).toHaveProperty(key);
    }
  });

  it('always includes status in every capability field', () => {
    const ctx = makeContext('vasp');

    const capabilityFields = [
      'languageDescription',
      'sectionKeywordSchema',
      'examples',
      'nextTokenGuidance',
      'diagnostics',
      'codeActions',
    ] as const;

    const validStatuses: CapabilityStatus[] = ['available', 'unavailable', 'unknown'];

    for (const field of capabilityFields) {
      const cap = ctx[field] as { status: CapabilityStatus };
      expect(validStatuses).toContain(cap.status);
    }
  });

  it('produces valid JSON output', () => {
    const ctx = makeContext('gaussian');
    const serialized = JSON.stringify(ctx);
    const parsed = JSON.parse(serialized);
    expect(parsed.serverId).toBe('gaussian-lsp');
    expect(parsed.languageDescription.status).toBe('available');
  });

  it('produces the same structure for json and markdown format requests', () => {
    const jsonCtx = makeContext('vasp', false, 'json');
    const mdCtx = makeContext('vasp', false, 'markdown');

    // Both return DSLAuthoringContext; markdown formatter is separate.
    expect(jsonCtx.serverId).toBe(mdCtx.serverId);
    expect(jsonCtx.languageDescription.status).toBe(mdCtx.languageDescription.status);
  });
});

// ---------------------------------------------------------------------------
// Markdown formatting
// ---------------------------------------------------------------------------

describe('Markdown formatting', () => {
  it('includes language name and server info', () => {
    const ctx = makeContext('vasp');
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('# DSL Authoring Context: VASP');
    expect(md).toContain('**Language**: vasp');
    expect(md).toContain('**Server**: vasp-lsp');
  });

  it('includes language description section', () => {
    const ctx = makeContext('gaussian');
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('## Language Description');
    expect(md).toContain('Gaussian Input');
  });

  it('shows unavailable marker for missing capabilities', () => {
    const ctx = makeContext('lammps');
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('## Section / Keyword Schema');
    expect(md).toContain('*unavailable*');
  });

  it('includes examples when available', () => {
    const ctx = makeContext('vasp');
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('## Examples');
    expect(md).toContain('### Basic relaxation');
  });

  it('includes next-token guidance candidates', () => {
    const ctx = makeContext('vasp');
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('## Next-Token Guidance');
    expect(md).toContain('ENCUT');
  });

  it('renders diagnostics section with unknown status', () => {
    const ctx = makeContext('vasp', false);
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('## Diagnostics');
    expect(md).toContain('*unknown*');
  });

  it('renders code actions section', () => {
    const ctx = makeContext('vasp', false);
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('## Code Actions');
  });

  it('handles completely unsupported language gracefully', () => {
    const ctx = makeContext('some-unknown-lang');
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('## Language Description');
    expect(md).toContain('*unavailable*');
    expect(md).toContain('## Section / Keyword Schema');
    expect(md).toContain('## Examples');
    expect(md).toContain('## Next-Token Guidance');
    expect(md).toContain('## Diagnostics');
    expect(md).toContain('## Code Actions');
  });

  it('renders diagnostics items when present', () => {
    const ctx: DSLAuthoringContext = {
      ...makeContext('vasp', true),
      diagnostics: {
        status: 'available',
        items: [
          { line: 5, character: 10, message: 'Unknown keyword XYZ', severity: 'error' },
          { line: 12, message: 'Deprecated parameter', severity: 'warning' },
        ],
      },
    };
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('[error] line 5: Unknown keyword XYZ');
    expect(md).toContain('[warning] line 12: Deprecated parameter');
  });

  it('renders code action items when present', () => {
    const ctx: DSLAuthoringContext = {
      ...makeContext('gaussian', true),
      codeActions: {
        status: 'available',
        items: [
          { title: 'Add missing basis set', kind: 'quickfix', isPreferred: true },
          { title: 'Convert to tight convergence', kind: 'refactor' },
        ],
      },
    };
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('Add missing basis set (quickfix)');
    expect(md).toContain('Convert to tight convergence (refactor)');
  });

  it('renders schema with type, default, and allowed values', () => {
    const ctx: DSLAuthoringContext = {
      ...makeContext('vasp'),
      sectionKeywordSchema: {
        status: 'available',
        data: {
          name: 'ISMEAR',
          description: 'Smearing method',
          type: 'integer',
          defaultValue: '1',
          allowedValues: ['-5', '0', '1'],
        },
      },
    };
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('**ISMEAR**');
    expect(md).toContain('Type: integer');
    expect(md).toContain('Default: 1');
    expect(md).toContain('Allowed: -5, 0, 1');
  });

  it('renders schema without type/default/allowedValues', () => {
    const ctx: DSLAuthoringContext = {
      ...makeContext('gaussian'),
      sectionKeywordSchema: {
        status: 'available',
        data: {
          name: 'opt',
          description: 'Requests geometry optimization.',
        },
      },
    };
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('**opt**');
    expect(md).toContain('Requests geometry optimization');
    expect(md).not.toContain('Type:');
    expect(md).not.toContain('Default:');
  });

  it('renders language description without documentation URI', () => {
    const ctx: DSLAuthoringContext = {
      ...makeContext('lammps'),
      languageDescription: {
        status: 'available',
        data: {
          name: 'LAMMPS Input',
          summary: 'LAMMPS input script format.',
        },
      },
    };
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('**LAMMPS Input**');
    expect(md).toContain('LAMMPS input script format');
    expect(md).not.toContain('Documentation:');
  });

  it('renders code action without kind', () => {
    const ctx: DSLAuthoringContext = {
      ...makeContext('vasp', true),
      codeActions: {
        status: 'available',
        items: [{ title: 'Fix indentation' }],
      },
    };
    const md = formatDSLAuthoringContextMarkdown(ctx);

    expect(md).toContain('- Fix indentation');
    expect(md).not.toContain('Fix indentation (');
  });
});

// ---------------------------------------------------------------------------
// All registered languages produce valid contexts
// ---------------------------------------------------------------------------

describe('All registered languages produce valid contexts', () => {
  const languageIds = [
    'abacus',
    'abinit',
    'cif',
    'cp2k',
    'vasp',
    'gaussian',
    'orca',
    'qe',
    'gamess',
    'nwchem',
    'gpumd',
    'gromacs',
    'lammps',
    'mlip',
    'pyatb',
    'pyscf',
  ];

  it.each(languageIds)('produces a valid context for language "%s"', languageId => {
    const ctx = makeContext(languageId);

    expect(ctx).toBeDefined();
    expect(ctx.languageId).toBe(languageId);
    expect(ctx.documentUri).toContain(languageId);
    expect(ctx.assembledAt).toBeTruthy();

    // Every context must have languageDescription with a valid status.
    const validStatuses: CapabilityStatus[] = ['available', 'unavailable', 'unknown'];
    expect(validStatuses).toContain(ctx.languageDescription.status);

    // If language description is available, it must have a name.
    if (ctx.languageDescription.status === 'available') {
      expect(ctx.languageDescription.data!.name).toBeTruthy();
    }
  });
});
