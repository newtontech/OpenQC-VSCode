/**
 * k-Point Grid Migration Tests
 */

import {
  parseVASP_KPOINTS,
  generateVASP_KPOINTS,
  parseQE_KPOINTS,
  generateQE_KPOINTS,
  monkhorstToGamma,
  gammaToMonkhorst,
  detectKpointFormat,
  parseKpoints,
  generateKpoints,
  KpointMigration,
  KpointGrid,
} from '../../../src/utils/migration/kpoints';

describe('VASP KPOINTS Parser', () => {
  it('should parse Gamma-centered grid', () => {
    const content = 'Gamma grid\n0\nGamma\n4 4 4\n0 0 0';
    const grid = parseVASP_KPOINTS(content);
    expect(grid.grid).toEqual([4, 4, 4]);
    expect(grid.type).toBe('Gamma-centered');
    expect(grid.shift).toEqual([0, 0, 0]);
    expect(grid.comment).toBe('Gamma grid');
  });

  it('should parse Monkhorst-Pack grid', () => {
    const content = 'MP grid\n0\nMonkhorst-Pack\n6 6 6\n0 0 0';
    const grid = parseVASP_KPOINTS(content);
    expect(grid.grid).toEqual([6, 6, 6]);
    expect(grid.type).toBe('Monkhorst-Pack');
  });

  it('should throw on invalid content', () => {
    expect(() => parseVASP_KPOINTS('invalid')).toThrow();
  });
});

describe('Quantum ESPRESSO K_POINTS Parser', () => {
  it('should parse automatic k-points', () => {
    const content = 'K_POINTS automatic\n4 4 4\n0 0 0';
    const grid = parseQE_KPOINTS(content);
    expect(grid.grid).toEqual([4, 4, 4]);
    expect(grid.type).toBe('Monkhorst-Pack');
  });

  it('should parse gamma point', () => {
    const content = 'K_POINTS gamma';
    const grid = parseQE_KPOINTS(content);
    expect(grid.grid).toEqual([1, 1, 1]);
    expect(grid.type).toBe('Gamma-centered');
  });
});

describe('Grid Type Conversion', () => {
  it('should convert Monkhorst-Pack to Gamma-centered', () => {
    expect(monkhorstToGamma([4, 4, 4])).toEqual([8, 8, 8]);
  });

  it('should convert Gamma-centered to Monkhorst-Pack', () => {
    expect(gammaToMonkhorst([8, 8, 8])).toEqual([4, 4, 4]);
  });
});

describe('KpointMigration', () => {
  const migration = new KpointMigration({} as any);

  it('should check supported formats', () => {
    expect(migration.isSupported('vasp', 'qe')).toBe(true);
    expect(migration.isSupported('vasp', 'unknown')).toBe(false);
  });
});
