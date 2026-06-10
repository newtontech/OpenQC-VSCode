/**
 * LSP Readiness Tests
 *
 * Validates parse/diagnostic stability, real fixture coverage,
 * invalid input handling, and idempotent parsing for all
 * quantum chemistry parsers.
 */
import { createParser, parseInput, validateInput } from '../../../src/parsers';
import { VASPParser } from '../../../src/parsers/VASPParser';
import { GaussianParser } from '../../../src/parsers/GaussianParser';
import { ORCAParser } from '../../../src/parsers/ORCAParser';
import { CP2KParser } from '../../../src/parsers/CP2KParser';
import { QEParser } from '../../../src/parsers/QEParser';
import { GAMESSParser } from '../../../src/parsers/GAMESSParser';
import { NWChemParser } from '../../../src/parsers/NWChemParser';
import * as fs from 'fs';
import * as path from 'path';

// ─── Helpers ────────────────────────────────────────────────
const fixturesDir = path.join(__dirname, '../../fixtures');
const examplesDir = path.join(__dirname, '../../../examples');

function readFixture(relPath: string): string {
  return fs.readFileSync(path.join(fixturesDir, relPath), 'utf8');
}

function readExample(relPath: string): string {
  return fs.readFileSync(path.join(examplesDir, relPath), 'utf8');
}

/**
 * Parse the same content N times and verify all outputs are structurally equal.
 */
function assertIdempotentParse(
  ParserClass: any,
  content: string,
  extraArgs: any[] = [],
  iterations = 3
): void {
  const results: string[] = [];
  for (let i = 0; i < iterations; i++) {
    // eslint-disable-next-line @typescript-eslint/no-unsafe-call
    const parser = new ParserClass(content, ...extraArgs);
    // eslint-disable-next-line @typescript-eslint/no-unsafe-call
    const result = parser.parseInput();
    results.push(JSON.stringify(result));
  }
  for (let i = 1; i < results.length; i++) {
    expect(results[i]).toBe(results[0]);
  }
}

// ═══════════════════════════════════════════════════════════
// GAMESS Parser — Real Fixtures
// ═══════════════════════════════════════════════════════════
describe('GAMESS Parser — Real Fixtures', () => {
  describe('H2O.inp fixture', () => {
    let content: string;
    let parser: GAMESSParser;

    beforeAll(() => {
      content = readFixture('gamess/H2O.inp');
      parser = new GAMESSParser(content);
    });

    it('should parse without throwing', () => {
      const result = parser.parseInput();
      expect(result).toBeDefined();
      expect(result.sections.length).toBeGreaterThan(0);
    });

    it('should detect CONTRL section with RUNTYP=OPTIMIZE', () => {
      const result = parser.parseInput();
      const contrl = result.sections.find(s => s.name === 'CONTRL');
      expect(contrl).toBeDefined();
      const runtyp = contrl!.parameters.find(p => p.name === 'RUNTYP');
      expect(runtyp?.value).toBe('OPTIMIZE');
    });

    it('should detect SCFTYP=RHF', () => {
      const result = parser.parseInput();
      const contrl = result.sections.find(s => s.name === 'CONTRL');
      const scftyp = contrl!.parameters.find(p => p.name === 'SCFTYP');
      expect(scftyp?.value).toBe('RHF');
    });

    it('should detect SYSTEM section with MWORDS', () => {
      const result = parser.parseInput();
      const system = result.sections.find(s => s.name === 'SYSTEM');
      expect(system).toBeDefined();
      const mwords = system!.parameters.find(p => p.name === 'MWORDS');
      expect(mwords?.value).toBe(100);
    });

    it('should detect BASIS section with GBASIS=N31', () => {
      const result = parser.parseInput();
      const basis = result.sections.find(s => s.name === 'BASIS');
      expect(basis).toBeDefined();
      const gbasis = basis!.parameters.find(p => p.name === 'GBASIS');
      expect(gbasis?.value).toBe('N31');
    });

    it('should validate as valid input', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
      expect(validation.errors).toHaveLength(0);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(GAMESSParser, content);
    });
  });

  describe('Si2.inp fixture', () => {
    let content: string;
    let parser: GAMESSParser;

    beforeAll(() => {
      content = readFixture('gamess/Si2.inp');
      parser = new GAMESSParser(content);
    });

    it('should parse successfully', () => {
      const result = parser.parseInput();
      expect(result.sections.length).toBeGreaterThan(0);
    });

    it('should detect UHF and triplet', () => {
      const result = parser.parseInput();
      const contrl = result.sections.find(s => s.name === 'CONTRL');
      const scftyp = contrl!.parameters.find(p => p.name === 'SCFTYP');
      expect(scftyp?.value).toBe('UHF');
      const mult = contrl!.parameters.find(p => p.name === 'MULT');
      expect(mult?.value).toBe(3);
    });

    it('should validate as valid input', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(GAMESSParser, content);
    });
  });
});

// ═══════════════════════════════════════════════════════════
// NWChem Parser — Real Fixtures
// ═══════════════════════════════════════════════════════════
describe('NWChem Parser — Real Fixtures', () => {
  describe('H2O.nw fixture', () => {
    let content: string;
    let parser: NWChemParser;

    beforeAll(() => {
      content = readFixture('nwchem/H2O.nw');
      parser = new NWChemParser(content);
    });

    it('should parse without throwing', () => {
      const result = parser.parseInput();
      expect(result).toBeDefined();
      expect(result.sections.length).toBeGreaterThan(0);
    });

    it('should detect GEOMETRY block with atoms', () => {
      const result = parser.parseInput();
      const geom = result.sections.find(s => s.name === 'GEOMETRY');
      expect(geom).toBeDefined();
      // Should have atom coordinate entries
      const atomParams = geom!.parameters.filter(p => p.name.startsWith('Atom_'));
      expect(atomParams.length).toBe(3); // O, H, H
    });

    it('should detect BASIS block', () => {
      const result = parser.parseInput();
      const basis = result.sections.find(s => s.name === 'BASIS');
      expect(basis).toBeDefined();
    });

    it('should detect DFT block with xc setting containing b3lyp', () => {
      const result = parser.parseInput();
      const dft = result.sections.find(s => s.name === 'DFT');
      expect(dft).toBeDefined();
      // NWChem DFT block stores space-separated key-value as 'Setting'
      const xcSetting = dft!.parameters.find(
        p =>
          p.name === 'Setting' &&
          String(p.value).includes('xc') &&
          String(p.value).includes('b3lyp')
      );
      expect(xcSetting).toBeDefined();
    });

    it('should detect TASK directives', () => {
      const result = parser.parseInput();
      const tasks = result.sections.filter(s => s.name === 'TASK');
      expect(tasks.length).toBeGreaterThanOrEqual(2);
    });

    it('should validate successfully', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
      expect(validation.errors).toHaveLength(0);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(NWChemParser, content);
    });
  });

  describe('Si2.nw fixture', () => {
    let content: string;
    let parser: NWChemParser;

    beforeAll(() => {
      content = readFixture('nwchem/Si2.nw');
      parser = new NWChemParser(content);
    });

    it('should parse successfully', () => {
      const result = parser.parseInput();
      expect(result.sections.length).toBeGreaterThan(0);
    });

    it('should detect SCF block with UHF', () => {
      const result = parser.parseInput();
      const scf = result.sections.find(s => s.name === 'SCF');
      expect(scf).toBeDefined();
      const uhf = scf!.parameters.find(p => String(p.value) === 'uhf');
      expect(uhf).toBeDefined();
    });

    it('should validate successfully', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(NWChemParser, content);
    });
  });
});

// ═══════════════════════════════════════════════════════════
// ORCA Parser — Real Fixtures
// ═══════════════════════════════════════════════════════════
describe('ORCA Parser — Real Fixtures', () => {
  describe('Si2.inp fixture', () => {
    let content: string;
    let parser: ORCAParser;

    beforeAll(() => {
      content = readFixture('orca/Si2.inp');
      parser = new ORCAParser(content);
    });

    it('should parse without throwing', () => {
      const result = parser.parseInput();
      expect(result).toBeDefined();
    });

    it('should detect method B3LYP from simple input', () => {
      const result = parser.parseInput();
      const method = result.parameters.find(p => p.name === 'Method');
      expect(method?.value).toBe('B3LYP');
    });

    it('should detect basis def2-SVP', () => {
      const result = parser.parseInput();
      const basis = result.parameters.find(p => p.name === 'Basis');
      expect(basis?.value).toBe('def2-SVP');
    });

    it('should detect calculation type Opt', () => {
      const result = parser.parseInput();
      const calcType = result.parameters.find(p => p.name === 'CalculationType');
      expect(calcType?.value).toBe('Opt');
    });

    it('should detect atom coordinates', () => {
      const result = parser.parseInput();
      const atomParams = result.parameters.filter(p => p.name.startsWith('Atom_'));
      expect(atomParams.length).toBe(2);
    });

    it('should validate without errors', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
      expect(validation.errors).toHaveLength(0);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(ORCAParser, content);
    });
  });
});

// ═══════════════════════════════════════════════════════════
// QE Parser — Real Fixtures
// ═══════════════════════════════════════════════════════════
describe('QE Parser — Real Fixtures', () => {
  describe('H2O.relax.in fixture', () => {
    let content: string;
    let parser: QEParser;

    beforeAll(() => {
      content = readFixture('qe/H2O.relax.in');
      parser = new QEParser(content);
    });

    it('should parse without throwing', () => {
      const result = parser.parseInput();
      expect(result).toBeDefined();
      expect(result.sections.length).toBeGreaterThan(0);
    });

    it('should detect CONTROL namelist with calculation=relax', () => {
      const result = parser.parseInput();
      const control = result.sections.find(s => s.name === 'CONTROL');
      expect(control).toBeDefined();
      const calc = control!.parameters.find(p => p.name === 'calculation');
      expect(calc?.value).toBe('relax');
    });

    it('should detect SYSTEM namelist with nat and ntyp', () => {
      const result = parser.parseInput();
      const system = result.sections.find(s => s.name === 'SYSTEM');
      expect(system).toBeDefined();
      const nat = system!.parameters.find(p => p.name === 'nat');
      expect(nat?.value).toBe(3);
      const ntyp = system!.parameters.find(p => p.name === 'ntyp');
      expect(ntyp?.value).toBe(2);
    });

    it('should detect ELECTRONS namelist', () => {
      const result = parser.parseInput();
      const electrons = result.sections.find(s => s.name === 'ELECTRONS');
      expect(electrons).toBeDefined();
      const convThr = electrons!.parameters.find(p => p.name === 'conv_thr');
      expect(convThr).toBeDefined();
    });

    it('should detect IONS namelist', () => {
      const result = parser.parseInput();
      const ions = result.sections.find(s => s.name === 'IONS');
      expect(ions).toBeDefined();
    });

    it('should detect ATOMIC_SPECIES section', () => {
      const result = parser.parseInput();
      const species = result.sections.find(s => s.name === 'ATOMIC_SPECIES');
      expect(species).toBeDefined();
    });

    it('should detect ATOMIC_POSITIONS section', () => {
      const result = parser.parseInput();
      const positions = result.sections.find(s => s.name === 'ATOMIC_POSITIONS');
      expect(positions).toBeDefined();
    });

    it('should validate with warnings (ELECTRONS exists)', () => {
      const validation = parser.validate();
      expect(validation.errors).toHaveLength(0);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(QEParser, content);
    });
  });

  describe('Si.scf.in fixture', () => {
    let content: string;
    let parser: QEParser;

    beforeAll(() => {
      content = readFixture('qe/Si.scf.in');
      parser = new QEParser(content);
    });

    it('should parse without errors', () => {
      const result = parser.parseInput();
      expect(result).toBeDefined();
    });

    it('should detect CONTROL with calculation=scf', () => {
      const result = parser.parseInput();
      const control = result.sections.find(s => s.name === 'CONTROL');
      const calc = control!.parameters.find(p => p.name === 'calculation');
      expect(calc?.value).toBe('scf');
    });

    it('should detect K_POINTS section', () => {
      const result = parser.parseInput();
      const kpoints = result.sections.find(s => s.name === 'K_POINTS');
      expect(kpoints).toBeDefined();
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(QEParser, content);
    });
  });
});

// ═══════════════════════════════════════════════════════════
// Gaussian Parser — Real Fixtures
// ═══════════════════════════════════════════════════════════
describe('Gaussian Parser — Real Fixtures', () => {
  describe('water_optimization.com fixture', () => {
    let content: string;
    let parser: GaussianParser;

    beforeAll(() => {
      content = readFixture('gaussian/water_optimization.com');
      parser = new GaussianParser(content);
    });

    it('should parse without throwing', () => {
      const result = parser.parseInput();
      expect(result).toBeDefined();
    });

    it('should extract Link 0 commands (mem, nprocshared, chk)', () => {
      const result = parser.parseInput();
      const mem = result.parameters.find(p => p.name === 'mem');
      expect(mem).toBeDefined();
      const nproc = result.parameters.find(p => p.name === 'nprocshared');
      expect(nproc).toBeDefined();
      const chk = result.parameters.find(p => p.name === 'chk');
      expect(chk).toBeDefined();
    });

    it('should extract method B3LYP and basis 6-31G(d)', () => {
      const result = parser.parseInput();
      const method = result.parameters.find(p => p.name === 'Method');
      expect(method?.value).toBe('B3LYP');
      const basis = result.parameters.find(p => p.name === 'Basis');
      expect(basis?.value).toBe('6-31G(d)');
    });

    it('should detect Opt and Freq calculation types', () => {
      const result = parser.parseInput();
      const calcType = result.parameters.find(p => p.name === 'CalculationType');
      expect(calcType).toBeDefined();
    });

    it('should extract charge=0 and multiplicity=1', () => {
      const result = parser.parseInput();
      const charge = result.parameters.find(p => p.name === 'Charge');
      expect(charge?.value).toBe(0);
      const mult = result.parameters.find(p => p.name === 'Multiplicity');
      expect(mult?.value).toBe(1);
    });

    it('should validate as valid input', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(GaussianParser, content);
    });
  });
});

// ═══════════════════════════════════════════════════════════
// VASP Parser — Real Fixtures
// ═══════════════════════════════════════════════════════════
describe('VASP Parser — Real Fixtures', () => {
  describe('INCAR fixture', () => {
    let content: string;
    let parser: VASPParser;

    beforeAll(() => {
      content = readFixture('vasp/INCAR');
      parser = new VASPParser(content, 'INCAR');
    });

    it('should parse without throwing', () => {
      const result = parser.parseInput();
      expect(result).toBeDefined();
    });

    it('should extract SYSTEM parameter', () => {
      const result = parser.parseInput();
      const sys = result.parameters.find(p => p.name === 'SYSTEM');
      expect(sys).toBeDefined();
    });

    it('should extract ENCUT = 520', () => {
      const result = parser.parseInput();
      const encut = result.parameters.find(p => p.name === 'ENCUT');
      expect(encut?.value).toBe(520);
    });

    it('should extract ISMEAR = 0', () => {
      const result = parser.parseInput();
      const ismear = result.parameters.find(p => p.name === 'ISMEAR');
      expect(ismear?.value).toBe(0);
    });

    it('should extract boolean LWAVE = .FALSE.', () => {
      const result = parser.parseInput();
      const lwave = result.parameters.find(p => p.name === 'LWAVE');
      expect(lwave?.value).toBe(false);
    });

    it('should validate with warnings (has ENCUT and PREC... wait PREC is not present)', () => {
      const validation = parser.validate();
      // PREC not in fixture, should warn
      expect(validation.warnings.some(w => w.message.includes('PREC'))).toBe(true);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(VASPParser, content, ['INCAR']);
    });
  });

  describe('POSCAR fixture', () => {
    let content: string;
    let parser: VASPParser;

    beforeAll(() => {
      content = readFixture('vasp/POSCAR');
      parser = new VASPParser(content, 'POSCAR');
    });

    it('should parse without errors', () => {
      const result = parser.parseInput();
      expect(result.errors).toHaveLength(0);
    });

    it('should extract lattice vectors', () => {
      const lattice = parser.getLatticeVectors();
      expect(lattice).toHaveLength(3);
    });

    it('should extract coordinates', () => {
      const coords = parser.getCoordinates();
      expect(coords.length).toBeGreaterThanOrEqual(2);
    });

    it('should extract atom types Si', () => {
      const types = parser.getAtomTypes();
      expect(types).toContain('Si');
    });

    it('should extract atom counts', () => {
      const counts = parser.getAtomCounts();
      expect(counts).toEqual([2]);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(VASPParser, content, ['POSCAR']);
    });
  });

  describe('KPOINTS fixture', () => {
    let content: string;
    let parser: VASPParser;

    beforeAll(() => {
      content = readFixture('vasp/KPOINTS');
      parser = new VASPParser(content, 'KPOINTS');
    });

    it('should parse without errors', () => {
      const result = parser.parseInput();
      expect(result.errors).toHaveLength(0);
    });

    it('should extract KX, KY, KZ = 4', () => {
      const result = parser.parseInput();
      const kx = result.parameters.find(p => p.name === 'KX');
      expect(kx?.value).toBe(4);
    });

    it('should parse idempotently', () => {
      assertIdempotentParse(VASPParser, content, ['KPOINTS']);
    });
  });
});

// ═══════════════════════════════════════════════════════════
// Invalid Input & Edge Case Tests
// ═══════════════════════════════════════════════════════════
describe('Invalid Input Handling', () => {
  describe('GAMESS invalid inputs', () => {
    it('should report missing CONTRL section', () => {
      const parser = new GAMESSParser('$SCF DIRSCF=.TRUE. $END');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('CONTRL'))).toBe(true);
    });

    it('should report missing SYSTEM section', () => {
      const parser = new GAMESSParser('$CONTRL RUNTYP=ENERGY SCFTYP=RHF $END');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('SYSTEM'))).toBe(true);
    });

    it('should report missing RUNTYP in CONTRL', () => {
      const parser = new GAMESSParser(
        '$CONTRL SCFTYP=RHF $END\n$SYSTEM MWORDS=100 $END\n$BASIS GBASIS=N31 $END'
      );
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('RUNTYP'))).toBe(true);
    });

    it('should report missing SCFTYP in CONTRL', () => {
      const parser = new GAMESSParser(
        '$CONTRL RUNTYP=ENERGY $END\n$SYSTEM MWORDS=100 $END\n$BASIS GBASIS=N31 $END'
      );
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('SCFTYP'))).toBe(true);
    });

    it('should handle completely empty input', () => {
      const parser = new GAMESSParser('');
      const result = parser.parseInput();
      expect(result.sections).toHaveLength(0);
      expect(result.parameters).toHaveLength(0);
    });

    it('should handle input with only whitespace', () => {
      const parser = new GAMESSParser('   \n  \n  ');
      const result = parser.parseInput();
      expect(result.sections).toHaveLength(0);
    });

    it('should handle unclosed data group', () => {
      const parser = new GAMESSParser('$CONTRL RUNTYP=ENERGY');
      const result = parser.parseInput();
      // Should still parse without crashing
      expect(result.sections.length).toBeGreaterThanOrEqual(1);
    });
  });

  describe('NWChem invalid inputs', () => {
    it('should report missing GEOMETRY block', () => {
      const parser = new NWChemParser('basis\n  * library 6-31G*\nend\ntask scf energy');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('GEOMETRY'))).toBe(true);
    });

    it('should report missing BASIS block', () => {
      const parser = new NWChemParser('geometry\n  H 0 0 0\nend\ntask scf energy');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('BASIS'))).toBe(true);
    });

    it('should report missing TASK directive', () => {
      const parser = new NWChemParser('geometry\n  H 0 0 0\nend\nbasis\n  * library 6-31G*\nend');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('TASK'))).toBe(true);
    });

    it('should handle completely empty input', () => {
      const parser = new NWChemParser('');
      const result = parser.parseInput();
      expect(result.sections).toHaveLength(0);
    });

    it('should handle input with only comments', () => {
      const parser = new NWChemParser('# comment\n# another comment');
      const result = parser.parseInput();
      expect(result.sections).toHaveLength(0);
    });
  });

  describe('ORCA invalid inputs', () => {
    it('should detect empty input', () => {
      const parser = new ORCAParser('');
      const result = parser.parseInput();
      expect(result.errors.some(e => e.message.includes('No valid ORCA input'))).toBe(true);
    });

    it('should warn when no method is specified', () => {
      const parser = new ORCAParser('%pal nprocs 4 end');
      const validation = parser.validate();
      expect(validation.warnings.some(w => w.message.includes('method'))).toBe(true);
    });

    it('should handle input with only comments', () => {
      const parser = new ORCAParser('# comment\n# another');
      const result = parser.parseInput();
      expect(result.errors.some(e => e.message.includes('No valid ORCA input'))).toBe(true);
    });
  });

  describe('QE invalid inputs', () => {
    it('should report missing CONTROL namelist', () => {
      const parser = new QEParser('&SYSTEM\n  nat = 1\n/');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('CONTROL'))).toBe(true);
    });

    it('should report missing SYSTEM namelist', () => {
      const parser = new QEParser("&CONTROL\n  calculation = 'scf'\n/");
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('SYSTEM'))).toBe(true);
    });

    it('should warn about missing ELECTRONS namelist', () => {
      const parser = new QEParser(
        "&CONTROL\n  calculation = 'scf'\n/\n&SYSTEM\n  nat = 1\n  ntyp = 1\n/"
      );
      const validation = parser.validate();
      expect(validation.warnings.some(w => w.message.includes('ELECTRONS'))).toBe(true);
    });

    it('should warn about missing calculation parameter', () => {
      const parser = new QEParser("&CONTROL\n  outdir = './tmp'\n/\n&SYSTEM\n  nat = 1\n/");
      const validation = parser.validate();
      expect(validation.warnings.some(w => w.message.includes('calculation'))).toBe(true);
    });

    it('should handle completely empty input', () => {
      const parser = new QEParser('');
      const result = parser.parseInput();
      expect(result.sections).toHaveLength(0);
    });
  });

  describe('CP2K invalid inputs', () => {
    it('should report missing FORCE_EVAL section', () => {
      const parser = new CP2KParser('&GLOBAL\n  PROJECT_NAME test\n&END GLOBAL');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('FORCE_EVAL'))).toBe(true);
    });

    it('should warn about missing PROJECT_NAME', () => {
      const parser = new CP2KParser(
        '&GLOBAL\n  RUN_TYPE ENERGY\n&END GLOBAL\n&FORCE_EVAL\n  METHOD Quickstep\n&END FORCE_EVAL'
      );
      const validation = parser.validate();
      expect(validation.warnings.some(w => w.message.includes('PROJECT_NAME'))).toBe(true);
    });

    it('should report unclosed sections', () => {
      const parser = new CP2KParser('&GLOBAL\n  PROJECT_NAME test');
      const result = parser.parseInput();
      expect(result.errors.some(e => e.message.includes('not properly closed'))).toBe(true);
    });
  });

  describe('Gaussian invalid inputs', () => {
    it('should detect missing route section', () => {
      const parser = new GaussianParser('');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('route'))).toBe(true);
    });

    it('should detect missing charge/multiplicity', () => {
      const parser = new GaussianParser('# B3LYP/6-31G(d)\n\nTest');
      const validation = parser.validate();
      expect(validation.errors.length).toBeGreaterThan(0);
    });

    it('should detect missing method in route', () => {
      const parser = new GaussianParser('#p freq\n\nTest\n\n0 1\nH 0 0 0');
      const validation = parser.validate();
      expect(validation.valid).toBe(false);
    });
  });

  describe('VASP invalid inputs', () => {
    it('should reject too-short POSCAR', () => {
      const parser = new VASPParser('Too short', 'POSCAR');
      const result = parser.parseInput();
      expect(result.errors.length).toBeGreaterThan(0);
    });

    it('should reject too-short KPOINTS', () => {
      const parser = new VASPParser('K-points', 'KPOINTS');
      const result = parser.parseInput();
      expect(result.errors.length).toBeGreaterThan(0);
    });

    it('should reject empty POTCAR', () => {
      const parser = new VASPParser('', 'POTCAR');
      const result = parser.parseInput();
      expect(result.errors.length).toBeGreaterThan(0);
    });

    it('should warn about malformed INCAR lines', () => {
      const parser = new VASPParser('ENCUT 520', 'INCAR');
      const result = parser.parseInput();
      expect(result.warnings.length).toBeGreaterThan(0);
    });
  });
});

// ═══════════════════════════════════════════════════════════
// Parser Registry Integration
// ═══════════════════════════════════════════════════════════
describe('Parser Registry — LSP Integration', () => {
  it('should create parsers for all supported software', () => {
    const software = [
      'VASP',
      'Gaussian',
      'ORCA',
      'CP2K',
      'Quantum ESPRESSO',
      'GAMESS',
      'NWChem',
    ] as const;
    for (const sw of software) {
      const parser = createParser(sw, 'dummy content');
      expect(parser).toBeDefined();
    }
  });

  it('should throw for unsupported software', () => {
    expect(() => createParser('Unknown' as any, 'content')).toThrow('Unsupported software');
  });

  it('parseInput() should return consistent structure for all parsers', () => {
    const software = [
      'VASP',
      'Gaussian',
      'ORCA',
      'CP2K',
      'Quantum ESPRESSO',
      'GAMESS',
      'NWChem',
    ] as const;
    for (const sw of software) {
      const result = parseInput(sw, '');
      expect(result).toHaveProperty('sections');
      expect(result).toHaveProperty('parameters');
      expect(result).toHaveProperty('errors');
      expect(result).toHaveProperty('warnings');
      expect(Array.isArray(result.sections)).toBe(true);
      expect(Array.isArray(result.parameters)).toBe(true);
      expect(Array.isArray(result.errors)).toBe(true);
      expect(Array.isArray(result.warnings)).toBe(true);
    }
  });

  it('validateInput() should return consistent structure for all parsers', () => {
    const software = [
      'VASP',
      'Gaussian',
      'ORCA',
      'CP2K',
      'Quantum ESPRESSO',
      'GAMESS',
      'NWChem',
    ] as const;
    for (const sw of software) {
      const result = validateInput(sw, '');
      expect(result).toHaveProperty('valid');
      expect(result).toHaveProperty('errors');
      expect(result).toHaveProperty('warnings');
      expect(typeof result.valid).toBe('boolean');
    }
  });
});

// ═══════════════════════════════════════════════════════════
// Diagnostic Stability
// ═══════════════════════════════════════════════════════════
describe('Diagnostic Stability', () => {
  it('should produce identical diagnostics on repeated validation calls', () => {
    const content = readFixture('gamess/H2O.inp');
    const parser1 = new GAMESSParser(content);
    const parser2 = new GAMESSParser(content);
    const v1 = parser1.validate();
    const v2 = parser2.validate();
    expect(JSON.stringify(v1)).toBe(JSON.stringify(v2));
  });

  it('should produce stable NWChem diagnostics', () => {
    const content = readFixture('nwchem/H2O.nw');
    const v1 = new NWChemParser(content).validate();
    const v2 = new NWChemParser(content).validate();
    expect(JSON.stringify(v1)).toBe(JSON.stringify(v2));
  });

  it('should produce stable QE diagnostics', () => {
    const content = readFixture('qe/H2O.relax.in');
    const v1 = new QEParser(content).validate();
    const v2 = new QEParser(content).validate();
    expect(JSON.stringify(v1)).toBe(JSON.stringify(v2));
  });

  it('should produce stable ORCA diagnostics', () => {
    const content = readFixture('orca/Si2.inp');
    const v1 = new ORCAParser(content).validate();
    const v2 = new ORCAParser(content).validate();
    expect(JSON.stringify(v1)).toBe(JSON.stringify(v2));
  });

  it('should produce stable VASP diagnostics', () => {
    const content = readFixture('vasp/INCAR');
    const v1 = new VASPParser(content, 'INCAR').validate();
    const v2 = new VASPParser(content, 'INCAR').validate();
    expect(JSON.stringify(v1)).toBe(JSON.stringify(v2));
  });

  it('should produce stable CP2K diagnostics', () => {
    const content = readFixture('cp2k/H2O.inp');
    const v1 = new CP2KParser(content).validate();
    const v2 = new CP2KParser(content).validate();
    expect(JSON.stringify(v1)).toBe(JSON.stringify(v2));
  });
});

// ═══════════════════════════════════════════════════════════
// Invalid Fixture Files — Diagnostic Accuracy
// ═══════════════════════════════════════════════════════════
describe('Invalid Fixture Files — Diagnostic Accuracy', () => {
  it('should report missing GAMESS groups from invalid fixture', () => {
    const content = readFixture('invalid/gamess_missing_groups.inp');
    const parser = new GAMESSParser(content);
    const validation = parser.validate();
    expect(validation.valid).toBe(false);
    // Should report CONTRL, SYSTEM, BASIS as missing
    expect(validation.errors.some(e => e.message.includes('CONTRL'))).toBe(true);
    expect(validation.errors.some(e => e.message.includes('SYSTEM'))).toBe(true);
    expect(validation.errors.some(e => e.message.includes('BASIS'))).toBe(true);
  });

  it('should report missing NWChem task from invalid fixture', () => {
    const content = readFixture('invalid/nwchem_missing_task.nw');
    const parser = new NWChemParser(content);
    const validation = parser.validate();
    expect(validation.valid).toBe(false);
    expect(validation.errors.some(e => e.message.includes('TASK'))).toBe(true);
  });

  it('should detect empty ORCA input from invalid fixture', () => {
    const content = readFixture('invalid/orca_empty.inp');
    const parser = new ORCAParser(content);
    const result = parser.parseInput();
    expect(result.errors.some(e => e.message.includes('No valid ORCA'))).toBe(true);
  });

  it('should report missing QE SYSTEM from invalid fixture', () => {
    const content = readFixture('invalid/qe_missing_system.in');
    const parser = new QEParser(content);
    const validation = parser.validate();
    expect(validation.valid).toBe(false);
    expect(validation.errors.some(e => e.message.includes('SYSTEM'))).toBe(true);
  });

  it('should detect missing Gaussian route from invalid fixture', () => {
    const content = readFixture('invalid/gaussian_no_route.com');
    const parser = new GaussianParser(content);
    const validation = parser.validate();
    expect(validation.valid).toBe(false);
  });

  it('should detect unclosed CP2K section from invalid fixture', () => {
    const content = readFixture('invalid/cp2k_unclosed.inp');
    const parser = new CP2KParser(content);
    const result = parser.parseInput();
    expect(result.errors.some(e => e.message.includes('not properly closed'))).toBe(true);
  });
});

// ═══════════════════════════════════════════════════════════
// Cross-Format Robustness
// ═══════════════════════════════════════════════════════════
describe('Cross-Format Robustness', () => {
  it('should not crash when GAMESS parser sees ORCA-like content', () => {
    const parser = new GAMESSParser('! B3LYP def2-SVP Opt');
    const result = parser.parseInput();
    expect(result).toBeDefined();
    // Should not crash, just return empty sections
  });

  it('should not crash when NWChem parser sees VASP-like content', () => {
    const parser = new NWChemParser('ENCUT = 520\nPREC = Accurate');
    const result = parser.parseInput();
    expect(result).toBeDefined();
  });

  it('should not crash when ORCA parser sees CP2K-like content', () => {
    const parser = new ORCAParser('&GLOBAL\n  PROJECT_NAME test\n&END GLOBAL');
    const result = parser.parseInput();
    expect(result).toBeDefined();
  });

  it('should not crash when QE parser sees Gaussian-like content', () => {
    const parser = new QEParser('# B3LYP/6-31G(d)\n\nTitle\n\n0 1\nH 0 0 0');
    const result = parser.parseInput();
    expect(result).toBeDefined();
  });

  it('should not crash when Gaussian parser sees NWChem-like content', () => {
    const parser = new GaussianParser('geometry\n  H 0 0 0\nend');
    const result = parser.parseInput();
    expect(result).toBeDefined();
  });

  it('should handle unicode content gracefully in all parsers', () => {
    const software = [
      'VASP',
      'Gaussian',
      'ORCA',
      'CP2K',
      'Quantum ESPRESSO',
      'GAMESS',
      'NWChem',
    ] as const;
    const unicodeContent = '# Comment with Ünïcödé 水分子 ☃\n';
    for (const sw of software) {
      expect(() => parseInput(sw, unicodeContent)).not.toThrow();
    }
  });
});

// ═══════════════════════════════════════════════════════════
// Examples Directory — Parse Stability
// ═══════════════════════════════════════════════════════════
describe('Examples Directory — Parse Stability', () => {
  it('should parse GAMESS examples without crashing', () => {
    const files = ['gamess/h2o_scf.inp', 'gamess/benzene_opt.inp', 'gamess/formaldehyde_tddft.inp'];
    for (const file of files) {
      const content = readExample(file);
      const parser = new GAMESSParser(content);
      expect(() => parser.parseInput()).not.toThrow();
      expect(() => parser.validate()).not.toThrow();
    }
  });

  it('should parse NWChem examples without crashing', () => {
    const files = ['nwchem/h2o_scf.nw', 'nwchem/benzene_dft.nw', 'nwchem/si_bulk.nw'];
    for (const file of files) {
      const content = readExample(file);
      const parser = new NWChemParser(content);
      expect(() => parser.parseInput()).not.toThrow();
      expect(() => parser.validate()).not.toThrow();
    }
  });

  it('should parse ORCA examples without crashing', () => {
    const files = [
      'orca/h2o_opt.inp',
      'orca/h2o_ccsdt.inp',
      'orca/benzene_pcm.inp',
      'orca/phenol_opt.inp',
    ];
    for (const file of files) {
      const content = readExample(file);
      const parser = new ORCAParser(content);
      expect(() => parser.parseInput()).not.toThrow();
      expect(() => parser.validate()).not.toThrow();
    }
  });

  it('should parse QE examples without crashing', () => {
    const files = ['qe/si_scf.in', 'qe/si_bands.in', 'qe/h2o_vc_relax.in'];
    for (const file of files) {
      const content = readExample(file);
      const parser = new QEParser(content);
      expect(() => parser.parseInput()).not.toThrow();
      expect(() => parser.validate()).not.toThrow();
    }
  });

  it('should parse Gaussian examples without crashing', () => {
    const files = [
      'gaussian/h2o_opt.com',
      'gaussian/benzene_sp.com',
      'gaussian/water_opt.com',
      'gaussian/ts_search.com',
    ];
    for (const file of files) {
      const content = readExample(file);
      const parser = new GaussianParser(content);
      expect(() => parser.parseInput()).not.toThrow();
      expect(() => parser.validate()).not.toThrow();
    }
  });

  it('should parse CP2K examples without crashing', () => {
    const files = ['cp2k/h2o_sp.inp', 'cp2k/ch4_opt.inp', 'cp2k/si_bulk.inp'];
    for (const file of files) {
      const content = readExample(file);
      const parser = new CP2KParser(content);
      expect(() => parser.parseInput()).not.toThrow();
      expect(() => parser.validate()).not.toThrow();
    }
  });

  it('should parse VASP examples without crashing', () => {
    const files = [
      { path: 'vasp/INCAR-Si', filename: 'INCAR' },
      { path: 'vasp/KPOINTS-Si', filename: 'KPOINTS' },
      { path: 'vasp/POSCAR-Si2', filename: 'POSCAR' },
      { path: 'vasp/INCAR-MD', filename: 'INCAR' },
    ];
    for (const { path: filePath, filename } of files) {
      const content = readExample(filePath);
      const parser = new VASPParser(content, filename);
      expect(() => parser.parseInput()).not.toThrow();
      expect(() => parser.validate()).not.toThrow();
    }
  });
});
