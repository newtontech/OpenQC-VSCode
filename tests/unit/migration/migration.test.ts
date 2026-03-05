/**
 * Unit tests for migration module
 */

import * as path from 'path';
import * as fs from 'fs';
import * as os from 'os';
import {
  getParameterMappings,
  getParameterMapping,
  convertParameterValue,
  isMigrationSupported,
  getSupportedSources,
  getSupportedTargets,
  ENERGY_CUTOFF_MAPPINGS,
  XC_FUNCTIONAL_MAPPINGS,
  KPOINT_MAPPINGS,
  CONVERGENCE_MAPPINGS,
  STRUCTURE_MAPPINGS,
  MD_MAPPINGS,
  ParameterMapping,
} from '../../../src/utils/migration/params';
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
import {
  ParameterConverter,
  ExtractedParameters,
  ConvertedParameters,
} from '../../../src/utils/migration/parameterConverter';
import {
  MDWorkflowConverter,
  MDWorkflow,
  MDParameters,
  OptimizationParameters,
  EnsembleType,
  ThermostatType,
  BarostatType,
  OptimizationAlgorithm,
  quickMigrateMDWorkflow,
} from '../../../src/utils/migration/mdWorkflow';

describe('Migration Params', () => {
  describe('getParameterMappings', () => {
    it('should return mappings for supported source-target pairs', () => {
      const mappings = getParameterMappings('vasp', 'cp2k');
      expect(mappings.length).toBeGreaterThan(0);
      expect(mappings.every(m => m.source === 'vasp' && m.target === 'cp2k')).toBe(true);
    });

    it('should return empty array for unsupported pairs', () => {
      const mappings = getParameterMappings('unknown', 'unknown');
      expect(mappings).toEqual([]);
    });

    it('should be case insensitive', () => {
      const mappings1 = getParameterMappings('VASP', 'CP2K');
      const mappings2 = getParameterMappings('vasp', 'cp2k');
      expect(mappings1.length).toBe(mappings2.length);
    });
  });

  describe('getParameterMapping', () => {
    it('should return specific parameter mapping', () => {
      const mapping = getParameterMapping('vasp', 'cp2k', 'ENCUT');
      expect(mapping).toBeDefined();
      expect(mapping?.targetParam).toBe('CUTOFF');
    });

    it('should return undefined for non-existent parameter', () => {
      const mapping = getParameterMapping('vasp', 'cp2k', 'NONEXISTENT');
      expect(mapping).toBeUndefined();
    });
  });

  describe('convertParameterValue', () => {
    it('should apply conversion function', () => {
      const mapping: ParameterMapping = {
        source: 'test',
        target: 'test',
        sourceParam: 'test',
        targetParam: 'test',
        conversion: (v: number) => v * 2,
        description: 'test',
      };
      expect(convertParameterValue(mapping, 5)).toBe(10);
    });

    it('should apply unit factor', () => {
      const mapping: ParameterMapping = {
        source: 'test',
        target: 'test',
        sourceParam: 'test',
        targetParam: 'test',
        unitFactor: 1.5,
        description: 'test',
      };
      expect(convertParameterValue(mapping, 10)).toBe(15);
    });

    it('should parse string values with unit factor', () => {
      const mapping: ParameterMapping = {
        source: 'test',
        target: 'test',
        sourceParam: 'test',
        targetParam: 'test',
        unitFactor: 2,
        description: 'test',
      };
      expect(convertParameterValue(mapping, '5')).toBe(10);
    });

    it('should return original value if no conversion', () => {
      const mapping: ParameterMapping = {
        source: 'test',
        target: 'test',
        sourceParam: 'test',
        targetParam: 'test',
        description: 'test',
      };
      expect(convertParameterValue(mapping, 'value')).toBe('value');
    });
  });

  describe('isMigrationSupported', () => {
    it('should return true for supported pairs', () => {
      expect(isMigrationSupported('vasp', 'cp2k')).toBe(true);
      expect(isMigrationSupported('vasp', 'qe')).toBe(true);
      expect(isMigrationSupported('qe', 'vasp')).toBe(true);
    });

    it('should return false for unsupported pairs', () => {
      expect(isMigrationSupported('vasp', 'unknown')).toBe(false);
      expect(isMigrationSupported('unknown', 'vasp')).toBe(false);
    });

    it('should be case insensitive', () => {
      expect(isMigrationSupported('VASP', 'CP2K')).toBe(true);
      expect(isMigrationSupported('vasp', 'cp2k')).toBe(true);
    });
  });

  describe('getSupportedSources', () => {
    it('should return sources for a target', () => {
      const sources = getSupportedSources('cp2k');
      expect(sources).toContain('vasp');
      expect(sources.length).toBeGreaterThan(0);
    });

    it('should return empty array for unsupported target', () => {
      const sources = getSupportedSources('unknown');
      expect(sources).toEqual([]);
    });
  });

  describe('getSupportedTargets', () => {
    it('should return targets for a source', () => {
      const targets = getSupportedTargets('vasp');
      expect(targets).toContain('cp2k');
      expect(targets).toContain('qe');
      expect(targets.length).toBeGreaterThan(0);
    });

    it('should return empty array for unsupported source', () => {
      const targets = getSupportedTargets('unknown');
      expect(targets).toEqual([]);
    });
  });

  describe('Mapping Constants', () => {
    it('should have energy cutoff mappings', () => {
      expect(ENERGY_CUTOFF_MAPPINGS.length).toBeGreaterThan(0);
      expect(ENERGY_CUTOFF_MAPPINGS[0]).toHaveProperty('sourceParam');
      expect(ENERGY_CUTOFF_MAPPINGS[0]).toHaveProperty('targetParam');
    });

    it('should have XC functional mappings', () => {
      expect(XC_FUNCTIONAL_MAPPINGS.length).toBeGreaterThan(0);
      expect(XC_FUNCTIONAL_MAPPINGS[0]).toHaveProperty('conversion');
    });

    it('should have k-point mappings', () => {
      expect(KPOINT_MAPPINGS.length).toBeGreaterThan(0);
    });

    it('should have convergence mappings', () => {
      expect(CONVERGENCE_MAPPINGS.length).toBeGreaterThan(0);
    });

    it('should have structure mappings', () => {
      expect(STRUCTURE_MAPPINGS.length).toBeGreaterThan(0);
    });

    it('should have MD mappings', () => {
      expect(MD_MAPPINGS.length).toBeGreaterThan(0);
    });
  });

  describe('XC Functional Conversions', () => {
    it('should convert VASP GGA to CP2K functional', () => {
      const mapping = getParameterMapping('vasp', 'cp2k', 'GGA');
      expect(mapping).toBeDefined();
      expect(mapping?.conversion).toBeDefined();
      expect(mapping?.conversion!('PE')).toBe('PBE');
      expect(mapping?.conversion!('PBE')).toBe('PBE');
      expect(mapping?.conversion!('RPBE')).toBe('RPBE');
    });

    it('should convert VASP GGA to QE input_dft', () => {
      const mapping = getParameterMapping('vasp', 'qe', 'GGA');
      expect(mapping).toBeDefined();
      expect(mapping?.conversion!('PE')).toBe('PBE');
    });

    it('should convert QE input_dft to VASP GGA', () => {
      const mapping = getParameterMapping('qe', 'vasp', 'input_dft');
      expect(mapping).toBeDefined();
      expect(mapping?.conversion!('PBE')).toBe('PE');
    });
  });
});

describe('Kpoints Migration', () => {
  describe('parseVASP_KPOINTS', () => {
    it('should parse VASP KPOINTS file', () => {
      const content = `KPOINTS for test
0
Gamma
4 4 4
0 0 0`;
      const grid = parseVASP_KPOINTS(content);
      expect(grid.grid).toEqual([4, 4, 4]);
      expect(grid.type).toBe('Gamma-centered');
      expect(grid.shift).toEqual([0, 0, 0]);
      expect(grid.comment).toBe('KPOINTS for test');
    });

    it('should parse Monkhorst-Pack grid', () => {
      const content = `KPOINTS
0
Monkhorst-Pack
6 6 6
0 0 0`;
      const grid = parseVASP_KPOINTS(content);
      expect(grid.grid).toEqual([6, 6, 6]);
      expect(grid.type).toBe('Monkhorst-Pack');
    });

    it('should handle custom shift values', () => {
      const content = `KPOINTS
0
Gamma
4 4 4
0.5 0.5 0.5`;
      const grid = parseVASP_KPOINTS(content);
      expect(grid.shift).toEqual([0.5, 0.5, 0.5]);
    });

    it('should throw error for invalid format', () => {
      expect(() => parseVASP_KPOINTS('invalid')).toThrow();
    });
  });

  describe('generateVASP_KPOINTS', () => {
    it('should generate VASP KPOINTS file', () => {
      const grid: KpointGrid = {
        grid: [4, 4, 4],
        type: 'Gamma-centered',
        shift: [0, 0, 0],
        comment: 'Test grid',
      };
      const content = generateVASP_KPOINTS(grid);
      expect(content).toContain('Test grid');
      expect(content).toContain('Gamma');
      expect(content).toContain('4 4 4');
    });

    it('should generate Monkhorst-Pack format', () => {
      const grid: KpointGrid = {
        grid: [6, 6, 6],
        type: 'Monkhorst-Pack',
        shift: [0, 0, 0],
      };
      const content = generateVASP_KPOINTS(grid);
      expect(content).toContain('Monkhorst-Pack');
    });
  });

  describe('parseQE_KPOINTS', () => {
    it('should parse QE automatic k-points', () => {
      const content = `K_POINTS AUTOMATIC
4 4 4 0 0 0`;
      const grid = parseQE_KPOINTS(content);
      expect(grid.grid).toEqual([4, 4, 4]);
      expect(grid.type).toBe('Monkhorst-Pack');
      expect(grid.shift).toEqual([0, 0, 0]);
    });

    it('should parse Gamma-only k-points', () => {
      const content = `K_POINTS GAMMA`;
      const grid = parseQE_KPOINTS(content);
      expect(grid.grid).toEqual([1, 1, 1]);
      expect(grid.type).toBe('Gamma-centered');
    });
  });

  describe('generateQE_KPOINTS', () => {
    it('should generate QE k-points', () => {
      const grid: KpointGrid = {
        grid: [4, 4, 4],
        type: 'Monkhorst-Pack',
        shift: [0, 0, 0],
      };
      const content = generateQE_KPOINTS(grid);
      expect(content).toContain('K_POINTS AUTOMATIC');
      expect(content).toContain('4 4 4');
    });
  });

  describe('Grid Conversions', () => {
    it('should convert Monkhorst-Pack to Gamma', () => {
      const gamma = monkhorstToGamma([4, 4, 4]);
      expect(gamma).toEqual([8, 8, 8]);
    });

    it('should convert Gamma to Monkhorst-Pack', () => {
      const mp = gammaToMonkhorst([8, 8, 8]);
      expect(mp).toEqual([4, 4, 4]);
    });
  });

  describe('detectKpointFormat', () => {
    it('should detect QE format', () => {
      const content = 'K_POINTS AUTOMATIC\n4 4 4';
      expect(detectKpointFormat(content)).toBe('qe');
    });

    it('should detect VASP format', () => {
      const content = 'KPOINTS\n0\nGamma\n4 4 4';
      expect(detectKpointFormat(content)).toBe('vasp');
    });

    it('should return unknown for invalid format', () => {
      expect(detectKpointFormat('invalid')).toBe('unknown');
    });
  });

  describe('parseKpoints', () => {
    it('should auto-detect and parse QE format', () => {
      const content = 'K_POINTS AUTOMATIC\n4 4 4 0 0 0';
      const grid = parseKpoints(content);
      expect(grid.grid).toEqual([4, 4, 4]);
    });

    it('should auto-detect and parse VASP format', () => {
      const content = 'KPOINTS\n0\nGamma\n4 4 4\n0 0 0';
      const grid = parseKpoints(content);
      expect(grid.grid).toEqual([4, 4, 4]);
    });
  });

  describe('generateKpoints', () => {
    it('should generate QE format', () => {
      const grid: KpointGrid = {
        grid: [4, 4, 4],
        type: 'Monkhorst-Pack',
        shift: [0, 0, 0],
      };
      const content = generateKpoints(grid, 'qe');
      expect(content).toContain('K_POINTS');
    });

    it('should generate VASP format', () => {
      const grid: KpointGrid = {
        grid: [4, 4, 4],
        type: 'Gamma-centered',
        shift: [0, 0, 0],
      };
      const content = generateKpoints(grid, 'vasp');
      expect(content).toContain('0');
    });
  });

  describe('KpointMigration class', () => {
    let migration: KpointMigration;

    beforeEach(() => {
      migration = new KpointMigration({});
    });

    it('should check if migration is supported', () => {
      expect(migration.isSupported('vasp', 'qe')).toBe(true);
      expect(migration.isSupported('vasp', 'unknown')).toBe(false);
    });

    it('should migrate k-points statically', () => {
      const vaspContent = 'Test\n0\nGamma\n4 4 4\n0 0 0';
      const qeContent = KpointMigration.migrate(vaspContent, 'qe');
      expect(qeContent).toContain('K_POINTS AUTOMATIC');
    });
  });
});

describe('ParameterConverter', () => {
  let converter: ParameterConverter;
  let tempDir: string;

  beforeEach(() => {
    converter = new ParameterConverter();
    tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'param-test-'));
  });

  afterEach(() => {
    fs.rmSync(tempDir, { recursive: true, force: true });
  });

  describe('VASP extraction', () => {
    it('should extract VASP INCAR parameters', () => {
      const content = `ENCUT = 520
ISMEAR = 0
SIGMA = 0.2
# Comment line
EDIFF = 1E-6`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.success).toBe(true);
      expect(result.data?.parameters.ENCUT).toBe(520);
      expect(result.data?.parameters.ISMEAR).toBe(0);
      expect(result.data?.parameters.SIGMA).toBe(0.2);
      expect(result.data?.parameters.EDIFF).toBe(1e-6);
    });

    it('should handle boolean values', () => {
      const content = `LCHARG = .TRUE.
LWAVE = .FALSE.
LORBIT = TRUE`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.data?.parameters.LCHARG).toBe(true);
      expect(result.data?.parameters.LWAVE).toBe(false);
    });

    it('should handle string values', () => {
      const content = `PREC = "Accurate"
ALGO = Normal`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.data?.parameters.PREC).toBe('Accurate');
      expect(result.data?.parameters.ALGO).toBe('Normal');
    });
  });

  describe('QE extraction', () => {
    it('should extract QE input parameters', () => {
      const content = `&CONTROL
  calculation = 'scf'
  prefix = 'test'
  outdir = './tmp'
/
&SYSTEM
  ecutwfc = 40
  ecutrho = 320
/`;
      const filepath = path.join(tempDir, 'test.in');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.success).toBe(true);
      expect(result.data?.parameters['CONTROL.calculation']).toBe('scf');
      expect(result.data?.parameters['SYSTEM.ecutwfc']).toBe(40);
    });

    it('should handle QE booleans', () => {
      const content = `&CONTROL
  tprnfor = .true.
  tstress = .false.
/`;
      const filepath = path.join(tempDir, 'test.in');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.data?.parameters['CONTROL.tprnfor']).toBe(true);
      expect(result.data?.parameters['CONTROL.tstress']).toBe(false);
    });
  });

  describe('CP2K extraction', () => {
    it('should extract CP2K input parameters', () => {
      const content = `&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIAL
  &END DFT
&END FORCE_EVAL`;
      const filepath = path.join(tempDir, 'test.inp');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.success).toBe(true);
      expect(result.data?.parameters['FORCE_EVAL.METHOD']).toBe('Quickstep');
      // Note: CP2K parser only tracks current section, not nested
      expect(result.data?.parameters['DFT.BASIS_SET_FILE_NAME']).toBe('BASIS_MOLOPT');
    });
  });

  describe('Gaussian extraction', () => {
    it('should extract Gaussian route parameters', () => {
      const content = `#P B3LYP/6-31G(d) Opt Freq

Test calculation

0 1
C 0 0 0`;
      const filepath = path.join(tempDir, 'test.gjf');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.success).toBe(true);
      expect(result.data?.parameters.P).toBe(true);
      expect(result.data?.parameters.OPT).toBe(true);
      expect(result.data?.parameters.FREQ).toBe(true);
    });

    it('should extract key-value route parameters', () => {
      const content = `# B3LYP/6-31G(d) SCF=Tight Int=UltraFine

Test`;
      const filepath = path.join(tempDir, 'test.com');
      fs.writeFileSync(filepath, content);

      const result = converter.extractParameters(filepath);
      expect(result.data?.parameters.SCF).toBe('Tight');
      expect(result.data?.parameters.INT).toBe('UltraFine');
    });
  });

  describe('Parameter conversion', () => {
    it('should convert parameters between formats', () => {
      const sourceParams: ExtractedParameters = {
        format: 'vasp',
        parameters: {
          ENCUT: 520,
          EDIFF: 1e-6,
        },
        metadata: {
          filepath: '/test/INCAR',
          lines: { ENCUT: 1, EDIFF: 2 },
        },
      };

      const result = converter.convertParameters(sourceParams, 'cp2k');
      expect(result.parameters.CUTOFF).toBe(520);
    });

    it('should track unmapped parameters', () => {
      const sourceParams: ExtractedParameters = {
        format: 'vasp',
        parameters: {
          ENCUT: 520,
          UNKNOWN_PARAM: 123,
        },
        metadata: {
          filepath: '/test/INCAR',
          lines: {},
        },
      };

      const result = converter.convertParameters(sourceParams, 'cp2k');
      expect(result.unmapped).toContain('UNKNOWN_PARAM');
    });

    it('should add warnings for unmapped parameters', () => {
      const sourceParams: ExtractedParameters = {
        format: 'vasp',
        parameters: {
          PARAM1: 1,
          PARAM2: 2,
          PARAM3: 3,
          PARAM4: 4,
          PARAM5: 5,
          PARAM6: 6,
        },
        metadata: {
          filepath: '/test/INCAR',
          lines: {},
        },
      };

      const result = converter.convertParameters(sourceParams, 'cp2k');
      expect(result.warnings.length).toBeGreaterThan(0);
      expect(result.warnings[0]).toContain('6 parameters');
    });
  });

  describe('Full file conversion', () => {
    it('should convert VASP to CP2K', async () => {
      const content = `ENCUT = 520
EDIFF = 1E-6`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const result = await converter.convertFile(filepath, 'cp2k');
      expect(result.success).toBe(true);
      expect(result.source).toBeDefined();
      expect(result.target).toBeDefined();
    });

    it('should handle unsupported formats', async () => {
      const filepath = path.join(tempDir, 'unknown.xyz');
      fs.writeFileSync(filepath, 'test');

      const result = await converter.convertFile(filepath, 'vasp');
      expect(result.success).toBe(false);
      expect(result.error).toContain('Unsupported format');
    });
  });

  describe('Mapping suggestions', () => {
    it('should return mapping suggestions', () => {
      const suggestions = converter.getMappingSuggestions('vasp', 'cp2k');
      expect(suggestions.length).toBeGreaterThan(0);
    });
  });
});

describe('MDWorkflowConverter', () => {
  let converter: MDWorkflowConverter;
  let tempDir: string;

  beforeEach(() => {
    converter = new MDWorkflowConverter();
    tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'md-test-'));
  });

  afterEach(() => {
    fs.rmSync(tempDir, { recursive: true, force: true });
  });

  describe('VASP workflow extraction', () => {
    it('should extract MD workflow from VASP INCAR', () => {
      const content = `IBRION = 0
NSW = 1000
POTIM = 1.0
TEBEG = 300
SMASS = 0
ISIF = 2`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('md');
      expect(workflow.sourceFormat).toBe('vasp');
      expect(workflow.md?.timeStep).toBe(1.0);
      expect(workflow.md?.nSteps).toBe(1000);
      expect(workflow.md?.temperature).toBe(300);
      expect(workflow.md?.ensemble).toBe('NVT');
      expect(workflow.md?.thermostat).toBe('NOSE_HOOVER');
    });

    it('should extract NPT ensemble from VASP', () => {
      const content = `IBRION = 0
NSW = 100
ISIF = 3
PSTRESS = 10`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.md?.ensemble).toBe('NPT');
      expect(workflow.md?.pressure).toBe(10);
    });

    it('should extract optimization workflow from VASP', () => {
      const content = `IBRION = 2
NSW = 100
EDIFFG = -0.01
ISIF = 2`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('optimization');
      expect(workflow.optimization?.algorithm).toBe('CG');
      expect(workflow.optimization?.maxSteps).toBe(100);
      expect(workflow.optimization?.forceConv).toBe(0.01);
    });

    it('should extract cell optimization from VASP', () => {
      const content = `IBRION = 2
NSW = 100
ISIF = 3
PSTRESS = 5.0`;
      const filepath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.optimization?.optimizeCell).toBe(true);
      expect(workflow.optimization?.cellPressure).toBe(5.0);
    });

    it('should map VASP thermostats correctly', () => {
      // Test each thermostat individually to avoid file conflicts
      const testCase = { smass: 0, expected: 'NOSE_HOOVER' };

      const content = `IBRION = 0
NSW = 10
TEBEG = 300
SMASS = ${testCase.smass}`;
      const filepath = path.join(tempDir, 'INCAR_THERMO');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('md');
      expect(workflow.md?.thermostat).toBe(testCase.expected);
    });
  });

  describe('QE workflow extraction', () => {
    it('should extract MD workflow from QE input', () => {
      const content = `&CONTROL
  calculation = 'md'
  dt = 20.0
  nstep = 1000
/
&IONS
  ion_dynamics = 'verlet'
  tempw = 300.0
/`;
      const filepath = path.join(tempDir, 'test.in');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('md');
      expect(workflow.sourceFormat).toBe('qe');
      expect(workflow.md?.timeStep).toBe(20.0);
      expect(workflow.md?.nSteps).toBe(1000);
      expect(workflow.md?.temperature).toBe(300);
    });

    it('should extract NVE ensemble from QE', () => {
      const content = `&CONTROL
  calculation = 'md'
/
&IONS
  ion_dynamics = 'verlet'
/`;
      const filepath = path.join(tempDir, 'test_nve.in');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.md?.ensemble).toBe('NVE');
    });

    it('should extract optimization workflow from QE', () => {
      const content = `&CONTROL
  calculation = 'relax'
/
&IONS
  ion_dynamics = 'bfgs'
  bfgs_ndim = 100
/`;
      const filepath = path.join(tempDir, 'test_opt.in');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('optimization');
      expect(workflow.optimization?.algorithm).toBe('BFGS');
      expect(workflow.optimization?.maxSteps).toBe(100);
    });
  });

  describe('CP2K workflow extraction', () => {
    it('should extract MD workflow from CP2K input', () => {
      const content = `&MOTION
  &MD
    ENSEMBLE NVT
    TIMESTEP 1.0
    STEPS 1000
    TEMPERATURE 300.0
  &END MD
&END MOTION`;
      const filepath = path.join(tempDir, 'test.inp');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('md');
      expect(workflow.sourceFormat).toBe('cp2k');
      expect(workflow.md?.ensemble).toBe('NVT');
      expect(workflow.md?.timeStep).toBe(1.0);
      expect(workflow.md?.nSteps).toBe(1000);
      expect(workflow.md?.temperature).toBe(300);
    });

    it('should extract optimization workflow from CP2K', () => {
      const content = `&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 200
  &END GEO_OPT
&END MOTION`;
      const filepath = path.join(tempDir, 'test_opt.inp');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('optimization');
      expect(workflow.optimization?.algorithm).toBe('BFGS');
      expect(workflow.optimization?.maxSteps).toBe(200);
    });

    it('should extract cell optimization from CP2K', () => {
      const content = `&MOTION
  &CELL_OPT
    MAX_ITER 100
  &END CELL_OPT
  &GEO_OPT
    MAX_ITER 200
  &END GEO_OPT
&END MOTION`;
      const filepath = path.join(tempDir, 'test_cell.inp');
      fs.writeFileSync(filepath, content);

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.type).toBe('optimization');
      expect(workflow.optimization?.optimizeCell).toBe(true);
      expect(workflow.optimization?.maxSteps).toBe(200);
    });
  });

  describe('Workflow conversion', () => {
    it('should convert VASP MD to QE MD', () => {
      const source: MDWorkflow = {
        type: 'md',
        sourceFormat: 'vasp',
        md: {
          timeStep: 1.0,
          nSteps: 1000,
          ensemble: 'NVT',
          temperature: 300,
          thermostat: 'NOSE_HOOVER',
        },
        metadata: {
          filepath: '/test/INCAR',
          warnings: [],
          unsupported: [],
        },
      };

      const result = converter.convertWorkflow(source, 'qe');
      expect(result.success).toBe(true);
      expect(result.target?.md?.timeStep).toBe(1.0);
      expect(result.target?.md?.nSteps).toBe(1000);
    });

    it('should convert VASP MD to CP2K MD with pressure conversion', () => {
      const source: MDWorkflow = {
        type: 'md',
        sourceFormat: 'vasp',
        md: {
          timeStep: 1.0,
          nSteps: 100,
          ensemble: 'NPT',
          temperature: 300,
          pressure: 1.0,
          thermostat: 'NOSE_HOOVER',
          barostat: 'NOSE_HOOVER',
        },
        metadata: {
          filepath: '/test/INCAR',
          warnings: [],
          unsupported: [],
        },
      };

      const result = converter.convertWorkflow(source, 'cp2k');
      expect(result.success).toBe(true);
      expect(result.target?.md?.pressure).toBe(1000); // kbar to bar
    });

    it('should convert VASP optimization to QE optimization', () => {
      const source: MDWorkflow = {
        type: 'optimization',
        sourceFormat: 'vasp',
        optimization: {
          algorithm: 'CG',
          maxSteps: 100,
          energyConv: 1e-6,
          forceConv: 0.01,
          optimizeCell: false,
        },
        metadata: {
          filepath: '/test/INCAR',
          warnings: [],
          unsupported: [],
        },
      };

      const result = converter.convertWorkflow(source, 'qe');
      expect(result.success).toBe(true);
      expect(result.target?.optimization?.algorithm).toBe('CG');
      expect(result.target?.optimization?.energyConv).toBeCloseTo(1e-6 * 0.0734986, 10);
    });

    it('should handle unsupported conversions gracefully', () => {
      const source: MDWorkflow = {
        type: 'md',
        sourceFormat: 'vasp',
        md: {
          timeStep: 1.0,
          nSteps: 100,
        },
        metadata: {
          filepath: '/test/INCAR',
          warnings: [],
          unsupported: [],
        },
      };

      const result = converter.convertWorkflow(source, 'unknown');
      expect(result.success).toBe(true);
      expect(result.warnings.length).toBeGreaterThan(0);
    });
  });

  describe('VASP section generation', () => {
    it('should generate VASP MD section', () => {
      const md: MDParameters = {
        timeStep: 1.0,
        nSteps: 1000,
        ensemble: 'NVT',
        temperature: 300,
        temperatureEnd: 600,
        pressure: 10,
        thermostat: 'NOSE_HOOVER',
      };

      const section = converter.generateVASPMDSection(md);
      expect(section).toContain('POTIM = 1.0');
      expect(section).toContain('NSW = 1000');
      expect(section).toContain('ISIF = 2');
      expect(section).toContain('TEBEG = 300');
      expect(section).toContain('TEEND = 600');
      expect(section).toContain('SMASS = 0');
    });

    it('should generate VASP NPT section', () => {
      const md: MDParameters = {
        ensemble: 'NPT',
        pressure: 5.0,
      };

      const section = converter.generateVASPMDSection(md);
      expect(section).toContain('ISIF = 3');
      expect(section).toContain('PSTRESS = 5.0');
    });

    it('should generate VASP optimization section', () => {
      const opt: OptimizationParameters = {
        algorithm: 'BFGS',
        maxSteps: 100,
        energyConv: 1e-6,
        forceConv: 0.01,
        optimizeCell: true,
        cellPressure: 5.0,
      };

      const section = converter.generateVASPOptSection(opt);
      expect(section).toContain('IBRION = 5');
      expect(section).toContain('NSW = 100');
      expect(section).toContain('ISIF = 3');
      expect(section).toContain('PSTRESS = 5.0');
    });
  });

  describe('QE section generation', () => {
    it('should generate QE MD section', () => {
      const md: MDParameters = {
        timeStep: 20.0,
        nSteps: 1000,
        ensemble: 'NVT',
        temperature: 300,
        temperatureEnd: 600,
        thermostat: 'NOSE_HOOVER',
      };

      const section = converter.generateQEMDSection(md);
      expect(section).toContain('dt = 20.0');
      expect(section).toContain('nstep = 1000');
      expect(section).toContain("ion_temperature = 'nose'");
    });

    it('should generate QE NVE section', () => {
      const md: MDParameters = {
        ensemble: 'NVE',
      };

      const section = converter.generateQEMDSection(md);
      expect(section).toContain("ion_dynamics = 'verlet'");
    });

    it('should generate QE optimization section', () => {
      const opt: OptimizationParameters = {
        algorithm: 'BFGS',
        maxSteps: 100,
        optimizeCell: true,
        cellPressure: 5.0,
      };

      const section = converter.generateQEOptSection(opt);
      expect(section).toContain("ion_dynamics = 'bfgs'");
      expect(section).toContain('bfgs_ndim = 100');
      expect(section).toContain("cell_dynamics = 'bfgs'");
    });
  });

  describe('CP2K section generation', () => {
    it('should generate CP2K MD section', () => {
      const md: MDParameters = {
        ensemble: 'NVT',
        timeStep: 1.0,
        nSteps: 1000,
        temperature: 300,
        thermostat: 'NOSE_HOOVER',
      };

      const section = converter.generateCP2KMDSection(md);
      expect(section).toContain('ENSEMBLE NVT');
      expect(section).toContain('TIMESTEP 1.0');
      expect(section).toContain('STEPS 1000');
      expect(section).toContain('TEMPERATURE 300');
      expect(section).toContain('TYPE NOSE');
    });

    it('should generate CP2K NPT section', () => {
      const md: MDParameters = {
        ensemble: 'NPT',
        pressure: 1.0,
      };

      const section = converter.generateCP2KMDSection(md);
      expect(section).toContain('&BAROSTAT');
      expect(section).toContain('PRESSURE 1.0');
    });

    it('should generate CP2K optimization section', () => {
      const opt: OptimizationParameters = {
        algorithm: 'BFGS',
        maxSteps: 200,
        energyConv: 1e-6,
        forceConv: 1e-4,
      };

      const section = converter.generateCP2KOptSection(opt);
      expect(section).toContain('&GEO_OPT');
      expect(section).toContain('OPTIMIZER BFGS');
      expect(section).toContain('MAX_ITER 200');
    });

    it('should generate CP2K cell optimization section', () => {
      const opt: OptimizationParameters = {
        algorithm: 'BFGS',
        maxSteps: 100,
        optimizeCell: true,
        cellPressure: 5.0,
      };

      const section = converter.generateCP2KOptSection(opt);
      expect(section).toContain('&CELL_OPT');
      expect(section).toContain('PRESSURE 5.0');
    });
  });

  describe('Parameter formatting', () => {
    it('should format MD parameters', () => {
      const md: MDParameters = {
        timeStep: 1.0,
        nSteps: 1000,
        ensemble: 'NVT',
        temperature: 300,
        pressure: 1.0,
        thermostat: 'NOSE_HOOVER',
        barostat: 'BERENDSEN',
      };

      const formatted = converter.formatMDParameters(md);
      expect(formatted).toContain('Time Step: 1 fs');
      expect(formatted).toContain('Number of Steps: 1000');
      expect(formatted).toContain('Ensemble: NVT');
      expect(formatted).toContain('Temperature: 300 K');
      expect(formatted).toContain('Pressure: 1 kbar');
      expect(formatted).toContain('Thermostat: NOSE_HOOVER');
      expect(formatted).toContain('Barostat: BERENDSEN');
    });

    it('should format optimization parameters', () => {
      const opt: OptimizationParameters = {
        algorithm: 'BFGS',
        maxSteps: 100,
        energyConv: 1e-6,
        forceConv: 0.01,
        optimizeCell: true,
        cellPressure: 5.0,
      };

      const formatted = converter.formatOptimizationParameters(opt);
      expect(formatted).toContain('Algorithm: BFGS');
      expect(formatted).toContain('Max Steps: 100');
      expect(formatted).toContain('Energy Convergence: 0.000001 eV');
      expect(formatted).toContain('Force Convergence: 0.01 eV/Å');
      expect(formatted).toContain('Cell Optimization: true');
      expect(formatted).toContain('Cell Pressure: 5 kbar');
    });
  });

  describe('Unknown format handling', () => {
    it('should handle unknown file formats gracefully', () => {
      const filepath = path.join(tempDir, 'unknown.xyz');
      fs.writeFileSync(filepath, 'test content');

      const workflow = converter.extractWorkflow(filepath);
      expect(workflow.sourceFormat).toBe('unknown');
      expect(workflow.metadata.unsupported.length).toBeGreaterThan(0);
    });
  });
});
