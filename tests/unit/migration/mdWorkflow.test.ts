/**
 * MD Workflow Migration Tests
 */

import * as path from 'path';
import * as fs from 'fs';
import * as os from 'os';
import {
  MDWorkflowConverter,
  EnsembleType,
  ThermostatType,
  OptimizationAlgorithm,
} from '../../../src/utils/migration/mdWorkflow';

describe('MDWorkflowConverter', () => {
  let converter: MDWorkflowConverter;
  let tempDir: string;

  beforeEach(() => {
    converter = new MDWorkflowConverter();
    tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-md-test-'));
  });

  afterEach(() => {
    if (fs.existsSync(tempDir)) {
      fs.rmSync(tempDir, { recursive: true, force: true });
    }
  });

  describe('VASP MD Workflow Extraction', () => {
    it('should extract MD parameters from VASP INCAR', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `# VASP MD INCAR
ENCUT = 520
POTIM = 1.0
NSW = 1000
ISIF = 2
TEBEG = 300
TEEND = 400
PSTRESS = 0.0
SMASS = 0
`
      );

      const workflow = converter.extractWorkflow(incarPath);

      expect(workflow.type).toBe('md');
      expect(workflow.sourceFormat).toBe('vasp');
      expect(workflow.md).toBeDefined();
      expect(workflow.md?.timeStep).toBe(1.0);
      expect(workflow.md?.nSteps).toBe(1000);
      expect(workflow.md?.ensemble).toBe('NVT');
      expect(workflow.md?.temperature).toBe(300);
      expect(workflow.md?.temperatureEnd).toBe(400);
      expect(workflow.md?.pressure).toBe(0.0);
      expect(workflow.md?.thermostat).toBe('NOSE_HOOVER');
    });

    it('should detect NPT ensemble from ISIF=3', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `POTIM = 1.0
NSW = 1000
ISIF = 3
TEBEG = 300
PSTRESS = 1.0
`
      );

      const workflow = converter.extractWorkflow(incarPath);

      expect(workflow.type).toBe('md');
      expect(workflow.md?.ensemble).toBe('NPT');
      expect(workflow.md?.pressure).toBe(1.0);
    });

    it('should detect Berendsen thermostat from SMASS=2', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `POTIM = 1.0
NSW = 1000
IBRION = 0
SMASS = 2
`
      );

      const workflow = converter.extractWorkflow(incarPath);

      expect(workflow.md?.thermostat).toBe('BERENDSEN');
    });

    it('should extract optimization parameters from VASP INCAR', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `# VASP Optimization
ENCUT = 520
IBRION = 2
NSW = 100
EDIFF = 1E-5
EDIFFG = -0.02
ISIF = 2
`
      );

      const workflow = converter.extractWorkflow(incarPath);

      expect(workflow.type).toBe('optimization');
      expect(workflow.sourceFormat).toBe('vasp');
      expect(workflow.optimization).toBeDefined();
      expect(workflow.optimization?.algorithm).toBe('CG');
      expect(workflow.optimization?.maxSteps).toBe(100);
      expect(workflow.optimization?.energyConv).toBe(1e-5);
      expect(workflow.optimization?.forceConv).toBe(0.02);
    });

    it('should detect BFGS algorithm from IBRION=5', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `IBRION = 5
NSW = 100
`
      );

      const workflow = converter.extractWorkflow(incarPath);

      expect(workflow.optimization?.algorithm).toBe('BFGS');
    });
  });

  describe('QE MD Workflow Extraction', () => {
    it('should extract MD parameters from QE input', () => {
      const qePath = path.join(tempDir, 'md.in');
      fs.writeFileSync(
        qePath,
        `&CONTROL
  dt = 1.0
  nstep = 1000
  calculation = 'md'
/
&IONS
  ion_dynamics = 'verlet'
  tempw = 300
/
&CELL
  press = 0.0
/
`
      );

      const workflow = converter.extractWorkflow(qePath);

      expect(workflow.type).toBe('md');
      expect(workflow.sourceFormat).toBe('qe');
      expect(workflow.md?.timeStep).toBe(1.0);
      expect(workflow.md?.nSteps).toBe(1000);
      expect(workflow.md?.ensemble).toBe('NVE');
      expect(workflow.md?.temperature).toBe(300);
    });

    it('should detect NVT ensemble with Nose thermostat', () => {
      const qePath = path.join(tempDir, 'md.in');
      fs.writeFileSync(
        qePath,
        `&CONTROL
  dt = 1.0
  nstep = 1000
/
&IONS
  ion_dynamics = 'verlet'
  ion_temperature = 'nose'
  tempw = 300
/
`
      );

      const workflow = converter.extractWorkflow(qePath);

      expect(workflow.md?.ensemble).toBe('NVT');
      expect(workflow.md?.thermostat).toBe('NOSE_HOOVER');
    });

    it('should detect Langevin thermostat', () => {
      const qePath = path.join(tempDir, 'md.in');
      fs.writeFileSync(
        qePath,
        `&IONS
  ion_dynamics = 'langevin'
  ion_temperature = 'langevin'
  tempw = 300
/
`
      );

      const workflow = converter.extractWorkflow(qePath);

      expect(workflow.md?.thermostat).toBe('LANGEVIN');
    });

    it('should extract optimization parameters from QE input', () => {
      const qePath = path.join(tempDir, 'relax.in');
      fs.writeFileSync(
        qePath,
        `&CONTROL
  calculation = 'relax'
/
&IONS
  ion_dynamics = 'bfgs'
  bfgs_ndim = 100
/
&ELECTRONS
  conv_thr = 1E-6
/
`
      );

      const workflow = converter.extractWorkflow(qePath);

      expect(workflow.type).toBe('optimization');
      expect(workflow.optimization?.algorithm).toBe('BFGS');
      expect(workflow.optimization?.maxSteps).toBe(100);
      expect(workflow.optimization?.energyConv).toBe(1e-6);
    });

    it('should detect CG algorithm', () => {
      const qePath = path.join(tempDir, 'relax.in');
      fs.writeFileSync(
        qePath,
        `&IONS
  ion_dynamics = 'cg'
/
`
      );

      const workflow = converter.extractWorkflow(qePath);

      expect(workflow.optimization?.algorithm).toBe('CG');
    });
  });

  describe('CP2K MD Workflow Extraction', () => {
    it('should extract MD parameters from CP2K input', () => {
      const cp2kPath = path.join(tempDir, 'md.inp');
      fs.writeFileSync(
        cp2kPath,
        `&MOTION
  &MD
    ENSEMBLE NVT
    TIMESTEP 1.0
    STEPS 1000
    TEMPERATURE 300
    &THERMOSTAT
      TYPE NOSE
    &END THERMOSTAT
  &END MD
&END MOTION
`
      );

      const workflow = converter.extractWorkflow(cp2kPath);

      expect(workflow.type).toBe('md');
      expect(workflow.sourceFormat).toBe('cp2k');
      expect(workflow.md?.timeStep).toBe(1.0);
      expect(workflow.md?.nSteps).toBe(1000);
      expect(workflow.md?.ensemble).toBe('NVT');
      expect(workflow.md?.temperature).toBe(300);
      expect(workflow.md?.thermostat).toBe('NOSE_HOOVER');
    });

    it('should extract NPT ensemble', () => {
      const cp2kPath = path.join(tempDir, 'md.inp');
      fs.writeFileSync(
        cp2kPath,
        `&MOTION
  &MD
    ENSEMBLE NPT
    TIMESTEP 1.0
    STEPS 1000
    PRESSURE 1.0
    &BAROSTAT
      TYPE NOSE
    &END BAROSTAT
  &END MD
&END MOTION
`
      );

      const workflow = converter.extractWorkflow(cp2kPath);

      expect(workflow.md?.ensemble).toBe('NPT');
      expect(workflow.md?.pressure).toBe(1.0);
      expect(workflow.md?.barostat).toBe('NOSE_HOOVER');
    });

    it('should extract optimization parameters from CP2K input', () => {
      const cp2kPath = path.join(tempDir, 'opt.inp');
      fs.writeFileSync(
        cp2kPath,
        `&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 100
    &CONVERGENCE
      EPS_ENERGY 1E-6
      EPS_FORCE 1E-4
    &END CONVERGENCE
  &END GEO_OPT
&END MOTION
`
      );

      const workflow = converter.extractWorkflow(cp2kPath);

      expect(workflow.type).toBe('optimization');
      expect(workflow.optimization?.algorithm).toBe('BFGS');
      expect(workflow.optimization?.maxSteps).toBe(100);
      expect(workflow.optimization?.energyConv).toBe(1e-6);
      expect(workflow.optimization?.forceConv).toBe(1e-4);
    });
  });

  describe('MD Workflow Conversion', () => {
    it('should convert VASP MD to QE MD', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `POTIM = 1.0
NSW = 1000
ISIF = 2
TEBEG = 300
SMASS = 0
`
      );

      const workflow = converter.extractWorkflow(incarPath);
      const result = converter.convertWorkflow(workflow, 'qe');

      expect(result.success).toBe(true);
      expect(result.target?.type).toBe('md');
      expect(result.target?.md?.timeStep).toBe(1.0);
      expect(result.target?.md?.nSteps).toBe(1000);
      expect(result.target?.md?.ensemble).toBe('NVT');
      expect(result.target?.md?.temperature).toBe(300);
    });

    it('should convert VASP MD to CP2K MD with unit conversion', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `POTIM = 1.0
NSW = 1000
ISIF = 3
TEBEG = 300
PSTRESS = 1.0  # kbar
`
      );

      const workflow = converter.extractWorkflow(incarPath);
      const result = converter.convertWorkflow(workflow, 'cp2k');

      expect(result.success).toBe(true);
      expect(result.target?.md?.pressure).toBe(1000); // kbar to bar
    });

    it('should convert QE MD to VASP MD', () => {
      const qePath = path.join(tempDir, 'md.in');
      fs.writeFileSync(
        qePath,
        `&CONTROL
  dt = 1.0
  nstep = 1000
/
&IONS
  ion_dynamics = 'verlet'
  ion_temperature = 'nose'
  tempw = 300
/
`
      );

      const workflow = converter.extractWorkflow(qePath);
      const result = converter.convertWorkflow(workflow, 'vasp');

      expect(result.success).toBe(true);
      expect(result.target?.md?.timeStep).toBe(1.0);
      expect(result.target?.md?.nSteps).toBe(1000);
      expect(result.target?.md?.temperature).toBe(300);
      expect(result.target?.md?.thermostat).toBe('NOSE_HOOVER');
    });

    it('should convert VASP optimization to QE optimization', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `IBRION = 2
NSW = 100
EDIFF = 1E-5
`
      );

      const workflow = converter.extractWorkflow(incarPath);
      const result = converter.convertWorkflow(workflow, 'qe');

      expect(result.success).toBe(true);
      expect(result.target?.type).toBe('optimization');
      expect(result.target?.optimization?.algorithm).toBe('CG');
      expect(result.target?.optimization?.maxSteps).toBe(100);
      // eV to Ry conversion
      expect(result.target?.optimization?.energyConv).toBeCloseTo(1e-5 * 0.0734986, 10);
    });

    it('should convert VASP optimization to CP2K optimization with unit conversion', () => {
      const incarPath = path.join(tempDir, 'INCAR');
      fs.writeFileSync(
        incarPath,
        `IBRION = 5
NSW = 100
EDIFF = 1E-5
EDIFFG = -0.02
ISIF = 3
PSTRESS = 1.0
`
      );

      const workflow = converter.extractWorkflow(incarPath);
      const result = converter.convertWorkflow(workflow, 'cp2k');

      expect(result.success).toBe(true);
      expect(result.target?.optimization?.algorithm).toBe('BFGS');
      // eV to hartree
      expect(result.target?.optimization?.energyConv).toBeCloseTo(1e-5 * 0.0367493, 10);
      // eV/Å to hartree/bohr
      expect(result.target?.optimization?.forceConv).toBeCloseTo(0.02 * 0.01943, 10);
      // kbar to bar
      expect(result.target?.optimization?.cellPressure).toBe(1000);
    });
  });

  describe('VASP Output Generation', () => {
    it('should generate VASP MD section', () => {
      const mdParams = {
        timeStep: 1.0,
        nSteps: 1000,
        ensemble: 'NVT' as EnsembleType,
        temperature: 300,
        thermostat: 'NOSE_HOOVER' as ThermostatType,
      };

      const output = converter.generateVASPMDSection(mdParams);

      expect(output).toContain('POTIM = 1.0');
      expect(output).toContain('NSW = 1000');
      expect(output).toContain('IBRION = 0');
      expect(output).toContain('ISIF = 2');
      expect(output).toContain('TEBEG = 300');
      expect(output).toContain('SMASS = 0');
    });

    it('should generate VASP NPT MD section', () => {
      const mdParams = {
        timeStep: 1.0,
        ensemble: 'NPT' as EnsembleType,
        temperature: 300,
        pressure: 1.0,
        thermostat: 'NOSE_HOOVER' as ThermostatType,
      };

      const output = converter.generateVASPMDSection(mdParams);

      expect(output).toContain('ISIF = 3');
      expect(output).toContain('PSTRESS = 1.0');
    });

    it('should generate VASP optimization section', () => {
      const optParams = {
        algorithm: 'BFGS' as OptimizationAlgorithm,
        maxSteps: 100,
        energyConv: 1e-5,
        forceConv: 0.02,
        optimizeCell: false,
      };

      const output = converter.generateVASPOptSection(optParams);

      expect(output).toContain('IBRION = 5');
      expect(output).toContain('NSW = 100');
      expect(output).toContain('EDIFF = 1e-05');
      expect(output).toContain('EDIFFG = -0.02');
      expect(output).toContain('ISIF = 2');
    });

    it('should generate VASP variable-cell optimization', () => {
      const optParams = {
        algorithm: 'CG' as OptimizationAlgorithm,
        maxSteps: 100,
        energyConv: 1e-5,
        optimizeCell: true,
        cellPressure: 1.0,
      };

      const output = converter.generateVASPOptSection(optParams);

      expect(output).toContain('ISIF = 3');
      expect(output).toContain('PSTRESS = 1.0');
    });
  });

  describe('QE Output Generation', () => {
    it('should generate QE MD section', () => {
      const mdParams = {
        timeStep: 1.0,
        nSteps: 1000,
        ensemble: 'NVT' as EnsembleType,
        temperature: 300,
        thermostat: 'NOSE_HOOVER' as ThermostatType,
      };

      const output = converter.generateQEMDSection(mdParams);

      expect(output).toContain('&CONTROL');
      expect(output).toContain('dt = 1.0');
      expect(output).toContain('nstep = 1000');
      expect(output).toContain('&IONS');
      expect(output).toContain("ion_dynamics = 'verlet'");
      expect(output).toContain('tempw = 300');
      expect(output).toContain("ion_temperature = 'nose'");
    });

    it('should generate QE NPT MD section', () => {
      const mdParams = {
        timeStep: 1.0,
        ensemble: 'NPT' as EnsembleType,
        temperature: 300,
        pressure: 1.0,
      };

      const output = converter.generateQEMDSection(mdParams);

      expect(output).toContain('&CELL');
      expect(output).toContain('press = 1.0');
      expect(output).toContain("cell_dynamics = 'pr'");
    });

    it('should generate QE optimization section', () => {
      const optParams = {
        algorithm: 'BFGS' as OptimizationAlgorithm,
        maxSteps: 100,
        energyConv: 1e-6,
        optimizeCell: false,
      };

      const output = converter.generateQEOptSection(optParams);

      expect(output).toContain('&IONS');
      expect(output).toContain("ion_dynamics = 'bfgs'");
      expect(output).toContain('bfgs_ndim = 100');
    });

    it('should generate QE variable-cell optimization', () => {
      const optParams = {
        algorithm: 'BFGS' as OptimizationAlgorithm,
        optimizeCell: true,
        cellPressure: 1.0,
      };

      const output = converter.generateQEOptSection(optParams);

      expect(output).toContain('&CELL');
      expect(output).toContain("cell_dynamics = 'bfgs'");
      expect(output).toContain('press = 1.0');
    });
  });

  describe('CP2K Output Generation', () => {
    it('should generate CP2K MD section', () => {
      const mdParams = {
        timeStep: 1.0,
        nSteps: 1000,
        ensemble: 'NVT' as EnsembleType,
        temperature: 300,
        thermostat: 'NOSE_HOOVER' as ThermostatType,
      };

      const output = converter.generateCP2KMDSection(mdParams);

      expect(output).toContain('&MOTION');
      expect(output).toContain('&MD');
      expect(output).toContain('ENSEMBLE NVT');
      expect(output).toContain('TIMESTEP 1.0');
      expect(output).toContain('STEPS 1000');
      expect(output).toContain('TEMPERATURE 300');
      expect(output).toContain('&THERMOSTAT');
      expect(output).toContain('TYPE NOSE');
    });

    it('should generate CP2K NPT MD section', () => {
      const mdParams = {
        timeStep: 1.0,
        ensemble: 'NPT' as EnsembleType,
        temperature: 300,
        pressure: 1.0,
        thermostat: 'NOSE_HOOVER' as ThermostatType,
      };

      const output = converter.generateCP2KMDSection(mdParams);

      expect(output).toContain('&BAROSTAT');
      expect(output).toContain('TYPE NOSE');
      expect(output).toContain('PRESSURE 1.0');
    });

    it('should generate CP2K optimization section', () => {
      const optParams = {
        algorithm: 'BFGS' as OptimizationAlgorithm,
        maxSteps: 100,
        energyConv: 1e-6,
        forceConv: 1e-4,
        optimizeCell: false,
      };

      const output = converter.generateCP2KOptSection(optParams);

      expect(output).toContain('&MOTION');
      expect(output).toContain('&GEO_OPT');
      expect(output).toContain('OPTIMIZER BFGS');
      expect(output).toContain('MAX_ITER 100');
      expect(output).toContain('&CONVERGENCE');
      expect(output).toContain('EPS_ENERGY 1e-06');
      expect(output).toContain('EPS_FORCE 1e-04');
    });

    it('should generate CP2K cell optimization', () => {
      const optParams = {
        algorithm: 'BFGS' as OptimizationAlgorithm,
        maxSteps: 100,
        optimizeCell: true,
        cellPressure: 1.0,
      };

      const output = converter.generateCP2KOptSection(optParams);

      expect(output).toContain('&CELL_OPT');
      expect(output).toContain('PRESSURE 1.0');
    });
  });

  describe('Parameter Formatting', () => {
    it('should format MD parameters for display', () => {
      const mdParams = {
        timeStep: 1.0,
        nSteps: 1000,
        ensemble: 'NVT' as EnsembleType,
        temperature: 300,
        temperatureEnd: 400,
        pressure: 1.0,
        thermostat: 'NOSE_HOOVER' as ThermostatType,
      };

      const output = converter.formatMDParameters(mdParams);

      expect(output).toContain('Time Step: 1 fs');
      expect(output).toContain('Number of Steps: 1000');
      expect(output).toContain('Ensemble: NVT');
      expect(output).toContain('Temperature: 300 K');
      expect(output).toContain('Final Temperature: 400 K');
      expect(output).toContain('Pressure: 1 kbar');
      expect(output).toContain('Thermostat: NOSE_HOOVER');
    });

    it('should format optimization parameters for display', () => {
      const optParams = {
        algorithm: 'BFGS' as OptimizationAlgorithm,
        maxSteps: 100,
        energyConv: 1e-5,
        forceConv: 0.02,
        optimizeCell: false,
      };

      const output = converter.formatOptimizationParameters(optParams);

      expect(output).toContain('Algorithm: BFGS');
      expect(output).toContain('Max Steps: 100');
      expect(output).toContain('Energy Convergence: 0.00001 eV');
      expect(output).toContain('Force Convergence: 0.02 eV/Å');
      expect(output).toContain('Cell Optimization: false');
    });
  });

  describe('Error Handling', () => {
    it('should handle unsupported format', () => {
      const filePath = path.join(tempDir, 'unknown.xyz');
      fs.writeFileSync(filePath, 'C 0 0 0');

      const workflow = converter.extractWorkflow(filePath);

      expect(workflow.sourceFormat).toBe('unknown');
      expect(workflow.metadata.unsupported.length).toBeGreaterThan(0);
    });

    it('should handle missing file', () => {
      const workflow = converter.extractWorkflow('/nonexistent/file');

      expect(workflow.sourceFormat).toBe('unknown');
      expect(workflow.metadata.unsupported.length).toBeGreaterThan(0);
    });
  });
});
