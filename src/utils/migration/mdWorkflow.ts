/**
 * MD/Optimization Workflow Migration Module
 *
 * Provides comprehensive migration of MD and optimization parameters
 * across quantum chemistry codes (VASP, CP2K, QE, LAMMPS)
 *
 * @module mdWorkflow
 * @see Issue #12 - Phase 3.4
 */

import * as vscode from 'vscode';
import * as fs from 'fs';
import * as path from 'path';

/**
 * MD ensemble types
 */
export type EnsembleType = 'NVE' | 'NVT' | 'NPT' | 'NPH';

/**
 * Thermostat types
 */
export type ThermostatType = 'NOSE_HOOVER' | 'BERENDSEN' | 'LANGEVIN' | 'VELOCITY_SCALING' | 'ANDERSEN';

/**
 * Barostat types
 */
export type BarostatType = 'NOSE_HOOVER' | 'BERENDSEN' | 'PARRINELLO_RAHMAN' | 'MTTK';

/**
 * Optimization algorithm types
 */
export type OptimizationAlgorithm = 'CG' | 'BFGS' | 'LBFGS' | 'FIRE' | 'MD' | 'RMMDIIS' | 'DAMPED_MD';

/**
 * MD Parameters for a specific code
 */
export interface MDParameters {
  /** Time step (fs) */
  timeStep?: number;
  /** Number of MD steps */
  nSteps?: number;
  /** Ensemble type */
  ensemble?: EnsembleType;
  /** Initial temperature (K) */
  temperature?: number;
  /** Final temperature (K) - for temperature ramping */
  temperatureEnd?: number;
  /** Target pressure (kbar or GPa depending on code) */
  pressure?: number;
  /** Thermostat type */
  thermostat?: ThermostatType;
  /** Thermostat parameters */
  thermostatParams?: Record<string, any>;
  /** Barostat type */
  barostat?: BarostatType;
  /** Barostat parameters */
  barostatParams?: Record<string, any>;
  /** Temperature damping (fs) - for Berendsen */
  temperatureDamping?: number;
  /** Pressure damping (fs) - for Berendsen */
  pressureDamping?: number;
  /** Chain length for Nose-Hoover */
  chainLength?: number;
  /** Print frequency for MD output */
  printFreq?: number;
  /** Trajectory save frequency */
  trajFreq?: number;
  /** Random seed for stochastic dynamics */
  randomSeed?: number;
}

/**
 * Optimization Parameters for a specific code
 */
export interface OptimizationParameters {
  /** Optimization algorithm */
  algorithm?: OptimizationAlgorithm;
  /** Maximum number of ionic steps */
  maxSteps?: number;
  /** Energy convergence criterion (eV) */
  energyConv?: number;
  /** Force convergence criterion (eV/Å) */
  forceConv?: number;
  /** Stress convergence criterion (kbar) */
  stressConv?: number;
  /** Step size for line search */
  stepSize?: number;
  /** Maximum displacement per step (Å) */
  maxDisplacement?: number;
  /** Line minimization precision */
  lineMinPrecision?: number;
  /** Use cell optimization */
  optimizeCell?: boolean;
  /** Cell pressure for variable-cell optimization (kbar) */
  cellPressure?: number;
  /** Constrained atoms */
  constraints?: number[];
}

/**
 * Complete MD/Optimization workflow
 */
export interface MDWorkflow {
  /** Workflow type: 'md' or 'optimization' */
  type: 'md' | 'optimization';
  /** Source format */
  sourceFormat: string;
  /** MD parameters */
  md?: MDParameters;
  /** Optimization parameters */
  optimization?: OptimizationParameters;
  /** Additional metadata */
  metadata: {
    filepath: string;
    warnings: string[];
    unsupported: string[];
  };
}

/**
 * MD Workflow conversion result
 */
export interface MDWorkflowConversionResult {
  success: boolean;
  source?: MDWorkflow;
  target?: MDWorkflow;
  errors: string[];
  warnings: string[];
}

/**
 * MD Workflow Converter class
 */
export class MDWorkflowConverter {
  private outputChannel: vscode.OutputChannel;

  constructor() {
    this.outputChannel = vscode.window.createOutputChannel('OpenQC MD Workflow');
  }

  /**
   * Parse a value string into appropriate type
   */
  private parseValue(value: string): any {
    // Boolean
    if (value.toUpperCase() === '.TRUE.' || value.toUpperCase() === 'TRUE') return true;
    if (value.toUpperCase() === '.FALSE.' || value.toUpperCase() === 'FALSE') return false;

    // Number
    const num = parseFloat(value);
    if (!isNaN(num)) return num;

    // String
    return value.replace(/^["']|["']$/g, '');
  }

  /**
   * Extract MD/optimization workflow from VASP INCAR
   */
  private extractVASPWorkflow(filepath: string): MDWorkflow {
    const content = fs.readFileSync(filepath, 'utf-8');
    const params: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    // Parse INCAR
    content.split('\n').forEach(line => {
      const trimmed = line.trim();
      if (!trimmed || trimmed.startsWith('#')) return;

      const match = trimmed.match(/^(\w+)\s*=\s*(.+?)(?:\s*#.*)?$/);
      if (match) {
        const [, key, value] = match;
        params[key.toUpperCase()] = this.parseValue(value.trim());
      }
    });

    // Detect workflow type
    const ibrion = params.IBRION || 0;
    const nsw = params.NSW || 0;
    const isMD = (ibrion === 0 || ibrion === 1 || ibrion === 3) && nsw > 0;
    const isOptimization = (ibrion >= 1 && ibrion <= 3) || ibrion === 5 || ibrion === 7;

    const workflow: MDWorkflow = {
      type: isMD ? 'md' : 'optimization',
      sourceFormat: 'vasp',
      metadata: {
        filepath,
        warnings,
        unsupported
      }
    };

    if (isMD) {
      workflow.md = this.extractVASPMD(params);
    } else if (isOptimization) {
      workflow.optimization = this.extractVASPOpt(params);
    }

    return workflow;
  }

  /**
   * Extract VASP MD parameters
   */
  private extractVASPMD(params: Record<string, any>): MDParameters {
    const md: MDParameters = {};

    // Time step
    if (params.POTIM) md.timeStep = params.POTIM;
    if (params.NSW) md.nSteps = params.NSW;

    // Ensemble
    if (params.ISIF !== undefined) {
      if (params.ISIF === 2) md.ensemble = 'NVT';
      else if (params.ISIF === 3) md.ensemble = 'NPT';
      else if (params.ISIF === 7) md.ensemble = 'NPH';
      else md.ensemble = 'NVT';
    }

    // Temperature
    if (params.TEBEG) md.temperature = params.TEBEG;
    if (params.TEEND) md.temperatureEnd = params.TEEND;

    // Pressure
    if (params.PSTRESS) md.pressure = params.PSTRESS;

    // Thermostat
    if (params.SMASS !== undefined) {
      if (params.SMASS === 0) md.thermostat = 'NOSE_HOOVER';
      else if (params.SMASS === 1) md.thermostat = 'VELOCITY_SCALING';
      else if (params.SMASS === 2) md.thermostat = 'BERENDSEN';
      else if (params.SMASS === 3) md.thermostat = 'ANDERSEN';
    }

    // Random seed
    if (params.RANDOM_SEED) md.randomSeed = params.RANDOM_SEED;

    // Print frequency
    if (params.NWRITE) md.printFreq = params.NWRITE;
    if (params.NBLOCK) md.trajFreq = params.NBLOCK;

    return md;
  }

  /**
   * Extract VASP optimization parameters
   */
  private extractVASPOpt(params: Record<string, any>): OptimizationParameters {
    const opt: OptimizationParameters = {};

    // Algorithm
    if (params.IBRION !== undefined) {
      if (params.IBRION === 1) opt.algorithm = 'RMMDIIS';
      else if (params.IBRION === 2) opt.algorithm = 'CG';
      else if (params.IBRION === 3) opt.algorithm = 'DAMPED_MD';
      else if (params.IBRION === 5) opt.algorithm = 'BFGS';
      else if (params.IBRION === 7) opt.algorithm = 'MD';
      else opt.algorithm = 'CG';
    }

    // Max steps
    if (params.NSW) opt.maxSteps = params.NSW;

    // Convergence
    if (params.EDIFF) opt.energyConv = params.EDIFF;
    if (params.EDIFFG) opt.forceConv = Math.abs(params.EDIFFG);

    // Step size
    if (params.POTIM) opt.stepSize = params.POTIM;

    // Cell optimization
    if (params.ISIF !== undefined) {
      opt.optimizeCell = params.ISIF >= 3 && params.ISIF <= 7;
      if (params.PSTRESS) opt.cellPressure = params.PSTRESS;
    }

    return opt;
  }

  /**
   * Extract MD/optimization workflow from QE input
   */
  private extractQEWorkflow(filepath: string): MDWorkflow {
    const content = fs.readFileSync(filepath, 'utf-8');
    const params: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    let currentSection = '';
    content.split('\n').forEach(line => {
      const trimmed = line.trim();

      if (trimmed.startsWith('&')) {
        currentSection = trimmed.substring(1).toUpperCase();
        return;
      }
      if (trimmed === '/') {
        currentSection = '';
        return;
      }
      if (!trimmed || trimmed.startsWith('!')) return;

      if (currentSection && trimmed.includes('=')) {
        const match = trimmed.match(/^(\w+)\s*=\s*(.+?)(?:\s*!.*)?$/);
        if (match) {
          const [, key, value] = match;
          params[`${currentSection}.${key.toLowerCase()}`] = this.parseValue(value.trim());
        }
      }
    });

    // Detect workflow type
    const ionDynamics = params['IONS.ion_dynamics'] || '';
    const isMD = ionDynamics.toLowerCase().includes('verlet') ||
                 ionDynamics.toLowerCase().includes('langevin');
    const isOptimization = ionDynamics.toLowerCase().includes('bfgs') ||
                          ionDynamics.toLowerCase().includes('cg');

    const workflow: MDWorkflow = {
      type: isMD ? 'md' : 'optimization',
      sourceFormat: 'qe',
      metadata: { filepath, warnings, unsupported }
    };

    if (isMD) {
      workflow.md = this.extractQEMD(params);
    } else if (isOptimization) {
      workflow.optimization = this.extractQEOpt(params);
    }

    return workflow;
  }

  /**
   * Extract QE MD parameters
   */
  private extractQEMD(params: Record<string, any>): MDParameters {
    const md: MDParameters = {};

    // Time step
    if (params['CONTROL.dt']) md.timeStep = params['CONTROL.dt'];
    if (params['CONTROL.nstep']) md.nSteps = params['CONTROL.nstep'];

    // Ensemble
    if (params['IONS.ion_dynamics']) {
      const dyn = params['IONS.ion_dynamics'].toLowerCase();
      if (dyn === 'verlet') md.ensemble = 'NVE';
      else if (dyn === 'langevin') md.ensemble = 'NVT';
      else md.ensemble = 'NVT';
    }

    // Temperature
    if (params['IONS.tempw']) md.temperature = params['IONS.tempw'];
    if (params['IONS.templ']) md.temperatureEnd = params['IONS.templ'];

    // Pressure
    if (params['CELL.press']) md.pressure = params['CELL.press'];

    // Thermostat
    if (params['IONS.ion_temperature']) {
      const therm = params['IONS.ion_temperature'].toLowerCase();
      if (therm.includes('nose')) md.thermostat = 'NOSE_HOOVER';
      else if (therm.includes('berendsen')) md.thermostat = 'BERENDSEN';
      else if (therm.includes('langevin')) md.thermostat = 'LANGEVIN';
      else if (therm.includes('rescaling')) md.thermostat = 'VELOCITY_SCALING';
    }

    return md;
  }

  /**
   * Extract QE optimization parameters
   */
  private extractQEOpt(params: Record<string, any>): OptimizationParameters {
    const opt: OptimizationParameters = {};

    // Algorithm
    if (params['IONS.ion_dynamics']) {
      const dyn = params['IONS.ion_dynamics'].toLowerCase();
      if (dyn.includes('bfgs')) opt.algorithm = 'BFGS';
      else if (dyn.includes('lbfgs')) opt.algorithm = 'LBFGS';
      else if (dyn.includes('cg')) opt.algorithm = 'CG';
      else if (dyn.includes('damp')) opt.algorithm = 'DAMPED_MD';
      else if (dyn.includes('verlet')) opt.algorithm = 'MD';
      else opt.algorithm = 'BFGS';
    }

    // Max steps
    if (params['IONS.bfgs_ndim']) opt.maxSteps = params['IONS.bfgs_ndim'];

    // Convergence
    if (params['ELECTRONS.conv_thr']) opt.energyConv = params['ELECTRONS.conv_thr'];

    // Cell optimization
    if (params['CELL.cell_dynamics']) {
      opt.optimizeCell = params['CELL.cell_dynamics'].toLowerCase() !== 'none';
      if (params['CELL.press']) opt.cellPressure = params['CELL.press'];
    }

    return opt;
  }

  /**
   * Extract MD/optimization workflow from CP2K input
   */
  private extractCP2KWorkflow(filepath: string): MDWorkflow {
    const content = fs.readFileSync(filepath, 'utf-8');
    const params: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    let currentSection = '';
    let sectionPath: string[] = [];

    content.split('\n').forEach(line => {
      const trimmed = line.trim();

      if (trimmed.startsWith('&') && !trimmed.startsWith('&END')) {
        const sectionName = trimmed.substring(1).split(/\s+/)[0].toUpperCase();
        sectionPath.push(sectionName);
        currentSection = sectionPath.join('.');
        return;
      }
      if (trimmed.startsWith('&END')) {
        sectionPath.pop();
        currentSection = sectionPath.join('.');
        return;
      }
      if (!trimmed || trimmed.startsWith('#') || trimmed.startsWith('!')) return;

      const parts = trimmed.split(/\s+/);
      if (parts.length >= 2) {
        const key = parts[0].toUpperCase();
        const value = parts.slice(1).join(' ');
        params[`${currentSection}.${key}`] = this.parseValue(value);
      }
    });

    // Detect workflow type
    const hasMD = params['MOTION.MD.TIMESTEP'] || params['MOTION.MD.STEPS'];
    const hasOpt = params['MOTION.GEO_OPT.OPTIMIZER'] || params['MOTION.GEO_OPT.MAX_ITER'];

    const workflow: MDWorkflow = {
      type: hasMD ? 'md' : (hasOpt ? 'optimization' : 'md'),
      sourceFormat: 'cp2k',
      metadata: { filepath, warnings, unsupported }
    };

    if (hasMD || !hasOpt) {
      workflow.md = this.extractCP2KMD(params);
    }
    if (hasOpt) {
      workflow.optimization = this.extractCP2KOpt(params);
    }

    return workflow;
  }

  /**
   * Extract CP2K MD parameters
   */
  private extractCP2KMD(params: Record<string, any>): MDParameters {
    const md: MDParameters = {};

    // Time step
    if (params['MOTION.MD.TIMESTEP']) md.timeStep = params['MOTION.MD.TIMESTEP'];
    if (params['MOTION.MD.STEPS']) md.nSteps = params['MOTION.MD.STEPS'];

    // Ensemble
    if (params['MOTION.MD.ENSEMBLE']) {
      const ens = params['MOTION.MD.ENSEMBLE'].toUpperCase();
      if (ens.includes('NVE')) md.ensemble = 'NVE';
      else if (ens.includes('NVT')) md.ensemble = 'NVT';
      else if (ens.includes('NPT')) md.ensemble = 'NPT';
      else if (ens.includes('NPH')) md.ensemble = 'NPH';
      else md.ensemble = 'NVT';
    }

    // Temperature
    if (params['MOTION.MD.TEMPERATURE']) md.temperature = params['MOTION.MD.TEMPERATURE'];
    if (params['MOTION.MD.TEMPTOL']) md.temperatureEnd = params['MOTION.MD.TEMPTOL'];

    // Thermostat
    if (params['MOTION.MD.THERMOSTAT.TYPE']) {
      const therm = params['MOTION.MD.THERMOSTAT.TYPE'].toUpperCase();
      if (therm.includes('NOSE')) md.thermostat = 'NOSE_HOOVER';
      else if (therm.includes('BERENDSEN')) md.thermostat = 'BERENDSEN';
      else if (therm.includes('LANGEVIN')) md.thermostat = 'LANGEVIN';
      else md.thermostat = 'NOSE_HOOVER';
    }

    // Pressure
    if (params['MOTION.MD.PRESSURE']) md.pressure = params['MOTION.MD.PRESSURE'];
    if (params['MOTION.MD.BAROSTAT.TYPE']) {
      const baro = params['MOTION.MD.BAROSTAT.TYPE'].toUpperCase();
      if (baro.includes('NOSE')) md.barostat = 'NOSE_HOOVER';
      else if (baro.includes('BERENDSEN')) md.barostat = 'BERENDSEN';
      else if (baro.includes('MTTK')) md.barostat = 'MTTK';
    }

    return md;
  }

  /**
   * Extract CP2K optimization parameters
   */
  private extractCP2KOpt(params: Record<string, any>): OptimizationParameters {
    const opt: OptimizationParameters = {};

    // Algorithm
    if (params['MOTION.GEO_OPT.OPTIMIZER']) {
      const algo = params['MOTION.GEO_OPT.OPTIMIZER'].toUpperCase();
      if (algo.includes('BFGS')) opt.algorithm = 'BFGS';
      else if (algo.includes('LBFGS')) opt.algorithm = 'LBFGS';
      else if (algo.includes('CG')) opt.algorithm = 'CG';
      else opt.algorithm = 'BFGS';
    }

    // Max steps
    if (params['MOTION.GEO_OPT.MAX_ITER']) opt.maxSteps = params['MOTION.GEO_OPT.MAX_ITER'];

    // Convergence
    if (params['MOTION.GEO_OPT.CONVERGENCE.EPS_ENERGY']) {
      opt.energyConv = params['MOTION.GEO_OPT.CONVERGENCE.EPS_ENERGY'];
    }
    if (params['MOTION.GEO_OPT.CONVERGENCE.EPS_FORCE']) {
      opt.forceConv = params['MOTION.GEO_OPT.CONVERGENCE.EPS_FORCE'];
    }

    // Cell optimization
    if (params['MOTION.CELL_OPT.MAX_ITER']) {
      opt.optimizeCell = true;
      if (params['MOTION.CELL_OPT.PRESSURE']) {
        opt.cellPressure = params['MOTION.CELL_OPT.PRESSURE'];
      }
    }

    return opt;
  }

  /**
   * Extract MD/optimization workflow from file (auto-detect format)
   */
  public extractWorkflow(filepath: string): MDWorkflow {
    const ext = path.extname(filepath).toLowerCase();
    const basename = path.basename(filepath).toUpperCase();

    // Detect format
    let format = 'unknown';
    if (basename === 'INCAR' || basename.startsWith('INCAR-')) {
      format = 'vasp';
    } else if (ext === '.in' || ext === '.pw') {
      format = 'qe';
    } else if (ext === '.inp' && !basename.includes('ORCA')) {
      format = 'cp2k';
    }

    // Extract based on format
    switch (format) {
      case 'vasp':
        return this.extractVASPWorkflow(filepath);
      case 'qe':
        return this.extractQEWorkflow(filepath);
      case 'cp2k':
        return this.extractCP2KWorkflow(filepath);
      default:
        return {
          type: 'md',
          sourceFormat: 'unknown',
          metadata: {
            filepath,
            warnings: [],
            unsupported: [`Unsupported format: ${format}`]
          }
        };
    }
  }

  /**
   * Convert MD parameters from VASP to QE
   */
  private convertVASPtoQEMD(md: MDParameters): MDParameters {
    const target: MDParameters = {};

    if (md.timeStep) target.timeStep = md.timeStep;
    if (md.nSteps) target.nSteps = md.nSteps;
    if (md.ensemble) target.ensemble = md.ensemble;
    if (md.temperature) target.temperature = md.temperature;
    if (md.temperatureEnd) target.temperatureEnd = md.temperatureEnd;
    if (md.pressure) target.pressure = md.pressure;
    if (md.thermostat) target.thermostat = md.thermostat;
    if (md.randomSeed) target.randomSeed = md.randomSeed;

    return target;
  }

  /**
   * Convert MD parameters from VASP to CP2K
   */
  private convertVASPtoCP2KMD(md: MDParameters): MDParameters {
    const target: MDParameters = {};

    if (md.timeStep) target.timeStep = md.timeStep;
    if (md.nSteps) target.nSteps = md.nSteps;
    if (md.ensemble) target.ensemble = md.ensemble;
    if (md.temperature) target.temperature = md.temperature;
    if (md.temperatureEnd) target.temperatureEnd = md.temperatureEnd;
    // kbar to bar
    if (md.pressure) target.pressure = md.pressure * 1000;
    if (md.thermostat) target.thermostat = md.thermostat;
    if (md.barostat) target.barostat = md.barostat;

    return target;
  }

  /**
   * Convert MD parameters from QE to VASP
   */
  private convertQEtoVASPMD(md: MDParameters): MDParameters {
    const target: MDParameters = {};

    if (md.timeStep) target.timeStep = md.timeStep;
    if (md.nSteps) target.nSteps = md.nSteps;
    if (md.ensemble) target.ensemble = md.ensemble;
    if (md.temperature) target.temperature = md.temperature;
    if (md.temperatureEnd) target.temperatureEnd = md.temperatureEnd;
    if (md.pressure) target.pressure = md.pressure;
    if (md.thermostat) {
      // QE Langevin -> VASP Nose-Hoover (no direct equivalent)
      if (md.thermostat === 'LANGEVIN') target.thermostat = 'NOSE_HOOVER';
      else target.thermostat = md.thermostat;
    }

    return target;
  }

  /**
   * Convert optimization parameters from VASP to QE
   */
  private convertVASPtoQEOpt(opt: OptimizationParameters): OptimizationParameters {
    const target: OptimizationParameters = {};

    // Algorithm mapping
    if (opt.algorithm) {
      const algoMap: Record<OptimizationAlgorithm, OptimizationAlgorithm> = {
        'CG': 'CG',
        'BFGS': 'BFGS',
        'LBFGS': 'LBFGS',
        'FIRE': 'BFGS',
        'MD': 'MD',
        'RMMDIIS': 'BFGS',
        'DAMPED_MD': 'DAMPED_MD'
      };
      target.algorithm = algoMap[opt.algorithm] || 'BFGS';
    }

    if (opt.maxSteps) target.maxSteps = opt.maxSteps;
    // eV to Ry
    if (opt.energyConv) target.energyConv = opt.energyConv * 0.0734986;
    if (opt.forceConv) target.forceConv = opt.forceConv;
    if (opt.optimizeCell !== undefined) target.optimizeCell = opt.optimizeCell;
    if (opt.cellPressure) target.cellPressure = opt.cellPressure;

    return target;
  }

  /**
   * Convert optimization parameters from VASP to CP2K
   */
  private convertVASPtoCP2KOpt(opt: OptimizationParameters): OptimizationParameters {
    const target: OptimizationParameters = {};

    if (opt.algorithm) {
      const algoMap: Record<OptimizationAlgorithm, OptimizationAlgorithm> = {
        'CG': 'CG',
        'BFGS': 'BFGS',
        'LBFGS': 'LBFGS',
        'FIRE': 'BFGS',
        'MD': 'CG',
        'RMMDIIS': 'BFGS',
        'DAMPED_MD': 'CG'
      };
      target.algorithm = algoMap[opt.algorithm] || 'BFGS';
    }

    if (opt.maxSteps) target.maxSteps = opt.maxSteps;
    // eV to hartree
    if (opt.energyConv) target.energyConv = opt.energyConv * 0.0367493;
    // eV/Å to hartree/bohr
    if (opt.forceConv) target.forceConv = opt.forceConv * 0.01943;
    if (opt.optimizeCell !== undefined) target.optimizeCell = opt.optimizeCell;
    // kbar to bar
    if (opt.cellPressure) target.cellPressure = opt.cellPressure * 1000;

    return target;
  }

  /**
   * Convert optimization parameters from QE to VASP
   */
  private convertQEtoVASPOpt(opt: OptimizationParameters): OptimizationParameters {
    const target: OptimizationParameters = {};

    if (opt.algorithm) {
      const algoMap: Record<OptimizationAlgorithm, OptimizationAlgorithm> = {
        'CG': 'CG',
        'BFGS': 'BFGS',
        'LBFGS': 'BFGS',
        'FIRE': 'CG',
        'MD': 'MD',
        'RMMDIIS': 'RMMDIIS',
        'DAMPED_MD': 'DAMPED_MD'
      };
      target.algorithm = algoMap[opt.algorithm] || 'CG';
    }

    if (opt.maxSteps) target.maxSteps = opt.maxSteps;
    // Ry to eV
    if (opt.energyConv) target.energyConv = opt.energyConv * 13.6057;
    if (opt.forceConv) target.forceConv = opt.forceConv;
    if (opt.optimizeCell !== undefined) target.optimizeCell = opt.optimizeCell;
    if (opt.cellPressure) target.cellPressure = opt.cellPressure;

    return target;
  }

  /**
   * Convert workflow from source to target format
   */
  public convertWorkflow(
    source: MDWorkflow,
    targetFormat: string
  ): MDWorkflowConversionResult {
    const result: MDWorkflowConversionResult = {
      success: false,
      errors: [],
      warnings: []
    };

    result.source = source;

    // Create target workflow
    const target: MDWorkflow = {
      type: source.type,
      sourceFormat: targetFormat,
      metadata: {
        filepath: '',
        warnings: [],
        unsupported: []
      }
    };

    // Convert based on type
    if (source.type === 'md' && source.md) {
      switch (`${source.sourceFormat}->${targetFormat}`) {
        case 'vasp->qe':
          target.md = this.convertVASPtoQEMD(source.md);
          break;
        case 'vasp->cp2k':
          target.md = this.convertVASPtoCP2KMD(source.md);
          break;
        case 'qe->vasp':
          target.md = this.convertQEtoVASPMD(source.md);
          break;
        case 'cp2k->vasp':
          target.md = { ...source.md };
          if (target.md.pressure) target.md.pressure = target.md.pressure / 1000;
          break;
        default:
          result.warnings.push(`Direct conversion ${source.sourceFormat}->${targetFormat} not fully supported, using generic mapping`);
          target.md = { ...source.md };
      }
    } else if (source.type === 'optimization' && source.optimization) {
      switch (`${source.sourceFormat}->${targetFormat}`) {
        case 'vasp->qe':
          target.optimization = this.convertVASPtoQEOpt(source.optimization);
          break;
        case 'vasp->cp2k':
          target.optimization = this.convertVASPtoCP2KOpt(source.optimization);
          break;
        case 'qe->vasp':
          target.optimization = this.convertQEtoVASPOpt(source.optimization);
          break;
        default:
          result.warnings.push(`Direct conversion ${source.sourceFormat}->${targetFormat} not fully supported, using generic mapping`);
          target.optimization = { ...source.optimization };
      }
    }

    result.target = target;
    result.success = true;

    return result;
  }

  /**
   * Generate VASP INCAR MD section from MD parameters
   */
  public generateVASPMDSection(md: MDParameters): string {
    const lines: string[] = [];

    lines.push('# MD Parameters (migrated from OpenQC)');
    if (md.timeStep) lines.push(`POTIM = ${md.timeStep.toFixed(1)}`);
    if (md.nSteps) lines.push(`NSW = ${md.nSteps}`);

    lines.push('IBRION = 0');  // MD

    // ISIF for ensemble
    if (md.ensemble === 'NVT') {
      lines.push('ISIF = 2');
    } else if (md.ensemble === 'NPT') {
      lines.push('ISIF = 3');
    }

    if (md.temperature) lines.push(`TEBEG = ${md.temperature}`);
    if (md.temperatureEnd) lines.push(`TEEND = ${md.temperatureEnd}`);
    if (md.pressure) lines.push(`PSTRESS = ${md.pressure.toFixed(1)}`);

    if (md.thermostat) {
      switch (md.thermostat) {
        case 'NOSE_HOOVER':
          lines.push('SMASS = 0');
          break;
        case 'BERENDSEN':
          lines.push('SMASS = 2');
          break;
        case 'ANDERSEN':
          lines.push('SMASS = 3');
          break;
      }
    }

    return lines.join('\n');
  }

  /**
   * Generate QE MD section from MD parameters
   */
  public generateQEMDSection(md: MDParameters): string {
    const lines: string[] = [];

    lines.push('&CONTROL');
    if (md.timeStep) lines.push(`  dt = ${md.timeStep}`);
    if (md.nSteps) lines.push(`  nstep = ${md.nSteps}`);
    lines.push('/');

    lines.push('&IONS');
    if (md.ensemble === 'NVE') {
      lines.push('  ion_dynamics = \'verlet\'');
    } else if (md.thermostat === 'LANGEVIN') {
      lines.push('  ion_dynamics = \'langevin\'');
    } else {
      lines.push('  ion_dynamics = \'verlet\'');
    }

    if (md.temperature) lines.push(`  tempw = ${md.temperature}`);
    if (md.temperatureEnd) lines.push(`  templ = ${md.temperatureEnd}`);

    if (md.thermostat && md.ensemble !== 'NVE') {
      switch (md.thermostat) {
        case 'NOSE_HOOVER':
          lines.push('  ion_temperature = \'nose\'');
          break;
        case 'BERENDSEN':
          lines.push('  ion_temperature = \'berendsen\'');
          break;
        case 'LANGEVIN':
          lines.push('  ion_temperature = \'langevin\'');
          break;
      }
    }
    lines.push('/');

    if (md.ensemble === 'NPT' && md.pressure) {
      lines.push('&CELL');
      lines.push(`  press = ${md.pressure}`);
      lines.push('  cell_dynamics = \'pr\'');
      lines.push('/');
    }

    return lines.join('\n');
  }

  /**
   * Generate CP2K MD section from MD parameters
   */
  public generateCP2KMDSection(md: MDParameters): string {
    const lines: string[] = [];

    lines.push('&MOTION');
    lines.push('  &MD');

    if (md.ensemble) lines.push(`    ENSEMBLE ${md.ensemble}`);
    if (md.timeStep) lines.push(`    TIMESTEP ${md.timeStep}`);
    if (md.nSteps) lines.push(`    STEPS ${md.nSteps}`);
    if (md.temperature) lines.push(`    TEMPERATURE ${md.temperature}`);

    if (md.thermostat && md.ensemble !== 'NVE') {
      lines.push('    &THERMOSTAT');
      switch (md.thermostat) {
        case 'NOSE_HOOVER':
          lines.push('      TYPE NOSE');
          break;
        case 'BERENDSEN':
          lines.push('      TYPE BERENDSEN');
          break;
        case 'LANGEVIN':
          lines.push('      TYPE LANGEVIN');
          break;
      }
      lines.push('    &END THERMOSTAT');
    }

    if (md.ensemble === 'NPT') {
      lines.push('    &BAROSTAT');
      lines.push('      TYPE NOSE');
      if (md.pressure) lines.push(`      PRESSURE ${md.pressure}`);
      lines.push('    &END BAROSTAT');
    }

    lines.push('  &END MD');
    lines.push('&END MOTION');

    return lines.join('\n');
  }

  /**
   * Generate VASP INCAR optimization section
   */
  public generateVASPOptSection(opt: OptimizationParameters): string {
    const lines: string[] = [];

    lines.push('# Optimization Parameters (migrated from OpenQC)');

    if (opt.algorithm) {
      const algoMap: Record<OptimizationAlgorithm, number> = {
        'CG': 2,
        'BFGS': 5,
        'LBFGS': 5,
        'FIRE': 2,
        'MD': 7,
        'RMMDIIS': 1,
        'DAMPED_MD': 3
      };
      lines.push(`IBRION = ${algoMap[opt.algorithm] || 2}`);
    }

    if (opt.maxSteps) lines.push(`NSW = ${opt.maxSteps}`);
    if (opt.energyConv) lines.push(`EDIFF = ${opt.energyConv.toExponential(1)}`);
    if (opt.forceConv) lines.push(`EDIFFG = -${opt.forceConv.toFixed(2)}`);

    if (opt.optimizeCell) {
      lines.push('ISIF = 3');
      if (opt.cellPressure) lines.push(`PSTRESS = ${opt.cellPressure}`);
    } else {
      lines.push('ISIF = 2');
    }

    return lines.join('\n');
  }

  /**
   * Generate QE optimization section
   */
  public generateQEOptSection(opt: OptimizationParameters): string {
    const lines: string[] = [];

    lines.push('&IONS');
    if (opt.algorithm) {
      const algoMap: Record<OptimizationAlgorithm, string> = {
        'CG': 'cg',
        'BFGS': 'bfgs',
        'LBFGS': 'lbfgs',
        'FIRE': 'fire',
        'MD': 'verlet',
        'RMMDIIS': 'bfgs',
        'DAMPED_MD': 'damp'
      };
      lines.push(`  ion_dynamics = '${algoMap[opt.algorithm] || 'bfgs'}'`);
    }
    if (opt.maxSteps) lines.push(`  bfgs_ndim = ${opt.maxSteps}`);
    lines.push('/');

    if (opt.optimizeCell) {
      lines.push('&CELL');
      lines.push('  cell_dynamics = \'bfgs\'');
      if (opt.cellPressure) lines.push(`  press = ${opt.cellPressure}`);
      lines.push('/');
    }

    return lines.join('\n');
  }

  /**
   * Generate CP2K optimization section
   */
  public generateCP2KOptSection(opt: OptimizationParameters): string {
    const lines: string[] = [];

    if (opt.optimizeCell) {
      lines.push('&MOTION');
      lines.push('  &CELL_OPT');
      if (opt.algorithm) lines.push(`    OPTIMIZER ${opt.algorithm}`);
      if (opt.maxSteps) lines.push(`    MAX_ITER ${opt.maxSteps}`);
      if (opt.cellPressure) lines.push(`    PRESSURE ${opt.cellPressure}`);
      lines.push('  &END CELL_OPT');
      lines.push('&END MOTION');
    } else {
      lines.push('&MOTION');
      lines.push('  &GEO_OPT');
      if (opt.algorithm) lines.push(`    OPTIMIZER ${opt.algorithm}`);
      if (opt.maxSteps) lines.push(`    MAX_ITER ${opt.maxSteps}`);
      lines.push('    &CONVERGENCE');
      if (opt.energyConv) lines.push(`      EPS_ENERGY ${opt.energyConv}`);
      if (opt.forceConv) lines.push(`      EPS_FORCE ${opt.forceConv}`);
      lines.push('    &END CONVERGENCE');
      lines.push('  &END GEO_OPT');
      lines.push('&END MOTION');
    }

    return lines.join('\n');
  }

  /**
   * Format MD parameters for display
   */
  public formatMDParameters(md: MDParameters): string {
    const lines: string[] = [];

    if (md.timeStep) lines.push(`Time Step: ${md.timeStep} fs`);
    if (md.nSteps) lines.push(`Number of Steps: ${md.nSteps}`);
    if (md.ensemble) lines.push(`Ensemble: ${md.ensemble}`);
    if (md.temperature) lines.push(`Temperature: ${md.temperature} K`);
    if (md.temperatureEnd) lines.push(`Final Temperature: ${md.temperatureEnd} K`);
    if (md.pressure) lines.push(`Pressure: ${md.pressure} kbar`);
    if (md.thermostat) lines.push(`Thermostat: ${md.thermostat}`);
    if (md.barostat) lines.push(`Barostat: ${md.barostat}`);

    return lines.join('\n');
  }

  /**
   * Format optimization parameters for display
   */
  public formatOptimizationParameters(opt: OptimizationParameters): string {
    const lines: string[] = [];

    if (opt.algorithm) lines.push(`Algorithm: ${opt.algorithm}`);
    if (opt.maxSteps) lines.push(`Max Steps: ${opt.maxSteps}`);
    if (opt.energyConv) lines.push(`Energy Convergence: ${opt.energyConv} eV`);
    if (opt.forceConv) lines.push(`Force Convergence: ${opt.forceConv} eV/Å`);
    if (opt.optimizeCell !== undefined) lines.push(`Cell Optimization: ${opt.optimizeCell}`);
    if (opt.cellPressure) lines.push(`Cell Pressure: ${opt.cellPressure} kbar`);

    return lines.join('\n');
  }
}

/**
 * Quick MD workflow migration command handler
 */
export async function quickMigrateMDWorkflow(
  context: vscode.ExtensionContext
): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const sourcePath = editor.document.uri.fsPath;
  const converter = new MDWorkflowConverter();

  // Extract workflow
  const workflow = converter.extractWorkflow(sourcePath);

  if (workflow.sourceFormat === 'unknown') {
    vscode.window.showErrorMessage('Could not detect file format');
    return;
  }

  // Select target format
  const targetFormat = await vscode.window.showQuickPick(
    ['vasp', 'cp2k', 'qe', 'lammps'].map(format => ({
      label: format.toUpperCase(),
      description: `Migrate to ${format.toUpperCase()}`,
      value: format,
    })),
    {
      placeHolder: 'Select target format',
    }
  );

  if (!targetFormat) {
    return;
  }

  // Convert
  const result = converter.convertWorkflow(workflow, targetFormat.value);

  if (!result.success || !result.target) {
    vscode.window.showErrorMessage(`Migration failed: ${result.errors.join(', ')}`);
    return;
  }

  // Generate output
  let output = '';
  if (result.target.type === 'md' && result.target.md) {
    switch (targetFormat.value) {
      case 'vasp':
        output = converter.generateVASPMDSection(result.target.md);
        break;
      case 'qe':
        output = converter.generateQEMDSection(result.target.md);
        break;
      case 'cp2k':
        output = converter.generateCP2KMDSection(result.target.md);
        break;
    }
  } else if (result.target.type === 'optimization' && result.target.optimization) {
    switch (targetFormat.value) {
      case 'vasp':
        output = converter.generateVASPOptSection(result.target.optimization);
        break;
      case 'qe':
        output = converter.generateQEOptSection(result.target.optimization);
        break;
      case 'cp2k':
        output = converter.generateCP2KOptSection(result.target.optimization);
        break;
    }
  }

  // Show result
  const action = await vscode.window.showInformationMessage(
    `Migration successful! Type: ${result.target.type}`,
    'Copy to Clipboard',
    'View Output',
    'Open Preview'
  );

  if (action === 'Copy to Clipboard') {
    vscode.env.clipboard.writeText(output);
    vscode.window.showInformationMessage('MD/Optimization parameters copied to clipboard');
  } else if (action === 'View Output' || action === 'Open Preview') {
    const doc = await vscode.workspace.openTextDocument({
      content: output,
      language: targetFormat.value === 'vasp' ? 'vasp' : 'plaintext'
    });
    await vscode.window.showTextDocument(doc);
  }
}
