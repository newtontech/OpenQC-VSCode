/**
 * Complex Property Mapper
 *
 * Maps complex quantum chemistry properties between different codes including:
 * - Hubbard U parameter mapping
 * - Constraints and frozen atoms
 * - Excited state methods
 * - Pseudopotential and basis set strategies
 */

import * as vscode from 'vscode';
import { CalculatorType } from './ASECalculator';

/**
 * Hubbard U parameters for a single element
 */
export interface HubbardUParams {
  element: string;
  l: number;  // Angular momentum (0=s, 1=p, 2=d, 3=f)
  u: number;  // U value in eV
  j?: number; // J value in eV (for DFT+U+J)
}

/**
 * Constraint definition
 */
export interface AtomConstraint {
  atomIndex: number;
  fixX?: boolean;
  fixY?: boolean;
  fixZ?: boolean;
  fixAll?: boolean;
}

/**
 * Excited state method configuration
 */
export interface ExcitedStateConfig {
  method: 'TDDFT' | 'GW' | 'BSE' | 'CIS' | 'CISD' | 'EOM-CCSD';
  nStates: number;
  spinState?: 'singlet' | 'triplet' | 'both';
  tda?: boolean;  // Tamm-Dancoff approximation
}

/**
 * Pseudopotential strategy
 */
export interface PseudopotentialStrategy {
  type: 'normconserving' | 'ultrasoft' | 'paw' | 'hartwigner';
  functional: string;
  accuracy: 'low' | 'medium' | 'high' | 'very_high';
  relativistic?: boolean;
}

/**
 * Basis set strategy
 */
export interface BasisSetStrategy {
  type: 'plane_wave' | 'gaussian' | 'numerical' | 'augmented';
  quality: 'minimal' | 'dzvp' | 'tzvp' | 'qzvp' | 'custom';
  customBasis?: string;
}

/**
 * Complex property mapping result
 */
export interface PropertyMappingResult {
  success: boolean;
  mappedProperties: Record<string, any>;
  warnings: string[];
  unsupported: string[];
  error?: string;
}

/**
 * Complex Property Mapper class
 *
 * Handles mapping of complex quantum chemistry properties between codes
 */
export class ComplexPropertyMapper {
  private sourceType: CalculatorType;
  private targetType: CalculatorType;

  constructor(sourceType: CalculatorType, targetType: CalculatorType) {
    this.sourceType = sourceType;
    this.targetType = targetType;
  }

  /**
   * Map Hubbard U parameters between codes
   */
  public mapHubbardU(params: HubbardUParams[]): PropertyMappingResult {
    const mapped: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    switch (this.targetType) {
      case CalculatorType.VASP:
        mapped.ldau = true;
        mapped.ldautype = 2;  // Dudarev scheme
        mapped.ldau_luj = {} as Record<string, { L: number; U: number; J: number }>;

        for (const p of params) {
          const symbol = p.element;
          mapped.ldau_luj[symbol] = {
            L: p.l,
            U: p.u,
            J: p.j || 0,
          };
        }
        break;

      case CalculatorType.CP2K:
        mapped.dft_plus_u = true;
        mapped.plus_u_parms = [] as Array<{ element: string; l: number; u: number }>;

        for (const p of params) {
          mapped.plus_u_parms.push({
            element: p.element,
            l: p.l,
            u: p.u,
          });

          if (p.j !== undefined) {
            warnings.push('CP2K does not support J parameter in DFT+U');
          }
        }
        break;

      case CalculatorType.QE:
        mapped.hubbard_u = {} as Record<string, number>;

        for (const p of params) {
          const hubbardL = ['s', 'p', 'd', 'f'][p.l];
          const key = `hubbard_u(${hubbardL}-${p.element})`;
          mapped.hubbard_u[key] = p.u;

          if (p.j !== undefined) {
            warnings.push('Quantum ESPRESSO J parameter requires special handling');
          }
        }
        break;

      default:
        unsupported.push(`Hubbard U for ${this.targetType}`);
    }

    return {
      success: unsupported.length === 0,
      mappedProperties: mapped,
      warnings,
      unsupported,
    };
  }

  /**
   * Map atom constraints between codes
   */
  public mapConstraints(constraints: AtomConstraint[]): PropertyMappingResult {
    const mapped: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    switch (this.targetType) {
      case CalculatorType.VASP:
        // VASP uses selective dynamics in POSCAR
        mapped.selective_dynamics = true;
        mapped.fixed_atoms = [] as number[];
        mapped.selective_flags = [] as boolean[][];

        for (const c of constraints) {
          if (c.fixAll) {
            mapped.fixed_atoms.push(c.atomIndex);
            mapped.selective_flags.push([false, false, false]);
          } else {
            mapped.selective_flags.push([
              !c.fixX,
              !c.fixY,
              !c.fixZ,
            ]);
          }
        }
        break;

      case CalculatorType.CP2K:
        // CP2K uses CONSTRAINT or FIXED_ATOMS section
        mapped.fixed_atoms = [] as number[];

        for (const c of constraints) {
          if (c.fixAll) {
            mapped.fixed_atoms.push(c.atomIndex + 1);  // CP2K uses 1-based indexing
          } else if (c.fixX || c.fixY || c.fixZ) {
            warnings.push('CP2K partial coordinate constraints require COLVAR');
          }
        }

        if (mapped.fixed_atoms.length > 0) {
          mapped.constraints_section = `&CONSTRAINT
  &FIXED_ATOMS
    LIST ${mapped.fixed_atoms.join(' ')}
  &END FIXED_ATOMS
&END CONSTRAINT`;
        }
        break;

      case CalculatorType.QE:
        // QE uses constraints in &IONS namelist
        mapped.if_pos = [] as number[][];

        for (const c of constraints) {
          if (c.fixAll) {
            mapped.if_pos.push([0, 0, 0]);
          } else {
            mapped.if_pos.push([
              c.fixX ? 0 : 1,
              c.fixY ? 0 : 1,
              c.fixZ ? 0 : 1,
            ]);
          }
        }
        break;

      default:
        unsupported.push(`Constraints for ${this.targetType}`);
    }

    return {
      success: unsupported.length === 0,
      mappedProperties: mapped,
      warnings,
      unsupported,
    };
  }

  /**
   * Map excited state methods between codes
   */
  public mapExcitedStateMethod(config: ExcitedStateConfig): PropertyMappingResult {
    const mapped: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    // Check if source supports the method
    const sourceSupport = this.checkExcitedStateSupport(this.sourceType, config.method);
    if (!sourceSupport.supported) {
      warnings.push(`${config.method} may not be fully supported in ${this.sourceType}`);
    }

    switch (this.targetType) {
      case CalculatorType.VASP:
        if (config.method === 'TDDFT') {
          mapped.ismeartype = 'TDDFT';
          mapped.nedos = 5000;

          // VASP uses different implementations
          if (config.spinState === 'singlet') {
            mapped.lralpha = true;
          } else if (config.spinState === 'triplet') {
            mapped.lbeta = true;
          }

          mapped.nomega = config.nStates;
        } else {
          unsupported.push(`${config.method} in VASP`);
        }
        break;

      case CalculatorType.CP2K:
        if (config.method === 'TDDFT') {
          mapped.run_type = 'tddft';
          mapped.nstates = config.nStates;
          mapped.tda = config.tda !== false;  // Default to TDA

          if (config.spinState === 'both') {
            mapped.rks_triplet = true;
          }
        } else if (config.method === 'GW') {
          mapped.run_type = 'gw';
          warnings.push('GW in CP2K requires special compilation flags');
        } else {
          unsupported.push(`${config.method} in CP2K`);
        }
        break;

      case CalculatorType.QE:
        if (config.method === 'TDDFT' || config.method === 'BSE') {
          mapped.calculation = 'turbo_davidson';
          mapped.nstates = config.nStates;
          mapped.ipol = config.spinState === 'both' ? 4 : 1;
        } else {
          unsupported.push(`${config.method} in Quantum ESPRESSO`);
        }
        break;

      default:
        unsupported.push(`Excited state methods for ${this.targetType}`);
    }

    return {
      success: unsupported.length === 0,
      mappedProperties: mapped,
      warnings,
      unsupported,
    };
  }

  /**
   * Map pseudopotential strategy between codes
   */
  public mapPseudopotentialStrategy(
    strategy: PseudopotentialStrategy,
    elements: string[]
  ): PropertyMappingResult {
    const mapped: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    const pseudos: Record<string, string> = {};

    for (const element of elements) {
      switch (this.targetType) {
        case CalculatorType.VASP:
          // VASP POTCAR naming convention
          pseudos[element] = this.getVASPPseudopotential(element, strategy);
          break;

        case CalculatorType.CP2K:
          // CP2K GTH pseudopotentials
          pseudos[element] = this.getCP2KPseudopotential(element, strategy);
          break;

        case CalculatorType.QE:
          // QE UPF format
          pseudos[element] = this.getQEPseudopotential(element, strategy);
          break;

        default:
          unsupported.push(`Pseudopotentials for ${this.targetType}`);
      }
    }

    mapped.pseudopotentials = pseudos;

    return {
      success: unsupported.length === 0,
      mappedProperties: mapped,
      warnings,
      unsupported,
    };
  }

  /**
   * Map basis set strategy between codes
   */
  public mapBasisSetStrategy(
    strategy: BasisSetStrategy,
    elements: string[]
  ): PropertyMappingResult {
    const mapped: Record<string, any> = {};
    const warnings: string[] = [];
    const unsupported: string[] = [];

    const basisSets: Record<string, string> = {};

    for (const element of elements) {
      switch (this.targetType) {
        case CalculatorType.VASP:
          // VASP uses plane waves, no per-element basis
          if (strategy.type !== 'plane_wave') {
            warnings.push('VASP only supports plane wave basis');
          }
          mapped.encut = this.getVASPENCut(strategy.quality);
          break;

        case CalculatorType.CP2K:
          basisSets[element] = this.getCP2KBasisSet(element, strategy);
          break;

        case CalculatorType.QE:
          // QE uses plane waves
          if (strategy.type !== 'plane_wave') {
            warnings.push('Quantum ESPRESSO only supports plane wave basis');
          }
          mapped.ecutwfc = this.getQEENCut(strategy.quality);
          break;

        default:
          unsupported.push(`Basis sets for ${this.targetType}`);
      }
    }

    if (Object.keys(basisSets).length > 0) {
      mapped.basis_sets = basisSets;
    }

    return {
      success: unsupported.length === 0,
      mappedProperties: mapped,
      warnings,
      unsupported,
    };
  }

  /**
   * Check if calculator supports excited state method
   */
  private checkExcitedStateSupport(
    calculatorType: CalculatorType,
    method: string
  ): { supported: boolean; notes?: string } {
    const support: Record<CalculatorType, Record<string, { supported: boolean; notes?: string }>> = {
      [CalculatorType.VASP]: {
        TDDFT: { supported: true },
        GW: { supported: true },
        BSE: { supported: true },
        CIS: { supported: false, notes: 'Use TDDFT instead' },
        CISD: { supported: false },
        'EOM-CCSD': { supported: false },
      },
      [CalculatorType.CP2K]: {
        TDDFT: { supported: true },
        GW: { supported: true, notes: 'Requires specific compilation' },
        BSE: { supported: false },
        CIS: { supported: true },
        CISD: { supported: false },
        'EOM-CCSD': { supported: false },
      },
      [CalculatorType.QE]: {
        TDDFT: { supported: true },
        GW: { supported: true },
        BSE: { supported: true },
        CIS: { supported: false },
        CISD: { supported: false },
        'EOM-CCSD': { supported: false },
      },
    };

    return support[calculatorType]?.[method] || { supported: false };
  }

  /**
   * Get VASP pseudopotential name
   */
  private getVASPPseudopotential(
    element: string,
    strategy: PseudopotentialStrategy
  ): string {
    const suffixes: Record<PseudopotentialStrategy['accuracy'], string> = {
      low: '',
      medium: '_pv',
      high: '_sv',
      very_high: '_sv_GW',
    };

    const suffix = suffixes[strategy.accuracy];
    return `${element}${suffix}`;
  }

  /**
   * Get CP2K pseudopotential name
   */
  private getCP2KPseudopotential(
    element: string,
    strategy: PseudopotentialStrategy
  ): string {
    const func = strategy.functional.toUpperCase();
    return `GTH-${func}`;
  }

  /**
   * Get Quantum ESPRESSO pseudopotential name
   */
  private getQEPseudopotential(
    element: string,
    strategy: PseudopotentialStrategy
  ): string {
    const types: Record<PseudopotentialStrategy['type'], string> = {
      normconserving: '_nc',
      ultrasoft: '_us',
      paw: '_paw',
      hartwigner: '_hgh',
    };

    const typeSuffix = types[strategy.type] || '';
    const func = strategy.functional.toLowerCase();
    return `${element}${typeSuffix}.${func}.UPF`;
  }

  /**
   * Get VASP ENCUT based on quality
   */
  private getVASPENCut(quality: BasisSetStrategy['quality']): number {
    const encuts: Record<BasisSetStrategy['quality'], number> = {
      minimal: 250,
      dzvp: 400,
      tzvp: 520,
      qzvp: 650,
      custom: 400,
    };

    return encuts[quality];
  }

  /**
   * Get QE energy cutoff based on quality
   */
  private getQEENCut(quality: BasisSetStrategy['quality']): number {
    const encuts: Record<BasisSetStrategy['quality'], number> = {
      minimal: 20,
      dzvp: 40,
      tzvp: 60,
      qzvp: 80,
      custom: 40,
    };

    return encuts[quality];
  }

  /**
   * Get CP2K basis set name
   */
  private getCP2KBasisSet(
    element: string,
    strategy: BasisSetStrategy
  ): string {
    const bases: Record<BasisSetStrategy['quality'], string> = {
      minimal: 'SZV-MOLOPT-SR-GTH',
      dzvp: 'DZVP-MOLOPT-SR-GTH',
      tzvp: 'TZVP-MOLOPT-GTH',
      qzvp: 'QZVP-MOLOPT-GTH',
      custom: strategy.customBasis || 'DZVP-MOLOPT-SR-GTH',
    };

    return bases[strategy.quality];
  }
}

/**
 * Factory for creating property mappers
 */
export class PropertyMapperFactory {
  /**
   * Create a property mapper
   */
  createMapper(
    sourceType: CalculatorType,
    targetType: CalculatorType
  ): ComplexPropertyMapper {
    return new ComplexPropertyMapper(sourceType, targetType);
  }
}

/**
 * Get supported complex properties for a calculator
 */
export function getSupportedComplexProperties(
  calculatorType: CalculatorType
): {
  hubbardU: boolean;
  constraints: boolean;
  excitedStates: string[];
  pseudopotentials: boolean;
  basisSets: boolean;
} {
  const support: Record<CalculatorType, any> = {
    [CalculatorType.VASP]: {
      hubbardU: true,
      constraints: true,
      excitedStates: ['TDDFT', 'GW', 'BSE'],
      pseudopotentials: true,
      basisSets: true,
    },
    [CalculatorType.CP2K]: {
      hubbardU: true,
      constraints: true,
      excitedStates: ['TDDFT', 'GW'],
      pseudopotentials: true,
      basisSets: true,
    },
    [CalculatorType.QE]: {
      hubbardU: true,
      constraints: true,
      excitedStates: ['TDDFT', 'BSE'],
      pseudopotentials: true,
      basisSets: true,
    },
  };

  return (
    support[calculatorType] || {
      hubbardU: false,
      constraints: false,
      excitedStates: [],
      pseudopotentials: false,
      basisSets: false,
    }
  );
}
