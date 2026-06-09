/**
 * OpenQCResults - Canonical DTO for calculation output results.
 *
 * @module results/OpenQCResults
 */

export const OPENQC_RESULTS_SCHEMA_VERSION = 'openqc.results.v1' as const;

export interface EnergyValue {
  value: number;
  unit: 'eV' | 'hartree' | 'kJ/mol' | 'kcal/mol';
}

export interface OptimizationInfo {
  energies?: number[];
  frames?: any[]; // OpenQCStructure[] when available
  converged?: boolean;
  stepCount?: number;
}

export interface OrbitalInfo {
  homo?: number;
  lumo?: number;
  energies?: number[];
}

export interface OpenQCResults {
  schemaVersion: typeof OPENQC_RESULTS_SCHEMA_VERSION;
  sourceFile?: string;
  software?: string;
  success?: boolean;
  finalEnergy?: EnergyValue;
  scfEnergies?: number[];
  optimization?: OptimizationInfo;
  orbitals?: OrbitalInfo;
  frequencies?: number[];
  charges?: Record<string, number[]>;
  dipole?: [number, number, number];
  warnings?: string[];
  errors?: string[];
  cclibAvailable?: boolean;
  finalStructure?: any;
  metadata?: Record<string, unknown>;
}
