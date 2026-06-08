/**
 * OpenQC Structures module.
 *
 * Re-exports the canonical DTO types, validation, and conversion utilities.
 *
 * @module structures
 */

export {
  OPENQC_STRUCTURE_SCHEMA_VERSION,
  type OpenQCStructure,
  type OpenQCStructureKind,
  type OpenQCAtom,
  type OpenQCCell,
  type OpenQCBond,
  type OpenQCStructureMetadata,
  type OpenQCStructureSchemaVersion,
} from './OpenQCStructure';

export {
  validateOpenQCStructure,
  type ValidationResult,
} from './validation';

export {
  openQCStructureToXYZ,
  xyzToOpenQCStructure,
  molecularStructureToOpenQCStructure,
  openQCStructureToMolecularStructure,
  poscarToOpenQCStructure,
  createOpenQCStructure,
} from './converters';
