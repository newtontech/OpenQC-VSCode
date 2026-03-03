/**
 * ASE Integration Module
 *
 * Exports for ASE (Atomic Simulation Environment) integration
 * providing cross-code workflow migration and structure conversion.
 */

export { ASEConverter, ASEFormat, ASEAtoms, ConversionResult } from './ASEConverter';
export { registerASECommands } from './commands';
export {
  ASECalculator,
  CalculatorFactory,
  CalculatorType,
  CalculatorConfig,
  CalculationResult,
  InputGenerationResult,
  getCalculatorConfiguration,
} from './ASECalculator';
export { registerCalculatorCommands } from './calculatorCommands';
export {
  ComplexPropertyMapper,
  PropertyMapperFactory,
  HubbardUParams,
  AtomConstraint,
  ExcitedStateConfig,
  PseudopotentialStrategy,
  BasisSetStrategy,
  PropertyMappingResult,
  getSupportedComplexProperties,
} from './ComplexPropertyMapper';
