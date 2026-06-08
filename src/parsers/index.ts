/**
 * OpenQC-VSCode Parsers
 * Input file parsers for quantum chemistry software
 *
 * Uses a registry pattern to satisfy the Open/Closed Principle:
 * adding a new parser only requires calling `registerParser()`,
 * not modifying this file's switch statement.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/38
 */

export * from './base';
export { VASPParser } from './VASPParser';
export { GaussianParser } from './GaussianParser';
export { ORCAParser } from './ORCAParser';
export { CP2KParser } from './CP2KParser';
export { QEParser } from './QEParser';
export { GAMESSParser } from './GAMESSParser';
export { NWChemParser } from './NWChemParser';

import { BaseParser } from './base';
import { VASPParser } from './VASPParser';
import { GaussianParser } from './GaussianParser';
import { ORCAParser } from './ORCAParser';
import { CP2KParser } from './CP2KParser';
import { QEParser } from './QEParser';
import { GAMESSParser } from './GAMESSParser';
import { NWChemParser } from './NWChemParser';
import { QuantumChemistrySoftware } from '../managers/FileTypeDetector';

/**
 * Parser factory function signature.
 * Each registered parser provides a factory that accepts (content, filename?) and returns a BaseParser.
 */
type ParserFactory = (content: string, filename?: string) => BaseParser;

/**
 * Internal registry mapping software names to parser factories.
 * Satisfies Open/Closed Principle: extend by registering new entries, not by editing switch cases.
 */
const parserRegistry = new Map<QuantumChemistrySoftware, ParserFactory>();

/**
 * Register a parser factory for a given software type.
 *
 * @param software - Quantum chemistry software identifier
 * @param factory - Factory function that creates the parser instance
 */
export function registerParser(software: QuantumChemistrySoftware, factory: ParserFactory): void {
  parserRegistry.set(software, factory);
}

// Register built-in parsers
registerParser('VASP', (content, filename) => new VASPParser(content, filename || 'INCAR'));
registerParser('Gaussian', content => new GaussianParser(content));
registerParser('ORCA', content => new ORCAParser(content));
registerParser('CP2K', content => new CP2KParser(content));
registerParser('Quantum ESPRESSO', content => new QEParser(content));
registerParser('GAMESS', content => new GAMESSParser(content));
registerParser('NWChem', content => new NWChemParser(content));

/**
 * Create appropriate parser for the given software and content
 *
 * Uses the parser registry to look up the factory function.
 * Throws if no parser has been registered for the given software.
 *
 * @param software - Quantum chemistry software identifier
 * @param content - File content to parse
 * @param filename - Optional filename hint
 * @returns Parser instance for the given software
 */
export function createParser(
  software: QuantumChemistrySoftware,
  content: string,
  filename?: string
): BaseParser {
  const factory = parserRegistry.get(software);
  if (!factory) {
    throw new Error(`Unsupported software: ${software}`);
  }
  return factory(content, filename);
}

/**
 * Parse input file and return structured result
 */
export function parseInput(software: QuantumChemistrySoftware, content: string, filename?: string) {
  const parser = createParser(software, content, filename);
  return parser.parseInput();
}

/**
 * Validate input file syntax
 */
export function validateInput(
  software: QuantumChemistrySoftware,
  content: string,
  filename?: string
) {
  const parser = createParser(software, content, filename);
  return parser.validate();
}
