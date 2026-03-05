/**
 * Electronic Structure Parameter Converter
 *
 * Extracts, converts, and applies parameters across quantum chemistry codes
 */

import * as vscode from 'vscode';
import * as fs from 'fs';
import * as path from 'path';
import { ParameterMapping, getParameterMappings, convertParameterValue } from './params';

/**
 * Extracted parameters from a quantum chemistry input file
 */
export interface ExtractedParameters {
  format: string;
  parameters: Record<string, any>;
  metadata: {
    filepath: string;
    lines: Record<string, number>; // parameter name -> line number
  };
}

/**
 * Converted parameters for target format
 */
export interface ConvertedParameters {
  format: string;
  parameters: Record<string, any>;
  warnings: string[];
  errors: string[];
  unmapped: string[]; // parameters that couldn't be mapped
}

/**
 * Parameter extraction result
 */
export interface ExtractionResult {
  success: boolean;
  data?: ExtractedParameters;
  error?: string;
  warnings: string[];
}

/**
 * Parameter conversion result
 */
export interface ConversionResult {
  success: boolean;
  source?: ExtractedParameters;
  target?: ConvertedParameters;
  error?: string;
  warnings: string[];
}

/**
 * Parameter Converter class
 */
export class ParameterConverter {
  /**
   * Extract parameters from VASP INCAR file
   */
  private extractVASP(filepath: string): ExtractionResult {
    const result: ExtractionResult = {
      success: false,
      warnings: [],
    };

    try {
      const content = fs.readFileSync(filepath, 'utf-8');
      const lines = content.split('\n');
      const parameters: Record<string, any> = {};
      const lineNumbers: Record<string, number> = {};

      lines.forEach((line, idx) => {
        // Skip comments and empty lines
        const trimmed = line.trim();
        if (!trimmed || trimmed.startsWith('#')) {
          return;
        }

        // Parse parameter = value
        const match = trimmed.match(/^(\w+)\s*=\s*(.+?)(?:\s*#.*)?$/);
        if (match) {
          const [, key, value] = match;
          parameters[key.toUpperCase()] = this.parseValue(value.trim());
          lineNumbers[key.toUpperCase()] = idx + 1;
        }
      });

      result.success = true;
      result.data = {
        format: 'vasp',
        parameters,
        metadata: {
          filepath,
          lines: lineNumbers,
        },
      };
    } catch (error: any) {
      result.error = `Failed to extract VASP parameters: ${error.message}`;
    }

    return result;
  }

  /**
   * Extract parameters from Quantum ESPRESSO input file
   */
  private extractQE(filepath: string): ExtractionResult {
    const result: ExtractionResult = {
      success: false,
      warnings: [],
    };

    try {
      const content = fs.readFileSync(filepath, 'utf-8');
      const parameters: Record<string, any> = {};
      const lineNumbers: Record<string, number> = {};
      const lines = content.split('\n');

      let currentSection = '';
      let lineIdx = 0;

      for (const line of lines) {
        lineIdx++;
        const trimmed = line.trim();

        // Detect section
        if (trimmed.startsWith('&')) {
          currentSection = trimmed.substring(1).toUpperCase();
          continue;
        }

        if (trimmed === '/') {
          currentSection = '';
          continue;
        }

        // Skip comments and empty lines
        if (!trimmed || trimmed.startsWith('!') || trimmed.startsWith('#')) {
          continue;
        }

        // Parse namelist parameters
        if (currentSection && trimmed.includes('=')) {
          const match = trimmed.match(/^(\w+)\s*=\s*(.+?)(?:\s*!.*)?$/);
          if (match) {
            const [, key, value] = match;
            const fullKey = `${currentSection}.${key.toLowerCase()}`;
            parameters[fullKey] = this.parseValue(value.trim());
            lineNumbers[fullKey] = lineIdx;
          }
        }
      }

      result.success = true;
      result.data = {
        format: 'qe',
        parameters,
        metadata: {
          filepath,
          lines: lineNumbers,
        },
      };
    } catch (error: any) {
      result.error = `Failed to extract QE parameters: ${error.message}`;
    }

    return result;
  }

  /**
   * Extract parameters from CP2K input file
   */
  private extractCP2K(filepath: string): ExtractionResult {
    const result: ExtractionResult = {
      success: false,
      warnings: [],
    };

    try {
      const content = fs.readFileSync(filepath, 'utf-8');
      const parameters: Record<string, any> = {};
      const lineNumbers: Record<string, number> = {};
      const lines = content.split('\n');

      let currentSection = '';
      let lineIdx = 0;

      for (const line of lines) {
        lineIdx++;
        const trimmed = line.trim();

        // Detect section
        if (trimmed.startsWith('&') && !trimmed.startsWith('&END')) {
          currentSection = trimmed.substring(1).toUpperCase();
          continue;
        }

        if (trimmed.startsWith('&END')) {
          currentSection = '';
          continue;
        }

        // Skip comments and empty lines
        if (!trimmed || trimmed.startsWith('#')) {
          continue;
        }

        // Parse keyword parameters
        if (currentSection && trimmed.includes(' ')) {
          const parts = trimmed.split(/\s+/);
          if (parts.length >= 2) {
            const key = parts[0].toUpperCase();
            const value = parts.slice(1).join(' ');
            const fullKey = `${currentSection}.${key}`;
            parameters[fullKey] = this.parseValue(value);
            lineNumbers[fullKey] = lineIdx;
          }
        }
      }

      result.success = true;
      result.data = {
        format: 'cp2k',
        parameters,
        metadata: {
          filepath,
          lines: lineNumbers,
        },
      };
    } catch (error: any) {
      result.error = `Failed to extract CP2K parameters: ${error.message}`;
    }

    return result;
  }

  /**
   * Extract parameters from Gaussian input file
   */
  private extractGaussian(filepath: string): ExtractionResult {
    const result: ExtractionResult = {
      success: false,
      warnings: [],
    };

    try {
      const content = fs.readFileSync(filepath, 'utf-8');
      const parameters: Record<string, any> = {};
      const lineNumbers: Record<string, number> = {};
      const lines = content.split('\n');

      let lineIdx = 0;
      let inRouteSection = false;

      for (const line of lines) {
        lineIdx++;
        const trimmed = line.trim();

        // Detect route section (starts with #)
        if (trimmed.startsWith('#')) {
          inRouteSection = true;
          // Parse route line
          const routeParams = this.parseGaussianRoute(trimmed);
          Object.entries(routeParams).forEach(([key, value]) => {
            parameters[key] = value;
            lineNumbers[key] = lineIdx;
          });
          continue;
        }

        // Empty line ends route section
        if (inRouteSection && trimmed === '') {
          inRouteSection = false;
          continue;
        }

        // Skip blank lines and charge/multiplicity line
        if (!trimmed || /^\d+\s+\d+$/.test(trimmed)) {
          continue;
        }
      }

      result.success = true;
      result.data = {
        format: 'gaussian',
        parameters,
        metadata: {
          filepath,
          lines: lineNumbers,
        },
      };
    } catch (error: any) {
      result.error = `Failed to extract Gaussian parameters: ${error.message}`;
    }

    return result;
  }

  /**
   * Parse Gaussian route section
   */
  private parseGaussianRoute(route: string): Record<string, any> {
    const params: Record<string, any> = {};

    // Remove leading #
    const cleanRoute = route.replace(/^#+\s*/, '');

    // Split by spaces
    const tokens = cleanRoute.split(/\s+/);

    tokens.forEach(token => {
      // Parse key=value pairs
      if (token.includes('=')) {
        const [key, value] = token.split('=');
        params[key.toUpperCase()] = value;
      } else if (token) {
        // Single keyword
        params[token.toUpperCase()] = true;
      }
    });

    return params;
  }

  /**
   * Parse a value string into appropriate type
   */
  private parseValue(value: string): any {
    // Boolean
    if (value.toUpperCase() === '.TRUE.' || value.toUpperCase() === 'TRUE') {
      return true;
    }
    if (value.toUpperCase() === '.FALSE.' || value.toUpperCase() === 'FALSE') {
      return false;
    }

    // Number
    const num = parseFloat(value);
    if (!isNaN(num)) {
      return num;
    }

    // String (remove quotes if present)
    return value.replace(/^["']|["']$/g, '');
  }

  /**
   * Extract parameters from file (auto-detect format)
   */
  public extractParameters(filepath: string): ExtractionResult {
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
    } else if (ext === '.gjf' || ext === '.com') {
      format = 'gaussian';
    }

    // Extract based on format
    switch (format) {
      case 'vasp':
        return this.extractVASP(filepath);
      case 'qe':
        return this.extractQE(filepath);
      case 'cp2k':
        return this.extractCP2K(filepath);
      case 'gaussian':
        return this.extractGaussian(filepath);
      default:
        return {
          success: false,
          error: `Unsupported format: ${format}`,
          warnings: [],
        };
    }
  }

  /**
   * Get the base parameter name (without section prefix for QE/CP2K)
   */
  private getBaseParamName(key: string, format: string): string {
    if (format === 'qe' || format === 'cp2k') {
      // For SECTION.param format, return just 'param'
      const parts = key.split('.');
      return parts.length > 1 ? parts[parts.length - 1] : key;
    }
    return key;
  }

  /**
   * Convert parameters from source to target format
   */
  public convertParameters(
    sourceParams: ExtractedParameters,
    targetFormat: string
  ): ConvertedParameters {
    const result: ConvertedParameters = {
      format: targetFormat,
      parameters: {},
      warnings: [],
      errors: [],
      unmapped: [],
    };

    // Get mappings
    const mappings = getParameterMappings(sourceParams.format, targetFormat);

    // Create mapping lookup (handle both full and base parameter names)
    const mappingLookup = new Map<string, ParameterMapping>();
    mappings.forEach(m => {
      // Store both uppercase versions
      mappingLookup.set(m.sourceParam.toUpperCase(), m);
    });

    // Convert each parameter
    Object.entries(sourceParams.parameters).forEach(([key, value]) => {
      // Try full key first
      let mapping = mappingLookup.get(key.toUpperCase());

      // If not found, try base parameter name (for QE/CP2K)
      if (!mapping) {
        const baseKey = this.getBaseParamName(key, sourceParams.format);
        mapping = mappingLookup.get(baseKey.toUpperCase());
      }

      if (mapping) {
        try {
          const convertedValue = convertParameterValue(mapping, value);
          result.parameters[mapping.targetParam] = convertedValue;
        } catch (error: any) {
          result.warnings.push(`Failed to convert ${key}: ${error.message}`);
        }
      } else {
        result.unmapped.push(key);
      }
    });

    // Add warnings for unmapped parameters
    if (result.unmapped.length > 0) {
      result.warnings.push(
        `${result.unmapped.length} parameters could not be mapped: ${result.unmapped
          .slice(0, 5)
          .join(', ')}${result.unmapped.length > 5 ? '...' : ''}`
      );
    }

    return result;
  }

  /**
   * Full conversion workflow
   */
  public async convertFile(sourcePath: string, targetFormat: string): Promise<ConversionResult> {
    const result: ConversionResult = {
      success: false,
      warnings: [],
    };

    // Extract source parameters
    const extraction = this.extractParameters(sourcePath);

    if (!extraction.success || !extraction.data) {
      result.error = extraction.error || 'Failed to extract parameters';
      result.warnings = extraction.warnings;
      return result;
    }

    result.source = extraction.data;
    result.warnings.push(...extraction.warnings);

    // Convert parameters
    const conversion = this.convertParameters(extraction.data, targetFormat);

    result.target = conversion;
    result.warnings.push(...conversion.warnings);
    result.success = true;

    return result;
  }

  /**
   * Get parameter mapping suggestions
   */
  public getMappingSuggestions(sourceFormat: string, targetFormat: string): ParameterMapping[] {
    return getParameterMappings(sourceFormat, targetFormat);
  }
}
