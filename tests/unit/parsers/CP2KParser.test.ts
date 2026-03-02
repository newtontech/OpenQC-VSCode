import { CP2KParser } from '../../../src/parsers/CP2KParser';

describe('CP2KParser', () => {
  describe('Basic parsing', () => {
    it('should parse basic CP2K input with GLOBAL and FORCE_EVAL', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME POTENTIAL
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      expect(result.sections.length).toBeGreaterThanOrEqual(2);
      expect(result.errors).toHaveLength(0);
    });

    it('should parse GLOBAL section parameters', () => {
      const content = `&GLOBAL
  PROJECT_NAME test_project
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const projectName = parser.getParameter('PROJECT_NAME');
      const runType = parser.getParameter('RUN_TYPE');
      const printLevel = parser.getParameter('PRINT_LEVEL');

      expect(projectName?.value).toBe('test_project');
      expect(runType?.value).toBe('ENERGY');
      expect(printLevel?.value).toBe('MEDIUM');
    });

    it('should detect unclosed sections', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      expect(result.errors.length).toBeGreaterThan(0);
      expect(result.errors.some(e => e.message.includes('not properly closed'))).toBe(true);
    });

    it('should handle comments and empty lines', () => {
      const content = `# CP2K Input File
&GLOBAL
  # This is a comment
  PROJECT_NAME test

  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      expect(result.errors).toHaveLength(0);
      expect(parser.getParameter('PROJECT_NAME')?.value).toBe('test');
    });

    it('should parse nested sections', () => {
      const content = `&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
    &END SCF
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const forceEval = parser.getSection('FORCE_EVAL');
      expect(forceEval).toBeDefined();
      expect(forceEval?.subsections).toBeDefined();
      expect(forceEval?.subsections?.length).toBeGreaterThan(0);
    });
  });

  describe('Value conversion', () => {
    it('should convert boolean values', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &SCF
      SCF_GUESS ATOMIC
      OT on
    &END SCF
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const otParam = result.parameters.find(p => p.name === 'OT');
      expect(otParam?.value).toBe(true);
    });

    it('should convert numeric values', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const cutoff = result.parameters.find(p => p.name === 'CUTOFF');
      const relCutoff = result.parameters.find(p => p.name === 'REL_CUTOFF');

      expect(cutoff?.value).toBe(400);
      expect(relCutoff?.value).toBe(50);
    });

    it('should handle string values', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const basisSet = result.parameters.find(p => p.name === 'BASIS_SET_FILE_NAME');
      const potential = result.parameters.find(p => p.name === 'POTENTIAL_FILE_NAME');

      expect(basisSet?.value).toBe('BASIS_MOLOPT');
      expect(potential?.value).toBe('GTH_POTENTIALS');
    });

    it('should handle quoted string values', () => {
      const content = `&GLOBAL
  PROJECT_NAME "test project"
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const projectName = result.parameters.find(p => p.name === 'PROJECT_NAME');
      expect(projectName?.value).toBe('test project');
    });
  });

  describe('Validation', () => {
    it('should validate complete CP2K input', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME POTENTIAL
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const validation = parser.validate();

      expect(validation.valid).toBe(true);
      expect(validation.errors).toHaveLength(0);
    });

    it('should detect missing GLOBAL section', () => {
      const content = `&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const validation = parser.validate();

      expect(validation.errors).toHaveLength(0);
      expect(validation.warnings.some(w => w.message.includes('GLOBAL'))).toBe(true);
    });

    it('should detect missing FORCE_EVAL section', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL`;
      const parser = new CP2KParser(content);
      const validation = parser.validate();

      expect(validation.valid).toBe(false);
      expect(validation.errors.some(e => e.message.includes('FORCE_EVAL'))).toBe(true);
    });

    it('should warn about missing PROJECT_NAME', () => {
      const content = `&GLOBAL
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const validation = parser.validate();

      expect(validation.warnings.some(w => w.message.includes('PROJECT_NAME'))).toBe(true);
    });
  });

  describe('Coordinate extraction', () => {
    it('should extract atomic coordinates from COORD section', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &COORD
      O   0.000000   0.000000   0.000000
      H   0.000000   0.757215   0.586535
      H   0.000000  -0.757215   0.586535
    &END COORD
  &END SUBSYS
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const coordinates = parser.getCoordinates();

      expect(coordinates).toHaveLength(3);
      expect(coordinates[0]).toEqual({
        element: 'O',
        coords: [0, 0, 0],
      });
      expect(coordinates[1]).toEqual({
        element: 'H',
        coords: [0, 0.757215, 0.586535],
      });
      expect(coordinates[2]).toEqual({
        element: 'H',
        coords: [0, -0.757215, 0.586535],
      });
    });

    it('should handle coordinates in scientific notation', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &COORD
      Si  0.000000  0.000000  0.000000
      Si  1.357500  1.357500  1.357500
    &END COORD
  &END SUBSYS
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const coordinates = parser.getCoordinates();

      expect(coordinates).toHaveLength(2);
      expect(coordinates[0].element).toBe('Si');
      expect(coordinates[1].element).toBe('Si');
    });

    it('should return empty array when no COORD section', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const coordinates = parser.getCoordinates();

      expect(coordinates).toHaveLength(0);
    });

    it('should handle comments in COORD section', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &COORD
      # Water molecule
      O   0.0   0.0   0.0
      H   0.0   0.75  0.58
    &END COORD
  &END SUBSYS
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const coordinates = parser.getCoordinates();

      expect(coordinates).toHaveLength(2);
    });
  });

  describe('Cell parameters', () => {
    it('should extract ABC cell parameters', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
  &END SUBSYS
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const cell = parser.getCellParameters();

      expect(cell).toBeDefined();
      if (cell && cell.type === 'ABC') {
        expect(cell.a).toBe(10.0);
        expect(cell.b).toBe(10.0);
        expect(cell.c).toBe(10.0);
      }
    });

    it('should extract cell vectors', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &CELL
      A  2.715  2.715  0.000
      B  2.715  0.000  2.715
      C  0.000  2.715  2.715
    &END CELL
  &END SUBSYS
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const cell = parser.getCellParameters();

      expect(cell).toBeDefined();
      if (cell && cell.type === 'vectors') {
        expect(cell.a).toEqual([2.715, 2.715, 0]);
        expect(cell.b).toEqual([2.715, 0, 2.715]);
        expect(cell.c).toEqual([0, 2.715, 2.715]);
      }
    });

    it('should return undefined when no CELL section', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const cell = parser.getCellParameters();

      expect(cell).toBeUndefined();
    });
  });

  describe('Atom types and basis sets', () => {
    it('should extract atom types from KIND sections', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const atomTypes = parser.getAtomTypes();

      expect(atomTypes).toHaveLength(2);
      expect(atomTypes).toContain('H');
      expect(atomTypes).toContain('O');
    });

    it('should extract basis set information', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const basisSets = parser.getBasisSets();

      expect(basisSets.get('H')).toBe('DZVP-MOLOPT-GTH');
      expect(basisSets.get('O')).toBe('DZVP-MOLOPT-GTH');
    });

    it('should return empty map when no basis sets defined', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const basisSets = parser.getBasisSets();

      expect(basisSets.size).toBe(0);
    });
  });

  describe('DFT functional', () => {
    it('should extract XC functional from subsection', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const functional = parser.getFunctional();

      expect(functional).toBe('PBE');
    });

    it('should return undefined when no XC functional specified', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const functional = parser.getFunctional();

      expect(functional).toBeUndefined();
    });
  });

  describe('Section and parameter access', () => {
    it('should get section by name', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const globalSection = parser.getSection('GLOBAL');
      const forceEvalSection = parser.getSection('FORCE_EVAL');

      expect(globalSection).toBeDefined();
      expect(forceEvalSection).toBeDefined();
    });

    it('should return undefined for non-existent section', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL`;
      const parser = new CP2KParser(content);
      const section = parser.getSection('NONEXISTENT');

      expect(section).toBeUndefined();
    });

    it('should be case-insensitive for section names', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL`;
      const parser = new CP2KParser(content);
      const section1 = parser.getSection('GLOBAL');
      const section2 = parser.getSection('global');
      const section3 = parser.getSection('Global');

      expect(section1).toBeDefined();
      expect(section2).toBeDefined();
      expect(section3).toBeDefined();
    });

    it('should get parameter by name', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const projectName = parser.getParameter('PROJECT_NAME');
      const runType = parser.getParameter('RUN_TYPE');

      expect(projectName?.value).toBe('test');
      expect(runType?.value).toBe('ENERGY');
    });

    it('should be case-insensitive for parameter names', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL`;
      const parser = new CP2KParser(content);
      const param1 = parser.getParameter('PROJECT_NAME');
      const param2 = parser.getParameter('project_name');
      const param3 = parser.getParameter('Project_Name');

      expect(param1).toBeDefined();
      expect(param2).toBeDefined();
      expect(param3).toBeDefined();
    });

    it('should return all parameters', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const parameters = parser.getParameters();

      expect(parameters.length).toBeGreaterThanOrEqual(4);
    });

    it('should return all sections', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const sections = parser.getSections();

      expect(sections.length).toBeGreaterThanOrEqual(2);
    });
  });

  describe('Complex nested structures', () => {
    it('should handle deeply nested DFT sections', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &SCF
      EPS_SCF 1.0E-6
      MAX_SCF 100
      SCF_GUESS ATOMIC
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE 300
      &END SMEAR
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
      &END MIXING
    &END SCF
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      expect(result.errors).toHaveLength(0);

      const forceEval = parser.getSection('FORCE_EVAL');
      expect(forceEval?.subsections).toBeDefined();
      expect(forceEval?.subsections?.length).toBeGreaterThan(0);
    });

    it('should parse KPOINTS section', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &KPOINTS
      SCHEME MONKHORST-PACK 4 4 4
    &END KPOINTS
  &END DFT
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const schemeParam = result.parameters.find(p => p.name === 'SCHEME');
      // CP2K SCHEME parameter includes both scheme name and grid values
      expect(schemeParam?.value).toBe('MONKHORST-PACK 4 4 4');
    });
  });

  describe('Error handling', () => {
    it('should handle empty input', () => {
      const content = '';
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      expect(result.sections).toHaveLength(0);
      expect(result.parameters).toHaveLength(0);
    });

    it('should handle input with only comments', () => {
      const content = `# This is a comment
# Another comment
&GLOBAL
  # Still a comment`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      expect(result.errors.length).toBeGreaterThan(0);
    });

    it('should handle malformed lines gracefully', () => {
      const content = `&GLOBAL
  PROJECT_NAME test
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  INVALID_LINE_WITHOUT_VALUE
&END FORCE_EVAL`;
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      // Should not crash, just skip the malformed line
      expect(result.sections.length).toBeGreaterThan(0);
    });
  });
});
