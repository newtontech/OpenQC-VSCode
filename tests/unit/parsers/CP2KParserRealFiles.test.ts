import { CP2KParser } from '../../../src/parsers/CP2KParser';
import * as fs from 'fs';
import * as path from 'path';

describe('CP2KParser with Real Files', () => {
  // Navigate from tests/unit/parsers to tests/fixtures
  const fixturesDir = path.join(__dirname, '../../fixtures/cp2k');

  describe('H2O.inp - Water molecule geometry optimization', () => {
    let content: string;
    let parser: CP2KParser;

    beforeAll(() => {
      const filePath = path.join(fixturesDir, 'H2O.inp');
      content = fs.readFileSync(filePath, 'utf8');
      parser = new CP2KParser(content);
    });

    it('should parse without errors', () => {
      const result = parser.parseInput();
      expect(result.errors).toHaveLength(0);
    });

    it('should extract project name', () => {
      const projectParam = parser.getParameter('PROJECT');
      expect(projectParam).toBeDefined();
      expect(projectParam?.value).toBe('H2O_Optimization');
    });

    it('should detect geometry optimization run type', () => {
      const runType = parser.getParameter('RUN_TYPE');
      expect(runType?.value).toBe('GEO_OPT');
    });

    it('should extract atomic coordinates', () => {
      const coordinates = parser.getCoordinates();
      expect(coordinates.length).toBe(3);
      expect(coordinates[0]).toEqual({
        element: 'O',
        coords: [0, 0, 0],
      });
      expect(coordinates[1].element).toBe('H');
      expect(coordinates[2].element).toBe('H');
    });

    it('should extract cell parameters', () => {
      const cell = parser.getCellParameters();
      expect(cell).toBeDefined();
      if (cell && cell.type === 'ABC') {
        expect(cell.a).toBe(10.0);
        expect(cell.b).toBe(10.0);
        expect(cell.c).toBe(10.0);
      }
    });

    it('should extract atom types', () => {
      const atomTypes = parser.getAtomTypes();
      expect(atomTypes).toContain('H');
      expect(atomTypes).toContain('O');
    });

    it('should extract basis sets', () => {
      const basisSets = parser.getBasisSets();
      expect(basisSets.get('H')).toBeDefined();
      expect(basisSets.get('O')).toBeDefined();
    });

    it('should pass validation', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
    });
  });

  describe('Si_bulk.inp - Silicon bulk calculation', () => {
    let content: string;
    let parser: CP2KParser;

    beforeAll(() => {
      const filePath = path.join(fixturesDir, 'Si_bulk.inp');
      content = fs.readFileSync(filePath, 'utf8');
      parser = new CP2KParser(content);
    });

    it('should parse without errors', () => {
      const result = parser.parseInput();
      expect(result.errors).toHaveLength(0);
    });

    it('should extract project name', () => {
      const projectParam = parser.getParameter('PROJECT');
      expect(projectParam).toBeDefined();
      expect(projectParam?.value).toBe('Si_bulk');
    });

    it('should detect energy run type', () => {
      const runType = parser.getParameter('RUN_TYPE');
      expect(runType?.value).toBe('ENERGY');
    });

    it('should extract atomic coordinates', () => {
      const coordinates = parser.getCoordinates();
      expect(coordinates.length).toBe(2);
      expect(coordinates.every(c => c.element === 'Si')).toBe(true);
    });

    it('should extract cell vectors', () => {
      const cell = parser.getCellParameters();
      expect(cell).toBeDefined();
      if (cell && cell.type === 'vectors') {
        expect(cell.a.length).toBe(3);
        expect(cell.b.length).toBe(3);
        expect(cell.c.length).toBe(3);
      }
    });

    it('should extract atom types', () => {
      const atomTypes = parser.getAtomTypes();
      expect(atomTypes).toContain('Si');
    });

    it('should pass validation', () => {
      const validation = parser.validate();
      expect(validation.valid).toBe(true);
    });
  });

  describe('MGRID parameters', () => {
    it('should extract MGRID cutoff from H2O file', () => {
      const filePath = path.join(fixturesDir, 'H2O.inp');
      const content = fs.readFileSync(filePath, 'utf8');
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const cutoffParam = result.parameters.find(p => p.name === 'CUTOFF');
      expect(cutoffParam).toBeDefined();
      expect(cutoffParam?.value).toBe(400);
    });

    it('should extract MGRID cutoff from Si_bulk file', () => {
      const filePath = path.join(fixturesDir, 'Si_bulk.inp');
      const content = fs.readFileSync(filePath, 'utf8');
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const cutoffParam = result.parameters.find(p => p.name === 'CUTOFF');
      expect(cutoffParam).toBeDefined();
      expect(cutoffParam?.value).toBe(500);
    });
  });

  describe('SCF parameters', () => {
    it('should extract SCF parameters from H2O file', () => {
      const filePath = path.join(fixturesDir, 'H2O.inp');
      const content = fs.readFileSync(filePath, 'utf8');
      const parser = new CP2KParser(content);
      const result = parser.parseInput();

      const epsScf = result.parameters.find(p => p.name === 'EPS_SCF');
      const maxScf = result.parameters.find(p => p.name === 'MAX_SCF');
      const scfGuess = result.parameters.find(p => p.name === 'SCF_GUESS');

      expect(epsScf).toBeDefined();
      expect(epsScf?.value).toBe(1.0E-6);
      expect(maxScf?.value).toBe(100);
      expect(scfGuess?.value).toBe('ATOMIC');
    });
  });

  describe('Complex nested structure', () => {
    it('should handle SMEAR section in Si_bulk', () => {
      const filePath = path.join(fixturesDir, 'Si_bulk.inp');
      const content = fs.readFileSync(filePath, 'utf8');
      const parser = new CP2KParser(content);

      const forceEval = parser.getSection('FORCE_EVAL');
      expect(forceEval?.subsections).toBeDefined();
      expect(forceEval?.subsections?.length).toBeGreaterThan(0);
    });
  });
});
