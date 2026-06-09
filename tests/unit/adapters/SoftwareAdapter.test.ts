/**
 * Tests for the software adapter registry (issue #79).
 * @module tests/unit/adapters/SoftwareAdapter.test
 */

import { adapterRegistry } from '../../../src/adapters/SoftwareAdapter';

describe('SoftwareAdapter Registry', () => {
  describe('existing software coverage', () => {
    const existingIds = ['cp2k', 'vasp', 'gaussian', 'orca', 'qe', 'gamess', 'nwchem'];

    for (const id of existingIds) {
      it(`has adapter for existing software: ${id}`, () => {
        const adapter = adapterRegistry.get(id);
        expect(adapter).toBeDefined();
        expect(adapter!.id).toBe(id);
        expect(adapter!.displayName).toBeTruthy();
        expect(adapter!.filePatterns.length).toBeGreaterThan(0);
      });
    }
  });

  describe('coming soon software (#79)', () => {
    const comingSoonIds = [
      'psi4',
      'molpro',
      'openmolcas',
      'dalton',
      'turbomole',
      'castep',
      'abinit',
      'crystal',
    ];

    for (const id of comingSoonIds) {
      it(`has adapter for: ${id}`, () => {
        const adapter = adapterRegistry.get(id);
        expect(adapter).toBeDefined();
        expect(adapter!.id).toBe(id);
        expect(adapter!.displayName).toBeTruthy();
        expect(adapter!.category).toBeDefined();
        expect(adapter!.filePatterns.length).toBeGreaterThan(0);
      });
    }
  });

  describe('category coverage', () => {
    it('has molecular-qc adapters', () => {
      const adapters = adapterRegistry.getByCategory('molecular-qc');
      expect(adapters.length).toBeGreaterThanOrEqual(5);
    });

    it('has periodic-dft adapters', () => {
      const adapters = adapterRegistry.getByCategory('periodic-dft');
      expect(adapters.length).toBeGreaterThanOrEqual(3);
    });

    it('has multireference adapters', () => {
      const adapters = adapterRegistry.getByCategory('multireference');
      expect(adapters.length).toBeGreaterThanOrEqual(1);
    });

    it('has properties adapters', () => {
      const adapters = adapterRegistry.getByCategory('properties');
      expect(adapters.length).toBeGreaterThanOrEqual(1);
    });
  });

  describe('detection', () => {
    it('detects VASP POSCAR from filename', () => {
      const adapter = adapterRegistry.detect('', '/path/to/POSCAR');
      expect(adapter?.id).toBe('vasp');
    });

    it('detects VASP CONTCAR from filename', () => {
      const adapter = adapterRegistry.detect('', 'CONTCAR');
      expect(adapter?.id).toBe('vasp');
    });

    it('detects CP2K from content', () => {
      const content = '&GLOBAL\n  PROJECT_NAME test\n&END GLOBAL\n&FORCE_EVAL\n&END FORCE_EVAL';
      const adapter = adapterRegistry.detect(content, 'input.inp');
      expect(adapter?.id).toBe('cp2k');
    });

    it('detects Gaussian from content', () => {
      const content = '%chk=test.chk\n# B3LYP/6-31G(d) Opt\n\nTitle\n\n0 1\nH 0 0 0\nH 0 0 0.74\n';
      const adapter = adapterRegistry.detect(content, 'test.gjf');
      expect(adapter?.id).toBe('gaussian');
    });

    it('detects ORCA from content', () => {
      const content = '! B3LYP def2-SVP Opt\n\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*';
      const adapter = adapterRegistry.detect(content, 'test.inp');
      expect(adapter?.id).toBe('orca');
    });

    it('detects Psi4 from content', () => {
      const content =
        'import psi4\n\nmolecule {\n  H 0 0 0\n  H 0 0 0.74\n}\n\nset { basis cc-pvdz }\nenergy("scf")';
      const adapter = adapterRegistry.detect(content, 'psi4.in');
      expect(adapter?.id).toBe('psi4');
    });

    it('detects CASTEP from content', () => {
      const content =
        '%block lattice_cart\n  5.43 0.00 0.00\n%endblock lattice_cart\n%block positions_abs\n  Si 0.0 0.0 0.0\n%endblock positions_abs';
      const adapter = adapterRegistry.detect(content, 'Si.cell');
      expect(adapter?.id).toBe('castep');
    });

    it('detects Crystal from content', () => {
      const content = 'CRYSTAL\n0 0 0\n5.43\n2\n0.0 0.0 0.0\n0.25 0.25 0.25\nEND';
      const adapter = adapterRegistry.detect(content, 'Si.d12');
      expect(adapter?.id).toBe('crystal');
    });

    it('detects Abinit from content', () => {
      const content = 'ndtset 1\necut 30\nnband 10\n';
      const adapter = adapterRegistry.detect(content, 'run.abi');
      expect(adapter?.id).toBe('abinit');
    });

    it('returns undefined for unknown content', () => {
      const adapter = adapterRegistry.detect('random text file content', 'readme.txt');
      expect(adapter).toBeUndefined();
    });
  });

  describe('adapter properties', () => {
    it('cclib-supported adapters have cclibSupport=true', () => {
      const cclibAdapters = adapterRegistry.getAll().filter(a => a.cclibSupport);
      expect(cclibAdapters.length).toBeGreaterThanOrEqual(5);
    });

    it('all adapters have display names', () => {
      for (const adapter of adapterRegistry.getAll()) {
        expect(adapter.displayName.length).toBeGreaterThan(0);
      }
    });

    it('total adapter count covers all software', () => {
      expect(adapterRegistry.getAll().length).toBeGreaterThanOrEqual(15);
    });
  });
});
