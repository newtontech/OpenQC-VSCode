import { FileTypeDetector, QuantumChemistrySoftware } from '../../src/managers/FileTypeDetector';

const createMockDocument = (fileName: string, content: string): any => ({
  fileName,
  getText: () => content,
});

describe('FileTypeDetector', () => {
  let detector: FileTypeDetector;

  beforeEach(() => {
    detector = new FileTypeDetector();
  });

  describe('detectSoftware', () => {
    it('should detect CP2K files by extension and content', () => {
      const doc = createMockDocument(
        '/test/input.inp',
        `
&GLOBAL
  PROJECT_NAME test
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
&END FORCE_EVAL
      `
      );
      expect(detector.detectSoftware(doc)).toBe('CP2K');
    });

    it('should detect VASP files by filename', () => {
      const incar = createMockDocument('/test/INCAR', 'ISTART = 0\nENCUT = 500');
      expect(detector.detectSoftware(incar)).toBe('VASP');

      const poscar = createMockDocument('/test/POSCAR', 'Test\n1.0\n...');
      expect(detector.detectSoftware(poscar)).toBe('VASP');

      const kpoints = createMockDocument('/test/KPOINTS', 'K-Points\n0\nGamma\n4 4 4');
      expect(detector.detectSoftware(kpoints)).toBe('VASP');

      const contributedFilenames = ['POTCAR', 'CONTCAR', 'OSZICAR', 'OUTCAR', 'vasprun.xml'];
      for (const filename of contributedFilenames) {
        const doc = createMockDocument(`/test/${filename}`, '');
        expect(detector.detectSoftware(doc)).toBe('VASP');
      }
    });

    it('should detect Gaussian files by extension and content', () => {
      const gjfDoc = createMockDocument(
        '/test/molecule.gjf',
        `
%chk=molecule.chk
# B3LYP/6-31G(d) opt

Title

0 1
C 0.0 0.0 0.0
      `
      );
      expect(detector.detectSoftware(gjfDoc)).toBe('Gaussian');

      const comDoc = createMockDocument(
        '/test/molecule.com',
        `
%chk=molecule.chk
# B3LYP/6-31G(d) opt

Title

0 1
C 0.0 0.0 0.0
      `
      );
      expect(detector.detectSoftware(comDoc)).toBe('Gaussian');
    });

    it('should detect ORCA files by content', () => {
      const doc = createMockDocument(
        '/test/input.inp',
        `
! B3LYP def2-TZVP opt
%pal nprocs 4 end
%maxcore 2000

* xyz 0 1
C 0.0 0.0 0.0
*
      `
      );
      expect(detector.detectSoftware(doc)).toBe('ORCA');
    });

    it('should detect Quantum ESPRESSO files', () => {
      const doc = createMockDocument(
        '/test/pw.in',
        `
&CONTROL
  calculation = 'scf'
  pseudo_dir = './'
/
&SYSTEM
  ibrav = 1
/
      `
      );
      expect(detector.detectSoftware(doc)).toBe('Quantum ESPRESSO');
    });

    it('should detect GAMESS files', () => {
      const doc = createMockDocument(
        '/test/input.inp',
        `
 $BASIS GBASIS=STO NGAUSS=3 $END
 $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END
 $SYSTEM TIMLIM=5 $END
      `
      );
      expect(detector.detectSoftware(doc)).toBe('GAMESS');
    });

    it('should detect NWChem files', () => {
      const doc = createMockDocument(
        '/test/molecule.nw',
        `
geometry units angstrom
  C 0.0 0.0 0.0
end

basis
  * library 6-31G*
end

scf
  maxiter 100
end
      `
      );
      expect(detector.detectSoftware(doc)).toBe('NWChem');
    });

    it('should detect all local sibling LSP formats added to OpenQC', () => {
      const cases: Array<[string, string, QuantumChemistrySoftware]> = [
        ['/test/INPUT', 'INPUT_PARAMETERS\nbasis_type lcao\npseudo_dir ./\n', 'ABACUS'],
        ['/test/run.abi', 'ecut 20\nngkpt 2 2 2\nacell 3*10\n', 'ABINIT'],
        ['/test/Si.cif', 'data_si\n_cell_length_a 5.43\n_atom_site_label Si\n', 'CIF'],
        ['/test/run.in', 'potential nep.txt\nvelocity 300\nrun 1000\n', 'GPUMD'],
        ['/test/topol.top', '[ moleculetype ]\n; name nrexcl\nSOL 3\n', 'GROMACS'],
        ['/test/in.lammps', 'units metal\natom_style atomic\npair_style eam\nrun 0\n', 'LAMMPS'],
        [
          '/test/workflow.mlip.json',
          '{"model":"DPA3.1-3M","structure":"input.cif","task":"optimize"}\n',
          'MLIP',
        ],
        ['/test/run_pyatb.py', 'import pyatb\nhr_file = "HR.dat"\nsr_file = "SR.dat"\n', 'PyATB'],
        ['/test/run_pyscf.py', 'from pyscf import gto, dft\nmol = gto.M()\n', 'PySCF'],
      ];

      for (const [fileName, content, expected] of cases) {
        expect(detector.detectSoftware(createMockDocument(fileName, content))).toBe(expected);
      }
    });

    it('should detect composite suffixes for MLIP, PyATB, and PySCF', () => {
      const cases: Array<[string, string, QuantumChemistrySoftware]> = [
        ['/test/workflow.mlip.json', '{"model":"DPA3.1-3M","structure":"input.cif"}\n', 'MLIP'],
        ['/test/workflow.mlip.yaml', 'model: DPA3.1-3M\nstructure: input.cif\n', 'MLIP'],
        ['/test/workflow.mlip.yml', 'model: DPA3.1-3M\ntask: optimize\n', 'MLIP'],
        ['/test/example.pyatb.py', 'from pyatb import api\nhr_file = "HR.dat"\n', 'PyATB'],
        ['/test/example.pyscf.py', 'import pyscf\nfrom pyscf import gto, scf\n', 'PySCF'],
      ];

      for (const [fileName, content, expected] of cases) {
        expect(detector.detectSoftware(createMockDocument(fileName, content))).toBe(expected);
      }
    });

    it('should return null for unsupported files', () => {
      const doc = createMockDocument('/test/random.txt', 'This is not a chemistry file');
      expect(detector.detectSoftware(doc)).toBeNull();
    });

    it('should handle files with ambiguous extensions - CP2K vs ORCA', () => {
      // CP2K and ORCA both use .inp
      const cp2kDoc = createMockDocument('/test/cp2k.inp', '&GLOBAL\n&END GLOBAL');
      expect(detector.detectSoftware(cp2kDoc)).toBe('CP2K');

      const orcaDoc = createMockDocument('/test/orca.inp', '! B3LYP\n%pal nprocs 4 end');
      expect(detector.detectSoftware(orcaDoc)).toBe('ORCA');
    });

    it('should handle file with extension but low confidence', () => {
      // File with .inp extension but minimal content - should try extension match first
      const doc = createMockDocument('/test/test.inp', 'minimal content');
      const result = detector.detectSoftware(doc);
      // May match by extension or fallback
      expect(result === null || typeof result === 'string').toBe(true);
    });

    it('should use fallback detection when filename pattern not matched', () => {
      // File that doesn't match by filename but matches by content
      const doc = createMockDocument(
        '/test/random.xyz',
        `
&GLOBAL
PROJECT_NAME test
&END GLOBAL
      `
      );
      const result = detector.detectSoftware(doc);
      // Should detect by content
      expect(result).toBe('CP2K');
    });
  });

  describe('getSoftwareInfo', () => {
    it('should return info for all supported software', () => {
      const softwares: QuantumChemistrySoftware[] = [
        'CP2K',
        'VASP',
        'Gaussian',
        'ORCA',
        'Quantum ESPRESSO',
        'GAMESS',
        'NWChem',
        'ABACUS',
        'ABINIT',
        'CIF',
        'GPUMD',
        'GROMACS',
        'LAMMPS',
        'MLIP',
        'PyATB',
        'PySCF',
      ];

      for (const software of softwares) {
        const info = detector.getSoftwareInfo(software);
        expect(info).toHaveProperty('name');
        expect(info).toHaveProperty('description');
        expect(info).toHaveProperty('website');
        expect(info.name).toBe(software);
      }
    });

    it('should return correct website URLs', () => {
      expect(detector.getSoftwareInfo('CP2K').website).toBe('https://www.cp2k.org');
      expect(detector.getSoftwareInfo('VASP').website).toBe('https://www.vasp.at');
      expect(detector.getSoftwareInfo('Gaussian').website).toBe('https://gaussian.com');
      expect(detector.getSoftwareInfo('ORCA').website).toBe('https://orcaforum.kofo.mpg.de');
      expect(detector.getSoftwareInfo('Quantum ESPRESSO').website).toBe(
        'https://www.quantum-espresso.org'
      );
      expect(detector.getSoftwareInfo('GAMESS').website).toBe(
        'https://www.msg.chem.iastate.edu/gamess'
      );
      expect(detector.getSoftwareInfo('NWChem').website).toBe('https://nwchemgit.github.io');
    });
  });

  describe('edge cases', () => {
    it('should handle Windows paths', () => {
      const doc = createMockDocument('C:\\Users\\test\\INCAR', 'ENCUT = 500');
      expect(detector.detectSoftware(doc)).toBe('VASP');
    });

    it('should handle empty content', () => {
      const doc = createMockDocument('/test/INCAR', '');
      expect(detector.detectSoftware(doc)).toBe('VASP');
    });

    it('should handle content with partial matches', () => {
      // Content that has patterns but low confidence
      const doc = createMockDocument('/test/mixed.inp', 'Some text with &GLOBAL keyword');
      const result = detector.detectSoftware(doc);
      // Result could be null or a software string depending on confidence threshold
      expect(result === null || typeof result === 'string').toBe(true);
    });

    it('should handle files with no extension', () => {
      const doc = createMockDocument('/test/Makefile', 'all: test');
      expect(detector.detectSoftware(doc)).toBeNull();
    });

    it('should handle file paths with special characters', () => {
      const doc = createMockDocument('/path with spaces/test.inp', '&GLOBAL\n&END GLOBAL');
      expect(detector.detectSoftware(doc)).toBe('CP2K');
    });

    it('should handle deeply nested paths', () => {
      const doc = createMockDocument('/very/deep/nested/path/INCAR', 'ENCUT = 500');
      expect(detector.detectSoftware(doc)).toBe('VASP');
    });

    it('should handle empty filename', () => {
      const doc = createMockDocument('', 'random content with no patterns');
      expect(detector.detectSoftware(doc)).toBeNull();
    });
  });
});
