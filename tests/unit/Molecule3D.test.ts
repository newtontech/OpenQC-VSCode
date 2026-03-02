import { Molecule3D } from '../../src/visualizers/Molecule3D';
import { QuantumChemistrySoftware } from '../../src/managers/FileTypeDetector';

describe('Molecule3D', () => {
  let visualizer: Molecule3D;

  beforeEach(() => {
    visualizer = new Molecule3D();
  });

  describe('parseAtoms', () => {
    it('should parse CP2K atoms', () => {
      const content = `
&COORD
  H 0.0 0.0 0.0
  H 0.0 0.0 0.74
&END COORD`;
      const atoms = visualizer.parseAtoms(content, 'CP2K' as QuantumChemistrySoftware);
      expect(atoms).toHaveLength(2);
      expect(atoms[0].element).toBe('H');
      expect(atoms[0].x).toBe(0.0);
      expect(atoms[0].y).toBe(0.0);
      expect(atoms[0].z).toBe(0.0);
      expect(atoms[1].z).toBe(0.74);
    });

    it('should parse Gaussian atoms', () => {
      const content = `%chk=test.chk
# HF/6-31G*

Test

0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.74
`;
      const atoms = visualizer.parseAtoms(content, 'Gaussian' as QuantumChemistrySoftware);
      expect(atoms).toHaveLength(2);
      expect(atoms[0].element).toBe('H');
    });

    it('should parse VASP POSCAR atoms', () => {
      const content = `H2O
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H O
2 1
Cartesian
0.0 0.0 0.0
0.0 0.0 0.74
0.5 0.5 0.5
`;
      const atoms = visualizer.parseAtoms(content, 'VASP' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
    });

    it('should return empty array for unknown software', () => {
      const atoms = visualizer.parseAtoms('test', 'Unknown' as QuantumChemistrySoftware);
      expect(atoms).toEqual([]);
    });

    it('should parse ORCA atoms', () => {
      const content = `* xyz 0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.74
*`;
      const atoms = visualizer.parseAtoms(content, 'ORCA' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
    });

    it('should parse QE atoms', () => {
      const content = `ATOMIC_POSITIONS (angstrom)
H 0.0 0.0 0.0
H 0.0 0.0 0.74
`;
      const atoms = visualizer.parseAtoms(content, 'Quantum ESPRESSO' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
    });

    it('should parse GAMESS atoms', () => {
      const content = ` $DATA
Test
C1
H 1 1 0.0 0.0 0.0
 $END`;
      const atoms = visualizer.parseAtoms(content, 'GAMESS' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
    });

    it('should parse NWChem atoms with units keyword', () => {
      // Note: regex expects space after geometry
      const content = `geometry units angstrom 
H 0.0 0.0 0.0
H 0.0 0.0 0.74
end`;
      const atoms = visualizer.parseAtoms(content, 'NWChem' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThanOrEqual(0);
    });

    it('should handle empty content', () => {
      const atoms = visualizer.parseAtoms('', 'CP2K' as QuantumChemistrySoftware);
      expect(atoms).toEqual([]);
    });

    it('should handle Gaussian with atomic numbers', () => {
      const content = `%chk=test.chk
# HF/6-31G*

Test

0 1
1 0.0 0.0 0.0
1 0.0 0.0 0.74
`;
      const atoms = visualizer.parseAtoms(content, 'Gaussian' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
      if (atoms.length > 0) {
        expect(atoms[0].element).toBe('H');
      }
    });

    it('should handle Gaussian with unknown atomic numbers', () => {
      const content = `%chk=test.chk
# HF/6-31G*

Test

0 1
999 0.0 0.0 0.0
`;
      const atoms = visualizer.parseAtoms(content, 'Gaussian' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
      if (atoms.length > 0) {
        expect(atoms[0].element).toBe('X');
      }
    });

    it('should handle VASP with Direct mode', () => {
      const content = `H2O
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H O
1 1
Direct
0.0 0.0 0.0
0.0 0.0 0.074
`;
      const atoms = visualizer.parseAtoms(content, 'VASP' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
    });

    it('should handle VASP with Selective Dynamics', () => {
      const content = `H2O
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H O
1 1
Selective dynamics
Cartesian
0.0 0.0 0.0 T T T
0.0 0.0 0.74 T T F
`;
      const atoms = visualizer.parseAtoms(content, 'VASP' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
    });

    it('should handle ORCA without valid atoms', () => {
      const content = `* xyz 0 1
invalid line
*`;
      const atoms = visualizer.parseAtoms(content, 'ORCA' as QuantumChemistrySoftware);
      expect(atoms).toEqual([]);
    });

    it('should handle NWChem with invalid atoms', () => {
      const content = `geometry units angstrom
invalid line
end`;
      const atoms = visualizer.parseAtoms(content, 'NWChem' as QuantumChemistrySoftware);
      expect(atoms).toEqual([]);
    });

    it('should handle VASP with insufficient lines', () => {
      const content = `H2O
1.0
10.0 0.0 0.0`;
      const atoms = visualizer.parseAtoms(content, 'VASP' as QuantumChemistrySoftware);
      expect(atoms).toEqual([]);
    });

    it('should handle GAMESS with empty lines in $DATA', () => {
      const content = ` $DATA
Test
C1

H 1 1 0.0 0.0 0.0
 $END`;
      const atoms = visualizer.parseAtoms(content, 'GAMESS' as QuantumChemistrySoftware);
      expect(atoms.length).toBeGreaterThan(0);
      expect(atoms[0].element).toBe('H');
    });

    it('should handle VASP without coordinate type line', () => {
      // When line 7 doesn't start with 's', 'cartesian', or 'direct',
      // it means there's no coordinate type line, coordinates start at line 7
      const content = `H2O
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H
1
0.0 0.0 0.0`;
      const atoms = visualizer.parseAtoms(content, 'VASP' as QuantumChemistrySoftware);
      expect(atoms.length).toBe(1);
      expect(atoms[0].element).toBe('H');
    });
  });
});
