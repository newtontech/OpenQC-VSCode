import { QuantumChemistrySoftware } from '../managers/FileTypeDetector';

interface Atom {
    elem: string;
    x: number;
    y: number;
    z: number;
}

export class Molecule3D {
    parseAtoms(content: string, software: QuantumChemistrySoftware): Atom[] {
        switch (software) {
            case 'CP2K':
                return this.parseCp2kAtoms(content);
            case 'VASP':
                return this.parseVaspAtoms(content);
            case 'Gaussian':
                return this.parseGaussianAtoms(content);
            case 'ORCA':
                return this.parseOrcaAtoms(content);
            case 'Quantum ESPRESSO':
                return this.parseQeAtoms(content);
            case 'GAMESS':
                return this.parseGamessAtoms(content);
            case 'NWChem':
                return this.parseNwchemAtoms(content);
            default:
                return [];
        }
    }

    private parseCp2kAtoms(content: string): Atom[] {
        const atoms: Atom[] = [];
        // Look for &COORD section
        const coordMatch = content.match(/&COORD[\s\S]*?&END\s*COORD/i);
        if (coordMatch) {
            const lines = coordMatch[0].split('\n').slice(1, -1);
            for (const line of lines) {
                const parts = line.trim().split(/\s+/);
                if (parts.length >= 4) {
                    atoms.push({
                        elem: parts[0],
                        x: parseFloat(parts[1]),
                        y: parseFloat(parts[2]),
                        z: parseFloat(parts[3])
                    });
                }
            }
        }
        return atoms;
    }

    private parseVaspAtoms(content: string): Atom[] {
        const atoms: Atom[] = [];
        const lines = content.split('\n');
        
        // POSCAR format
        // Line 1: Comment
        // Line 2: Scale factor
        // Lines 3-5: Lattice vectors
        // Line 6: Element names
        // Line 7: Atom counts
        
        if (lines.length < 7) return atoms;
        
        const elementNames = lines[5].trim().split(/\s+/);
        const atomCounts = lines[6].trim().split(/\s+/).map(Number);
        
        // Check if it's selective dynamics or Cartesian
        const modeLine = lines[7].trim().toLowerCase();
        const startLine = (modeLine.startsWith('s') || modeLine === 'cartesian' || modeLine === 'direct') ? 8 : 7;
        
        let atomIndex = 0;
        for (let i = 0; i < elementNames.length; i++) {
            const elem = elementNames[i];
            const count = atomCounts[i];
            for (let j = 0; j < count; j++) {
                const line = lines[startLine + atomIndex];
                if (line) {
                    const parts = line.trim().split(/\s+/);
                    if (parts.length >= 3) {
                        atoms.push({
                            elem,
                            x: parseFloat(parts[0]),
                            y: parseFloat(parts[1]),
                            z: parseFloat(parts[2])
                        });
                    }
                }
                atomIndex++;
            }
        }
        
        return atoms;
    }

    private parseGaussianAtoms(content: string): Atom[] {
        const atoms: Atom[] = [];
        const lines = content.split('\n');
        let inAtoms = false;
        let chargeMultFound = false;
        
        for (const line of lines) {
            // Look for charge/multiplicity line (e.g., "0 1")
            if (!chargeMultFound && /^-?\d+\s+-?\d+\s*$/.test(line.trim())) {
                chargeMultFound = true;
                inAtoms = true;
                continue;
            }
            
            // Stop at blank line after atoms
            if (inAtoms && line.trim() === '') {
                break;
            }
            
            // Parse atom lines
            if (inAtoms) {
                const parts = line.trim().split(/\s+/);
                if (parts.length >= 4) {
                    // Check if first part is element symbol or atomic number
                    const firstPart = parts[0];
                    const isElement = /^[A-Za-z]+$/.test(firstPart);
                    
                    if (isElement) {
                        atoms.push({
                            elem: firstPart.charAt(0).toUpperCase() + firstPart.slice(1).toLowerCase(),
                            x: parseFloat(parts[1]),
                            y: parseFloat(parts[2]),
                            z: parseFloat(parts[3])
                        });
                    } else {
                        // Atomic number - convert to element
                        const atomicNumber = parseInt(firstPart);
                        const element = this.atomicNumberToElement(atomicNumber);
                        atoms.push({
                            elem: element,
                            x: parseFloat(parts[1]),
                            y: parseFloat(parts[2]),
                            z: parseFloat(parts[3])
                        });
                    }
                }
            }
        }
        
        return atoms;
    }

    private parseOrcaAtoms(content: string): Atom[] {
        const atoms: Atom[] = [];
        const lines = content.split('\n');
        let inAtoms = false;
        
        for (const line of lines) {
            // Look for coordinate section
            if (line.includes('* xyz') || line.includes('* xyzfile')) {
                inAtoms = true;
                continue;
            }
            
            // End of coordinate section
            if (inAtoms && line.trim() === '*') {
                inAtoms = false;
                continue;
            }
            
            // Parse atom lines
            if (inAtoms) {
                const parts = line.trim().split(/\s+/);
                if (parts.length >= 4) {
                    const elem = parts[0];
                    if (/^[A-Za-z]+$/.test(elem)) {
                        atoms.push({
                            elem: elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase(),
                            x: parseFloat(parts[1]),
                            y: parseFloat(parts[2]),
                            z: parseFloat(parts[3])
                        });
                    }
                }
            }
        }
        
        return atoms;
    }

    private parseQeAtoms(content: string): Atom[] {
        const atoms: Atom[] = [];
        // Look for ATOMIC_POSITIONS section
        const posMatch = content.match(/ATOMIC_POSITIONS[^\n]*\n([\s\S]*?)(?=ATOMIC|CELL|K_POINTS|&|$)/i);
        
        if (posMatch) {
            const lines = posMatch[1].trim().split('\n');
            for (const line of lines) {
                const parts = line.trim().split(/\s+/);
                if (parts.length >= 4) {
                    atoms.push({
                        elem: parts[0],
                        x: parseFloat(parts[1]),
                        y: parseFloat(parts[2]),
                        z: parseFloat(parts[3])
                    });
                }
            }
        }
        
        return atoms;
    }

    private parseGamessAtoms(content: string): Atom[] {
        const atoms: Atom[] = [];
        const lines = content.split('\n');
        let inData = false;
        
        for (const line of lines) {
            // Look for $DATA section
            if (line.includes('$DATA')) {
                inData = true;
                continue;
            }
            
            // End of $DATA
            if (inData && line.includes('$END')) {
                break;
            }
            
            // Skip title and symmetry lines
            if (inAtoms && line.trim().length === 0) continue;
            
            // Parse atom lines in GAMESS format
            if (inData) {
                const parts = line.trim().split(/\s+/);
                // GAMESS format: ELEMENT ATOM# ATOMIC# X Y Z
                if (parts.length >= 6 && /^[A-Za-z]+$/.test(parts[0])) {
                    atoms.push({
                        elem: parts[0],
                        x: parseFloat(parts[3]),
                        y: parseFloat(parts[4]),
                        z: parseFloat(parts[5])
                    });
                }
            }
        }
        
        return atoms;
    }

    private parseNwchemAtoms(content: string): Atom[] {
        const atoms: Atom[] = [];
        // Look for geometry section
        const geomMatch = content.match(/geometry\s+(?:units\s+\w+\s+)?(?:\n[\s\S]*?)end/i);
        
        if (geomMatch) {
            const lines = geomMatch[0].split('\n').slice(1, -1);
            for (const line of lines) {
                const parts = line.trim().split(/\s+/);
                if (parts.length >= 4 && /^[A-Za-z]+$/.test(parts[0])) {
                    atoms.push({
                        elem: parts[0],
                        x: parseFloat(parts[1]),
                        y: parseFloat(parts[2]),
                        z: parseFloat(parts[3])
                    });
                }
            }
        }
        
        return atoms;
    }

    private atomicNumberToElement(atomicNumber: number): string {
        const elements: Record<number, string> = {
            1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
            9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
            16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 26: 'Fe', 29: 'Cu',
            30: 'Zn', 35: 'Br', 47: 'Ag', 53: 'I', 79: 'Au', 80: 'Hg', 82: 'Pb'
        };
        return elements[atomicNumber] || 'X';
    }
}