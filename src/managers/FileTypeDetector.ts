import * as vscode from 'vscode';

export type QuantumChemistrySoftware = 
    | 'CP2K' 
    | 'VASP' 
    | 'Gaussian' 
    | 'ORCA' 
    | 'Quantum ESPRESSO' 
    | 'GAMESS' 
    | 'NWChem';

interface FilePattern {
    software: QuantumChemistrySoftware;
    extensions?: string[];
    filenames?: string[];
    contentPatterns?: RegExp[];
}

export class FileTypeDetector {
    private readonly patterns: FilePattern[] = [
        {
            software: 'CP2K',
            extensions: ['.inp'],
            contentPatterns: [
                /&GLOBAL/i,
                /&FORCE_EVAL/i,
                /PROJECT_NAME/i
            ]
        },
        {
            software: 'VASP',
            filenames: ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR'],
            contentPatterns: [
                /ISTART|ICHARG|ENCUT|PREC/i,
                /^\s*\d+\.\d+/m  // POSCAR coordinate pattern
            ]
        },
        {
            software: 'Gaussian',
            extensions: ['.gjf', '.com'],
            contentPatterns: [
                /^%chk=/i,
                /^#.*(?:HF|DFT|MP2|CCSD|B3LYP)/i,
                /^(?:\s*\d+\s+\d+\s*)$/m
            ]
        },
        {
            software: 'ORCA',
            extensions: ['.inp'],
            contentPatterns: [
                /^!.*(?:HF|DFT|MP2|CCSD)/i,
                /%pal/i,
                /%maxcore/i
            ]
        },
        {
            software: 'Quantum ESPRESSO',
            extensions: ['.in', '.pw.in', '.relax.in', '.vc-relax.in'],
            contentPatterns: [
                /&CONTROL/i,
                /&SYSTEM/i,
                /calculation\s*=/i,
                /pseudo_dir/i
            ]
        },
        {
            software: 'GAMESS',
            extensions: ['.inp'],
            contentPatterns: [
                /^\s*\$BASIS/i,
                /^\s*\$CONTRL/i,
                /^\s*\$SYSTEM/i,
                /runtyp/i
            ]
        },
        {
            software: 'NWChem',
            extensions: ['.nw', '.nwinp'],
            contentPatterns: [
                /^(?:geometry|basis|scf|dft|mp2|ccsd)/i,
                /(?:geometry|basis|scf|dft)\s+\w+/i
            ]
        }
    ];

    detectSoftware(document: vscode.TextDocument): QuantumChemistrySoftware | null {
        const filename = document.fileName;
        const basename = filename.split('/').pop()?.split('\\').pop() || '';
        const extension = basename.includes('.') ? basename.slice(basename.lastIndexOf('.')) : '';
        const content = document.getText();

        // Check filename matches first
        for (const pattern of this.patterns) {
            if (pattern.filenames?.includes(basename)) {
                return pattern.software;
            }
        }

        // Check extension matches
        for (const pattern of this.patterns) {
            if (pattern.extensions?.includes(extension)) {
                // For ambiguous extensions, check content
                if (pattern.contentPatterns) {
                    const confidence = this.calculateConfidence(content, pattern.contentPatterns);
                    if (confidence > 0.5) {
                        return pattern.software;
                    }
                } else {
                    return pattern.software;
                }
            }
        }

        // Fallback: check content patterns for all
        let bestMatch: QuantumChemistrySoftware | null = null;
        let bestConfidence = 0;

        for (const pattern of this.patterns) {
            if (pattern.contentPatterns) {
                const confidence = this.calculateConfidence(content, pattern.contentPatterns);
                if (confidence > bestConfidence && confidence > 0.3) {
                    bestConfidence = confidence;
                    bestMatch = pattern.software;
                }
            }
        }

        return bestMatch;
    }

    private calculateConfidence(content: string, patterns: RegExp[]): number {
        let matches = 0;
        for (const pattern of patterns) {
            if (pattern.test(content)) {
                matches++;
            }
        }
        return matches / patterns.length;
    }

    getSoftwareInfo(software: QuantumChemistrySoftware): {
        name: string;
        description: string;
        website: string;
    } {
        const info: Record<QuantumChemistrySoftware, { name: string; description: string; website: string }> = {
            'CP2K': {
                name: 'CP2K',
                description: 'Quantum chemistry and solid state physics software package',
                website: 'https://www.cp2k.org'
            },
            'VASP': {
                name: 'VASP',
                description: 'Vienna Ab initio Simulation Package',
                website: 'https://www.vasp.at'
            },
            'Gaussian': {
                name: 'Gaussian',
                description: 'Comprehensive computational chemistry software',
                website: 'https://gaussian.com'
            },
            'ORCA': {
                name: 'ORCA',
                description: 'Quantum chemistry program package',
                website: 'https://orcaforum.kofo.mpg.de'
            },
            'Quantum ESPRESSO': {
                name: 'Quantum ESPRESSO',
                description: 'Integrated suite of open-source computer codes',
                website: 'https://www.quantum-espresso.org'
            },
            'GAMESS': {
                name: 'GAMESS',
                description: 'General Atomic and Molecular Electronic Structure System',
                website: 'https://www.msg.chem.iastate.edu/gamess'
            },
            'NWChem': {
                name: 'NWChem',
                description: 'Open Source High-Performance Computational Chemistry',
                website: 'https://nwchemgit.github.io'
            }
        };
        return info[software];
    }
}