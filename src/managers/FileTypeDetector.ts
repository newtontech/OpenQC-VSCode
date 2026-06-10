import * as path from 'path';
import * as vscode from 'vscode';

export type QuantumChemistrySoftware =
  | 'ABACUS'
  | 'ABINIT'
  | 'CIF'
  | 'CP2K'
  | 'VASP'
  | 'Gaussian'
  | 'ORCA'
  | 'Quantum ESPRESSO'
  | 'GAMESS'
  | 'NWChem'
  | 'GPUMD'
  | 'GROMACS'
  | 'LAMMPS'
  | 'MLIP'
  | 'PyATB'
  | 'PySCF';

interface FilePattern {
  software: QuantumChemistrySoftware;
  extensions?: string[];
  filenames?: string[];
  contentPatterns?: RegExp[];
}

export class FileTypeDetector {
  private readonly patterns: FilePattern[] = [
    {
      software: 'ABACUS',
      filenames: ['INPUT', 'STRU', 'KPT'],
      contentPatterns: [/INPUT_PARAMETERS|ATOMIC_POSITIONS|K_POINTS/i, /basis_type|pseudo_dir/i],
    },
    {
      software: 'ABINIT',
      extensions: ['.abi', '.abinit'],
      contentPatterns: [/^\s*(ecut|ngkpt|nstep|toldfe|ixc)\b/im, /\bacell\b|\brprim\b|\bxred\b/i],
    },
    {
      software: 'CIF',
      extensions: ['.cif'],
      contentPatterns: [/^data_/im, /_cell_length_[abc]/i, /_atom_site_/i],
    },
    {
      software: 'CP2K',
      extensions: ['.inp'],
      contentPatterns: [/&GLOBAL/i, /&FORCE_EVAL/i, /PROJECT_NAME/i],
    },
    {
      software: 'VASP',
      filenames: [
        'INCAR',
        'POSCAR',
        'KPOINTS',
        'POTCAR',
        'CONTCAR',
        'OSZICAR',
        'OUTCAR',
        'vasprun.xml',
      ],
      contentPatterns: [
        /ISTART|ICHARG|ENCUT|PREC/i,
        /^\s*\d+\.\d+/m, // POSCAR coordinate pattern
      ],
    },
    {
      software: 'Gaussian',
      extensions: ['.gjf', '.com'],
      contentPatterns: [/^%chk=/i, /^#.*(?:HF|DFT|MP2|CCSD|B3LYP)/i, /^(?:\s*\d+\s+\d+\s*)$/m],
    },
    {
      software: 'ORCA',
      extensions: ['.inp'],
      contentPatterns: [/^!.*(?:HF|DFT|MP2|CCSD)/i, /%pal/i, /%maxcore/i],
    },
    {
      software: 'Quantum ESPRESSO',
      extensions: ['.in', '.pw.in', '.relax.in', '.vc-relax.in'],
      contentPatterns: [/&CONTROL/i, /&SYSTEM/i, /calculation\s*=/i, /pseudo_dir/i],
    },
    {
      software: 'GAMESS',
      extensions: ['.inp'],
      contentPatterns: [/^\s*\$BASIS/i, /^\s*\$CONTRL/i, /^\s*\$SYSTEM/i, /runtyp/i],
    },
    {
      software: 'NWChem',
      extensions: ['.nw', '.nwinp'],
      contentPatterns: [
        /^(?:geometry|basis|scf|dft|mp2|ccsd)/i,
        /(?:geometry|basis|scf|dft)\s+\w+/i,
      ],
    },
    {
      software: 'GPUMD',
      filenames: ['run.in', 'nep.in'],
      contentPatterns: [
        /^\s*(potential|velocity|ensemble|dump_thermo|run)\b/im,
        /\b(nep|gpumd)\b/i,
      ],
    },
    {
      software: 'GROMACS',
      extensions: ['.top', '.itp', '.mdp', '.gro'],
      contentPatterns: [
        /\[\s*(moleculetype|atoms|system|molecules|defaults)\s*\]/i,
        /^\s*(integrator|nsteps|dt)\s*=/im,
      ],
    },
    {
      software: 'LAMMPS',
      extensions: ['.lmp', '.lammps', '.lmps'],
      filenames: ['in.lammps'],
      contentPatterns: [/^\s*(units|atom_style|pair_style|read_data|run)\b/im, /\bLAMMPS\b/i],
    },
    {
      software: 'MLIP',
      extensions: ['.mlip.json', '.mlip.yaml', '.mlip.yml'],
      filenames: ['mlip.json', 'mlip.yaml', 'mlip.yml'],
      contentPatterns: [
        /"model"\s*:|model\s*:/i,
        /"structure"\s*:|structure\s*:/i,
        /"task"\s*:|task\s*:/i,
      ],
    },
    {
      software: 'PyATB',
      extensions: ['.pyatb.py'],
      filenames: ['run_pyatb.py'],
      contentPatterns: [/import\s+pyatb|from\s+pyatb\s+import/i, /\b(hr_file|sr_file|HR\.dat)\b/i],
    },
    {
      software: 'PySCF',
      extensions: ['.pyscf.py'],
      filenames: ['run_pyscf.py'],
      contentPatterns: [/from\s+pyscf\s+import|import\s+pyscf/i, /\b(gto|scf|dft|cc|mp)\b/i],
    },
  ];

  /**
   * Detect the quantum chemistry software represented by a VS Code document.
   *
   * The detector prefers exact filename matches, then validates ambiguous file
   * extensions with content patterns, and finally falls back to content-only matching.
   *
   * @param document - VS Code document to inspect.
   * @returns Detected software name, or null when no supported format matches.
   */
  detectSoftware(document: vscode.TextDocument): QuantumChemistrySoftware | null {
    const basename = path.basename(document.fileName);
    const lowerBasename = basename.toLowerCase();
    const extension = path.extname(basename);
    const content = document.getText();

    // Check filename matches first
    for (const pattern of this.patterns) {
      if (pattern.filenames?.includes(basename)) {
        return pattern.software;
      }
    }

    // Check extension matches
    for (const pattern of this.patterns) {
      if (
        this.matchesExtension(lowerBasename, extension, pattern.extensions) &&
        pattern.contentPatterns
      ) {
        // For ambiguous extensions, check content
        const confidence = this.calculateConfidence(content, pattern.contentPatterns);
        if (confidence > 0.5) {
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

  private matchesExtension(
    lowerBasename: string,
    extension: string,
    extensions?: string[]
  ): boolean {
    return Boolean(
      extensions?.some(candidate => {
        const lowered = candidate.toLowerCase();
        return extension.toLowerCase() === lowered || lowerBasename.endsWith(lowered);
      })
    );
  }

  /**
   * Return display metadata for a supported quantum chemistry package.
   *
   * @param software - Supported software identifier.
   * @returns Human-readable name, description, and website metadata.
   */
  getSoftwareInfo(software: QuantumChemistrySoftware): {
    name: string;
    description: string;
    website: string;
  } {
    const info: Record<
      QuantumChemistrySoftware,
      { name: string; description: string; website: string }
    > = {
      CP2K: {
        name: 'CP2K',
        description: 'Quantum chemistry and solid state physics software package',
        website: 'https://www.cp2k.org',
      },
      VASP: {
        name: 'VASP',
        description: 'Vienna Ab initio Simulation Package',
        website: 'https://www.vasp.at',
      },
      Gaussian: {
        name: 'Gaussian',
        description: 'Comprehensive computational chemistry software',
        website: 'https://gaussian.com',
      },
      ORCA: {
        name: 'ORCA',
        description: 'Quantum chemistry program package',
        website: 'https://orcaforum.kofo.mpg.de',
      },
      'Quantum ESPRESSO': {
        name: 'Quantum ESPRESSO',
        description: 'Integrated suite of open-source computer codes',
        website: 'https://www.quantum-espresso.org',
      },
      GAMESS: {
        name: 'GAMESS',
        description: 'General Atomic and Molecular Electronic Structure System',
        website: 'https://www.msg.chem.iastate.edu/gamess',
      },
      NWChem: {
        name: 'NWChem',
        description: 'Open Source High-Performance Computational Chemistry',
        website: 'https://nwchemgit.github.io',
      },
      ABACUS: {
        name: 'ABACUS',
        description: 'Atomic-orbital based ab initio computation at USTC',
        website: 'https://abacus.ustc.edu.cn',
      },
      ABINIT: {
        name: 'ABINIT',
        description: 'First-principles materials and nanostructure simulation suite',
        website: 'https://www.abinit.org',
      },
      CIF: {
        name: 'CIF',
        description: 'Crystallographic Information File format',
        website: 'https://www.iucr.org/resources/cif',
      },
      GPUMD: {
        name: 'GPUMD',
        description: 'Graphics Processing Units Molecular Dynamics',
        website: 'https://gpumd.org',
      },
      GROMACS: {
        name: 'GROMACS',
        description: 'Molecular dynamics package for biomolecular simulation',
        website: 'https://www.gromacs.org',
      },
      LAMMPS: {
        name: 'LAMMPS',
        description: 'Large-scale Atomic/Molecular Massively Parallel Simulator',
        website: 'https://www.lammps.org',
      },
      MLIP: {
        name: 'MLIP',
        description: 'Machine-learning interatomic potential workflow files',
        website: 'https://github.com/newtontech/mlip-lsp',
      },
      PyATB: {
        name: 'PyATB',
        description: 'Python workflows for Atomic Tight-Binding analysis',
        website: 'https://github.com/newtontech/pyatb-lsp',
      },
      PySCF: {
        name: 'PySCF',
        description: 'Python-based Simulations of Chemistry Framework',
        website: 'https://pyscf.org',
      },
    };
    return info[software];
  }
}
