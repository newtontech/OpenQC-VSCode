/**
 * Software adapter registry for expandable quantum chemistry support.
 *
 * Provides a registry-style interface for adding new software packages
 * without editing multiple switch statements.
 *
 * @module adapters/SoftwareAdapter
 */

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

export type SoftwareCategory =
  | 'molecular-qc'
  | 'periodic-dft'
  | 'multireference'
  | 'properties'
  | 'ml-md';

export interface OpenQCFilePattern {
  extensions?: string[];
  filenames?: string[];
  contentPatterns?: RegExp[];
}

export interface OpenQCSoftwareAdapter {
  id: string;
  displayName: string;
  category: SoftwareCategory;
  languageId?: string;
  filePatterns: OpenQCFilePattern[];
  description?: string;
  cclibSupport?: boolean;
  nativeParse?: boolean;
  pythonBridge?: boolean;
}

// ---------------------------------------------------------------------------
// Registry
// ---------------------------------------------------------------------------

class AdapterRegistry {
  private adapters: Map<string, OpenQCSoftwareAdapter> = new Map();

  register(adapter: OpenQCSoftwareAdapter): void {
    this.adapters.set(adapter.id, adapter);
  }

  get(id: string): OpenQCSoftwareAdapter | undefined {
    return this.adapters.get(id);
  }

  getAll(): OpenQCSoftwareAdapter[] {
    return Array.from(this.adapters.values());
  }

  getByCategory(category: SoftwareCategory): OpenQCSoftwareAdapter[] {
    return this.getAll().filter(a => a.category === category);
  }

  /** Detect software from file content and name. */
  detect(content: string, fileName: string): OpenQCSoftwareAdapter | undefined {
    const basename = fileName.split(/[/\\]/).pop()?.toUpperCase() ?? '';

    for (const adapter of this.getAll()) {
      // Check filename patterns
      for (const pattern of adapter.filePatterns) {
        if (pattern.filenames) {
          if (pattern.filenames.some(f => basename === f.toUpperCase())) {
            return adapter;
          }
        }
      }

      // Check content patterns
      for (const pattern of adapter.filePatterns) {
        if (pattern.contentPatterns) {
          if (pattern.contentPatterns.some(r => r.test(content))) {
            return adapter;
          }
        }
      }
    }

    return undefined;
  }
}

// ---------------------------------------------------------------------------
// Singleton registry with all adapters
// ---------------------------------------------------------------------------

export const adapterRegistry = new AdapterRegistry();

// -- Existing supported software --

adapterRegistry.register({
  id: 'cp2k',
  displayName: 'CP2K',
  category: 'periodic-dft',
  languageId: 'cp2k',
  filePatterns: [
    {
      extensions: ['.inp'],
      contentPatterns: [/&GLOBAL/i, /&FORCE_EVAL/i],
    },
  ],
  nativeParse: true,
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'vasp',
  displayName: 'VASP',
  category: 'periodic-dft',
  languageId: 'vasp',
  filePatterns: [
    {
      filenames: ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'CONTCAR', 'OSZICAR', 'OUTCAR'],
    },
    {
      extensions: ['.vasp'],
    },
  ],
  nativeParse: true,
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'gaussian',
  displayName: 'Gaussian',
  category: 'molecular-qc',
  languageId: 'gaussian',
  filePatterns: [
    { extensions: ['.gjf', '.com'] },
    { contentPatterns: [/^%chk=/m, /^#\s*(?:HF|B3LYP|MP2)/m] },
  ],
  cclibSupport: true,
  nativeParse: true,
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'orca',
  displayName: 'ORCA',
  category: 'molecular-qc',
  languageId: 'orca',
  filePatterns: [
    {
      extensions: ['.inp'],
      contentPatterns: [/^!\s*(?:HF|DFT|MP2|CCSD)/i, /%pal/i, /\* xyz/i],
    },
  ],
  cclibSupport: true,
  nativeParse: true,
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'qe',
  displayName: 'Quantum ESPRESSO',
  category: 'periodic-dft',
  languageId: 'qe',
  filePatterns: [
    {
      extensions: ['.in', '.pw.in', '.relax.in', '.vc-relax.in'],
      contentPatterns: [/&CONTROL/i, /&SYSTEM/i],
    },
  ],
  nativeParse: true,
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'gamess',
  displayName: 'GAMESS',
  category: 'molecular-qc',
  languageId: 'gamess',
  filePatterns: [
    {
      extensions: ['.inp'],
      contentPatterns: [/^\s*\$BASIS/i, /^\s*\$CONTRL/i],
    },
  ],
  cclibSupport: true,
  nativeParse: true,
});

adapterRegistry.register({
  id: 'nwchem',
  displayName: 'NWChem',
  category: 'molecular-qc',
  languageId: 'nwchem',
  filePatterns: [
    {
      extensions: ['.nw', '.nwinp'],
      contentPatterns: [/^(?:geometry|basis|scf|dft|task)/i],
    },
  ],
  cclibSupport: true,
  nativeParse: true,
});

// -- Coming Soon: Molecular QC (#79) --

adapterRegistry.register({
  id: 'psi4',
  displayName: 'Psi4',
  category: 'molecular-qc',
  description: 'Open-source quantum chemistry with Python API',
  filePatterns: [
    {
      extensions: ['.in', '.dat'],
      contentPatterns: [/set\s+\{/i, /energy\(/i, /molecule\s+\{/i, /^import\s+psi4/m],
    },
  ],
  cclibSupport: true,
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'molpro',
  displayName: 'Molpro',
  category: 'molecular-qc',
  description: 'High-accuracy quantum chemistry with explicitly correlated methods',
  filePatterns: [
    {
      extensions: ['.com', '.inp'],
      contentPatterns: [/^\s*\{/, /rhf\s*;/i, /geometry\s*=/i],
    },
  ],
  cclibSupport: true,
});

adapterRegistry.register({
  id: 'openmolcas',
  displayName: 'OpenMolcas',
  category: 'multireference',
  description: 'Multiconfigurational quantum chemistry (CASSCF, CASPT2)',
  filePatterns: [
    {
      extensions: ['.input', '.inp'],
      contentPatterns: [/&SEWARD/i, /&SCF/i, /&RASSCF/i, /&CASPT2/i],
    },
  ],
  cclibSupport: true,
});

adapterRegistry.register({
  id: 'dalton',
  displayName: 'DALTON',
  category: 'properties',
  description: 'Molecular properties: polarizabilities, NMR, EPR, UV/vis',
  filePatterns: [
    {
      extensions: ['.dal'],
      contentPatterns: [/^\*\*DALTON/i, /\.RUN\s/i],
    },
  ],
  cclibSupport: true,
});

adapterRegistry.register({
  id: 'turbomole',
  displayName: 'Turbomole',
  category: 'molecular-qc',
  description: 'Efficient DFT and MP2 for molecules',
  filePatterns: [
    {
      filenames: ['CONTROL', 'COORD', 'BASIS', 'ENERGIES'],
      contentPatterns: [/\$end/i],
    },
  ],
  cclibSupport: true,
});

// -- Coming Soon: Periodic DFT (#79) --

adapterRegistry.register({
  id: 'castep',
  displayName: 'CASTEP',
  category: 'periodic-dft',
  description: 'Plane-wave DFT for materials science',
  filePatterns: [
    {
      extensions: ['.cell', '.param'],
      contentPatterns: [/^\s*%block\s+lattice_cart/i, /%block\s+positions_abs/i],
    },
    {
      extensions: ['.castep'],
    },
  ],
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'abinit',
  displayName: 'Abinit',
  category: 'periodic-dft',
  description: 'Plane-wave DFT with GW and TDDFT',
  filePatterns: [
    {
      extensions: ['.abi', '.abo'],
      contentPatterns: [/ndtset\s+/i, /ecut\s+/i, /nband\s+/i],
    },
  ],
  pythonBridge: true,
});

adapterRegistry.register({
  id: 'crystal',
  displayName: 'Crystal',
  category: 'periodic-dft',
  description: 'Gaussian basis periodic DFT',
  filePatterns: [
    {
      extensions: ['.d12', '.gui'],
      contentPatterns: [/CRYSTAL/i, /SLAB/i, /POLYMER/i, /MOLECULE/i],
    },
  ],
  pythonBridge: true,
});
