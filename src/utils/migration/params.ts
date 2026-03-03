/**
 * Parameter Mapping Database for Quantum Chemistry Code Migration
 *
 * Maps common parameters across VASP, CP2K, QE, Gaussian, ORCA, NWChem, GAMESS, LAMMPS
 */

/**
 * Parameter mapping for a specific property
 */
export interface ParameterMapping {
  /** Source code (e.g., 'vasp', 'cp2k', 'qe') */
  source: string;
  /** Target code */
  target: string;
  /** Parameter name in source code */
  sourceParam: string;
  /** Parameter name in target code */
  targetParam: string;
  /** Conversion function (if needed) */
  conversion?: (value: any) => any;
  /** Unit conversion factor */
  unitFactor?: number;
  /** Description */
  description: string;
  /** Is this parameter required? */
  required?: boolean;
  /** Default value if not specified */
  defaultValue?: any;
}

/**
 * Energy cutoff mappings
 */
export const ENERGY_CUTOFF_MAPPINGS: ParameterMapping[] = [
  {
    source: 'vasp',
    target: 'cp2k',
    sourceParam: 'ENCUT',
    targetParam: 'CUTOFF',
    unitFactor: 1.0,
    description: 'Plane-wave cutoff energy (eV)',
    required: true,
  },
  {
    source: 'vasp',
    target: 'qe',
    sourceParam: 'ENCUT',
    targetParam: 'ecutwfc',
    unitFactor: 1.0,
    description: 'Kinetic energy cutoff for wavefunctions (Ry)',
    required: true,
  },
  {
    source: 'qe',
    target: 'vasp',
    sourceParam: 'ecutwfc',
    targetParam: 'ENCUT',
    unitFactor: 1.0,
    description: 'Kinetic energy cutoff (eV)',
    required: true,
  },
  {
    source: 'cp2k',
    target: 'vasp',
    sourceParam: 'CUTOFF',
    targetParam: 'ENCUT',
    unitFactor: 1.0,
    description: 'Plane-wave cutoff energy (eV)',
    required: true,
  },
];

/**
 * Exchange-correlation functional mappings
 */
export const XC_FUNCTIONAL_MAPPINGS: ParameterMapping[] = [
  {
    source: 'vasp',
    target: 'cp2k',
    sourceParam: 'GGA',
    targetParam: 'FUNCTIONAL',
    conversion: (vaspGGA: string) => {
      const mapping: Record<string, string> = {
        PE: 'PBE',
        PE91: 'PBE',
        RPBE: 'RPBE',
        PBE: 'PBE',
        PS: 'PBE',
        AM: 'Pade',
      };
      return mapping[vaspGGA.toUpperCase()] || vaspGGA;
    },
    description: 'Exchange-correlation functional',
    required: true,
  },
  {
    source: 'vasp',
    target: 'qe',
    sourceParam: 'GGA',
    targetParam: 'input_dft',
    conversion: (vaspGGA: string) => {
      const mapping: Record<string, string> = {
        PE: 'PBE',
        PE91: 'PBE',
        RPBE: 'RPBE',
        PBE: 'PBE',
      };
      return mapping[vaspGGA.toUpperCase()] || vaspGGA;
    },
    description: 'Exchange-correlation functional',
    required: true,
  },
  {
    source: 'qe',
    target: 'vasp',
    sourceParam: 'input_dft',
    targetParam: 'GGA',
    conversion: (qeDFT: string) => {
      const mapping: Record<string, string> = {
        PBE: 'PE',
        RPBE: 'RPBE',
        BLYP: 'BLYP',
      };
      return mapping[qeDFT.toUpperCase()] || qeDFT;
    },
    description: 'Exchange-correlation functional',
  },
];

/**
 * k-Point grid mappings
 */
export const KPOINT_MAPPINGS: ParameterMapping[] = [
  {
    source: 'vasp',
    target: 'qe',
    sourceParam: 'KPOINTS',
    targetParam: 'K_POINTS',
    description: 'k-point grid specification',
    required: true,
  },
  {
    source: 'qe',
    target: 'vasp',
    sourceParam: 'K_POINTS',
    targetParam: 'KPOINTS',
    description: 'k-point grid specification',
    required: true,
  },
  {
    source: 'vasp',
    target: 'cp2k',
    sourceParam: 'KPOINTS',
    targetParam: 'KPOINTS',
    description: 'k-point grid specification (for periodic calculations)',
  },
];

/**
 * Electronic convergence mappings
 */
export const CONVERGENCE_MAPPINGS: ParameterMapping[] = [
  {
    source: 'vasp',
    target: 'cp2k',
    sourceParam: 'EDIFF',
    targetParam: 'EPS_SCF',
    conversion: (vaspEDIFF: number) => vaspEDIFF,
    description: 'Electronic convergence criterion',
    required: true,
  },
  {
    source: 'vasp',
    target: 'qe',
    sourceParam: 'EDIFF',
    targetParam: 'conv_thr',
    conversion: (vaspEDIFF: number) => vaspEDIFF,
    description: 'Convergence threshold for self-consistency',
    required: true,
  },
  {
    source: 'qe',
    target: 'vasp',
    sourceParam: 'conv_thr',
    targetParam: 'EDIFF',
    description: 'Electronic convergence criterion',
    required: true,
  },
];

/**
 * Structure parameter mappings (lattice, positions, constraints)
 */
export const STRUCTURE_MAPPINGS: ParameterMapping[] = [
  {
    source: 'vasp',
    target: 'cp2k',
    sourceParam: 'POSCAR',
    targetParam: 'COORD',
    description: 'Atomic structure and coordinates',
    required: true,
  },
  {
    source: 'vasp',
    target: 'qe',
    sourceParam: 'POSCAR',
    targetParam: 'ATOMIC_POSITIONS',
    description: 'Atomic positions and cell',
    required: true,
  },
  {
    source: 'qe',
    target: 'vasp',
    sourceParam: 'ATOMIC_POSITIONS',
    targetParam: 'POSCAR',
    description: 'Atomic structure',
    required: true,
  },
  {
    source: 'cp2k',
    target: 'vasp',
    sourceParam: 'COORD',
    targetParam: 'POSCAR',
    description: 'Atomic structure',
    required: true,
  },
];

/**
 * Molecular dynamics mappings
 */
export const MD_MAPPINGS: ParameterMapping[] = [
  {
    source: 'vasp',
    target: 'cp2k',
    sourceParam: 'POTIM',
    targetParam: 'TIMESTEP',
    unitFactor: 1.0,
    description: 'MD time step (fs)',
  },
  {
    source: 'vasp',
    target: 'qe',
    sourceParam: 'POTIM',
    targetParam: 'dt',
    unitFactor: 1.0,
    description: 'MD time step (fs)',
  },
  {
    source: 'qe',
    target: 'vasp',
    sourceParam: 'dt',
    targetParam: 'POTIM',
    description: 'MD time step (fs)',
  },
];

/**
 * Get all parameter mappings for a source-target pair
 */
export function getParameterMappings(source: string, target: string): ParameterMapping[] {
  const allMappings = [
    ...ENERGY_CUTOFF_MAPPINGS,
    ...XC_FUNCTIONAL_MAPPINGS,
    ...KPOINT_MAPPINGS,
    ...CONVERGENCE_MAPPINGS,
    ...STRUCTURE_MAPPINGS,
    ...MD_MAPPINGS,
  ];

  return allMappings.filter(
    m =>
      m.source.toLowerCase() === source.toLowerCase() &&
      m.target.toLowerCase() === target.toLowerCase()
  );
}

/**
 * Get mapping for a specific parameter
 */
export function getParameterMapping(
  source: string,
  target: string,
  sourceParam: string
): ParameterMapping | undefined {
  const mappings = getParameterMappings(source, target);
  return mappings.find(m => m.sourceParam.toLowerCase() === sourceParam.toLowerCase());
}

/**
 * Convert parameter value using mapping
 */
export function convertParameterValue(mapping: ParameterMapping, value: any): any {
  if (mapping.conversion) {
    return mapping.conversion(value);
  }
  if (mapping.unitFactor) {
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    return numValue * mapping.unitFactor;
  }
  return value;
}

/**
 * Validate if migration from source to target is supported
 */
export function isMigrationSupported(source: string, target: string): boolean {
  const supportedPairs = [
    ['vasp', 'cp2k'],
    ['vasp', 'qe'],
    ['vasp', 'gaussian'],
    ['vasp', 'orca'],
    ['cp2k', 'vasp'],
    ['cp2k', 'qe'],
    ['qe', 'vasp'],
    ['qe', 'cp2k'],
    ['qe', 'gaussian'],
    ['gaussian', 'vasp'],
    ['gaussian', 'qe'],
    ['gaussian', 'orca'],
    ['orca', 'gaussian'],
    ['xyz', 'vasp'],
    ['xyz', 'cp2k'],
    ['xyz', 'qe'],
    ['cif', 'vasp'],
    ['cif', 'qe'],
  ];

  return supportedPairs.some(([s, t]) => s === source.toLowerCase() && t === target.toLowerCase());
}

/**
 * Get list of supported source formats for a target
 */
export function getSupportedSources(target: string): string[] {
  const sources = new Set<string>();
  const supportedPairs = [
    ['vasp', 'cp2k'],
    ['vasp', 'qe'],
    ['vasp', 'gaussian'],
    ['vasp', 'orca'],
    ['cp2k', 'vasp'],
    ['cp2k', 'qe'],
    ['qe', 'vasp'],
    ['qe', 'cp2k'],
    ['qe', 'gaussian'],
    ['gaussian', 'vasp'],
    ['gaussian', 'qe'],
    ['gaussian', 'orca'],
    ['orca', 'gaussian'],
    ['xyz', 'vasp'],
    ['xyz', 'cp2k'],
    ['xyz', 'qe'],
    ['cif', 'vasp'],
    ['cif', 'qe'],
  ];

  supportedPairs.forEach(([s, t]) => {
    if (t === target.toLowerCase()) {
      sources.add(s);
    }
  });

  return Array.from(sources);
}

/**
 * Get list of supported target formats for a source
 */
export function getSupportedTargets(source: string): string[] {
  const targets = new Set<string>();
  const supportedPairs = [
    ['vasp', 'cp2k'],
    ['vasp', 'qe'],
    ['vasp', 'gaussian'],
    ['vasp', 'orca'],
    ['cp2k', 'vasp'],
    ['cp2k', 'qe'],
    ['qe', 'vasp'],
    ['qe', 'cp2k'],
    ['qe', 'gaussian'],
    ['gaussian', 'vasp'],
    ['gaussian', 'qe'],
    ['gaussian', 'orca'],
    ['orca', 'gaussian'],
    ['xyz', 'vasp'],
    ['xyz', 'cp2k'],
    ['xyz', 'qe'],
    ['cif', 'vasp'],
    ['cif', 'qe'],
  ];

  supportedPairs.forEach(([s, t]) => {
    if (s === source.toLowerCase()) {
      targets.add(t);
    }
  });

  return Array.from(targets);
}
