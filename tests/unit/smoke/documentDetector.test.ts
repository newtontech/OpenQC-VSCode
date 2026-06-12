/**
 * Document Detector Tests
 *
 * Tests for output/log document detection. Note that the detector prioritizes
 * input file detection (Phase 1) over output detection (Phase 2). Files that
 * appear in the LSP registry as input fileNames or extensions will always be
 * detected as "input" even if they also appear in output patterns.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/159
 */

import {
  detectDocument,
  isOutputOrLogFile,
  getOutputPatternsForLanguage,
  getLogPatternsForLanguage,
} from '../../../src/smoke/documentDetector';

// ---------------------------------------------------------------------------
// Input file detection
// ---------------------------------------------------------------------------

describe('detectDocument - input files', () => {
  it('detects VASP INCAR as input', () => {
    const result = detectDocument('INCAR');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('vasp');
    expect(result.confidence).toBe(1.0);
  });

  it('detects VASP POSCAR as input', () => {
    const result = detectDocument('POSCAR');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('vasp');
  });

  it('detects VASP KPOINTS as input', () => {
    const result = detectDocument('KPOINTS');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('vasp');
  });

  it('detects VASP POTCAR as input', () => {
    const result = detectDocument('POTCAR');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('vasp');
  });

  // Note: OUTCAR and CONTCAR are in the VASP registry's fileNames list,
  // so they are detected as "input" in Phase 1 even though they are
  // conceptually output files. This is correct behavior for the registry
  // alignment -- the registry declares them as files the VASP LSP handles.

  it('detects Gaussian .gjf files as input', () => {
    const result = detectDocument('water.gjf');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('gaussian');
  });

  it('detects Gaussian .com files as input', () => {
    const result = detectDocument('opt.com');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('gaussian');
  });

  it('detects LAMMPS .lmp files as input', () => {
    const result = detectDocument('relax.lmp');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('lammps');
  });

  it('detects LAMMPS in.lammps files as input', () => {
    const result = detectDocument('in.lammps');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('lammps');
  });

  it('detects CIF .cif files as input', () => {
    const result = detectDocument('structure.cif');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('cif');
  });

  it('detects ABACUS INPUT file as input', () => {
    const result = detectDocument('INPUT');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('abacus');
  });

  it('detects ABACUS STRU file as input', () => {
    const result = detectDocument('STRU');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('abacus');
  });

  it('detects NWChem .nw files as input', () => {
    const result = detectDocument('calc.nw');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('nwchem');
  });

  it('detects GROMACS .top files as input', () => {
    const result = detectDocument('system.top');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('gromacs');
  });

  it('detects GROMACS .mdp files as input', () => {
    const result = detectDocument('grompp.mdp');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('gromacs');
  });

  it('detects GPUMD run.in files as input via fileNames match', () => {
    // run.in is in GPUMD's fileNames in the registry
    const result = detectDocument('run.in');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('gpumd');
  });

  it('detects GPUMD nep.in files as input via fileNames match', () => {
    // nep.in is in GPUMD's fileNames in the registry
    const result = detectDocument('nep.in');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('gpumd');
  });

  it('detects QE .scf.in files as input', () => {
    const result = detectDocument('pw.scf.in');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('qe');
  });

  it('detects ABINIT .abi files as input', () => {
    const result = detectDocument('t1.abi');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('abinit');
  });
});

// ---------------------------------------------------------------------------
// Ambiguous extension tests
// ---------------------------------------------------------------------------

describe('detectDocument - ambiguous extensions', () => {
  it('detects .inp files as input (first registry match wins)', () => {
    // .inp is shared by CP2K, ORCA, GAMESS. The registry order determines priority.
    const result = detectDocument('calc.inp');
    expect(result.kind).toBe('input');
    // CP2K comes first in the registry for .inp extension
    expect(result.languageId).toBe('cp2k');
  });

  it('detects .lammps files as input via extension', () => {
    // .lammps extension is registered for LAMMPS
    const result = detectDocument('log.lammps');
    expect(result.kind).toBe('input');
    expect(result.languageId).toBe('lammps');
  });
});

// ---------------------------------------------------------------------------
// Output file detection
// ---------------------------------------------------------------------------

describe('detectDocument - output files', () => {
  it('detects ORCA prefixed .out as output', () => {
    const result = detectDocument('ORCA_job.out');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('orca');
  });

  it('detects GPUMD thermo.out as output', () => {
    const result = detectDocument('thermo.out');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('gpumd');
  });

  it('detects CP2K .ener as output', () => {
    const result = detectDocument('md.ener');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('cp2k');
  });

  it('detects GROMACS .xvg as output', () => {
    const result = detectDocument('energy.xvg');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('gromacs');
  });

  it('detects ABACUS running_*.log as output', () => {
    const result = detectDocument('running_scf.log');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('abacus');
  });

  it('detects NWChem .nwout as output', () => {
    const result = detectDocument('calc.nwout');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('nwchem');
  });

  it('detects VASP DOSCAR as output', () => {
    // DOSCAR is NOT in the registry fileNames, so it hits output patterns
    const result = detectDocument('DOSCAR');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('vasp');
  });

  it('detects VASP WAVECAR as output', () => {
    const result = detectDocument('WAVECAR');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('vasp');
  });

  it('detects LAMMPS dump files as output', () => {
    const result = detectDocument('dump.atoms');
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('lammps');
  });
});

// ---------------------------------------------------------------------------
// Content-based detection
// ---------------------------------------------------------------------------

describe('detectDocument - content detection', () => {
  it('detects Gaussian output from content', () => {
    const content = ' Entering Gaussian System ... some output ... Normal termination of Gaussian.';
    const result = detectDocument('unknown_file.txt', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('gaussian');
  });

  it('detects VASP output content', () => {
    const content = 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM\nE0= -.12345E+02';
    const result = detectDocument('output.dat', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('vasp');
  });

  it('detects LAMMPS output from content', () => {
    const content = 'LAMMPS (29 Oct 2020)\nTotal wall time: 0:00:01';
    const result = detectDocument('out.txt', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('lammps');
  });

  it('detects ORCA output from content', () => {
    const content = '* O   R   C   A *\nORCA TERMINATED NORMALLY';
    const result = detectDocument('result.txt', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('orca');
  });

  it('detects CP2K output from content', () => {
    const content = 'CP2K| Program started at 2026-01-01\nENERGY| Total energy: -123.45';
    const result = detectDocument('output.txt', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('cp2k');
  });

  it('detects QE output from content', () => {
    const content = 'PROGRAM PWSCF ENDED\n some data';
    const result = detectDocument('output.dat', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('qe');
  });

  it('detects GAMESS output from content', () => {
    const content = 'GAMESS VERSION 2022\n some data';
    const result = detectDocument('output.dat', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('gamess');
  });

  it('detects GROMACS output from content', () => {
    const content = 'GROMACS mdrun finished\n some data';
    const result = detectDocument('output.dat', content);
    expect(result.kind).toBe('output');
    expect(result.languageId).toBe('gromacs');
  });
});

// ---------------------------------------------------------------------------
// Unknown files
// ---------------------------------------------------------------------------

describe('detectDocument - unknown files', () => {
  it('returns unknown for unrecognized file names without content', () => {
    const result = detectDocument('random_file.dat');
    expect(result.kind).toBe('unknown');
    expect(result.confidence).toBe(0);
  });

  it('returns unknown for empty string', () => {
    const result = detectDocument('');
    expect(result.kind).toBe('unknown');
  });

  it('returns unknown for unrecognized file with non-matching content', () => {
    const result = detectDocument('random_file.xyz', 'just some random text with no signatures');
    expect(result.kind).toBe('unknown');
  });
});

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

describe('isOutputOrLogFile', () => {
  it('returns true for known output files', () => {
    expect(isOutputOrLogFile('DOSCAR')).toBe(true);
    expect(isOutputOrLogFile('WAVECAR')).toBe(true);
    expect(isOutputOrLogFile('ORCA_job.out')).toBe(true);
    expect(isOutputOrLogFile('thermo.out')).toBe(true);
    expect(isOutputOrLogFile('dump.atoms')).toBe(true);
    expect(isOutputOrLogFile('energy.xvg')).toBe(true);
  });

  it('returns false for input files that are in the registry', () => {
    // OUTCAR is in VASP's fileNames, so it is detected as input
    expect(isOutputOrLogFile('INCAR')).toBe(false);
    expect(isOutputOrLogFile('water.gjf')).toBe(false);
  });

  it('returns false for unknown files', () => {
    expect(isOutputOrLogFile('random.txt')).toBe(false);
  });
});

describe('getOutputPatternsForLanguage', () => {
  it('returns patterns for known languages', () => {
    const patterns = getOutputPatternsForLanguage('vasp');
    expect(patterns.length).toBeGreaterThan(0);
  });

  it('returns empty array for unknown language', () => {
    const patterns = getOutputPatternsForLanguage('nonexistent');
    expect(patterns).toEqual([]);
  });
});

describe('getLogPatternsForLanguage', () => {
  it('returns log patterns for known languages', () => {
    const patterns = getLogPatternsForLanguage('gromacs');
    expect(patterns.length).toBeGreaterThan(0);
  });

  it('returns empty array for unknown language', () => {
    const patterns = getLogPatternsForLanguage('nonexistent');
    expect(patterns).toEqual([]);
  });
});
