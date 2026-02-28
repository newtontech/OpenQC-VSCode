# OpenQC-VSCode Examples

This directory contains example input files for all supported quantum chemistry packages.

## Directory Structure

```
examples/
├── cp2k/           # CP2K examples
├── vasp/           # VASP examples
├── gaussian/       # Gaussian examples
├── orca/           # ORCA examples
├── qe/             # Quantum ESPRESSO examples
├── gamess/         # GAMESS examples
└── nwchem/         # NWChem examples
```

## CP2K Examples

| File | Description |
|------|-------------|
| `h2o_sp.inp` | Single point energy calculation for water molecule |
| `si_bulk.inp` | Bulk silicon with k-points |
| `ch4_opt.inp` | Methane geometry optimization |

## VASP Examples

| File | Description |
|------|-------------|
| `INCAR-Si` | Band structure calculation setup |
| `INCAR-MD` | Molecular dynamics simulation |
| `POSCAR-Si2` | Silicon dimer structure |
| `KPOINTS-Si` | K-point grid definition |

## Gaussian Examples

| File | Description |
|------|-------------|
| `h2o_opt.com` | Water optimization and frequencies |
| `benzene_sp.com` | Benzene single point with solvation |
| `ts_search.com` | Transition state search |

## ORCA Examples

| File | Description |
|------|-------------|
| `h2o_opt.inp` | Water geometry optimization |
| `h2o_ccsdt.inp` | CCSD(T) high-accuracy calculation |
| `benzene_pcm.inp` | Solvent effects with CPCM |

## Quantum ESPRESSO Examples

| File | Description |
|------|-------------|
| `si_scf.in` | Silicon bulk SCF calculation |
| `h2o_vc_relax.in` | Water cell optimization |
| `si_bands.in` | Silicon band structure |

## GAMESS Examples

| File | Description |
|------|-------------|
| `h2o_scf.inp` | Water SCF with 6-31G* basis |
| `benzene_opt.inp` | Benzene geometry optimization |
| `formaldehyde_tddft.inp` | TDDFT excited states |

## NWChem Examples

| File | Description |
|------|-------------|
| `h2o_scf.nw` | Water SCF calculation |
| `benzene_dft.nw` | DFT optimization with B3LYP |
| `si_bulk.nw` | Periodic boundary conditions |

## Usage

1. Open any example file in VSCode
2. OpenQC extension will automatically detect the file type
3. Use command palette (Ctrl+Shift+P) to access OpenQC features:
   - `OpenQC: Visualize Molecular Structure`
   - `OpenQC: Plot Calculation Data`
   - `OpenQC: Preview Input File`

## Notes

- These examples are designed for educational purposes
- Some calculations may require significant computational resources
- Check pseudopotential/basis set paths before running
- Adjust memory and CPU settings based on your system
