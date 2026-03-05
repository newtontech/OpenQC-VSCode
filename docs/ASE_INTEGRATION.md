# ASE Integration

## Overview

OpenQC-VSCode now supports [ASE (Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/) integration for cross-code workflow migration. This enables seamless conversion between different quantum chemistry packages.

## Features

### Format Conversion

Convert between the following formats:

- **VASP** (POSCAR, CONTCAR)
- **CP2K** (.inp)
- **Quantum ESPRESSO** (.in, .pw.in)
- **Gaussian** (.com, .gjf)
- **ORCA** (.inp)
- **NWChem** (.nw, .nwinp)
- **GAMESS** (.inp)
- **LAMMPS** (.lmp, .data)
- **XYZ** (.xyz)
- **PDB** (.pdb)
- **CIF** (.cif)

### Commands

- **OpenQC: Convert to ASE Atoms** - Convert current file to ASE Atoms format
- **OpenQC: Convert from ASE Atoms** - Convert ASE Atoms to target format
- **OpenQC: Migrate Format (ASE)** - Migrate from one format to another with format selection
- **OpenQC: Quick Convert (ASE)** - Quick convert using common presets

### Quick Convert Presets

- VASP → CP2K
- VASP → Quantum ESPRESSO
- CP2K → VASP
- Gaussian → ORCA
- ORCA → Gaussian
- XYZ → VASP
- CIF → VASP

## Installation

### Prerequisites

1. Install Python 3.8 or later
2. Install ASE and dependencies:

```bash
cd python
pip install -r requirements.txt
```

### Configuration

Configure Python path in VSCode settings:

```json
{
  "openqc.pythonPath": "/usr/bin/python3"
}
```

## Usage

### Basic Conversion

1. Open a quantum chemistry input file
2. Press `F1` (or `Ctrl+Shift+P`) to open Command Palette
3. Type "ASE" or select from Quick Convert
4. Choose target format
5. Specify output file name
6. Converted file opens automatically

### Example: Convert VASP to CP2K

1. Open `POSCAR` file
2. Run command `OpenQC: Migrate Format (ASE)`
3. Select "CP2K" as target format
4. Enter output filename (e.g., `cp2k.inp`)
5. File is converted and opened

### Quick Convert Presets

Use presets for common conversions:

1. Open file
2. Run command `OpenQC: Quick Convert (ASE)`
3. Select preset (e.g., "VASP → CP2K")
4. File is converted with suggested filename

## Architecture

### Python Backend

- `python/openqc/ase/converter.py` - Core conversion logic
- `python/openqc/ase/utils.py` - Validation and utility functions

### TypeScript Interface

- `src/ase/ASEConverter.ts` - TypeScript wrapper for Python backend
- `src/ase/commands.ts` - VSCode command handlers
- `src/ase/index.ts` - Module exports

### Integration Points

- **File → ASE Atoms**: Read quantum chemistry files into ASE Atoms
- **ASE Atoms → File**: Write ASE Atoms to quantum chemistry formats
- **Format → Format**: Direct conversion via ASE as intermediate

## Testing

### Python Tests

```bash
cd python
pytest tests/test_ase_converter.py -v
```

### TypeScript Tests

```bash
npm test -- ASEConverter
```

## Known Limitations

1. **Complex Features**: Hubbard U, special constraints, and excited state methods require per-code handling
2. **Pseudopotentials**: Automatic mapping may need manual adjustment
3. **Code-Specific Keywords**: Not all keywords transfer perfectly between codes

## Future Work

- [ ] Calculator integration for running jobs via ASE
- [ ] k-point grid migration with density preservation
- [ ] Electronic structure parameter mapping
- [ ] MD/Optimization workflow migration
- [ ] Hubbert U parameter mapping
- [ ] Advanced constraint handling

## Contributing

To add support for a new format:

1. Update `python/openqc/ase/converter.py` with new format
2. Add format to `get_supported_formats()` in `utils.py`
3. Update `ASEFormat` enum in `src/ase/ASEConverter.ts`
4. Add tests for the new format

## References

- [ASE Documentation](https://wiki.fysik.dtu.dk/ase/)
- [ASE Calculators](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html)
- [dpdata - DeepMD Data Format](https://github.com/deepmodeling/dpdata)

## License

MIT License - See [LICENSE](../../LICENSE) for details.
