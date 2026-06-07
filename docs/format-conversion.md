# OpenQC Format Converter

## Overview

The Format Converter provides seamless conversion between various quantum chemistry file formats using a Python backend powered by `dpdata`.

## Supported Formats

### Input Formats
- VASP (POSCAR, CONTCAR)
- Gaussian (.gjf, .com)
- ORCA (.inp)
- XYZ
- PDB
- CIF
- CP2K
- Quantum ESPRESSO

### Output Formats
- VASP (POSCAR)
- Gaussian (.gjf, .com)
- ORCA (.inp)
- XYZ
- PDB
- CIF

## Installation

### Prerequisites
- Python 3.7 or higher
- pip (Python package installer)

### Install Python Dependencies

```bash
pip install -r python/requirements.txt
```

Or install the core dependency:

```bash
pip install dpdata
```

ASE is not required for the current converter. Install it only for future
advanced structure manipulation or migration experiments:

```bash
pip install -r python/requirements-ase.txt
```

## Usage

### Command Line

#### Convert a single file:
```bash
python python/format_converter.py input.xyz output.POSCAR
```

#### Convert with explicit format specification:
```bash
python python/format_converter.py input.gjf output.inp --from gaussian --to orca
```

#### Batch convert multiple files:
```bash
python python/format_converter.py --batch file1.xyz file2.xyz --to vasp --output-dir ./converted
```

#### List supported formats:
```bash
python python/format_converter.py --list-formats
```

### VSCode Commands

- `OpenQC: Convert File Format` - Open format picker to convert current file
- `OpenQC: Convert to XYZ` - Quick convert current file to XYZ format
- `OpenQC: Convert to PDB` - Quick convert current file to PDB format
- `OpenQC: Convert to VASP` - Quick convert current file to VASP format
- `OpenQC: Convert to Gaussian` - Quick convert current file to Gaussian format
- `OpenQC: Batch Convert Files` - Batch convert multiple files
- `OpenQC: Check Format Converter Backend` - Verify Python backend is available

## Configuration

Add to your VSCode `settings.json`:

```json
{
  "openqc.converter.pythonPath": "python3",
  "openqc.converter.preserveMetadata": true
}
```

## Examples

### VASP to XYZ
```bash
python python/format_converter.py POSCAR output.xyz
```

### Gaussian to ORCA
```bash
python python/format_converter.py input.gjf output.inp
```

### XYZ to PDB
```bash
python python/format_converter.py molecule.xyz protein.pdb
```

## Testing

```bash
# Run unit tests
npm run test:unit -- --testPathPattern=converters

# Run integration tests (requires Python backend)
npm run test:integration -- --testPathPattern=formatConversion

# Quick test conversion
make -f Makefile.converter test-vasp-to-xyz
```

## Architecture

- **Python Backend** (`python/format_converter.py`): Uses `dpdata` library for format conversion
- **TypeScript Adapter** (`src/converters/FormatConverter.ts`): VSCode extension integration
- **Command Handlers** (`src/commands/formatConversionCommands.ts`): User interface commands

## Metadata Preservation

The converter preserves metadata during conversion:
- Number of atoms
- Element types and counts
- Number of frames (for trajectories)
- Original format information

## Error Handling

- Invalid input files: Returns descriptive error message
- Unsupported formats: Falls back to auto-detection
- Missing Python backend: Prompts user to install dependencies
- Conversion errors: Reports specific error type and message

## Future Enhancements

- Support for additional formats (CASTEP, VASP OUTCAR, etc.)
- Optional ASE integration for advanced structure operations; see
  [ASE Integration Architecture Evaluation](./architecture/ASE_INTEGRATION_EVALUATION.md)
- Preview of conversion results before saving
- Custom format templates
