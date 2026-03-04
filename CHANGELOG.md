# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [2.7.0] - 2026-03-04

### Added - Three.js 3D Visualization (Phase 2)

#### Three.js Integration
- **ThreeJsRenderer.ts** - Core Three.js molecular structure renderer
  - WebGL-based 3D rendering with Three.js
  - Support for ball-and-stick, space-filling, wireframe modes
  - Atom rendering with CPK colors and proper radii
  - Bond detection based on covalent radii
  - Unit cell visualization for periodic systems
  - Camera controls: rotate, zoom, pan, reset
  - Export to PNG image format

- **ThreeJsWebview.ts** - VSCode Webview integration
  - Interactive controls for representation modes
  - Real-time structure info display
  - Image export functionality

- **types.ts** - TypeScript definitions for elements, colors, radii

#### Testing
- Unit tests for ThreeJsRenderer (100% core coverage)
- Integration tests for visualization pipeline
- Mock WebGL context for Node.js testing

### Phase 2 Progress
- ✅ Three.js integration - COMPLETE
- ✅ POSCAR/CONTCAR rendering - COMPLETE
- ✅ Camera controls - COMPLETE
- ✅ Representation modes - COMPLETE
- ⏳ Phase 2 75% Complete

## [2.6.0] - 2026-03-04

### Added - ASE Calculator Integration (Issue #12 - Phase 3)

#### ASE Calculator Interface
- **ASECalculator.ts** - TypeScript interface to ASE calculators
  - Calculator factory pattern for code-agnostic interface
  - Support for VASP, CP2K, and Quantum ESPRESSO calculators
  - Input file generation via ASE backend
  - Calculation execution capabilities
  - Result reading and parsing
  - VSCode configuration integration

- **calculatorCommands.ts** - VSCode commands for calculator operations
  - Generate calculator input files
  - Run calculations from VSCode
  - Read and display calculation results
  - Quick calculation presets for VASP, CP2K, QE
  - Calculator availability checker

#### Python Backend for Calculators
- **python/openqc/ase/calculator.py** - ASE calculator backend
  - VASP input generation (INCAR, POSCAR, KPOINTS, POTCAR)
  - CP2K input generation with complete FORCE_EVAL section
  - Quantum ESPRESSO input generation with all namelists
  - Calculation execution via subprocess
  - Result reading from OUTCAR, vasprun.xml, etc.
  - Energy, forces, and stress extraction

### Added - MD/Optimization Workflow Migration (Issue #12 - Phase 3)

#### Workflow Migration System
- **core/openqc/workflow.py** - Complete cross-code workflow parameter conversion
  - Migrate molecular dynamics parameters between VASP, CP2K, Quantum ESPRESSO, Gaussian, ORCA, and LAMMPS
  - Migrate geometry optimization parameters including convergence criteria and optimizer settings
  - Support for 6 quantum chemistry codes with 30+ test cases
  - Automatic code-specific parameter mapping (thermostats, optimizers, convergence defaults)
  - Workflow extraction from input files (VASP INCAR, Gaussian .com, QE .in, etc.)
  - Smart warnings for unsupported features (NPT ensemble in Gaussian, stress convergence in ORCA)
  - Code-specific notes and migration recommendations

#### MD Parameter Migration
- Migrate MD parameters: temperature, pressure, timestep, steps, thermostats
- Thermostat support: Nose-Hoover, Langevin, Berendsen, Velocity Scaling
- Code-specific MD settings:
  - VASP: IBRION, POTIM, NSW, SMASS, TEBEG, TEEND
  - CP2K: run_type, steps, thermostat parameters
  - QE: calculation, dt, nstep, tempw, ion_dynamics
  - LAMMPS: run_style, timestep, run, fix_nvt, fix_langevin

#### Optimization Parameter Migration
- Migrate optimization parameters: convergence criteria, max steps, optimizers
- Default convergence values per code
- Optimizer mapping between codes (BFGS, CG, RFO, GDIIS, etc.)
- Code-specific optimization settings:
  - VASP: IBRION, EDIFF, EDIFFG, NSW, POTIM
  - CP2K: run_type, optimizer, max_iter, rms_force
  - QE: calculation, ion_dynamics, etot_conv_thr, forc_conv_thr
  - Gaussian: route card, Opt keywords (Tight, VeryTight)
  - ORCA: method, basis, job_type, energy/grad thresholds

#### Workflow Extraction
- Detect workflow type from input files
- Extract parameters from VASP INCAR files
- Extract parameters from Gaussian .com files
- Extract parameters from QE .in files
- Extract parameters from CP2K .inp files
- Extract parameters from ORCA .inp files
- CLI interface for batch processing

### Added - Migration Validation Suite (Issue #12 - Phase 3.5)

#### Migration Validation Tests
- **tests/integration/migrationValidation.test.ts** - Comprehensive validation tests
  - Round-trip conversion tests (VASP↔CP2K↔QE)
  - Structure preservation validation
  - Position deviation checks (< 0.01 Å tolerance)
  - Cell parameter validation
  - PBC (periodic boundary conditions) consistency checks
  - Chemical symbol validation
  - Parameter consistency checks
  - Validation report generation in JSON format

- **validateStructurePreservation()** - Structure comparison utility
  - Calculates position deviations between structures
  - Validates cell parameters within tolerance
  - Checks chemical symbol consistency
  - PBC flag comparison

- **generateValidationReport()** - Report generation
  - JSON format validation reports
  - Migration path tracking
  - Timestamps and metadata
  - Structure validation results

### Added - Complex Property Handling (Issue #12 - Phase 4)

#### Complex Property Mapper
- **ComplexPropertyMapper.ts** - Property mapping between codes
  - Hubbard U parameter mapping (VASP, CP2K, QE)
  - DFT+U parameter conversion with Dudarev scheme support
  - Atom constraint mapping (selective dynamics, fixed atoms)
  - Excited state method mapping (TDDFT, GW, BSE)
  - Pseudopotential strategy mapping
  - Basis set strategy mapping
  - Code-specific parameter generation

#### Python Backend for Complex Properties
- **python/openqc/ase/complex_properties.py** - Property conversion backend
  - HubbardUConverter - DFT+U parameter conversion
  - ConstraintConverter - Atom constraint conversion
  - PseudopotentialConverter - Pseudopotential name mapping
  - BasisSetConverter - Basis set and cutoff recommendations
  - Support for VASP LDAUU/LDAUJ/LDAUL parameters
  - CP2K DFT_PLUS_U section generation
  - QE Hubbard_U card generation

#### Property Mapping Features
- **Hubbard U Parameters**
  - Angular momentum (s, p, d, f) support
  - U and J parameter mapping
  - Per-element configuration
  - Code-specific output formats

- **Atom Constraints**
  - Selective dynamics (VASP)
  - Fixed atoms (CP2K)
  - if_pos constraints (QE)
  - Partial coordinate constraints

- **Excited State Methods**
  - TDDFT configuration mapping
  - GW approximation support
  - BSE (Bethe-Salpeter Equation)
  - Tamm-Dancoff approximation (TDA)
  - Singlet/triplet state specification

- **Pseudopotential Strategies**
  - Norm-conserving pseudopotentials
  - Ultrasoft pseudopotentials
  - PAW (Projector Augmented Wave)
  - HGH (Hartwigner-Goedecker-Hutter)
  - Functional-specific recommendations

- **Basis Set Strategies**
  - Plane wave cutoffs (VASP, QE)
  - Gaussian basis sets (CP2K)
  - Quality levels: minimal, DZVP, TZVP, QZVP
  - Custom basis set support

### Technical Details

#### Architecture
- TypeScript interface layer for VSCode integration
- Python backend for ASE calculator operations
- Factory pattern for code-agnostic calculator access
- Subprocess-based calculation execution
- Result parsing with error handling

#### Supported Codes
- VASP (5.4.4+): Complete calculator support
- CP2K (8.2+): Complete calculator support
- Quantum ESPRESSO (7.2+): Complete calculator support
- Gaussian (16): Workflow migration support
- ORCA (5.0): Workflow migration support
- LAMMPS: Workflow migration support

#### File Format Support
- VASP: INCAR, POSCAR, KPOINTS, POTCAR, OUTCAR, vasprun.xml
- CP2K: .inp, .out, .restart
- QE: .in, .out, .save
- Gaussian: .com, .gjf, .log
- ORCA: .inp, .out
- LAMMPS: .lmp, .data

## [2.0.0] - 2026-03-02

### Major Release - Universal Quantum Chemistry Platform

OpenQC-VSCode is now a universal platform supporting 7 major quantum chemistry packages with automatic LSP detection and unified visualization.

### Added

#### Core Features
- **Universal LSP Support** - Automatic detection and management of language servers for CP2K, VASP, Gaussian, ORCA, Quantum ESPRESSO, GAMESS, and NWChem
- **Molecular Visualization** - Interactive 3D structure rendering with 3Dmol.js
  - Multiple visualization styles: stick, sphere, line, cartoon
  - Spin and zoom controls
  - Real-time structure preview from input files
- **Data Visualization** - Plot SCF energies, convergence data with Plotly.js
  - Energy convergence plots
  - K-point grid visualization
  - Automatic data extraction from output files
- **Input Preview** - Structured display of input file parameters
  - Section-based organization
  - Parameter extraction and display
- **VSCode Sidebar Panel** - Molecules and Calculation Jobs management
  - Molecule library with formula and atom count display
