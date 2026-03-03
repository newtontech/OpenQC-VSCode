# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- TypeScript interfaces for all property types
- Factory pattern for property mappers
- Code-agnostic property definitions
- Extensible mapping framework

#### Testing
- Unit tests for property mapping functions
- Integration tests for round-trip conversions
- Structure validation with tolerance checks
- Report generation verification

#### Documentation
- Comprehensive JSDoc comments
- Type definitions for all interfaces
- Usage examples in test files

## [2.5.1] - 2026-03-03

### Fixed - MD Workflow Test Failures

#### Number Formatting Fixes
- **CP2K Output Generation** - Fixed number formatting for consistency
  - TIMESTEP now uses `toFixed(1)` to ensure decimal display (e.g., `1.0`)
  - PRESSURE now uses `toFixed(1)` for consistent formatting
  - EPS_ENERGY now uses `toExponential(0).replace('e-', 'e-0')` for scientific notation
  - EPS_FORCE now uses `toExponential(0).replace('e-', 'e-0')` for scientific notation

#### Extraction Logic Fixes
- **VASP Pressure Extraction** - Fixed extraction of `PSTRESS = 0.0`
  - Changed from truthy check to `!== undefined` check
  - Now correctly extracts zero pressure values

- **QE Ensemble Detection** - Fixed NVT ensemble detection
  - Added check for `ion_temperature` parameter
  - Correctly identifies NVT when `ion_temperature = 'nose'` with `ion_dynamics = 'verlet'`

#### Test Fixes
- Fixed test expectations for thermostat type (`NOOSE_HOOVER` → `NOSE_HOOVER` typo)

#### Code Quality
- All 593 tests now passing (34 test suites)
- Prettier formatting applied
- ESLint warnings reviewed (pre-existing)

## [2.5.0] - 2026-03-03


### Added - MD/Optimization Workflow Migration (Issue #12 - Phase 3.4)

#### MD/Optimization Workflow Migration Module
- **Comprehensive MD Parameter Extraction** - Extract MD parameters from VASP, QE, and CP2K inputs
  - Time step and number of steps
  - Ensemble types: NVE, NVT, NPT, NPH
  - Temperature and temperature ramping
  - Pressure for NPT ensembles
  - Thermostat types: Nose-Hoover, Berendsen, Langevin, Velocity Scaling, Andersen
  - Barostat types: Nose-Hoover, Berendsen, Parrinello-Rahman, MTTK
  - Damping times for thermostats/barostats
  - Chain length for Nose-Hoover thermostats
  - Print and trajectory frequencies

#### Optimization Parameter Extraction
- **Optimization Algorithm Detection** - Automatic detection of optimization algorithms
  - VASP IBRION mapping (0=MD, 1=RMMDIIS, 2=CG, 3=DampedMD, 5=BFGS, 7=MD)
  - QE ion_dynamics mapping (bfgs, cg, lbfgs, fire, verlet, damp)
  - CP2K optimizer selection (BFGS, LBFGS, CG)

- **Convergence Parameters** - Extract convergence criteria
  - Energy convergence (eV, Ry, hartree)
  - Force convergence (eV/Å, hartree/bohr)
  - Stress convergence (kbar)

- **Cell Optimization** - Detect and extract cell optimization parameters
  - VASP ISIF configuration (2=ions only, 3=ions+cell, 7=cell+stress)
  - QE CELL dynamics
  - CP2K CELL_OPT and GEO_OPT sections
  - Cell pressure for variable-cell optimization

#### Workflow Conversion
- **Cross-Code MD Parameter Conversion** - Convert MD parameters between codes
  - VASP → QE: NVT/NPT, thermostats, pressures
  - VASP → CP2K: Pressure unit conversion (kbar→bar)
  - QE → VASP: Thermostat mapping (Langevin→Nose-Hoover fallback)
  - QE → CP2K: Ensemble and thermostat conversion
  - CP2K → VASP: Pressure unit conversion (bar→kbar)
  - CP2K → QE: Ensemble and thermostat conversion

- **Optimization Algorithm Mapping** - Map algorithms across codes
  - VASP ↔ QE: CG↔cg, BFGS↔bfgs, LBFGS↔lbfgs
  - VASP ↔ CP2K: Algorithm and convergence unit conversion
  - QE ↔ VASP: Algorithm mapping and Ry→eV conversion

- **Unit Conversions**
  - Energy: eV ↔ Ry (1 eV = 0.0734986 Ry)
  - Energy: eV ↔ hartree (1 eV = 0.0367493 hartree)
  - Force: eV/Å ↔ hartree/bohr (1 eV/Å = 0.01943 hartree/bohr)
  - Pressure: kbar ↔ bar (1 kbar = 1000 bar)

#### Output Generation
- **VASP Output Generation** - Generate VASP INCAR MD and optimization sections
  - MD: POTIM, NSW, IBRION=0, ISIF, TEBEG, TEEND, PSTRESS, SMASS
  - Optimization: IBRION, NSW, EDIFF, EDIFFG, ISIF, PSTRESS
  - Number formatting: toFixed(1) for floats, toExponential(1) for EDIFF

- **QE Output Generation** - Generate QE namelists for MD and optimization
  - CONTROL section: dt, nstep
  - IONS section: ion_dynamics, tempw, templ, ion_temperature
  - CELL section: press, cell_dynamics
  - Algorithm mapping: CG, BFGS, LBFGS → ion_dynamics

- **CP2K Output Generation** - Generate CP2K input sections
  - MOTION/MD section: ENSEMBLE, TIMESTEP, STEPS, TEMPERATURE
  - MOTION/MD/THERMOSTAT section: TYPE, TIMECON, NOSE.LENGTH
  - MOTION/MD/BAROSTAT section: TYPE, PRESSURE, TIMECON
  - MOTION/GEO_OPT section: OPTIMIZER, MAX_ITER, CONVERGENCE
  - MOTION/CELL_OPT section: OPTIMIZER, MAX_ITER, PRESSURE

#### VSCode Integration
- **New Command**: `OpenQC: Migrate MD/Optimization Workflow`
  - Comprehensive MD and optimization workflow migration
  - Auto-detect workflow type (MD vs optimization)
  - Extract ensemble, thermostat, barostat, algorithm parameters
  - Convert between VASP, QE, CP2K, LAMMPS
  - Generate formatted output for target code
  - Copy to clipboard or open preview

- **Enhanced Legacy Command**: `OpenQC: Migrate MD Parameters`
  - Improved with workflow migration suggestion
  - Better user guidance for comprehensive migration

#### Technical Details
- Modular architecture with separate extractors for each code
- Type-safe TypeScript interfaces for all parameters
- Comprehensive error handling and validation
- 71.96% code coverage for mdWorkflow.ts
- 34 unit tests covering extraction, conversion, and generation

### Phase 3 Progress Update
- ✅ Week 12-13 Task 1: Structure Migration Tool - COMPLETE
- ✅ Week 12-13 Task 2: k-Point Grid Migration - COMPLETE
- ✅ Week 12-13 Task 3: Electronic Structure Parameter Migration - COMPLETE
- ✅ Week 12-13 Task 4: MD/Optimization Workflow Migration - COMPLETE

#### Next Steps (Phase 3.5)
- Migration Validation Suite
- Automated round-trip conversion tests
- Energy/force consistency checks
- Integration test with real research workflows

### Changed
- Extended migrationCommands.ts with MD workflow migration command
- Improved legacy MD parameters migration with better user guidance
- Enhanced parameter conversion with workflow-level context

### Fixed
- VASP IBRION workflow type detection (added IBRION=5,7 support)
- CP2K parameter extraction (THERMOSTAT.TYPE, BAROSTAT.TYPE)
- Number formatting in output generation
- TypeScript compilation errors

### Known Limitations
- Some ensemble/thermostat combinations may not have direct equivalents
- Complex parameter combinations may require manual review
- Integration testing with real files recommended



## [2.0.0] - 2026-03-02

## [2.2.0] - 2026-03-03

### Added - Phase 3: Workflow Migration Tools (Issue #12)

#### Structure Migration Tool
- **Parameter Mapping Database** - Comprehensive parameter mappings across quantum chemistry codes
  - Energy cutoff mappings (ENCUT/CUTOFF/ecutwfc)
  - Exchange-correlation functional mappings (PBE, RPBE, etc.)
  - k-Point grid mappings (Monkhorst-Pack ↔ Gamma-centered)
  - Convergence parameter mappings (EDIFF/EPS_SCF/conv_thr)
  - Structure parameter mappings (POSCAR/COORD/ATOMIC_POSITIONS)
  - MD parameter mappings (POTIM/dt/TIMESTEP)
- **Structure Migration Engine** - Full-featured migration with validation
  - Bidirectional format conversion (VASP ↔ CP2K ↔ QE ↔ Gaussian ↔ ORCA)
  - Automatic format detection from file extension and content
  - Structure validation (atom count, elements, cell vectors)
  - Round-trip validation for accuracy
  - Warnings for constraints and structural changes
  - Metadata extraction (atom count, formula, conversions)
- **VSCode Commands** - User-friendly migration commands
  - `OpenQC: Migrate Structure Format` - Quick migration with format selection
  - `OpenQC: Migrate Structure Format (Advanced)` - Full options (constraints, validation, output location)
  - `OpenQC: Migrate k-Point Grid` - k-point migration (UI placeholder for Phase 3.2)
  - `OpenQC: Migrate Electronic Parameters` - Parameter mapping viewer
  - `OpenQC: Migrate MD Parameters` - MD parameters migration (UI placeholder for Phase 3.4)
  - `OpenQC: Show Migration Validation` - Validation report display

#### Supported Migration Paths
- VASP → CP2K, QE, Gaussian, ORCA, XYZ
- CP2K → VASP, QE
- QE → VASP, CP2K, Gaussian
- Gaussian → VASP, QE, ORCA
- ORCA → Gaussian
- XYZ → VASP, CP2K, QE
- CIF → VASP, QE

#### Technical Details
- Parameter mapping database with automatic conversion functions
- Unit conversion support for numerical parameters
- Extensible architecture for adding new mappings
- Comprehensive test suite for migration validation
- Integration with existing ASE converter backend

#### Testing
- Unit tests for structure migration engine
- Parameter mapping tests for all supported codes
- Round-trip conversion validation tests
- Structure preservation tests (atom count, elements, cell)

### Changed
- Extended ASE integration with migration capabilities
- Enhanced VSCode command palette with migration commands
- Improved error reporting for migration failures

### Phase 3 Progress
- ✅ Week 12-13 Task 1: Structure Migration Tool - COMPLETE
- ⏳ Week 12-13 Task 2: k-Point Grid Migration - IN PROGRESS
- ⏳ Week 12-13 Task 3: Electronic Structure Parameter Migration - PARTIAL (parameter mappings complete, conversion tool pending)
- ⏳ Week 12-13 Task 4: MD/Optimization Workflow Migration - PENDING

### Next Steps (Phase 3.2)
- Implement automatic k-point grid generation
- Add Monkhorst-Pack ↔ Gamma-centered conversion
- Support k-point density preservation across codes


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
  - Job tracking with status icons and progress indicators
  - Job management: run, cancel, restart, view results
  - Data export to JSON/CSV format


## [2.3.0] - 2026-03-03

### Added - k-Point Grid Migration (Issue #12 - Phase 3.2)

#### k-Point Grid Migration Module
- **Format Support** - Full support for VASP, Quantum ESPRESSO, and CP2K k-point formats
  - VASP KPOINTS format (Automatic, Gamma-centered, Monkhorst-Pack, Explicit)
  - Quantum ESPRESSO K_POINTS section (automatic, gamma, crystal, tpiba)
  - CP2K &KPOINTS section (Gamma and Monkhorst-Pack schemes)
- **Grid Type Conversion** - Monkhorst-Pack ↔ Gamma-centered conversion
  - Automatic grid doubling when converting MP → Gamma
  - Automatic grid halving when converting Gamma → MP
  - Preserves k-point sampling density
- **k-Point Density Calculations**
  - Calculate density from grid and lattice vectors (kpoints per Å⁻¹)
  - Calculate optimal grid from target density
  - Volume-based density preservation across migrations
- **Auto-Format Detection** - Automatically detect k-point format from content
- **VSCode Command** - `OpenQC: Migrate k-Point Grid`
  - Interactive format selection
  - Grid type conversion options
  - Progress notifications
  - Success/error reporting

#### Technical Details
- Modular architecture with separate parsers/generators for each format
- Comprehensive error handling and validation
- Full TypeScript type safety
- Lattice vector support for density calculations

#### Testing
- Unit tests for all parsers (VASP, QE, CP2K)
- Unit tests for generators (VASP, QE, CP2K)
- Grid conversion tests (MP ↔ Gamma)
- Density calculation tests
- Format detection tests
- Full migration workflow tests
- Error handling and edge case tests
- **Code Coverage: 90%+**

### Phase 3 Progress Update
- ✅ Week 12-13 Task 1: Structure Migration Tool - COMPLETE
- ✅ Week 12-13 Task 2: k-Point Grid Migration - COMPLETE
- ⏳ Week 12-13 Task 3: Electronic Structure Parameter Migration - PARTIAL
- ⏳ Week 12-13 Task 4: MD/Optimization Workflow Migration - PENDING

#### Language Support
- Syntax highlighting for all 7 quantum chemistry formats
- File type auto-detection
- Language server management (start/stop/restart)
- Error diagnostics from LSPs

#### Developer Tools
- Completion provider for all supported formats
- Hover provider for parameter documentation
- Definition provider for navigation
- Diagnostics provider for validation
- Auto-start LSP on document open

#### Testing
- Comprehensive unit tests for all major components
- LSP Manager tests with full lifecycle coverage
- Parser tests for VASP with TDD approach
- Visualization flow tests
- DataPlotter and StructureViewer tests

### Changed
- Improved parser performance with better error handling
- Enhanced file type detection with multi-layer confidence scoring
- Better error messages throughout the extension
- Refactored visualization architecture for better maintainability

### Fixed
- Memory leaks in visualization panels
- Parser edge cases for various quantum chemistry formats
- LSP client cleanup on document close
- Sidebar state persistence

## [0.1.0] - TBD

### Added
- Basic VSCode extension scaffold
- VASP file parsers (INCAR, POSCAR, KPOINTS)
- Gaussian input file parser
- ORCA input file parser
- Syntax highlighting for supported formats
- Basic validation and error reporting
- Unit test framework (Jest + pytest)
- CI/CD pipeline with GitHub Actions
- Documentation website

### Security
- Secure API key storage
- Input validation for external processes

## [0.2.0] - TBD

### Added
- 3D molecular visualization (NGL Viewer)
- Interactive atom selection
- Multiple representation modes
- Real-time structure editing

### Changed
- Improved parser performance
- Better error messages

## [0.3.0] - TBD

### Added
- Format conversion via dpdata
- ASE integration
- PyMOL export
- Batch conversion support

### Fixed
- Memory leaks in visualization
- Parser edge cases

## [0.4.0] - TBD

### Added
- AI-powered parameter suggestions
- Natural language input generation
- Error debugging assistant
- Parameter explanations

### Changed
- Refactored AI service architecture

## [1.0.0] - TBD

### Added
- Full feature set complete
- Production-ready release
- Comprehensive documentation
- Video tutorials
- Example gallery

### Performance
- Optimized for large systems (>1000 atoms)
- Lazy loading and WebWorkers
- Memory optimization

---

## Version History

| Version | Date | Description |
|---------|------|-------------|
| 0.1.0 | TBD | Foundation - Basic parsing and highlighting |
| 0.2.0 | TBD | Visualization - 3D rendering |
| 0.3.0 | TBD | Conversion - Format conversion |
| 0.4.0 | TBD | AI - Smart features |
| 1.0.0 | TBD | Production - Full release |

---

[Unreleased]: https://github.com/newtontech/OpenQC-VSCode/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/newtontech/OpenQC-VSCode/releases/tag/v0.1.0
[0.2.0]: https://github.com/newtontech/OpenQC-VSCode/releases/tag/v0.2.0
[0.3.0]: https://github.com/newtontech/OpenQC-VSCode/releases/tag/v0.3.0
[0.4.0]: https://github.com/newtontech/OpenQC-VSCode/releases/tag/v0.4.0
[1.0.0]: https://github.com/newtontech/OpenQC-VSCode/releases/tag/v1.0.0

## [2.2.0] - 2026-03-03 (Dynamic LSP Discovery - Issue #13)

### Added

#### Dynamic LSP Discovery (Issue #13 - Phase 2.5)
- **LSPDiscovery Module** (`src/utils/LSPDiscovery.ts`)
  - GitHub API integration for fetching LSP repositories from OpenQuantumChemistry org
  - Automatic discovery of 7+ LSP servers: VASP, Gaussian, ORCA, CP2K, Quantum ESPRESSO, GAMESS, NWChem
  - Smart caching with 1-hour TTL to avoid API rate limits
  - Graceful error handling with hardcoded fallback definitions
  - Cache persistence via VSCode global state

- **LSPManager Integration**
  - `discoverAvailableLSPs()` - async method to fetch available LSPs dynamically
  - `refreshLSPList()` - force refresh LSP list from GitHub API
  - Constructor now accepts optional VSCode ExtensionContext for cache storage

- **Update Script** (`scripts/update-lsp-list.sh`)
  - Manual refresh script for debugging and validation
  - Supports JSON, raw, and pretty output formats

- **Testing**
  - Comprehensive Jest test suite (`tests/unit/utils/LSPDiscovery.test.ts`)
  - 8 tests covering API failures, caching, and fallback behavior
  - ~81% code coverage for LSPDiscovery module

### Technical Details
- GitHub API: `https://api.github.com/orgs/OpenQuantumChemistry/repos`
- Cache TTL: 60 minutes
- Fallback definitions for offline operation
- Rate limit handling with informative error messages


## [2.1.0] - 2026-03-03 (ASE Integration Phase 1)

### Added

#### ASE Integration (Issue #12 - Phase 3)
- **ASE Converter Backend** - Python module for format conversion
  - Support for 11 quantum chemistry formats
  - ASE Atoms object as universal intermediate layer
  - Validation and utility functions
- **TypeScript Interface** - VSCode integration for ASE backend
  - ASEConverter class with subprocess communication
  - 4 new VSCode commands:
    - `OpenQC: Convert to ASE Atoms`
    - `OpenQC: Convert from ASE Atoms`
    - `OpenQC: Migrate Format (ASE)`
    - `OpenQC: Quick Convert (ASE)`
- **Quick Convert Presets** - Common format conversions
  - VASP → CP2K, VASP → QE, CP2K → VASP
  - Gaussian → ORCA, ORCA → Gaussian
  - XYZ → VASP, CIF → VASP
- **Configuration** - Python path setting for ASE backend
- **Documentation** - Comprehensive ASE_INTEGRATION.md guide

### Changed
- Updated extension.ts to register ASE commands
- Enhanced package.json with ASE command categories

### Testing
- Python test suite for ASE converter (pytest)
- TypeScript test suite for VSCode integration (Jest)
- Test fixtures: H2O.xyz, Si.POSCAR, CH4.gjf

### Technical Details
- Backend: Python 3.8+ with ASE>=3.22.0
- Format support: VASP, CP2K, QE, Gaussian, ORCA, NWChem, GAMESS, LAMMPS, XYZ, PDB, CIF
- Bidirectional conversion via ASE Atoms objects
- Error handling and progress notifications

### Phase 3 Progress
- ✅ Week 10-11: ASE Core Integration - COMPLETE
  - ASE Atoms Converter Module
  - Format conversions for all 8 quantum chemistry codes
  - Calculator interface wrapper
  - Structure validation utilities

### Next Steps (Phase 2 - Week 12-13)
- Structure Migration Tool with code selection
- k-Point Grid Migration (Monkhorst-Pack ↔ Gamma-centered)
- Electronic Structure Parameter Mapping
- MD/Optimization Workflow Migration

### Known Limitations
- Calculator execution (Phase 3) not yet implemented
- Complex features (Hubbard U, constraints) require manual adjustment
- Pseudopotential mapping is code-specific

---


## [2.4.0] - 2026-03-03

### Added - Electronic Structure Parameter Migration (Issue #12 - Phase 3.3)

#### Parameter Conversion Module
- **Parameter Converter Engine** - Full-featured parameter extraction and conversion
  - VASP INCAR parameter extraction with full parsing support
  - Quantum ESPRESSO namelist parameter extraction (CONTROL, SYSTEM, ELECTRONS sections)
  - CP2K parameter extraction with section-aware parsing
  - Gaussian route section parameter extraction
  - Automatic format detection from file extension and content
  - Comprehensive error handling and validation

#### Parameter Conversion Features
- **Bi-directional Conversion** - Convert parameters between quantum chemistry codes
  - VASP ↔ CP2K parameter conversion
  - VASP ↔ Quantum ESPRESSO parameter conversion
  - VASP ↔ Gaussian parameter conversion
  - CP2K ↔ VASP, QE parameter conversion
  - QE ↔ VASP, CP2K, Gaussian parameter conversion

- **Intelligent Parameter Mapping**
  - Energy cutoff mappings (ENCUT/CUTOFF/ecutwfc)
  - Exchange-correlation functional mappings (PBE, RPBE, etc.)
  - k-Point grid mappings
  - Convergence parameter mappings (EDIFF/EPS_SCF/conv_thr)
  - MD parameter mappings (POTIM/dt/TIMESTEP)

- **Section-Aware Parsing**
  - Handles QE/CP2K section prefixes (CONTROL., SYSTEM., ELECTRONS.)
  - Supports nested parameter lookup for format conversion
  - Preserves parameter metadata and line numbers

#### VSCode Integration
- **Parameter Migration Command** (`OpenQC: Migrate Parameters`)
  - Extract parameters from active editor
  - Show extracted parameters summary
  - Convert to target format with interactive UI
  - Display converted parameters with unmapped warnings
  - Copy parameters to clipboard
  - View parameters in output channel

- **MD Parameter Migration** (`OpenQC: Migrate MD Parameters`)
  - Extract MD-related parameters (POTIM, TEBEG, PRESS, etc.)
  - Convert MD parameters between codes
  - Filter and display only MD-specific parameters

#### Technical Details
- Full TypeScript type safety with comprehensive interfaces
- Modular architecture for easy extension
- Error handling with detailed error messages
- Warnings for unmapped parameters
- Extensive test coverage with 16 test cases

#### Testing
- Unit tests for VASP parameter extraction
- Unit tests for QE parameter extraction
- Unit tests for CP2K parameter extraction
- Unit tests for Gaussian parameter extraction
- Parameter conversion tests (VASP→QE, VASP→CP2K, QE→VASP)
- Functional conversion tests (GGA→input_dft)
- Unmapped parameter warning tests
- Error handling tests
- **Code Coverage: 74.86% for parameterConverter.ts**

### Phase 3 Progress Update
- ✅ Week 12-13 Task 1: Structure Migration Tool - COMPLETE
- ✅ Week 12-13 Task 2: k-Point Grid Migration - COMPLETE
- ✅ Week 12-13 Task 3: Electronic Structure Parameter Migration - COMPLETE
- ⏳ Week 12-13 Task 4: MD/Optimization Workflow Migration - IN PROGRESS

#### Next Steps (Phase 3.4)
- Complete MD/Optimization parameter conversion with UI
- Add automatic MD workflow generation
- Implement optimization criteria migration
- Add MD trajectory analysis tools

### Changed
- Enhanced migrationCommands.ts with real parameter conversion implementation
- Improved error reporting for parameter conversion failures
- Better user feedback during parameter extraction

### Fixed
- Section-prefix handling for QE/CP2K parameters
- Test failures related to nested parameter names
- Temporary file handling in tests
