# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

