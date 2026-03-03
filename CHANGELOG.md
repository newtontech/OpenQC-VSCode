# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
  - Job tracking with status icons and progress indicators
  - Job management: run, cancel, restart, view results
  - Data export to JSON/CSV format

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

