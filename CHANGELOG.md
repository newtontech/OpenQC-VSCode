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
