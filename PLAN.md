# OpenQC-VSCode - Project Roadmap

> VSCode Extension for Quantum Chemistry Input Files: Visualization, Editing, Conversion, and AI-Assisted Modification

## Project Vision

OpenQC-VSCode aims to become the definitive VSCode extension for quantum chemistry researchers and computational chemists. By providing seamless integration with popular quantum chemistry software (VASP, Gaussian, ORCA, etc.), powerful visualization capabilities, and AI-assisted input file modification, we enable researchers to focus on science rather than file format wrestling.

### Core Goals

1. **Universal Compatibility** - Support all major quantum chemistry input file formats
2. **Intuitive Visualization** - 3D molecular structure rendering directly in VSCode
3. **Seamless Conversion** - Convert between different formats with zero data loss
4. **AI-Powered Assistance** - Leverage AI to suggest modifications and optimizations
5. **Extensible Architecture** - Plugin system for custom workflows and integrations

---

## Phase 1: Foundation (Weeks 1-4)

### Objectives
- Establish core architecture
- Basic file parsing and syntax highlighting
- Initial visualization capabilities

### Deliverables

#### Week 1-2: Project Setup
- [x] VSCode extension scaffolding
- [x] Build system configuration (TypeScript compiler)
- [x] Test framework setup (Jest + ts-jest, 184 tests, 97.97% coverage)
- [x] CI/CD pipeline (GitHub Actions: ci.yml, test.yml, release.yml, prod-build.yml)
- [x] Documentation framework (Vitepress)

#### Week 3-4: Core Parsing Engine
- [x] VASP input file parser (INCAR, POSCAR, KPOINTS, POTCAR) - PR: feat/vasp-parser
- [x] Gaussian input file parser (.gjf, .com) - Already implemented
- [x] ORCA input file parser (.inp) - Already implemented
- [x] Abstract parser interface for extensibility - BaseParser abstract class
- [x] Syntax highlighting for all supported formats - 7 tmLanguage.json files
- [x] Basic validation and error reporting - ValidationResult interface

### Technology Stack
- **Language**: TypeScript (extension), Python (backend utilities)
- **Parser**: Chevrotain (TypeScript), Pyparsing (Python fallback)
- **UI**: VSCode Webview API + React
- **Testing**: Jest (TypeScript), pytest (Python)
- **Build**: esbuild (fast bundling)

### Success Metrics
- Parse 95% of sample VASP/Gaussian/ORCA files without errors
- Syntax highlighting covers all major keywords
- Test coverage ≥ 80%

---

## Phase 2: Visualization (Weeks 5-8)

### Objectives
- 3D molecular structure visualization
- Interactive editing capabilities
- Visual feedback for structure modifications

### Deliverables

#### Week 5-6: 3D Rendering
- [ ] Integrate 3D visualization library (Three.js / Mol* / NGL Viewer)
- [ ] POSCAR/CONTCAR structure rendering
- [ ] Gaussian/ORCA geometry visualization
- [ ] Atom/bond rendering with proper coloring
- [ ] Camera controls (rotate, zoom, pan)
- [ ] Multiple representation modes (ball-stick, space-filling, wireframe)

#### Week 7-8: Interactive Editing
- [ ] Click-to-select atoms
- [ ] Real-time coordinate editing in 3D view
- [ ] Bond distance/angle measurements
- [ ] Unit cell visualization for periodic systems
- [ ] Export modified structures back to input format

### Technology Stack
- **3D Engine**: Three.js or NGL Viewer
- **State Management**: Zustand or Redux
- **Integration**: VSCode Webview with bidirectional messaging

### Success Metrics
- Render structures at 60 FPS for systems up to 1000 atoms
- Sub-second load time for typical input files
- Intuitive editing with undo/redo support

---

## Phase 2.5: Dynamic LSP Discovery (Week 8-9)

### Objectives
- 动态获取 OpenQuantumChemistry 组织下的 LSP 仓库列表
- 消除硬编码的 LSP 配置，实现自动发现新 LSP
- 更新 cron job 和自动化脚本以支持动态列表

### GitHub Issue
**#13**: [Feature] 动态获取组织下的 LSP 仓库列表 (https://github.com/newtontech/OpenQC-VSCode/issues/13)

### Background
Currently, OpenQC-VSCode hardcodes the list of supported LSPs in \`package.json\` and \`LSPManager.ts\` (CP2K, VASP, Gaussian, ORCA, QE, GAMESS, NWChem). This requires manual code updates whenever OpenQuantumChemistry organization adds new LSP repositories.

### Available LSP Repositories
From https://github.com/orgs/OpenQuantumChemistry/repositories:
| Repository | Description | Status |
|------------|-------------|--------|
| orca-lsp | ORCA quantum chemistry software LSP | ✅ Active |
| gamess-lsp | GAMESS (US) input files LSP | ✅ Active |
| qe-lsp | Quantum ESPRESSO LSP | ✅ Active |
| cp2k-lsp-enhanced | CP2K input file tools | ✅ Active |
| gaussian-lsp | Gaussian LSP | ✅ Active |
| vasp-lsp | VASP input/output files LSP | ✅ Active |

### Deliverables

#### Week 8: Dynamic Discovery Implementation
- [ ] **GitHub API Integration** (Issue #13 - Phase 1)
  - [ ] Create \`src/utils/LSPDiscovery.ts\` module
  - [ ] Implement \`fetchLSPRepositories()\` using GitHub API
  - [ ] Cache repository list with TTL (e.g., 1 hour)
  - [ ] Handle API rate limiting and errors gracefully
- [ ] **Dynamic Configuration** (Issue #13 - Phase 2)
  - [ ] Modify \`LSPManager.ts\` to use discovered LSP list
  - [ ] Auto-generate VSCode configuration contributions
  - [ ] Support runtime LSP registration/unregistration

#### Week 9: Automation & Tooling Updates
- [ ] **Cron Job Updates** (Issue #13 - Phase 3)
  - [ ] Update \`scripts/qclsp-dev.sh\` to fetch LSP list dynamically
  - [ ] Modify HEARTBEAT.md to discover projects from GitHub API
  - [ ] Create \`scripts/update-lsp-list.sh\` for manual refresh
- [ ] **CI/CD Integration** (Issue #13 - Phase 4)
  - [ ] GitHub Action to validate LSP list daily
  - [ ] Auto-generate PR when new LSP is detected
  - [ ] Update documentation when LSP list changes

### Implementation Details

#### GitHub API Endpoint
\`\`\`bash
# List org repos with 'lsp' in name
curl -s "https://api.github.com/orgs/OpenQuantumChemistry/repos?per_page=100" | \\
  jq '.[] | select(.name | contains("lsp")) | .name'
\`\`\`

#### Dynamic Configuration Schema
\`\`\`typescript
interface LSPServerDefinition {
  id: string;           // e.g., "vasp-lsp"
  name: string;         // e.g., "VASP"
  repository: string;   // e.g., "OpenQuantumChemistry/vasp-lsp"
  executable: string;   // e.g., "vasp-lsp"
  languageId: string;   // e.g., "vasp"
  fileExtensions: string[];  // e.g., ["INCAR", "POSCAR"]
  enabled: boolean;
}
\`\`\`

### Success Metrics
- LSP list auto-updates within 1 hour of new repo creation
- Zero manual code changes required to support new LSP
- Cron jobs automatically include new LSP projects
- Backward compatibility maintained for existing configurations

---

## Phase 3: ASE Integration & Cross-Code Migration (Weeks 10-15)

### Objectives
- **Primary**: Full ASE integration as the "universal intermediate layer"
- Cross-code workflow migration (VASP ↔ CP2K ↔ QE ↔ Gaussian ↔ ORCA)
- Seamless format conversion with ASE Atoms as the bridge
- Integration with external tools and automated workflows

### Why ASE?
ASE (Atomic Simulation Environment) provides:
- **Unified Atoms object**: Single representation for all chemical structures
- **Calculator abstraction**: Code-agnostic interface to DFT/MD/force-field engines
- **Multi-backend support**: VASP, CP2K, QE, Gaussian, ORCA, NWChem, GAMESS, LAMMPS, etc.
- **Workflow capabilities**: Structure manipulation, running jobs, reading results

### GitHub Issue
**#12**: [feat(ase): Integrate ASE for cross-code workflow migration](https://github.com/newtontech/OpenQC-VSCode/issues/12)

### Deliverables

#### Week 10-11: ASE Core Integration
- [ ] **ASE Atoms Converter Module** (Issue #12 - Phase 1)
  - [ ] VASP POSCAR ↔ ASE Atoms
  - [ ] CP2K input ↔ ASE Atoms
  - [ ] QE input ↔ ASE Atoms
  - [ ] Gaussian input ↔ ASE Atoms
  - [ ] ORCA input ↔ ASE Atoms
  - [ ] NWChem input ↔ ASE Atoms
  - [ ] GAMESS input ↔ ASE Atoms
  - [ ] LAMMPS data ↔ ASE Atoms
- [ ] ASE Calculator interface wrapper
- [ ] Unified structure validation via ASE

#### Week 12-13: Cross-Code Migration Tools
- [ ] **Structure Migration Tool** (Issue #12 - Phase 2)
  - VASP→CP2K, QE→Gaussian, etc.
  - Preserve unit cells, atomic positions, constraints
  - Handle periodic boundary conditions
- [ ] **k-Point Grid Migration**
  - Convert between Monkhorst-Pack and Gamma-centered
  - Maintain density of k-points across codes
- [ ] **Electronic Structure Parameter Migration**
  - Map common parameters (ENCUT ↔ cutoff, etc.)
- [ ] **MD/Optimization Workflow Migration**
  - Convert MD parameters, optimization criteria

#### Week 14-15: Advanced ASE Features
- [ ] **ASE Calculator Integration** (Issue #12 - Phase 3)
  - Direct job execution via ASE
  - Generate inputs and run calculations
  - Read results back into OpenQC-VSCode
- [ ] **Complex Property Handling** (Issue #12 - Phase 4)
  - Hubbard U parameters with per-code mapping
  - Special constraints and frozen atoms
  - Excited state methods
  - Pseudopotential/basis set strategies
- [ ] **Migration Validation Suite**
  - Automated tests for round-trip conversions
  - Energy/force consistency checks

### Technology Stack
- **Conversion**: dpdata (Python backend)
- **Visualization**: PyMOL (external), NGL (internal)
- **Integration**: Python subprocess + VSCode tasks

### Success Metrics
- 99% conversion accuracy (validated against reference datasets)
- Support for 10+ file format pairs
- Integration tests with real research workflows

---

## Phase 4: AI Assistance (Weeks 16-19)

### Objectives
- AI-powered input file optimization
- Intelligent suggestions and validation
- Natural language interface for modifications

### Deliverables

#### Week 16-17: AI Core
- [ ] Integrate LLM API (OpenAI/Anthropic/local models)
- [ ] Context-aware prompt engineering for quantum chemistry
- [ ] Input file analysis and optimization suggestions
- [ ] Automatic parameter tuning recommendations

#### Week 18-19: Smart Features
- [ ] Natural language to input file generation
- [ ] "Fix this structure" command
- [ ] Explain input parameters in plain English
- [ ] Predict calculation outcomes (time, resources)
- [ ] AI-assisted debugging of failed calculations

### Technology Stack
- **AI Backend**: OpenAI API / Ollama (local)
- **Context**: RAG with quantum chemistry documentation
- **UI**: Inline suggestions + chat sidebar

### Success Metrics
- AI suggestions accepted by users in ≥ 70% of cases
- 50% reduction in input file creation time
- User satisfaction score ≥ 4.5/5

---

## Phase 5: Advanced Features & Polish (Weeks 20-23)

### Objectives
- Performance optimization
- User experience refinement
- Community features

### Deliverables

#### Week 20-21: Performance
- [ ] Lazy loading for large structures
- [ ] WebWorker for heavy computations
- [ ] Caching and incremental parsing
- [ ] Memory optimization for 10k+ atom systems

#### Week 22-23: Community & Documentation
- [ ] Extension marketplace listing
- [ ] Comprehensive documentation website
- [ ] Video tutorials and quick-start guide
- [ ] Example gallery with common workflows
- [ ] Community templates repository

### Success Metrics
- 10,000+ marketplace downloads
- 4.5+ star rating
- Active community (GitHub Discussions, Discord)

---

## Phase 6: Enterprise & Ecosystem (Future)

### Potential Features
- [ ] HPC cluster integration (Slurm, PBS job submission)
- [ ] Collaborative editing (real-time sync)
- [ ] Database integration (Materials Project, NOMAD)
- [ ] Custom plugin API for proprietary formats
- [ ] Enterprise support and training

---

## Technology Stack Summary

### Core Technologies
- **Extension**: TypeScript + VSCode Extension API
- **Frontend**: React + Tailwind CSS
- **3D Graphics**: Three.js / NGL Viewer
- **Parsing**: Chevrotain (custom DSL)

### Backend & Integrations
- **Python Tools**: dpdata, ASE, PyMOL
- **Conversion**: dpdata + custom adapters
- **AI**: OpenAI API / Ollama / Local models
- **GitHub API**: For dynamic LSP discovery

### Testing & Quality
- **Unit Tests**: Jest (TS), pytest (Python)
- **E2E Tests**: VSCode extension testing framework
- **CI/CD**: GitHub Actions
- **Code Quality**: ESLint, Prettier, Black, MyPy

### Documentation
- **Framework**: Vitepress or Docusaurus
- **API Docs**: TypeDoc (TS), Sphinx (Python)
- **Examples**: Interactive playground

---

## Integration Plans

### dpdata Integration
- **Purpose**: Format conversion and structure manipulation
- **Method**: Python subprocess or native port
- **Timeline**: Phase 3
- **Dependencies**: Python 3.8+

### ASE Integration
- **Purpose**: Advanced structure operations, trajectory analysis
- **Method**: Python backend + VSCode tasks
- **Timeline**: Phase 3
- **Dependencies**: ASE package

### PyMOL Integration
- **Purpose**: Publication-quality molecular images
- **Method**: External process or PyMOL API
- **Timeline**: Phase 3
- **Dependencies**: PyMOL installation (optional)

### AI Integration
- **Purpose**: Smart suggestions and natural language interface
- **Method**: API calls to LLM providers or local models
- **Timeline**: Phase 4
- **Dependencies**: API keys or local model setup

### Dynamic LSP Discovery
- **Purpose**: Auto-discover and integrate new LSP servers from GitHub org
- **Method**: GitHub API + dynamic configuration
- **Timeline**: Phase 2.5
- **Dependencies**: GitHub API access

---

## Milestones & Deliverables

| Phase | Duration | Key Deliverable | Target Date |
|-------|----------|----------------|-------------|
| Phase 1 | 4 weeks | Basic parsing & highlighting | Week 4 |
| Phase 2 | 4 weeks | 3D visualization working | Week 8 |
| Phase 2.5 | 2 weeks | Dynamic LSP discovery | Week 9 |
| Phase 3 | 6 weeks | ASE integration complete | Week 15 |
| Phase 4 | 4 weeks | AI features live | Week 19 |
| Phase 5 | 4 weeks | Production-ready release | Week 23 |
| Phase 6 | Ongoing | Enterprise features | Future |

---

## Risk Assessment

### Technical Risks
- **Large file performance** - Mitigate with lazy loading and WebWorkers
- **3D rendering complexity** - Use established libraries (Three.js, NGL)
- **Format compatibility** - Extensive test suite with real-world files
- **GitHub API rate limits** - Implement caching and graceful degradation

### Resource Risks
- **Development time** - Prioritize core features, defer nice-to-haves
- **AI API costs** - Offer local model alternatives

### Adoption Risks
- **User learning curve** - Invest in documentation and tutorials
- **Competition** - Focus on unique AI features and integration depth

---

## Success Criteria

### Technical Excellence
- ✅ Parse all major quantum chemistry formats
- ✅ 60 FPS rendering for 1000+ atom systems
- ✅ 99% conversion accuracy
- ✅ Test coverage ≥ 80%
- ✅ Dynamic LSP discovery working

### User Impact
- ✅ 10,000+ downloads in first year
- ✅ 4.5+ star rating on marketplace
- ✅ Active community with contributions

### Research Impact
- ✅ Cited in research papers
- ✅ Adopted by research groups
- ✅ Positive testimonials from users

---

## Contributing

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - See [LICENSE](LICENSE) for details.

---

**Last Updated**: 2026-03-03
**Version**: 1.1.0
**Status**: In Progress
