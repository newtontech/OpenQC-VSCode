# Frequently Asked Questions (FAQ)

## General

### Q: What quantum chemistry packages are supported?

A: OpenQC-VSCode supports 7 major packages:
- CP2K
- VASP
- Gaussian
- ORCA
- Quantum ESPRESSO
- GAMESS
- NWChem

### Q: Is this extension free?

A: Yes, OpenQC-VSCode is open-source and free to use under the MIT License.

### Q: Do I need to install the quantum chemistry software?

A: The extension itself works without the software installed, but to run the language servers, you need to install the corresponding LSP packages.

## Installation

### Q: How do I install the extension?

A: Search for "OpenQC-VSCode" in the VSCode Marketplace, or install from source:
```bash
git clone https://github.com/newtontech/OpenQC-VSCode.git
cd OpenQC-VSCode
npm install
npx vsce package
code --install-extension openqc-vscode-*.vsix
```

### Q: Do I need to configure anything after installation?

A: The extension works out of the box. Optional configuration for LSP paths can be done in VSCode settings.

## Usage

### Q: How do I open the 3D structure viewer?

A: Use any of these methods:
- Click the structure icon in the editor title bar
- Press `Ctrl+Shift+P` and run "OpenQC: Visualize Molecular Structure"
- Right-click in the editor and select "OpenQC: Visualize Molecular Structure"

### Q: Why isn't my file being recognized?

A: Check that:
1. The file extension is supported
2. The file content matches expected patterns
3. The file is saved (not unsaved/untitled)

### Q: Can I visualize output files?

A: Currently, visualization works best with input files. Output file support is planned for future releases.

## Language Servers

### Q: What are language servers?

A: Language servers provide intelligent features like auto-completion, hover information, and error diagnostics. They run as separate processes.

### Q: How do I install language servers?

A: Each quantum chemistry package has its own LSP:
- CP2K: `cp2k-lsp-enhanced`
- VASP: `vasp-lsp`
- Gaussian: `gaussian-lsp`
- ORCA: `orca-lsp`
- QE: `qe-lsp`
- GAMESS: `gamess-lsp`
- NWChem: `nwchem-lsp`

### Q: How do I configure LSP paths?

A: Open VSCode settings (File > Preferences > Settings) and search for "openqc.lsp". Set the paths for each package.

### Q: The language server isn't starting. What should I do?

A: Check:
1. The LSP executable is installed and in your PATH
2. The path in settings is correct
3. The LSP is enabled in settings
4. Check the Output panel (View > Output) for error messages

## Troubleshooting

### Q: The 3D viewer shows "No atomic coordinates found"

A: Make sure:
1. The file has valid atomic coordinates in the expected format
2. The file type was detected correctly (check the status bar)
3. Try reloading the window (Ctrl+Shift+P > Developer: Reload Window)

### Q: The plot shows no data

A: Currently, plotting works with input files for basic parameters. Full output file parsing is planned for future releases.

### Q: Syntax highlighting isn't working

A: Try:
1. Reload the window (Ctrl+Shift+P > Developer: Reload Window)
2. Check that the file extension is correct
3. Verify the language mode in the status bar

## Development

### Q: How can I contribute?

A: See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines. Areas include:
- Bug fixes
- New features
- Documentation
- Tests
- New quantum chemistry package support

### Q: How do I report a bug?

A: Open an issue on GitHub with:
- Description of the problem
- Steps to reproduce
- Expected vs actual behavior
- VSCode version and extension version
- Sample file (if applicable)

### Q: Can I request a feature?

A: Yes! Open a feature request on GitHub. Check existing issues first to avoid duplicates.

## Future Plans

### Q: What features are planned?

A: See our [roadmap](PLAN.md). Highlights include:
- Format conversion between packages
- Output file parsing and visualization
- Real-time calculation monitoring
- AI-powered features

### Q: Will you support [other quantum chemistry package]?

A: We're open to adding support for additional packages. Open a feature request to discuss.
