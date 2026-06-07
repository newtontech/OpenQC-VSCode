import * as vscode from 'vscode';
import { createParser } from '../../parsers';
import { FileTypeDetector, QuantumChemistrySoftware } from '../../managers/FileTypeDetector';

interface CompletionItem {
  label: string;
  kind: vscode.CompletionItemKind;
  detail?: string;
  documentation?: string;
  insertText?: string;
}

export class CompletionProvider implements vscode.CompletionItemProvider {
  private fileTypeDetector: FileTypeDetector;

  constructor() {
    this.fileTypeDetector = new FileTypeDetector();
  }

  /**
   * Provide context-aware completion items for supported quantum chemistry inputs.
   *
   * Implements the VS Code completion provider contract by detecting the active
   * software format and returning completions for the current line and document.
   *
   * @param document - Document requesting completions.
   * @param position - Cursor position where completion was triggered.
   * @param token - Cancellation token supplied by VS Code.
   * @param context - Completion trigger context supplied by VS Code.
   * @returns Completion items, a completion list, or an empty list when unsupported.
   */
  provideCompletionItems(
    document: vscode.TextDocument,
    position: vscode.Position,
    token: vscode.CancellationToken,
    context: vscode.CompletionContext
  ): vscode.ProviderResult<vscode.CompletionItem[] | vscode.CompletionList> {
    const software = this.fileTypeDetector.detectSoftware(document);
    if (!software) {
      return [];
    }

    const line = document.lineAt(position).text.substring(0, position.character);
    const content = document.getText();

    return this.getCompletions(software, line, content, position);
  }

  private getCompletions(
    software: QuantumChemistrySoftware,
    line: string,
    content: string,
    position: vscode.Position
  ): vscode.CompletionItem[] {
    const completions: vscode.CompletionItem[] = [];

    switch (software) {
      case 'VASP':
        completions.push(...this.getVASPCompletions(line));
        break;
      case 'Gaussian':
        completions.push(...this.getGaussianCompletions(line));
        break;
      case 'ORCA':
        completions.push(...this.getORCACompletions(line));
        break;
      case 'CP2K':
        completions.push(...this.getCP2KCompletions(line, content));
        break;
      case 'Quantum ESPRESSO':
        completions.push(...this.getQECompletions(line, content));
        break;
      case 'GAMESS':
        completions.push(...this.getGAMESSCompletions(line, content));
        break;
      case 'NWChem':
        completions.push(...this.getNWChemCompletions(line, content));
        break;
    }

    return completions;
  }

  private getVASPCompletions(line: string): vscode.CompletionItem[] {
    const items: CompletionItem[] = [
      {
        label: 'ENCUT',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Energy cutoff (eV)',
        documentation: 'Specifies the cutoff energy for the plane wave basis set',
      },
      {
        label: 'PREC',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Precision',
        documentation: 'Determines the precision of the calculation',
      },
      {
        label: 'EDIFF',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Electronic convergence',
        documentation: 'Convergence criterion for electronic relaxation',
      },
      {
        label: 'NELM',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Max electronic steps',
        documentation: 'Maximum number of electronic self-consistency steps',
      },
      {
        label: 'ISMEAR',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Smearing method',
        documentation: 'Determines how the partial occupancies are set for each orbital',
      },
      {
        label: 'SIGMA',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Smearing width',
        documentation: 'The width of the smearing in eV',
      },
      {
        label: 'IBRION',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Ionic relaxation',
        documentation: 'Determines how the ions are updated and moved',
      },
      {
        label: 'NSW',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Max ionic steps',
        documentation: 'Maximum number of ionic steps',
      },
      {
        label: 'ISIF',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Stress/relaxation',
        documentation: 'Determines whether the stress tensor is calculated',
      },
      {
        label: 'ISPIN',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Spin polarization',
        documentation: 'Specifies whether spin-polarized calculation is performed',
      },
      {
        label: 'MAGMOM',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Magnetic moments',
        documentation: 'Initial magnetic moments for each atom',
      },
      {
        label: 'LREAL',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Real space projection',
        documentation: 'Determines whether the projection operators are evaluated in real space',
      },
      {
        label: 'ALGO',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Algorithm',
        documentation: 'Selects the algorithm for electron optimization',
      },
      {
        label: 'LWAVE',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Write WAVECAR',
        documentation: 'Determines whether the wavefunctions are written to WAVECAR',
      },
      {
        label: 'LCHARG',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Write CHGCAR',
        documentation: 'Determines whether the charge densities are written to CHGCAR',
      },
    ];

    return items.map(item => {
      const completion = new vscode.CompletionItem(item.label, item.kind);
      completion.detail = item.detail;
      completion.documentation = item.documentation;
      completion.insertText = new vscode.SnippetString(`${item.label} = \${1:value}`);
      return completion;
    });
  }

  private getGaussianCompletions(line: string): vscode.CompletionItem[] {
    const items: CompletionItem[] = [
      { label: 'HF', kind: vscode.CompletionItemKind.Method, detail: 'Hartree-Fock' },
      { label: 'B3LYP', kind: vscode.CompletionItemKind.Method, detail: 'Hybrid DFT functional' },
      { label: 'PBE', kind: vscode.CompletionItemKind.Method, detail: 'GGA DFT functional' },
      {
        label: 'MP2',
        kind: vscode.CompletionItemKind.Method,
        detail: 'Second-order Møller-Plesset',
      },
      {
        label: 'CCSD',
        kind: vscode.CompletionItemKind.Method,
        detail: 'Coupled Cluster Singles and Doubles',
      },
      { label: 'opt', kind: vscode.CompletionItemKind.Keyword, detail: 'Geometry optimization' },
      { label: 'freq', kind: vscode.CompletionItemKind.Keyword, detail: 'Frequency calculation' },
      { label: 'sp', kind: vscode.CompletionItemKind.Keyword, detail: 'Single point' },
      { label: 'td', kind: vscode.CompletionItemKind.Keyword, detail: 'Time-dependent DFT' },
      { label: '6-31G(d)', kind: vscode.CompletionItemKind.Struct, detail: 'Pople basis set' },
      {
        label: '6-311G(d,p)',
        kind: vscode.CompletionItemKind.Struct,
        detail: 'Pople basis set with polarization',
      },
      { label: 'def2-TZVP', kind: vscode.CompletionItemKind.Struct, detail: 'Ahlrichs basis set' },
      {
        label: 'cc-pVTZ',
        kind: vscode.CompletionItemKind.Struct,
        detail: 'Dunning correlation consistent',
      },
    ];

    return items.map(item => {
      const completion = new vscode.CompletionItem(item.label, item.kind);
      completion.detail = item.detail;
      return completion;
    });
  }

  private getORCACompletions(line: string): vscode.CompletionItem[] {
    const items: CompletionItem[] = [
      { label: 'RHF', kind: vscode.CompletionItemKind.Method, detail: 'Restricted Hartree-Fock' },
      { label: 'UHF', kind: vscode.CompletionItemKind.Method, detail: 'Unrestricted Hartree-Fock' },
      { label: 'DFT', kind: vscode.CompletionItemKind.Method, detail: 'Density Functional Theory' },
      {
        label: 'MP2',
        kind: vscode.CompletionItemKind.Method,
        detail: 'Second-order Møller-Plesset',
      },
      { label: 'CCSD', kind: vscode.CompletionItemKind.Method, detail: 'Coupled Cluster' },
      { label: 'def2-SVP', kind: vscode.CompletionItemKind.Struct, detail: 'Basis set' },
      { label: 'def2-TZVP', kind: vscode.CompletionItemKind.Struct, detail: 'Basis set' },
      { label: 'def2-TZVPP', kind: vscode.CompletionItemKind.Struct, detail: 'Basis set' },
      { label: 'OPT', kind: vscode.CompletionItemKind.Keyword, detail: 'Geometry optimization' },
      { label: 'FREQ', kind: vscode.CompletionItemKind.Keyword, detail: 'Frequency calculation' },
      {
        label: 'TightSCF',
        kind: vscode.CompletionItemKind.Keyword,
        detail: 'Tight SCF convergence',
      },
      { label: 'Grid5', kind: vscode.CompletionItemKind.Keyword, detail: 'Integration grid' },
    ];

    return items.map(item => {
      const completion = new vscode.CompletionItem(item.label, item.kind);
      completion.detail = item.detail;
      return completion;
    });
  }

  private getCP2KCompletions(line: string, content: string): vscode.CompletionItem[] {
    const items: CompletionItem[] = [];

    if (line.trim().startsWith('&')) {
      items.push(
        { label: 'GLOBAL', kind: vscode.CompletionItemKind.Struct, detail: 'Global settings' },
        { label: 'FORCE_EVAL', kind: vscode.CompletionItemKind.Struct, detail: 'Force evaluation' },
        { label: 'DFT', kind: vscode.CompletionItemKind.Struct, detail: 'DFT settings' },
        { label: 'QS', kind: vscode.CompletionItemKind.Struct, detail: 'Quickstep settings' },
        { label: 'SCF', kind: vscode.CompletionItemKind.Struct, detail: 'SCF settings' },
        { label: 'MOTION', kind: vscode.CompletionItemKind.Struct, detail: 'Motion/MD settings' },
        {
          label: 'GEO_OPT',
          kind: vscode.CompletionItemKind.Struct,
          detail: 'Geometry optimization',
        }
      );
    } else {
      items.push(
        { label: 'PROJECT_NAME', kind: vscode.CompletionItemKind.Property, detail: 'Project name' },
        { label: 'RUN_TYPE', kind: vscode.CompletionItemKind.Property, detail: 'Calculation type' },
        {
          label: 'METHOD',
          kind: vscode.CompletionItemKind.Property,
          detail: 'Electronic structure method',
        },
        {
          label: 'BASIS_SET_FILE_NAME',
          kind: vscode.CompletionItemKind.Property,
          detail: 'Basis set file',
        },
        {
          label: 'POTENTIAL_FILE_NAME',
          kind: vscode.CompletionItemKind.Property,
          detail: 'Potential file',
        }
      );
    }

    return items.map(item => {
      const completion = new vscode.CompletionItem(item.label, item.kind);
      completion.detail = item.detail;
      return completion;
    });
  }

  private getQECompletions(line: string, content: string): vscode.CompletionItem[] {
    const items: CompletionItem[] = [
      {
        label: 'calculation',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Type of calculation',
      },
      { label: 'restart_mode', kind: vscode.CompletionItemKind.Property, detail: 'Restart mode' },
      {
        label: 'pseudo_dir',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Pseudopotential directory',
      },
      { label: 'outdir', kind: vscode.CompletionItemKind.Property, detail: 'Output directory' },
      { label: 'ibrav', kind: vscode.CompletionItemKind.Property, detail: 'Bravais lattice index' },
      {
        label: 'ecutwfc',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Wavefunction cutoff (Ry)',
      },
      {
        label: 'ecutrho',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Charge density cutoff (Ry)',
      },
      {
        label: 'conv_thr',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Convergence threshold',
      },
      {
        label: 'diagonalization',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Diagonalization method',
      },
      {
        label: 'mixing_mode',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Charge mixing mode',
      },
    ];

    return items.map(item => {
      const completion = new vscode.CompletionItem(item.label, item.kind);
      completion.detail = item.detail;
      completion.insertText = new vscode.SnippetString(`${item.label} = \${1:value}`);
      return completion;
    });
  }

  private getGAMESSCompletions(line: string, content: string): vscode.CompletionItem[] {
    const items: CompletionItem[] = [
      { label: 'RUNTYP', kind: vscode.CompletionItemKind.Property, detail: 'Run type' },
      {
        label: 'SCFTYP',
        kind: vscode.CompletionItemKind.Property,
        detail: 'SCF type (RHF/UHF/ROHF)',
      },
      { label: 'DFTTYP', kind: vscode.CompletionItemKind.Property, detail: 'DFT functional' },
      { label: 'MAXIT', kind: vscode.CompletionItemKind.Property, detail: 'Maximum iterations' },
      { label: 'CONV', kind: vscode.CompletionItemKind.Property, detail: 'Convergence criterion' },
      { label: 'GBASIS', kind: vscode.CompletionItemKind.Property, detail: 'Gaussian basis set' },
      { label: 'NGAUSS', kind: vscode.CompletionItemKind.Property, detail: 'Number of Gaussians' },
      {
        label: 'NDFUNC',
        kind: vscode.CompletionItemKind.Property,
        detail: 'd functions on heavy atoms',
      },
      { label: 'NPFUNC', kind: vscode.CompletionItemKind.Property, detail: 'p functions on H' },
      {
        label: 'DIFFSP',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Diffuse functions on heavy atoms',
      },
      {
        label: 'DIFFS',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Diffuse functions on H',
      },
    ];

    return items.map(item => {
      const completion = new vscode.CompletionItem(item.label, item.kind);
      completion.detail = item.detail;
      completion.insertText = new vscode.SnippetString(`${item.label}=\${1:value}`);
      return completion;
    });
  }

  private getNWChemCompletions(line: string, content: string): vscode.CompletionItem[] {
    const items: CompletionItem[] = [
      { label: 'geometry', kind: vscode.CompletionItemKind.Struct, detail: 'Geometry block' },
      { label: 'basis', kind: vscode.CompletionItemKind.Struct, detail: 'Basis set block' },
      { label: 'scf', kind: vscode.CompletionItemKind.Struct, detail: 'SCF block' },
      { label: 'dft', kind: vscode.CompletionItemKind.Struct, detail: 'DFT block' },
      { label: 'task', kind: vscode.CompletionItemKind.Keyword, detail: 'Task directive' },
      {
        label: 'xc',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Exchange-correlation functional',
      },
      {
        label: 'convergence',
        kind: vscode.CompletionItemKind.Property,
        detail: 'Convergence criteria',
      },
      { label: 'maxiter', kind: vscode.CompletionItemKind.Property, detail: 'Maximum iterations' },
      { label: 'thresh', kind: vscode.CompletionItemKind.Property, detail: 'Threshold' },
    ];

    return items.map(item => {
      const completion = new vscode.CompletionItem(item.label, item.kind);
      completion.detail = item.detail;
      return completion;
    });
  }
}
