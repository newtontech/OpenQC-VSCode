import * as vscode from 'vscode';
import { FileTypeDetector, QuantumChemistrySoftware } from '../../managers/FileTypeDetector';

interface ParameterInfo {
  name: string;
  description: string;
  defaultValue?: string;
  type?: string;
  allowedValues?: string[];
}

export class HoverProvider implements vscode.HoverProvider {
  private fileTypeDetector: FileTypeDetector;

  private parameterDocs: Partial<Record<QuantumChemistrySoftware, Record<string, ParameterInfo>>> =
    {
      VASP: {
        ENCUT: {
          name: 'ENCUT',
          description: 'Specifies the cutoff energy for the plane wave basis set in eV.',
          defaultValue: 'largest ENMAX from POTCAR',
          type: 'real',
        },
        PREC: {
          name: 'PREC',
          description: 'Determines the precision of the calculation.',
          allowedValues: ['Low', 'Normal', 'Accurate', 'Single'],
          defaultValue: 'Normal',
        },
        EDIFF: {
          name: 'EDIFF',
          description: 'Convergence criterion for electronic relaxation.',
          defaultValue: '1E-4',
          type: 'real',
        },
        ISMEAR: {
          name: 'ISMEAR',
          description: 'Determines how the partial occupancies are set for each orbital.',
          allowedValues: ['-5', '-4', '-3', '-2', '-1', '0', '1', 'N'],
          defaultValue: '1 (metals) or 0 (insulators)',
        },
        IBRION: {
          name: 'IBRION',
          description: 'Determines how the ions are updated and moved.',
          allowedValues: ['-1', '0', '1', '2', '3', '5', '6', '7', '8'],
          defaultValue: '-1 (no update)',
        },
        NSW: {
          name: 'NSW',
          description: 'Maximum number of ionic steps.',
          defaultValue: '0',
          type: 'integer',
        },
        ISIF: {
          name: 'ISIF',
          description:
            'Determines whether the stress tensor is calculated and which degrees of freedom are allowed to change.',
          allowedValues: ['0', '1', '2', '3', '4', '5', '6', '7'],
          defaultValue: '2',
        },
      },
      Gaussian: {
        HF: {
          name: 'HF',
          description: 'Hartree-Fock method. Performs a single-determinant SCF calculation.',
        },
        B3LYP: {
          name: 'B3LYP',
          description:
            'Hybrid DFT functional combining Becke 3-parameter exchange with Lee-Yang-Parr correlation.',
        },
        opt: {
          name: 'opt',
          description: 'Requests geometry optimization to a stationary point.',
        },
        freq: {
          name: 'freq',
          description: 'Requests vibrational frequency and thermochemical analysis.',
        },
      },
      ORCA: {
        RHF: {
          name: 'RHF',
          description: 'Restricted Hartree-Fock for closed-shell systems.',
        },
        UHF: {
          name: 'UHF',
          description: 'Unrestricted Hartree-Fock for open-shell systems.',
        },
        'def2-SVP': {
          name: 'def2-SVP',
          description: 'Split-valence polarized basis set of Ahlrichs.',
        },
        'def2-TZVP': {
          name: 'def2-TZVP',
          description: 'Triple-zeta valence polarized basis set of Ahlrichs.',
        },
      },
      CP2K: {
        RUN_TYPE: {
          name: 'RUN_TYPE',
          description: 'Specifies the type of calculation to perform.',
          allowedValues: ['ENERGY', 'GEO_OPT', 'MD', 'BAND', 'MONTECARLO', 'EP'],
          defaultValue: 'ENERGY',
        },
        PROJECT_NAME: {
          name: 'PROJECT_NAME',
          description: 'Name of the project, used for output files.',
          defaultValue: 'PROJECT',
        },
      },
      'Quantum ESPRESSO': {
        calculation: {
          name: 'calculation',
          description: 'Type of calculation to be performed.',
          allowedValues: ['scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md'],
          defaultValue: 'scf',
        },
        ecutwfc: {
          name: 'ecutwfc',
          description: 'Kinetic energy cutoff for wavefunctions in Ry.',
          type: 'real',
        },
        ibrav: {
          name: 'ibrav',
          description: 'Bravais lattice index.',
          allowedValues: ['0-14'],
        },
      },
      GAMESS: {
        RUNTYP: {
          name: 'RUNTYP',
          description: 'Type of calculation to perform.',
          allowedValues: ['ENERGY', 'OPTIMIZE', 'SADPOINT', 'HESSIAN', 'IRC'],
          defaultValue: 'ENERGY',
        },
        SCFTYP: {
          name: 'SCFTYP',
          description: 'Type of SCF calculation.',
          allowedValues: ['RHF', 'UHF', 'ROHF', 'GVB', 'MCSCF'],
          defaultValue: 'RHF',
        },
      },
      NWChem: {
        task: {
          name: 'task',
          description: 'Specifies the theory and operation to perform.',
          allowedValues: ['scf', 'dft', 'mp2', 'ccsd', 'tce'],
        },
        xc: {
          name: 'xc',
          description: 'Exchange-correlation functional for DFT.',
          allowedValues: ['b3lyp', 'pbe', 'blyp', 'bp86'],
        },
      },
    };

  constructor() {
    this.fileTypeDetector = new FileTypeDetector();
  }

  /**
   * Provide hover documentation for recognized quantum chemistry parameters.
   *
   * @param document - Document containing the hover target.
   * @param position - Cursor position for the hover request.
   * @param token - Cancellation token supplied by VS Code.
   * @returns Hover markdown for a known parameter, or null when no match is found.
   */
  provideHover(
    document: vscode.TextDocument,
    position: vscode.Position,
    token: vscode.CancellationToken
  ): vscode.ProviderResult<vscode.Hover> {
    const software = this.fileTypeDetector.detectSoftware(document);
    if (!software) {
      return null;
    }

    const wordRange = document.getWordRangeAtPosition(position);
    if (!wordRange) {
      return null;
    }

    const word = document.getText(wordRange);
    const docs = this.parameterDocs[software];

    if (!docs || !docs[word]) {
      return null;
    }

    const info = docs[word];
    const contents = this.formatHoverContent(info);

    return new vscode.Hover(contents, wordRange);
  }

  private formatHoverContent(info: ParameterInfo): vscode.MarkdownString {
    let content = `**${info.name}**\n\n`;
    content += `${info.description}\n\n`;

    if (info.defaultValue) {
      content += `*Default: ${info.defaultValue}*\n\n`;
    }

    if (info.type) {
      content += `*Type: ${info.type}*\n\n`;
    }

    if (info.allowedValues) {
      content += `*Allowed values: ${info.allowedValues.join(', ')}*\n\n`;
    }

    const markdown = new vscode.MarkdownString(content);
    markdown.isTrusted = true;
    return markdown;
  }
}
