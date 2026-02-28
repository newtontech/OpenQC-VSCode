/**
 * Completion Provider - Auto-completion for quantum chemistry keywords
 */

import * as vscode from 'vscode';

interface CompletionItem {
    label: string;
    kind: vscode.CompletionItemKind;
    documentation: string;
    insertText: string;
}

// Gaussian route section keywords
const gaussianKeywords: CompletionItem[] = [
    // Functionals
    { label: 'B3LYP', kind: vscode.CompletionItemKind.Function, documentation: 'Becke, 3-parameter, Lee-Yang-Parr hybrid functional', insertText: 'B3LYP' },
    { label: 'PBE0', kind: vscode.CompletionItemKind.Function, documentation: 'Perdew-Burke-Ernzerhof hybrid functional', insertText: 'PBE0' },
    { label: 'M06-2X', kind: vscode.CompletionItemKind.Function, documentation: 'Minnesota 2006 meta-hybrid functional', insertText: 'M06-2X' },
    { label: 'wB97X-D', kind: vscode.CompletionItemKind.Function, documentation: 'Long-range corrected hybrid with dispersion', insertText: 'wB97X-D' },
    { label: 'CAM-B3LYP', kind: vscode.CompletionItemKind.Function, documentation: 'Coulomb attenuated B3LYP', insertText: 'CAM-B3LYP' },
    { label: 'B2PLYP', kind: vscode.CompletionItemKind.Function, documentation: 'Double hybrid functional', insertText: 'B2PLYP' },
    
    // Basis sets
    { label: '6-31G(d)', kind: vscode.CompletionItemKind.Constant, documentation: 'Pople basis set with polarization on heavy atoms', insertText: '6-31G(d)' },
    { label: '6-311+G(d,p)', kind: vscode.CompletionItemKind.Constant, documentation: 'Extended Pople basis with diffuse and polarization', insertText: '6-311+G(d,p)' },
    { label: 'def2-SVP', kind: vscode.CompletionItemKind.Constant, documentation: 'Ahlrichs split-valence basis', insertText: 'def2-SVP' },
    { label: 'def2-TZVP', kind: vscode.CompletionItemKind.Constant, documentation: 'Ahlrichs triple-zeta basis', insertText: 'def2-TZVP' },
    { label: 'cc-pVDZ', kind: vscode.CompletionItemKind.Constant, documentation: 'Dunning correlation-consistent double-zeta', insertText: 'cc-pVDZ' },
    { label: 'cc-pVTZ', kind: vscode.CompletionItemKind.Constant, documentation: 'Dunning correlation-consistent triple-zeta', insertText: 'cc-pVTZ' },
    { label: 'aug-cc-pVTZ', kind: vscode.CompletionItemKind.Constant, documentation: 'Dunning basis with diffuse functions', insertText: 'aug-cc-pVTZ' },
    
    // Job types
    { label: 'opt', kind: vscode.CompletionItemKind.Method, documentation: 'Geometry optimization', insertText: 'opt' },
    { label: 'freq', kind: vscode.CompletionItemKind.Method, documentation: 'Frequency calculation', insertText: 'freq' },
    { label: 'sp', kind: vscode.CompletionItemKind.Method, documentation: 'Single point energy', insertText: 'sp' },
    { label: 'td', kind: vscode.CompletionItemKind.Method, documentation: 'Time-dependent DFT', insertText: 'td' },
    { label: 'irc', kind: vscode.CompletionItemKind.Method, documentation: 'Intrinsic reaction coordinate', insertText: 'irc' },
    { label: 'scan', kind: vscode.CompletionItemKind.Method, documentation: 'Potential energy scan', insertText: 'scan' },
    { label: 'nmr', kind: vscode.CompletionItemKind.Method, documentation: 'NMR calculation', insertText: 'nmr' },
    
    // Additional keywords
    { label: 'empiricaldispersion=GD3BJ', kind: vscode.CompletionItemKind.Property, documentation: 'Grimme D3 dispersion with Becke-Johnson damping', insertText: 'empiricaldispersion=GD3BJ' },
    { label: 'scrf=(pcm,solvent=water)', kind: vscode.CompletionItemKind.Property, documentation: 'Polarizable continuum model solvation', insertText: 'scrf=(pcm,solvent=water)' },
    { label: 'int=ultrafine', kind: vscode.CompletionItemKind.Property, documentation: 'Ultrafine integration grid', insertText: 'int=ultrafine' },
    { label: 'scf=tight', kind: vscode.CompletionItemKind.Property, documentation: 'Tight SCF convergence', insertText: 'scf=tight' },
    { label: 'pop=full', kind: vscode.CompletionItemKind.Property, documentation: 'Full population analysis', insertText: 'pop=full' },
    { label: 'output=wfn', kind: vscode.CompletionItemKind.Property, documentation: 'Output wavefunction file', insertText: 'output=wfn' }
];

// VASP keywords
const vaspKeywords: CompletionItem[] = [
    { label: 'SYSTEM', kind: vscode.CompletionItemKind.Property, documentation: 'System description', insertText: 'SYSTEM = ' },
    { label: 'ENCUT', kind: vscode.CompletionItemKind.Property, documentation: 'Cutoff energy (eV)', insertText: 'ENCUT = 520' },
    { label: 'EDIFF', kind: vscode.CompletionItemKind.Property, documentation: 'Electronic convergence', insertText: 'EDIFF = 1E-6' },
    { label: 'EDIFFG', kind: vscode.CompletionItemKind.Property, documentation: 'Ionic convergence', insertText: 'EDIFFG = -0.01' },
    { label: 'IBRION', kind: vscode.CompletionItemKind.Property, documentation: 'Ionic relaxation method', insertText: 'IBRION = 2' },
    { label: 'ISIF', kind: vscode.CompletionItemKind.Property, documentation: 'Stress/strain relaxation', insertText: 'ISIF = 3' },
    { label: 'NSW', kind: vscode.CompletionItemKind.Property, documentation: 'Number of ionic steps', insertText: 'NSW = 100' },
    { label: 'ISPIN', kind: vscode.CompletionItemKind.Property, documentation: 'Spin polarization', insertText: 'ISPIN = 2' },
    { label: 'LCHARG', kind: vscode.CompletionItemKind.Property, documentation: 'Write charge density', insertText: 'LCHARG = .TRUE.' },
    { label: 'LWAVE', kind: vscode.CompletionItemKind.Property, documentation: 'Write wavefunction', insertText: 'LWAVE = .FALSE.' }
];

// ORCA keywords
const orcaKeywords: CompletionItem[] = [
    { label: '! B3LYP def2-SVP', kind: vscode.CompletionItemKind.Property, documentation: 'Simple method line', insertText: '! B3LYP def2-SVP' },
    { label: '%pal nprocs 8 end', kind: vscode.CompletionItemKind.Property, documentation: 'Parallel processing', insertText: '%pal nprocs 8 end' },
    { label: '%maxcore 4000', kind: vscode.CompletionItemKind.Property, documentation: 'Memory per core (MB)', insertText: '%maxcore 4000' },
    { label: '%scf convergence tight end', kind: vscode.CompletionItemKind.Property, documentation: 'Tight SCF convergence', insertText: '%scf\n  convergence tight\nend' },
    { label: '%geom constraints end', kind: vscode.CompletionItemKind.Property, documentation: 'Geometry constraints', insertText: '%geom\n  constraints\n    { ... }\n  end\nend' }
];

export class GaussianCompletionProvider implements vscode.CompletionItemProvider {
    public provideCompletionItems(
        document: vscode.TextDocument,
        position: vscode.Position,
        token: vscode.CancellationToken,
        context: vscode.CompletionContext
    ): vscode.ProviderResult<vscode.CompletionItem[] | vscode.CompletionList> {
        const items: vscode.CompletionItem[] = [];
        const languageId = document.languageId;
        
        let keywords: CompletionItem[] = [];
        
        switch (languageId) {
            case 'gaussian':
                keywords = gaussianKeywords;
                break;
            case 'vasp':
                keywords = vaspKeywords;
                break;
            case 'orca':
                keywords = orcaKeywords;
                break;
            default:
                keywords = gaussianKeywords;
        }
        
        for (const kw of keywords) {
            const item = new vscode.CompletionItem(kw.label, kw.kind);
            item.documentation = new vscode.MarkdownString(kw.documentation);
            item.insertText = kw.insertText;
            items.push(item);
        }
        
        return items;
    }
}
