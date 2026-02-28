/**
 * Hover Provider - Show documentation on hover
 */

import * as vscode from 'vscode';

const documentation: Record<string, Record<string, string>> = {
    gaussian: {
        'B3LYP': '**B3LYP Functional**\n\nBecke 3-parameter hybrid with Lee-Yang-Parr correlation.\n\n- Type: Hybrid GGA\n- Recommended for: General organic molecules\n- Reference: Becke, JCP 98, 5648 (1993)',
        'opt': '**Geometry Optimization**\n\nOptimizes molecular geometry to find local minimum.\n\nOptions:\n- `opt=ts` - Transition state search\n- `opt=calcfc` - Calculate initial Hessian\n- `opt=modredundant` - Constrained optimization',
        'freq': '**Frequency Calculation**\n\nCalculates vibrational frequencies and thermochemistry.\n\n- Confirms minimum (all positive frequencies)\n- Provides ZPE, enthalpy, entropy\n- Required for thermochemical analysis',
        'scrf': '**Solvation Model**\n\nImplicit solvation via PCM or SMD.\n\nCommon options:\n- `pcm` - Polarizable Continuum Model\n- `smd` - Solvation Model based on Density\n- `solvent=water` - Specify solvent'
    },
    vasp: {
        'ENCUT': '**Cutoff Energy**\n\nPlane-wave cutoff energy in eV.\n\n- Default: Usually 1.3x ENMAX from POTCAR\n- Typical values: 400-600 eV\n- Higher = more accurate but slower',
        'IBRION': '**Ionic Relaxation**\n\nControls how ions are moved during relaxation.\n\n- 2: Conjugate gradient (safe)\n- 5: Molecular dynamics\n- 7: Damped MD',
        'ISIF': **Stress/Strain**\n\nControls what degrees of freedom are relaxed.\n\n- 2: Ions only\n- 4: Ions + volume\n- 3: Full (ions + cell shape + volume)'
    },
    orca: {
        '!': '**Simple Input Line**\n\nSpecifies method, basis, and keywords.\n\nExample: `! B3LYP def2-SVP TightSCF Opt`',
        '%pal': '**Parallel Processing**\n\nControls parallel execution.\n\nExample:\n```\n%pal nprocs 8\nend\n```'
    }
};

export class QCHoverProvider implements vscode.HoverProvider {
    constructor(private languageId: string) {}
    
    public provideHover(
        document: vscode.TextDocument,
        position: vscode.Position,
        token: vscode.CancellationToken
    ): vscode.ProviderResult<vscode.Hover> {
        const range = document.getWordRangeAtPosition(position);
        if (!range) {
            return undefined;
        }
        
        const word = document.getText(range);
        const docs = documentation[this.languageId];
        
        if (docs && docs[word]) {
            return new vscode.Hover(
                new vscode.MarkdownString(docs[word]),
                range
            );
        }
        
        return undefined;
    }
}
