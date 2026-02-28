/**
 * Document Symbol Provider - Outline view for QC input files
 */

import * as vscode from 'vscode';

export class QCDocumentSymbolProvider implements vscode.DocumentSymbolProvider {
    constructor(private languageId: string) {}
    
    public provideDocumentSymbols(
        document: vscode.TextDocument,
        token: vscode.CancellationToken
    ): vscode.ProviderResult<vscode.SymbolInformation[] | vscode.DocumentSymbol[]> {
        const symbols: vscode.DocumentSymbol[] = [];
        const text = document.getText();
        const lines = text.split('\n');
        
        switch (this.languageId) {
            case 'gaussian':
                return this.parseGaussianSymbols(document, lines);
            case 'vasp':
                return this.parseVASPSymbols(document, lines);
            case 'quantumespresso':
                return this.parseQESymbols(document, lines);
            case 'orca':
                return this.parseORCASymbols(document, lines);
            default:
                return symbols;
        }
    }
    
    private parseGaussianSymbols(
        document: vscode.TextDocument,
        lines: string[]
    ): vscode.DocumentSymbol[] {
        const symbols: vscode.DocumentSymbol[] = [];
        
        let routeStart = -1;
        let titleStart = -1;
        let chargeStart = -1;
        let geometryStart = -1;
        
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];
            
            // Route section (starts with #)
            if (line.trim().startsWith('#') && routeStart === -1) {
                routeStart = i;
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    'Route Section',
                    line.trim(),
                    vscode.SymbolKind.Property,
                    range,
                    range
                ));
            }
            
            // Title line (after route)
            if (routeStart !== -1 && i > routeStart && line.trim() && titleStart === -1) {
                titleStart = i;
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    'Title',
                    line.trim(),
                    vscode.SymbolKind.String,
                    range,
                    range
                ));
            }
            
            // Charge/multiplicity
            if (titleStart !== -1 && i > titleStart && line.match(/^\s*[\d-]+\s+\d+\s*$/)) {
                chargeStart = i;
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    'Charge/Multiplicity',
                    line.trim(),
                    vscode.SymbolKind.Number,
                    range,
                    range
                ));
            }
            
            // Geometry
            if (chargeStart !== -1 && i > chargeStart && line.trim().match(/^[A-Za-z]+\s+/)) {
                if (geometryStart === -1) {
                    geometryStart = i;
                }
            }
        }
        
        // Add geometry section if found
        if (geometryStart !== -1) {
            let geometryEnd = lines.length - 1;
            for (let i = geometryStart; i < lines.length; i++) {
                if (lines[i].trim() === '' || lines[i].trim().match(/^--$/)) {
                    geometryEnd = i - 1;
                    break;
                }
            }
            const range = new vscode.Range(geometryStart, 0, geometryEnd, lines[geometryEnd].length);
            symbols.push(new vscode.DocumentSymbol(
                'Molecular Geometry',
                `${geometryEnd - geometryStart + 1} atoms`,
                vscode.SymbolKind.Struct,
                range,
                range
            ));
        }
        
        return symbols;
    }
    
    private parseVASPSymbols(
        document: vscode.TextDocument,
        lines: string[]
    ): vscode.DocumentSymbol[] {
        const symbols: vscode.DocumentSymbol[] = [];
        
        if (lines.length > 0) {
            // Comment line
            const range0 = new vscode.Range(0, 0, 0, lines[0].length);
            symbols.push(new vscode.DocumentSymbol(
                'Comment',
                lines[0].trim(),
                vscode.SymbolKind.String,
                range0,
                range0
            ));
        }
        
        if (lines.length > 1) {
            // Scaling factor
            const range1 = new vscode.Range(1, 0, 1, lines[1].length);
            symbols.push(new vscode.DocumentSymbol(
                'Scaling Factor',
                lines[1].trim(),
                vscode.SymbolKind.Number,
                range1,
                range1
            ));
        }
        
        if (lines.length > 4) {
            // Lattice vectors
            for (let i = 2; i < 5; i++) {
                const range = new vscode.Range(i, 0, i, lines[i].length);
                symbols.push(new vscode.DocumentSymbol(
                    `Lattice Vector ${i - 1}`,
                    lines[i].trim(),
                    vscode.SymbolKind.Array,
                    range,
                    range
                ));
            }
        }
        
        return symbols;
    }
    
    private parseQESymbols(
        document: vscode.TextDocument,
        lines: string[]
    ): vscode.DocumentSymbol[] {
        const symbols: vscode.DocumentSymbol[] = [];
        
        let currentSection = '';
        let sectionStart = 0;
        
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();
            
            // Check for namelist
            const namelistMatch = line.match(/^&(\w+)/);
            if (namelistMatch) {
                currentSection = namelistMatch[1].toUpperCase();
                sectionStart = i;
                
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    `${currentSection} Namelist`,
                    '',
                    vscode.SymbolKind.Namespace,
                    range,
                    range
                ));
            }
            
            // Check for cards
            const cardMatch = line.match(/^(ATOMIC_POSITIONS|CELL_PARAMETERS|K_POINTS|ATOMIC_SPECIES)/i);
            if (cardMatch) {
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    cardMatch[1].toUpperCase(),
                    '',
                    vscode.SymbolKind.Struct,
                    range,
                    range
                ));
            }
        }
        
        return symbols;
    }
    
    private parseORCASymbols(
        document: vscode.TextDocument,
        lines: string[]
    ): vscode.DocumentSymbol[] {
        const symbols: vscode.DocumentSymbol[] = [];
        
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();
            
            // Simple input line
            if (line.startsWith('!')) {
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    'Simple Input',
                    line.substring(1).trim(),
                    vscode.SymbolKind.Property,
                    range,
                    range
                ));
            }
            
            // Block
            const blockMatch = line.match(/^%(\w+)/);
            if (blockMatch) {
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    blockMatch[1].toUpperCase() + ' Block',
                    '',
                    vscode.SymbolKind.Class,
                    range,
                    range
                ));
            }
            
            // Coordinates
            if (line.match(/^\*\s*xyz/)) {
                const range = new vscode.Range(i, 0, i, line.length);
                symbols.push(new vscode.DocumentSymbol(
                    'Coordinates',
                    '',
                    vscode.SymbolKind.Struct,
                    range,
                    range
                ));
            }
        }
        
        return symbols;
    }
}
