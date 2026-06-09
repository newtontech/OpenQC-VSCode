/**
 * VS Code commands for scientific Python bridges.
 * @module commands/scientificBridgeCommands
 */

import * as vscode from 'vscode';

export function registerScientificBridgeCommands(context: vscode.ExtensionContext): void {
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.parseStructurePython', async () => {
      vscode.window.showInformationMessage(
        'Parse structure with Python backend (requires Python bridge)'
      );
    }),
    vscode.commands.registerCommand('openqc.generateSupercell', async () => {
      vscode.window.showInformationMessage('Generate supercell (requires periodic structure)');
    }),
    vscode.commands.registerCommand('openqc.parseCalculationOutput', async () => {
      vscode.window.showInformationMessage('Parse calculation output (requires cclib)');
    }),
    vscode.commands.registerCommand('openqc.showSCFConvergence', async () => {
      vscode.window.showInformationMessage('Show SCF convergence (parse output first)');
    })
  );
}
