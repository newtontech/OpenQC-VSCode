/**
 * Export Commands - VSCode commands for exporting structures
 *
 * Phase 2: Week 7-8 - Export modified structures
 *
 * Provides VSCode commands to export molecular structures that have been
 * modified in the 3D viewer back to quantum chemistry input file formats.
 */

import * as vscode from 'vscode';
import { StructureExporter, ExportFormat, ExportOptions } from '../utils/structureExporter';
import { Atom } from '../visualizers/types';

/**
 * Register export commands
 */
export function registerExportCommands(
  context: vscode.ExtensionContext,
  getStructureData?: () => {
    atoms: Atom[];
    cell?: { a: number; b: number; c: number; alpha: number; beta: number; gamma: number };
    pbc?: [boolean, boolean, boolean];
  }
): void {
  const exporter = new StructureExporter(context);

  // Export current structure command
  const exportStructureCommand = vscode.commands.registerCommand(
    'openqc.exportStructure',
    async () => {
      if (!getStructureData) {
        vscode.window.showErrorMessage('No structure data available');
        return;
      }

      const structureData = getStructureData();
      if (!structureData || structureData.atoms.length === 0) {
        vscode.window.showErrorMessage('No atoms to export');
        return;
      }

      // Show format picker
      const formats = StructureExporter.getSupportedFormats();
      const formatItems = formats.map(format => ({
        label: StructureExporter.getFormatDisplayName(format),
        description: format,
      }));

      const selected = await vscode.window.showQuickPick(formatItems, {
        placeHolder: 'Select export format',
      });

      if (!selected) {
        return;
      }

      // Export structure
      const result = await exporter.exportStructure(structureData, {
        format: selected.description as ExportFormat,
      });

      if (!result.success) {
        vscode.window.showErrorMessage(`Export failed: ${result.error}`);
      }
    }
  );

  // Export with format picker
  const exportWithPickerCommand = vscode.commands.registerCommand(
    'openqc.exportStructureWithPicker',
    async () => {
      if (!getStructureData) {
        vscode.window.showErrorMessage('No structure data available');
        return;
      }

      const structureData = getStructureData();
      if (!structureData || structureData.atoms.length === 0) {
        vscode.window.showErrorMessage('No atoms to export');
        return;
      }

      // Show format picker with recent formats
      const formats = StructureExporter.getSupportedFormats();
      const formatItems = formats.map(format => ({
        label: StructureExporter.getFormatDisplayName(format),
        description: format,
      }));

      const selected = await vscode.window.showQuickPick(formatItems, {
        placeHolder: 'Select export format',
      });

      if (!selected) {
        return;
      }

      // Prompt for output path
      const extensions = getFormatExtensions(selected.description as ExportFormat);
      const defaultFilename = `structure.${extensions[0]}`;

      const uri = await vscode.window.showSaveDialog({
        defaultUri: vscode.Uri.file(defaultFilename),
        saveLabel: 'Export Structure',
        filters: {
          [selected.label]: extensions,
        },
      });

      if (!uri) {
        return;
      }

      // Export structure
      const result = await exporter.exportStructure(structureData, {
        format: selected.description as ExportFormat,
        outputPath: uri.fsPath,
      });

      if (!result.success) {
        vscode.window.showErrorMessage(`Export failed: ${result.error}`);
      }
    }
  );

  context.subscriptions.push(
    exportStructureCommand,
    exportWithPickerCommand
  );
}

/**
 * Capitalize first letter
 */
function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}

/**
 * Get file extensions for format
 */
function getFormatExtensions(format: string): string[] {
  const extensionMap: Record<string, string[]> = {
    vasp: ['POSCAR', 'CONTCAR'],
    cp2k: ['inp'],
    qe: ['in', 'pw.in'],
    gaussian: ['com', 'gjf'],
    orca: ['inp'],
    nwchem: ['nw', 'nwinp'],
    gamess: ['inp'],
    lammps: ['lmp', 'data'],
    xyz: ['xyz'],
    extxyz: ['extxyz'],
    pdb: ['pdb'],
    cif: ['cif'],
  };
  return extensionMap[format] || ['txt'];
}
