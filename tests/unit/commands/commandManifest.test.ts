import * as fs from 'fs';
import * as path from 'path';

describe('command manifest contract', () => {
  const repoRoot = path.resolve(__dirname, '../../..');
  const unadvertisedRuntimeCommands = new Set(['openqc.getDslAuthoringContext']);

  it('keeps contributed and user-facing runtime commands synchronized', () => {
    const manifest = JSON.parse(
      fs.readFileSync(path.join(repoRoot, 'package.json'), 'utf8')
    ) as PackageManifest;
    const contributed = new Set(manifest.contributes.commands.map(entry => entry.command));
    const registered = collectRegisteredCommands(path.join(repoRoot, 'src'));

    expect([...contributed].filter(command => command.startsWith('OpenQC.'))).toEqual([]);
    expect([...contributed].filter(command => !registered.has(command))).toEqual([]);
    expect(
      [...registered]
        .filter(command => command.startsWith('openqc.'))
        .filter(command => !unadvertisedRuntimeCommands.has(command))
        .filter(command => !contributed.has(command))
        .sort()
    ).toEqual([]);
  });

  it('activates ASE, calculator, and export command registrations from the extension entrypoint', () => {
    const extensionSource = fs.readFileSync(path.join(repoRoot, 'src/extension.ts'), 'utf8');

    expect(extensionSource).toContain('registerASECommands(context)');
    expect(extensionSource).toContain('registerCalculatorCommands(context)');
    expect(extensionSource).toContain('registerExportCommands(');
    expect(extensionSource).toContain('() => OpenQCViewerPanel.getCurrentStructureData()');
    expect(extensionSource).toContain('() => OpenQCViewerPanel.saveCurrentStructureToSource()');
  });

  it('contributes calculator settings consumed by the ASE calculator wrapper', () => {
    const manifest = JSON.parse(
      fs.readFileSync(path.join(repoRoot, 'package.json'), 'utf8')
    ) as PackageManifest;
    const properties = manifest.contributes.configuration.properties;

    for (const calculator of ['vasp', 'cp2k', 'qe']) {
      expect(hasOwn(properties, `openqc.calculators.${calculator}.command`)).toBe(true);
      expect(hasOwn(properties, `openqc.calculators.${calculator}.defaultParams`)).toBe(true);
    }
  });
});

interface PackageManifest {
  contributes: {
    commands: Array<{ command: string }>;
    configuration: {
      properties: Record<string, unknown>;
    };
  };
}

function collectRegisteredCommands(srcRoot: string): Set<string> {
  const commands = new Set<string>();

  for (const filePath of walkTsFiles(srcRoot)) {
    const source = fs.readFileSync(filePath, 'utf8');
    for (const match of source.matchAll(/registerCommand\(\s*['"]([^'"]+)['"]/g)) {
      commands.add(match[1]);
    }
    for (const match of source.matchAll(
      /registerCommandWithLegacyAlias\(\s*context,\s*['"]([^'"]+)['"]/g
    )) {
      commands.add(match[1]);
    }
  }

  return commands;
}

function walkTsFiles(root: string): string[] {
  const files: string[] = [];
  for (const entry of fs.readdirSync(root, { withFileTypes: true })) {
    const fullPath = path.join(root, entry.name);
    if (entry.isDirectory()) {
      files.push(...walkTsFiles(fullPath));
    } else if (entry.isFile() && entry.name.endsWith('.ts')) {
      files.push(fullPath);
    }
  }
  return files;
}

function hasOwn(object: Record<string, unknown>, key: string): boolean {
  return Object.prototype.hasOwnProperty.call(object, key);
}
