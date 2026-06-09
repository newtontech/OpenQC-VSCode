/**
 * Manifest hygiene tests.
 *
 * Ensures package.json has no duplicate command contributions and
 * that activation events are properly scoped.
 */

import * as fs from 'fs';
import * as path from 'path';

interface PackageJson {
  contributes?: {
    commands?: Array<{ command: string }>;
  };
  activationEvents?: string[];
}

const pkgPath = path.resolve(__dirname, '../../package.json');
const pkg: PackageJson = JSON.parse(fs.readFileSync(pkgPath, 'utf-8'));

describe('Manifest hygiene', () => {
  test('no duplicate command contributions', () => {
    const commands = pkg.contributes?.commands ?? [];
    const ids = commands.map(c => c.command);
    const seen = new Set<string>();
    const dupes: string[] = [];

    for (const id of ids) {
      if (seen.has(id)) {
        dupes.push(id);
      }
      seen.add(id);
    }

    expect(dupes).toEqual([]);
  });

  test('activation events include language triggers', () => {
    const events = pkg.activationEvents ?? [];
    const languageEvents = events.filter(e => e.startsWith('onLanguage:'));

    expect(languageEvents.length).toBeGreaterThan(0);
  });

  test('all commands follow naming convention', () => {
    const commands = pkg.contributes?.commands ?? [];

    for (const cmd of commands) {
      expect(cmd.command).toMatch(/^openqc\./);
    }
  });
});
