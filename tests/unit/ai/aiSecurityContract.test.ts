import * as fs from 'fs';
import * as path from 'path';

describe('AI manifest security contract', () => {
  const root = path.resolve(__dirname, '../../..');
  const manifest = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));

  it('does not contribute a plaintext API-key setting', () => {
    const properties = manifest.contributes.configuration.properties;
    expect(properties).not.toHaveProperty('openqc.ai.apiKey');
  });

  it('contributes and activates SecretStorage set and clear commands', () => {
    const commands = manifest.contributes.commands.map((item: { command: string }) => item.command);
    for (const command of ['openqc.aiSetApiKey', 'openqc.aiClearApiKey']) {
      expect(commands).toContain(command);
      expect(manifest.activationEvents).toContain(`onCommand:${command}`);
    }
  });
});
