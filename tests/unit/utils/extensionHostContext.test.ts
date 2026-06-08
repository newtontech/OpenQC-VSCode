import {
  describeExtensionHostContext,
  getExtensionHostContext,
} from '../../../src/utils/extensionHostContext';

jest.mock('vscode', () => ({
  env: {
    remoteName: undefined,
  },
}));

describe('extensionHostContext', () => {
  it('describes local extension host command resolution', () => {
    const context = getExtensionHostContext({ remoteName: undefined });

    expect(context).toEqual({
      isRemote: false,
      remoteName: undefined,
      label: 'local extension host',
      commandResolution: 'LSP commands and paths resolve on this machine.',
    });
    expect(describeExtensionHostContext(context)).toContain('local extension host');
  });

  it('labels common remote extension host names', () => {
    const context = getExtensionHostContext({ remoteName: 'ssh-remote' });

    expect(context.isRemote).toBe(true);
    expect(context.remoteName).toBe('ssh-remote');
    expect(context.label).toBe('remote extension host (Remote SSH)');
    expect(describeExtensionHostContext(context)).toContain('remoteName=ssh-remote');
    expect(describeExtensionHostContext(context)).toContain('remote workspace environment');
  });

  it('falls back to the raw remoteName for unknown remote hosts', () => {
    const context = getExtensionHostContext({ remoteName: 'custom-remote' });

    expect(context.label).toBe('remote extension host (custom-remote)');
  });
});
