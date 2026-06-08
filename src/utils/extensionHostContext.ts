import * as vscode from 'vscode';

export interface ExtensionHostContext {
  isRemote: boolean;
  remoteName: string | undefined;
  label: string;
  commandResolution: string;
}

const REMOTE_LABELS: Record<string, string> = {
  'ssh-remote': 'Remote SSH',
  wsl: 'WSL',
  'dev-container': 'Dev Container',
  'attached-container': 'Dev Container',
  codespaces: 'Codespaces',
};

export function getExtensionHostContext(
  env: Pick<typeof vscode.env, 'remoteName'> = vscode.env
): ExtensionHostContext {
  const remoteName = env.remoteName;
  if (!remoteName) {
    return {
      isRemote: false,
      remoteName: undefined,
      label: 'local extension host',
      commandResolution: 'LSP commands and paths resolve on this machine.',
    };
  }

  const remoteLabel = REMOTE_LABELS[remoteName] || remoteName;
  return {
    isRemote: true,
    remoteName,
    label: `remote extension host (${remoteLabel})`,
    commandResolution: 'LSP commands and paths resolve in the remote workspace environment.',
  };
}

export function describeExtensionHostContext(context = getExtensionHostContext()): string {
  if (!context.isRemote) {
    return `${context.label}; ${context.commandResolution}`;
  }

  return `${context.label}, remoteName=${context.remoteName}; ${context.commandResolution}`;
}
