import * as vscode from 'vscode';

export type LspLifecycleState = 'idle' | 'starting' | 'running' | 'stopped' | 'failed';

export interface LspStatusSnapshot {
  key: string;
  state: LspLifecycleState;
  software: string;
  languageId: string;
  workspaceLabel: string;
  workspaceIdentity: string;
  commandSummary: string;
  hostContext: string;
  message?: string;
  updatedAt: Date;
}

export interface LspLifecycleDetails {
  key: string;
  software: string;
  languageId: string;
  workspaceLabel?: string;
  workspaceIdentity: string;
  commandSummary: string;
  hostContext: string;
  message?: string;
}

export class LspStatusService implements vscode.Disposable {
  private readonly output: vscode.OutputChannel;
  private readonly statusBarItem: vscode.StatusBarItem;
  private readonly statuses = new Map<string, LspStatusSnapshot>();

  constructor() {
    this.output = vscode.window.createOutputChannel('OpenQC LSP');
    this.statusBarItem = vscode.window.createStatusBarItem(vscode.StatusBarAlignment.Left, 95);
    this.statusBarItem.command = 'openqc.lsp.showStatus';
    this.statusBarItem.name = 'OpenQC LSP Status';
    this.setIdleStatus();
  }

  starting(details: LspLifecycleDetails): void {
    this.record('starting', details, 'Starting');
  }

  running(details: LspLifecycleDetails): void {
    this.record('running', details, 'Started');
  }

  stopped(details: LspLifecycleDetails): void {
    this.record('stopped', details, 'Stopped');
  }

  failed(details: LspLifecycleDetails): void {
    this.record('failed', details, 'Failed');
  }

  showLogs(): void {
    this.output.show(true);
  }

  showStatus(): void {
    this.output.appendLine(this.formatStatusReport());
    this.output.show(true);
  }

  appendCompatibilityReport(report: string): void {
    this.output.appendLine(report);
    this.output.show(true);
  }

  getStatusReport(): string {
    return this.formatStatusReport();
  }

  dispose(): void {
    this.statusBarItem.dispose();
    this.output.dispose();
  }

  private record(state: LspLifecycleState, details: LspLifecycleDetails, eventLabel: string): void {
    const snapshot: LspStatusSnapshot = {
      key: details.key,
      state,
      software: details.software,
      languageId: details.languageId,
      workspaceLabel: details.workspaceLabel || 'single file',
      workspaceIdentity: details.workspaceIdentity,
      commandSummary: details.commandSummary,
      hostContext: details.hostContext,
      message: details.message,
      updatedAt: new Date(),
    };

    this.statuses.set(details.key, snapshot);
    this.output.appendLine(this.formatLogLine(eventLabel, snapshot));
    this.updateStatusBar();
  }

  private formatLogLine(eventLabel: string, snapshot: LspStatusSnapshot): string {
    const message = snapshot.message ? ` message="${snapshot.message}"` : '';
    return (
      `[${snapshot.updatedAt.toISOString()}] ${eventLabel} ${snapshot.software} LSP` +
      ` state=${snapshot.state}` +
      ` languageId=${snapshot.languageId}` +
      ` workspace="${snapshot.workspaceLabel}"` +
      ` identity="${snapshot.workspaceIdentity}"` +
      ` command="${snapshot.commandSummary}"` +
      ` context="${snapshot.hostContext}"` +
      message
    );
  }

  private formatStatusReport(): string {
    const snapshots = Array.from(this.statuses.values()).sort((a, b) => a.key.localeCompare(b.key));

    if (snapshots.length === 0) {
      return 'OpenQC LSP status: no language servers have been started in this session.';
    }

    const lines = ['OpenQC LSP status:'];
    for (const snapshot of snapshots) {
      lines.push(
        [
          `- ${snapshot.software} (${snapshot.languageId})`,
          `state=${snapshot.state}`,
          `workspace=${snapshot.workspaceLabel}`,
          `identity=${snapshot.workspaceIdentity}`,
          `command=${snapshot.commandSummary}`,
          `context=${snapshot.hostContext}`,
          `updated=${snapshot.updatedAt.toISOString()}`,
          snapshot.message ? `message=${snapshot.message}` : undefined,
        ]
          .filter(Boolean)
          .join('; ')
      );
    }

    return lines.join('\n');
  }

  private updateStatusBar(): void {
    const snapshots = Array.from(this.statuses.values());
    const active =
      snapshots.find(snapshot => snapshot.state === 'starting') ||
      snapshots.find(snapshot => snapshot.state === 'failed') ||
      snapshots.find(snapshot => snapshot.state === 'running') ||
      snapshots[0];

    if (!active) {
      this.setIdleStatus();
      return;
    }

    const icon = this.getIcon(active.state);
    this.statusBarItem.text = `${icon} OpenQC LSP: ${active.software} ${active.state}`;
    this.statusBarItem.tooltip = this.formatStatusReport();
    this.statusBarItem.show();
  }

  private setIdleStatus(): void {
    this.statusBarItem.text = '$(circle-outline) OpenQC LSP: Idle';
    this.statusBarItem.tooltip = 'OpenQC LSP status: no active language servers.';
    this.statusBarItem.show();
  }

  private getIcon(state: LspLifecycleState): string {
    switch (state) {
      case 'starting':
        return '$(sync~spin)';
      case 'running':
        return '$(check)';
      case 'failed':
        return '$(error)';
      case 'stopped':
        return '$(debug-stop)';
      default:
        return '$(circle-outline)';
    }
  }
}
