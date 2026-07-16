import * as vscode from 'vscode';
import {
  CalculatorConfig,
  CalculatorFactory,
  CalculatorType,
  getCalculatorConfiguration,
} from '../ase/ASECalculator';
import { JobItem, JobResultData, JobTreeProvider } from '../sidebar/JobTreeProvider';
import { renderResultsWebviewHtml } from '../webviews/resultsWebview';

interface CalculatorPickItem {
  label: string;
  description: CalculatorType;
}

const CALCULATOR_ITEMS: CalculatorPickItem[] = [
  { label: 'VASP', description: CalculatorType.VASP },
  { label: 'CP2K', description: CalculatorType.CP2K },
  { label: 'Quantum ESPRESSO', description: CalculatorType.QE },
];

export function registerJobCommands(
  context: vscode.ExtensionContext,
  jobProvider: JobTreeProvider
): void {
  const runningJobs = new Map<string, AbortController>();

  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.sidebar.refreshJobs', () => {
      jobProvider.refresh();
      vscode.window.showInformationMessage('Jobs refreshed');
    }),

    vscode.commands.registerCommand('openqc.sidebar.runCalculation', async () => {
      await runCalculationFromSidebar(context, jobProvider, runningJobs);
    }),

    vscode.commands.registerCommand('openqc.sidebar.viewResults', async (item: JobItem) => {
      await viewJobResults(item);
    }),

    vscode.commands.registerCommand('openqc.sidebar.exportData', async (item: JobItem) => {
      await exportJobData(item);
    }),

    vscode.commands.registerCommand('openqc.sidebar.cancelJob', (item: JobItem) => {
      const controller = runningJobs.get(item.id);
      controller?.abort();
      const cancelled = jobProvider.cancelJob(item.id);
      if (cancelled) {
        vscode.window.showInformationMessage(`Cancelled job: ${item.label}`);
      } else {
        vscode.window.showInformationMessage(`Job is not running: ${item.label}`);
      }
    }),

    vscode.commands.registerCommand('openqc.sidebar.restartJob', async (item: JobItem) => {
      await restartJobFromSidebar(context, jobProvider, item, runningJobs);
    })
  );
}

async function runCalculationFromSidebar(
  context: vscode.ExtensionContext,
  jobProvider: JobTreeProvider,
  runningJobs: Map<string, AbortController>
): Promise<void> {
  const selectedCalculator = await vscode.window.showQuickPick(CALCULATOR_ITEMS, {
    placeHolder: 'Select calculator for an existing input directory',
  });
  if (!selectedCalculator) {
    return;
  }

  const directoryUris = await vscode.window.showOpenDialog({
    canSelectFiles: false,
    canSelectFolders: true,
    canSelectMany: false,
    openLabel: 'Select Input Directory',
  });
  const inputDir = directoryUris?.[0]?.fsPath;
  if (!inputDir) {
    return;
  }

  const defaultName = `${selectedCalculator.label} calculation`;
  const nameInput = await vscode.window.showInputBox({
    prompt: 'Enter calculation name',
    placeHolder: defaultName,
    value: defaultName,
  });
  if (nameInput === undefined) {
    return;
  }
  const name = nameInput.trim() || defaultName;

  const job = new JobItem(
    `job-${Date.now()}`,
    name,
    'running',
    5,
    selectedCalculator.label,
    new Date(),
    undefined,
    inputDir
  );
  jobProvider.addJob(job);
  vscode.window.showInformationMessage(`Started calculation: ${name}`);

  await executeJob(context, jobProvider, job, selectedCalculator.description, runningJobs);
}

async function restartJobFromSidebar(
  context: vscode.ExtensionContext,
  jobProvider: JobTreeProvider,
  item: JobItem,
  runningJobs: Map<string, AbortController>
): Promise<void> {
  if (!item.workingDirectory) {
    vscode.window.showWarningMessage(`Cannot restart ${item.label}: no working directory recorded`);
    return;
  }

  const calculatorType = getCalculatorTypeForSoftware(item.software);
  if (!calculatorType) {
    vscode.window.showWarningMessage(
      `Cannot restart ${item.label}: unsupported calculator ${item.software}`
    );
    return;
  }

  const job = jobProvider.restartJob(item.id);
  if (!job) {
    vscode.window.showWarningMessage(`Cannot restart ${item.label}: job was not found`);
    return;
  }

  vscode.window.showInformationMessage(`Restarted job: ${item.label}`);
  await executeJob(context, jobProvider, job, calculatorType, runningJobs);
}

async function executeJob(
  context: vscode.ExtensionContext,
  jobProvider: JobTreeProvider,
  job: JobItem,
  calculatorType: CalculatorType,
  runningJobs: Map<string, AbortController>
): Promise<void> {
  const inputDir = job.workingDirectory;
  if (!inputDir) {
    jobProvider.updateJobResult(
      job.id,
      'failed',
      {
        success: false,
        error: 'No working directory recorded for calculation job',
        warnings: [],
        metadata: {},
        outputFiles: [],
      },
      100
    );
    return;
  }

  const calculator = new CalculatorFactory(context).createCalculator(
    createCalculatorConfig(calculatorType)
  );

  const controller = new AbortController();
  runningJobs.set(job.id, controller);
  jobProvider.updateJobStatus(job.id, 'running', 25);
  try {
    const result = await calculator.runCalculation(inputDir, undefined, controller.signal);
    const latest = jobProvider.getJob(job.id);
    if (controller.signal.aborted || latest?.status === 'cancelled') {
      return;
    }
    const resultData = normalizeJobResult(result);

    if (result.success) {
      jobProvider.updateJobResult(job.id, 'completed', resultData, 100);
      vscode.window.showInformationMessage(`Completed calculation: ${job.label}`);
    } else {
      jobProvider.updateJobResult(job.id, 'failed', resultData, 100);
      vscode.window.showErrorMessage(`Calculation failed: ${result.error ?? 'unknown error'}`);
    }
  } finally {
    runningJobs.delete(job.id);
  }
}

function getCalculatorTypeForSoftware(software: string): CalculatorType | undefined {
  return CALCULATOR_ITEMS.find(
    item => item.label === software || item.description === software.toLowerCase()
  )?.description;
}

function createCalculatorConfig(type: CalculatorType): CalculatorConfig {
  const calculatorConfig = getCalculatorConfiguration();
  const settings = calculatorConfig[type] ?? {};
  return {
    type,
    command: settings.command,
    parameters: settings.defaultParams ?? {},
  };
}

function normalizeJobResult(result: {
  success: boolean;
  energy?: number;
  forces?: number[][];
  stress?: number[][];
  error?: string;
  warnings?: string[];
  metadata?: Record<string, any>;
  outputFiles?: string[];
}): JobResultData {
  return {
    success: result.success,
    energy: result.energy,
    forces: result.forces,
    stress: result.stress,
    error: result.error,
    warnings: result.warnings ?? [],
    metadata: result.metadata ?? {},
    outputFiles: result.outputFiles ?? [],
  };
}

async function viewJobResults(item: JobItem): Promise<void> {
  if (item.status !== 'completed' && item.status !== 'failed') {
    vscode.window.showInformationMessage(`Results not available for ${item.status} jobs`);
    return;
  }
  if (!item.result) {
    vscode.window.showWarningMessage(`No captured result data is available for ${item.label}`);
    return;
  }

  const panel = vscode.window.createWebviewPanel(
    'openqc.results',
    `Results: ${item.label}`,
    vscode.ViewColumn.Two,
    {
      retainContextWhenHidden: true,
      localResourceRoots: [],
    }
  );

  panel.webview.html = renderResultsWebviewHtml(
    {
      jobId: item.id,
      jobName: item.label,
      software: item.software,
      status: item.status,
      startTime: item.startTime?.toISOString(),
      endTime: item.endTime?.toISOString(),
      duration: typeof item.tooltip === 'string' ? item.tooltip.split('Duration: ')[1] : undefined,
      output: formatJobResultOutput(item),
    },
    panel.webview.cspSource
  );
}

async function exportJobData(item: JobItem): Promise<void> {
  if (item.status !== 'completed') {
    vscode.window.showWarningMessage('Can only export data from completed jobs');
    return;
  }
  if (!item.result) {
    vscode.window.showWarningMessage(`No captured result data is available for ${item.label}`);
    return;
  }

  const uri = await vscode.window.showSaveDialog({
    filters: {
      JSON: ['json'],
      CSV: ['csv'],
      'All Files': ['*'],
    },
    defaultUri: vscode.Uri.file(`${item.label.replace(/\s+/g, '_')}_results.json`),
    saveLabel: 'Export Results',
  });
  if (!uri) {
    return;
  }

  const payload = createJobExportPayload(item);
  const content = uri.fsPath.toLowerCase().endsWith('.csv')
    ? renderJobResultCsv(payload)
    : JSON.stringify(payload, null, 2);
  await vscode.workspace.fs.writeFile(vscode.Uri.file(uri.fsPath), Buffer.from(content));
  vscode.window.showInformationMessage(`Exported ${item.label} to ${uri.fsPath}`);
}

export function createJobExportPayload(item: JobItem): Record<string, any> {
  return {
    jobId: item.id,
    jobName: item.label,
    software: item.software,
    status: item.status,
    workingDirectory: item.workingDirectory,
    startTime: item.startTime?.toISOString(),
    endTime: item.endTime?.toISOString(),
    timestamp: new Date().toISOString(),
    result: item.result,
  };
}

export function formatJobResultOutput(item: JobItem): string {
  const result = item.result;
  if (!result) {
    return 'No captured result data is available.';
  }

  const lines = [
    `Results for ${item.label}`,
    '',
    `Software: ${item.software}`,
    `Status: ${item.status}`,
  ];
  if (item.workingDirectory) {
    lines.push(`Directory: ${item.workingDirectory}`);
  }
  if (typeof result.energy === 'number') {
    lines.push(`Final energy: ${result.energy}`);
  }
  if (result.forces?.length) {
    lines.push(`Forces: ${result.forces.length} atoms`);
  }
  if (result.stress?.length) {
    lines.push(`Stress rows: ${result.stress.length}`);
  }
  if (result.outputFiles?.length) {
    lines.push(`Output files: ${result.outputFiles.join(', ')}`);
  }
  if (result.warnings?.length) {
    lines.push('', 'Warnings:', ...result.warnings.map(warning => `- ${warning}`));
  }
  if (result.error) {
    lines.push('', 'Error:', result.error);
  }
  if (result.metadata && Object.keys(result.metadata).length > 0) {
    lines.push('', 'Metadata:', JSON.stringify(result.metadata, null, 2));
  }

  return lines.join('\n');
}

function renderJobResultCsv(payload: Record<string, any>): string {
  const result = payload.result ?? {};
  return [
    'jobId,jobName,software,status,workingDirectory,energy,error,outputFiles',
    [
      payload.jobId,
      payload.jobName,
      payload.software,
      payload.status,
      payload.workingDirectory ?? '',
      result.energy ?? '',
      result.error ?? '',
      Array.isArray(result.outputFiles) ? result.outputFiles.join(';') : '',
    ]
      .map(csvCell)
      .join(','),
  ].join('\n');
}

function csvCell(value: unknown): string {
  const text = String(value);
  if (!/[",\n]/.test(text)) {
    return text;
  }
  return `"${text.replace(/"/g, '""')}"`;
}
