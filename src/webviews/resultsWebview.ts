import { generateNonce } from '../utils/nonce';

export interface ResultsWebviewData {
  jobId: string;
  jobName: string;
  software: string;
  status: string;
  startTime?: string;
  endTime?: string;
  duration?: string;
  output: string;
}

interface ResultsViewModel {
  csp: string;
  nonce: string;
  jobName: string;
  jobId: string;
  software: string;
  status: string;
  statusColor: string;
  duration: string;
  startTime: string;
  endTime: string;
  output: string;
}

export function renderResultsWebviewHtml(data: ResultsWebviewData, cspSource: string): string {
  const viewModel = createResultsViewModel(data, cspSource);

  return `<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta http-equiv="Content-Security-Policy" content="${viewModel.csp}">
    <title>Results: ${viewModel.jobName}</title>
    <style nonce="${viewModel.nonce}">
        body {
            margin: 0;
            padding: 0;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1e1e1e;
            color: #cccccc;
        }
        #header {
            padding: 15px 20px;
            background: #252526;
            border-bottom: 1px solid #3c3c3c;
        }
        #title {
            font-size: 16px;
            font-weight: 600;
            margin-bottom: 5px;
        }
        #status {
            display: inline-block;
            padding: 2px 8px;
            border-radius: 3px;
            font-size: 11px;
            text-transform: uppercase;
            background: ${viewModel.statusColor};
            color: white;
        }
        #content {
            padding: 20px;
        }
        .info-section {
            background: #252526;
            border: 1px solid #3c3c3c;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 15px;
        }
        .info-section h3 {
            margin-top: 0;
            margin-bottom: 10px;
            font-size: 14px;
            color: #0e639c;
        }
        .info-row {
            display: flex;
            padding: 5px 0;
            border-bottom: 1px solid #3c3c3c;
        }
        .info-row:last-child {
            border-bottom: none;
        }
        .info-label {
            width: 150px;
            color: #9cdcfe;
            font-weight: 500;
        }
        .info-value {
            flex: 1;
            color: #ce9178;
        }
        #output {
            background: #1e1e1e;
            border: 1px solid #3c3c3c;
            border-radius: 5px;
            padding: 15px;
            font-family: 'Consolas', monospace;
            font-size: 12px;
            white-space: pre-wrap;
            color: #cccccc;
        }
    </style>
</head>
<body>
    <div id="header">
        <div id="title">${viewModel.jobName}</div>
        <div><span id="status">${viewModel.status}</span> ${viewModel.software}</div>
    </div>
    <div id="content">
        <div class="info-section">
            <h3>Job Information</h3>
            <div class="info-row">
                <span class="info-label">Job ID:</span>
                <span class="info-value">${viewModel.jobId}</span>
            </div>
            <div class="info-row">
                <span class="info-label">Software:</span>
                <span class="info-value">${viewModel.software}</span>
            </div>
            <div class="info-row">
                <span class="info-label">Status:</span>
                <span class="info-value">${viewModel.status}</span>
            </div>
            <div class="info-row">
                <span class="info-label">Duration:</span>
                <span class="info-value">${viewModel.duration}</span>
            </div>
            ${renderOptionalInfoRow('Started:', viewModel.startTime)}
            ${renderOptionalInfoRow('Completed:', viewModel.endTime)}
        </div>
        <div class="info-section">
            <h3>Output</h3>
            <div id="output">${viewModel.output}</div>
        </div>
    </div>
</body>
</html>`;
}

function createResultsViewModel(data: ResultsWebviewData, cspSource: string): ResultsViewModel {
  const nonce = getNonce();

  return {
    csp: getResultsWebviewCsp(cspSource, nonce),
    nonce,
    jobName: escapeHtml(data.jobName),
    jobId: escapeHtml(data.jobId),
    software: escapeHtml(data.software),
    status: escapeHtml(data.status),
    statusColor: data.status === 'completed' ? '#3c9' : '#f66',
    duration: escapeHtml(data.duration || 'N/A'),
    startTime: formatOptionalDate(data.startTime),
    endTime: formatOptionalDate(data.endTime),
    output: escapeHtml(data.output),
  };
}

function renderOptionalInfoRow(label: string, value: string): string {
  if (!value) {
    return '';
  }

  return `
            <div class="info-row">
                <span class="info-label">${label}</span>
                <span class="info-value">${value}</span>
            </div>
            `;
}

function formatOptionalDate(value: string | undefined): string {
  return value ? escapeHtml(new Date(value).toLocaleString()) : '';
}

function getResultsWebviewCsp(cspSource: string, nonce: string): string {
  return [
    `default-src 'none';`,
    `script-src 'none';`,
    `style-src 'nonce-${nonce}';`,
    `img-src ${cspSource} data:;`,
    `font-src ${cspSource};`,
  ].join(' ');
}

function escapeHtml(value: unknown): string {
  return String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

const getNonce = generateNonce;
