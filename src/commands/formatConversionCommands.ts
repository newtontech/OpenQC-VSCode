/**
 * Format Conversion Commands
 *
 * VSCode commands for quantum chemistry file format conversion.
 */

import * as vscode from 'vscode';
import {
  FormatConverter,
  SupportedFormat,
  type ConversionResult,
  type BatchConversionResult,
} from '../converters';

/**
 * Register format conversion commands
 */
export function registerFormatConversionCommands(context: vscode.ExtensionContext): void {
  // Single file conversion
  const convertCommand = vscode.commands.registerCommand(
    'openqc.convertFormat',
    async (uri?: vscode.Uri) => {
      await handleConvertCommand(uri);
    }
  );

  // Quick convert commands for common formats
  const convertToXYZ = vscode.commands.registerCommand(
    'openqc.convertToXYZ',
    async (uri?: vscode.Uri) => {
      await handleQuickConvert(uri, SupportedFormat.XYZ);
    }
  );

  const convertToPDB = vscode.commands.registerCommand(
    'openqc.convertToPDB',
    async (uri?: vscode.Uri) => {
      await handleQuickConvert(uri, SupportedFormat.PDB);
    }
  );

  const convertToVASP = vscode.commands.registerCommand(
    'openqc.convertToVASP',
    async (uri?: vscode.Uri) => {
      await handleQuickConvert(uri, SupportedFormat.POSCAR);
    }
  );

  const convertToGaussian = vscode.commands.registerCommand(
    'openqc.convertToGaussian',
    async (uri?: vscode.Uri) => {
      await handleQuickConvert(uri, SupportedFormat.Gaussian);
    }
  );

  // Batch conversion
  const batchConvertCommand = vscode.commands.registerCommand('openqc.batchConvert', async () => {
    await handleBatchConvert();
  });

  // Check format converter availability
  const checkBackendCommand = vscode.commands.registerCommand(
    'openqc.checkConverterBackend',
    async () => {
      const converter = new FormatConverter();
      const available = await converter.checkBackend();

      if (available) {
        vscode.window.showInformationMessage('Format converter backend is available and ready.');
      } else {
        const action = await vscode.window.showWarningMessage(
          'Format converter backend is not available. Python and dpdata are required.',
          'Install Instructions',
          'Close'
        );

        if (action === 'Install Instructions') {
          await showInstallInstructions();
        }
      }
    }
  );

  context.subscriptions.push(
    convertCommand,
    convertToXYZ,
    convertToPDB,
    convertToVASP,
    convertToGaussian,
    batchConvertCommand,
    checkBackendCommand
  );
}

/**
 * Handle the convert format command
 */
async function handleConvertCommand(uri?: vscode.Uri): Promise<void> {
  const fileUri = uri || vscode.window.activeTextEditor?.document.uri;
  if (!fileUri) {
    vscode.window.showErrorMessage('No file selected for conversion');
    return;
  }

  // Show format picker
  const targetFormat = await vscode.window.showQuickPick(
    [
      { label: 'XYZ', description: 'Standard molecular coordinates', format: SupportedFormat.XYZ },
      { label: 'PDB', description: 'Protein Data Bank format', format: SupportedFormat.PDB },
      {
        label: 'CIF',
        description: 'Crystallographic Information File',
        format: SupportedFormat.CIF,
      },
      { label: 'VASP', description: 'VASP POSCAR format', format: SupportedFormat.POSCAR },
      { label: 'Gaussian', description: 'Gaussian input format', format: SupportedFormat.Gaussian },
      { label: 'ORCA', description: 'ORCA input format', format: SupportedFormat.ORCA },
    ],
    {
      placeHolder: 'Select target format',
    }
  );

  if (!targetFormat) {
    return;
  }

  await convertFile(fileUri.fsPath, targetFormat.format);
}

/**
 * Handle quick convert to a specific format
 */
async function handleQuickConvert(
  uri: vscode.Uri | undefined,
  format: SupportedFormat
): Promise<void> {
  const fileUri = uri || vscode.window.activeTextEditor?.document.uri;
  if (!fileUri) {
    vscode.window.showErrorMessage('No file selected for conversion');
    return;
  }

  await convertFile(fileUri.fsPath, format);
}

/**
 * Convert a file to the target format
 */
async function convertFile(filePath: string, targetFormat: SupportedFormat): Promise<void> {
  const converter = new FormatConverter();

  // Check backend availability
  if (!(await converter.checkBackend())) {
    const action = await vscode.window.showWarningMessage(
      'Format converter backend is not available. Install Python and dpdata?',
      'Install',
      'Cancel'
    );

    if (action === 'Install') {
      await showInstallInstructions();
    }
    return;
  }

  // Determine output path
  const inputBasename = filePath.replace(/\.[^/.]+$/, '');
  const extension = converter.getExtensionForFormat(targetFormat);
  const outputPath = `${inputBasename}.${extension}`;

  // Show progress
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Converting to ${targetFormat.toUpperCase()}...`,
      cancellable: false,
    },
    async () => {
      const result = await converter.convert(filePath, outputPath, undefined, targetFormat);

      if (result.success) {
        const action = await vscode.window.showInformationMessage(
          `Successfully converted to ${targetFormat.toUpperCase()}`,
          'Open File',
          'Copy to Clipboard',
          'Show in Folder'
        );

        if (action === 'Open File') {
          const doc = await vscode.workspace.openTextDocument(outputPath);
          await vscode.window.showTextDocument(doc);
        } else if (action === 'Copy to Clipboard') {
          const content = await vscode.workspace.fs.readFile(vscode.Uri.file(outputPath));
          await vscode.env.clipboard.writeText(Buffer.from(content).toString('utf8'));
          vscode.window.showInformationMessage('Content copied to clipboard');
        } else if (action === 'Show in Folder') {
          vscode.env.openExternal(vscode.Uri.file(outputPath));
        }
      } else {
        vscode.window.showErrorMessage(`Conversion failed: ${result.error}`);
      }
    }
  );
}

/**
 * Handle batch conversion command
 */
async function handleBatchConvert(): Promise<void> {
  // Select files
  const files = await vscode.window.showOpenDialog({
    canSelectFiles: true,
    canSelectFolders: false,
    canSelectMany: true,
    title: 'Select files to convert',
    filters: {
      'All Quantum Chemistry Files': [
        '*',
        'POSCAR',
        'CONTCAR',
        'INCAR',
        '*.xyz',
        '*.pdb',
        '*.cif',
        '*.gjf',
        '*.com',
        '*.inp',
      ],
    },
  });

  if (!files || files.length === 0) {
    return;
  }

  // Select target format
  const targetFormat = await vscode.window.showQuickPick(
    [
      { label: 'XYZ', format: SupportedFormat.XYZ },
      { label: 'PDB', format: SupportedFormat.PDB },
      { label: 'CIF', format: SupportedFormat.CIF },
      { label: 'VASP', format: SupportedFormat.POSCAR },
    ],
    { placeHolder: 'Select target format' }
  );

  if (!targetFormat) {
    return;
  }

  // Select output directory
  const outputDirUri = await vscode.window.showOpenDialog({
    canSelectFiles: false,
    canSelectFolders: true,
    canSelectMany: false,
    title: 'Select output directory',
  });

  if (!outputDirUri || outputDirUri.length === 0) {
    return;
  }

  const outputDir = outputDirUri[0].fsPath;
  const filePaths = files.map(f => f.fsPath);

  const converter = new FormatConverter();

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Converting ${filePaths.length} files...`,
      cancellable: false,
    },
    async () => {
      const result = await converter.batchConvert(filePaths, outputDir, targetFormat.format);

      if (result.success) {
        vscode.window
          .showInformationMessage(
            `Successfully converted ${result.successful} of ${result.total} files`,
            'Open Folder'
          )
          .then(action => {
            if (action === 'Open Folder') {
              vscode.env.openExternal(vscode.Uri.file(outputDir));
            }
          });
      } else {
        vscode.window.showErrorMessage(
          `Batch conversion completed with errors: ${result.failed} of ${result.total} failed`
        );
      }
    }
  );
}

/**
 * Show installation instructions for the Python backend
 */
async function showInstallInstructions(): Promise<void> {
  const panel = vscode.window.createWebviewPanel(
    'openqc.converter.setup',
    'OpenQC Format Converter Setup',
    vscode.ViewColumn.One,
    { enableScripts: true }
  );

  panel.webview.html = `
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8">
      <meta name="viewport" content="width=device-width, initial-scale=1.0">
      <title>Format Converter Setup</title>
      <style>
        body {
          font-family: var(--vscode-font-family);
          padding: 20px;
          color: var(--vscode-foreground);
          background-color: var(--vscode-editor-background);
        }
        h1 { color: var(--vscode-foreground); }
        h2 { color: var(--vscode-foreground); margin-top: 24px; }
        code {
          background-color: var(--vscode-textCodeBlock-background);
          padding: 2px 6px;
          border-radius: 3px;
          font-family: var(--vscode-editor-font-family);
        }
        pre {
          background-color: var(--vscode-textCodeBlock-background);
          padding: 12px;
          border-radius: 6px;
          overflow-x: auto;
        }
        ul { line-height: 1.6; }
      </style>
    </head>
    <body>
      <h1>OpenQC Format Converter Setup</h1>
      <p>The format converter requires Python and the dpdata library.</p>

      <h2>Prerequisites</h2>
      <ul>
        <li>Python 3.7 or higher</li>
        <li>pip (Python package installer)</li>
      </ul>

      <h2>Installation Steps</h2>
      <ol>
        <li>Install Python if not already installed:
          <pre><code># On Ubuntu/Debian
sudo apt-get install python3 python3-pip

# On macOS (with Homebrew)
brew install python3

# On Windows
# Download from https://www.python.org/downloads/</code></pre>
        </li>
        <li>Install dpdata:
          <pre><code>pip install dpdata</code></pre>
        </li>
      </ol>

      <h2>Verify Installation</h2>
      <p>Run the following command to verify:</p>
      <pre><code>python3 -c "import dpdata; print(dpdata.__version__)"</code></pre>

      <h2>Supported Formats</h2>
      <ul>
        <li><strong>Input:</strong> VASP, Gaussian, ORCA, XYZ, PDB, CIF, CP2K, Quantum ESPRESSO</li>
        <li><strong>Output:</strong> VASP, Gaussian, ORCA, XYZ, PDB, CIF</li>
      </ul>

      <p>After installation, click the "Check Backend" command in the command palette to verify.</p>
    </body>
    </html>
  `;
}
