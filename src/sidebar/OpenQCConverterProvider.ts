import * as vscode from 'vscode';
import * as path from 'path';

/**
 * OpenQC Sidebar Provider
 *
 * 提供VSCode侧边栏界面，包含：
 * - 格式转换按钮
 * - 快速转换预设
 * - 格式选择器
 * - 最近转换历史
 */
export class OpenQCConverterProvider implements vscode.WebviewViewProvider {
  public static readonly viewType = 'openqc.converterSidebar';

  private _view?: vscode.WebviewView;

  constructor(private readonly _extensionUri: vscode.Uri) {}

  public resolveWebviewView(
    webviewView: vscode.WebviewView,
    context: vscode.WebviewViewResolveContext,
    _token: vscode.CancellationToken
  ) {
    this._view = webviewView;

    webviewView.webview.options = {
      enableScripts: true,
      localResourceRoots: [this._extensionUri],
    };

    webviewView.webview.html = this._getHtmlForWebview(webviewView.webview);

    // 处理来自Webview的消息
    webviewView.webview.onDidReceiveMessage(async data => {
      switch (data.type) {
        case 'convertToASE':
          vscode.commands.executeCommand('openqc.convertToASE');
          break;
        case 'convertFromASE':
          vscode.commands.executeCommand('openqc.convertFromASE');
          break;
        case 'migrateFormat':
          vscode.commands.executeCommand('openqc.migrateFormat');
          break;
        case 'quickConvert':
          vscode.commands.executeCommand('openqc.quickConvert', data.from, data.to);
          break;
        case 'openSettings':
          vscode.commands.executeCommand('workbench.action.openSettings', 'openqc');
          break;
      }
    });
  }

  private _getHtmlForWebview(webview: vscode.Webview): string {
    const nonce = getNonce();

    return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src ${webview.cspSource} 'unsafe-inline'; script-src 'nonce-${nonce}';">
    <style>
        body {
            font-family: var(--vscode-font-family);
            font-size: var(--vscode-font-size);
            color: var(--vscode-foreground);
            padding: 10px;
        }
        
        .section {
            margin-bottom: 20px;
        }
        
        .section-title {
            font-weight: bold;
            margin-bottom: 10px;
            color: var(--vscode-descriptionForeground);
            text-transform: uppercase;
            font-size: 11px;
            letter-spacing: 0.5px;
        }
        
        .button {
            display: flex;
            align-items: center;
            padding: 8px 12px;
            margin-bottom: 6px;
            background: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 13px;
            width: 100%;
            transition: background 0.2s;
        }
        
        .button:hover {
            background: var(--vscode-button-hoverBackground);
        }
        
        .button .icon {
            margin-right: 8px;
            font-size: 14px;
        }
        
        .quick-convert {
            display: grid;
            grid-template-columns: 1fr auto 1fr auto;
            gap: 6px;
            align-items: center;
            margin-bottom: 6px;
        }
        
        select {
            background: var(--vscode-dropdown-background);
            color: var(--vscode-dropdown-foreground);
            border: 1px solid var(--vscode-dropdown-border);
            padding: 4px 8px;
            border-radius: 3px;
            font-size: 12px;
        }
        
        .convert-btn {
            background: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            padding: 4px 10px;
            border-radius: 3px;
            cursor: pointer;
            font-size: 11px;
        }
        
        .supported-formats {
            display: flex;
            flex-wrap: wrap;
            gap: 4px;
            margin-top: 8px;
        }
        
        .format-tag {
            background: var(--vscode-badge-background);
            color: var(--vscode-badge-foreground);
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 10px;
        }
        
        .info-box {
            background: var(--vscode-textBlockQuote-background);
            border-left: 3px solid var(--vscode-textBlockQuote-border);
            padding: 8px 12px;
            margin-top: 15px;
            font-size: 11px;
            line-height: 1.4;
        }
        
        .status-indicator {
            display: inline-block;
            width: 8px;
            height: 8px;
            border-radius: 50%;
            margin-right: 6px;
        }
        
        .status-ready { background: #4caf50; }
        .status-busy { background: #ff9800; }
    </style>
</head>
<body>
    <div class="section">
        <div class="section-title">
            <span class="status-indicator status-ready"></span>
            ASE Converter
        </div>
        
        <button class="button" onclick="convertToASE()">
            <span class="icon">→</span>
            Convert to ASE Atoms
        </button>
        
        <button class="button" onclick="convertFromASE()">
            <span class="icon">←</span>
            Convert from ASE Atoms
        </button>
        
        <button class="button" onclick="migrateFormat()">
            <span class="icon">⇄</span>
            Migrate Format
        </button>
    </div>
    
    <div class="section">
        <div class="section-title">Quick Convert</div>
        
        <div class="quick-convert">
            <select id="fromFormat">
                <option value="vasp">VASP</option>
                <option value="cp2k">CP2K</option>
                <option value="espresso">QE</option>
                <option value="gaussian">Gaussian</option>
                <option value="orca">ORCA</option>
                <option value="xyz">XYZ</option>
            </select>
            
            <span>→</span>
            
            <select id="toFormat">
                <option value="cp2k">CP2K</option>
                <option value="vasp">VASP</option>
                <option value="espresso">QE</option>
                <option value="gaussian">Gaussian</option>
                <option value="orca">ORCA</option>
                <option value="xyz">XYZ</option>
            </select>
            
            <button class="convert-btn" onclick="quickConvert()">Go</button>
        </div>
        
        <div class="supported-formats">
            <span class="format-tag">VASP</span>
            <span class="format-tag">CP2K</span>
            <span class="format-tag">QE</span>
            <span class="format-tag">Gaussian</span>
            <span class="format-tag">ORCA</span>
            <span class="format-tag">NWChem</span>
            <span class="format-tag">GAMESS</span>
            <span class="format-tag">LAMMPS</span>
            <span class="format-tag">XYZ</span>
            <span class="format-tag">PDB</span>
            <span class="format-tag">CIF</span>
        </div>
    </div>
    
    <div class="info-box">
        <strong>OpenQC ASE Integration</strong><br>
        Universal converter supporting 11 formats.<br>
        Powered by ASE (Atomic Simulation Environment).
    </div>
    
    <script nonce="${nonce}">
        const vscode = acquireVsCodeApi();
        
        function convertToASE() {
            vscode.postMessage({ type: 'convertToASE' });
        }
        
        function convertFromASE() {
            vscode.postMessage({ type: 'convertFromASE' });
        }
        
        function migrateFormat() {
            vscode.postMessage({ type: 'migrateFormat' });
        }
        
        function quickConvert() {
            const from = document.getElementById('fromFormat').value;
            const to = document.getElementById('toFormat').value;
            vscode.postMessage({ type: 'quickConvert', from, to });
        }
    </script>
</body>
</html>`;
  }
}

function getNonce(): string {
  let text = '';
  const possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
  for (let i = 0; i < 32; i++) {
    text += possible.charAt(Math.floor(Math.random() * possible.length));
  }
  return text;
}
