/**
 * Internationalization (i18n) utility for OpenQC extension
 *
 * Provides a simple key-based localization system.
 * Falls back to the default (English) string when a translation is missing.
 *
 * @module i18n
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/15
 */

/**
 * Supported locale codes.
 * Add new locales here as translations are contributed.
 */
export type SupportedLocale = 'en' | 'zh-CN';

/**
 * Flat key-value map for a single locale.
 * Keys use dot-notation for namespacing: "commands.visualize", "errors.noAtoms".
 */
export type LocaleMessages = Record<string, string>;

/**
 * All locale bundles keyed by locale code.
 */
export type LocaleBundle = Record<SupportedLocale, LocaleMessages>;

const DEFAULT_LOCALE: SupportedLocale = 'en';

/**
 * English message bundle (canonical source of truth).
 */
const en: LocaleMessages = {
  // Extension activation
  'extension.activated': 'OpenQC-VSCode extension activated',

  // Commands
  'command.visualizeStructure': 'OpenQC: Visualize Molecular Structure',
  'command.plotData': 'OpenQC: Plot Calculation Data',
  'command.previewInput': 'OpenQC: Preview Input File',
  'command.startLSP': 'OpenQC: Start Language Server',
  'command.stopLSP': 'OpenQC: Stop Language Server',
  'command.restartLSP': 'OpenQC: Restart Language Server',
  'command.validate': 'OpenQC: Validate Input File',
  'command.convertFormat': 'OpenQC: Convert File Format',
  'command.convertToXYZ': 'OpenQC: Convert to XYZ',
  'command.convertToPDB': 'OpenQC: Convert to PDB',
  'command.convertToVASP': 'OpenQC: Convert to VASP',
  'command.convertToGaussian': 'OpenQC: Convert to Gaussian',
  'command.batchConvert': 'OpenQC: Batch Convert Files',

  // Sidebar
  'sidebar.molecules': 'Molecules',
  'sidebar.jobs': 'Calculation Jobs',
  'sidebar.refreshMolecules': 'Refresh Molecules',
  'sidebar.refreshJobs': 'Refresh Jobs',
  'sidebar.runCalculation': 'Run Calculation',
  'sidebar.openMolecule': 'Open Molecule',
  'sidebar.deleteMolecule': 'Delete Molecule',
  'sidebar.viewResults': 'View Results',
  'sidebar.exportData': 'Export Data',
  'sidebar.cancelJob': 'Cancel Job',
  'sidebar.restartJob': 'Restart Job',

  // Messages
  'message.noActiveEditor': 'No active text editor',
  'message.unsupportedFileType': 'Unsupported file type for visualization',
  'message.noStructureFound': 'No molecular structure found in file',
  'message.visualizeFailed': 'Failed to visualize structure: {0}',
  'message.validationComplete': 'Input file validated',
  'message.moleculesRefreshed': 'Molecules refreshed',
  'message.jobsRefreshed': 'Jobs refreshed',
  'message.calculationStarted': 'Started calculation: {0}',
  'message.calculationNamePrompt': 'Enter calculation name',
  'message.calculationNamePlaceholder': 'e.g., Geometry Optimization',
  'message.removedMolecule': 'Removed molecule: {0}',
  'message.selectedMolecule': 'Selected molecule: {0} ({1})',
  'message.resultsNotAvailable': 'Results not available for {0} jobs',
  'message.cannotExportRunning': 'Can only export data from completed jobs',
  'message.exported': 'Exported {0} to {1}',
  'message.exportFailed': 'Failed to export data: {0}',
  'message.cancelledJob': 'Cancelled job: {0}',
  'message.restartedJob': 'Restarted job: {0}',
  'message.openFileFailed': 'Failed to open file: {0}',

  // LSP messages
  'lsp.started': '{0} Language Server started',
  'lsp.startFailed': 'Failed to start {0} Language Server: {1}',
  'lsp.stopFailed': 'Error stopping {0} Language Server: {1}',
  'lsp.restartFailed': 'Failed to restart {0} Language Server: {1}',
  'lsp.notDetected': 'Could not detect quantum chemistry software for this file',

  // Errors
  'error.unknown': 'An unknown error occurred',
  'error.loadJobs': 'Failed to load jobs',
  'error.saveJobs': 'Failed to save jobs',

  // Configuration
  'config.logging.level': 'Set the logging level for OpenQC extension',
  'config.logging.showUserMessages': 'Show error/warning messages to user via VS Code UI',
};

/**
 * Simplified Chinese message bundle.
 * Add translations below as they are contributed.
 */
const zhCN: LocaleMessages = {
  // Extension activation
  'extension.activated': 'OpenQC-VSCode 扩展已激活',

  // Commands
  'command.visualizeStructure': 'OpenQC: 可视化分子结构',
  'command.plotData': 'OpenQC: 绘制计算数据',
  'command.previewInput': 'OpenQC: 预览输入文件',
  'command.startLSP': 'OpenQC: 启动语言服务器',
  'command.stopLSP': 'OpenQC: 停止语言服务器',
  'command.restartLSP': 'OpenQC: 重启语言服务器',
  'command.validate': 'OpenQC: 验证输入文件',

  // Messages
  'message.noActiveEditor': '没有活动的文本编辑器',
  'message.unsupportedFileType': '不支持的可视化文件类型',
  'message.noStructureFound': '文件中未找到分子结构',
  'message.validationComplete': '输入文件验证完成',
  'message.moleculesRefreshed': '分子列表已刷新',
  'message.jobsRefreshed': '计算任务已刷新',
  'message.calculationStarted': '已开始计算: {0}',
  'message.cancelledJob': '已取消任务: {0}',
  'message.restartedJob': '已重启任务: {0}',

  // LSP messages
  'lsp.started': '{0} 语言服务器已启动',
  'lsp.startFailed': '启动 {0} 语言服务器失败: {1}',
  'lsp.notDetected': '无法检测此文件的量子化学软件类型',

  // Errors
  'error.loadJobs': '加载任务失败',
  'error.saveJobs': '保存任务失败',

  // Configuration
  'config.logging.level': '设置 OpenQC 扩展的日志级别',
  'config.logging.showUserMessages': '通过 VS Code UI 显示错误/警告消息',
};

const bundles: LocaleBundle = { en, 'zh-CN': zhCN };

let currentLocale: SupportedLocale = DEFAULT_LOCALE;

/**
 * Set the active locale for message lookups.
 *
 * @param locale - Locale code to activate
 */
export function setLocale(locale: SupportedLocale): void {
  currentLocale = locale;
}

/**
 * Get the currently active locale.
 *
 * @returns Current locale code
 */
export function getLocale(): SupportedLocale {
  return currentLocale;
}

/**
 * Look up a localized message by key.
 *
 * Substitutes `{0}`, `{1}`, etc. with the provided positional args.
 * Falls back to the English bundle when the key is missing from the current locale.
 * Returns the raw key when not found in any bundle.
 *
 * @param key - Dot-notation message key (e.g. "lsp.started")
 * @param args - Positional substitution values
 * @returns Localized (or fallback) message string
 */
export function t(key: string, ...args: (string | number)[]): string {
  const bundle = bundles[currentLocale] || bundles[DEFAULT_LOCALE];
  let message = bundle[key] || bundles[DEFAULT_LOCALE][key] || key;

  args.forEach((arg, index) => {
    message = message.replace(`{${index}}`, String(arg));
  });

  return message;
}
