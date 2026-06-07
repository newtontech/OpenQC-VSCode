/**
 * OpenQC Logger Utility
 *
 * Provides centralized logging with configurable log levels and output to VS Code Output panel.
 * Replaces inconsistent console.log/error patterns throughout the codebase.
 *
 * @module Logger
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/39
 */

import * as vscode from 'vscode';

/**
 * Log levels following standard severity ordering
 */
export enum LogLevel {
  DEBUG = 0,
  INFO = 1,
  WARN = 2,
  ERROR = 3,
}

/**
 * Configuration for logger behavior
 */
export interface LoggerConfig {
  /** Current log level - only messages at this level or higher will be logged */
  level: LogLevel;
  /** Whether to show error/warning messages to user via VS Code UI */
  showUserMessages: boolean;
  /** Optional custom output channel name */
  channelName?: string;
}

/**
 * Logger interface for component loggers
 */
export interface ILogger {
  debug(message: string, ...args: any[]): void;
  info(message: string, ...args: any[]): void;
  warn(message: string, showUser?: boolean, ...args: any[]): void;
  error(message: string, error?: Error, showUser?: boolean): void;
  show(): void;
  clear(): void;
  setConfig(config: Partial<LoggerConfig>): void;
  getConfig(): LoggerConfig;
}

/**
 * Centralized logger for OpenQC extension
 *
 * Provides consistent logging interface with:
 * - Configurable log levels
 * - Output to VS Code Output panel
 * - Optional user notifications for warnings/errors
 * - Structured log format with timestamps
 * - Stack trace preservation for errors
 *
 * @example
 * ```typescript
 * // Basic usage
 * Logger.getInstance().info('LSP started');
 * Logger.getInstance().error('Failed to load jobs', error, true);
 *
 * // With configuration
 * const logger = Logger.getInstance();
 * logger.setConfig({ level: LogLevel.DEBUG, showUserMessages: true });
 * logger.debug('Detailed debug info');
 * ```
 */
export class Logger implements ILogger {
  private static instance: Logger;
  private config: LoggerConfig;
  private outputChannel: vscode.OutputChannel;

  /**
   * Private constructor to enforce singleton pattern
   */
  private constructor(config?: Partial<LoggerConfig>) {
    this.config = {
      level: LogLevel.INFO,
      showUserMessages: true,
      channelName: config?.channelName || 'OpenQC',
    };
    this.outputChannel = vscode.window.createOutputChannel(this.config.channelName || 'OpenQC');
  }

  /**
   * Get the singleton Logger instance
   *
   * @param config - Optional configuration to apply on first access
   * @returns The singleton Logger instance
   */
  static getInstance(config?: Partial<LoggerConfig>): Logger {
    if (!Logger.instance) {
      Logger.instance = new Logger(config);
    }
    return Logger.instance;
  }

  /**
   * Update logger configuration
   *
   * @param config - New configuration values
   */
  setConfig(config: Partial<LoggerConfig>): void {
    this.config = { ...this.config, ...config };
    if (config.channelName && config.channelName !== this.config.channelName) {
      this.outputChannel.dispose();
      this.config.channelName = config.channelName;
      this.outputChannel = vscode.window.createOutputChannel(config.channelName);
    }
  }

  /**
   * Get current logger configuration
   *
   * @returns Current configuration
   */
  getConfig(): LoggerConfig {
    return { ...this.config };
  }

  /**
   * Log debug-level message
   *
   * Use for detailed diagnostic information that users typically don't need to see
   *
   * @param message - Log message
   * @param args - Additional arguments to log
   */
  debug(message: string, ...args: any[]): void {
    if (this.config.level <= LogLevel.DEBUG) {
      this.log('DEBUG', message, ...args);
    }
  }

  /**
   * Log info-level message
   *
   * Use for general informational messages about normal operation
   *
   * @param message - Log message
   * @param args - Additional arguments to log
   */
  info(message: string, ...args: any[]): void {
    if (this.config.level <= LogLevel.INFO) {
      this.log('INFO', message, ...args);
    }
  }

  /**
   * Log warning message
   *
   * Use for potentially harmful situations that don't prevent execution
   *
   * @param message - Log message
   * @param showUser - Whether to show warning to user (default: false)
   * @param args - Additional arguments to log
   */
  warn(message: string, showUser = false, ...args: any[]): void {
    if (this.config.level <= LogLevel.WARN) {
      this.log('WARN', message, ...args);
    }
    if (showUser && this.config.showUserMessages) {
      vscode.window.showWarningMessage(message);
    }
  }

  /**
   * Log error message
   *
   * Use for error events that might still allow the application to continue
   *
   * @param message - Log message
   * @param error - Optional Error object with stack trace
   * @param showUser - Whether to show error to user (default: false)
   */
  error(message: string, error?: Error, showUser = false): void {
    if (this.config.level <= LogLevel.ERROR) {
      this.log('ERROR', message, error?.message || '');
      if (error?.stack) {
        this.outputChannel.appendLine(error.stack);
      }
    }
    if (showUser && this.config.showUserMessages) {
      vscode.window.showErrorMessage(message);
    }
  }

  /**
   * Internal logging method that formats and writes messages
   *
   * @param level - Log level string
   * @param message - Log message
   * @param args - Additional arguments to serialize and append
   */
  private log(level: string, message: string, ...args: any[]): void {
    const timestamp = new Date().toISOString();
    const argsStr = args.length > 0 ? ' ' + args.map(a => {
      try {
        return JSON.stringify(a, null, 2);
      } catch {
        return String(a);
      }
    }).join(' ') : '';
    this.outputChannel.appendLine(`[${timestamp}] [${level}] ${message}${argsStr}`);
  }

  /**
   * Show the Output panel to display logs to user
   */
  show(): void {
    this.outputChannel.show();
  }

  /**
   * Clear all log content from the Output panel
   */
  clear(): void {
    this.outputChannel.clear();
  }

  /**
   * Dispose the Output panel
   *
   * Call this when the extension deactivates
   */
  dispose(): void {
    this.outputChannel.dispose();
  }

  /**
   * Convert string log level to enum
   *
   * @param level - String log level ('debug', 'info', 'warn', 'error')
   * @returns LogLevel enum value
   */
  static parseLogLevel(level: string): LogLevel {
    switch (level.toLowerCase()) {
      case 'debug': return LogLevel.DEBUG;
      case 'info': return LogLevel.INFO;
      case 'warn': return LogLevel.WARN;
      case 'error': return LogLevel.ERROR;
      default: return LogLevel.INFO;
    }
  }
}

/**
 * Create a component-specific logger with automatic prefixing
 *
 * @param component - Component name for log prefixes
 * @returns Logger instance that automatically prefixes messages
 *
 * @example
 * ```typescript
 * const logger = createComponentLogger('LSPDiscovery');
 * logger.info('Fetching repositories'); // Logs: [INFO] [LSPDiscovery] Fetching repositories
 * ```
 */
export function createComponentLogger(component: string): ILogger {
  const baseLogger = Logger.getInstance();
  return {
    debug: (message: string, ...args: any[]) => baseLogger.debug(`[${component}] ${message}`, ...args),
    info: (message: string, ...args: any[]) => baseLogger.info(`[${component}] ${message}`, ...args),
    warn: (message: string, showUser = false, ...args: any[]) => baseLogger.warn(`[${component}] ${message}`, showUser, ...args),
    error: (message: string, error?: Error, showUser = false) => baseLogger.error(`[${component}] ${message}`, error, showUser),
    show: () => baseLogger.show(),
    clear: () => baseLogger.clear(),
    setConfig: (config: Partial<LoggerConfig>) => baseLogger.setConfig(config),
    getConfig: () => baseLogger.getConfig(),
  };
}
