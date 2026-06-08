import { Logger, LogLevel, createComponentLogger } from '../../src/utils/Logger';

// Mock vscode module
jest.mock('vscode', () => ({
  window: {
    createOutputChannel: jest.fn(() => ({
      appendLine: jest.fn(),
      append: jest.fn(),
      clear: jest.fn(),
      show: jest.fn(),
      hide: jest.fn(),
      dispose: jest.fn(),
    })),
    showInformationMessage: jest.fn(),
    showWarningMessage: jest.fn(),
    showErrorMessage: jest.fn(),
  },
  workspace: {
    getConfiguration: jest.fn(() => ({
      get: jest.fn((key: string, defaultValue: any) => defaultValue),
    })),
  },
}));

describe('Logger', () => {
  beforeEach(() => {
    Logger.resetInstance();
  });

  afterEach(() => {
    Logger.resetInstance();
  });

  describe('singleton', () => {
    it('should return the same instance', () => {
      const instance1 = Logger.getInstance();
      const instance2 = Logger.getInstance();
      expect(instance1).toBe(instance2);
    });

    it('should create a new instance after reset', () => {
      const instance1 = Logger.getInstance();
      Logger.resetInstance();
      const instance2 = Logger.getInstance();
      expect(instance1).not.toBe(instance2);
    });
  });

  describe('log levels', () => {
    it('should parse log level strings', () => {
      expect(Logger.parseLogLevel('debug')).toBe(LogLevel.DEBUG);
      expect(Logger.parseLogLevel('info')).toBe(LogLevel.INFO);
      expect(Logger.parseLogLevel('warn')).toBe(LogLevel.WARN);
      expect(Logger.parseLogLevel('error')).toBe(LogLevel.ERROR);
      expect(Logger.parseLogLevel('unknown')).toBe(LogLevel.INFO);
    });

    it('should respect log level filtering', () => {
      const logger = Logger.getInstance();
      logger.setConfig({ level: LogLevel.WARN });

      const channel = {
        appendLine: jest.fn(),
        append: jest.fn(),
        clear: jest.fn(),
        show: jest.fn(),
        hide: jest.fn(),
        dispose: jest.fn(),
      };

      // Access internal output channel via reflection
      (logger as any).outputChannel = channel;

      logger.debug('should not log');
      logger.info('should not log');
      logger.warn('should log');
      logger.error('should log');

      const debugCalls = channel.appendLine.mock.calls.filter((call: string[]) =>
        call[0].includes('DEBUG')
      );
      const infoCalls = channel.appendLine.mock.calls.filter((call: string[]) =>
        call[0].includes('INFO')
      );
      expect(debugCalls).toHaveLength(0);
      expect(infoCalls).toHaveLength(0);
    });
  });

  describe('configuration', () => {
    it('should get current config', () => {
      const logger = Logger.getInstance();
      const config = logger.getConfig();
      expect(config).toHaveProperty('level');
      expect(config).toHaveProperty('showUserMessages');
    });

    it('should update config', () => {
      const logger = Logger.getInstance();
      logger.setConfig({ level: LogLevel.DEBUG, showUserMessages: false });
      const config = logger.getConfig();
      expect(config.level).toBe(LogLevel.DEBUG);
      expect(config.showUserMessages).toBe(false);
    });
  });

  describe('logging methods', () => {
    let logger: Logger;
    let channel: {
      appendLine: jest.Mock;
      append: jest.Mock;
      clear: jest.Mock;
      show: jest.Mock;
      hide: jest.Mock;
      dispose: jest.Mock;
    };

    beforeEach(() => {
      logger = Logger.getInstance();
      channel = {
        appendLine: jest.fn(),
        append: jest.fn(),
        clear: jest.fn(),
        show: jest.fn(),
        hide: jest.fn(),
        dispose: jest.fn(),
      };
      (logger as any).outputChannel = channel;
      logger.setConfig({ level: LogLevel.DEBUG, showUserMessages: false });
    });

    it('should log debug messages', () => {
      logger.debug('test debug message');
      expect(channel.appendLine).toHaveBeenCalledWith(
        expect.stringContaining('[DEBUG] test debug message')
      );
    });

    it('should log info messages', () => {
      logger.info('test info message');
      expect(channel.appendLine).toHaveBeenCalledWith(
        expect.stringContaining('[INFO] test info message')
      );
    });

    it('should log warn messages', () => {
      logger.warn('test warning');
      expect(channel.appendLine).toHaveBeenCalledWith(
        expect.stringContaining('[WARN] test warning')
      );
    });

    it('should log error messages with stack trace', () => {
      const error = new Error('test error');
      error.stack = 'Error: test error\n  at test';
      logger.error('something failed', error);
      expect(channel.appendLine).toHaveBeenCalledWith(
        expect.stringContaining('[ERROR] something failed')
      );
      expect(channel.appendLine).toHaveBeenCalledWith(error.stack);
    });

    it('should log error messages without error object', () => {
      logger.error('something failed');
      expect(channel.appendLine).toHaveBeenCalledWith(
        expect.stringContaining('[ERROR] something failed')
      );
    });

    it('should include timestamp in log messages', () => {
      logger.info('test');
      const call = channel.appendLine.mock.calls[0][0] as string;
      expect(call).toMatch(/^\[\d{4}-\d{2}-\d{2}T/);
    });

    it('should serialize extra arguments', () => {
      logger.info('message', { key: 'value' });
      expect(channel.appendLine).toHaveBeenCalledWith(expect.stringContaining('"key": "value"'));
    });

    it('should handle non-serializable arguments gracefully', () => {
      const circular: any = { ref: null };
      circular.ref = circular;
      logger.info('message', circular);
      expect(channel.appendLine).toHaveBeenCalled();
    });
  });

  describe('show and clear', () => {
    it('should show output channel', () => {
      const logger = Logger.getInstance();
      const channel = {
        appendLine: jest.fn(),
        append: jest.fn(),
        clear: jest.fn(),
        show: jest.fn(),
        hide: jest.fn(),
        dispose: jest.fn(),
      };
      (logger as any).outputChannel = channel;

      logger.show();
      expect(channel.show).toHaveBeenCalled();
    });

    it('should clear output channel', () => {
      const logger = Logger.getInstance();
      const channel = {
        appendLine: jest.fn(),
        append: jest.fn(),
        clear: jest.fn(),
        show: jest.fn(),
        hide: jest.fn(),
        dispose: jest.fn(),
      };
      (logger as any).outputChannel = channel;

      logger.clear();
      expect(channel.clear).toHaveBeenCalled();
    });
  });

  describe('createComponentLogger', () => {
    it('should prefix messages with component name', () => {
      const logger = Logger.getInstance();
      const channel = {
        appendLine: jest.fn(),
        append: jest.fn(),
        clear: jest.fn(),
        show: jest.fn(),
        hide: jest.fn(),
        dispose: jest.fn(),
      };
      (logger as any).outputChannel = channel;
      logger.setConfig({ level: LogLevel.DEBUG, showUserMessages: false });

      const componentLogger = createComponentLogger('TestComponent');
      componentLogger.info('test message');

      expect(channel.appendLine).toHaveBeenCalledWith(
        expect.stringContaining('[TestComponent] test message')
      );
    });

    it('should implement ILogger interface', () => {
      const componentLogger = createComponentLogger('Test');
      expect(typeof componentLogger.debug).toBe('function');
      expect(typeof componentLogger.info).toBe('function');
      expect(typeof componentLogger.warn).toBe('function');
      expect(typeof componentLogger.error).toBe('function');
      expect(typeof componentLogger.show).toBe('function');
      expect(typeof componentLogger.clear).toBe('function');
      expect(typeof componentLogger.setConfig).toBe('function');
      expect(typeof componentLogger.getConfig).toBe('function');
    });
  });

  describe('user notifications', () => {
    it('should show warning to user when showUser is true', () => {
      const vscode = require('vscode');
      const logger = Logger.getInstance();
      logger.setConfig({ level: LogLevel.WARN, showUserMessages: true });

      logger.warn('user warning', true);
      expect(vscode.window.showWarningMessage).toHaveBeenCalledWith('user warning');
    });

    it('should show error to user when showUser is true', () => {
      const vscode = require('vscode');
      const logger = Logger.getInstance();
      logger.setConfig({ level: LogLevel.ERROR, showUserMessages: true });

      logger.error('user error', undefined, true);
      expect(vscode.window.showErrorMessage).toHaveBeenCalledWith('user error');
    });

    it('should not show user messages when showUserMessages is false', () => {
      const vscode = require('vscode');
      const logger = Logger.getInstance();
      logger.setConfig({ level: LogLevel.ERROR, showUserMessages: false });

      logger.error('hidden error', undefined, true);
      expect(vscode.window.showErrorMessage).not.toHaveBeenCalled();
    });
  });
});
