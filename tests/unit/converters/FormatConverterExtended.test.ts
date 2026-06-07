/**
 * Extended tests for FormatConverter to improve coverage
 * Tests for quickConvert, showSetupInstructions, and error paths
 */

jest.setTimeout(15000);

const mockSpawn = jest.fn();

// Mock vscode module FIRST before imports
jest.mock(
  'vscode',
  () => {
    const mockShowWarningMessage = jest.fn();
    const mockCreateTerminal = jest.fn();
    const mockShowInformationMessage = jest.fn();
    const mockShowErrorMessage = jest.fn();
    const mockOpenTextDocument = jest.fn();
    const mockShowTextDocument = jest.fn();
    const mockReadFile = jest.fn();
    const mockWriteText = jest.fn();

    return {
      window: {
        createOutputChannel: jest.fn(() => ({
          appendLine: jest.fn(),
          dispose: jest.fn(),
        })),
        showWarningMessage: mockShowWarningMessage,
        showInformationMessage: mockShowInformationMessage,
        showErrorMessage: mockShowErrorMessage,
        createTerminal: mockCreateTerminal,
        showTextDocument: mockShowTextDocument,
        activeTextEditor: {
          document: {
            uri: { fsPath: '/test/input.xyz' },
          },
        },
      },
      workspace: {
        openTextDocument: mockOpenTextDocument,
        fs: {
          readFile: mockReadFile,
        },
      },
      env: {
        clipboard: {
          writeText: mockWriteText,
        },
      },
      Uri: {
        file: jest.fn(path => ({ fsPath: path })),
      },
      // Export mocks for access in tests
      __mocks__: {
        mockShowWarningMessage,
        mockCreateTerminal,
        mockShowInformationMessage,
        mockShowErrorMessage,
        mockOpenTextDocument,
        mockShowTextDocument,
        mockReadFile,
        mockWriteText,
      },
    };
  },
  { virtual: true }
);

// Mock child_process
jest.mock('child_process', () => ({
  spawn: mockSpawn,
}));

import { FormatConverter, SupportedFormat, quickConvert } from '../../../src/converters';

// Get mocks from the vscode module
const vscode = require('vscode');
const {
  mockShowWarningMessage,
  mockCreateTerminal,
  mockShowInformationMessage,
  mockShowErrorMessage,
  mockOpenTextDocument,
  mockShowTextDocument,
  mockReadFile,
  mockWriteText,
} = vscode.__mocks__;

describe('FormatConverter Extended Tests', () => {
  let converter: FormatConverter;

  function mockPythonResult(stdout = '', stderr = '', code = 0): void {
    mockSpawn.mockImplementationOnce(() => {
      const { EventEmitter } = require('events');
      const child = new EventEmitter();
      child.stdout = new EventEmitter();
      child.stderr = new EventEmitter();

      process.nextTick(() => {
        if (stdout) {
          child.stdout.emit('data', Buffer.from(stdout));
        }
        if (stderr) {
          child.stderr.emit('data', Buffer.from(stderr));
        }
        child.emit('close', code);
      });

      return child;
    });
  }

  function mockPythonError(error: Error): void {
    mockSpawn.mockImplementationOnce(() => {
      const { EventEmitter } = require('events');
      const child = new EventEmitter();
      child.stdout = new EventEmitter();
      child.stderr = new EventEmitter();

      process.nextTick(() => {
        child.emit('error', error);
      });

      return child;
    });
  }

  beforeEach(() => {
    converter = new FormatConverter();
    jest.clearAllMocks();
  });

  afterEach(() => {
    converter.dispose();
  });

  describe('showSetupInstructions', () => {
    it('should show warning and create terminal on Install action', async () => {
      const mockTerminal = {
        sendText: jest.fn(),
        show: jest.fn(),
      };
      mockShowWarningMessage.mockResolvedValue('Install');
      mockCreateTerminal.mockReturnValue(mockTerminal);

      await converter.showSetupInstructions();

      expect(mockShowWarningMessage).toHaveBeenCalledWith(
        'OpenQC Format Converter requires Python and dpdata. Install now?',
        'Install',
        'Cancel'
      );
      expect(mockCreateTerminal).toHaveBeenCalledWith('OpenQC Setup');
      expect(mockTerminal.sendText).toHaveBeenCalledWith('pip install dpdata');
      expect(mockTerminal.show).toHaveBeenCalled();
    });

    it('should not create terminal on Cancel action', async () => {
      mockShowWarningMessage.mockResolvedValue('Cancel');

      await converter.showSetupInstructions();

      expect(mockCreateTerminal).not.toHaveBeenCalled();
    });

    it('should not create terminal when undefined action', async () => {
      mockShowWarningMessage.mockResolvedValue(undefined);

      await converter.showSetupInstructions();

      expect(mockCreateTerminal).not.toHaveBeenCalled();
    });
  });

  describe('quickConvert helper function', () => {
    it('should show setup instructions when backend not available', async () => {
      mockPythonError(new Error('Python not found'));

      await quickConvert(SupportedFormat.XYZ);

      expect(mockShowWarningMessage).toHaveBeenCalled();
    });

    it('should show error message when conversion fails', async () => {
      // checkBackend passes
      mockPythonResult('Python 3.9');
      mockPythonResult();
      // conversion fails
      mockPythonResult(JSON.stringify({ success: false, error: 'Test error' }));

      await quickConvert(SupportedFormat.XYZ);

      expect(mockShowErrorMessage).toHaveBeenCalledWith('Conversion failed: Test error');
    });

    it('should handle successful conversion with Open File action', async () => {
      // checkBackend passes
      mockPythonResult('Python 3.9');
      mockPythonResult();
      // conversion succeeds
      mockPythonResult(JSON.stringify({ success: true, output_file: '/test/output.xyz' }));

      mockShowInformationMessage.mockResolvedValue('Open File');
      mockOpenTextDocument.mockResolvedValue({});
      mockShowTextDocument.mockResolvedValue({});

      await quickConvert(SupportedFormat.XYZ);

      expect(mockShowInformationMessage).toHaveBeenCalled();
      expect(mockOpenTextDocument).toHaveBeenCalled();
      expect(mockShowTextDocument).toHaveBeenCalled();
    });

    it('should handle successful conversion with Copy to Clipboard action', async () => {
      // checkBackend passes
      mockPythonResult('Python 3.9');
      mockPythonResult();
      // conversion succeeds
      mockPythonResult(JSON.stringify({ success: true, output_file: '/test/output.xyz' }));

      mockShowInformationMessage.mockResolvedValue('Copy to Clipboard');
      mockReadFile.mockResolvedValue(Buffer.from('test content'));

      await quickConvert(SupportedFormat.XYZ);

      expect(mockWriteText).toHaveBeenCalledWith('test content');
    });

    it('should handle successful conversion with no action (dismiss)', async () => {
      // checkBackend passes
      mockPythonResult('Python 3.9');
      mockPythonResult();
      // conversion succeeds
      mockPythonResult(JSON.stringify({ success: true, output_file: '/test/output.xyz' }));

      mockShowInformationMessage.mockResolvedValue(undefined);

      await quickConvert(SupportedFormat.XYZ);

      expect(mockShowInformationMessage).toHaveBeenCalled();
    });
  });

  describe('Error handling paths', () => {
    it('should handle JSON parse errors in convert', async () => {
      mockPythonResult('invalid json');

      const result = await converter.convert('/test/input', '/test/output');

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should handle empty output in batch convert', async () => {
      mockPythonResult();

      const result = await converter.batchConvert(
        ['/test/input.xyz'],
        '/test/output',
        SupportedFormat.VASP
      );

      expect(result.success).toBe(false);
    });
  });

  describe('Configuration options', () => {
    it('should use custom python path', () => {
      const customConverter = new FormatConverter({ pythonPath: '/custom/python' });
      expect(customConverter).toBeDefined();
      customConverter.dispose();
    });

    it('should handle preserveMetadata option', () => {
      const noMetaConverter = new FormatConverter({ preserveMetadata: false });
      expect(noMetaConverter).toBeDefined();
      noMetaConverter.dispose();
    });
  });
});
