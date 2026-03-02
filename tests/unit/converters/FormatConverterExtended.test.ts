/**
 * Extended tests for FormatConverter to improve coverage
 * Tests for quickConvert, showSetupInstructions, and error paths
 */

jest.setTimeout(15000);

// Mock util.promisify - must be before imports
const mockExecAsync = jest.fn();
jest.mock('util', () => {
  const actualUtil = jest.requireActual('util');
  return {
    ...actualUtil,
    promisify: jest.fn(() => mockExecAsync),
  };
});

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
  exec: jest.fn(),
  spawn: jest.fn(),
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
      mockExecAsync.mockRejectedValueOnce(new Error('Python not found'));

      await quickConvert(SupportedFormat.XYZ);

      expect(mockShowWarningMessage).toHaveBeenCalled();
    });

    it('should show error message when conversion fails', async () => {
      // checkBackend passes
      mockExecAsync.mockResolvedValueOnce({ stdout: 'Python 3.9', stderr: '' });
      mockExecAsync.mockResolvedValueOnce({ stdout: '', stderr: '' });
      // conversion fails
      mockExecAsync.mockResolvedValueOnce({
        stdout: JSON.stringify({ success: false, error: 'Test error' }),
        stderr: '',
      });

      await quickConvert(SupportedFormat.XYZ);

      expect(mockShowErrorMessage).toHaveBeenCalledWith('Conversion failed: Test error');
    });

    it('should handle successful conversion with Open File action', async () => {
      // checkBackend passes
      mockExecAsync.mockResolvedValueOnce({ stdout: 'Python 3.9', stderr: '' });
      mockExecAsync.mockResolvedValueOnce({ stdout: '', stderr: '' });
      // conversion succeeds
      mockExecAsync.mockResolvedValueOnce({
        stdout: JSON.stringify({ success: true, output_file: '/test/output.xyz' }),
        stderr: '',
      });

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
      mockExecAsync.mockResolvedValueOnce({ stdout: 'Python 3.9', stderr: '' });
      mockExecAsync.mockResolvedValueOnce({ stdout: '', stderr: '' });
      // conversion succeeds
      mockExecAsync.mockResolvedValueOnce({
        stdout: JSON.stringify({ success: true, output_file: '/test/output.xyz' }),
        stderr: '',
      });

      mockShowInformationMessage.mockResolvedValue('Copy to Clipboard');
      mockReadFile.mockResolvedValue(Buffer.from('test content'));

      await quickConvert(SupportedFormat.XYZ);

      expect(mockWriteText).toHaveBeenCalledWith('test content');
    });

    it('should handle successful conversion with no action (dismiss)', async () => {
      // checkBackend passes
      mockExecAsync.mockResolvedValueOnce({ stdout: 'Python 3.9', stderr: '' });
      mockExecAsync.mockResolvedValueOnce({ stdout: '', stderr: '' });
      // conversion succeeds
      mockExecAsync.mockResolvedValueOnce({
        stdout: JSON.stringify({ success: true, output_file: '/test/output.xyz' }),
        stderr: '',
      });

      mockShowInformationMessage.mockResolvedValue(undefined);

      await quickConvert(SupportedFormat.XYZ);

      expect(mockShowInformationMessage).toHaveBeenCalled();
    });
  });

  describe('Error handling paths', () => {
    it('should handle JSON parse errors in convert', async () => {
      mockExecAsync.mockResolvedValueOnce({
        stdout: 'invalid json',
        stderr: '',
      });

      const result = await converter.convert('/test/input', '/test/output');

      expect(result.success).toBe(false);
      expect(result.error).toBeDefined();
    });

    it('should handle empty output in batch convert', async () => {
      mockExecAsync.mockResolvedValueOnce({
        stdout: '',
        stderr: '',
      });

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
