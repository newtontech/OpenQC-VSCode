import * as vscode from 'vscode';
import {
  getLanguageFeaturePolicy,
  normalizeLanguageFeatureMode,
  readLanguageFeatureMode,
} from '../../src/languageFeatures';

describe('language feature mode', () => {
  it('normalizes supported modes', () => {
    expect(normalizeLanguageFeatureMode('lsp')).toBe('lsp');
    expect(normalizeLanguageFeatureMode('local')).toBe('local');
    expect(normalizeLanguageFeatureMode('hybrid')).toBe('hybrid');
  });

  it('falls back to lsp for invalid or missing modes', () => {
    expect(normalizeLanguageFeatureMode('invalid')).toBe('lsp');
    expect(normalizeLanguageFeatureMode(undefined)).toBe('lsp');
    expect(normalizeLanguageFeatureMode(42)).toBe('lsp');
  });

  it('reads openqc.languageFeatures.mode from VS Code configuration', () => {
    const get = jest.fn().mockReturnValue('local');
    (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({ get });

    expect(readLanguageFeatureMode()).toBe('local');
    expect(vscode.workspace.getConfiguration).toHaveBeenCalledWith('openqc.languageFeatures');
    expect(get).toHaveBeenCalledWith('mode', 'lsp');
  });

  it('defaults configuration reads to lsp when the configured value is invalid', () => {
    (vscode.workspace.getConfiguration as jest.Mock).mockReturnValue({
      get: jest.fn().mockReturnValue('always-local'),
    });

    expect(readLanguageFeatureMode()).toBe('lsp');
  });
});

describe('language feature registration policy', () => {
  it('makes LSP the default owner in lsp mode', () => {
    expect(getLanguageFeaturePolicy('lsp')).toEqual({
      mode: 'lsp',
      autoStartLsp: true,
      registerLocalProviders: false,
      registerLocalDiagnostics: false,
      validateWithLocalDiagnostics: false,
    });
  });

  it('preserves local provider behavior in local mode without LSP auto-start', () => {
    expect(getLanguageFeaturePolicy('local')).toEqual({
      mode: 'local',
      autoStartLsp: false,
      registerLocalProviders: true,
      registerLocalDiagnostics: true,
      validateWithLocalDiagnostics: true,
    });
  });

  it('keeps hybrid conservative to avoid duplicate diagnostics and providers', () => {
    expect(getLanguageFeaturePolicy('hybrid')).toEqual({
      mode: 'hybrid',
      autoStartLsp: true,
      registerLocalProviders: false,
      registerLocalDiagnostics: false,
      validateWithLocalDiagnostics: false,
    });
  });
});
