import * as vscode from 'vscode';

export const LANGUAGE_FEATURE_MODES = ['lsp', 'local', 'hybrid'] as const;

export type LanguageFeatureMode = (typeof LANGUAGE_FEATURE_MODES)[number];

export interface LanguageFeaturePolicy {
  mode: LanguageFeatureMode;
  autoStartLsp: boolean;
  registerLocalProviders: boolean;
  registerLocalDiagnostics: boolean;
  validateWithLocalDiagnostics: boolean;
}

const DEFAULT_LANGUAGE_FEATURE_MODE: LanguageFeatureMode = 'lsp';

export function isLanguageFeatureMode(value: unknown): value is LanguageFeatureMode {
  return typeof value === 'string' && (LANGUAGE_FEATURE_MODES as readonly string[]).includes(value);
}

export function normalizeLanguageFeatureMode(value: unknown): LanguageFeatureMode {
  return isLanguageFeatureMode(value) ? value : DEFAULT_LANGUAGE_FEATURE_MODE;
}

export function readLanguageFeatureMode(): LanguageFeatureMode {
  const config = vscode.workspace.getConfiguration('openqc.languageFeatures');
  return normalizeLanguageFeatureMode(config.get<unknown>('mode', DEFAULT_LANGUAGE_FEATURE_MODE));
}

export function getLanguageFeaturePolicy(mode: LanguageFeatureMode): LanguageFeaturePolicy {
  if (mode === 'local') {
    return {
      mode,
      autoStartLsp: false,
      registerLocalProviders: true,
      registerLocalDiagnostics: true,
      validateWithLocalDiagnostics: true,
    };
  }

  return {
    mode,
    autoStartLsp: true,
    registerLocalProviders: false,
    registerLocalDiagnostics: false,
    validateWithLocalDiagnostics: false,
  };
}
