/**
 * Compatibility Report Generator
 *
 * Generates a comprehensive compatibility report that checks alignment between
 * the LSP registry, package.json contributions, docs, and rule manifests.
 * Extends the basic registry alignment checks with rule manifest integration.
 *
 * @module smoke/compatibilityReport
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/163
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/167
 */

import * as fs from 'fs';
import * as path from 'path';
import { listBundledLspServers, getLspDiagnosticReadiness } from '../lsp/registry';
import type { LSPServerRegistryEntry } from '../lsp/types';
import type {
  CompatibilityCheck,
  CompatibilityReport,
  CheckStatus,
  LspCompatibilityEntry,
  RuleManifest,
} from './types';
import { tryReadManifest } from './ruleManifestReader';

// ---------------------------------------------------------------------------
// Required dependencies for VSIX verification
// ---------------------------------------------------------------------------

const REQUIRED_VSIX_DEPENDENCIES = ['vscode-languageclient', 'vscode-languageserver-protocol'];

// ---------------------------------------------------------------------------
// Report generation
// ---------------------------------------------------------------------------

/**
 * Generate a full compatibility report checking all LSP server integrations.
 *
 * @param repoRoot - Root directory of the OpenQC-VSCode repository.
 * @param manifestDir - Optional directory containing rule manifest JSON files.
 * @returns A complete compatibility report.
 */
export function generateCompatibilityReport(
  repoRoot: string,
  manifestDir?: string
): CompatibilityReport {
  const servers = listBundledLspServers();
  const pkgPath = path.join(repoRoot, 'package.json');
  const docPath = path.join(repoRoot, 'docs', 'LSP_COMPATIBILITY.md');

  const packageJson = loadJsonSafe(pkgPath);
  const compatDoc = loadTextSafe(docPath);

  const entries: LspCompatibilityEntry[] = servers.map(server => {
    const checks: CompatibilityCheck[] = [
      checkRegistryEntryExists(server),
      checkPackageJsonLanguage(server, packageJson),
      checkPackageJsonConfiguration(server, packageJson),
      checkGrammarExists(server, repoRoot),
      checkLanguageConfigExists(server, repoRoot),
      checkDocsAlignment(server, compatDoc),
      checkDiagnosticReadiness(server),
      checkRuleManifest(server, manifestDir),
      checkLocalLaunchConfig(server),
      checkCapabilityManifest(server, repoRoot),
    ];

    const passed = checks.every(c => c.status === 'pass' || c.status === 'skip');

    return {
      serverId: server.id,
      serverName: server.name,
      languageId: server.languageId,
      stability: server.stability,
      checks,
      passed,
    };
  });

  const passedServers = entries.filter(e => e.passed).length;

  return {
    generatedAt: new Date().toISOString(),
    totalServers: servers.length,
    passedServers,
    failedServers: entries.length - passedServers,
    entries,
  };
}

// ---------------------------------------------------------------------------
// Individual checks
// ---------------------------------------------------------------------------

function checkRegistryEntryExists(server: LSPServerRegistryEntry): CompatibilityCheck {
  return {
    name: 'registry-entry-exists',
    description: `Registry entry exists for ${server.name}`,
    status: server.id && server.name ? 'pass' : 'fail',
    detail: server.id ? undefined : 'Missing required registry fields',
  };
}

function checkPackageJsonLanguage(
  server: LSPServerRegistryEntry,
  packageJson: Record<string, unknown> | null
): CompatibilityCheck {
  if (!packageJson) {
    return {
      name: 'package-json-language',
      description: `package.json language contribution for ${server.name}`,
      status: 'fail',
      detail: 'package.json could not be loaded',
    };
  }

  const languages = getPackageJsonLanguages(packageJson);
  const match = languages.find((l: Record<string, unknown>) => l.id === server.languageId);

  if (!match) {
    return {
      name: 'package-json-language',
      description: `package.json language contribution for ${server.name}`,
      status: 'fail',
      detail: `No language contribution found for languageId "${server.languageId}"`,
    };
  }

  return {
    name: 'package-json-language',
    description: `package.json language contribution for ${server.name}`,
    status: 'pass',
  };
}

function checkPackageJsonConfiguration(
  server: LSPServerRegistryEntry,
  packageJson: Record<string, unknown> | null
): CompatibilityCheck {
  if (!packageJson) {
    return {
      name: 'package-json-configuration',
      description: `package.json configuration settings for ${server.name}`,
      status: 'fail',
      detail: 'package.json could not be loaded',
    };
  }

  const properties = getPackageJsonProperties(packageJson);
  const requiredSettings = [
    `openqc.lsp.${server.languageId}.enabled`,
    `openqc.lsp.${server.languageId}.path`,
    `openqc.lsp.${server.languageId}.command`,
  ];

  const missing = requiredSettings.filter(key => !(key in properties));

  if (missing.length > 0) {
    return {
      name: 'package-json-configuration',
      description: `package.json configuration settings for ${server.name}`,
      status: 'warn',
      detail: `Missing configuration keys: ${missing.join(', ')}`,
    };
  }

  return {
    name: 'package-json-configuration',
    description: `package.json configuration settings for ${server.name}`,
    status: 'pass',
  };
}

function checkGrammarExists(server: LSPServerRegistryEntry, repoRoot: string): CompatibilityCheck {
  const grammarDir = path.join(repoRoot, 'syntaxes');
  const expectedName = `${server.languageId}.tmLanguage.json`;
  const grammarPath = path.join(grammarDir, expectedName);

  if (fs.existsSync(grammarPath)) {
    return {
      name: 'grammar-exists',
      description: `TextMate grammar for ${server.name}`,
      status: 'pass',
    };
  }

  // Check for YAML variant
  const yamlPath = path.join(grammarDir, `${server.languageId}.tmLanguage.yaml`);
  if (fs.existsSync(yamlPath)) {
    return {
      name: 'grammar-exists',
      description: `TextMate grammar for ${server.name}`,
      status: 'pass',
    };
  }

  // Check for plist variant
  const plistPath = path.join(grammarDir, `${server.languageId}.tmLanguage`);
  if (fs.existsSync(plistPath)) {
    return {
      name: 'grammar-exists',
      description: `TextMate grammar for ${server.name}`,
      status: 'pass',
    };
  }

  return {
    name: 'grammar-exists',
    description: `TextMate grammar for ${server.name}`,
    status: 'warn',
    detail: `No grammar found at syntaxes/${expectedName} (checked .json, .yaml, .plist)`,
  };
}

function checkLanguageConfigExists(
  server: LSPServerRegistryEntry,
  repoRoot: string
): CompatibilityCheck {
  const configDir = path.join(repoRoot, 'language-configurations');
  const expectedName = `${server.languageId}.json`;
  const configPath = path.join(configDir, expectedName);

  if (fs.existsSync(configPath)) {
    return {
      name: 'language-config-exists',
      description: `Language configuration for ${server.name}`,
      status: 'pass',
    };
  }

  return {
    name: 'language-config-exists',
    description: `Language configuration for ${server.name}`,
    status: 'warn',
    detail: `No language config at language-configurations/${expectedName}`,
  };
}

function checkDocsAlignment(
  server: LSPServerRegistryEntry,
  docContent: string | null
): CompatibilityCheck {
  if (!docContent) {
    return {
      name: 'docs-alignment',
      description: `LSP_COMPATIBILITY.md entry for ${server.name}`,
      status: 'fail',
      detail: 'docs/LSP_COMPATIBILITY.md could not be read',
    };
  }

  const hasRepo = docContent.includes(server.repository);
  const hasLanguageId = docContent.includes(`\`${server.languageId}\``);
  const hasBranch = docContent.includes(`\`${server.defaultBranch}\``);

  if (!hasRepo || !hasLanguageId || !hasBranch) {
    const missing: string[] = [];
    if (!hasRepo) {
      missing.push(`repository "${server.repository}"`);
    }
    if (!hasLanguageId) {
      missing.push(`languageId \`${server.languageId}\``);
    }
    if (!hasBranch) {
      missing.push(`branch \`${server.defaultBranch}\``);
    }

    return {
      name: 'docs-alignment',
      description: `LSP_COMPATIBILITY.md entry for ${server.name}`,
      status: 'fail',
      detail: `Missing in docs: ${missing.join(', ')}`,
    };
  }

  return {
    name: 'docs-alignment',
    description: `LSP_COMPATIBILITY.md entry for ${server.name}`,
    status: 'pass',
  };
}

function checkDiagnosticReadiness(server: LSPServerRegistryEntry): CompatibilityCheck {
  const readiness = getLspDiagnosticReadiness(server.id);

  if (!readiness) {
    return {
      name: 'diagnostic-readiness',
      description: `Diagnostic readiness for ${server.name}`,
      status: 'warn',
      detail: 'No diagnostic readiness metadata in registry',
    };
  }

  return {
    name: 'diagnostic-readiness',
    description: `Diagnostic readiness for ${server.name}`,
    status: readiness.diagnosticEngine === 'v1' ? 'pass' : 'fail',
    detail: `Engine: ${readiness.diagnosticEngine}, Rich: ${readiness.richDiagnostics}, Closed-loop: ${readiness.closedLoop}`,
  };
}

function checkRuleManifest(
  server: LSPServerRegistryEntry,
  manifestDir?: string
): CompatibilityCheck {
  if (!manifestDir) {
    return {
      name: 'rule-manifest',
      description: `Rule manifest for ${server.name}`,
      status: 'skip',
      detail: 'No manifest directory provided',
    };
  }

  const manifestPath = path.join(manifestDir, `${server.id}-rules.json`);
  const result = tryReadManifest(manifestPath);

  if (result.error) {
    return {
      name: 'rule-manifest',
      description: `Rule manifest for ${server.name}`,
      status: 'warn',
      detail: result.error,
    };
  }

  if (!result.manifest) {
    return {
      name: 'rule-manifest',
      description: `Rule manifest for ${server.name}`,
      status: 'fail',
      detail: 'Manifest was not loaded',
    };
  }
  const ruleCount = result.manifest.rules.length;
  return {
    name: 'rule-manifest',
    description: `Rule manifest for ${server.name}`,
    status: 'pass',
    detail: `${ruleCount} rules loaded from manifest`,
  };
}

function checkLocalLaunchConfig(server: LSPServerRegistryEntry): CompatibilityCheck {
  if (!server.localLaunch) {
    return {
      name: 'local-launch-config',
      description: `Local launch configuration for ${server.name}`,
      status: 'fail',
      detail: 'No localLaunch defined in registry entry',
    };
  }

  if (!server.localLaunch.repoName) {
    return {
      name: 'local-launch-config',
      description: `Local launch configuration for ${server.name}`,
      status: 'fail',
      detail: 'localLaunch missing repoName',
    };
  }

  return {
    name: 'local-launch-config',
    description: `Local launch configuration for ${server.name}`,
    status: 'pass',
    detail: `kind: ${server.localLaunch.kind}, repo: ${server.localLaunch.repoName}`,
  };
}

/**
 * Check for lsp-capabilities.json manifest in sibling repo.
 *
 * Verifies the manifest exists, has correct languageId, and reports
 * capability gaps vs. the registry's known feature set.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/187
 */
function checkCapabilityManifest(
  server: LSPServerRegistryEntry,
  repoRoot: string
): CompatibilityCheck {
  const codeRoot = path.resolve(repoRoot, '..');
  const repoName = server.repository.split('/')[1];

  // Search multiple checkout locations
  const searchPaths = [
    path.join(codeRoot, repoName, 'lsp-capabilities.json'),
    path.join(codeRoot, '.worktrees-lsp-latest', repoName, 'lsp-capabilities.json'),
  ];

  let manifestPath: string | null = null;
  for (const p of searchPaths) {
    if (fs.existsSync(p)) {
      manifestPath = p;
      break;
    }
  }

  if (!manifestPath) {
    return {
      name: 'capability-manifest',
      description: `Capability manifest for ${server.name}`,
      status: 'warn',
      detail: 'lsp-capabilities.json not found in sibling checkout',
    };
  }

  try {
    const content = fs.readFileSync(manifestPath, 'utf8');
    const manifest = JSON.parse(content);

    // Validate languageId match
    if (manifest.languageId && manifest.languageId !== server.languageId) {
      return {
        name: 'capability-manifest',
        description: `Capability manifest for ${server.name}`,
        status: 'fail',
        detail: `languageId mismatch: manifest="${manifest.languageId}" registry="${server.languageId}"`,
      };
    }

    // Report capability coverage
    const caps = manifest.capabilities || {};
    const capabilityNames = [
      'diagnostics',
      'completion',
      'hover',
      'documentSymbol',
      'semanticTokens',
      'formatting',
      'codeAction',
      'references',
      'rename',
      'definition',
      'inlayHint',
    ];
    const supported = capabilityNames.filter(k => caps[k] === true);
    const missing = capabilityNames.filter(k => caps[k] === false);
    const unknown = capabilityNames.filter(k => caps[k] === undefined);

    return {
      name: 'capability-manifest',
      description: `Capability manifest for ${server.name}`,
      status: 'pass',
      detail: `${supported.length}/${capabilityNames.length} capabilities declared, ${missing.length} explicitly unsupported, ${unknown.length} unknown`,
    };
  } catch (e) {
    return {
      name: 'capability-manifest',
      description: `Capability manifest for ${server.name}`,
      status: 'fail',
      detail: `Failed to parse manifest: ${(e as Error).message}`,
    };
  }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function loadJsonSafe(filePath: string): Record<string, unknown> | null {
  try {
    return JSON.parse(fs.readFileSync(filePath, 'utf8')) as Record<string, unknown>;
  } catch {
    return null;
  }
}

function loadTextSafe(filePath: string): string | null {
  try {
    return fs.readFileSync(filePath, 'utf8');
  } catch {
    return null;
  }
}

function getPackageJsonLanguages(pkg: Record<string, unknown>): Record<string, unknown>[] {
  const contributes = pkg.contributes as Record<string, unknown> | undefined;
  if (!contributes) {
    return [];
  }
  return (contributes.languages as Record<string, unknown>[]) ?? [];
}

function getPackageJsonProperties(pkg: Record<string, unknown>): Record<string, unknown> {
  const contributes = pkg.contributes as Record<string, unknown> | undefined;
  if (!contributes) {
    return {};
  }
  const configuration = contributes.configuration as Record<string, unknown> | undefined;
  if (!configuration) {
    return {};
  }
  return (configuration.properties as Record<string, unknown>) ?? {};
}

// ---------------------------------------------------------------------------
// Report formatting
// ---------------------------------------------------------------------------

/**
 * Format a compatibility report as a human-readable markdown string.
 *
 * @param report - The compatibility report to format.
 * @returns Markdown-formatted report string.
 */
export function formatReportAsMarkdown(report: CompatibilityReport): string {
  const lines: string[] = [
    '# OpenQC LSP Compatibility Report',
    '',
    `Generated: ${report.generatedAt}`,
    '',
    `**Total:** ${report.totalServers} | **Passed:** ${report.passedServers} | **Failed:** ${report.failedServers}`,
    '',
    '## Per-Server Results',
    '',
  ];

  for (const entry of report.entries) {
    const icon = entry.passed ? '[PASS]' : '[FAIL]';
    lines.push(`### ${icon} ${entry.serverName} (${entry.serverId})`);
    lines.push('');
    lines.push(`- Language: \`${entry.languageId}\``);
    lines.push(`- Stability: ${entry.stability}`);

    for (const check of entry.checks) {
      const checkIcon =
        check.status === 'pass'
          ? '[PASS]'
          : check.status === 'fail'
            ? '[FAIL]'
            : check.status === 'warn'
              ? '[WARN]'
              : '[SKIP]';
      lines.push(`- ${checkIcon} ${check.name}: ${check.detail ?? 'OK'}`);
    }

    lines.push('');
  }

  return lines.join('\n');
}
