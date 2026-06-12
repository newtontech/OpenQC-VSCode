/**
 * Smoke Test Infrastructure
 *
 * Public API surface for the OpenQC smoke test and verification system.
 *
 * @module smoke
 */

export type {
  RuleSeverity,
  RuleManifestEntry,
  RuleManifest,
  CheckStatus,
  CompatibilityCheck,
  LspCompatibilityEntry,
  CompatibilityReport,
  DocumentKind,
  DocumentDetectionResult,
  SmokeTestResult,
  SmokeTestSummary,
  VsixVerificationResult,
  RepoGitHubStatus,
  RepoGateResult,
  RepoCleanStatus,
  CampaignReport,
} from './types';

export {
  validateManifest,
  readManifestFromFile,
  parseManifestString,
  tryReadManifest,
  getRulesBySeverity,
  getRuleByCode,
  getManifestCategories,
  countRulesBySeverity,
} from './ruleManifestReader';

export {
  detectDocument,
  isOutputOrLogFile,
  getOutputPatternsForLanguage,
  getLogPatternsForLanguage,
} from './documentDetector';

export { generateCompatibilityReport, formatReportAsMarkdown } from './compatibilityReport';

export {
  runSmokeTestForLsp,
  runAllSmokeTests,
  verifyVsixPackage,
  checkGitHubStatus,
  runLocalGates,
  checkRepoCleanliness,
  generateCampaignReport,
  formatCampaignReportAsMarkdown,
} from './smokeRunner';
