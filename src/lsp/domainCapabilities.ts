/**
 * Domain Language Description Capability Detection and API Consumer
 *
 * Discovers and consumes domain-specific language description APIs from
 * standalone LSP servers. Uses LSP capability negotiation (ServerCapabilities
 * extensions) to detect support and sends custom LSP requests when available.
 *
 * Integration is fully optional: LSPs that do not advertise these capabilities
 * are tracked as unsupported and all consumer methods return `undefined`
 * gracefully.
 *
 * @module lsp/domainCapabilities
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/145
 */

import type { LanguageClient } from 'vscode-languageclient/node';
import type {
  DomainLSPCapabilities,
  LanguageDescriptionParams,
  LanguageDescriptionResponse,
  MinimalExampleParams,
  MinimalExampleResponse,
  NextTokenGuidanceParams,
  NextTokenGuidanceResponse,
  OpenQCServerCapabilities,
  SectionKeywordParams,
  SectionKeywordResponse,
} from './types';
import { createComponentLogger } from '../utils/Logger';

const logger = createComponentLogger('DomainCapabilities');

// ---------------------------------------------------------------------------
// Capability defaults
// ---------------------------------------------------------------------------

/** All domain capabilities disabled. */
const NO_CAPABILITIES: Readonly<DomainLSPCapabilities> = {
  languageDescription: false,
  sectionKeywordLookup: false,
  minimalExamples: false,
  nextTokenGuidance: false,
} as const;

/** All domain capabilities enabled. */
const ALL_CAPABILITIES: Readonly<DomainLSPCapabilities> = {
  languageDescription: true,
  sectionKeywordLookup: true,
  minimalExamples: true,
  nextTokenGuidance: true,
} as const;

// ---------------------------------------------------------------------------
// Capability store
// ---------------------------------------------------------------------------

/**
 * Tracks which domain capabilities each running LSP server supports.
 *
 * Keyed by the client identity string used by `LSPManager`.
 */
const capabilityStore = new Map<string, DomainLSPCapabilities>();

// ---------------------------------------------------------------------------
// Capability detection
// ---------------------------------------------------------------------------

/**
 * Extract domain capabilities from raw server capabilities returned during
 * the LSP `initialize` handshake.
 *
 * @param rawCapabilities - The `capabilities` field from `InitializeResult`.
 * @returns A normalized `DomainLSPCapabilities` object.
 */
export function detectCapabilities(rawCapabilities: unknown): DomainLSPCapabilities {
  if (!rawCapabilities || typeof rawCapabilities !== 'object') {
    return { ...NO_CAPABILITIES };
  }

  const caps = rawCapabilities as OpenQCServerCapabilities;
  const domain = caps?.openqc?.domainCapabilities;

  if (!domain) {
    return { ...NO_CAPABILITIES };
  }

  return {
    languageDescription: domain.languageDescription === true,
    sectionKeywordLookup: domain.sectionKeywordLookup === true,
    minimalExamples: domain.minimalExamples === true,
    nextTokenGuidance: domain.nextTokenGuidance === true,
  };
}

/**
 * Record the detected capabilities for a running LSP client.
 *
 * @param clientKey - Unique key identifying the LSP client.
 * @param capabilities - Detected domain capabilities.
 */
export function recordCapabilities(clientKey: string, capabilities: DomainLSPCapabilities): void {
  capabilityStore.set(clientKey, capabilities);
  const flagList = Object.entries(capabilities)
    .filter(([, v]) => v)
    .map(([k]) => k);
  if (flagList.length > 0) {
    logger.info(`Client "${clientKey}" supports domain capabilities: ${flagList.join(', ')}`);
  } else {
    logger.info(`Client "${clientKey}" does not advertise domain capabilities`);
  }
}

/**
 * Remove capability tracking for a stopped LSP client.
 *
 * @param clientKey - Unique key identifying the LSP client.
 */
export function removeCapabilities(clientKey: string): void {
  capabilityStore.delete(clientKey);
}

/**
 * Retrieve the recorded capabilities for a running LSP client.
 *
 * @param clientKey - Unique key identifying the LSP client.
 * @returns The recorded capabilities, or `NO_CAPABILITIES` if unknown.
 */
export function getCapabilities(clientKey: string): DomainLSPCapabilities {
  return capabilityStore.get(clientKey) ?? { ...NO_CAPABILITIES };
}

/**
 * Check whether a specific domain capability is supported by a client.
 *
 * @param clientKey - Unique key identifying the LSP client.
 * @param capability - Name of the capability to check.
 * @returns `true` if the capability is advertised and the client is known.
 */
export function hasCapability(clientKey: string, capability: keyof DomainLSPCapabilities): boolean {
  const caps = capabilityStore.get(clientKey);
  return caps?.[capability] === true;
}

/**
 * Return all client keys that support at least one domain capability.
 */
export function clientsWithDomainSupport(): string[] {
  return Array.from(capabilityStore.entries())
    .filter(
      ([, caps]) =>
        caps.languageDescription ||
        caps.sectionKeywordLookup ||
        caps.minimalExamples ||
        caps.nextTokenGuidance
    )
    .map(([key]) => key);
}

// ---------------------------------------------------------------------------
// Client lifecycle integration
// ---------------------------------------------------------------------------

/**
 * Inspect a running `LanguageClient` and record its domain capabilities.
 *
 * Safe to call even when the client has not completed initialization --
 * returns `NO_CAPABILITIES` in that case.
 *
 * @param clientKey - Unique key identifying the LSP client.
 * @param client - The running `LanguageClient` instance.
 */
export function inspectAndRecordClient(
  clientKey: string,
  client: LanguageClient
): DomainLSPCapabilities {
  let rawCaps: unknown;
  try {
    rawCaps = client.initializeResult?.capabilities;
  } catch {
    rawCaps = undefined;
  }

  const capabilities = detectCapabilities(rawCaps);
  recordCapabilities(clientKey, capabilities);
  return capabilities;
}

// ---------------------------------------------------------------------------
// API consumers
// ---------------------------------------------------------------------------

/**
 * Request a language description from the LSP server.
 *
 * @param client - The running `LanguageClient` to send the request through.
 * @param params - Request parameters.
 * @returns The language description, or `undefined` if unavailable.
 */
export async function requestLanguageDescription(
  client: LanguageClient,
  params: LanguageDescriptionParams
): Promise<LanguageDescriptionResponse | undefined> {
  return sendDomainRequest<LanguageDescriptionResponse>(
    client,
    'openqc/languageDescription',
    params
  );
}

/**
 * Request section/keyword lookup from the LSP server.
 *
 * @param client - The running `LanguageClient` to send the request through.
 * @param params - Request parameters.
 * @returns Matching entries, or `undefined` if unavailable.
 */
export async function requestSectionKeywordLookup(
  client: LanguageClient,
  params: SectionKeywordParams
): Promise<SectionKeywordResponse | undefined> {
  return sendDomainRequest<SectionKeywordResponse>(client, 'openqc/sectionKeywordLookup', params);
}

/**
 * Request a minimal usage example from the LSP server.
 *
 * @param client - The running `LanguageClient` to send the request through.
 * @param params - Request parameters.
 * @returns The example, or `undefined` if unavailable.
 */
export async function requestMinimalExample(
  client: LanguageClient,
  params: MinimalExampleParams
): Promise<MinimalExampleResponse | undefined> {
  return sendDomainRequest<MinimalExampleResponse>(client, 'openqc/minimalExample', params);
}

/**
 * Request next-token guidance from the LSP server.
 *
 * @param client - The running `LanguageClient` to send the request through.
 * @param params - Request parameters.
 * @returns Candidate tokens, or `undefined` if unavailable.
 */
export async function requestNextTokenGuidance(
  client: LanguageClient,
  params: NextTokenGuidanceParams
): Promise<NextTokenGuidanceResponse | undefined> {
  return sendDomainRequest<NextTokenGuidanceResponse>(client, 'openqc/nextTokenGuidance', params);
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/**
 * Send a custom LSP request, catching errors from servers that have not
 * implemented the endpoint.
 */
async function sendDomainRequest<T>(
  client: LanguageClient,
  method: string,
  params: unknown
): Promise<T | undefined> {
  try {
    const result = await client.sendRequest(method, params);
    return result as T;
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    // MethodNotFound or similar LSP errors are expected for unsupported servers
    logger.info(`Domain request "${method}" not handled: ${message}`);
    return undefined;
  }
}

// ---------------------------------------------------------------------------
// Test-only reset
// ---------------------------------------------------------------------------

/**
 * Clear all recorded capabilities. Exported for test teardown only.
 */
export function resetCapabilityStore(): void {
  capabilityStore.clear();
}

/**
 * Return a frozen copy of `ALL_CAPABILITIES` for use in tests.
 */
export function allCapabilitiesSupported(): Readonly<DomainLSPCapabilities> {
  return ALL_CAPABILITIES;
}

/**
 * Return a frozen copy of `NO_CAPABILITIES` for use in tests.
 */
export function noCapabilitiesSupported(): Readonly<DomainLSPCapabilities> {
  return NO_CAPABILITIES;
}
