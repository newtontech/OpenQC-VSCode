/**
 * Tests for domain language description capability detection and API consumers.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/145
 */

import type { LanguageClient } from 'vscode-languageclient/node';
import {
  allCapabilitiesSupported,
  clientsWithDomainSupport,
  detectCapabilities,
  getCapabilities,
  hasCapability,
  inspectAndRecordClient,
  noCapabilitiesSupported,
  recordCapabilities,
  removeCapabilities,
  requestLanguageDescription,
  requestMinimalExample,
  requestNextTokenGuidance,
  requestSectionKeywordLookup,
  resetCapabilityStore,
} from '../../../src/lsp/domainCapabilities';
import type {
  DomainLSPCapabilities,
  LanguageDescriptionResponse,
  MinimalExampleResponse,
  NextTokenGuidanceResponse,
  SectionKeywordResponse,
} from '../../../src/lsp/types';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Build a mock LanguageClient with configurable initializeResult. */
function mockClient(
  capabilities?: unknown,
  requestHandler?: (method: string, params: unknown) => unknown
): LanguageClient {
  return {
    initializeResult: capabilities ? { capabilities } : undefined,
    sendRequest: jest.fn(async (method: string, params: unknown) => {
      if (requestHandler) {
        return requestHandler(method, params);
      }
      throw new Error(`Method not found: ${method}`);
    }),
  } as unknown as LanguageClient;
}

/** Full domain capabilities as returned by a fully-supporting server. */
function fullServerCapabilities(): object {
  return {
    completionProvider: {},
    hoverProvider: true,
    openqc: {
      domainCapabilities: {
        languageDescription: true,
        sectionKeywordLookup: true,
        minimalExamples: true,
        nextTokenGuidance: true,
      },
    },
  };
}

/** Partial domain capabilities (only languageDescription and minimalExamples). */
function partialServerCapabilities(): object {
  return {
    completionProvider: {},
    openqc: {
      domainCapabilities: {
        languageDescription: true,
        minimalExamples: true,
      },
    },
  };
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

describe('domainCapabilities', () => {
  beforeEach(() => {
    resetCapabilityStore();
  });

  // -------------------------------------------------------------------------
  // detectCapabilities
  // -------------------------------------------------------------------------

  describe('detectCapabilities', () => {
    it('detects full domain capabilities from server capabilities', () => {
      const caps = detectCapabilities(fullServerCapabilities());
      expect(caps).toEqual({
        languageDescription: true,
        sectionKeywordLookup: true,
        minimalExamples: true,
        nextTokenGuidance: true,
      });
    });

    it('detects partial domain capabilities from server capabilities', () => {
      const caps = detectCapabilities(partialServerCapabilities());
      expect(caps).toEqual({
        languageDescription: true,
        sectionKeywordLookup: false,
        minimalExamples: true,
        nextTokenGuidance: false,
      });
    });

    it('returns no capabilities for a server without custom extensions', () => {
      const caps = detectCapabilities({
        completionProvider: {},
        hoverProvider: true,
      });
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('returns no capabilities for null input', () => {
      const caps = detectCapabilities(null);
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('returns no capabilities for undefined input', () => {
      const caps = detectCapabilities(undefined);
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('returns no capabilities for non-object input', () => {
      const caps = detectCapabilities('not an object');
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('treats non-boolean capability values as false', () => {
      const caps = detectCapabilities({
        openqc: {
          domainCapabilities: {
            languageDescription: 'yes',
            sectionKeywordLookup: 1,
            minimalExamples: null,
            nextTokenGuidance: undefined,
          },
        },
      });
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('handles empty openqc object', () => {
      const caps = detectCapabilities({ openqc: {} });
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('handles openqc with missing domainCapabilities', () => {
      const caps = detectCapabilities({ openqc: { otherFeature: true } });
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });
  });

  // -------------------------------------------------------------------------
  // recordCapabilities / getCapabilities / removeCapabilities
  // -------------------------------------------------------------------------

  describe('capability store', () => {
    it('records and retrieves capabilities', () => {
      const caps: DomainLSPCapabilities = {
        languageDescription: true,
        sectionKeywordLookup: false,
        minimalExamples: true,
        nextTokenGuidance: false,
      };
      recordCapabilities('test-client', caps);
      expect(getCapabilities('test-client')).toEqual(caps);
    });

    it('returns no capabilities for unknown client keys', () => {
      expect(getCapabilities('unknown')).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('overwrites capabilities on re-recording', () => {
      recordCapabilities('client-a', allCapabilitiesSupported());
      recordCapabilities('client-a', noCapabilitiesSupported());
      expect(getCapabilities('client-a')).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('removes capabilities for a client', () => {
      recordCapabilities('client-b', allCapabilitiesSupported());
      removeCapabilities('client-b');
      expect(getCapabilities('client-b')).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('does not throw when removing a non-existent client', () => {
      expect(() => removeCapabilities('nonexistent')).not.toThrow();
    });
  });

  // -------------------------------------------------------------------------
  // hasCapability
  // -------------------------------------------------------------------------

  describe('hasCapability', () => {
    it('returns true for an advertised capability', () => {
      recordCapabilities('client-c', allCapabilitiesSupported());
      expect(hasCapability('client-c', 'languageDescription')).toBe(true);
      expect(hasCapability('client-c', 'sectionKeywordLookup')).toBe(true);
      expect(hasCapability('client-c', 'minimalExamples')).toBe(true);
      expect(hasCapability('client-c', 'nextTokenGuidance')).toBe(true);
    });

    it('returns false for a capability not advertised', () => {
      recordCapabilities('client-d', noCapabilitiesSupported());
      expect(hasCapability('client-d', 'languageDescription')).toBe(false);
    });

    it('returns false for unknown client', () => {
      expect(hasCapability('unknown', 'languageDescription')).toBe(false);
    });
  });

  // -------------------------------------------------------------------------
  // clientsWithDomainSupport
  // -------------------------------------------------------------------------

  describe('clientsWithDomainSupport', () => {
    it('returns only clients that support at least one domain capability', () => {
      recordCapabilities('with-support', allCapabilitiesSupported());
      recordCapabilities('no-support', noCapabilitiesSupported());
      recordCapabilities('partial', {
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: true,
        nextTokenGuidance: false,
      });

      const supported = clientsWithDomainSupport();
      expect(supported).toContain('with-support');
      expect(supported).toContain('partial');
      expect(supported).not.toContain('no-support');
    });

    it('returns empty array when no clients have domain support', () => {
      recordCapabilities('client-e', noCapabilitiesSupported());
      expect(clientsWithDomainSupport()).toEqual([]);
    });
  });

  // -------------------------------------------------------------------------
  // inspectAndRecordClient
  // -------------------------------------------------------------------------

  describe('inspectAndRecordClient', () => {
    it('inspects a fully-supporting client and records capabilities', () => {
      const client = mockClient(fullServerCapabilities());
      const caps = inspectAndRecordClient('inspect-full', client);

      expect(caps).toEqual({
        languageDescription: true,
        sectionKeywordLookup: true,
        minimalExamples: true,
        nextTokenGuidance: true,
      });
      expect(getCapabilities('inspect-full')).toEqual(caps);
    });

    it('inspects a partially-supporting client and records capabilities', () => {
      const client = mockClient(partialServerCapabilities());
      const caps = inspectAndRecordClient('inspect-partial', client);

      expect(caps).toEqual({
        languageDescription: true,
        sectionKeywordLookup: false,
        minimalExamples: true,
        nextTokenGuidance: false,
      });
    });

    it('inspects an unsupported client (no openqc extensions)', () => {
      const client = mockClient({ completionProvider: {} });
      const caps = inspectAndRecordClient('inspect-unsupported', client);

      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('handles client with no initializeResult', () => {
      const client = mockClient(undefined);
      const caps = inspectAndRecordClient('inspect-no-init', client);

      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });

    it('handles client where initializeResult throws', () => {
      const client = {
        get initializeResult() {
          throw new Error('not initialized');
        },
        sendRequest: jest.fn(),
      } as unknown as LanguageClient;

      const caps = inspectAndRecordClient('inspect-throw', client);
      expect(caps).toEqual({
        languageDescription: false,
        sectionKeywordLookup: false,
        minimalExamples: false,
        nextTokenGuidance: false,
      });
    });
  });

  // -------------------------------------------------------------------------
  // API consumers
  // -------------------------------------------------------------------------

  describe('requestLanguageDescription', () => {
    it('returns a language description from a supporting server', async () => {
      const mockResponse: LanguageDescriptionResponse = {
        name: 'Gaussian',
        description: 'Quantum chemistry input format',
        sections: ['route', 'title', 'charge_and_multiplicity', 'geometry'],
      };
      const client = mockClient(fullServerCapabilities(), (_method, _params) => mockResponse);

      const result = await requestLanguageDescription(client, { languageId: 'gaussian' });
      expect(result).toEqual(mockResponse);
    });

    it('returns undefined when the server does not implement the endpoint', async () => {
      const client = mockClient(fullServerCapabilities());
      const result = await requestLanguageDescription(client, { languageId: 'gaussian' });
      expect(result).toBeUndefined();
    });
  });

  describe('requestSectionKeywordLookup', () => {
    it('returns section/keyword entries from a supporting server', async () => {
      const mockResponse: SectionKeywordResponse = {
        entries: [
          { name: '#P', documentation: 'Normal processing output', section: 'route' },
          { name: '#N', documentation: 'No output', section: 'route' },
        ],
      };
      const client = mockClient(fullServerCapabilities(), (_method, _params) => mockResponse);

      const result = await requestSectionKeywordLookup(client, {
        languageId: 'gaussian',
        query: '#',
      });
      expect(result).toEqual(mockResponse);
    });

    it('returns undefined when the server does not implement the endpoint', async () => {
      const client = mockClient(fullServerCapabilities());
      const result = await requestSectionKeywordLookup(client, {
        languageId: 'gaussian',
        query: '#',
      });
      expect(result).toBeUndefined();
    });
  });

  describe('requestMinimalExample', () => {
    it('returns a minimal example from a supporting server', async () => {
      const mockResponse: MinimalExampleResponse = {
        example: '#P B3LYP/6-31G(d) opt\n\nTitle\n\n0 1\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n',
        description: 'Geometry optimization of H2 at B3LYP/6-31G(d)',
      };
      const client = mockClient(fullServerCapabilities(), (_method, _params) => mockResponse);

      const result = await requestMinimalExample(client, { languageId: 'gaussian' });
      expect(result).toEqual(mockResponse);
    });

    it('returns undefined when the server does not implement the endpoint', async () => {
      const client = mockClient(fullServerCapabilities());
      const result = await requestMinimalExample(client, { languageId: 'gaussian' });
      expect(result).toBeUndefined();
    });
  });

  describe('requestNextTokenGuidance', () => {
    it('returns next-token candidates from a supporting server', async () => {
      const mockResponse: NextTokenGuidanceResponse = {
        candidates: [
          { token: 'B3LYP', documentation: 'Becke three-parameter hybrid functional' },
          { token: 'HF', documentation: 'Hartree-Fock' },
        ],
      };
      const client = mockClient(fullServerCapabilities(), (_method, _params) => mockResponse);

      const result = await requestNextTokenGuidance(client, {
        languageId: 'gaussian',
        textDocument: { uri: 'file:///test.gjf' },
        position: { line: 0, character: 3 },
      });
      expect(result).toEqual(mockResponse);
    });

    it('returns undefined when the server does not implement the endpoint', async () => {
      const client = mockClient(fullServerCapabilities());
      const result = await requestNextTokenGuidance(client, {
        languageId: 'gaussian',
        textDocument: { uri: 'file:///test.gjf' },
        position: { line: 0, character: 3 },
      });
      expect(result).toBeUndefined();
    });
  });

  // -------------------------------------------------------------------------
  // Graceful degradation
  // -------------------------------------------------------------------------

  describe('graceful degradation', () => {
    it('does not throw when a request is sent to a server that drops the connection', async () => {
      const client = {
        initializeResult: { capabilities: fullServerCapabilities() },
        sendRequest: jest.fn().mockRejectedValue(new Error('Connection dropped')),
      } as unknown as LanguageClient;

      await expect(
        requestLanguageDescription(client, { languageId: 'gaussian' })
      ).resolves.toBeUndefined();
    });

    it('all consumers return undefined for unsupported servers', async () => {
      const client = mockClient({ completionProvider: {} });

      await expect(
        requestLanguageDescription(client, { languageId: 'gaussian' })
      ).resolves.toBeUndefined();

      await expect(
        requestSectionKeywordLookup(client, { languageId: 'gaussian', query: '#' })
      ).resolves.toBeUndefined();

      await expect(
        requestMinimalExample(client, { languageId: 'gaussian' })
      ).resolves.toBeUndefined();

      await expect(
        requestNextTokenGuidance(client, {
          languageId: 'gaussian',
          textDocument: { uri: 'file:///test.gjf' },
          position: { line: 0, character: 0 },
        })
      ).resolves.toBeUndefined();
    });
  });
});
