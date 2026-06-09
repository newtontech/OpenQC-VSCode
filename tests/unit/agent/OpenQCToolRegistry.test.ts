/**
 * Tests for the OpenQC tool registry and MCP tools (issue #80).
 * @module tests/unit/agent/OpenQCToolRegistry.test
 */

import {
  OPENQC_TOOLS,
  getToolByName,
  getToolsByCategory,
  getToolList,
  requiresConfirmation,
  isReadOnly,
} from '../../../src/agent/OpenQCToolRegistry';

describe('OpenQCToolRegistry', () => {
  describe('tool definitions', () => {
    it('has at least 9 tools defined', () => {
      expect(OPENQC_TOOLS.length).toBeGreaterThanOrEqual(9);
    });

    it('all tools have required fields', () => {
      for (const tool of OPENQC_TOOLS) {
        expect(tool.name).toBeTruthy();
        expect(tool.description).toBeTruthy();
        expect(tool.category).toBeDefined();
        expect(tool.inputSchema).toBeDefined();
        expect(tool.handler).toBeTruthy();
      }
    });

    it('tools have valid categories', () => {
      const validCategories = ['read-only', 'preview-only', 'write', 'external'];
      for (const tool of OPENQC_TOOLS) {
        expect(validCategories).toContain(tool.category);
      }
    });
  });

  describe('getToolByName', () => {
    it('finds parseStructure tool', () => {
      const tool = getToolByName('openqc.parseStructure');
      expect(tool).toBeDefined();
      expect(tool!.category).toBe('read-only');
    });

    it('finds parseCalculationOutput tool', () => {
      const tool = getToolByName('openqc.parseCalculationOutput');
      expect(tool).toBeDefined();
      expect(tool!.category).toBe('read-only');
    });

    it('finds checkPythonBackend tool', () => {
      const tool = getToolByName('openqc.checkPythonBackend');
      expect(tool).toBeDefined();
      expect(tool!.category).toBe('read-only');
    });

    it('finds generateSupercell tool', () => {
      const tool = getToolByName('openqc.generateSupercell');
      expect(tool).toBeDefined();
      expect(tool!.category).toBe('preview-only');
    });

    it('returns undefined for unknown tool', () => {
      expect(getToolByName('openqc.nonexistent')).toBeUndefined();
    });
  });

  describe('getToolsByCategory', () => {
    it('returns read-only tools', () => {
      const tools = getToolsByCategory('read-only');
      expect(tools.length).toBeGreaterThanOrEqual(3);
      for (const tool of tools) {
        expect(tool.category).toBe('read-only');
      }
    });

    it('returns preview-only tools', () => {
      const tools = getToolsByCategory('preview-only');
      expect(tools.length).toBeGreaterThanOrEqual(2);
    });

    it('returns write tools', () => {
      const tools = getToolsByCategory('write');
      expect(tools.length).toBeGreaterThanOrEqual(1);
    });
  });

  describe('safety classification', () => {
    it('read-only tools do not require confirmation', () => {
      expect(requiresConfirmation('openqc.parseStructure')).toBe(false);
      expect(requiresConfirmation('openqc.checkPythonBackend')).toBe(false);
    });

    it('write tools require confirmation', () => {
      expect(requiresConfirmation('openqc.convertStructure')).toBe(true);
      expect(requiresConfirmation('openqc.generateInput')).toBe(true);
    });

    it('unknown tools require confirmation', () => {
      expect(requiresConfirmation('openqc.unknown')).toBe(true);
    });

    it('isReadOnly returns true for read-only tools', () => {
      expect(isReadOnly('openqc.parseStructure')).toBe(true);
      expect(isReadOnly('openqc.generateSupercell')).toBe(false);
    });
  });

  describe('getToolList', () => {
    it('returns list with name, description, category', () => {
      const list = getToolList();
      expect(list.length).toBe(OPENQC_TOOLS.length);
      for (const item of list) {
        expect(item.name).toBeTruthy();
        expect(item.description).toBeTruthy();
        expect(item.category).toBeDefined();
      }
    });
  });

  describe('input schemas', () => {
    it('parseStructure has required path parameter', () => {
      const tool = getToolByName('openqc.parseStructure')!;
      expect(tool.inputSchema.properties).toHaveProperty('path');
    });

    it('generateSupercell has na, nb, nc parameters', () => {
      const tool = getToolByName('openqc.generateSupercell')!;
      expect(tool.inputSchema.properties).toHaveProperty('na');
      expect(tool.inputSchema.properties).toHaveProperty('nb');
      expect(tool.inputSchema.properties).toHaveProperty('nc');
    });
  });
});
