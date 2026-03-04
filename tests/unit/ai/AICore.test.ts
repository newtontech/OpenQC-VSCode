/**
 * AICore Unit Tests
 */

import { AICore, AICoreFactory, AIProvider, AIRequestType } from '../../../src/ai/AICore';

// Mock vscode
const mockGetConfiguration = jest.fn();
const mockContext = {
  extensionPath: '/test/extension',
} as any;

jest.mock('vscode', () => ({
  workspace: {
    getConfiguration: () => mockGetConfiguration(),
  },
  window: {
    showErrorMessage: jest.fn(),
    showInformationMessage: jest.fn(),
  },
}));

// Mock child_process
const mockSpawn = jest.fn();
jest.mock('child_process', () => ({
  spawn: (...args: any[]) => mockSpawn(...args),
}));

describe('AICore', () => {
  let aiCore: AICore;

  beforeEach(() => {
    jest.clearAllMocks();
    AICoreFactory.reset();
    
    // Default mock configuration
    mockGetConfiguration.mockReturnValue({
      get: jest.fn((key: string, defaultValue?: any) => {
        const config: Record<string, any> = {
          'enabled': true,
          'provider': AIProvider.Ollama,
          'apiKey': undefined,
          'model': 'llama2',
          'ollamaUrl': 'http://localhost:11434',
          'maxTokens': 2048,
          'temperature': 0.7,
          'pythonPath': 'python3',
        };
        return config[key] ?? defaultValue;
      }),
    });
    
    aiCore = new AICore(mockContext);
  });

  describe('Configuration', () => {
    it('should load configuration correctly', () => {
      const config = aiCore.getConfig();
      
      expect(config.enabled).toBe(true);
      expect(config.provider).toBe(AIProvider.Ollama);
      expect(config.model).toBe('llama2');
      expect(config.ollamaUrl).toBe('http://localhost:11434');
      expect(config.maxTokens).toBe(2048);
      expect(config.temperature).toBe(0.7);
    });

    it('should refresh configuration', () => {
      mockGetConfiguration.mockReturnValue({
        get: jest.fn((key: string, defaultValue?: any) => {
          const config: Record<string, any> = {
            'enabled': false,
            'provider': AIProvider.OpenAI,
            'apiKey': 'test-key',
            'model': 'gpt-4',
            'ollamaUrl': 'http://localhost:11434',
            'maxTokens': 4096,
            'temperature': 0.5,
            'pythonPath': 'python3',
          };
          return config[key] ?? defaultValue;
        }),
      });
      
      aiCore.refreshConfig();
      const config = aiCore.getConfig();
      
      expect(config.enabled).toBe(false);
      expect(config.provider).toBe(AIProvider.OpenAI);
      expect(config.apiKey).toBe('test-key');
      expect(config.model).toBe('gpt-4');
    });

    it('should check if AI is enabled', () => {
      expect(aiCore.isEnabled()).toBe(true);
      
      mockGetConfiguration.mockReturnValue({
        get: jest.fn(() => false),
      });
      aiCore.refreshConfig();
      
      expect(aiCore.isEnabled()).toBe(false);
    });
  });

  describe('Validation', () => {
    it('should validate configuration when enabled', () => {
      const validation = aiCore.validateConfig();
      
      expect(validation.valid).toBe(true);
      expect(validation.errors).toHaveLength(0);
    });

    it('should fail validation when disabled', () => {
      mockGetConfiguration.mockReturnValue({
        get: jest.fn((key: string) => {
          if (key === 'enabled') return false;
          return undefined;
        }),
      });
      aiCore.refreshConfig();
      
      const validation = aiCore.validateConfig();
      
      expect(validation.valid).toBe(false);
      expect(validation.errors).toContain('AI features are disabled');
    });

    it('should require API key for OpenAI', () => {
      mockGetConfiguration.mockReturnValue({
        get: jest.fn((key: string, defaultValue?: any) => {
          const config: Record<string, any> = {
            'enabled': true,
            'provider': AIProvider.OpenAI,
            'apiKey': undefined,
            'model': 'gpt-4',
            'ollamaUrl': 'http://localhost:11434',
            'maxTokens': 2048,
            'temperature': 0.7,
          };
          return config[key] ?? defaultValue;
        }),
      });
      aiCore.refreshConfig();
      
      const validation = aiCore.validateConfig();
      
      expect(validation.valid).toBe(false);
      expect(validation.errors).toContain('OpenAI API key is required');
    });

    it('should validate temperature range', () => {
      mockGetConfiguration.mockReturnValue({
        get: jest.fn((key: string, defaultValue?: any) => {
          const config: Record<string, any> = {
            'enabled': true,
            'provider': AIProvider.Ollama,
            'model': 'llama2',
            'temperature': 3.0, // Invalid: > 2
          };
          return config[key] ?? defaultValue;
        }),
      });
      aiCore.refreshConfig();
      
      const validation = aiCore.validateConfig();
      
      expect(validation.valid).toBe(false);
      expect(validation.errors).toContain('Temperature must be between 0 and 2');
    });
  });

  describe('AI Operations', () => {
    const mockSuccessResponse = {
      success: true,
      content: 'Test response',
      suggestions: [],
    };

    beforeEach(() => {
      // Mock successful spawn
      mockSpawn.mockReturnValue({
        stdin: { write: jest.fn(), end: jest.fn() },
        stdout: { on: jest.fn((event: string, callback: Function) => {
          if (event === 'data') {
            callback(JSON.stringify(mockSuccessResponse));
          }
        })},
        stderr: { on: jest.fn() },
        on: jest.fn((event: string, callback: Function) => {
          if (event === 'close') {
            callback(0);
          }
        }),
      });
    });

    it('should optimize input', async () => {
      const result = await aiCore.optimizeInput('test input', 'vasp');
      
      expect(result.success).toBe(true);
      expect(mockSpawn).toHaveBeenCalled();
    });

    it('should generate input', async () => {
      const result = await aiCore.generateInput('test description', 'cp2k');
      
      expect(result.success).toBe(true);
      expect(mockSpawn).toHaveBeenCalled();
    });

    it('should explain parameters', async () => {
      const result = await aiCore.explainParameters('test input', 'gaussian');
      
      expect(result.success).toBe(true);
      expect(mockSpawn).toHaveBeenCalled();
    });

    it('should debug calculation', async () => {
      const result = await aiCore.debugCalculation('input', 'output', 'qe');
      
      expect(result.success).toBe(true);
      expect(mockSpawn).toHaveBeenCalled();
    });

    it('should handle Python errors', async () => {
      mockSpawn.mockReturnValue({
        stdin: { write: jest.fn(), end: jest.fn() },
        stdout: { on: jest.fn() },
        stderr: { on: jest.fn((event: string, callback: Function) => {
          if (event === 'data') {
            callback('Python error');
          }
        })},
        on: jest.fn((event: string, callback: Function) => {
          if (event === 'close') {
            callback(1);
          }
        }),
      });
      
      const result = await aiCore.optimizeInput('test', 'vasp');
      
      expect(result.success).toBe(false);
      expect(result.error).toContain('AI backend failed');
    });

    it('should handle invalid JSON response', async () => {
      mockSpawn.mockReturnValue({
        stdin: { write: jest.fn(), end: jest.fn() },
        stdout: { on: jest.fn((event: string, callback: Function) => {
          if (event === 'data') {
            callback('invalid json');
          }
        })},
        stderr: { on: jest.fn() },
        on: jest.fn((event: string, callback: Function) => {
          if (event === 'close') {
            callback(0);
          }
        }),
      });
      
      const result = await aiCore.optimizeInput('test', 'vasp');
      
      expect(result.success).toBe(false);
      expect(result.error).toContain('Failed to parse AI response');
    });
  });

  describe('Availability Check', () => {
    it('should return false when disabled', async () => {
      mockGetConfiguration.mockReturnValue({
        get: jest.fn(() => false),
      });
      aiCore.refreshConfig();
      
      const available = await aiCore.isAvailable();
      expect(available).toBe(false);
    });

    it('should check availability via Python backend', async () => {
      mockSpawn.mockReturnValue({
        stdin: { write: jest.fn(), end: jest.fn() },
        stdout: { on: jest.fn((event: string, callback: Function) => {
          if (event === 'data') {
            callback(JSON.stringify({ success: true }));
          }
        })},
        stderr: { on: jest.fn() },
        on: jest.fn((event: string, callback: Function) => {
          if (event === 'close') {
            callback(0);
          }
        }),
      });
      
      const available = await aiCore.isAvailable();
      expect(available).toBe(true);
    });
  });

  describe('Available Models', () => {
    it('should return OpenAI models when provider is OpenAI', async () => {
      // Mock OpenAI provider
      mockGetConfiguration.mockReturnValue({
        get: jest.fn((key: string, defaultValue?: any) => {
          const config: Record<string, any> = {
            'enabled': true,
            'provider': AIProvider.OpenAI,
            'apiKey': 'test-key',
            'model': 'gpt-4',
            'ollamaUrl': 'http://localhost:11434',
            'maxTokens': 2048,
            'temperature': 0.7,
            'pythonPath': 'python3',
          };
          return config[key] ?? defaultValue;
        }),
      });
      aiCore.refreshConfig();
      
      const models = await aiCore.getAvailableModels();
      
      expect(models).toContain('gpt-4');
      expect(models).toContain('gpt-3.5-turbo');
    });

    it('should return Ollama models', async () => {
      global.fetch = jest.fn().mockResolvedValue({
        ok: true,
        json: () => Promise.resolve({ models: [{ name: 'llama2' }, { name: 'mistral' }] }),
      });
      
      const models = await aiCore.getAvailableModels();
      
      expect(models).toContain('llama2');
      expect(models).toContain('mistral');
    });

    it('should handle Ollama fetch error', async () => {
      global.fetch = jest.fn().mockRejectedValue(new Error('Connection failed'));
      
      const models = await aiCore.getAvailableModels();
      
      expect(models).toContain('llama2');
      expect(models).toContain('codellama');
    });
  });

  describe('Factory', () => {
    it('should return singleton instance', () => {
      const instance1 = AICoreFactory.getInstance(mockContext);
      const instance2 = AICoreFactory.getInstance(mockContext);
      
      expect(instance1).toBe(instance2);
    });

    it('should reset instance', () => {
      const instance1 = AICoreFactory.getInstance(mockContext);
      AICoreFactory.reset();
      const instance2 = AICoreFactory.getInstance(mockContext);
      
      expect(instance1).not.toBe(instance2);
    });
  });
});
