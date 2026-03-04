/**
 * AI Core Module - TypeScript Interface
 *
 * Provides TypeScript interface to the Python AI backend for
 * AI-powered input file optimization, generation, and debugging.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { spawn } from 'child_process';

/**
 * Supported AI providers
 */
export enum AIProvider {
  OpenAI = 'openai',
  Ollama = 'ollama',
}

/**
 * AI Configuration interface
 */
export interface AIConfig {
  enabled: boolean;
  provider: AIProvider;
  apiKey?: string;
  model: string;
  ollamaUrl: string;
  maxTokens: number;
  temperature: number;
}

/**
 * AI Request types
 */
export enum AIRequestType {
  OptimizeInput = 'optimize_input',
  GenerateInput = 'generate_input',
  ExplainParameters = 'explain_parameters',
  DebugCalculation = 'debug_calculation',
}

/**
 * AI Request interface
 */
export interface AIRequest {
  type: AIRequestType;
  content: string;
  software?: string;
  context?: Record<string, any>;
}

/**
 * AI Response interface
 */
export interface AIResponse {
  success: boolean;
  content?: string;
  suggestions?: Array<{
    type: 'optimization' | 'warning' | 'info';
    message: string;
    parameter?: string;
    currentValue?: string;
    suggestedValue?: string;
    explanation?: string;
  }>;
  generatedInput?: string;
  error?: string;
  metadata?: Record<string, any>;
}

/**
 * AI Core class
 *
 * Interfaces with Python AI backend for LLM operations
 */
export class AICore {
  private config: AIConfig;
  private pythonPath: string;
  private aiScript: string;
  private context: vscode.ExtensionContext;

  constructor(context: vscode.ExtensionContext) {
    this.context = context;
    this.pythonPath = this.getPythonPath();
    this.aiScript = path.join(
      context.extensionPath,
      'python',
      'openqc',
      'ai',
      'client.py'
    );
    this.config = this.loadConfig();
  }

  /**
   * Get Python executable path from configuration
   */
  private getPythonPath(): string {
    const config = vscode.workspace.getConfiguration('openqc');
    return config.get<string>('pythonPath', 'python3');
  }

  /**
   * Load AI configuration from VSCode settings
   */
  private loadConfig(): AIConfig {
    const config = vscode.workspace.getConfiguration('openqc.ai');
    return {
      enabled: config.get<boolean>('enabled', false),
      provider: config.get<AIProvider>('provider', AIProvider.Ollama),
      apiKey: config.get<string>('apiKey'),
      model: config.get<string>('model', 'llama2'),
      ollamaUrl: config.get<string>('ollamaUrl', 'http://localhost:11434'),
      maxTokens: config.get<number>('maxTokens', 2048),
      temperature: config.get<number>('temperature', 0.7),
    };
  }

  /**
   * Refresh configuration (call when settings change)
   */
  public refreshConfig(): void {
    this.config = this.loadConfig();
  }

  /**
   * Check if AI features are enabled and available
   */
  public isEnabled(): boolean {
    return this.config.enabled;
  }

  /**
   * Check if AI backend is available
   */
  public async isAvailable(): Promise<boolean> {
    if (!this.config.enabled) {
      return false;
    }

    try {
      const result = await this.executeAI(['check'], '');
      return result.success;
    } catch {
      return false;
    }
  }

  /**
   * Get current configuration
   */
  public getConfig(): AIConfig {
    return { ...this.config };
  }

  /**
   * Optimize input file parameters
   */
  public async optimizeInput(
    inputContent: string,
    software: string
  ): Promise<AIResponse> {
    const request: AIRequest = {
      type: AIRequestType.OptimizeInput,
      content: inputContent,
      software,
    };

    return this.executeAI(['optimize'], JSON.stringify(request));
  }

  /**
   * Generate input file from natural language description
   */
  public async generateInput(
    description: string,
    software: string,
    context?: Record<string, any>
  ): Promise<AIResponse> {
    const request: AIRequest = {
      type: AIRequestType.GenerateInput,
      content: description,
      software,
      context,
    };

    return this.executeAI(['generate'], JSON.stringify(request));
  }

  /**
   * Explain input file parameters
   */
  public async explainParameters(
    inputContent: string,
    software: string
  ): Promise<AIResponse> {
    const request: AIRequest = {
      type: AIRequestType.ExplainParameters,
      content: inputContent,
      software,
    };

    return this.executeAI(['explain'], JSON.stringify(request));
  }

  /**
   * Debug failed calculation
   */
  public async debugCalculation(
    inputContent: string,
    outputContent: string,
    software: string
  ): Promise<AIResponse> {
    const request: AIRequest = {
      type: AIRequestType.DebugCalculation,
      content: inputContent,
      software,
      context: { output: outputContent },
    };

    return this.executeAI(['debug'], JSON.stringify(request));
  }

  /**
   * Execute AI Python script
   */
  private async executeAI(args: string[], input: string): Promise<AIResponse> {
    return new Promise((resolve, reject) => {
      const env = {
        ...process.env,
        OPENQC_AI_PROVIDER: this.config.provider,
        OPENQC_AI_MODEL: this.config.model,
        OPENQC_AI_API_KEY: this.config.apiKey || '',
        OPENQC_AI_OLLAMA_URL: this.config.ollamaUrl,
        OPENQC_AI_MAX_TOKENS: this.config.maxTokens.toString(),
        OPENQC_AI_TEMPERATURE: this.config.temperature.toString(),
      };

      const process_ = spawn(this.pythonPath, [this.aiScript, ...args], { env });

      let stdout = '';
      let stderr = '';

      if (input) {
        process_.stdin.write(input);
        process_.stdin.end();
      }

      process_.stdout.on('data', data => {
        stdout += data.toString();
      });

      process_.stderr.on('data', data => {
        stderr += data.toString();
      });

      process_.on('close', code => {
        if (code !== 0) {
          resolve({
            success: false,
            error: `AI backend failed with code ${code}: ${stderr}`,
          });
          return;
        }

        try {
          const result = JSON.parse(stdout) as AIResponse;
          resolve(result);
        } catch (error) {
          resolve({
            success: false,
            error: `Failed to parse AI response: ${error}`,
            metadata: { raw_output: stdout },
          });
        }
      });

      process_.on('error', error => {
        reject({
          success: false,
          error: `Failed to execute AI backend: ${error.message}`,
        });
      });
    });
  }

  /**
   * Get available models for current provider
   */
  public async getAvailableModels(): Promise<string[]> {
    if (this.config.provider === AIProvider.OpenAI) {
      return ['gpt-4', 'gpt-4-turbo', 'gpt-3.5-turbo'];
    } else if (this.config.provider === AIProvider.Ollama) {
      try {
        const response = await fetch(`${this.config.ollamaUrl}/api/tags`);
        if (!response.ok) {
          return ['llama2', 'codellama', 'mistral'];
        }
        const data = await response.json();
        return data.models?.map((m: any) => m.name) || ['llama2'];
      } catch {
        return ['llama2', 'codellama', 'mistral'];
      }
    }
    return [];
  }

  /**
   * Validate AI configuration
   */
  public validateConfig(): { valid: boolean; errors: string[] } {
    const errors: string[] = [];

    if (!this.config.enabled) {
      errors.push('AI features are disabled');
    }

    if (this.config.provider === AIProvider.OpenAI && !this.config.apiKey) {
      errors.push('OpenAI API key is required');
    }

    if (!this.config.model) {
      errors.push('Model name is required');
    }

    if (this.config.temperature < 0 || this.config.temperature > 2) {
      errors.push('Temperature must be between 0 and 2');
    }

    return {
      valid: errors.length === 0,
      errors,
    };
  }
}

/**
 * AI Core factory for creating instances
 */
export class AICoreFactory {
  private static instance: AICore | null = null;

  static getInstance(context: vscode.ExtensionContext): AICore {
    if (!AICoreFactory.instance) {
      AICoreFactory.instance = new AICore(context);
    }
    return AICoreFactory.instance;
  }

  static reset(): void {
    AICoreFactory.instance = null;
  }
}
