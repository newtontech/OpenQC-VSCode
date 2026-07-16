/** Secure TypeScript bridge to the Python AI module. */

import * as path from 'path';
import { spawn } from 'child_process';
import * as vscode from 'vscode';

export const OPENAI_API_KEY_SECRET = 'openqc.ai.openaiApiKey';
const STDERR_LIMIT = 8_192;

export enum AIProvider {
  OpenAI = 'openai',
  Ollama = 'ollama',
}

export interface AIConfig {
  enabled: boolean;
  provider: AIProvider;
  model: string;
  openaiBaseUrl: string;
  ollamaUrl: string;
  maxTokens: number;
  maxOutputChars: number;
  timeoutSeconds: number;
  temperature: number;
}

export enum AIRequestType {
  OptimizeInput = 'optimize_input',
  GenerateInput = 'generate_input',
  ExplainParameters = 'explain_parameters',
  DebugCalculation = 'debug_calculation',
}

export interface AIRequest {
  type: AIRequestType;
  content: string;
  software?: string;
  context?: Record<string, unknown>;
}

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
  metadata?: Record<string, unknown>;
}

export class AICore {
  private config: AIConfig;
  private pythonPath: string;
  private readonly pythonRoot: string;

  constructor(private readonly context: vscode.ExtensionContext) {
    this.pythonPath = this.getPythonPath();
    this.pythonRoot = path.join(context.extensionPath, 'python');
    this.config = this.loadConfig();
  }

  private getPythonPath(): string {
    return vscode.workspace.getConfiguration('openqc').get<string>('pythonPath', 'python3');
  }

  private loadConfig(): AIConfig {
    const config = vscode.workspace.getConfiguration('openqc.ai');
    return {
      enabled: config.get<boolean>('enabled', false),
      provider: config.get<AIProvider>('provider', AIProvider.Ollama),
      model: config.get<string>('model', 'llama2'),
      openaiBaseUrl: config.get<string>('openaiBaseUrl', 'https://api.openai.com/v1'),
      ollamaUrl: config.get<string>('ollamaUrl', 'http://localhost:11434'),
      maxTokens: config.get<number>('maxTokens', 2048),
      maxOutputChars: config.get<number>('maxOutputChars', 100_000),
      timeoutSeconds: config.get<number>('timeoutSeconds', 120),
      temperature: config.get<number>('temperature', 0.7),
    };
  }

  public refreshConfig(): void {
    this.pythonPath = this.getPythonPath();
    this.config = this.loadConfig();
  }

  public isEnabled(): boolean {
    return this.config.enabled;
  }

  public async setOpenAIApiKey(apiKey: string): Promise<void> {
    const value = apiKey.trim();
    if (!value) {
      throw new Error('OpenAI API key cannot be empty');
    }
    await this.context.secrets.store(OPENAI_API_KEY_SECRET, value);
  }

  public async clearOpenAIApiKey(): Promise<void> {
    await this.context.secrets.delete(OPENAI_API_KEY_SECRET);
  }

  public async hasOpenAIApiKey(): Promise<boolean> {
    return Boolean(await this.context.secrets.get(OPENAI_API_KEY_SECRET));
  }

  public async isAvailable(token?: vscode.CancellationToken): Promise<boolean> {
    if (!this.config.enabled) {
      return false;
    }
    const result = await this.executeAI(['check'], '', token);
    return result.success;
  }

  public getConfig(): AIConfig {
    return { ...this.config };
  }

  public optimizeInput(
    inputContent: string,
    software: string,
    token?: vscode.CancellationToken
  ): Promise<AIResponse> {
    return this.executeAI(
      ['optimize'],
      JSON.stringify({
        type: AIRequestType.OptimizeInput,
        content: inputContent,
        software,
      } satisfies AIRequest),
      token
    );
  }

  public generateInput(
    description: string,
    software: string,
    context?: Record<string, unknown>,
    token?: vscode.CancellationToken
  ): Promise<AIResponse> {
    return this.executeAI(
      ['generate'],
      JSON.stringify({
        type: AIRequestType.GenerateInput,
        content: description,
        software,
        context,
      } satisfies AIRequest),
      token
    );
  }

  public explainParameters(
    inputContent: string,
    software: string,
    token?: vscode.CancellationToken
  ): Promise<AIResponse> {
    return this.executeAI(
      ['explain'],
      JSON.stringify({
        type: AIRequestType.ExplainParameters,
        content: inputContent,
        software,
      } satisfies AIRequest),
      token
    );
  }

  public debugCalculation(
    inputContent: string,
    outputContent: string,
    software: string,
    token?: vscode.CancellationToken
  ): Promise<AIResponse> {
    return this.executeAI(
      ['debug'],
      JSON.stringify({
        type: AIRequestType.DebugCalculation,
        content: inputContent,
        software,
        context: { output: outputContent },
      } satisfies AIRequest),
      token
    );
  }

  private async executeAI(
    args: string[],
    input: string,
    token?: vscode.CancellationToken
  ): Promise<AIResponse> {
    const apiKey =
      this.config.provider === AIProvider.OpenAI
        ? await this.context.secrets.get(OPENAI_API_KEY_SECRET)
        : undefined;
    if (this.config.provider === AIProvider.OpenAI && !apiKey) {
      return { success: false, error: 'OpenAI API key is not configured' };
    }
    if (token?.isCancellationRequested) {
      return { success: false, error: 'AI request was cancelled' };
    }

    return new Promise(resolve => {
      const pathSeparator = process.platform === 'win32' ? ';' : ':';
      const currentPythonPath = process.env.PYTHONPATH;
      const child = spawn(this.pythonPath, ['-m', 'openqc.ai.client', ...args], {
        cwd: this.pythonRoot,
        shell: false,
        env: {
          ...process.env,
          PYTHONPATH: currentPythonPath
            ? `${this.pythonRoot}${pathSeparator}${currentPythonPath}`
            : this.pythonRoot,
          OPENQC_AI_PROVIDER: this.config.provider,
          OPENQC_AI_MODEL: this.config.model,
          OPENQC_AI_API_KEY: apiKey ?? '',
          OPENQC_AI_OPENAI_URL: this.config.openaiBaseUrl,
          OPENQC_AI_OLLAMA_URL: this.config.ollamaUrl,
          OPENQC_AI_MAX_TOKENS: String(this.config.maxTokens),
          OPENQC_AI_MAX_OUTPUT_CHARS: String(this.config.maxOutputChars),
          OPENQC_AI_TIMEOUT_SECONDS: String(this.config.timeoutSeconds),
          OPENQC_AI_TEMPERATURE: String(this.config.temperature),
        },
      });

      let stdout = '';
      let stderr = '';
      let settled = false;
      let exceededOutputLimit = false;
      const outputLimit = Math.max(4096, this.config.maxOutputChars * 4);

      const finish = (response: AIResponse): void => {
        if (settled) {
          return;
        }
        settled = true;
        clearTimeout(timeout);
        cancellation.dispose();
        resolve(response);
      };

      const terminate = (message: string): void => {
        child.kill();
        finish({ success: false, error: message });
      };

      const timeout = setTimeout(
        () => terminate(`AI request timed out after ${this.config.timeoutSeconds} seconds`),
        Math.max(100, this.config.timeoutSeconds * 1000)
      );
      const cancellation = token?.onCancellationRequested(() =>
        terminate('AI request was cancelled')
      ) ?? {
        dispose: () => undefined,
      };

      if (input) {
        child.stdin.write(input);
      }
      child.stdin.end();

      child.stdout.on('data', data => {
        if (settled || exceededOutputLimit) {
          return;
        }
        stdout += data.toString();
        if (Buffer.byteLength(stdout, 'utf8') > outputLimit) {
          exceededOutputLimit = true;
          terminate('AI backend output exceeded the configured safe limit');
        }
      });

      child.stderr.on('data', data => {
        if (stderr.length < STDERR_LIMIT) {
          stderr += data.toString().slice(0, STDERR_LIMIT - stderr.length);
        }
      });

      child.on('close', code => {
        if (settled) {
          return;
        }
        if (code !== 0) {
          const detail = this.sanitizeError(stderr, apiKey);
          finish({
            success: false,
            error: detail
              ? `AI backend failed with code ${code}: ${detail}`
              : `AI backend failed with code ${code}`,
          });
          return;
        }
        try {
          finish(JSON.parse(stdout) as AIResponse);
        } catch {
          finish({ success: false, error: 'AI backend returned an invalid response' });
        }
      });

      child.on('error', error => {
        finish({
          success: false,
          error: `Failed to execute AI backend: ${this.sanitizeError(error.message, apiKey)}`,
        });
      });
    });
  }

  private sanitizeError(message: string, apiKey?: string): string {
    let safe = message;
    if (apiKey) {
      safe = safe.split(apiKey).join('[REDACTED]');
    }
    return safe
      .replace(/bearer\s+[^\s,;]+/gi, 'Bearer [REDACTED]')
      .replace(/\bsk-[A-Za-z0-9_-]{8,}\b/g, '[REDACTED]')
      .trim()
      .slice(0, 500);
  }

  public async getAvailableModels(): Promise<string[]> {
    if (this.config.provider === AIProvider.OpenAI) {
      return this.config.model ? [this.config.model] : [];
    }
    const controller = new AbortController();
    const timeout = setTimeout(() => controller.abort(), this.config.timeoutSeconds * 1000);
    try {
      const response = await fetch(`${this.config.ollamaUrl}/api/tags`, {
        signal: controller.signal,
      });
      if (!response.ok) {
        return [];
      }
      const data = (await response.json()) as { models?: Array<{ name?: string }> };
      return data.models?.flatMap(model => (model.name ? [model.name] : [])) ?? [];
    } catch {
      return [];
    } finally {
      clearTimeout(timeout);
    }
  }

  public validateConfig(): { valid: boolean; errors: string[] } {
    const errors: string[] = [];
    if (!this.config.enabled) {
      errors.push('AI features are disabled');
    }
    if (!this.config.model) {
      errors.push('Model name is required');
    }
    if (this.config.temperature < 0 || this.config.temperature > 2) {
      errors.push('Temperature must be between 0 and 2');
    }
    if (this.config.timeoutSeconds <= 0 || this.config.timeoutSeconds > 600) {
      errors.push('Timeout must be between 0 and 600 seconds');
    }
    if (this.config.maxOutputChars < 256 || this.config.maxOutputChars > 1_000_000) {
      errors.push('Maximum output characters must be between 256 and 1000000');
    }
    return { valid: errors.length === 0, errors };
  }
}

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
