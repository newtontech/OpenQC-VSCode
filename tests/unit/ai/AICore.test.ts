import { AICore, AICoreFactory, AIProvider, OPENAI_API_KEY_SECRET } from '../../../src/ai/AICore';

const mockGetConfiguration = jest.fn();
const mockSpawn = jest.fn();
const secretGet = jest.fn();
const secretStore = jest.fn();
const secretDelete = jest.fn();

jest.mock('vscode', () => ({
  workspace: { getConfiguration: () => mockGetConfiguration() },
}));

jest.mock('child_process', () => ({
  spawn: (...args: unknown[]) => mockSpawn(...args),
}));

function config(overrides: Record<string, unknown> = {}): void {
  const values: Record<string, unknown> = {
    enabled: true,
    provider: AIProvider.Ollama,
    model: 'llama3.2',
    openaiBaseUrl: 'https://api.openai.invalid/v1',
    ollamaUrl: 'http://localhost:11434',
    maxTokens: 128,
    maxOutputChars: 1024,
    timeoutSeconds: 5,
    temperature: 0.2,
    pythonPath: 'python-test',
    ...overrides,
  };
  mockGetConfiguration.mockReturnValue({
    get: jest.fn((key: string, defaultValue?: unknown) => values[key] ?? defaultValue),
  });
}

function childProcess(
  options: {
    stdout?: string;
    stderr?: string;
    code?: number;
    neverClose?: boolean;
  } = {}
) {
  const kill = jest.fn();
  return {
    stdin: { write: jest.fn(), end: jest.fn() },
    stdout: {
      on: jest.fn((event: string, callback: (data: string) => void) => {
        if (event === 'data' && options.stdout !== undefined) {
          callback(options.stdout);
        }
      }),
    },
    stderr: {
      on: jest.fn((event: string, callback: (data: string) => void) => {
        if (event === 'data' && options.stderr !== undefined) {
          callback(options.stderr);
        }
      }),
    },
    on: jest.fn((event: string, callback: (value: number) => void) => {
      if (event === 'close' && !options.neverClose) {
        callback(options.code ?? 0);
      }
    }),
    kill,
  };
}

describe('AICore secure bridge', () => {
  const context = {
    extensionPath: '/extension',
    secrets: { get: secretGet, store: secretStore, delete: secretDelete },
  } as any;

  beforeEach(() => {
    jest.clearAllMocks();
    jest.useRealTimers();
    AICoreFactory.reset();
    config();
    secretGet.mockResolvedValue(undefined);
    secretStore.mockResolvedValue(undefined);
    secretDelete.mockResolvedValue(undefined);
  });

  it('stores and clears the OpenAI credential only through SecretStorage', async () => {
    const core = new AICore(context);
    await core.setOpenAIApiKey('  test-secret  ');
    await core.clearOpenAIApiKey();

    expect(secretStore).toHaveBeenCalledWith(OPENAI_API_KEY_SECRET, 'test-secret');
    expect(secretDelete).toHaveBeenCalledWith(OPENAI_API_KEY_SECRET);
    expect(core.getConfig()).not.toHaveProperty('apiKey');
  });

  it('refreshes settings and validates unsafe bounds', () => {
    const core = new AICore(context);
    expect(core.isEnabled()).toBe(true);
    config({ enabled: false, model: '', temperature: 3, timeoutSeconds: 0, maxOutputChars: 10 });

    core.refreshConfig();

    expect(core.isEnabled()).toBe(false);
    expect(core.validateConfig()).toEqual({
      valid: false,
      errors: [
        'AI features are disabled',
        'Model name is required',
        'Temperature must be between 0 and 2',
        'Timeout must be between 0 and 600 seconds',
        'Maximum output characters must be between 256 and 1000000',
      ],
    });
  });

  it('reports whether SecretStorage contains an OpenAI key', async () => {
    secretGet.mockResolvedValueOnce('credential').mockResolvedValueOnce(undefined);
    const core = new AICore(context);

    await expect(core.hasOpenAIApiKey()).resolves.toBe(true);
    await expect(core.hasOpenAIApiKey()).resolves.toBe(false);
    await expect(core.setOpenAIApiKey('   ')).rejects.toThrow('cannot be empty');
  });

  it('rejects OpenAI requests with a missing SecretStorage credential', async () => {
    config({ provider: AIProvider.OpenAI, model: 'gpt-test' });
    const result = await new AICore(context).explainParameters('input', 'vasp');

    expect(result).toEqual({ success: false, error: 'OpenAI API key is not configured' });
    expect(mockSpawn).not.toHaveBeenCalled();
  });

  it('runs the Python client as a module and passes the secret only in the child environment', async () => {
    config({ provider: AIProvider.OpenAI, model: 'gpt-test' });
    secretGet.mockResolvedValue('secret-value');
    mockSpawn.mockReturnValue(
      childProcess({ stdout: JSON.stringify({ success: true, content: 'ok' }) })
    );

    const result = await new AICore(context).explainParameters('input', 'vasp');

    expect(result).toMatchObject({ success: true, content: 'ok' });
    const [executable, args, options] = mockSpawn.mock.calls[0];
    expect(executable).toBe('python-test');
    expect(args).toEqual(['-m', 'openqc.ai.client', 'explain']);
    expect(options.env.OPENQC_AI_API_KEY).toBe('secret-value');
    expect(JSON.stringify(args)).not.toContain('secret-value');
  });

  it('routes generate, debug, and availability operations to their module commands', async () => {
    mockSpawn.mockImplementation(() =>
      childProcess({ stdout: JSON.stringify({ success: true, generatedInput: 'input' }) })
    );
    const core = new AICore(context);

    await expect(core.generateInput('water', 'cp2k', { structure: 'xyz' })).resolves.toMatchObject({
      success: true,
    });
    await expect(core.debugCalculation('input', 'output', 'orca')).resolves.toMatchObject({
      success: true,
    });
    await expect(core.isAvailable()).resolves.toBe(true);

    expect(mockSpawn.mock.calls.map(call => call[1][2])).toEqual(['generate', 'debug', 'check']);
  });

  it('does not check availability while AI is disabled', async () => {
    config({ enabled: false });
    await expect(new AICore(context).isAvailable()).resolves.toBe(false);
    expect(mockSpawn).not.toHaveBeenCalled();
  });

  it('short-circuits an already-cancelled request', async () => {
    const token = { isCancellationRequested: true } as any;
    await expect(new AICore(context).optimizeInput('input', 'vasp', token)).resolves.toEqual({
      success: false,
      error: 'AI request was cancelled',
    });
    expect(mockSpawn).not.toHaveBeenCalled();
  });

  it('kills a timed-out bridge process', async () => {
    jest.useFakeTimers();
    config({ timeoutSeconds: 0.01 });
    const child = childProcess({ neverClose: true });
    mockSpawn.mockReturnValue(child);

    const pending = new AICore(context).optimizeInput('input', 'vasp');
    await jest.advanceTimersByTimeAsync(100);

    await expect(pending).resolves.toEqual({
      success: false,
      error: 'AI request timed out after 0.01 seconds',
    });
    expect(child.kill).toHaveBeenCalled();
  });

  it('kills the process when VS Code cancellation is requested', async () => {
    let cancel: (() => void) | undefined;
    const token = {
      isCancellationRequested: false,
      onCancellationRequested: jest.fn((callback: () => void) => {
        cancel = callback;
        return { dispose: jest.fn() };
      }),
    } as any;
    const child = childProcess({ neverClose: true });
    mockSpawn.mockReturnValue(child);

    const pending = new AICore(context).optimizeInput('input', 'vasp', token);
    await Promise.resolve();
    cancel?.();

    await expect(pending).resolves.toEqual({ success: false, error: 'AI request was cancelled' });
    expect(child.kill).toHaveBeenCalled();
  });

  it('bounds stdout and does not echo invalid output', async () => {
    config({ maxOutputChars: 256 });
    const child = childProcess({ stdout: 'x'.repeat(5000), neverClose: true });
    mockSpawn.mockReturnValue(child);

    const result = await new AICore(context).optimizeInput('input', 'vasp');

    expect(result.error).toBe('AI backend output exceeded the configured safe limit');
    expect(result).not.toHaveProperty('metadata.raw_output');
    expect(child.kill).toHaveBeenCalled();
  });

  it('redacts the credential from child-process errors', async () => {
    config({ provider: AIProvider.OpenAI });
    secretGet.mockResolvedValue('example-credential');
    mockSpawn.mockReturnValue(
      childProcess({ code: 1, stderr: 'Authorization: Bearer example-credential' })
    );

    const result = await new AICore(context).optimizeInput('input', 'vasp');

    expect(result.error).toContain('[REDACTED]');
    expect(result.error).not.toContain('example-credential');
  });

  it('returns configured OpenAI model and discovers Ollama models', async () => {
    config({ provider: AIProvider.OpenAI, model: 'gpt-test' });
    await expect(new AICore(context).getAvailableModels()).resolves.toEqual(['gpt-test']);

    config({ provider: AIProvider.Ollama });
    const fetchMock = jest
      .spyOn(global, 'fetch')
      .mockResolvedValueOnce({
        ok: true,
        json: async () => ({ models: [{ name: 'local-a' }, {}, { name: 'local-b' }] }),
      } as any)
      .mockRejectedValueOnce(new Error('offline'));
    const core = new AICore(context);

    await expect(core.getAvailableModels()).resolves.toEqual(['local-a', 'local-b']);
    await expect(core.getAvailableModels()).resolves.toEqual([]);
    fetchMock.mockRestore();
  });

  it('reuses and resets the factory singleton', () => {
    const first = AICoreFactory.getInstance(context);
    expect(AICoreFactory.getInstance(context)).toBe(first);
    AICoreFactory.reset();
    expect(AICoreFactory.getInstance(context)).not.toBe(first);
  });
});
