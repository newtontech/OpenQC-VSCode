import * as path from 'path';

export interface VisualSmokeOptions {
  workspaceRoot: string;
  extensionRoot: string;
  outputDir: string;
  outputPrefix: string;
}

export interface VisualSmokePaths {
  htmlPath: string;
  screenshotPath: string;
  narrowScreenshotPath: string;
  reportPath: string;
}

function readOptionValue(args: string[], index: number, flag: string): string {
  const value = args[index + 1];
  if (!value || value.startsWith('--')) {
    throw new Error(`Missing value for ${flag}`);
  }
  return value;
}

export function parseVisualSmokeArgs(args: string[], workspaceRoot: string): VisualSmokeOptions {
  let extensionRoot = workspaceRoot;
  let outputDir = path.join(workspaceRoot, 'output', 'playwright');
  let outputPrefix = 'openqc-viewer-smoke';

  for (let index = 0; index < args.length; index += 1) {
    const arg = args[index];
    switch (arg) {
      case '--extension-root':
        extensionRoot = readOptionValue(args, index, arg);
        index += 1;
        break;
      case '--output-dir':
        outputDir = readOptionValue(args, index, arg);
        index += 1;
        break;
      case '--output-prefix':
        outputPrefix = readOptionValue(args, index, arg);
        index += 1;
        break;
      default:
        throw new Error(`Unknown visual smoke option: ${arg}`);
    }
  }

  return {
    workspaceRoot,
    extensionRoot: path.resolve(extensionRoot),
    outputDir: path.resolve(outputDir),
    outputPrefix,
  };
}

export function buildVisualSmokePaths(options: VisualSmokeOptions): VisualSmokePaths {
  return {
    htmlPath: path.join(options.outputDir, `${options.outputPrefix}.html`),
    screenshotPath: path.join(options.outputDir, `${options.outputPrefix}.png`),
    narrowScreenshotPath: path.join(options.outputDir, `${options.outputPrefix}-narrow.png`),
    reportPath: path.join(options.outputDir, `${options.outputPrefix}.json`),
  };
}
