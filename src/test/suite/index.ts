import * as path from 'path';
import { glob } from 'glob';

interface TestCase {
  suite: string;
  name: string;
  run: () => void | Promise<void>;
}

export function run(): Promise<void> {
  const testsRoot = path.resolve(__dirname, '..');
  const tests: TestCase[] = [];
  let currentSuite = '';

  const previousSuite = (globalThis as any).suite;
  const previousTest = (globalThis as any).test;

  (globalThis as any).suite = (name: string, callback: () => void) => {
    const parentSuite = currentSuite;
    currentSuite = parentSuite ? `${parentSuite} ${name}` : name;
    callback();
    currentSuite = parentSuite;
  };

  (globalThis as any).test = (name: string, callback: () => void | Promise<void>) => {
    tests.push({ suite: currentSuite, name, run: callback });
  };

  return new Promise((c, e) => {
    glob('**/*.test.js', { cwd: testsRoot })
      .then(async files => {
        for (const file of files) {
          require(path.resolve(testsRoot, file));
        }

        let failures = 0;
        for (const testCase of tests) {
          try {
            await testCase.run();
            console.log(`PASS ${testCase.suite} - ${testCase.name}`);
          } catch (err) {
            failures += 1;
            console.error(`FAIL ${testCase.suite} - ${testCase.name}`);
            console.error(err);
          }
        }

        if (failures > 0) {
          e(new Error(`${failures} tests failed.`));
        } else {
          c();
        }
      })
      .catch(err => {
        console.error(err);
        e(err);
      })
      .finally(() => {
        (globalThis as any).suite = previousSuite;
        (globalThis as any).test = previousTest;
      });
  });
}
