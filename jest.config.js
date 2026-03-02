/** @type {import('ts-jest').JestConfigWithTsJest} */
module.exports = {
  preset: 'ts-jest',
  testEnvironment: 'node',
  roots: ['<rootDir>'],
  testMatch: ['**/tests/unit/**/*.test.ts'],
  moduleNameMapper: {
    '^vscode$': '<rootDir>/tests/mocks/vscode.ts',
  },
  collectCoverageFrom: [
    'src/**/*.ts',
    '!src/**/*.d.ts',
    '!src/test/**',
    '!src/**/index.ts',
    '!src/parsers/**',
    '!src/providers/lsp/**',
    '!src/server.ts',
    '!src/features/**',
    '!src/data/**',
    '!src/extension.ts',
    '!src/managers/LSPManager.ts',
    '!src/providers/DataPlotter.ts',
    '!src/providers/StructureViewer.ts',
    '!src/visualizers/MoleculeViewerPanel.ts',
  ],
  coverageDirectory: 'coverage',
  coverageReporters: ['text', 'text-summary', 'lcov', 'html'],
  coverageThreshold: {
    global: {
      branches: 85,
      functions: 90,
      lines: 90,
      statements: 90,
    },
  },
  moduleFileExtensions: ['ts', 'tsx', 'js', 'jsx', 'json', 'node'],
  setupFilesAfterEnv: ['<rootDir>/tests/setup.ts'],
  verbose: true,
};
