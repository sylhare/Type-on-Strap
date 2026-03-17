import type { Config } from 'jest';

const config: Config = {
  testEnvironment: 'node',
  testMatch: ['<rootDir>/test/**/*.test.ts'],
  testPathIgnorePatterns: ['/node_modules/', '/test/browser/', '/test/integration/'],
  transform: {
    '^.+\\.tsx?$': ['ts-jest', { tsconfig: '<rootDir>/tsconfig.json', diagnostics: false }],
  },
  roots: ['<rootDir>/test'],
  clearMocks: true,
  resetMocks: true,
  restoreMocks: true,
  verbose: true,
};

export default config;
