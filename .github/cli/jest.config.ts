import type { Config } from 'jest';

const config: Config = {
  projects: [
    {
      displayName: 'unit',
      testEnvironment: 'node',
      testMatch: ['<rootDir>/test/**/*.test.ts'],
      testPathIgnorePatterns: ['/node_modules/', '/test/browser/', '/test/integration/'],
      transform: {
        '^.+\\.tsx?$': ['ts-jest', { tsconfig: '<rootDir>/tsconfig.json', diagnostics: false }],
      },
      roots: ['<rootDir>/test'],
    },
    {
      displayName: 'browser',
      testEnvironment: 'jsdom',
      testMatch: ['<rootDir>/test/browser/**/*.test.js'],
      roots: ['<rootDir>/test/browser'],
      transform: {},
    },
    {
      displayName: 'integration',
      testEnvironment: 'node',
      testMatch: ['<rootDir>/test/integration/**/*.test.js'],
      roots: ['<rootDir>/test/integration'],
      transform: {
        '^.+\\.tsx?$': ['ts-jest', { tsconfig: '<rootDir>/tsconfig.json', diagnostics: false }],
      },
    },
  ],
  clearMocks: true,
  resetMocks: true,
  restoreMocks: true,
};

export default config;
