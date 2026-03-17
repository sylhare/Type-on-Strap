import type { Config } from 'jest';

/**
 * Jest configuration for Type-on-Strap integration tests.
 * Uses real file I/O and real npm dependencies — requires `npm install` at project root.
 */
const config: Config = {
  testEnvironment: 'node',
  testMatch: ['**/test/integration/**/*.test.js'],
  roots: ['<rootDir>/test/integration'],
  clearMocks: true,
  verbose: true,
};

export default config;
