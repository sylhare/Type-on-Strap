import type { Config } from 'jest';

const config: Config = {
  testEnvironment: 'jsdom',
  testMatch: ['**/test/browser/**/*.test.js'],
  roots: ['<rootDir>/test/browser'],
  clearMocks: true,
  resetMocks: true,
  restoreMocks: true,
  verbose: true,
};

export default config;
