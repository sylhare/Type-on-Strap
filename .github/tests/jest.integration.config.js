/**
 * Jest configuration for Type-on-Strap integration tests.
 * Uses real file I/O and real npm dependencies — requires `npm install` at project root.
 */
module.exports = {
  testEnvironment: 'node',
  testMatch: ['**/integration/**/*.test.js'],
  roots: ['<rootDir>/integration'],
  clearMocks: true,
  verbose: true,
};
