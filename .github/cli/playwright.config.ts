import { defineConfig, devices } from '@playwright/test';
import * as path from 'node:path';

/**
 * Playwright configuration for Type-on-Strap theme e2e tests
 * @see https://playwright.dev/docs/test-configuration
 */
export default defineConfig({
  testDir: './e2e',
  timeout: 30 * 1000,
  fullyParallel: true,
  forbidOnly: !!process.env['CI'],
  retries: process.env['CI'] ? 2 : 0,
  workers: process.env['CI'] ? 1 : undefined,
  reporter: [
    ['html', { outputFolder: 'playwright-report' }],
    ['list'],
  ],
  use: {
    // Base URL for tests
    // Note: Site is served from _site/ root, but links contain /Type-on-Strap/ prefix
    baseURL: process.env['BASE_URL'] ?? 'http://localhost:4000',
    trace: 'on-first-retry',
    screenshot: 'only-on-failure',
    video: 'retain-on-failure',
  },
  expect: {
    toHaveScreenshot: {
      maxDiffPixelRatio: 0.02,
    },
  },
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] },
      grep: /@desktop/,
    },
    {
      name: 'firefox',
      use: { ...devices['Desktop Firefox'] },
      grep: /@desktop/,
    },
    {
      name: 'webkit',
      use: { ...devices['Desktop Safari'] },
      grep: /@desktop/,
    },
    // Mobile tests
    {
      name: 'mobile-chrome',
      use: { ...devices['Pixel 5'] },
      grep: /@mobile/,
    },
    {
      name: 'mobile-safari',
      use: { ...devices['iPhone 12'] },
      grep: /@mobile/,
    },
    {
      name: 'visual-chromium',
      use: { ...devices['Desktop Chrome'] },
      grep: /@visual/,
    },
    {
      name: 'visual-firefox',
      use: { ...devices['Desktop Firefox'] },
      grep: /@visual/,
    },
    {
      name: 'visual-webkit',
      use: { ...devices['Desktop Safari'] },
      grep: /@visual/,
    },
    {
      name: 'visual-mobile-chrome',
      use: { ...devices['Pixel 5'] },
      grep: /@visual/,
    },
    {
      name: 'visual-mobile-safari',
      use: { ...devices['iPhone 12'] },
      grep: /@visual/,
    },
  ],
  webServer: {
    command: 'bundle exec jekyll build --baseurl "" --quiet && npm run server',
    cwd: path.resolve(__dirname, '../..'),
    url: 'http://localhost:4000',
    reuseExistingServer: !process.env['CI'],
    timeout: 120 * 1000,
  },
});
